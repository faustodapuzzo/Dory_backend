
# coding: utf-8

# In[10]:

"""Class handles Snowy Range communications.
"""

import serial
import time
import numpy as np

ON, OFF = 'ON', 'OFF'
OPEN, CLOSE = 'OPEN', 'CLOSE'
ID, RASTER, REF, PWR, SHUTTER, DATA = \
            'ID', 'RASTER', 'REFERENCING', 'LASERPWR', 'SHUTTER', 'READ'

# different dalays (sec) for different actions
mech, non_mech = 0.4, 0.4


class serial_dummy(object):
    def __init__(self):
        pass

    def read(self, *args, **kwargs):
        return ''

    def write(self, *args, **kwargs):
        print(args)

    def readline(self, *args, **kwargs):
        return ''

    
class snri(object):
    """This class does controls the Snowy Range Natinal Instrument
    """
    
    def __init__(self, ctrl_port, debug=False, timeout=20, pbar=None):
        """Handle the communication with Snowy Range.

        Args: 
            sett: the dictionary with default settings

            debug: a flag to debug the serial communication. Defaults
                to False.

            timeout: the timeout for the serial read in s. Defaults to
                20.

            pbar: a pointer to a progress bar

        """
        self.debug = debug
        self.timeout = 1 if self.debug else timeout
        self.coll_time = 1
        self.pbar = pbar
        self.enabled = False
        
        # default settings 
        self.setts = {'INT': 0.3,  'AVE': 1, 'RAS': OFF, 'REF': ON,
                      'BAS': OFF, 'PLV': 3, 'PWR': OFF, 'SHU': CLOSE}

        # attempt to open the COM port
        if pbar is not None:
            pbar.set_status('Scanning for a valid spectrometer COM port...')
        self.find_COM(b'id\n', b'\r\n Sierra 2',ctrl_port)

        # initialize the spectrometer
        if pbar is not None:
            pbar.set_status('Initializing the spectrometer ...')
        self.settings(**self.setts)
        
 
    def en_flag(self):
        return self.enabled
    
    def settings(self, **opts):
        """Updates the setting in the settings dictionary "self.setts"
            and in the hardware
            
        """
        if not self.enabled:
            return
        
        for key in opts.keys():
            # ignore the power on/off for now
            if key == 'PWR':
                continue
            
            # update the hardware settings
            if key == 'INT':
                self.int_time(opts[key])
            elif key == 'RAS':
                self.raster(opts[key])
            elif key == 'REF':
                self.ref(opts[key])
            elif key == 'BAS':
                self.baseline(opts[key])
            elif key == 'PLV':
                self.pwr(None, opts[key])
            elif key == 'SHU':
                self.shutter(opts[key])
            elif key == 'AVE':
                self.ave(opts[key])
                
            # update the self.setts dictionary
            self.setts[key] = opts[key]

        # turn the laser on at the very end
        if 'PWR' in opts.keys():
            self.pwr(opts['PWR'], None)
            self.setts['PWR'] = opts['PWR']

            
    def int_time(self, sec= ''):
        """Set the integration time in seconds
        
        Args:
            sec: a floating point number for seconds of integration.
                 default: '' to quarry the int_time status
            
        Returns:
            True|False if the integration time was valid.
            
        """
        if not self.enabled:
            return
        
        if sec == '':
            self.ser.write(b'int_time' + bytes(str(sec), 'utf-8') + b'\n')
            time.sleep(mech)
            msg = self.read()
            if msg != b'':
                print('Timing:' + str(msg))   
            return True
        elif (sec >= 0.1 and sec <= 30):
            try:
                self.ser.write(b'int_time set ' + bytes(str(sec), 'utf-8') \
                               + b'\n')
                time.sleep(mech)
                msg = self.read()

            except:
                pass

            return True
        
        return False

    
    def ave(self, ave=1):
        """Set averages number
        """
        if not self.enabled:
            return

        self.ser.write(b'averages ' + bytes(str(ave), 'utf-8') + b'\n')
        time.sleep(non_mech)
        try:
            msg = self.read()
            if msg != b'':
                print('Averages: ' + str(msg))
        except:
            pass
        
    def raster(self, status=''):
        """Control the Raster.
        Args:
            status:  turn it ON/OFF with True/False
            default: '' to quarry the raster status
            
        """
        if not self.enabled:
            return

        self.ser.write(b'raster ' + bytes(status, 'utf-8') + b'\n')
        time.sleep(mech)
        msg = self.read()
        if msg != b'':
            print('Raster:' + str(msg))       

            
    def pwr(self, status=None, level=None):
        """Set the laser power level and ON/OFF state
        """
        
        if not self.enabled or (status is None and level is None):
            return

        # get the power level requested
        plev = status if status is not None else str(level)

        # send the request to the device
        self.ser.write(b'laserpwr ' + bytes(plev, 'utf-8') + b'\n')

        # read the return
        msg = self.read()
        if msg != b'' and self.debug:
            print('Power Level:' + str(msg) + '['+plev+']')

            
    def ref(self, status=''):
        """Set the referencing status ON/OFF
            default = '' to quarry referencing status
        """
        if not self.enabled:
            return
        
        self.ser.write(b'referencing ' + bytes(status, 'utf-8') + b'\n')
        time.sleep(non_mech)
        msg = self.read()
        if msg != b'':
            print('Referencing:' + str(msg))
        
        return
    
    def baseline(self, status=''):
        """Set the baseline mode ON/OFF
            default = '' to quarry baselining status
        """
        if not self.enabled:
            return

        self.ser.write(b'baselining ' + bytes(status, 'utf-8') + b'\n')
        time.sleep(non_mech)
        msg = self.read()
        if msg != b'':
            print('Baselining:' + str(msg))
        
        return
    
    def shutter(self, status=''):
        """Set the shutter mode OPEN/CLOSE
            default = '' to quarry shutter status
        """
        if not self.enabled:
            return

        try:
            self.ser.write(b'shutter ' + bytes(status, 'utf-8') + b'\n')
            time.sleep(2*mech)
            msg = self.read()
            if msg != b'':
                print('Shutter:' + str(msg))
        except:
            pass
        
        return

    
    def read(self):                     
        """Reads from the spectrometer the status of the parameter 
            Parameter: ID, RASTER, PWR, REF, SHUTTER, SPECTRA
        """
        if not self.enabled:
            return
        
        complete = False
        attempts = 1
        line = []
        length =0
        lettura = 'empty'

        start_time = time.time()
        
        while not complete and time.time() < start_time+self.timeout:
            s=self.ser.readline()     # read until line empty
            if s==b'':                      
                complete = True
            else:
                line.append(s)
                complete = False
                attempts += 1

        if attempts > 0:             # remove useless start and end characters
            if len(line) and (line[0] == b'\x06\r\n' or line[0] == b'\r\n'):
                if line[-1] == b'>' and line[-2] == b'\r\n':
                    lettura = b''.join(line[1:-2])
                    length = len(line[1:-2])
            else:
                lettura = b''.join(line)
                length = len(line)

        return lettura
    
        
    def acquire(self, flush=False, count=1801, data_seg_len=3607):
        """starts acquisition
           Parameter: ID, RASTER, PWR, REF, SHUTTER, SPECTRA
        """
        
        if not self.enabled:
            return np.zeros(count)

        self.ser.write(b'start\n')

        # wait for reference to be measured and read
        if self.setts['REF'] == ON:
            time.sleep(self.setts['AVE']*self.setts['INT'])
            time.sleep(self.non_mech)
            
        # wait for integration to finish
        self.coll_time = self.setts['AVE']*self.setts['INT']
        time.sleep(self.coll_time)
        self.ser.read(100)

        # return the spectrum
        return self.read_spectrum(count, data_seg_len, self.timeout)

    
    def read_spectrum(self, count=1801, data_seg_len=3607, timeout=20):
        # request the data from the device
        self.ser.write(b'read\n')

        # accumulate the data
        data, start_time  = b'', time.time()
        while time.time() < start_time+timeout:
            data = self.ser.read(5000)

            if len(data) < data_seg_len:
                self.ser.write(b'read\n')
                time.sleep(non_mech)            
            else:
                break

        # get the values from the buffer
        if len(data) == data_seg_len:
            try:
                return np.array([(data[2*i] << 8) + data[(2*i)+1]
                                 for i in range(count)])
            except:
                pass

        # if all fails, return zeros
        return np.zeros(count)

        
    def find_COM(self, cmd, name, ctrl_port):
        """Find COM port where instruments responds to command 'cmd' with
        'name' string

        """

        for self.port in [ctrl_port] + ['COM%s' % (i+1) for i in range(256)]:

            try:
                # attempt to open a COM port
                self.ser = serial.Serial(self.port, baudrate=9600, timeout=0.2)

                # try to get the instrument ID
                self.ser.write(cmd)     
                time.sleep(0.2)

                # get the response
                instrument = self.ser.read(100)

                if instrument[:len(name)] == name:
                    self.enabled = True
                    break
                else:
                    self.ser.close()
                
            except Exception as e:
                pass


        
    def close(self):
        """Handle shuttind down the device."""
    
        if not self.enabled or self.debug:
            return
        
        closing_par= {'RAS': OFF, 'PLV': 0, 'PWR': OFF, 'SHU': CLOSE} 
        
        # just update with instance.settings(INT=3, BAS=ON,... )
        self.settings(**closing_par)
        
        self.ser.close()
        self.enabled = False

