"""Class handles xyz stage communications.  Commands are written on
   the serial port in binary (e.g. b"x movr 5\n") Use requires
   installation of "Pyserial" module ("pip install pyserial")

"""

import serial
import time

X, Y, Z = 'X', 'Y', 'Z'
REL, ABS = 'REL', 'ABS'

ZMIN = 0.2
ZMAX = 9.8

class serial_dummy(object):
    def __init__(self):
        pass

    def read(self, *args, **kwargs):
        return ''

    def write(self, *args, **kwargs):
        pass

    def readline(self, *args, **kwargs):
        return ''


class d3000(object):
    """This class controls the xyz stage motors
    """
    
    def __init__(self, ctrl_port, debug=False, pbar=None):
        """Initialize communication with stage
        """
        self.enabled = False
        
        if pbar:
            pbar.set_status('Scanning for a valid stage COM port...')

        # find the correct COM port
        self.find_COM(b'X PRINT DN\n', b'"X"\r\n', ctrl_port, pbar=pbar)
        
    def limit_z(self, pos, span=0):
        return min(ZMAX-span/2, max(ZMIN+span/2, pos))

        
    def en_flag(self):
        return self.enabled

    def move(self, axis, x, move_type='ABS', wait=True, timeout=10):
        """Args:

            axis: the axis id (X, Y, Z)

            x: the distance to move or the final destination in mm

            move_type: absolute ('ABS') or relative ('REL')

            wait: a flag to wait for move complete. Defaults to True.

            timeout: maximum wait time in s, set to -1 to never
                expire. Defaults to 10.

        """
        if not self.enabled:
            return
        
        # Motors are enabled
        cmd = 'MOVR' if move_type == 'REL' else 'MOVA'
        self.send_cmd(axis, cmd+' {:.3f}'.format(x))

        # Wait for motor to stop and check timeout
        if wait:
            self.wait(axis, timeout)

            
    def wait(self, axis, timeout=10):
        """Waits for the motor relative to 'axis' to stop moving
           Returns a flag: TRUE/FALSE if the waiting timed-out
           
        """
        if not self.enabled:
            return

        # Wait for the motor to stop moving
        moving = True
        seconds = int(round(time.time() * 1000))
                
        # check moving flag
        while moving:
            time.sleep(0.01)
            flags = self.status(axis)
            if (flags[0] and flags[1])==True:
                moving = False
                return False
            else:                 # Timeout
                moving = True
                if timeout == -1:
                    pass
                elif (int(round(time.time() * 1000))-seconds)/1000 > timeout:
                    return True
        


    def home(self, hard=False, pbar=None):
        """Home the motor"""
        if not self.enabled:
            return

        for name in list(['X', 'Y', 'Z']):
            if pbar is not None:
                pbar.set_status('Homing axis '+name)
            
            self.home_axis(name)
            self.set_position(name, 0)
    
    
    def get_position(self, axis):
        """Get the current position"""
        if not self.enabled:
            return

        self.send_cmd(axis, ' PRINT POS')
        return float(self.ser.read(100))
    

    def set_position(self, axis, x):
        """Set the absolute reference frame"""

        if not self.enabled:
            return

        self.send_cmd(axis, ' POS={:.3f}'.format(x))
        return float(self.get_position(axis))
    
    
    def status(self,axis):
        """check the current status"""
        
        if not self.enabled:
            return (False,False)
            
        enabled = True
        self.send_cmd(axis, ' PRINT MVG')

        flag = self.ser.read(100)
        moving = True
        
        if flag[:4] == b'FALS':                             
            moving = False
        elif flag[:4] == b'TRUE':
            moving = True

        non_moving = not moving
        return (enabled, non_moving)
        

    def initialize(self, pbar=None):
        """Initialize myself"""
        self.enabled = False

        try:                                              
            self.enabled = True

            if pbar is not None:
                pbar.set_status('Adjusting stage acceleration ...')

            self.set_accel(slow=False)

            if pbar is not None:
                pbar.set_status('Homing the stage ...')

            self.home(hard=True, pbar=pbar)
            
        except:
            pass
                

    def set_accel(self, slow=False):
        """Set the stage acceleration factors. For the X and Y axis the
        default value is 2540. For the Z axis, it is 762.

        Args:
            slow: choose to slow down the sage 10X. Default is False.

        """
        
        if self.enabled:
            for axis,accel in zip(['X', 'Y', 'Z'], [2540, 2540, 762]):
                sf = (12 if axis == 'Z' else 10) if slow else 2
                for key,val in zip(['ACCL', 'DECL', 'ACLT', 'DCLT', 'LDCLT'],
                                   [accel/sf, accel/sf, 1, 1, 1]):
                    self.send_cmd(axis, ' {:s}={:.0f}'.format(key, val))
                    self.ser.read(100)

    
    def find_COM(self, cmd, name, ctrl_port, pbar=None):
        """Find COM port where instruments responds to command 'cmd' with
        'name' string

        """
        
        for self.port in [ctrl_port] + ['COM%s' % (i + 1) for i in range(256)]:

            try:
                # open the port
                self.ser = serial.Serial(self.port, baudrate=9600, timeout=0.2)

                # send an identifier command
                self.ser.write(cmd)
                time.sleep(0.2)

                # read the response
                resp = self.ser.read(100)
                
                if resp[:len(name)] == name:
                    # run initialization if a valid port was open
                    self.initialize(pbar=pbar)
                    break
                
                else:
                    self.ser.close()
                
            except:
                pass
    
    
    def close(self, pbar=None):
        """Handle shuttind down the device."""

        if not self.enabled:
            return
        
        # home the axis
        self.home(pbar=pbar)
        
        # close the communication and free the port
        self.ser.close()
        

    def send_cmd(self, axis, cmd, delay=0.1):
        self.ser.write(bytes(axis, 'utf-8'))
        self.ser.write(bytes('{:s}\n'.format(cmd), 'utf-8'))
        time.sleep(delay)

    def home_axis(self, axis, timeout=30):
        self.send_cmd(axis, 'IOS 21=12')
        self.send_cmd(axis, 'FIOS -25,+0.5,21')
        
        moving, stime = True, time.time()
        while(moving):
            self.send_cmd(axis, ' PRINT MVG')
            
            if self.ser.read(100)[:4] == b'FALS' \
               or time.time() > stime+timeout:
                break
