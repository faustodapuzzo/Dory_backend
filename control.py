"""Controls the SIERRA spectrometer and the XYZ stage to perform
   integrated function such as autofocus, samples finding and scanning
   
   Requires: loading of the 'Spectrometer' and 'stage_xyz' classes
             calls the functions: hexLattice, gauss

"""

import numpy as np
from numpy import linspace, array, exp, cos, sin, pi, zeros, linalg
from numpy.fft import fftshift, ifftshift, fftfreq, fft, ifft
from numpy import real, exp, pi
from scipy.optimize import curve_fit

from logan.stage import d3000
from logan.spectrometer import snri
from ui_plot import update_plot, ui_plot_graph
from suzy_constants import STAGE_CMD_DELAY, WAVENUMBER_RANGE
from prog_bar import prog_window

import rampy
import time

ON, OFF, OPEN, CLOSE = 'ON', 'OFF', 'OPEN', 'CLOSE'
X, Y, Z, ABS, REL  = 'X', 'Y', 'Z', 'ABS', 'REL'


#########################################################################
##                       MATH FUNCTIONS                                ##
#########################################################################

def hexLattice(center, area, spot = 2, err = 0.1):
    """Make lattice center: center of area of interest size: Return: list
    of positions in x,y that follow an hexagonal lattice in spiral
    order

    """
    # set up Hexagonal lattice constants
    L=3**0.5*spot/2
    h=3*spot/4
    r=spot/2
    
    # finds max value of miller index necessary to cover area
    n = int((area[0]-L/(3**0.5))/h + 1 - err*h) + 1
    m = int((area[1]-L/2)/L + 1 - err*L) + 1
    
    s = zeros((n*m+1,2))
    integral = 0
    
    for w in range(10):
        
        # keep track of previous cycle indexes
        if w==0:
            tot=0
        else:
            integral += -8*(w-1)
            tot = 2*(w)*n + 2*(w)*m + integral - 4*(w)
        
        # cycle through each edge of the perimeter
        # top edge:
        for k in range(n-2*w):
            # start from top right corner 
            if w==0 and k==0:
                # set up first point 
                s[0,0] = h*(n-1)
                s[0,1] = L*(m-1) + (L/4)*(1+(-1)**(n-1))
                
                # transform coordinates to be centered
                s[0,0] += center[0] - h*(n-1)/2
                s[0,1] += center[1] - (L*(m-1)/2 + L/4) 
                
            else:
                try:
                    s[k+tot,0] = s[k+tot-1,0] - h
                    s[k+tot,1] = s[k+tot-1,1] + (L/2)*((-1)**(k+n+w+1))
                except:
                    pass
        # left edge
        for l in range(m-2*w-1):
            try:
                s[l+n-2*w+tot,0] = s[l+n-2*w+tot-1,0] 
                s[l+n-2*w+tot,1] = s[l+n-2*w+tot-1,1] - L
            except:
                pass
                
        # bottom edge
        for k in range(n-2*w-1):
            try:
                s[k+n+m-4*w+tot-1,0] = s[k+n+m-4*w+tot-2,0] + h
                s[k+n+m-4*w+tot-1,1] = s[k+n+m-4*w+tot-2,1] \
                                       + (L/2)*((-1)**(k+w+1))   
            except:
                pass
                
        # right edge
        for l in range(m-2*w-2):
            try:
                s[l+2*n+m-6*w+tot-2,0] = s[l+2*n+m-6*w+tot-3,0] 
                s[l+2*n+m-6*w+tot-2,1] = s[l+2*n+m-6*w+tot-3,1] + L
            except:
                pass
    
    spiral = s[:]
    spiral = spiral[::-1]
    
    return  spiral[1:]

def gauss(x, *p):
    """Gaussian function with offset for fitting
       param: 
           B: vertical offset
           A: gaussian amplitude
           mu: mean
           sigma: s.d.
    """
    B, A, mu, sigma = p
    return B + A*exp(-(x-mu)**2/(2.*sigma**2)) 

################################################################
##                   scanSierra CLASS                         ## 
################################################################

class scanSierra(object):
    """This class controls both the SIERRA and XYZ stage
       and allows to perform autofocus, sample finding, 
       automated data collection over a whole caddy

    Example:
        hpl = scanSierra()
        if hpl.enabled:
           hpl.find_sample(A1)
    
    """
    
    def __init__(self, parent, origin, ai_ports, debug=False,
                 stage_cmd_delay=STAGE_CMD_DELAY, prog=None):
        """Starts communication with both instruments and homes the stage

        Args:
            origin: 3D xyz array for center of A1 caddy element

            prog: the pointer to a progress bar

        """
        # save the pointer to the parent
        self.parent = parent

        
        # delay between moving commands
        self.stage_cmd_delay = stage_cmd_delay
        
        # enabled flags: [spectr., stage]
        self.enabled = [False, False]
        self.ports = ai_ports

        if prog:
            prog.start(2)
            prog.set_status('Initializing the spectrometer')
        
        try:
            # Initialize Spectrometer
            self.snri = snri(self.ports[0], debug=debug, pbar=prog)
            if self.snri.en_flag():
                  self.enabled[0] = True
                  self.ports[0] = self.snri.port
                  
        except Exception as e:
            self.enabled[0] = False

        if prog:
            prog.update(1)
            prog.set_status('Initializing the stage')
            
        try:
            # Initialize Stage
            self.stage = d3000(self.ports[1], debug=debug, pbar=prog)
            if self.stage.en_flag():
                self.enabled[1] = True
                self.ports[1] = self.stage.port
        except Exception as e:                 
            self.enabled[1] = False

        # default settings for Spectrometer
        self.setts = {}        
        self.set_cell_pos(origin)

        
    def set_cell_pos(self, origin, spacing=12.75):
        # compute coordinate for 3x3 Caddy holder
        caddy_names = ['A1', 'A2', 'A3', 'B1', 'B2', 'B3', 'C1', 'C2', 'C3']
        
        self.coordinates = {}
        
        for j in range(3):
            for i in range(3):
                self.coordinates[caddy_names[i+3*j]] = [origin[0]-i*spacing,
                                                        origin[1]+j*spacing,
                                                        origin[2]]
        
        self.coordinates['Load'] = [5, 85, 1]
        self.coordinates['Close'] = [50, 0, 0]
        

        
    def goto_caddy(self, pos, wait_on_move=True, z_dim=False):
        """move caddy element under laser port
           or moves stage to loading position
           (for motor control use control.stage.move(...) function)

        Args: 
           pos: caddy element label 'A1', ... 'C3' or 'Load'
        
           wait_on_move: set wether to wait for the move to
               complete. Defaults to True.
           
           z_dim: a flag to move or not the z dimension. 

        Example:
           object.goto_caddy('B2')

        """
        if not self.enabled[1]:
            return
        
        self.stage.move(X, self.coordinates[pos][0], ABS, wait_on_move)

        time.sleep(self.stage_cmd_delay)
        self.stage.move(Y, self.coordinates[pos][1], ABS, wait_on_move)

        if z_dim:
            time.sleep(self.stage_cmd_delay)
            self.stage.move(Z, self.coordinates[pos][2], ABS, wait_on_move)

            
    def settings(self, **opts):
        """Update settings for Spectrometer

        Arg: PWR, PLV, RAS, BAS, REF, SHU

        Example:
             hpl.settings(PWR=ON, PLV=2, RAS=ON, SHU=OPEN)
        """

        if not self.enabled[0]:
            return

        nopts = {}
        for key in opts.keys():
            if (key not in self.setts or
                self.setts[key] != opts[key] or key == 'PWR'):
                nopts[key] = opts[key]
        
        # setup the instrument
        self.snri.settings(**nopts)

        # keep track of the current settings        
        self.setts.update(opts)

        
    def scan_square(self, element='A1', N=1, area=[10,10], spot=1.6,
                    err=0.1, thresh=1.15e-4, pbar=None):
        """Finds a number N of samples in the caddy element indicated

        Args:

            element:  element label for 3x3 caddy

            N:        number of samples present in caddy element

            area:     size of rectangle to span

            spot:     detection range with raster on (1.6-1.8 mm)

            err:      tollerance on coverage at the edge of the area

            thresh: trigger sample detection for fom above this value
                (bg fom depends on power Level)
            
        Returns:
            array of positions of samples found (up to N)
        
        Example: 
            hpl.scan_square('A3',1)

        """
        
        if (not self.enabled[0]) or (not self.enabled[1]):
            return
        
        # center of caddy element
        center = self.coordinates[element]
        self.goto_caddy(element)
        
        # Generate list of x,y  coordinates to probe
        pos = hexLattice(center, area, spot, err)
        loc, found  = zeros((N,2)), 0

        if pbar is not None:
            pbar.start(len(pos))
        
        # scan caddy locations
        for i in range(len(pos)):
            # compute the distance to all previous locations
            d = [linalg.norm(pos[i]-loc[k,:]) for k in range(found)]

            if len(d) and min(d) < 2:
                self.fom_snri(pos[i], False)
                continue
            
            fom_p, spec = self.fom_snri(pos[i])

            if pbar is not None:
                self.update_master_plot(pbar.parent, spec)
            
            if fom_p > thresh:
                
                if found and len(d) and min(d) < 2:
                    # refine the sample location vector
                    ind = d.index(min(d))
                    loc[ind,:] = (loc[ind,:]+pos[i])/2
                    
                elif found < N:
                    # new sample found
                    loc[found,:] = pos[i]
                    found += 1

                
            # stop here if all requested samples are found
            if found == N:
                break

            if pbar is not None and pbar.is_aborted():
                break
            elif pbar is not None:
                pbar.set_status('Looking for sample {:.0f} out of {:.0f}'\
                                .format(found+1, N))
                pbar.update(i+1)

            
        return loc[:found,:], found

    
    def tuning(self, pos_in, delta=0.3, leap=1.5):
        """Optimizing sample position before profile measurement
        """
        if (not self.enabled[0]) or (not self.enabled[1]):
            return
        
        # calculate gradient sampling points delta mm away from
        # starting pos
        delta = 0.3
        
        # move leap times delta in the gradient direction
        leap = 1.5  
        
        #============= PHASE II ================
        # compute neighbors points for gradient
        dxy = array([pos_in + delta*array([cos(pi*x/2), sin(pi*x/2)]) \
                     for x in range(4)])
        foms = zeros(4)

        for i in range(len(dxy)):
            fom_x = self.fom_snri(dxy[i])
            foms[i] = fom_x[0]

        # find gradient vector
        grad_x = -(foms[0]-foms[2])/(2*delta)
        grad_y = -(foms[1]-foms[3])/(2*delta)
        vector = array([grad_x, grad_y])
        versor = vector/linalg.norm(vector)

        # estimated best next point
        return pos_in + versor*delta*leap

    
    def find_samples(self, element='A1', N=1, area=[10,10], spot=1.6,
                     err=0.1, thresh=1.15e-4, span=3, delta=0.4,
                     z_span=1, pbar=None):
        """Given a caddy element label, it will find N samples and optimize
           signal and acquire Returns an array of spectra and 1D array
           of wavenumbers

        """
        if (not self.enabled[0]) or (not self.enabled[1]):
            return

        # find rough positions 
        pos_a, N_eff = self.scan_square(element, N, area, spot, err,
                                        thresh, pbar=pbar)

        spectra = []

        if pbar is not None:
            pbar.start(N_eff)
        
        # optimize the location of each sample found
        for i in range(N_eff):
            if pbar is not None:
                pbar.set_status('Optimizing sample {:d} of {:d}'\
                                .format(i+1, N_eff))
            
            # pick a place to start the profile
            cent_x, cent_y = self.tuning(pos_a[i])

            self.scan_along_axis(X, cent_x, span, delta, pbar=pbar)
            self.scan_along_axis(Y, cent_y, span, delta, pbar=pbar)

            # only scan z once
            if i == 0:
                z_in = self.stage.limit_z(self.stage.get_position('Z'),
                                          z_span)
                self.scan_along_axis(Z, z_in, z_span, delta/3, pbar)

            # we are now at the best position, grab a spectrum here
            self.parent.setup_device()
            spectra.append(self.acquire(flush=True))
            self.parent.setup_optimize()

            if pbar is not None:
                pbar.update(i+1)

        return spectra


    def scan_along_axis(self, axis, cent, span, delta, pbar=None):
        if pbar is not None:
            pbar.set_status('Scanning axis '+axis+' ...')

        self.stage.move(axis,
                        self.profile(axis, cent, span, delta),
                        ABS, True)
    
        
    def fom_snri(self, pos, measure_spectrum=True):
        """Measures the FOM at a particular position pos (2D or 3D list)
           Returns:
               FOM estimate and Spectrum
        """
        if not self.enabled[0]:
            return
        
        # move in 2D or 3D
        for i,axis in enumerate(['X', 'Y', 'Z']):
            if i < len(pos):
                self.stage.move(axis, pos[i], ABS, False) 
    
        # wait for motion complete
        self.stage.wait(X)
        self.stage.wait(Y)
        if len(pos) == 3:
            self.stage.wait(Z)

        if measure_spectrum:
            # MEASURE spectrum
            spectrum = self.snri.acquire()
    
            # compute FOM
            fom = self.fom(spectrum)

        else:
            fom, spectrum = 0, []
        
        return fom, spectrum


    def gfunc(self, x, a, sig):
        """This is a non-normalized zero-mean Gaussian function.
        
        Args:
            x: a NumPy array with the the x-axis values
            a: the amplitude of the max
            sig: the sigma of the Gaussian distribution

        Returns:
            A NumPy array of the Gaussian y-values per each x-value.

        """
        return a*np.exp(-1*(x**2)/(2*sig**2))

    def highPassFFT(self, y, fmin):
        """This function performs the high pass FFT filter on the input
        data. In this process, the data is zero-mean'ed and then
        passed to the FFT.

        Args:
            fmin: the cutoff frequency.

        """
        y = np.array(y)
        xr = np.arange(y.shape[0])
        
        # remove the DC component
        sfft = np.fft.fft(y-y.mean())
        
        # calculate the Amplitude and Phase
        A, P = np.abs(sfft), np.angle(sfft)
        fr = fftshift(fftfreq(xr.shape[0]))
        g_top = max(fftshift(A))
        
        # apply the FFT high pass filter
        yFFT_filtered = fftshift(A) \
                        * (1 - self.gfunc(fr, g_top, fmin)/g_top)

        # calculate the filtered result in real space
        return real(ifft(ifftshift(yFFT_filtered) * exp(1j*P)))

    
    def fom(self, raw, bounds=[400, 1700], ui=None):
        """Computes the FOM of the spectrum
           Removing contribution from caddy raman signal
        """

        # calculate the baseline:
        spec_bl = rampy.modPoly(np.array(raw))

        # define the x-range
        xr = WAVENUMBER_RANGE
        
        # get the spectrum
        spec = spec_bl.y-spec_bl.bl
        
        # calculate the noise
        noise = rampy.calc_sig_noise(xr, spec, bounds)

        # calculate bounds indecies
        bl, br = np.argmin(abs(xr-bounds[0])), np.argmin(abs(xr-bounds[-1]))

        if ui is not None:
            self.update_master_plot(ui, spec, xr)

        # return the SNR for the largest deviation (could be negative)
        return spec[bl+np.argmax(abs(spec[bl:br]))]/noise
    
        
    def profile(self, axis, x_in, span=3, delta=0.2, ui=None):
        """Find the best position along the 'axis'
           It will search around x_in in a range of 'span' mm 
           along direction 'axis' with a step 'delta' mm
           Then will FIT a Gaussian to the FOM profile
           
           Returns:
               x_out, FOM profile, positions array
           Example:
               hpl.profile(X, 10, 5, 0.1)
        """
        if (not self.enabled[0]) or (not self.enabled[1]):
            return
        
        # set positions to probe
        steps = round(span/delta)
        pos = linspace(x_in-span/2, x_in+span/2, steps)
        int_vals = []
        
        # iterate through positions
        for i in range(steps):
            self.stage.move(axis, pos[i], ABS, False)

            # get the spectrum and baseline-correct it:
            spec_raw = self.snri.acquire()
            
            # integrate the signal
            int_vals.append(self.fom(spec_raw, ui=ui))


        # initial guess
        p0 = [min(int_vals),
              max(int_vals)-min(int_vals),
              pos[np.argmax(int_vals)],
              3 if axis == 'Z' else 0.5]

        self.update_master_plot(ui, int_vals, pos)

        idx_max = np.argmax(int_vals)

        if idx_max > 0 and idx_max < len(pos)-1:
            if int_vals[idx_max-1] > int_vals[idx_max+1]:
                return 0.5*(pos[idx_max]+pos[idx_max-1])
            else:
                return 0.5*(pos[idx_max]+pos[idx_max+1])

        elif idx_max > 0:
            return 0.5*(pos[idx_max]+pos[idx_max-1])

        return 0.5*(pos[idx_max]+pos[idx_max+1])


    def acquire(self, flush=False):
        """Collect spectrum from SNRI

        Returns:
            An NumPy array of intensities read from the instrument.

        """

        if not self.enabled[0]:
            return
        
        intensities = self.snri.acquire(flush=flush)
        return intensities

    
    def autofocus(self, ui_ref, pbar=None, pbar_step_i=0):
        """autofocus around current z position
           Decide how many millimeters to scan
        """
        if (not self.enabled[0]) or (not self.enabled[1]):
            return

        for cur_z,span_z,delta_z in [[self.stage.get_position('Z'), 3, 0.7],
                                     [None, 1, 0.2]]:
            if cur_z is None:
                if pbar is not None:
                    pbar.update(pbar_step_i+1)
                    pbar.set_status('Running fine scan')
                
                # this is a second pass
                cur_z = pos[np.argmax(ivals)]

            elif pbar is not None:
                pbar.set_status('Running rough scan')
                
            #  ensure the z is within bounds
            cur_z = self.stage.limit_z(cur_z, span_z)

            # find rough position of focus
            z_max = self.profile(Z, cur_z, span_z, delta_z, ui_ref)
              
        if pbar is not None:
            pbar.update(pbar_step_i+2)
            pbar.set_status('Moving to best focus')
        
        # go to best position
        self.stage.move(Z, z_max, ABS, True)
        self.update_master_plot(ui_ref, self.snri.acquire(), WAVENUMBER_RANGE)
        
        return z_max
    

    def close(self, pbar=None):
        """Shut down instruments
           go to desider close position
           Close communication
        """

        if pbar is not None:
            pbar.set_status('Closing connection to the stage ...')

        if self.enabled[1]:
            self.stage.close(pbar=pbar)
            
        if pbar is not None:
            pbar.update(2)
            pbar.set_status('Closing connection to the spectrometer ...')

        if self.enabled[0]:
            self.snri.close()
            

    def update_master_plot(self, ui, y, x=None):
        """This function will plot something on the master GUI.

        Args:
            ui: the pointer ot the master ui
        
            y: the y-axis values

            x: the x-axis values, if not provided, just the indecies
                of the y-values. Defaults to None.

        """

        if ui is not None:
            # show the graph in the UI
            ui.mpl_ax.cla()
            update_plot(ui)
            ui_plot_graph(ui, x if x is not None else WAVENUMBER_RANGE, y)
            ui.figCanvas.show()
