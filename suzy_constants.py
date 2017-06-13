"""This file presets the constants that can be accessed from other
parts of this package.

"""

import sys
import os
import numpy as np

# setup the paths
RAMPY_PAR_DIR = os.path.abspath(os.path.join(sys.path[0], os.pardir))
ICON_PATH = os.path.abspath(os.path.join(sys.path[0], 'logo.ico'))

# make sure the rampy module is loaded
sys.path.append(RAMPY_PAR_DIR)

# system constants
VERSION = 2.0

DEBUG_SERIAL = False

# enable/disable feature sets
FEATURE_SET2 = 'normal'
FEATURE_SET3 = 'disabled'
FEATURE_SET4 = 'normal'

# controls constants
ORIGIN_POS = [73.6,63.6,1.0]
LOAD_POS = [5,5,1]

WAVENUMBER_RANGE = np.linspace(200, 2000, 1801)

AVE_OPTIONS = ['1', '2', '3', '5', '7', '10', '15', '20', '30', '50']
CELL_NAMES = ['A', 'B', 'C']
CELL_IDS = list(map(str, range(1,4)))
CADDY_SIZE = np.array([10,10])


CELL_MAP = dict(A1='C1', A2='B1', A3='A1',
                B1='C2', B2='B2', B3='A2',
                C1='C3', C2='B3', C3='A3')

CELL_RMAP = dict(A1='A3', A2='B3', A3='C3',
                 B1='A2', B2='B2', B3='C2',
                 C1='A1', C2='B1', C3='C1')


# controls optimization constrants
SETT_OPTIMIZE = {'AVE': 1, 'REF': 'OFF', 'BAS': 'OFF', 'RAS': 'ON',
                 'PWR': 'ON', 'SHU': 'OPEN'}
OPT_THRESHOLD = 20
OPT_SPAN = 3
OPT_DELTA = 0.4
OPT_SPOT = 1.3
OPT_ERR = 0.1
STAGE_CMD_DELAY = 0.2

# plotting constants
RANGE_NONE = 'None'
RANGE_FIXED = 'Fixed'
RANGE_RANGE = 'Range'
NORMAL_TYPE = ['None', 'Fixed', 'First', 'Last', 'Average', 'Minimum',
               'Maximum']
LINE_STYLES = ['None', '-', '--', '-.', ':']
MARKER_TYPES = ['None', 'o', 's', 'D', '^']

# define some global constants
DT_SHOW_ALL = 0
DT_SHOW_CORRECTED = 1
DT_SHOW_RAW = 2

# define the smoothing algorithm
SMOOTH_NONE = 0
SMOOTH_GAUSS = 1
SMOOTH_SAVGOL = 2

# define the baseline algorithm constants
BL_ALGOS = ['ALS', 'HPL', 'MP']
ALS, HPL, MP = 'ALS', 'HPL', 'MP'

ALS_KEYS = ['lam', 'p', 'niter']
ALS_KNAMES = ['Lambda', 'p', 'num. of iterations']

HPL_KEYS = ['fmin', 'errTol', 'rawKernel', 'sigKernel', 'sigGradKernel']
HPL_KNAMES = ['Roll off f.', 'Error tol.', 'Raw Kernel (1/cm)',
              'Sig Kernel (1/cm)', 'Grad. Kernel (1/cm)']

MP_KEYS = ['pwr', 'FOM', 'maxiter']
MP_KNAMES = ['Polynomial pwr.', 'FOM', 'Num. iterations']

BL_RANGE = {'fmin': [0.001, 0.5],
            'errTol': [0, 10000],
            'rawKernel': [1, 100],
            'sigKernel': [1, 100],
            'sigGradKernel': [1, 100],
            'lam': [1e2, 1e9],
            'p': [0.001, 0.1],
            'niter': [1, 1000],
            'pwr': [3, 6],
            'FOM': [0.00001, 1],
            'maxiter': [1, 1000]}


# default settings
PEAK_DET_SETT = {'fmin': 0.007, 'rawKernel': 6, 'minSNR': 1.5,
                 'noiseFWHM': 30, 'noise_slen': 300, 'range': [400, 1600]}
PEAK_DET_SETT_RANGES = [[0.001, 2], [2, 60], [1, 1000],
                        [2, 60], [50, 1000], [200, 1800], [200, 1800]]


BASELINE_SETT = {'fmin': 0.007, 'errTol': 0,
                 'rawKernel': 30, 'sigKernel': 10,
                 'sigGradKernel': 12,
                 'lam': 1e5, 'p': 0.01, 'niter': 10,
                 'pwr': 4, 'FOM': 0.001, 'maxiter': 100,
                 'algo': ALS}


# big strings
als_instructions = """
The assymetric least squares algorithm implementation is described
in "Asymmetric Least Squares Smoothing" by P. Eilers and H. Boelens 
in 2005. There are three parameters to configure:

* Lambda - takes values between 1e2 and 1e9. Default is 1e5.

* p - takes values between 0.001 and 0.1. Default is 0.01.

* # interations - default is 10."""
        
hpl_instructions = """    
The settings for the HP Labs baseline are shown below (with valid range 
indicated in the paranthesis):

* Roll-off frequency: the cutoff frequency for the high-pass 
      FFT filter (0.001, 0.5)

* Error threshold: the threshold where a negative amplitude 
      counts as noise (0, 10000)

* Raw smoothing FWHM: the first-pass smoothing kernel FWHM 
      in px (1, 100)

* Sig. smoothing FWHM: the second-pass smoothing kernel FWHM 
      in px (1, 100)

* Grad. smoothing FWHM: the smoothing kernel for the gradient 
      calc. (1, 100)"""

mp_instructions = """
The ModPoly baseline method attempts to find an polynomial 
estimate of the baseline. The algorithm fits a polynomial to 
the raw data and selectively subtracts it in each calculation
cycle. The residue is then used as the new input signal and
fitted again in the next cycle.

"Automated Autofluorescence Background Subtraction Algorithm
for Biomedical Raman Spectroscopy" by JIANHUA ZHAO, HARVEY LUI, 
DAVID I. MCLEAN, and HAISHAN ZENG (2007).

* Polynomial power: defaults to 4.

* Figure of Merit (FOM): defaults to 0.001.

* # iterations: defaults to 100."""

about_string = """Hi there! I'm Suzy. Are you interested in plotting SERS data?
Maybe I can help you with checking your data using a baseline
correction and cross-set data comparison?

Jokes aside, Suzy provides a platform to quickly browse through the
SERS data and run common data processing tasks.  The goal is to
provide a place for testing algorithms developed (or adopted) at HP
for data analysis. For advanced data massaging, please use a proper
software/development environment.

NOTE: THIS SOFTWARE HAS BUGS AND FUNDAMENTAL MISTAKES.
      YOU CAN HELP FIX THEM BY SENDING AND EMAIL TO 
      Alex Polyakov (polyakov@hp.com)

Enjoy!
"""

peak_det_inst = """The peak detection works by computing a derivative 
on the smoothed input signal. For each peak found, the SNR is calulated
and only the peaks exceeding the min. SNR threshold are kept. To 
estimate the noise level, a blurred (smoothed) input signal is 
subtracted from the raw input signal. This yields a measure of the high
frequency noise present in the raw input. By running a window of a set
number of points (noise # points) and reperatedly calculating the STD,
we find the minimum STD as the estimate of the noise.

The algorithm has the following settings (with valid range indicated 
in the parenthesis):

    Cutoff freq [1/(1/cm)]: the cut-off frequency for the initial input 
        signal high pass Fourier filter (0.001, 2).

    Signal blur (1/cm): the kernel used for smoothing the input in the 
        gradient calculation (2, 60)

    Minimum SNR: reject any peaks with low SNR (1, 1000)

    Noise blur (1/cm): the kernel used for smoothing the noise in the 
        noise STD calculation (2, 60)

    Noise # points: number of points to consider for estimating 
                    noise (50, 1000)

    Minimum wnum (1/cm): the minimum peak wave number (200, 1800)

    Minimum wnum (1/cm): the maximum peak wave number (200, 1800)

"""
    
