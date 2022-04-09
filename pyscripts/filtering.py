from scipy.signal import butter, filtfilt, freqz
import numpy as np

def xpass(var, frads, dt_obs, filt_type, order = 4):
    # Create digital filter for evenly spaced data
    Wn = makefreq( frads, dt_obs ) 
    b,a = butter( order, Wn, btype = filt_type, analog = False)
    #w,h = freqz(b,a)
    #plt.plot( w, np.log10(abs(h)))
    # Assume that longest dimension is time dimension (generally true)
    inputshape = var.shape;
    return filtfilt(b,a,var, axis = inputshape.index(max(inputshape)));

def makefreq(frads, dt_obs):
    # Create frequency input for xpass based on:
    # frads  - desired frequency in radians / second
    # dt_obs - sampling rate in datapoints / second
    
    # Scipy.signal tools take frequencies in units (half cycles / sample). 
    normT = 1/dt_obs; # seconds per sample
    Wn = frads / np.pi ; # convert to half-cycles / second. 
    Wn = Wn * normT ; # half-cycles per sample 
    return Wn     
    
def nearf(lat,width = [0.8,1.2]):
    # Find band-limit frequencies for the near-inertial band
    f = 4*np.pi*np.sin(lat/180*np.pi)/24/3600; # inertial freq in rads/second
    return f*np.array(width)
