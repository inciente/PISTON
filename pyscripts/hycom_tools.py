import numpy as np; import sys; 
import xarray as xr; import datetime;
import matplotlib.pyplot as plt;
import pandas as pd; 
sys.path.append('/mnt/sda1/PhysOc/modview/modview/')
import loader; import viztools; import timetools; import mapper
sys.path.append('/mnt/sda1/PhysOc/DataFinder/')
import filtering

def get_prep( xr_obj ):
    # First steps that a xr.DataArray() must go through before
    # it is compatible with all the functions below.
    xr_obj = mapper.standardize_coords( xr_obj ); 
    return xr_obj 

def get_dt(xr_obj):
    # Return mean dt
    time = pd.to_datetime( xr_obj.time.values )
    dt = ( time[1] - time[0] ).total_seconds()
    return dt

def latlon_deriv(xr_obj, func=None):
    # Return lat lon derivatives of xr object
    if func is not None:
        xr_obj = func(xr_obj); # Apply function
    var_x = xr_obj.differentiate(coord='longitude')/ \
            110e3/np.cos(xr_obj.latitude.mean()/180*np.pi)
    var_y = xr_obj.differentiate(coord='latitude')/110e3
    return var_x, var_y

def get_vort_div(flow):
    f = 4*np.pi*np.sin(flow.latitude/180*np.pi)/24/3600;
    u_x, u_y = latlon_deriv(np.real(flow)) # spatial derivatives
    v_x, v_y = latlon_deriv(np.imag(flow))
    vort = (v_x - u_y)/f; # vorticity
    div = (u_x + v_y)/f; # divergence
    return vort, div

def near_inertial_filter( xr_obj ):
    # Take in xr_object and run it through filtering.xpass
    # does it for a single value of latitude
    latitude = xr_obj.latitude.values; 
    if latitude.shape != () :
        print('I dont know how to use multiple latitudes yet. Using mean')
        latitude = np.mean( latitude )
    dt = get_dt( xr_obj )     
    my_filter = filtering.xpass; # in case i want to use a different one
    inertial_freq = 4*np.pi*np.sin( latitude /180*np.pi)/24/3600
    freq_band = np.array( [0.7, 1.3] ) * inertial_freq
    # ===== now run the filter
    fltrd = my_filter( xr_obj, freq_band, 1/dt, 'bandpass' );
    fltrd = xr.DataArray( data = fltrd, coords = xr_obj.coords )
    return fltrd

def sel_dat( xr_obj, x=None, y=None, t=None, p=None ): 
    # run sel at the given coordinates
    all_dims = xr_obj.dims;
    if x is not None:
        xr_obj = xr_obj.sel( longitude=x, method='nearest' ); 
    if y is not None:
        xr_obj = xr_obj.sel( latitude=y, method='nearest' ); 
    if t is not None:
        xr_obj = xr_obj.sel( time=t, method='nearest');
    if p is not None:
        xr_obj = xr_obj.sel( pressure=p, method='nearest'); 
    return xr_obj            

def filter_by_lat( xr_obj, filt=near_inertial_filter):
    # Run a series of filters for all latitude values in 
    pass
