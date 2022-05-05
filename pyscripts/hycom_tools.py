import numpy as np; import sys, socket; 
import xarray as xr; import datetime;
import matplotlib.pyplot as plt;
import pandas as pd; 

if socket.gethostname() == 'kuxan-suum':
    sys.path.append('/media/mydrive/PhysOc/modview/modview/');
    sys.path.append('/media/mydrive/PhysOc/DataFinder/');
elif socket.gethostname() == 'Chijpiri':
    sys.path.append('/mnt/sda1/PhysOc/modview/modview/')
    sys.path.append('/mnt/sda1/PhysOc/DataFinder/')
# import the home made software
import loader; import viztools; import timetools; import mapper
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
 
def scale_by_f( xr_obj, power=1 ):
    const = 4*np.pi/24/3600;
    f = const * np.sin( xr_obj.latitude / 180 * np.pi ); 
    scaled = xr_obj / f**power; 
    return scaled

class flow_field:
    def __init__(self, u, v, w=None):
        # takes in xr.DataArrays
        self.u = get_prep( u + 1j*v ); # standardize coords
        self.w = w; 
    def vort_div(self, scale=True):
        ux, uy = latlon_deriv( self.u ); 
        vort = np.imag( ux ) - np.real( uy ); 
        div = np.real( ux ) + np.imag( uy ); 
        if scale:
            vort = scale_by_f( vort, power=1 ); 
            div = scale_by_f ( div, power=1 ); 
        return vort, div
    
    def nonlinear_uv( self, scale=False ):
        # Compute nonlinear terms in the equations of motion
        ux, uy = latlon_deriv( self.u ); 
        accel_u = np.real(self.u) * np.real(ux) \
                + np.imag(self.u) * np.real(uy) ; # u*u_x + v*u_y
        accel_v = np.real(self.u) * np.imag(ux) \
                + np.imag(self.u) * np.imag(uy); # u*v_x + v*v_y
        if scale: 
            accel_u = scale_by_f( accel_u, power=1 );
            accel_v = scale_by_f( accel_v, power=1 ); 
        return accel_u, accel_v
    
    def nonlinear_vort_div( self, scale=False ):
        # Nonlinear contributions to dzeta/dt, dGamma/dt
        # following Nevir & Sommer (2009) in JAS
        vort, div = self.vort_div( scale=False ); 
        vort_x, vort_y = latlon_deriv( vort ); 
        vort_transport = np.real( self.u ) * vort_x \
                + 1j * np.imag( self.u ) * vort_y;
        trans_x, trans_y = latlon_deriv( vort_transport ); 
        # divergence of vorticity transport
        accel_vort = np.real( trans_x ) + np.imag( trans_y );
        # curl of vorticity transport
        accel_div = np.imag( trans_x ) - np.real( trans_y ); 
        if scale: 
            accel_vort = scale_by_f( accel_vort, power=2 ); 
            accel_div = scale_by_f( accel_div, power=2 ); 
        return accel_vort, accel_div
    
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
        xr_obj = xr_obj.sel( time=t);
    if p is not None:
        xr_obj = xr_obj.sel( pressure=p, method='nearest'); 
    return xr_obj            

def filter_by_lat( xr_obj, filt=near_inertial_filter):
    # Run a series of filters for all latitude values in 
    pass
