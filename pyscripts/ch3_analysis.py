''' Script holding the main functions necessary for analyses of mooring data
collected during PISTON '''

import xarray as xr; import numpy as np; import pandas as pd;
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
import sys, os 
sys.path.append('/mnt/sda1/PhysOc/modview/modview');
import loader, mapper
import cartopy as cart

# Some of the main functions that will need to be included are:
# - - - load argo climatology
# - - - load mooring data using loader
# - - - load wind/rinfall climatology
# - - - load ERA5 wind/rainfall/heatfluxes for 2018-2019
# - - - interpolate mooring ADCP data (nc may need reprocessing)

''' SAVE DIRECTORIES AND FILENAMES FOR MAIN DATASOURCES'''
data_dirs = dict(); 
data_dirs['wind'] =  '/mnt/sda1/SciData/Wind/'
data_dirs['argo'] = '/mnt/sda1/SciData/ARGO/'
data_dirs['precip'] = '/mnt/sda1/SciData/CMAP/'
data_dirs['era5'] = '/mnt/sda1/SciData/ERA5/';
data_dirs['SIO3'] = '/mnt/sda1/SciData/PISTON/Moorings/SIO3/';
data_dirs['SIO1'] = '/mnt/sda1/SciData/PISTON/Moorings/SIO1/';

# Set information necessary to access mooring data
def mooring_paths(num):
    ''' Dictionaries storing arguments that will be passed to instantiate
    loader.assemble '''
    moor_paths = dict(); 
    if num == 1:
        adcp_sn = [3160,4021,11181,14255];
        prel = 'SIO1';
    elif num == 3:
        adcp_sn = [13596,14435,22525];
        prel = 'SIO3'
    adcp_filenames = [prel+'_'+str(kk)+'.nc' for kk in adcp_sn];
    moor_paths['ADCP'] = [data_dirs[prel]+'ADCP/'+ jj for jj in adcp_filenames]
    moor_paths['chipod'] = data_dirs[prel] + prel + '_chipod.mat';
    moor_paths['chiflags'] = data_dirs[prel] + prel + '_chipod_finflag.csv'
    moor_paths['temp'] = data_dirs[prel] + prel + '_gridT_forxr.mat'
    
    moor_id = {'name':prel, 'project':'PISTON', 'type':'Mooring',
            'limits':{'t0':'2018-09-20','t1':'2019-11-01'} }
    return moor_paths, moor_id

# Climatology datasets
clim_paths = dict(); 
clim_paths['wind_clim'] = data_dirs['wind'] \
        + 'Monthly_global_winds.nc';
clim_paths['argo'] = data_dirs['argo'] \
        + 'RG_ArgoClim_33pfit_2019_annual.nc';
clim_paths['precip'] = data_dirs['precip'] \
        + 'precip_monthly_means.nc';
clim_paths['era5_monthly'] = data_dirs['era5'] \
        + 'indopacific_monthly_20yrs.nc';
clim_paths['era5_pentad'] = data_dirs['era5'] \
        + 'indopacific_pentad_20yrs.nc';

# Datasets corresponding to our period of observations
obs_paths = dict(); 
obs_paths['era5'] = data_dirs['era5'] + 'heat_wind_rain_PISTON_6hr.nc';

# ------ set geographical limits ---------
bg_lims = {'t0':3,'t1':10,'lon0':60,'lon1':180,'lat0':-20,'lat1':20,
        'p0_ocean':0,'p1_ocean':200,
        'p0_atmosphere':1100,'p1_atmosphere':800};

clim_id = {'name':'Background climate','project':'PISTON','type':'Hybrid',
        'limits':{'lon0':80,'lon1':160,'lat0':-25,'lat1':25}};

# Prepare data from reanalyses and remote sensing datasets

def get_argo():
    ARGO = xr.open_dataset(clim_paths['argo'], decode_times=False); 
    ARGO = mapper.standardize_coords(ARGO,'ARGO_TEMPERATURE_ANNUAL_ANOMALY'); 
    return ARGO

def get_era5(use='monthly', cut=True):
    if use == 'obs':
        ERA5 = xr.open_dataset( obs_paths['era5'] );
    else:
        key = 'era5_' + use;
        ERA5 = xr.open_dataset(clim_paths[key]); 
    # Reorder and slice data
    ERA5 = ERA5.reindex( {\
            'latitude':list(reversed(ERA5.latitude))})
    if cut:
        ERA5 = ERA5.sel(longitude= \
                slice(bg_lims['lon0'],bg_lims['lon1']));
        ERA5 = ERA5.sel(latitude= \
                slice(bg_lims['lat0'],bg_lims['lat1'])); 
    ERA5['u_comp'] = ERA5['u10'] + 1j*ERA5['v10']; 
    #ERA5['Q_rad'] = ERA5['msnswrf'] + ERA5['msnlwrf']; 
    #ERA5['Q_turb'] = ERA5['mslhf'] + ERA5['msshf']; 
    return ERA5

def get_precip(use=None):
    if use is None:
        # use NOAA CPC monthly climatology
        PRECIP = xr.open_dataset(clim_paths['precip']); 
        PRECIP = mapper.standardize_coords(PRECIP, 'precip');
        PRECIP = PRECIP.reindex({'latitude':list(reversed(PRECIP.latitude))});
        PRECIP = PRECIP['precip']; 
    elif isinstance(use,xr.Dataset):
        # use is the era5 datast
        PRECIP = use['mtpr']; # mean total precipitation rate [kg m^{-2} s^{-1}]
        PRECIP = PRECIP*24*3600; # convert to [mm day^{-1}]
    PRECIP = PRECIP.sel(longitude=slice(bg_lims['lon0'],bg_lims['lon1']));
    PRECIP = PRECIP.sel(latitude=slice(bg_lims['lat0'],bg_lims['lat1']))
    return PRECIP


# Create assemble object for mooring
def get_mooring(num=3, get_chipod=False):
    
    paths, ID = mooring_paths(num); # filepaths and settings
    moor = loader.assemble( ID, paths); # load methods 

    moor.store_nc(moor.paths['ADCP'], ['u','v','pressure']); # import adcp data
    # ----- Intermediate step necessary to treat discontinuities ---------
    moor.interp_grid( variables=['u','v'], sorter='z'); # unite ADCP files
    
    tgrid = loader.matfile(moor.paths['temp']); # temp grid from matfile
    #struct_name = os.path.basename( moor.paths['temp'] );
    #struct_name = struct_name[0:4]; # put grid into xr format
    tgrid.data['temp'] = tgrid.var_to_xr(['temp'],['depth','time'], 
                  in_struct=ID['name'])
    # standard time format
    tgrid.data['temp']['time'] = pd.to_datetime(tgrid.data['temp'].time.values); 
    moor.grids['temp'] = tgrid.data['temp'];
    # get chipod if requested
    if get_chipod:
        moor = add_chipod(moor);
    return moor


def add_chipod(assembler):
    # load data and access variables in structure
    chi_file = loader.matfile(assembler.paths['chipod']);
    vars2get = ['time','eps','Kt','Jq','SN','znom']
    struct_name = os.path.basename(assembler.paths['chiflags']);
    struct_name = struct_name[0:11];# name of structure to inspect
    chi_data = chi_file.read_struct(struct_name, vars2get); # read data
    # get flags for estimates on each instrument
    chi_flags = pd.read_csv(assembler.paths['chiflags']); 
    assembler.vars['chipod'] = dict(); 
    for item in ['SN','znom']:
        assembler.vars['chipod'][item] = chi_data[item]; # no flagging needed
    for item in ['eps','Kt','Jq']:
        assembler.vars['chipod'][item] = loader.clean_chivar( \
                chi_data, chi_flags, item, as_xr=True); # time extracted w/ each var
    # Before gridding, account for knockdown
    assembler.vars['knockdown'] = assembler.grids['pressure'][0] \
                - np.nanmin(assembler.grids['pressure'][0]); # p data from 1 adcp
    chi_time = pd.to_datetime(assembler.vars['chipod']['eps'].time.values); 
    z_offset = assembler.vars['knockdown'].interp( time=chi_time, 
                method='nearest'); 
    znom = assembler.vars['chipod']['znom'][0]; # nominal depths
    zchi = np.expand_dims(znom,axis=1) + z_offset.values; # add knockdown 
    # replace dict for xr.Dataset
    chi_data = xr.Dataset(data_vars={\
            'eps':(['time','pressure'], assembler.vars['chipod']['eps']),
            'Kt':(['time','pressure'], assembler.vars['chipod']['Kt']),
            'Jq':(['time','pressure'], assembler.vars['chipod']['Jq']),
            'znom':(['pressure'],assembler.vars['chipod']['znom'][0]),
            'SN':(['pressure'],assembler.vars['chipod']['SN'][0])},
            coords={'time':chi_time,'pressure':assembler.vars['chipod']['znom'][0],
                'depth':(['time','pressure'], zchi.transpose() ) } );
    # ------ sort instruments by increasing depth
    assembler.vars['chipod'] = chi_data.reindex( {'pressure':np.sort(znom)} )
    return assembler

def plot_chipod_grid( ax, assembler, variable, limits, settings):
    
    data = assembler.vars['chipod'][variable];
    data = data.sel(time=slice(limits['t0'],limits['t1']))
    data = data.resample(time='12H').mean()
    znom = assembler.vars['chipod']['znom'];
    img = data.values.transpose();
    if settings['log']:
        img = np.log10(img); # for plotting eps, kappa

    contour = ax.contourf( data.time, znom, img, cmap=settings['cmap'], 
                            levels=settings['levels'], extend='both')
    Ti = np.repeat( data.time[10].values, data.values.shape[1] );
    ax.scatter( Ti, znom, s=40, color='red', marker='*')
    return contour 
    

def at_mooring(gridded,lon=134.7,lat=15.7):
    gridded = gridded.sel(longitude=lon,method='nearest');
    gridded = gridded.sel(latitude=lat,method='nearest');
    return gridded

