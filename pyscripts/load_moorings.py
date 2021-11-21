# Pull data from database 
sio1_limits = {'t0':'2018 Sep 15','t1':'2019 Oct 31'};
sio1_data_dir = '/media/mydrive/PISTON/Moorings/SIO1/'; 
#sio1_data_dir = '/mnt/sda1/SciData/PISTON/Moorings/SIO1/';
sio1_adcp_paths = ['SIO1_4021.nc','SIO1_11181.nc','SIO1_3160.nc'];
sio1_paths = {'chipod':sio1_data_dir+'SIO1_chipod.mat','chiflags':sio1_data_dir+'SIO1_chipod_finflag.csv',
             'ADCP':[sio1_data_dir + 'ADCP/' + kk for kk in sio1_adcp_paths],
             'Temperature':sio1_data_dir + 'SIO1_gridT_lores.mat'};
sio1_id = {'name':'SIO1','project':'PISTON','type':'Mooring','limits':sio1_limits};

sio3_limits = {'t0':'2018 Sep 15','t1':'2019 Oct 31'};
sio3_data_dir = '/media/mydrive/PISTON/Moorings/SIO3/'; 
#sio1_data_dir = '/mnt/sda1/SciData/PISTON/Moorings/SIO1/';
sio3_adcp_paths = ['SIO3_13596.nc','SIO3_14435.nc','SIO3_22525.nc'];
sio3_paths = {'chipod':sio3_data_dir+'SIO3_chipod.mat','chiflags':sio3_data_dir+'SIO3_chipod_finflag.csv',
             'ADCP':[sio3_data_dir + 'ADCP/' + kk for kk in sio3_adcp_paths],
             'Temperature':sio3_data_dir + 'SIO3_gridT_forxr.mat'};
sio3_id = {'name':'SIO3','project':'PISTON','type':'Mooring','limits':sio1_limits};

def make_tempgrid(obj,namestr):
    temp_file = loader.matfile(obj.paths['Temperature']); 
    temp_file = temp_file.read_struct(namestr,['depth','time','temp']); 
    time = pd.to_datetime(temp_file['time']); 
    depth = np.squeeze(temp_file['depth']); 
    obj.vars['Temperature'] = xr.DataArray( data=temp_file['temp'], 
                     dims=['depth','time'], 
                     coords={'depth':depth,'time':time});
    return obj


# ------ -------- prepare chipod data
# now use the maftile class to load chipod data
def make_chipod(obj,namestr):
    chi_vars2get = ['time','eps','Kt','Jq','SN','znom'];
    chi_file = loader.matfile(obj.paths['chipod']); # get data
    chi_data = chi_file.read_struct(namestr+'_chipod',chi_vars2get); 
    
    chi_flags = pd.read_csv(obj.paths['chiflags']); # get flags
    
    obj.vars['chipod'] = dict(); # to store flagged data 
    for item in ['znom','SN']:
        obj.vars['chipod'][item] = chi_data[item] # save each chipods info
    for item in ['eps','Kt','Jq']: # vars to be interpolated
        obj.vars['chipod'][item] = loader.clean_chivar( \
                        chi_data, chi_flags,item, as_xr=True);
    del chi_data, chi_file
    # All data is now stored indexed by serial number and time. 
    # Make xr objects that facilitate plotting and accounting for knockdown
    obj.vars['knockdown'] = obj.grids['pressure'][0] \
                  - np.nanmin(obj.grids['pressure'][0]);
    # interpolate knockdown to chipod time
    raw_chi_time = pd.to_datetime(obj.vars['chipod']['eps'].time.values) 
    chi_offset = obj.vars['knockdown'].interp( time=raw_chi_time, 
                                           method='nearest')
    # knockdown + nominal depth
    zchi = np.expand_dims(obj.vars['chipod']['znom'][0],axis=1) \
            + chi_offset.values;
    obj.vars['chipod']['z'] = xr.DataArray( data=zchi.transpose(), \
                    dims=['time','SN'])
    # Replace separate grids for single dataset
    obj.vars['chipod'] = xr.Dataset(data_vars={\
         'eps':(['time','z'],obj.vars['chipod']['eps']),
        'Kt':(['time','z'],obj.vars['chipod']['Kt']),
        'Jq':(['time','z'],obj.vars['chipod']['Jq'])},
        coords={'time':raw_chi_time, 'z':obj.vars['chipod']['znom'][0],
               'depth':(['time','z'],zchi.transpose())}); 
    return obj


def make_mooring(strname='SIO1'):
    if strname=='SIO1':
        moor = loader.assemble(sio1_id,sio1_paths);
        limits = sio1_limits;
    elif strname=='SIO3':
        moor = loader.assemble(sio3_id,sio3_paths);
        limits = sio3_limits;
    # Dump data from ADCPs
    moor.store_nc(moor.paths['ADCP'],['u','v','pressure'],limits=limits);
    moor.interp_grid(['u','v'],'z'); 
    moor = make_tempgrid(moor,strname); 
    moor = make_chipod(moor,strname);
    return moor
