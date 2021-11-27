import xarray as xr

def tot_flux(era_data):
    total = era_data['mslhf'] + era_data['msshf'] \
               + era_data['msnswrf'] + era_data['msnlwrf']
    return total
    
def cut_resamp(data, limits, resample):
    data = data.sel(time=slice(limits['t0'],limits['t1']));
    if resample != 'none': 
        data = data.resample(time=resample).mean();
    return data

def flux_plot(era_data, ax, limits, resample='24h'):
    era_data = cut_resamp(era_data, limits, resample); 
    ax.plot(era_data.time, era_data['mslhf'],label='latent');
    ax.plot(era_data.time, era_data['msshf'],label='sensible'); 
    ax.plot(era_data.time, era_data['msnswrf'], label='shortwave');
    ax.plot(era_data.time, era_data['msnlwrf'], label='longwave'); 
    ax.plot(era_data.time, tot_flux(era_data), linewidth=2.5, label='total');
    ax.legend(loc='lower center', bbox_to_anchor=(0.5,1.05), fancybox=True,
              ncol=5);
    ax.set_ylim((-325,325))

def sst_plot(era_data, ax, limits, resample='24h',**kwargs):
    era_data = cut_resamp(era_data, limits, resample); 
    ax.plot( era_data.time, era_data['sst']-273.15, **kwargs); 
    #ax.set_ylim((27,30))

def tau_plot(era_data, ax, limits, resample='none', **kwargs):
    era_data = cut_resamp(era_data, limits, resample); 
    tau = 2e-3*1.22*(era_data['u10']**2 + era_data['v10']**2)
    ax.plot(era_data.time, tau);


