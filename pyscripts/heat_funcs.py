import xarray as xr

def tot_flux(era_data):
    total = era_data['mslhf'] + era_data['msshf'] \
               + era_data['msnswrf'] + era_data['msnlwrf']
    return total

def flux_plot(era_data, ax, limits, resample='24h'):
    era_data = era_data.sel(time=slice(limits['t0'],limits['t1']));
    era_data = era_data.resample(time=resample).mean();
    ax.plot(era_data.time, era_data['mslhf'],label='latent');
    ax.plot(era_data.time, era_data['msshf'],label='sensible'); 
    ax.plot(era_data.time, era_data['msnswrf'], label='shortwave');
    ax.plot(era_data.time, era_data['msnlwrf'], label='longwave'); 
    ax.plot(era_data.time, tot_flux(era_data), linewidth=3, label='total');
    ax.legend();

