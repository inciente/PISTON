import xarray as xr; import pandas as pd; import numpy as np;
import datetime; import scipy.io as sio

''' Assortment of functions that will help produce overviews of chipod data availability, trends, averages, and others. Initial processing is not handled here. '''

def get_SN(moor, sn):
    chipods = moor.vars['chipod'];
    sn_in_moor = chipods['SN'].values;
    nominal_depth = chipods['znom'].values; 
    # Find index of instrument requested (sn) in the pressure dimension
    which_chipod = np.array([ instrument == sn for instrument in sn_in_moor]); 
    print('getting ' + str(sn) )
    if sum(which_chipod) > 1 or sum(which_chipod)==0:
        print('Multiple or none chipods with serial number ' \
                + str(sn) )
        return
    else:
        nominal_depth = nominal_depth[which_chipod];
        # get data for that depth (corresponds to sn)
        chi_dat = moor.vars['chipod'].sel(pressure=nominal_depth); 
        return chi_dat

def monthly_availability(moor, sn):
    # Return percentage of data from instrument sn that is valid for different seasons
    chi_dat = get_SN(moor, sn); 
    time = pd.to_datetime(chi_dat.time.values); 
    good_ones = availability(chi_dat); # times of valid estimates
    perc_good = np.zeros(12); 
    for kk in range(1,13):
        check_month = time.month == kk # dates within month kk
        good_here = time[good_ones].month == kk; # those with valid estimates
        # percentage of estimates for a month that are valid
        perc_good[kk-1] = np.sum(good_here) / np.sum(check_month); 
    return perc_good

def availability_summary(moor):
    sn_here = moor.vars['chipod'] 
    summary = pd.DataFrame(columns=['Jan','Feb','Mar','Apr','May','Jun',
        'Jul','Aug','Sep','Oct','Nov','Dec']); 
    for sn in sn_here['SN'].values:
        perc_monthly = monthly_availability(moor, sn);
        index = len(summary.index)
        summary.loc[index] = [str(sn), *perc_monthly]; # add instrument to table
    return summary

def availability(grouped_chi):
    # Take in grouped dataset of chipod output and Return the fraction 
    # estimates that is valid for each group
    good_fraction = grouped_chi.count() # number of items that are not null
    good_fraction = good_fraction/good_fraction['check_time']
    return good_fraction


def group_summary(moor,group_by='season'):
    chipods = moor.vars['chipod'];
    chipods['check_time'] = chipods.time;
    grouped = chipods.groupby('time.'+group_by);
    good_fraction = availability(grouped);  
    return good_fraction



