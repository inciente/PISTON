import netCDF4, datetime
import pandas as pd
import matplotlib.pyplot as plt
import xarray as xr
import numpy as np
import scipy.interpolate as scinterp

''' Functions to import and plot data from Kerry Emanuel's TC database '''

Emanuel_folder = '/mnt/sda1/SciData/TCs/Emanuel_db/'
Emanuel_paths = ['attracks.nc','eptracks.nc','gbtracks.nc','iotracks.nc',
        'shtracks.nc','wptracks']; 

Emanuel_identifiers = ['yearic','stnamec'];
Emanuel_vars = ['longmc','latmc','vsmc','pcmc'];

class TC_collection(): 
    basins_loaded = [];
    files = []; # all instances will share files to avoid memory redundancies

    def __init__(self,title='default'):
        self.title = title;
        self.storms = []; # list of storm_info instances
        self.basins = []; # basins for which 
        
    def in_basin(self,query_basin):
        # Return all storms in collection within query_basin
        basin_list = [TC.basin for TC in self.storms];
        check_basin = basin_list == query_basin;
        TCs_in_basin = self.storms[check_basin];
        return TCs_in_basin
    
    def add_storms(self, basin, q_years=None, q_names=None):
        ''' Add an indefinite number of storms from a single 
        basin that is given as input. '''
        if basin in self.basins:
            # Has the data for this basin already been loaded?
            check_which = self.basins.index(basin)
            ncfile = self.files[check_which]
        else:
            ncfile = load_basin(basin);
            self.files.append(ncfile); 
            self.basins_loaded.append(basin); # keep track of what has been loaded

        # Find storms that match with q_years and q_names
        storm_check = find_storm(ncfile, q_years, q_names); 
        storm_indices = [i for i, x in enumerate(storm_check) if x]; # index if True
        for TC_id in storm_indices:
            # Instantite storm class with sourcefile and index within it
            self.storms.append( storm_info(ncfile, TC_id))
        return 

def load_basin(basin):
    ''' Get ncfile datastructures for all TCs in a given basin
    plus a list of identifying variables (from Emanuel_identifiers)
    that are determined by kwargs'''
    if basin in ['wp','sh','at','io']:
        track_file = Emanuel_folder + basin + 'tracks.nc'; 
    else:
        print('Basin ' + basin + ' not found');
        return
    # Import into pandas dataframe
    ncfile = xr.open_dataset(track_file) 
    # reformat storm names to strings
    ncfile['names'] = stname_to_str(ncfile['stnamec']);
    return ncfile

def stname_to_str(stnames):
    ''' Take in the variable ncfile['stnamec'] and manipulate
    to turn char into a list of strings.'''
    storm_names = [];
    for stn in range(stnames.shape[1]):
        char_name = stnames[:,stn].values;
        str_name = str();
        for let in char_name:
            str_name += str(let)[2]
        storm_names.append(str_name.strip()) # strip removes white spaces
    return storm_names

def find_storm(ncfile, q_years=None, q_names=None):
    names = ncfile['names']; 
    years = ncfile['yearic'];
    if q_years is not None:
        # Search which years match with items in list q_years
        check_years = [True if kk in q_years else False for kk in years]
        #names = names[check_years]; years= years[check_years]; 
    else:
        check_years = [True for kk in years]; 
    if q_names is not None:
        # Which names match with items in list q_names?
        check_names = np.array([0 for stn in names]); 
        for storm_name in q_names:
            check_names += np.array([1 if storm_name in stn else 0 \
                    for stn in names]);
        check_names = [True if ll > 0 else False for ll in check_names];
        #print(check_names)
    else: 
        check_names = [True for stn in names]; 
    # Combine both conditions
    check_both = [True if check_years[kk] & check_names[kk] else False \
            for kk in range(len(check_years))]
    return check_both

class storm_info:
    ''' This is the object that gathers all information from database 
    for any given storm.'''
    def __init__(self, ncfile, index):
        self.ncfile = ncfile; 
        self.index = index; 
        self.name = ncfile['names'].values[index];
        self.year = ncfile['yearic'].values[index]; 
        # Cut time dimension to find valid datapoints:
        self.valid = self.ncfile['vsmc'][:,self.index].values != 0.

    def time(self):
        ''' Use self.year, self.ncfile['monthmc','daymc','hourmc'] to 
        create np.array of datetime.datetime objects '''
        months = self.ncfile['monthmc'][:,self.index].values[self.valid];
        days = self.ncfile['daymc'][:,self.index].values[self.valid]; 
        hours = self.ncfile['hourmc'][:,self.index].values[self.valid]; 
        
        time_vector = [datetime.datetime(int(self.year), int(months[kk]), 
                  int(days[kk]), int(hours[kk])) for kk in range(len(months))];
        return np.array(time_vector)
    
    def coords(self, as_df=True):
        lon = self.ncfile['longmc'][:,self.index].values[self.valid];
        lat = self.ncfile['latmc'][:,self.index].values[self.valid]; 
        xy = {'lon':lon,'lat':lat};
        if as_df:
            xy = pd.DataFrame(data=xy, index=self.time()); 
        return xy

    def windspeed(self, as_df=True):
        wspd = self.ncfile['vsmc'][:,self.index].values[self.valid];
        if as_df:
            wspd = pd.DataFrame(data=wspd, index=self.time(), columns='Windspeed');
        else:
            pass
        return wspd

    def interp_track(self,dim,vals):
        ''' Interpolate the TC track onto the dimension dim 'lon','lat','time' using
        at the given values vals'''
        track = self.coords(); 
        base_time = datetime.datetime(1900,1,1,0); 
        numeric_time = track.index - base_time;
        numeric_time = numeric_time.total_seconds()/3600
        if dim == 'lon':
            interp_lat = scinterp.interp1d(track['lon'],track['lat'], 
                    kind='linear');
            interp_time = scinterp.interp1d(track['lon'],numeric_time,
                    kind='linear');
            interpolator = [interp_lat, interp_time]; 
        elif dim=='lat':
            interp_lon = scinterp.interp1d(track['lat'], track['lon'],
                    kind='linear'); 
            interp_time = scinterp.interp1d(track['lat'], numeric_time,
                    kind='linear'); 
            interpolator = [interp_lon, interp_time]; 
        # Get at query value. Times will be returned in hours since 1900, 01, 01.
        qvals = [interpolator[jj](vals) for jj in range(len(interpolator))];
        if dim in ['lat','lon']:
            # translate time to datetime
            qvals[1] = base_time + datetime.timedelta(hours=int(qvals[1]))

        return qvals


