import numpy as np; import xarray as xr
import pandas as pd; import scipy.io as sio
import sys, warnings
import importlib; import matplotlib.pyplot as plt
import cartopy as cart
import matplotlib.dates as mdates
import matplotlib as mpl
import matplotlib.patheffects as PathEffects
from matplotlib.patches import Rectangle
from datetime import datetime
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
sys.path.insert(0,'/mnt/sda1/PhysOc/modview/modview')
import phystools, timetools

mybox = [[133,137],[10,19]]; # for windwork visualization

def plot_range(df, limits, goodz, axis, reps=500, cut=True):
    ''' Take in data frame, compute mean and 95% CL, and return both. '''
    if cut:
        df_new = df.loc[limits['t0']:limits['t1']]
        df_new = df_new.loc[:,goodz]; 
    else:
        df_new = df; 
    df_bounds = bootstrap_matrix(df_new, reps, func=np.nanmean);
    with warnings.catch_warnings():
        warnings.simplefilter("ignore",category=RuntimeWarning);
        if axis==0:
            df_line = np.nanmean(df_new,axis=1); 
        elif axis==1:
            df_line = np.nanmean(df_new,axis=0);  
    return df_line, df_bounds
    

def bootstrap(data_vec,reps, ivals=[2.5,97.5], func=np.mean):
    ''' Resample a single vector data_vec reps number of times. 
    Then, apply func  to estimate confidence intervals of that 
    specific function, which can be std, mean, var, or other.
    '''
    dat_replicates = np.empty(reps); 
    with warnings.catch_warnings():
        warnings.simplefilter("ignore",category=RuntimeWarning);
        for kk in range(reps):
		    # artificial sample of data
            dat_sample = np.random.choice(data_vec,size=len(data_vec)); 
            dat_replicates[kk] = func(dat_sample); # save mean
    lo_bound = np.percentile( dat_replicates, [ivals[0]]); 
    hi_bound = np.percentile( dat_replicates, [ivals[1]]); 
    return lo_bound, hi_bound 

def bootstrap_matrix(data_mat,reps,ivals=[2.5,97.5],axis=1, func=np.mean):
    ''' Take in a matrix data_mat and cycle through its columns/rows 
    to calculate bootstrapping error bars via boot strapping on each column. 
    This way, a different error bar is estimated for each depth level in vertical profiles.
    '''
    N_est = data_mat.shape[axis];
    bounds = np.empty((N_est,2));
    
    for jj in range(N_est):
        if axis==0:
            dat_vector = data_mat.values[jj,:];  
        elif axis==1:
            dat_vector = data_mat.values[:,jj];
        bounds[jj,0], bounds[jj,1] = bootstrap(dat_vector, reps, ivals, func); 
    return bounds
    
    
def inset_ax(axis, locat):
    ''' Create box within a place a separate visual object there '''
    ans1 = inset_axes(axis, width="3%",  height="90%",  loc=locat,
                     bbox_to_anchor=(0.07,0.,1,1), bbox_transform=axis.transAxes, borderpad=0)
    return ans1

def SeaMap(fig_dict, in_gs, axind):
    ''' Make a map of coasts for the extent specified below ''' 
    figname = fig_dict['fig']; 
    gridplace = fig_dict['gs'][in_gs];
    axlist = fig_dict['axes']
    #set up map
    axlist[axind] = figname.add_subplot( gridplace, 
                projection=cart.crs.Orthographic(
                central_longitude=135, central_latitude=18))
    axlist[axind].set_extent([118,142,5,28]);
    axlist[axind].set_aspect('equal')
    axlist[axind].add_feature(cart.feature.LAND,zorder=100,
                    edgecolor='k', facecolor=[0.8,0.8,0.8]);
    return axlist
    
def dashedbox(fig_dict, axind):
    axlist = fig_dict['axes']; 
    for kk in axind:
        axlist[kk].plot([mybox[0][0],mybox[0][1],mybox[0][1],
                mybox[0][0], mybox[0][0]],[mybox[1][0],mybox[1][0],
                mybox[1][1],mybox[1][1],mybox[1][0]],
                linestyle='dashed',linewidth=2,color='blue',
                      transform=cart.crs.PlateCarree())

def box_avg(xr_obj, variable):
    line2plot = xr_obj[variable].sel(lon=slice(mybox[0][0],mybox[0][1]))
    line2plot = line2plot.sel(lat=slice(mybox[1][1],mybox[1][0]))
    line2plot = line2plot.mean(['lon','lat']);
    line2plot = line2plot.sel(time=slice('2018-08-11','2018-10-12'));
    return line2plot
    
# Create grid spec using figure information
def make_gs(fig_dict): 
    nx = len(fig_dict['widths']); ny = len(fig_dict['heights']); 
    gs = fig_dict['fig'].add_gridspec(ncols=nx, nrows=ny, 
                    width_ratios= fig_dict['widths'], 
                    height_ratios=fig_dict['heights']); 
    return gs

def make_panels(fig_dict):
    panlist = fig_dict['panels']; 
    axes = fig_dict['axes'];
    jj=0; 
    for kk in panlist:
        axes[jj] = fig_dict['fig'].add_subplot( fig_dict['gs'][kk[0],kk[1]]);    
        jj+=1
    return axes
    
def find_goodz(cham_obj, limits, threshold = 0.25):
    varhere = cham_obj[limits['t0']:limits['t1']]; 
    goodcount = np.sum( ~np.isnan(varhere), axis=0); 
    goodz = (goodcount/varhere.shape[0])>threshold;
    return goodz
    
def cham_fig_format(axes):
    axes[0].loglog([1e-5,1e-4],[1e2,2e2],color='black',linestyle='dashed',linewidth=2 , alpha=0.7,
               label='$c_g$ up')
    axes[0].loglog([1e-5,1e-4],[1e2,3e3],color='black',linestyle='solid',linewidth=3,label='$c_g$ down')
    axes[0].legend(); axes[1].legend()
    axes[0].grid(True, which='major',ls="-", color='0.65')
    axes[0].grid(True, which='minor',ls="-", color='0.8', alpha=0.5)
    axes[0].set_xlim((1e-3, 6e-2)); axes[0].set_ylim([0.85e-4,2.35e-3]);
    
    # Set special ticks
    axes[1].grid(True, which='both');
    ticks = np.array([0,5e-5,5e-4,2.5e-3,1e-2,5e-2])**(1/4)
    axes[1].set_yticks(ticks); axes[1].set_ylim((-(1e-6)**0.25, (5e-2)**0.25))
    axes[1].set_yticklabels(['0','5','50','250','1000','5000'])
    
    axes[2].set_yticks(np.arange(0,360,50));
    axes[2].set_yticklabels(['0','','100','','200','','300',''])
    axes[4].set_xticks([-5,-4,-3,-2])
    axes[5].set_xlim([-12.5,80])
    hf_ticks = np.array([0,25,50,75]);
    axes[5].set_xticks(hf_ticks )
    
    xlabs = [r'$k_z/2\pi$ [cpm]',r'$J_q$ [W m$^{-2}$]',
           r'$ \frac{\rho_0}{2} \left< \| \vec{U} \|^2\right>$ [J m$^{-3}$]',
           r'$ \left< \epsilon \right> $ [W kg$^{-1}$]',r'$\left< \kappa \right>$ [m$^2$ s$^{-1}$]',
           r'$\left< J_q \right>$ [W m$^{-2}$]']
    ylabs = [r'$\Phi_{sh}$ [m s$^{-2}$]',r'pdf [$10^{-5}$ W$^{-1}$ m$^2$]','Depth [m]','','','']

    for kk in range(len(axes)):
        axn = axes[kk]
        axn.set_xlabel(xlabs[kk])
        axn.set_ylabel(ylabs[kk])

    
def period_spectra(dataset, limits, ctd_dict, scaling=1, dz=10):
    # Take a slice of adcp data to compute spectra
    block = dataset['u'].sel(time=slice(limits['t0'],limits['t1']));
    block = block + 1j*dataset['v'].sel(time=slice(limits['t0'],limits['t1']));
    block = block * scaling; # wkb scaling
    block = block.sel(depth_cell=slice(limits['z0'],limits['z1']));
    # Prepare N2 
    block_n2 = np.interp(block['depth'], ctd_dict['p_mid'], ctd_dict['N2'])
    c_scale, z_wkb = phystools.wkb_stretch(block['depth'], block_n2, 30) # stretch coordinate
    z_hom = np.arange(min(z_wkb), max(z_wkb)+1, dz); # give it even spacing
    # Create xr version of block using wkb depth
    mod_block = xr.DataArray(data=block,dims=['time','depth'], coords={\
                  'time':block.time, 'depth':z_wkb})
    mod_block = mod_block.differentiate(coord='depth')
    mod_block = mod_block.sel(depth=z_hom,method='nearest'); # velocity on evenly spaced grid
    nan_count = np.sum(np.isnan(mod_block),axis=1); # how many nans per profile
    
    mod_power = timetools.spectrum_1D(mod_block, dt=dz, nseg=1, axis=1)
    wvnums = np.fft.fftfreq(n=mod_block.shape[1],d=dz)
    # Remove spectra whose data had more than 30% nans.
    mod_power[nan_count/mod_block.shape[0]>0.3, :] = np.nan; 
    mod_power = np.nanmean( mod_power, axis=0)
    # Separate frequencies
    pos_power = mod_power[wvnums>0];
    neg_power = mod_power[wvnums<0];
    mod_power = [pos_power, neg_power];
    wvnums = [wvnums[wvnums>0], wvnums[wvnums<0]];
    return mod_power, wvnums

def dateticks(fig_dict, axind, dts, combined=False):
    axlist = fig_dict['axes']
    # Major ticks.
    fmt_big = mdates.DayLocator(interval=dts[1])
    axlist[axind].xaxis.set_major_locator(fmt_big)
    # Minor ticks.
    if combined:
        fmt_small = mdates.HourLocator(byhour=dts[0]);
    else:
        fmt_small = mdates.DayLocator(interval=dts[0]);
    axlist[axind].xaxis.set_minor_locator(fmt_small)
    # FOrmat of text
    axlist[axind].xaxis.set_major_formatter(mdates.DateFormatter('%b %-d'))
    
def colgroup(df,ncols):
    # Take the average of contiguous groups of columns in a dataframe
    dfc = np.zeros( (df.shape[0],int(df.shape[1]/ncols)) )
    for kk in range(dfc.shape[1]):
        start = kk*ncols; 
        end = (kk+1)*ncols; 
        dfc[:,kk] = np.nanmean(df.loc[:,start:end], axis=1);
    
    dfc[dfc==0] = np.nan;
    dfc = pd.DataFrame( data=dfc, index=df.index)
    return dfc
    
def mapticks(fig_dict, axind, latlabels=True):
    # specify format of map ticks and their labels
    axlist = fig_dict['axes'];
    for kk in axind: 
        gl = axlist[kk].gridlines(crs=cart.crs.PlateCarree(), draw_labels=True,
                  linewidth=1, color='k', alpha=0.7)
        gl.top_labels = False
        gl.right_labels = False
        gl.left_labels = latlabels
        #gl.xlines = False
        gl.xlocator = mpl.ticker.FixedLocator([120, 130, 140])
        gl.xformatter = LATITUDE_FORMATTER
        gl.ylocator = mpl.ticker.FixedLocator([12, 18, 24])
        gl.yformatter = LATITUDE_FORMATTER
    
def mkxlim(fig_dict, axind,limits):
    axlist = fig_dict['axes'];
    startTime = datetime.strptime(limits['t0'],'%Y-%m-%d %H');
    endTime = datetime.strptime(limits['t1'],'%Y-%m-%d %H');
    startTime = mdates.date2num(startTime);
    endTime = mdates.date2num(endTime);
    axlist[axind].set_xlim((startTime, endTime));
    return
    
def shade_period(fig_dict, axind, limits, color='yellow'):
    axlist = fig_dict['axes']; 
    print(axlist)
    # convert to matplotlib date representation
    startTime = datetime.strptime(limits['t0'],'%Y-%m-%d %H');
    endTime = datetime.strptime(limits['t1'],'%Y-%m-%d %H');
    start = mdates.date2num(startTime)
    end = mdates.date2num(endTime)
    width = end - start;
    # Plot rectangle
    rect = Rectangle((start, 0), width, 1, color=color)
    axlist[axind].add_patch(rect)   
    
    return

def colbartop(bounds,im,fig_dict, barlegend,tickz):
    fig = fig_dict['fig']; 
    axins = fig.add_axes(bounds)
    cbar = fig.colorbar(im, cax=axins, orientation="horizontal", 
                 label=barlegend, ticklocation='top', ticks=tickz)
    print([str(kk) for kk in tickz])
    cbar.ax.set_xticklabels([str(kk) for kk in tickz])             
    return
    
def square_mean(arr, y, x):
    yy, xx = arr.shape
    vals = arr.reshape(y, yy//y, x, xx//x).mean((1,3))
    return vals
    
def format_direct_obs(fig_dict, axind, period, dt=[1,5], ylim=(350,0),yticks=False):
    axlist = fig_dict['axes']; # for access to axis
    # via other functions
    dateticks(fig_dict, axind, dt); 
    mkxlim(fig_dict,axind, period); 
    # directly on axis
    axlist[axind].set_ylim(ylim); 
    axlist[axind].set_yticks(np.arange(0,350,50)); 
    axlist[axind].set_yticklabels([]); 
    if yticks:
        axlist[axind].set_yticklabels(['0','','100','','200','','300'])
    
    
# EXPECT TO RESTRUCTURE THIS CONSIDERABLY 
def workseries(fig_dict, axind, work2plot, flux2plot1, flux2plot2):
    datcolor = [0.8, 0.26, 0.8];
    axlist = fig_dict['axes'];
    # plot turb fluxes
    newax = axlist[axind].twinx();
    newax.plot( flux2plot1.iloc[:,100:150].mean(axis=1),color=datcolor)
    newax.plot( flux2plot2.iloc[:,100:150].mean(axis=1),color=datcolor); 
    # plot windwork
    axlist[axind].plot( work2plot.time, work2plot, color='k',linewidth=2 , label='$\Pi$ dashed box')
    axlist[axind].set_ylabel('Windwork $\Pi$ [W m$^{-2}$]')
    axlist[axind].set_ylim((-0.1, 0.6)); 
    
    newax.set_ylim((-25,150)); dateticks(fig_dict,axind,[5,15])
    newax.grid(True)
    newax.grid(axis ='x', linewidth='1', color='black');
    newax.set_ylabel('$J_q$ [W m$^{-2}$]')
    mkxlim(fig_dict,axind,{'t0':'2018-08-10 00','t1':'2018-10-13 23'})
    newax.tick_params(axis='y', colors=datcolor)
    newax.yaxis.label.set_color(datcolor)
    newax.spines['right'].set_color(datcolor)
    return

    #axlist[axind].legend()

def leg2ylabel(legax):    
    # Turn off axis lines and ticks of the big subplot
    legax.spines['top'].set_color('none')
    legax.spines['bottom'].set_color('none')
    legax.spines['left'].set_color('none')
    legax.spines['right'].set_color('none')

    legax.tick_params(labelcolor='w', top=False, bottom=False, left=False, right=False)
    
def addlet(fig_dict, axind, x,y,letter):
    axlist = fig_dict['axes']; 
    for kk in range(len(axind)):
        txt = axlist[axind[kk]].text(x[kk], y[kk], letter[kk],
                transform=axlist[axind[kk]].transAxes, fontweight='bold',size=22)
        txt.set_path_effects([PathEffects.withStroke(linewidth=5, foreground='w')])
        plt.draw()

def cut_slow(p_dict, cruise_leg, var):
    #if var=='Jq':
    #    slovar = cruise_leg.turbflux('DTDZ',interval='1h')
    if var=='KAPPA':
        slovar = cruise_leg.diffusivity(interval='2h');  
    else:
        slovar = cruise_leg.slowcham(var, interval='2h'); 
    slovar = slovar[p_dict['t0']:p_dict['t1']]; 
    return slovar
     
def prepare_geost(p_dict,cruise_leg):
    p_dict['x'], p_dict['z'], p_dict['u'] = cruise_leg.cuttz(p_dict,'u','adcp'); 
    x,z,p_dict['v'] = cruise_leg.cuttz(p_dict,'v','adcp'); 
    vars2get = ['THETA','COORDS','SAL']; 
    for item in vars2get:
        p_dict[item] = cut_slow(p_dict, cruise_leg, item); 
    return p_dict

def makewave(timevec,dt,depth,kz,f,phi,amp):
    dtnum = np.arange(0, len(timevec),1)*dt;
    vertphase = np.expand_dims(kz * depth, axis=1);
    timephase = np.expand_dims(f*dtnum, axis=0);
    
    sine = amp*np.sin(vertphase - timephase + phi);
    return sine

# Calculate wkb-stretched depth for leg 2.
def wkbstretch(z,N2,MLD):
    belowML = z > MLD;
    N = np.sqrt(N2);
    sratio = np.nanmean(N[belowML])/N;
    # Rolling mean for smoothness
    sratio = pd.DataFrame(sratio); 
    sratio = sratio.rolling(10, min_periods=1).mean().to_numpy();

    # Scaling for velocities
    sfactor = np.sqrt(sratio);
    
    # Now stretch the depths
    tot_depth = z[~belowML];
    stretched = np.expand_dims( np.append( 0, np.diff(z)), axis=1);
    stretched = np.cumsum( stretched[belowML]/sratio[belowML]) + MLD;
    stretched = np.append(tot_depth, stretched);
    
    return sfactor, stretched

def make_cmaps():
    cvel = mpl.cm.get_cmap('seismic'); # velocity
    cvel = cvel(np.linspace(0,1,14));
    #cvel[6,:] = [0.9, 0.9, 0.9, 1];
    cvel = mpl.colors.ListedColormap(cvel)

    ceps = mpl.cm.get_cmap('magma_r',6); # epsilon
    ceps = ceps(np.linspace(0,1,6));
    ceps = np.vstack( ([0.85, 0.85, 0.85, 1], ceps[0:-1,:]))
    ceps = mpl.colors.ListedColormap(ceps)
    ceps.set_under([0.5, 0.5, 0.5,1]); ceps.set_over([0.25,0,0.51,1]);

    cwork = mpl.cm.get_cmap('YlGnBu'); # cumulative windwork 
    cwork = cwork(np.linspace(0,1,8));
    cwork = np.vstack( ([1, 1, 1, 1], cwork))
    cwork = mpl.colors.ListedColormap(cwork); 
    cwork.set_under([0.5, 0.5, 0.5, 1])
    return cvel, ceps, cwork
    
def make_cmaps2():
    KEmap = mpl.cm.get_cmap('viridis',10); 
    SHmap = mpl.cm.get_cmap('YlGnBu',9); 
    heatm = mpl.cm.get_cmap('OrRd'); 
    heatm = heatm(np.linspace(0.06,1,40)); 
    HFmap = mpl.colors.LinearSegmentedColormap.from_list('HFmap',heatm,N=10);
    return KEmap, SHmap, HFmap

def prepare_wwork(wwork):
    ww1 = wwork['windwork'].sel(time=slice('2018-08-11 00:00','2018-09-10 23:59'));
    ww1 = ww1.sum(dim='time')*3600*3;
    ww2 = wwork['windwork'].sel(time=slice('2018-09-11 00:00','2018-10-10 23:59'));
    ww2 = ww2.sum(dim='time')*3600*3;	
    return ww1, ww2

def gather_hycom(hycom_path):
    sst = xr.open_dataset(hycom_path+'sst_wakes.nc')
    ssh = xr.open_dataset(hycom_path+'ssh_wakes.nc');
    uv = xr.open_dataset(hycom_path+'uv_wakes.nc');
    wakes_tot = xr.merge([sst,ssh,uv]) 
    wakes = wakes_tot.sel(lon=slice(133,137)); 
    wakes = wakes.sel(lat=slice(12, 19)); # limits for visualization
    return wakes
    
def join_chamvars(leg1, leg2):
    Nsquared = colgroup( pd.concat( [leg1.cham['N2'],leg2.cham['N2']] ), 5)
    Epsilon = colgroup( pd.concat( [leg1.cham['EPSILON'],leg2.cham['EPSILON']] ), 5);
    Kappa = colgroup( pd.concat( [leg1.diffusivity(), leg2.diffusivity()]) , 5)
    HeatFlux = colgroup( pd.concat( [leg1.turbflux('DTDZ'), leg2.turbflux('DTDZ')]), 5);
    Theta = pd.concat( [leg1.cham['THETA'], leg2.cham['THETA']] );
    return Nsquared, Epsilon, Kappa, HeatFlux, Theta

def dense_grid(axobj, rangex, interval):
    xlocs = np.arange(rangex[0],rangex[1],interval);
    for kk in xlocs:
        axobj.plot([kk,kk], [50,0], color='black',alpha=0.5, linewidth=0.6); 

    
    
    
    
    
    
