import numpy as np; import xarray as xr; import pandas as pd; 
import matplotlib.pyplot as plt; import sys
sys.path.append('/mnt/sda1/PhysOc/modview/modview');
import viztools, mapper

import make_ch2_results.py as results


def make_F1():
    F1_dict = {'figsize':[10,6], 'widths':[1,0.1,2], 'heights':[0.15,2,1,1,2],
            'panels':( [slice(1,3,1), 0], [slice(3,5,1),0], 
                [slice(0,2,1),2], [slice(2,4,2), 2], [4,2], 
                [0,0] ), 
            'projections':['PlateCarree','PlateCarree', None, None, None,None], 
                'map_centers':([130, 18], [130, 18], None, None, None, None),
                'xlabels':(None,'Longitude',None, None, 'Date',None),
                'ylabels':('Latitude','Latitude',r'$\Pi$ [W m$^{-2}$]',
                    r'$\left< \| \vec{S} \| - \left< \| \vec{S} \| \right> \right>$'\
                            +'[10$^{-2}$ s$^{-1}$]',r'$\kappa$ [m s$^{-1}$]',None)}
    F1 = viztools.panel_plot(F1_dict);
    F1.draw();
    # ----------  Maps of cumulative windwork
    wwork, ww1, ww2 = import_windwork()
    F1.axes[0].set_title('12 Aug to 11 Sep')
    draw_int_wwork(F1, 0, ww1, lonlabels=False);
    draw_TC_tracks(F1, 0, L1_tracks); # soulik, jebi, cimaron
    draw_ship_track(F1, 0, [0]); # ship track for drifting period
    F1.axes[1].set_title('12 Sep to 11 Oct')
    workmap = draw_int_wwork(F1, 1, ww2, lonlabels=True); 
    draw_TC_tracks(F1, 1, L2_tracks); 
    draw_ship_track(F1, 1, [1,2]); 
    wwork_colbar(F1,-1,workmap); 
    # ---------- Time series
    plot_wwork_epsi(F1,2,wwork['windwork']) # \Pi and KE_{ni}
    draw_periods( F1.axes[2], 1.1, 0.2 )
    plot_shear_var(F1,3); # shear std
    plot_kap_jq(F1,4); # kappa and J_q
    plot_lims = {'t0':'2018-08-18 00','t1':'2018-10-15 00'};
    for panel in [2,3,4]:
        ch2funcs.dateticks(F1,panel,[5,15]);
        ch2funcs.mkxlim(F1,panel,plot_lims);
        F1.axes[panel].grid(True, which='minor',lw=0.75)
        F1.axes[panel].grid(True,which='major',lw=1.5);
        if panel in [2,3]:
            F1.axes[panel].set_xticklabels([])
    return F1

def make_F2():
    F2_dict = {'figsize':[8,8], 'widths':[1,1,1],'heights':[1,1,1],
            'panels':([0,0],[1,0],[2,0],[0,1],[1,1],[2,1],[0,2],[1,2],
                [2,2]), 'projections':[None]*9};
    F2 = viztools.panel_plot(F2_dict); 
    F2.draw(); 
    rot_list = [l1_rot, l2_rot, l2_rot]; # backrotated currents
    # Cycle through periods and columns
    for kk in [0,1,2]:
        limits = period_list[kk]; 
        # Get raw and backrotated (inertial) velocities
        v_cut = adcp.sel(time=slice(limits['t0'],limits['t1']));
        v_cut = v_cut.resample(time='1H').mean(); # time averaging
        v_rec = rot_list[kk].remake_data(comp='v');
        v_rec = v_rec.sel(time=slice(limits['t0'],limits['t1'])); 
        axinds = kk*3 + np.array([0,1,2]);  # where 2 plot for this period
        # Plot raw velocity
        F2.axes[axinds[0]].pcolormesh( v_cut.time, v_cut['depth'], 
                v_cut['v'].transpose(), cmap=cvel, vmin=-0.6, vmax=0.6);
        # Recreate complex-demodulated 
        F2.axes[axinds[1]].contourf( v_rec.time, v_rec.pressure, v_rec.values, 
                cmap=cvel, levels=np.arange(-0.25,0.28,0.05) );
        epsi = Epsilon[limits['t0']:limits['t1']]; 
        epsi = epsi.resample('1H').mean(); 
        F2.axes[axinds[2]].contourf( epsi.index, np.arange(2.5,401,5), 
                np.log10(epsi).transpose(), cmap=ceps, 
                levels=np.arange(-8.5,-5,0.5), extend='both'); 
    for kk in range(9):
        F2.axes[kk].set_ylim([350,0]);
        viztools.dateticks(F2.axes,kk,[1,4])

