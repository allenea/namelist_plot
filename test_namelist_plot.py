#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec  8 19:06:22 2018

Eric Allen, University of Delaware

Last Updated: 2/4/2019 11am

SHOWING SEVERAL CAPABILITIES OF namelist_plot class wps_info

Written in Python3(.6)
        Python/3.6
        Cartopy/0.17
        Matplotlib/3.0

@author: allenea

RUNS FROM THE MAIN namelist_plot.py class

### SAMPLE 1   - Plot All Domains
### SAMPLE 2   - Plot a Single Domain
### SAMPLE 3   - ZOOM IN ON A SINGLE DOMAIN
### SAMPLE 4   - ZOOM OUT ON A SINGLE DOMAIN
### SAMPLE 5   - Not sure where I am going with this but showing how you can use the data
### SAMPLE 6   - Using an invalid namelist.wps name  



Issues:
    - Id' rather have the x/y ticks/labels at certain increments like 10,20,30, or 15,20,25,30
    - Gridlines look like crap and don't cover the entire image
    
    

Sample 1: I am not sure why the gridlines do not cover teh entire image
Sample 2,4: I am not sure why the xticks and yticks aren't labeled.

Is it possible to have gridlines and x/y ticks on increments of 1,5,10,etc. instead of based on the extent of the projection
"""
import namelist_plot as nplt

#%% NOT USED DIRECTLY
import numpy as np
import cartopy
from cartopy.feature import OCEAN, LAKES, LAND, NaturalEarthFeature
import matplotlib.pyplot as plt

#%% SAMPLE 1        - Plot All Domains
print("Sample 1")

# Include the directory/path is file located outside the working directory
infile = 'namelist.wps'
wps = nplt.wps_info(infile)

wpsproj, latlonproj, corner_lat_full, corner_lon_full, length_x, length_y = wps.calc_wps_domain_info()
## SET UP PLOT
fig1 = plt.figure(figsize=(10,10))
ax1 = plt.axes(projection=wpsproj)

#PLOT ALL DOMAINS
nplt.wps_info.plot_all_domains(ax1,fig1,wps.max_dom, wpsproj, latlonproj, corner_lat_full, corner_lon_full, length_x, length_y)
ax1.set_title('WRF 3-Nested Domain', size=20)
fig1.savefig('WRF_DOMAIN_PLOT_COAWST.png', dpi=600)
plt.show()



#%% ## SAMPLE 2         - Plot a Single Domain
print("Sample 2")

## SET INPUT DATA
wpsfile = 'namelist.wps'
plot_domains= 3

## GET/STORE WPS Data, Print Data, Get Plot Domain, Calculate Domain Info
wps = nplt.wps_info(wpsfile)
#wps.print_info()
plt_idx = wps.plot_domain_number(plot_domains)
wpsproj, latlonproj, corner_lat_full, corner_lon_full, length_x, length_y = wps.calc_wps_domain_info()


## SET UP PLOT
fig2 = plt.figure(figsize=(10,10))
ax2 = plt.axes(projection=wpsproj)

## ADDING FEATURES
ax2.coastlines('10m','black')
ax2.add_feature(cartopy.feature.STATES.with_scale('10m'))

## REPORJECT
corner_x, corner_y = nplt.wps_info.reproject_corners(corner_lon_full[plt_idx,:], corner_lat_full[plt_idx,:], wpsproj, latlonproj)

### ZOOM FUNCTION
ns_zoom =0  
we_zoom = 0
corner_x, corner_y = nplt.wps_info.plot_zoom_in(corner_x, corner_y,ns_zoom,we_zoom)

########

## SET DOMAIN LIMITS TO ZOOM 
ax2.set_xlim([corner_x[0], corner_x[3]])
ax2.set_ylim([corner_y[0], corner_y[3]])

#print(list(ax2.get_extent(cartopy.crs.PlateCarree())))
nplt.wps_info.set_lambert_ticks(ax2,xskip=1.,yskip=1.,x_thickness=14,y_thickness=14)

ax2.set_title("WRF domain "+ str(plot_domains), size=20)
#fig2.savefig('WRF_SAMPLE_SINGLE_DOMAIN_PLOT.png', dpi=600)
plt.show()
#print("I'm not sure why the xtick and yticks aren't labeled.")

#%% SAMPLE 3   - ZOOM IN ON A SINGLE DOMAIN
print("Sample 3")

plt_idx = wps.plot_domain_number(plot_domains)
wpsproj, latlonproj, corner_lat_full, corner_lon_full, length_x, length_y = wps.calc_wps_domain_info()

## SET UP PLOT
fig3 = plt.figure(figsize=(10,10))
ax3 = plt.axes(projection=wpsproj)

## ADDING FEATURES
ax3.coastlines('10m','black')
ax3.add_feature(cartopy.feature.STATES.with_scale('10m'))

## REPORJECT
corner_x, corner_y = nplt.wps_info.reproject_corners(corner_lon_full[plt_idx,:], corner_lat_full[plt_idx,:], wpsproj, latlonproj)

### ZOOM FUNCTION
ns_zoom = 65  
we_zoom = 80
corner_x, corner_y = nplt.wps_info.plot_zoom_in(corner_x, corner_y,ns_zoom,we_zoom)
########

## SET DOMAIN LIMITS TO ZOOM 
ax3.set_xlim([corner_x[0], corner_x[3]])
ax3.set_ylim([corner_y[0], corner_y[3]])


ax3.set_title("WRF Domain "+ str(plot_domains)+" Zoom In", size=20)
#fig3.savefig('WRF_SAMPLE_SINGLE_DOMAIN_ZOOM_IN_PLOT.png', dpi=600)
plt.show()

#%% SAMPLE 4   - ZOOM OUT ON A SINGLE DOMAIN
print("Sample 4")

plt_idx = wps.plot_domain_number(plot_domains)
wpsproj, latlonproj, corner_lat_full, corner_lon_full, length_x, length_y = wps.calc_wps_domain_info()

## SET UP PLOT
fig4 = plt.figure(figsize=(10,10))
ax4 = plt.axes(projection=wpsproj)

## ADDING FEATURES
ax4.coastlines('10m','black')
ax4.add_feature(cartopy.feature.STATES.with_scale('10m'))

## REPORJECT
corner_x, corner_y = nplt.wps_info.reproject_corners(corner_lon_full[plt_idx,:], corner_lat_full[plt_idx,:], wpsproj, latlonproj)

### ZOOM FUNCTION - FULL ZOOM OUT
ns_zoom = 39 
we_zoom =39 
corner_x, corner_y = nplt.wps_info.plot_zoom_out(corner_x, corner_y,ns_zoom,we_zoom)

ax4 = nplt.wps_info.set_lambert_ticks(ax4,xskip=10.,yskip=10.,x_thickness=14,y_thickness=14)

## SET DOMAIN LIMITS TO ZOOM 
ax4.set_xlim([corner_x[0], corner_x[3]])
ax4.set_ylim([corner_y[0], corner_y[3]])

ax4.set_title("WRF Domain "+ str(plot_domains)+" Zoom Out", size=20)
#fig4.savefig('WRF_SAMPLE_SINGLE_DOMAIN_ZOOM_OUT_PLOT.png', dpi=600)
plt.show()

#%%  SAMPLE 5   - Not sure where I am going with this but showing how you can use the data
print("Sample 5")
wps = nplt.wps_info(wpsfile,True)
#dx = wps.dx #% Distance between 2 points in x direction
#dy = wps.dy #% Distance between 2 points in y direction
nx = np.empty(wps.max_dom)
ny = np.empty(wps.max_dom)
for i in range(wps.max_dom):
    nx[i] = wps.e_we[i] - wps.j_parent_start[i] #x direction #nx = 256 #% Number of grid points in the x direction
    ny[i] = wps.e_sn[i] - wps.i_parent_start[i] # y direction #ny = 256 #%*2; % Number of grid points in the y direction
    
    
print()  
print()
print("J Parents", wps.j_parent_start)
print("I Parents", wps.i_parent_start)
print("nx,ny: ",nx,ny)
print()
print()


#%% SAMPLE 6   - Using an invalid namelist.wps name  
print("Sample 6 - Failure")
infile = 'namelist.wps_failure'
wps = nplt.wps_info(infile)
wpsproj, latlonproj, corner_lat_full, corner_lon_full, length_x, length_y = wps.calc_wps_domain_info()

