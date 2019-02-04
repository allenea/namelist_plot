#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  6 19:21:36 2018

@author: allenea
Eric Allen, University of Delaware

Last Updated: 2/4/2019 11am

Written in Python3(.6)
        Python/3.6
        Cartopy/0.17
        Matplotlib/3.0
CLASS: wps_info ("namelist.wps") -- > Calls get_wps_info

METHODS:
    - get_wps_info(filename,plot_info=False) -- Default -- Intialized when calling wps_info(filename)... Optional to print info right away
    
    - print_info()
    
    - get_proj_LambertConformal()
        +(RETURNS: projection)
            
    - get_proj_Mercator()
        +(RETURNS: projection)
        
    - get_proj_Stereographic()
        +(RETURNS: projection)
        
    - get_proj_LatLon()
        +(RETURNS: projection)
        
    - plot_domain_number(domain_number)  -- default -- 1
        + (RETURNS: domain_number)
        
    - calc_wps_domain_info() **
        - (RETURNS: wpsproj, latlonproj, corner_lat_full, corner_lon_full, length_x, length_y)
        
    - calc_corner_point_latlon(center_lat, center_lon, e_we, e_ns, dx, dy, wpsproj, latlonproj, loc) **
       +(RETURNS: corner_lon, corner_lat)
       
    - calc_center_point_latlon(corner_lat_parent, corner_lon_parent, dx_parent, dy_parent, e_we, e_ns, dx, dy, i, j, wpsproj, latlonproj) **
        +(RETURNS: center_lon_child, center_lat_child)
        
    - reproject_corners(corner_lons, corner_lats, wpsproj, latlonproj) **
        +(RETURNS: return corner_x, corner_y)
        
    - plot_zoom_out(corner_x, corner_y, ns_zoom,we_zoom)
        +(RETURNS: return corner_x, corner_y)

    - plot_zoom_in(corner_x, corner_y, ns_zoom,we_zoom)
        +(RETURNS: return corner_x, corner_y)

    - plot_all_domains(ax1, fig1, wpsproj, latlonproj, corner_lat_full, corner_lon_full, length_x, length_y)
        +(RETURNS: ax1, fig1 )

    
OUTSIDE FUNCTION:
    - _ismissing(val, islat=True)  -- FROM WRF-PYTHON
        +(RETURNS: (bool) True/False)
        
        
** Adapted from Xiaodong Chen @lucas-uw on github 

Code heavily adapted from https://github.com/lucas-uw/WRF-tools/blob/master/WRF_input_tools/Visualize_WPS_domain.ipynb and combined with my previous
attempts at reading/using/plotting wrf's namelist.wps file-data.


He had a great calc_wps_domain_info function. I modified it to work with my wps_info class replacing my reprojection code.
    - Decided not to implement the lambert projection gridlines work-around
        - Wait until cartopy implements it in a future version or let the user implement at their own peril

I was having trouble recalculating the centerpoint and then getting the corners perfect in all nested domains for all projection types.
Instead of just plotting the bounds of each domain, I have made it  so that you can create a cartopy projection and use only that domain,
 then from there you can choose to zoom in and out as you wish.

You do not have to use this class for plotting, but that is one of the features... The other main feature is it hold namelist.wps data
    - You can also use it to recreate your grid for other uses such as spatial interpolation or plotting observation data.


#!! Problems/Improvements:
     - Gridlines
     - X/Y Ticks and Labels for geographic reference
     - Is it possible to improve zoom or shift center point of zoom?
     - No Testing
     - How plot_all_domains applies gridlines
     
"""

#IMPORTS
import numpy as np
import matplotlib as mpl
import cartopy.crs as ccrs
from matplotlib.font_manager import FontProperties
import math
from cartopy.feature import OCEAN, LAKES, LAND  #,NaturalEarthFeature
import cartopy
import os.path
import shapely.geometry as sgeom
from copy import copy
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.ticker as ticker

class wps_info(object):
    """
    Developed by Eric Allen, University of Delaware
    
    wps = wps_info(infile)
    plt_idx = wps.plot_domain_number(plot_domains)
    wpsproj, latlonproj, corner_lat_full, corner_lon_full, length_x, length_y = wps.calc_wps_domain_info()
    
    """
    def __init__(self,file,print_info=False):
        self.filelocname = file
        self._printinfo = print_info
        if os.path.exists(self.filelocname):
            wps_info.get_wps_info(self.filelocname)
            if self._printinfo == True:
                wps_info.print_info()
        else:
            print ("FILE NOT FOUND IN PROVIDED PATH AND/OR FILENAME")
            raise ValueError

        
    @classmethod
    def get_wps_info(cls,namelist_wps_filename="namelist.wps"):
        """
        Reads the namelist file and stores the important information. Default = namelist.wps
        
        """

        newdx = []
        newdy = []
        cls.pole_lat = 90.0 # these get changed if there is something to change it with
        cls.pole_lon = 0.0  # these get changed if there is something to change it with
        file = open(namelist_wps_filename,'r').readlines()
        for row in file:
            row = " ".join(row.split(","))
            row = " ".join(row.split("\n"))
        
            if 'wrf_core' in row:
                tmprow = row.split("=")[1]
                tmprow = " ".join(tmprow.split("'"))
                cls.wrf_core = str(tmprow.strip())
                
            elif 'max_dom' in row:
                tmprow = row.split("=")[1]
                cls.max_dom = int(tmprow)

            elif 'start_date' in row:
                tmprow = row.split("=")[1] #All Start/End Dates for each domain should be the same
                tmprow = " ".join(tmprow.split("'"))
                tmprow2 = tmprow.split(" ")
                start1 = [item for item in tmprow2 if item != '']
                cls.start_date = start1[0]
        
            elif 'end_date' in row:
                tmprow = row.split("=")[1] #All Start/End Dates for each domain should be the same
                tmprow = " ".join(tmprow.split("'"))
                tmprow2 = tmprow.split(" ")
                end1 = [item for item in tmprow2 if item != '']
                cls.end_date = end1[0]
                
            elif 'interval_seconds' in row:
                tmprow = row.split("=")[1]
                cls.interval_seconds = int(tmprow)
        
            elif 'parent_id' in row:
                tmprow = row.split("=")[1]
                tmprow2 = tmprow.split(" ")
                pid = [int(item) for item in tmprow2 if item != '']
                cls.parent_id = pid
                
            elif 'parent_grid_ratio' in row:
                tmprow = row.split("=")[1]
                tmprow2 = tmprow.split(" ")
                pgr = [int(item) for item in tmprow2 if item != '']
                cls.parent_grid_ratio = pgr
                
            elif 'i_parent_start' in row:
                tmprow = row.split("=")[1]
                tmprow2 = tmprow.split(" ")
                ipar = [int(item) for item in tmprow2 if item != '']
                cls.i_parent_start = ipar
                
                        
            elif 'j_parent_start' in row:
                tmprow = row.split("=")[1]
                tmprow2 = tmprow.split(" ")
                jpar = [int(item) for item in tmprow2 if item != '']
                cls.j_parent_start = jpar
                
                
            elif 'e_we' in row:
                tmprow = row.split("=")[1]
                tmprow2 = tmprow.split(" ")
                ewe = [int(item) for item in tmprow2 if item != '']
                cls.e_we = ewe
                
            elif 'e_sn' in row:
                tmprow = row.split("=")[1]
                tmprow2 = tmprow.split(" ")
                esn = [int(item) for item in tmprow2 if item != '']
                cls.e_sn = esn
                
            elif 'geog_data_res' in row:
                tmprow = row.split("=")[1]
                tmprow = " ".join(tmprow.split("'"))
                tmprow2 = tmprow.split(" ")
                gdr = [item for item in tmprow2 if item != '']
                cls.geog_data_res = gdr
                
            elif 'dx' in row:
                tmprow = row.split("=")[1]
                tmprow2 = tmprow.split(" ")
                dxval = [float(item) for item in tmprow2 if item != '']
                dx = dxval[0]
                for dom_idx in range(cls.max_dom):
                    if dom_idx == 0:
                        tmpdx = dx
                        newdx.append(tmpdx)
                    else:
                        tmpdx =tmpdx/cls.parent_grid_ratio[dom_idx]
                        newdx.append(tmpdx)
                cls.dx0 = newdx[0]
                cls.dx = newdx
        
                
            elif 'dy' in row:
                tmprow = row.split("=")[1]
                tmprow2 = tmprow.split(" ")
                dyval = [float(item) for item in tmprow2 if item != '']
                dy = dyval[0]
                for dom_idx in range(cls.max_dom):
                    if dom_idx == 0:
                        tmpdy = dy
                        newdy.append(tmpdy)
                    else:
                        tmpdy =tmpdy/cls.parent_grid_ratio[dom_idx]
                        newdy.append(tmpdy)
                cls.dy0 = newdy[0]
                cls.dy = newdy
                
            elif 'map_proj' in row:
                tmprow = row.split("=")[1]
                tmprow = "".join(tmprow.split("'"))
                cls.map_proj = str(tmprow.strip())
        
            elif 'ref_lat' in row:
                tmprow = row.split("=")[1]
                tmprow2 = tmprow.split(" ")
                rlat = [float(item) for item in tmprow2 if item != '']
                cls.ref_lat = rlat[0]
                
            elif 'ref_lon' in row:
                tmprow = row.split("=")[1]
                tmprow2 = tmprow.split(" ")
                rlon = [float(item) for item in tmprow2 if item != '']
                cls.ref_lon = rlon[0]
                
            elif 'truelat1' in row:
                tmprow = row.split("=")[1]
                tmprow2 = tmprow.split(" ")
                tl1 = [float(item) for item in tmprow2 if item != '']
                cls.truelat1 = tl1[0]
                
            elif 'truelat2' in row:
                tmprow = row.split("=")[1]
                tmprow2 = tmprow.split(" ")
                tl2 = [float(item) for item in tmprow2 if item != '']
                cls.truelat2 = tl2[0]
                
            elif 'stand_lon' in row:
                tmprow = row.split("=")[1]
                tmprow2 = tmprow.split(" ")
                stdlon = [float(item) for item in tmprow2 if item != '']
                cls.stand_lon = stdlon[0]
                
            elif 'ref_x' in row:
                tmprow = row.split("=")[1]
                tmprow2 = tmprow.split(" ")
                refx = [float(item) for item in tmprow2 if item != '']
                cls.ref_x = refx[0]
                
            elif 'ref_y' in row:
                tmprow = row.split("=")[1]
                tmprow2 = tmprow.split(" ")
                refy = [float(item) for item in tmprow2 if item != '']
                cls.ref_y = refy[0]
                
            elif 'pole_lat' in row:
                tmprow = row.split("=")[1]
                tmprow2 = tmprow.split(" ")
                polelat = [float(item) for item in tmprow2 if item != '']
                cls.pole_lat = polelat[0]
                
            elif 'pole_lon' in row:
                tmprow = row.split("=")[1]
                tmprow2 = tmprow.split(" ")
                polelon = [float(item) for item in tmprow2 if item != '']
                cls.pole_lon = polelon[0]
    
        #keys = ["wrf_core", "max_dom","start_date","end_date","interval_seconds","parent_id","parent_grid_ratio",\
        #             "i_parent_start","j_parent_start","e_we","e_sn","geog_data_res","dx","dy","map_proj","ref_lat","ref_lon",\
        #             "truelat1","truelat2","stand_lon","ref_x","ref_y","pole_lat","pole_lon","new_dx","new_dy"]
        # data = [ wrf_core, max_dom,start_date,end_date,interval_seconds,parent_id,parent_grid_ratio,\
        #             i_parent_start,j_parent_start,e_we,e_sn,geog_data_res,dx,dy,map_proj,ref_lat,ref_lon,\
        #             truelat1,truelat2,stand_lon,ref_x,ref_y,pole_lat,pole_lon,newdx,newdy]
        #dictionary = dict(zip(keys, data))
    
        #return dictionary
    

    @classmethod     
    def print_info(cls):
        print ("wrf_core:", cls.wrf_core)
        print ("max_dom:", cls.max_dom)
        print ("start_date:", cls.start_date)
        print ("end_date:", cls.end_date)   
        print ("interval_seconds:", cls.interval_seconds)
        print ("parent_id:", cls.parent_id)   
        print ("parent_grid_ratio:", cls.parent_grid_ratio)   
        print ("i_parent_start:", cls.i_parent_start)   
        print ("j_parent_start:", cls.j_parent_start)   
        print ("e_we:", cls.e_we)   
        print ("e_sn:", cls.e_sn)   
        print ("geog_data_res:", cls.geog_data_res)   
        print ("dx0 =",cls.dx0)
        print ("dx =",cls.dx)
        print ("dy0 =",cls.dy0)
        print ("dy =",cls.dy)
        print ("map_proj:", cls.map_proj)
        print ("ref_lat:", cls.ref_lat)   
        print ("ref_lon:", cls.ref_lon)   
        print ("truelat1:", cls.truelat1)   
        print ("truelat2:", cls.truelat2)   
        print ("stand_lon:", cls.stand_lon)   
        print ("ref_x:", cls.ref_x)   
        print ("ref_y:", cls.ref_y)   
        print ("pole_lat:", cls.pole_lat)         
        print ("pole_lon:", cls.pole_lon) 
 
    #%% Projection Class-Methods
    @classmethod
    def get_proj_LambertConformal(cls):
        lccproj = ccrs.LambertConformal(central_longitude=cls.ref_lon, central_latitude=cls.ref_lat,
                                        standard_parallels=(cls.truelat1, cls.truelat2), globe=None, cutoff=10)
        return lccproj
     
    @classmethod #UNTESTED
    def get_proj_Mercator(cls):
        merproj = ccrs.Mercator(central_longitude=cls.ref_lon, globe=None)
        return merproj
    

    @classmethod #UNTESTED
    def get_proj_Stereographic(cls):
        hemi = -90. if cls.truelat1 < 0 else 90.
        
        lat_ts = (None 
                  if wps_info._ismissing(cls.truelat1) 
                  else cls.truelat1)
                
        polarproj = ccrs.Stereographic(central_latitude=hemi, 
                                          central_longitude=cls.stand_lon, 
                                          true_scale_latitude=lat_ts, 
                                          globe=None)
        return polarproj
    
    
    @classmethod #UNTESTED
    def get_proj_LatLon(cls):
        latlonproj = ccrs.PlateCarree(central_longitude=cls.stand_lon,
                                            globe=None)
        return latlonproj
    
    
    #%% Calculating and Plotting WPS Domain    
    @classmethod
    def plot_domain_number(cls, number=1):
        """ Developed by Eric Allen to plot certain domains from namelist.wps file
        
        Optional Input - Default to coarsest domain, but configurable to the max_domain
        """
        if cls.max_dom < number:
            print ("Domain Number Out Of Range: Setting to coarsest domain")
            cls._plot_domain_number = 0
        else:
            cls._plot_domain_number = number - 1
            
        cls._maxdom_LIMIT = 10
        if cls.max_dom > cls._maxdom_LIMIT:
            print ("USING TOO MANY DOMAINS. Re-evaluate your model setup ")
            cls._plot_domain_number = 0
            
        return cls._plot_domain_number
    
    @classmethod
    def calc_wps_domain_info(cls):
         """ MODIFIED FROM Xiaodong Chen @lucas-uw on github"""
          
         # get WPS projection info
         if cls.map_proj=='lambert':
             cls.wpsproj = wps_info.get_proj_LambertConformal()
         elif cls.map_proj == "mercator":
             cls.wpsproj = wps_info.get_proj_Mercator()
         elif cls.map_proj == "polar":
             cls.wpsproj = wps_info.get_proj_Stereographic()
         elif cls.map_proj == "lat-lon":
             cls.wpsproj = wps_info.get_proj_LatLon()
         else:
            print("WARNING: PROJECTION UNKNOWN")
            raise ValueError

         # Geodetic, for lat/lon projection
         cls.latlonproj = ccrs.Geodetic()
         #cls.latlonproj = ccrs.PlateCarree()

         
         length_x = np.zeros((cls.max_dom, 1))
         length_y = np.zeros((cls.max_dom, 1))
         
         for i in np.arange(0,cls.max_dom):
            length_x[i] = cls.dx[i]*cls.e_we[i] #Units = Meters
            length_y[i] = cls.dy[i]*cls.e_sn[i] #Units = Meters
            
            
         center_lat_full = np.zeros((cls.max_dom, 1))
         center_lon_full = np.zeros((cls.max_dom, 1))
         corner_lat_full = np.zeros((cls.max_dom, 4)) # ll, lr, uf, ur
         corner_lon_full = np.zeros((cls.max_dom, 4)) # ll, lr, uf, ur

         center_lat_full[0] = cls.ref_lat
         center_lon_full[0] = cls.ref_lon

         corner_lon_full[0,0], corner_lat_full[0,0] = wps_info.calc_corner_point_latlon(float(center_lat_full[0]), float(center_lon_full[0]),
                                                                               cls.e_we[0], cls.e_sn[0], 
                                                                               float(cls.dx[0]), float(cls.dy[0]),
                                                                               cls.wpsproj, cls.latlonproj, 'll')
         corner_lon_full[0,1], corner_lat_full[0,1] = wps_info.calc_corner_point_latlon(center_lat_full[0], center_lon_full[0],
                                                                               cls.e_we[0], cls.e_sn[0], 
                                                                               cls.dx[0], cls.dy[0],
                                                                               cls.wpsproj, cls.latlonproj, 'lr')
         corner_lon_full[0,2], corner_lat_full[0,2] = wps_info.calc_corner_point_latlon(center_lat_full[0], center_lon_full[0],
                                                                               cls.e_we[0], cls.e_sn[0], 
                                                                               cls.dx[0], cls.dy[0],
                                                                               cls.wpsproj, cls.latlonproj, 'ul')
         corner_lon_full[0,3], corner_lat_full[0,3] = wps_info.calc_corner_point_latlon(center_lat_full[0], center_lon_full[0],
                                                                               cls.e_we[0], cls.e_sn[0], 
                                                                               cls.dx[0], cls.dy[0],
                                                                               cls.wpsproj, cls.latlonproj, 'ur')
         # Getting information for subsequent domains
         if cls.max_dom>1:
             for i in np.arange(1,cls.max_dom):
                 
                 center_lon_full[i], center_lat_full[i] = wps_info.calc_center_point_latlon(corner_lat_full[i-1,0], corner_lon_full[i-1,0],
                                                                                   cls.dx[i-1], cls.dy[i-1],
                                                                                   cls.e_we[i], cls.e_sn[i],
                                                                                   cls.dx[i], cls.dy[i],
                                                                                   cls.i_parent_start[i], cls.j_parent_start[i],
                                                                                   cls.wpsproj, cls.latlonproj)
                 
                 
                 corner_lon_full[i,0], corner_lat_full[i,0] = wps_info.calc_corner_point_latlon(center_lat_full[i], center_lon_full[i],
                                                                                   cls.e_we[i], cls.e_sn[i], 
                                                                                   cls.dx[i], cls.dy[i],
                                                                                   cls.wpsproj, cls.latlonproj, 'll')
                 corner_lon_full[i,1], corner_lat_full[i,1] = wps_info.calc_corner_point_latlon(center_lat_full[i], center_lon_full[i],
                                                                                   cls.e_we[i], cls.e_sn[i], 
                                                                                   cls.dx[i], cls.dy[i],
                                                                                   cls.wpsproj, cls.latlonproj, 'lr')
                 corner_lon_full[i,2], corner_lat_full[i,2] = wps_info.calc_corner_point_latlon(center_lat_full[i], center_lon_full[i],
                                                                                   cls.e_we[i], cls.e_sn[i], 
                                                                                   cls.dx[i], cls.dy[i],
                                                                                   cls.wpsproj, cls.latlonproj, 'ul')
                 corner_lon_full[i,3], corner_lat_full[i,3] = wps_info.calc_corner_point_latlon(center_lat_full[i], center_lon_full[i],
                                                                                   cls.e_we[i], cls.e_sn[i], 
                                                                                   cls.dx[i], cls.dy[i],
                                                                                   cls.wpsproj, cls.latlonproj, 'ur')
                 
         #DO I REALLY WANT THESE VARIABLES PUBLIC, Is there a use that this class doesn't cover for?
         #Making private since you get these variables in the return
         cls._corner_lat_full =corner_lat_full
         cls._corner_lon_full = corner_lon_full
         cls._length_x = length_x
         cls._length_y = length_y

         return cls.wpsproj, cls.latlonproj, cls._corner_lat_full, cls._corner_lon_full, cls._length_x, cls._length_y
         
         
   

    #%%  FROM Xiaodong Chen for finding accurate map information on inner domains to plot       

    def calc_corner_point_latlon(center_lat, center_lon, e_we, e_ns, dx, dy, wpsproj, latlonproj, loc):
        """ FROM Xiaodong Chen @lucas-uw on github"""

        center_x, center_y = wpsproj.transform_point(center_lon, center_lat, latlonproj)
        if loc=='ll':
            xpt = center_x - dx*e_we/2.0
            ypt = center_y - dy*e_ns/2.0
        elif loc=='lr':
            xpt = center_x + dx*e_we/2.0
            ypt = center_y - dy*e_ns/2.0
        elif loc=='ul':
            xpt = center_x - dx*e_we/2.0
            ypt = center_y + dy*e_ns/2.0
        elif loc=='ur':
            xpt = center_x + dx*e_we/2.0
            ypt = center_y + dy*e_ns/2.0
        else:
            print("Invalid corner location (VALID: ll, lr, ul, ur)")
            raise ValueError

        corner_lon, corner_lat = latlonproj.transform_point(xpt, ypt, wpsproj)
        # TRANSFORMS TO DD Lat, DD LON
        
        return corner_lon, corner_lat
    
    def calc_center_point_latlon(corner_lat_parent, corner_lon_parent, dx_parent, dy_parent, e_we, e_ns, dx, dy, i, j, wpsproj, latlonproj):
        """ FROM Xiaodong Chen @lucas-uw on github"""
        corner_x_parent, corner_y_parent = wpsproj.transform_point(corner_lon_parent, corner_lat_parent, latlonproj)
        center_x_child = corner_x_parent + dx_parent*i + dx*e_we/2.0
        center_y_child = corner_y_parent + dy_parent*j + dy*e_ns/2.0
        center_lon_child, center_lat_child = latlonproj.transform_point(center_x_child, center_y_child, wpsproj) # FROM #### to Lat/Lon
        return center_lon_child, center_lat_child

        
    def reproject_corners(corner_lons, corner_lats, wpsproj, latlonproj):
         """ FROM Xiaodong Chen @lucas-uw on github"""

         corner_x = np.zeros((4,1))
         corner_y = np.zeros((4,1))
         corner_x[0], corner_y[0] = wpsproj.transform_point(corner_lons[0], corner_lats[0], latlonproj)
         corner_x[1], corner_y[1] = wpsproj.transform_point(corner_lons[1], corner_lats[1], latlonproj)
         corner_x[2], corner_y[2] = wpsproj.transform_point(corner_lons[2], corner_lats[2], latlonproj)
         corner_x[3], corner_y[3] = wpsproj.transform_point(corner_lons[3], corner_lats[3], latlonproj)
         return corner_x, corner_y
     

#%%   PLOT ZOOM(s) and PLOT_ALL_DOMAINS
    def plot_zoom_out(corner_x, corner_y, ns_zoom,we_zoom):
        """NOT FULLY TESTED.... FOR ALL PROJECTION TYPES... trial and error"""
        if (ns_zoom >= 0 and ns_zoom <40) and (we_zoom >= 0 and we_zoom <40):
            corner_x = corner_x + corner_x * we_zoom
            corner_y = corner_y + corner_y * ns_zoom
            return corner_x, corner_y 
        else:
            print ("Invalid Zoom Level, No Zoom Applied")
            return corner_x, corner_y 

    def plot_zoom_in(corner_x, corner_y, ns_zoom,we_zoom):
        """HAVE HAD NO ISSUES BUT NOT FULLY TESTED FOR ALL PROJECTION TYPES... trial and error"""
        if (ns_zoom >= 0 and ns_zoom <100) and (we_zoom >= 0 and we_zoom <100):
            corner_x = corner_x - corner_x * we_zoom/100.
            corner_y = corner_y - corner_y * ns_zoom/100.
            return corner_x, corner_y 
        else:
            print ("Invalid Zoom Level, No Zoom Applied")
            return corner_x, corner_y 

    def plot_all_domains(ax1, fig1,max_dom, wpsproj, latlonproj, corner_lat_full, corner_lon_full, length_x, length_y):
        """
        Plot all the Domains
        
        Location of dNum string to be dynamic relative to domain size and figure size would be much better.
        
        """
        # d01
        font0 = FontProperties()
        font0.set_weight('bold')
        #print("Number of Domains: ", max_dom)
        colors = ['blue',"white","red", "cyan", "magenta","gold","black","green","yellow",'pink']
        dNum = ['D01','D02','D03','D04','D05','D06','D07','D08','D09','D10']
        ## ORIGINAL VALUES
        #xbuff = [0.05, 0.05, 0.1]
        #ybuff = [0.9, 1.1, 0.8]
        for i in range(max_dom):
            if i == 0:
                corner_x, corner_y = wps_info.reproject_corners(corner_lon_full[i,:], corner_lat_full[i,:], wpsproj, latlonproj)
                ax1.set_xlim([corner_x[0]-length_x[0]/15, corner_x[3]+length_x[0]/15]) ## In Geodetic Coordinates # Want Lat/Lon so I can add gridlines and labels on axis
                ax1.set_ylim([corner_y[0]-length_y[0]/15, corner_y[3]+length_y[0]/15]) ## In Geodetic Coordinates # Want Lat/Lon so I can add gridlines and labels on axis
            elif i < 10:
                corner_x, corner_y = wps_info.reproject_corners(corner_lon_full[i,:], corner_lat_full[i,:], wpsproj, latlonproj)
            else:
                print ("Maximum Domains To Plot is 9")
                raise IndexError
                
            ax1.add_patch(mpl.patches.Rectangle((corner_x[0], corner_y[0]),  length_x[i], length_y[i], 
                                        fill=None, lw=3, edgecolor=colors[i], zorder=10))
            
            ax1.text(corner_x[0]+length_x[i]*0.05, corner_y[0]+length_y[i]*0.85, dNum[i],
                     fontproperties=font0, size=15, color=colors[i], zorder=10) ## Can be improved
        
        fig1.canvas.draw()
        
        ax1.add_feature(LAND, edgecolor='k',facecolor='limegreen')
        ax1.add_feature(OCEAN, edgecolor='k',facecolor='deepskyblue')#, facecolor='deepskyblue')
        ax1.add_feature(LAKES, edgecolor='k',facecolor='deepskyblue')
        ax1.coastlines('10m','black')
        ax1.add_feature(cartopy.feature.STATES.with_scale('10m'))
        
        # I DO NOT LIKE THIS LAMBERT PROJ WORK AROUND. WAIT UNTIL CARTOPY SUPPORTS IT OR DON'T USE OR IMPLEMENT 
        if "Lambert" in str(wpsproj):
            # Add tick marks for lambert projection (only)
            #wps_info.set_lambert_ticks(ax1)
            wps_info.set_lambert_ticks(ax1,xskip=5.,yskip=5.,x_thickness=14,y_thickness=14)
        elif "Mercator" in str(wpsproj) or "PlateCarree" in str(wpsproj):
            ax1.gridlines(color='lightgrey', linestyle='-', draw_labels=True)
        elif "Polar" in str(wpsproj) or "Stereo" in str(wpsproj) or "lat-lon" in str(wpsproj):
            pass
        else:
            pass
        
        # Don't think this is needed
        return ax1, fig1  
    
    
    def _ismissing(val, islat=True):
        """
        From WRF-PYTHON....
        Return True if a value is None or out of bounds.
        This function is used to check for invalid latitude/longitude values.
        Args
            val (numeric): A numeric value.
            islat (:obj:`bool`): Set to False if checking for longitude values.
            
        Returns
            :obj:`bool`: True if the value is None, or an out of bounds value.
        
        """
        if islat:
            if val is None:
                return True 
            
            if math.fabs(val) > 90.:
                return True
        else:
            if val is None:
                return True 
            
            if math.fabs(val) > 360.:
                return True
        
        return False
#%%   Xiaodong Chen - Projection Features For Lambert Confromal Plots
      
    # all these functions are necessary only when LCC projection is used.
    def find_side(ls, side):
         """
         FROM Xiaodong Chen @lucas-uw on github

         Given a shapely LineString which is assumed to be rectangular, return the
         line corresponding to a given side of the rectangle.
         
         """
         minx, miny, maxx, maxy = ls.bounds
         points = {'left': [(minx, miny), (minx, maxy)],
                   'right': [(maxx, miny), (maxx, maxy)],
                   'bottom': [(minx, miny), (maxx, miny)],
                   'top': [(minx, maxy), (maxx, maxy)],}
         return sgeom.LineString(points[side])
     
     
    def lambert_xticks(ax, ticks, size):
         """FROM Xiaodong Chen @lucas-uw on github - Draw ticks on the bottom x-axis of a Lambert Conformal projection."""
         te = lambda xy: xy[0]
         lc = lambda t, n, b: np.vstack((np.zeros(n) + t, np.linspace(b[2], b[3], n))).T
         xticks, xticklabels = wps_info._lambert_ticks(ax, ticks, 'bottom', lc, te)
         ax.xaxis.tick_bottom()
         ax.set_xticks(xticks)
         ax.set_xticklabels([ax.xaxis.get_major_formatter()(xtick) for xtick in xticklabels], size=size)
         
     
    def lambert_yticks_left(ax, ticks, size):
         """FROM Xiaodong Chen @lucas-uw on github - Draw tricks on the left y-axis of a Lamber Conformal projection."""
         te = lambda xy: xy[1]
         lc = lambda t, n, b: np.vstack((np.linspace(b[0], b[1], n), np.zeros(n) + t)).T
         yticks, yticklabels = wps_info._lambert_ticks(ax, ticks, 'left', lc, te)
         ax.yaxis.tick_left()
         ax.set_yticks(yticks)
         ax.set_yticklabels([ax.yaxis.get_major_formatter()(ytick) for ytick in yticklabels], size=size)
     
         
    def lambert_yticks_right(ax, ticks, size):
         """FROM Xiaodong Chen @lucas-uw on github - Draw ricks on the left y-axis of a Lamber Conformal projection."""
         te = lambda xy: xy[1]
         lc = lambda t, n, b: np.vstack((np.linspace(b[0], b[1], n), np.zeros(n) + t)).T
         yticks, yticklabels = wps_info._lambert_ticks(ax, ticks, 'right', lc, te)
         ax.yaxis.tick_right()
         ax.set_yticks(yticks)
         ax.set_yticklabels([ax.yaxis.get_major_formatter()(ytick) for ytick in yticklabels], size=size)
     
    def _lambert_ticks(ax, ticks, tick_location, line_constructor, tick_extractor):
         """FROM Xiaodong Chen @lucas-uw on github - Get the tick locations and labels for an axis of a Lambert Conformal projection."""
         outline_patch = sgeom.LineString(ax.outline_patch.get_path().vertices.tolist())
         axis = wps_info.find_side(outline_patch, tick_location)
         n_steps = 30
         extent = ax.get_extent(ccrs.PlateCarree())
         _ticks = []
         for t in ticks:
             xy = line_constructor(t, n_steps, extent)
             proj_xyz = ax.projection.transform_points(ccrs.Geodetic(), xy[:, 0], xy[:, 1])
             xyt = proj_xyz[..., :2]
             ls = sgeom.LineString(xyt.tolist())
             locs = axis.intersection(ls)
             if not locs:
                 tick = [None]
             else:
                 tick = tick_extractor(locs.xy)
             _ticks.append(tick[0])
         # Remove ticks that aren't visible:    
         ticklabels = copy(ticks)
         #ticklabels.tolist()
         while True:
             """Eric: HAD TO MAKE A CHANGE HERE BECAUSE OFTEN ticklabels is a numpy array and not a list.
             np.arrays cannot pop a value like that"""
             try:
                 index = _ticks.index(None)
             except ValueError:
                 break
             #print (type(_ticks),type(ticklabels))
             _ticks.pop(index)
             try:
                 list(ticklabels).pop(index)
             except:
                 ticklabels.pop(index)
         return _ticks, ticklabels
 
    #!! PROBLEM
    def set_lambert_ticks(ax,xskip=10.,yskip=10.,x_thickness=14,y_thickness=14):
        """
        PROBLEMS HERE.... extent_rounded doesn't match the real extent but if the extent is less than a skip increment then it becomes all messed up.
        labels are inconsistent. Not sure why it plots sometimes and not othertimes. I want to be able to have gridlines on set increments like at 5N,10N,15N,20N or 10N,20N,30N.
        Gridlines sometimes don't cover the entire area. 
        """
        extent = list(ax.get_extent(ccrs.PlateCarree()))
    
        # Incase the extent is less than the x/y skip value. But then that messes with the only working gridlines
        extent_rounded = [round(val/10.0) * 10 for val in extent]
        #print("Rounded",extent_rounded)
        if extent_rounded[0] == extent_rounded[1]:
            extent_rounded[0] -=10
            extent_rounded[0] +=10
        if extent_rounded[0] == extent_rounded[1]:
            extent_rounded[0] -=10
            extent_rounded[0] +=10
            
        #print("Extent",extent)
        #print("Rounded",extent_rounded)
        
        xticks = np.arange(extent[0],extent[1]+1,xskip) # added one to get the full extent
        yticks = np.arange(extent[2],extent[3]+1,yskip) # added one to get the full extent
        #xticks = np.arange(extent_rounded[0],extent_rounded[1]+1,xskip)
        #yticks = np.arange(extent_rounded[2],extent_rounded[3]+1,yskip)
        
        
        #print("x_ticks ",xticks,"\t\t","y_ticks ", yticks)
        ax.gridlines(xlocs=xticks, ylocs=yticks)
        #ax.gridlines(color='lightgrey', linestyle='-', draw_labels=True)

        # OPTION 1: Use Cartopy's LAT/LON FORMATTER
        #ax.xaxis.set_major_formatter(LONGITUDE_FORMATTER) 
        #ax.yaxis.set_major_formatter(LATITUDE_FORMATTER) 
        
        # OPTION 2: SET PRECISION
        ax.xaxis.set_major_formatter(ticker.FormatStrFormatter("%3.2f"))
        ax.yaxis.set_major_formatter(ticker.FormatStrFormatter("%3.2f")) 

        wps_info.lambert_xticks(ax, xticks, x_thickness)
        wps_info.lambert_yticks_left(ax, yticks, y_thickness)
        return ax