#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  10 22:53:42 2020

@author: PAN Mengxin
"""

import numpy as np 
import pandas as pd
from geopy.distance import distance
import math
import xarray as xr
from scipy.stats import norm
import datetime


# the calculation of the area of each 0.25*0.25 grid in the study domain
# it will be used in the AR pathway area calculation
def Area_matrix_calculation(IVT_lat,IVT_long):
    Area_global=np.zeros([len(IVT_lat),len(IVT_long)])
    for i in range(len(IVT_lat)):
        lat=IVT_lat[i]
        for j in range(len(IVT_long)):
            long=IVT_long[j]
            long_dist=distance((lat,long-0.125),(lat,long+0.125)).km
            if lat+0.125>90:
                lat_dist=distance((lat-0.125,long),(lat,long)).km
            elif lat-0.125<-90:
                lat_dist=distance((lat,long),(lat+0.125,long)).km
            else:
                lat_dist=distance((lat-0.125,long),(lat+0.125,long)).km 
            Area_global[i,j]=long_dist*lat_dist
    
    return Area_global



# quantile_field: the 85th quantile of IVT field for each month
# bandwidth: the bandwidth in the gaussian kernel smoothing
def dual_threshold_calculation_GKS(quantile_field,IVT_lat,IVT_long,bandwidth):
    local_threshold=np.zeros(np.shape(quantile_field))
    for target_lat_index in range(np.shape(quantile_field)[0]):
        for target_long_index in range(np.shape(quantile_field)[1]):
            target_lat=IVT_lat[target_lat_index]
            target_long=IVT_long[target_long_index]
            bandwidth_km=distance((target_lat,target_long),(target_lat,target_long+bandwidth)).km  
            gaussian_kernel=norm(loc=0, scale=bandwidth_km)
            kernel_matrix=np.zeros([151,151])
            ivt_matrix=np.zeros([151,151])
            for a in range(0,151):
                shift_lat_index=a-75
                if target_lat_index+shift_lat_index < 0 or target_lat_index+shift_lat_index >= len(IVT_lat):
                    continue
                neighbor_lat=IVT_lat[target_lat_index+shift_lat_index]
                for b in range(0,151):
                    shift_long_index=b-75
                    if target_long_index+shift_long_index < 0 or target_long_index+shift_long_index >= len(IVT_long):
                        continue
                    neighbor_long=IVT_long[target_long_index+shift_long_index]
                    dist= distance((target_lat,target_long),(neighbor_lat,neighbor_long)).km  
                    pdf=gaussian_kernel.pdf(dist)
                    kernel_matrix[a,b]=pdf
                    ivt_matrix[a,b]=quantile_field[target_lat_index+shift_lat_index,target_long_index+shift_long_index]             
            
            kernel_matrix_standard=kernel_matrix/np.sum(kernel_matrix)
            target_ivt_smooth=np.sum(kernel_matrix_standard*ivt_matrix)
            local_threshold[target_lat_index,target_long_index]=target_ivt_smooth
    return local_threshold




def AR_pathway_detection(ivt_total,ivt_minus,IVT_lat,IVT_long):
    add=np.array([0, 1, 0, -1, -1, 0, 1,0,1,1,1,-1,-1,1,-1,-1 ])
    add=add.reshape([8,2])
    ar_pathway=np.zeros(ivt_minus.shape)
    ar_pathway[:]=np.nan
    ar_pathwaynew=np.zeros(ivt_minus.shape)
    ar_pathwaynew[:]=np.nan
    
    lat_stop_small_index=min(np.where(IVT_lat==-80)[0][0],np.where(IVT_lat==80)[0][0])
    lat_stop_large_index=max(np.where(IVT_lat==-80)[0][0],np.where(IVT_lat==80)[0][0])

# detect the AR pathway crossing the EA box only (60N-0,40E-180)    
    index_lat_small=np.where(IVT_lat==60)[0][0]
    index_lat_large=np.where(IVT_lat==0)[0][0]
    index_long_small=np.where(IVT_long==40)[0][0]
    index_long_large=np.where(IVT_long==179.75)[0][0]
    ivt_minus_EA=np.zeros(np.shape(ivt_minus))
    ivt_minus_EA[:]=-999
    ivt_minus_EA[index_lat_small:index_lat_large,index_long_small:index_long_large+1]=ivt_minus[index_lat_small:index_lat_large,index_long_small:index_long_large+1]
    ivt_minus_max_EA=np.max(ivt_minus_EA)
    if ivt_minus_max_EA>0:
        exist=2
        index_max_lat=np.where(ivt_minus==ivt_minus_max_EA)[0][0]
        index_max_long=np.where(ivt_minus==ivt_minus_max_EA)[1][0]
        
        ar_pathway[index_max_lat,index_max_long]=ivt_minus_max_EA
        ar_pathwaynew[index_max_lat,index_max_long]=ivt_minus_max_EA   
        for loop in range(0,10000):
            new_point=np.where((np.isnan(ar_pathwaynew)!=True ))        
            if len (new_point[0])==0:
                break
            new_point_lat=new_point[0]
            new_point_long=new_point[1]                     
            ar_pathwaynew[:]=np.nan
            
            for l in range(0,len(new_point_lat)):
                lati=new_point_lat[l]
                longi=new_point_long[l]
                for a in range(0,8):
                    if lati+add[a,0] < lat_stop_small_index or lati+add[a,0] > lat_stop_large_index:
                        continue
                    elif longi+add[a,1]<0:
                        if ivt_minus[lati+add[a,0],longi+add[a,1]+np.shape(ivt_minus)[1]]>0 and str(ar_pathway[lati+add[a,0],longi+add[a,1]+np.shape(ivt_minus)[1]])=='nan':
                            ar_pathway[lati+add[a,0],longi+add[a,1]+np.shape(ivt_minus)[1]]=ivt_total[lati+add[a,0],longi+add[a,1]+np.shape(ivt_minus)[1]]
                            ar_pathwaynew[lati+add[a,0],longi+add[a,1]+np.shape(ivt_minus)[1]]=ivt_total[lati+add[a,0],longi+add[a,1]+np.shape(ivt_minus)[1]]                         
                    
                    elif longi+add[a,1] >= np.shape(ivt_minus)[1]:
                        if ivt_minus[lati+add[a,0],longi+add[a,1]-np.shape(ivt_minus)[1]]>0 and str(ar_pathway[lati+add[a,0],longi+add[a,1]-np.shape(ivt_minus)[1]])=='nan':
                            ar_pathway[lati+add[a,0],longi+add[a,1]-np.shape(ivt_minus)[1]]=ivt_total[lati+add[a,0],longi+add[a,1]-np.shape(ivt_minus)[1]]
                            ar_pathwaynew[lati+add[a,0],longi+add[a,1]-np.shape(ivt_minus)[1]]=ivt_total[lati+add[a,0],longi+add[a,1]-np.shape(ivt_minus)[1]]                                            
                    else:
                        if ivt_minus[lati+add[a,0],longi+add[a,1]]>0 and str(ar_pathway[lati+add[a,0],longi+add[a,1]])=='nan':
                            ar_pathway[lati+add[a,0],longi+add[a,1]]=ivt_total[lati+add[a,0],longi+add[a,1]]
                            ar_pathwaynew[lati+add[a,0],longi+add[a,1]]=ivt_total[lati+add[a,0],longi+add[a,1]] 
    else:
        exist=0
        ar_pathway[:]=np.nan              
    return exist,ar_pathway


# from centroid to find next reference grid by moving 100 km
def move100km(x_centroid,y_centroid,angle_centroid,direction):
    
    if direction == 'forward':
        degree_large=2
        degree_small=0.5
    elif direction == 'backward':
        degree_large=-2
        degree_small=-0.5
        
    degree_initial=(degree_large+degree_small)/2
    for loop in range(100):

        y_centroid1=y_centroid+degree_initial*math.sin(angle_centroid)
        
        distance_geo=distance((y_centroid,x_centroid), (y_centroid1, \
                           x_centroid+degree_initial*math.cos(angle_centroid))).km
        if abs(distance_geo-100)<0.01:
            break
        
        if distance_geo>100:
            degree_large=degree_initial
        else:
            degree_small=degree_initial
            
        degree_initial=(degree_large+degree_small)/2
    y_referpoint=y_centroid1
    x_referpoint=x_centroid+degree_initial*math.cos(angle_centroid)

    return x_referpoint,y_referpoint


# find the nearest neighbors based on the geo-distance
def find_nearest_neighbors(xstart_center,ystart_center,NNratio,ar_pathway,ivt_east,ivt_north,IVT_lat,IVT_long,resolution_long,resolution_lat):

    ar_num_grids=len(np.where(np.isnan(ar_pathway)==False)[0])
    NN_num=int(ar_num_grids*NNratio)
    stop=0

    add=np.array([0, 1, 0, -1, -1, 0, 1,0,1,1,1,-1,-1,1,-1,-1 ])
    add=add.reshape([8,2])    
    
    ar_pathway1=np.zeros(ar_pathway.shape)
    ar_pathway1[:]=np.nan
    
    ar_pathwaynew1=np.zeros(ar_pathway.shape)
    ar_pathwaynew1[:]=np.nan


    round_lon=np.round(xstart_center/resolution_long)*resolution_long
    if round_lon<-180:
        round_lon=round_lon+360
    elif round_lon>180-resolution_long:
        round_lon=round_lon-360
        
    startgrid_long_index=np.where(IVT_long==round_lon)[0][0]
    startgrid_lat_index=np.where(IVT_lat==np.round(ystart_center/resolution_lat)*resolution_lat)[0][0]
    ar_pathway1[startgrid_lat_index,startgrid_long_index]=1
    ar_pathwaynew1[startgrid_lat_index,startgrid_long_index]=1

    lat_NN=[]
    long_NN=[]
    ivt_intensity_NN=[]
    ivteast_NN=[]
    ivtnorth_NN=[]
    d_NN=[]
    for loop in range(0,1000):
        new_point=np.where((np.isnan(ar_pathwaynew1)!=True ))        
        if len (new_point[0])==0:
            break 
        new_point_lat_index=new_point[0]
        new_point_long_index=new_point[1]
        dist_list=[]
        for l in range(0,len(new_point_lat_index)):
            dist=distance((IVT_lat[new_point_lat_index[l]],IVT_long[new_point_long_index[l]]), (ystart_center,xstart_center)).km
            dist_list.append(dist)
            
        dist_list=np.array(dist_list)
        sort_index = np.argsort(dist_list)  
        for l in range(0,len(sort_index)):
            if len(lat_NN) >= NN_num:
                stop=1
                break
            lat_NN.append(IVT_lat[new_point_lat_index[sort_index[l]]])
            long_NN.append(IVT_long[new_point_long_index[sort_index[l]]])
            ivt_e=ivt_east[new_point_lat_index[sort_index[l]],new_point_long_index[sort_index[l]]]
            ivt_n=ivt_north[new_point_lat_index[sort_index[l]],new_point_long_index[sort_index[l]]]
            ivteast_NN.append(ivt_e)
            ivtnorth_NN.append(ivt_n)
            ivt_intensity_NN.append(np.sqrt(ivt_e*ivt_e+ivt_n*ivt_n))
            d_NN.append(dist_list[sort_index[l]])
            
        if stop == 1:
            lat_NN=np.array(lat_NN)
            long_NN=np.array(long_NN)
            ivt_intensity_NN=np.array(ivt_intensity_NN)
            ivteast_NN=np.array(ivteast_NN)
            ivtnorth_NN=np.array(ivtnorth_NN)
            d_NN=np.array(d_NN)
            break            

                     
        ar_pathwaynew1[:]=np.nan

        if loop==0:
            lati=new_point_lat_index[0]
            longi=new_point_long_index[0]
            for latadd in np.arange(-2,3):
                for longadd in np.arange(-2,3):
                    lati_add=lati+latadd
                    longi_add=longi+longadd
                    
                    if lati_add>=np.shape(ar_pathway)[0] or lati_add<0:
                        continue
                    
                    if longi_add<0:
                        longi_add=longi_add+np.shape(ar_pathway)[1]
                    elif longi_add >= np.shape(ar_pathway)[1]:
                        longi_add=longi_add-np.shape(ar_pathway)[1]
                    else:
                        longi_add=longi_add
                        
                    if str(ar_pathway[lati_add,longi_add])!='nan' and str(ar_pathway1[lati_add,longi_add])=='nan':
                        ar_pathway1[lati_add,longi_add]=ar_pathway[lati_add,longi_add]
                        ar_pathwaynew1[lati_add,longi_add]=ar_pathway[lati_add,longi_add]
                        
# if the refer point fall into a hole after moving 100km
            for loop_search in np.arange(2,100,5):
                
                new_point1=np.where((np.isnan(ar_pathwaynew1)!=True )) 
                
                if len (new_point1[0])==0:
                    lati=new_point_lat_index[0]
                    longi=new_point_long_index[0]
                    for latadd in np.arange(-loop_search,loop_search+1):
                        for longadd in np.arange(-loop_search,loop_search+1):
                            lati_add=lati+latadd
                            longi_add=longi+longadd
                            
                            if lati_add>=np.shape(ar_pathway)[0] or lati_add<0:
                                continue        
                            
                            if longi_add<0:
                                longi_add=longi_add+np.shape(ar_pathway)[1]
                            elif longi_add >= np.shape(ar_pathway)[1]:
                                longi_add=longi_add-np.shape(ar_pathway)[1]
                            else:
                                longi_add=longi_add
                                
                            if str(ar_pathway[lati_add,longi_add])!='nan' and str(ar_pathway1[lati_add,longi_add])=='nan':
                                ar_pathway1[lati_add,longi_add]=ar_pathway[lati_add,longi_add]
                                ar_pathwaynew1[lati_add,longi_add]=ar_pathway[lati_add,longi_add]
                else:
                    break
        else:         
                            
            for l in range(0,len(new_point_lat_index)):
                lati=new_point_lat_index[l]
                longi=new_point_long_index[l]
                
                for a in range(0,8):
                    lati_add=lati+add[a,0]
                    longi_add=longi+add[a,1]
    
                    if lati_add>=np.shape(ar_pathway)[0] or lati_add<0:
                        continue
                    
                    if longi_add<0:
                        longi_add=longi_add+np.shape(ar_pathway)[1]
                    elif longi_add >= np.shape(ar_pathway)[1]:
                        longi_add=longi_add-np.shape(ar_pathway)[1]
                    else:
                        longi_add=longi+add[a,1]
                        
                    if str(ar_pathway[lati_add,longi_add])!='nan' and str(ar_pathway1[lati_add,longi_add])=='nan':
                        ar_pathway1[lati_add,longi_add]=ar_pathway[lati_add,longi_add]
                        ar_pathwaynew1[lati_add,longi_add]=ar_pathway[lati_add,longi_add]                                               
                        
    return long_NN,lat_NN,ivt_intensity_NN,ivteast_NN,ivtnorth_NN,d_NN


# find the nearest neighbors for the IDR (IVT coherence) calculation
def find_nearest_neighbors500km(cen_long,cen_lat,ar_pathway,ivteast,ivtnorth,IVT_lat,IVT_long,resolution_long,resolution_lat):

    add=np.array([0, 1, 0, -1, -1, 0, 1,0,1,1,1,-1,-1,1,-1,-1 ])
    add=add.reshape([8,2])

    ar_pathway1=np.zeros(ar_pathway.shape)
    ar_pathway1[:]=np.nan
    
    ar_pathwaynew1=np.zeros(ar_pathway.shape)
    ar_pathwaynew1[:]=np.nan

    round_lon=np.round(cen_long/resolution_long)*resolution_long
    if round_lon<-180:
        round_lon=round_lon+360
    elif round_lon>180-resolution_long:
        round_lon=round_lon-360
        
    startgrid_long_index=np.where(IVT_long==round_lon)[0][0]
    startgrid_lat_index=np.where(IVT_lat==np.round(cen_lat/resolution_lat)*resolution_lat)[0][0]
    
    
    ar_pathway1[startgrid_lat_index,startgrid_long_index]=1
    ar_pathwaynew1[startgrid_lat_index,startgrid_long_index]=1

    lat_NN=[]
    long_NN=[]
    ivt_intensity_NN=[]
    ivteast_NN=[]
    ivtnorth_NN=[]
    dist_list=[]
    for loop in range(0,1000):
        new_point=np.where((np.isnan(ar_pathwaynew1)!=True ))        
        if len (new_point[0])==0:
            break 
        new_point_lat_index=new_point[0]
        new_point_long_index=new_point[1]
        dist_1round=[]
        for l in range(0,len(new_point_lat_index)):
            
            dist=distance((IVT_lat[new_point_lat_index[l]],IVT_long[new_point_long_index[l]]), (cen_lat,cen_long)).km
            dist_1round.append(dist)
            if dist < 500:
                dist_list.append(dist)
                lat_NN.append(IVT_lat[new_point_lat_index[l]])
                long_NN.append(IVT_long[new_point_long_index[l]])
                ivt_intensity_NN.append(ar_pathway[new_point_lat_index[l],new_point_long_index[l]])
                ivteast_NN.append(ivteast[new_point_lat_index[l],new_point_long_index[l]])
                ivtnorth_NN.append(ivtnorth[new_point_lat_index[l],new_point_long_index[l]])
            
        dist_1round=np.array(dist_1round)
        dist_1round_min=np.min(dist_1round)
        
        if dist_1round_min > 500:
            break       
               
        ar_pathwaynew1[:]=np.nan
        
        for l in range(0,len(new_point_lat_index)):
            lati=new_point_lat_index[l]
            longi=new_point_long_index[l]
            
            for a in range(0,8):
                lati_add=lati+add[a,0]
                longi_add=longi+add[a,1]
                if longi_add<0:
                    longi_add=longi_add+np.shape(ar_pathway)[1]
                elif longi_add >= np.shape(ar_pathway)[1]:
                    longi_add=longi_add-np.shape(ar_pathway)[1]
                else:
                    longi_add=longi+add[a,1]
                if lati_add>=np.shape(ar_pathway)[0] or lati_add<0:
                    continue                    
                if str(ar_pathway[lati_add,longi_add])!='nan' and str(ar_pathway1[lati_add,longi_add])=='nan':
                    ar_pathway1[lati_add,longi_add]=ar_pathway[lati_add,longi_add]
                    ar_pathwaynew1[lati_add,longi_add]=ar_pathway[lati_add,longi_add] 

    dist_list=np.array(dist_list)
    lat_NN=np.array(lat_NN)
    long_NN=np.array(long_NN)
    ivt_intensity_NN=np.array(ivt_intensity_NN)
    ivteast_NN=np.array(ivteast_NN)
    ivtnorth_NN=np.array(ivtnorth_NN)


    return long_NN,lat_NN,ivt_intensity_NN,ivteast_NN,ivtnorth_NN


def shift(xstart_center,ystart_center,long_NN,lat_NN, angle): # clockwise
    
    """
    Rotate a point counterclockwise by a given angle around a given origin.

    The angle should be given in radians.
    """
    angle = -angle
    ox, oy = xstart_center,ystart_center
    px, py = long_NN,lat_NN

    qx = math.cos(angle) * (px - ox) - math.sin(angle) * (py - oy)
    qy = math.sin(angle) * (px - ox) + math.cos(angle) * (py - oy)

    
    return qx, qy   
 
# rotate the coordinate back
def shift_back(xstart_center,ystart_center,y_pre, angle):
    
    """
    Rotate a point counterclockwise by a given angle around a given origin.

    The angle should be given in radians.
    """
    
    ox, oy = xstart_center,ystart_center
    px, py = 0,y_pre

    qx =ox + math.cos(angle) * (px) - math.sin(angle) * (py)
    qy =oy + math.sin(angle) * (px) + math.cos(angle) * (py)

    
    return qx, qy    
    
def dotproduct(v1, v2):
  return sum((a*b) for a, b in zip(v1, v2))

def length_vector(v):
  return math.sqrt(dotproduct(v, v))

   
def angle_vector(v1, v2):
    if length_vector(v1) * length_vector(v2)==0:
        angle=0
    else:
        
        cos_value=dotproduct(v1, v2) / (length_vector(v1) * length_vector(v2))
        if cos_value >= 1:
            cos_value = 0.999
        elif cos_value <= -1:
            cos_value = -0.999
        angle=math.acos(cos_value)
    return angle 


#  AR trajectory generation
def tracking(direction,NNratio,ar_pathway,ivt_east,ivt_north,IVT_lat,IVT_long,resolution_long,resolution_lat):
    
    
    x_trajectory=[]
    y_trajectory=[]

# start from the grid with maximum IVT intensity
    ivt_intensity_max=np.nanmax(ar_pathway)
    index_max=np.where(ar_pathway==ivt_intensity_max)
    x_geo_new=IVT_long[index_max[1][0]]
    y_geo_new=IVT_lat[index_max[0][0]]
    Turning_Angle=0
    for loop in range(0,500):

        if loop==0:
            xstart_center=x_geo_new
            ystart_center=y_geo_new
        else:
            if y_centroid >80 or y_centroid < -80:
                break
            
            xstart_center,ystart_center=move100km(x_centroid,y_centroid,angle_centroid,direction)
             
# find the nearest neighbors
        
        if xstart_center>180-resolution_long*0.5:
            xstart_center=xstart_center-360
        elif xstart_center<-180-resolution_long*0.5:
            xstart_center=xstart_center+360
        
        
        long_NN,lat_NN,ivt_intensity_NN,ivteast_NN,ivtnorth_NN,d_NN=find_nearest_neighbors(xstart_center,ystart_center,NNratio,ar_pathway,ivt_east,ivt_north,IVT_lat,IVT_long,resolution_long,resolution_lat)

# get the IVT transport direction
        d_gap=np.max(d_NN)-np.min(d_NN)
        if d_gap==0:
            w_d=np.ones(len(d_NN))
            w_d=0.5*w_d
        else:
            w_d=((np.max(d_NN)-d_NN)/d_gap)
        
        IVT_east_NN_mean=sum(w_d*ivteast_NN)/sum(w_d)
        IVT_north_NN_mean=sum(w_d*ivtnorth_NN)/sum(w_d)
        angle_refer=math.atan2(IVT_north_NN_mean,IVT_east_NN_mean)

#  shift coordinate (IVT transport direction as X-axis, cross-section as Y-axis)
        long_NNlist_new=np.zeros(len(long_NN))
        if xstart_center > 0:
            long_NNlist_new[long_NN<(xstart_center-180)]=long_NN[long_NN<(xstart_center-180)]+360
            long_NNlist_new[long_NN>=(xstart_center-180)]=long_NN[long_NN>=(xstart_center-180)]
        else:
            long_NNlist_new[long_NN > (xstart_center+180-resolution_long)]=long_NN[long_NN > (xstart_center+180-resolution_long)]-360
            long_NNlist_new[long_NN <= (xstart_center+180-resolution_long)]=long_NN[long_NN <= (xstart_center+180-resolution_long)]          
        
        x_rotate,y_rotate=shift(xstart_center,ystart_center,long_NNlist_new,lat_NN, angle_refer)
        if direction == 'forward':        
# if more than 90% of the point behind the cross-section, break forward tracking         
            count_minus=len(np.where(x_rotate<0)[0])
            if len(long_NN)*0.9 < count_minus:
                break
        elif direction == 'backward':
# if less than 10% of the point behind the cross-section, break backward tracking                     
            count_minus=len(np.where(x_rotate<0)[0])
            if len(long_NN)*0.1 > count_minus:
                break      
# get the new centroid on the cross-section
                
        ivt_gap=np.max(ivt_intensity_NN)-np.min(ivt_intensity_NN)
        if ivt_gap==0:
            w_ivt=np.ones(len(ivt_intensity_NN))
            w_ivt=0.5*w_ivt
        else:
            w_ivt=((np.max(ivt_intensity_NN)-ivt_intensity_NN)/ivt_gap)
        
        
        w_both=((w_ivt)**2+(w_d)**2)**0.5        
        y_pre=np.sum((y_rotate*w_both)/np.sum(w_both))           
        x_centroid,y_centroid=shift_back(xstart_center,ystart_center,y_pre, angle_refer)

        if x_centroid >= 179.75:
            x_centroid=x_centroid-360
        if x_centroid < -180:
            x_centroid=x_centroid+360
        if y_centroid>80 or y_centroid<-80:
            break

        long_NN,lat_NN,ivt_intensity_NN,ivteast_NN,ivtnorth_NN,d_NN=find_nearest_neighbors(x_centroid,y_centroid,NNratio,ar_pathway,ivt_east,ivt_north,IVT_lat,IVT_long,resolution_long,resolution_lat)

        
        ivt_gap=np.max(ivt_intensity_NN)-np.min(ivt_intensity_NN)
        if ivt_gap==0:
            w_ivt=np.ones(len(ivt_intensity_NN))
            w_ivt=0.5*w_ivt
        else:
            w_ivt=((np.max(ivt_intensity_NN)-ivt_intensity_NN)/ivt_gap)
        
        d_gap=np.max(d_NN)-np.min(d_NN)
        if d_gap==0:
            w_d=np.ones(len(d_NN))
            w_d=0.5*w_d
        else:
            w_d=((np.max(d_NN)-d_NN)/d_gap)   
        IVT_east_NN_mean=sum((w_d*ivteast_NN)/sum(w_d))
        IVT_north_NN_mean=sum((w_d*ivtnorth_NN)/sum(w_d))
        angle_centroid=math.atan2(IVT_north_NN_mean,IVT_east_NN_mean)                        
       

        if loop > 3:
            x_trajectory_array=np.array(x_trajectory)
            y_trajectory_array=np.array(y_trajectory)
            
            if xstart_center > 0:
                
                x_trajectorylist_new=np.zeros(len(x_trajectory_array))
                x_trajectorylist_new[x_trajectory_array<(xstart_center-180)]=x_trajectory_array[x_trajectory_array<(xstart_center-180)]+360
                x_trajectorylist_new[x_trajectory_array>=(xstart_center-180)]=x_trajectory_array[x_trajectory_array>=(xstart_center-180)]
            else:
                x_trajectorylist_new=np.zeros(len(x_trajectory_array))
                x_trajectorylist_new[x_trajectory_array > xstart_center+180-resolution_long]=x_trajectory_array[x_trajectory_array > xstart_center+180-resolution_long]-360
                x_trajectorylist_new[x_trajectory_array <= xstart_center+180-resolution_long]=x_trajectory_array[x_trajectory_array <= xstart_center+180-resolution_long]   

            
            angle_turn=angle_vector([x_trajectorylist_new[-1]-x_trajectorylist_new[-2],y_trajectory_array[-1]-y_trajectory_array[-2]], \
                                    [x_centroid-x_trajectorylist_new[-1],y_centroid-y_trajectory_array[-1]])
            direction_turn=np.cross([x_trajectorylist_new[-1]-x_trajectorylist_new[-2],y_trajectory_array[-1]-y_trajectory_array[-2]], \
                                    [x_centroid-x_trajectorylist_new[-1],y_centroid-y_trajectory_array[-1]])
            if direction_turn<0: # clockwise
                angle_turn=-angle_turn*180/math.pi
            else:
                angle_turn=angle_turn*180/math.pi
            
            Turning_Angle=np.nansum([Turning_Angle,angle_turn])
            if abs(Turning_Angle)>540:
                break
        x_trajectory.append(x_centroid)
        y_trajectory.append(y_centroid)        
    return x_trajectory,y_trajectory






##  main function for the PanLu2.0 algorithm
#
## definition of input variables:
## timestep: the target timestep
## dual_threshold: the dual threshold (including local and regional threshold) field for the month of target timestep
## ivt_east: the zonal component of IVT field
## ivt_north: the meridional component of IVT field
## IVT_long: the longitude coordinate of the IVT field
## ivt_total: the IVT intensity field
## IVT_lat: the latitude coordinate of the IVT field
## resolution_long: the longitude resolution of the IVT field
## resolution_lat: the latitude resolution of the IVT field

def PanLU_algorithm(timestep,dual_threshold,ivt_east,ivt_north,ivt_total,IVT_long,IVT_lat,resolution_long,resolution_lat,Area_global):    

    ivt_minus=ivt_total-dual_threshold
    new_initial=0
    parameters=["Length0.1","Width0.1","LWratio0.1","Turning_Angle0.1", \
                "Length0.2","Width0.2","LWratio0.2","Turning_Angle0.2", \
                "Length0.3","Width0.3","LWratio0.3","Turning_Angle0.3", \
                "Length0.4","Width0.4","LWratio0.4","Turning_Angle0.4", \
                "Length0.5","Width0.5","LWratio0.5","Turning_Angle0.5", \
                "key","timestep","Area_total","IVT_mean","IVT_mean_east","IVT_mean_north", \
                "IVT_sum","IVT_direction","percentage_tropics",'percentage_ivtdirection15']        
    
    
    AR_Para_total=pd.DataFrame(columns=parameters)   

    jj=0
    for j in range(0,200):
        exist,ar_pathway=AR_pathway_detection(ivt_total,ivt_minus,IVT_lat,IVT_long)
        if exist==0:
            break
        ar_pathway_index=np.where(np.isnan(ar_pathway)!=True)
# to remove the current AR pathway in the ivt_minus field to the next detection            
        ivt_minus[ar_pathway_index]=-999   
# to remove the AR pathway with the number of grids is smaller 800 to save the computational time, since these ARs will be eliminated by the length criteria eventually
        if len(ar_pathway_index[0])<800:
            continue 
        
        


# if percentage in tropics is > 0.95, eliminating the AR pathway directly to save the computational time
        pathway_numgrid=len(ar_pathway_index[0])
        lat_tropics_small_index=min(np.where(IVT_lat==-20)[0][0],np.where(IVT_lat==20)[0][0])
        lat_tropics_large_index=max(np.where(IVT_lat==-20)[0][0],np.where(IVT_lat==20)[0][0])
        pathway_tropics=ar_pathway[lat_tropics_small_index:lat_tropics_large_index,:]
        pathway_tropics_numgrid=len(pathway_tropics[np.where(np.isnan(pathway_tropics)==False)])
        percentage_tropics=pathway_tropics_numgrid/pathway_numgrid  
        if percentage_tropics>0.95:
            continue
        
# Area_total calculation 
        ivt_intensity=np.array(ar_pathway[ar_pathway_index])
        index_lat=ar_pathway_index[0]
        y_geo=IVT_lat[index_lat]
        index_long=ar_pathway_index[1]
        x_geo=IVT_long[index_long]  
        ivtnorth_ar=np.array(ivt_north[ar_pathway_index])
        ivteast_ar=np.array(ivt_east[ar_pathway_index])
        Area_total=np.sum(Area_global[ar_pathway_index])
        
# percentage IVT direction calculation           
        count_15=0
        for l in range(0,len(y_geo)):                
            ivt_direction_grid=math.atan2(ivtnorth_ar[l],ivteast_ar[l])*180/np.pi
            if abs(ivt_direction_grid) < 15 or abs(ivt_direction_grid)>180-15:
                count_15=count_15+1   
        percentage_ivtdirection15=count_15/len(ivteast_ar)    
        if percentage_tropics > 0.5 and percentage_ivtdirection15 > 0.5:
            continue
         
# trajectory generation
        NNratio_list=[0.3,0.1,0.2,0.4,0.5]     
        ar_traj5=np.zeros([5,1000,3])
        ar_traj5[:]=np.nan
        length_lwratio_Turning_Angle_break=0
        para=np.zeros(20)
        for NNratio_index in range(0,5):
            NNratio=NNratio_list[NNratio_index]
        #  forward tracking
            forward_tracking_x,forward_tracking_y=tracking('forward',NNratio,ar_pathway,ivt_east,ivt_north,IVT_lat,IVT_long,resolution_long,resolution_lat)        
        #  backward tracking  
            backward_tracking_x,backward_tracking_y=tracking('backward',NNratio,ar_pathway,ivt_east,ivt_north,IVT_lat,IVT_long,resolution_long,resolution_lat)            
        #   combine forward and backward trajectory
            backward_tracking_x=backward_tracking_x[1:]
            backward_tracking_x=backward_tracking_x[::-1]   
            backward_tracking_y=backward_tracking_y[1:]
            backward_tracking_y=backward_tracking_y[::-1]
            
            trajectory_x=np.array(backward_tracking_x+forward_tracking_x)
            trajectory_y=np.array(backward_tracking_y+forward_tracking_y)    
            
            if len(trajectory_x)<=2:
                Length=0
                Turning_Angle=0                
                Width=0
                LWratio=0
            else:
                Length=0
                Turning_Angle=0
                turning_angle_list=[]
                turning_angle_list.append(0)
                
                trajectory_xlist_new=np.zeros(len(trajectory_x))
                x_median=trajectory_x[int(len(trajectory_x)/2)]
                if  x_median > 0:
                     trajectory_xlist_new[ trajectory_x<( x_median-180)]= trajectory_x[ trajectory_x<( x_median-180)]+360
                     trajectory_xlist_new[ trajectory_x>=( x_median-180)]= trajectory_x[ trajectory_x>=( x_median-180)]
                else:
                     trajectory_xlist_new[ trajectory_x >  x_median+180-resolution_long]= trajectory_x[ trajectory_x >  x_median+180-resolution_long]-360
                     trajectory_xlist_new[ trajectory_x <=  x_median+180-resolution_long]= trajectory_x[ trajectory_x <=  x_median+180-resolution_long]              

                for xi in range(1,len( trajectory_xlist_new)):
                    length=distance((trajectory_y[xi], trajectory_xlist_new[xi]), (trajectory_y[xi-1], trajectory_xlist_new[xi-1])).km
                    Length=length+Length
                    
                    if xi<len( trajectory_xlist_new)-1:
                        angle_turn=angle_vector([ trajectory_xlist_new[xi]- trajectory_xlist_new[xi-1],trajectory_y[xi]-trajectory_y[xi-1]], \
                                                [ trajectory_xlist_new[xi+1]- trajectory_xlist_new[xi],trajectory_y[xi+1]-trajectory_y[xi]])
                        direction_turn=np.cross([ trajectory_xlist_new[xi]- trajectory_xlist_new[xi-1],trajectory_y[xi]-trajectory_y[xi-1]], \
                                                [ trajectory_xlist_new[xi+1]- trajectory_xlist_new[xi],trajectory_y[xi+1]-trajectory_y[xi]])
                        if direction_turn<0: # clockwise
                            angle_turn=-angle_turn*180/math.pi
                        else:
                            angle_turn=angle_turn*180/math.pi
                        turning_angle_list.append(angle_turn)
                        Turning_Angle=np.nansum([Turning_Angle,angle_turn])
                turning_angle_list.append(0)
                    
            if Length==0:
                Width=0
                LWratio=0

            else:
                Width=Area_total/Length
                LWratio=Length/Width
            
            if NNratio_index==0 and (Length<2000 or LWratio<2):
                length_lwratio_Turning_Angle_break=1
                break
            if Turning_Angle>360:
                length_lwratio_Turning_Angle_break=1
                break                    
            para[int(NNratio*10-1)*4:int(NNratio*10-1)*4+4]=np.array([Length,Width,LWratio,Turning_Angle])
            ar_traj5[int(NNratio*10-1),0:len(trajectory_y),0]=trajectory_y
            ar_traj5[int(NNratio*10-1),0:len(trajectory_x),1]=trajectory_x
            ar_traj5[int(NNratio*10-1),0:len(turning_angle_list),2]=turning_angle_list        
    

        if length_lwratio_Turning_Angle_break==1:
            continue
            
# para  calculation
        IVT_sum=np.nansum(ivt_intensity)
        IVT_mean=np.nanmean(ivt_intensity)
        IVT_mean_east=np.nanmean(ivteast_ar)
        IVT_mean_north=np.nanmean(ivtnorth_ar)                    
        IVT_direction=math.atan2(IVT_mean_north,IVT_mean_east)*180/np.pi                    

        key=str(timestep)+'_'+str(jj)
        jj=jj+1   
        
# IVT coherence IDR calculation and segmentation
        ivt_coherence_idr=np.zeros(np.shape(ar_pathway))
        ivt_coherence_idr[:]=np.nan
        for l in range(0,len(y_geo)):
            cen_long=x_geo[l]
            cen_lat=y_geo[l]
            long_NN,lat_NN,ivt_intensity_NN,ivteast_NN,ivtnorth_NN=find_nearest_neighbors500km(cen_long,cen_lat,ar_pathway,ivt_east,ivt_north,IVT_lat,IVT_long,resolution_long,resolution_lat)
            ivt_direction_list=[]
            for k in range(0,len(long_NN)):
                ivt_direction_list.append(math.atan2(ivtnorth_NN[k],ivteast_NN[k])*180/np.pi)
            ivt_direction_list=np.array(ivt_direction_list)
            ivt_direction_median=np.median(ivt_direction_list)
            ivt_direction_list_shifted=ivt_direction_list-ivt_direction_median
            ivt_direction_list_shifted2=np.zeros(np.shape(ivt_direction_list_shifted))
            for k in range(len(ivt_direction_list_shifted)):
                ivt_direction_shifted1=ivt_direction_list_shifted[k]
                if ivt_direction_shifted1>180:
                    ivt_direction_shifted2=ivt_direction_shifted1-360
                elif ivt_direction_shifted1<-180:
                    ivt_direction_shifted2=ivt_direction_shifted1+360
                else:
                    ivt_direction_shifted2=ivt_direction_shifted1
                ivt_direction_list_shifted2[k]=ivt_direction_shifted2
            
            q90=np.nanpercentile(ivt_direction_list_shifted2,90)
            q10=np.nanpercentile(ivt_direction_list_shifted2,10)
            idr=q90-q10
            ivt_coherence_idr[index_lat[l],index_long[l]]=idr
            
        IVT_coherence_idr_xr=xr.DataArray(ivt_coherence_idr,dims=['latitude','longitude'], \
                                    coords={'latitude':IVT_lat, 'longitude':IVT_long})  
        IVT_coherence_idr_xr=IVT_coherence_idr_xr.assign_coords(AR_key=key)        
               
    
# parameter summary
        para2_list=[key,timestep,Area_total,IVT_mean,IVT_mean_east,IVT_mean_north,IVT_sum,IVT_direction, \
                   percentage_tropics,percentage_ivtdirection15]
        para_array=np.concatenate([para,np.array(para2_list)],axis=0)
        para_pd=pd.DataFrame(para_array.reshape([1,len(parameters)]),columns=parameters,index=[key])  
        
        AR_pathway_xr=xr.DataArray(ar_pathway,dims=['latitude','longitude'], \
                                    coords={'latitude':IVT_lat, 'longitude':IVT_long})
        AR_traj_xr=xr.DataArray(ar_traj5,dims=['NNratio','point','lat_long'], \
                                    coords={'NNratio': np.arange(0.1,0.6,0.1), \
                                            'point': np.arange(0,1000), 'lat_long': ['lat', 'long','angle']})
        AR_pathway_xr=AR_pathway_xr.assign_coords(AR_key=key)
        AR_traj_xr=AR_traj_xr.assign_coords(AR_key=key)       
        AR_Para_total=AR_Para_total.append(para_pd,ignore_index=False)
        if new_initial==0:
            AR_pathway_list_total=AR_pathway_xr
            AR_traj_list_total=AR_traj_xr
            IVT_coherence_idr_total=IVT_coherence_idr_xr
            new_initial=1
        else:
            AR_traj_list_total=xr.concat([AR_traj_list_total,AR_traj_xr],dim='AR_key')
            AR_pathway_list_total=xr.concat([AR_pathway_list_total,AR_pathway_xr],dim='AR_key')   
            IVT_coherence_idr_total=xr.concat([IVT_coherence_idr_total,IVT_coherence_idr_xr],dim='AR_key')          
    
    return timestep,AR_traj_list_total,AR_pathway_list_total,AR_Para_total,IVT_coherence_idr_total




# After identifying all the AR slice, the final step is to generation AR sequence/life cycle
# AR_Para_total: the dateframe of all the AR parameter in the study period
# date_max: the last day of the study period
def ARS_generation(AR_Para_total,date_max):
    
    AR_Para_total['overlap_ratio']=np.nan
    AR_Para_total['ARS_key']=np.nan
    AR_Para_ARS=pd.DataFrame(columns=AR_Para_total.columns)
    AR_keys_total=AR_Para_total['key']
    ARS_dict={}
    
    for i in range(0,len(AR_keys_total)):
        key1=AR_keys_total.iloc[i]
        if str(AR_Para_total.loc[key1,'ARS_key']) != 'nan':
            continue
    
        datetime1=pd.to_datetime(key1[0:19])
        AR_pathway_list=xr.open_dataset("path_of_AR_pathway_ncfile_in_the_datetime1")['AR_pathway']    
        pathway1=AR_pathway_list.loc[dict(AR_key=key1)].values
    
        keylist=[]
        keylist.append(key1)        
        overlap_list=[]
# in an ARS, the AR slice with overlap ratio = 2 indicate it is the genesis of AR sequence/event 
        overlap_list.append(2)    
        for j in range(0,1000):
            datetime2=datetime1+datetime.timedelta(hours=6)

# if out of the time range, save the AR para and break
            if datetime2>date_max:
                if len(keylist)>3:
                    ARS_dict[keylist[0]]=keylist
                    print(keylist[0])
                    for k in range(0,len(keylist)):
                        AR_Para_total.loc[keylist[k],'overlap_ratio']=overlap_list[k]
                        AR_Para_total.loc[keylist[k],'ARS_key']=keylist[0]            
                        AR_Para_ARS=AR_Para_ARS.append(AR_Para_total.loc[keylist[k],:],sort=False)                        
                break
            AR_pathway_list=xr.open_dataset("path_of_AR_pathway_ncfile_in_the_datetime2")['AR_pathway']    
            index2=np.where(AR_Para_total['timestep']==str(datetime2))[0]
            ar_para2=AR_Para_total.iloc[index2,:]
            
            overlap_rate_list=[]
            AR_key_list2=ar_para2['key']
            for p in range(0,len(ar_para2)):
                
                key2=AR_key_list2[p]
                
                if str(ar_para2.loc[key2,'ARS_key'])!='nan':                
                    overlap_rate_list.append(0)
                    continue
                pathway2=AR_pathway_list.loc[dict(AR_key=key2)].values
                pathway_sum=pathway1+pathway2
                num_sum=np.count_nonzero(~np.isnan(pathway_sum))
                num_1=np.count_nonzero(~np.isnan(pathway1))
                overlap_rate=num_sum/num_1
                overlap_rate_list.append(overlap_rate)
            if len(overlap_rate_list)==0:
                overlap_rate_max=0
            else:
                overlap_rate_max=max(np.array(overlap_rate_list))
                
            if overlap_rate_max>0.5: # pathway1 and pathway2 have overlap
                overlap_max_index=np.where(np.array(overlap_rate_list)==overlap_rate_max)[0][0]
                keylist.append(str(AR_key_list2[overlap_max_index]))
                overlap_list.append(overlap_rate_max)
                
                key1=str(AR_key_list2[overlap_max_index])
                datetime1=datetime2
                pathway1=AR_pathway_list.loc[dict(AR_key=key1)].values
                
            else:                
                datetime2=datetime1+datetime.timedelta(hours=12)
                if datetime2>date_max:
                    if len(keylist)>3:
                        ARS_dict[keylist[0]]=keylist
                        print(keylist[0])
                        for k in range(0,len(keylist)):
                            AR_Para_total.loc[keylist[k],'overlap_ratio']=overlap_list[k]
                            AR_Para_total.loc[keylist[k],'ARS_key']=keylist[0]            
                            AR_Para_ARS=AR_Para_ARS.append(AR_Para_total.loc[keylist[k],:],sort=False)                        
                    break
                
                AR_pathway_list=xr.open_dataset("path_of_AR_pathway_ncfile_in_the_datetime2")['AR_pathway']    
                index2=np.where(AR_Para_total['timestep']==str(datetime2))[0]
                ar_para2=AR_Para_total.iloc[index2,:]
                    
                overlap_rate_list=[]
                AR_key_list2=ar_para2['key']
                for p in range(0,len(ar_para2)):
                    
                    key2=AR_key_list2[p]
                    
                    if str(ar_para2.loc[key2,'ARS_key'])!='nan':                
                        overlap_rate_list.append(0)
                        continue
                    pathway2=AR_pathway_list.loc[dict(AR_key=key2)]
                    pathway_sum=pathway1+pathway2
                    num_sum=np.count_nonzero(~np.isnan(pathway_sum))
                    num_1=np.count_nonzero(~np.isnan(pathway1))
                    overlap_rate=num_sum/num_1
                    overlap_rate_list.append(overlap_rate)        
                if len(overlap_rate_list)==0:
                    overlap_rate_max=0
                else:
                    overlap_rate_max=max(np.array(overlap_rate_list))
                
                if overlap_rate_max>0.5: # pathway1 and pathway2 have overlap
                    overlap_max_index=np.where(np.array(overlap_rate_list)==overlap_rate_max)[0][0]
                    keylist.append(str(AR_key_list2[overlap_max_index]))
                    overlap_list.append(overlap_rate_max)
                    key1=str(AR_key_list2[overlap_max_index])
                    datetime1=datetime2
                    pathway1=AR_pathway_list.loc[dict(AR_key=key1)]
                    
                else:
                    if len(keylist)>3:
                        ARS_dict[keylist[0]]=keylist
                        print(keylist[0])
                        for k in range(0,len(keylist)):
                            AR_Para_total.loc[keylist[k],'overlap_ratio']=overlap_list[k]
                            AR_Para_total.loc[keylist[k],'ARS_key']=keylist[0]            
                            AR_Para_ARS=AR_Para_ARS.append(AR_Para_total.loc[keylist[k],:],sort=False)
                    break
    AR_pathway_list.close()      
    return AR_Para_ARS,ARS_dict