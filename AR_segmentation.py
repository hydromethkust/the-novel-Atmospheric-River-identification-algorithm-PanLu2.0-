#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 19 16:50:04 2020

@author: mengxinpan
"""

import numpy as np 
from geopy.distance import distance
from skimage.segmentation import felzenszwalb



# the input IVT_coherence_dir with 1 degree resolution
def AR_segmentation(IVT_coherence_idr,IVT_lat,IVT_long):
    

    ivt_coherence_max=np.nanmax(IVT_coherence_idr)
    ivt_coherence_min=np.nanmin(IVT_coherence_idr)
    ivt_coherence_image=((IVT_coherence_idr-ivt_coherence_min)*0.5/(ivt_coherence_max-ivt_coherence_min))
    ivt_coherence_image[np.where(np.isnan(ivt_coherence_image))]=1

    segments_felzenszwalb = felzenszwalb(ivt_coherence_image, scale=100, sigma=0, min_size=50)
    segments_felzenszwalb1=np.zeros(np.shape(IVT_coherence_idr))
    segments_felzenszwalb1[np.where(np.isnan(IVT_coherence_idr)==False)]=segments_felzenszwalb[np.where(np.isnan(IVT_coherence_idr)==False)]
    cluster_name=segments_felzenszwalb1[np.where(segments_felzenszwalb1!=0)]
    cluster_num=np.unique(cluster_name)
    segments_result=np.zeros(np.shape(IVT_coherence_idr))
    
    
    q50_list=[]
    q90_list=[]
    for i in range(len(cluster_num)):
        cluster=cluster_num[i]
        segments_result[np.where(segments_felzenszwalb1==cluster)]=i+1    
  
        segments_single=np.zeros(np.shape(IVT_coherence_idr))
        segments_single[:]=np.nan
        segments_single[np.where(segments_felzenszwalb1==cluster)]=IVT_coherence_idr[np.where(segments_felzenszwalb1==cluster)]
        segments_single_series=IVT_coherence_idr[np.where(segments_felzenszwalb1==cluster)]


        y_geo=np.array(IVT_lat[np.where(np.isnan(segments_single)==False)[0]])
        x_geo=np.array(IVT_long[np.where(np.isnan(segments_single)==False)[1]])  
        
        Area_total=0
        for l in range(0,len(x_geo)):
            long=x_geo[l]
            lat=y_geo[l]
            long_dist=distance((lat,long-0.5),(lat,long+0.5)).km

            if lat+0.5>90:
                continue
            elif lat-0.5<-90:
                continue
            else:
                lat_dist=distance((lat-0.5,long),(lat+0.5,long)).km
                area=long_dist*lat_dist
                Area_total=Area_total+area    

        q50=np.percentile(segments_single_series,50)
        q90=np.percentile(segments_single_series,90)
        q50_list.append(q50)
        q90_list.append(q90)
        
    return segments_result,q50_list,q90_list,Area_total
    