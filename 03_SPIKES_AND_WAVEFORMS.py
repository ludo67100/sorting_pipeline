# -*- coding: utf-8 -*-
"""
Created on Wed Apr  3 18:11:04 2019

This script allows to extract spike times and waveform from tridesclous catalogue
in excel sheet
+ figure plot 

This script should run in a tridesclous environement where TDC is installed 

Take a look at the DATA INFO section, fill the info, run and let the magic happen. 

@author: ludovic
"""

#------------------------------------------------------------------------------
#-----------------------------DATA INFO----------------------------------------
#----------------------------FILL BELOW----------------------------------------

#The path of the TDC catalogue file - must be STRING format
path ='C:/Users/ludov/Documents/Fede/2535-1300-P3/tdc_concatenate/'

#Name of the experiment, protocol... anything to discriminate the date. Will be used
#for datasheets/figure labeling - must be STRING format  
name = '2535-1300-P3'

#Where to save datasheets and figures. If None, nothing will be saved  - must be STRING format
savedir = 'C:/Users/ludov/Documents/Fede/2535-1300-P3'

sampling_rate = 20000 #in Hz

#Stim time ad stim duration in seconds
stim_time = 4.0
stim_duration = 1.0

#Number of files concatenated to built the unique file 
conc = 10

#Specify the channel group to explore as [#]. Feel free to do them all : [0,1,2,3]
channel_groups=[0,1]

#If True : close figure automatically (avoids to overhelm screen when looping and debug)
closefig = False

#------------------------------------------------------------------------------
#-----------------------------THE SCRIPT---------------------------------------
#---------------------------DO NOT MODIFY--------------------------------------
import tridesclous as tdc 
import numpy as np 
from matplotlib import pyplot as plt 
import pandas as pd 

#Load the catalogue
dataio = tdc.DataIO(path)

for chan_grp in channel_groups:
    #Define the constructor and the channel group 
    cc = tdc.CatalogueConstructor(dataio, chan_grp=chan_grp)
    print ('--- Experiment : {} ---'.format(name))
    print ('Catalogue loaded from {}'.format(path))
    print ('----Channel group {}----'.format(chan_grp))
    
    #The cluster list for the spikes 
    clust_id = cc.all_peaks['cluster_label']
    
    #The spike times
    sampling_period = 1.0/sampling_rate 
    spike_times = cc.all_peaks['index']*sampling_period
    
    #Time vector of the whole trace
    len_trace = cc.info['processed_length']
    time_vector = np.arange(0,len_trace,1)*sampling_period
    
    #Stim vector for the whole trace
    stim_vector = np.arange(stim_time,time_vector[-1],float(conc))
    
    #The cluster label for median waveform
    waveforms_label = cc.clusters['cluster_label']
    
    #The median waveforms 
    waveforms = cc.centroids_median
    
    #The probe geometry and specs 
    probe_geometry = cc.geometry
    probe_channels = cc.nb_channel
    
    #Figure for waveforms------------------------------------------------------
    fig1 = plt.figure(figsize=(4,6))
    plt.title('{} Average Waveforms (ch_group = {})'.format(name,chan_grp))
    plt.xlabel('Probe location (micrometers)')
    plt.ylabel('Probe location (micrometers)')
          
    for cluster in np.unique(clust_id):
        for loc, prob_loc in zip(range(len(probe_geometry)), probe_geometry): 
            x_offset, y_offset = prob_loc[0], prob_loc[1]
            #base_x = np.arange(0,len(waveforms[1,:,loc]),1)  
            base_x = np.linspace(-15,15,num=len(waveforms[1,:,loc])) #Basic x-array for plot, centered
         
            if cluster == -1:# First array of waveforms, contains the zero line 
                
                if y_offset!=0: #Top and down zero line
                    plt.plot(base_x,waveforms[0,:,loc]+y_offset, color='0.8')
                
                if x_offset!=0: #Left and right zeeo_line
                    plt.plot(base_x+2*x_offset,waveforms[0,:,loc]+y_offset,color='0.8')
                    
            else:
                                
                clust_color = 'C{}'.format(cluster+1)
    
                if y_offset!=0: #Top and down probe 
                    if loc == 0 : #to avoid fucking legend redundancy
                       median = plt.plot(base_x,waveforms[cluster+1,:,loc]+y_offset,color=clust_color,label='Cluster {}'.format(cluster+1))
                    else :
                       median = plt.plot(base_x,waveforms[cluster+1,:,loc]+y_offset,color=clust_color)
    
                plt.legend()
                   
                if x_offset!=0: #Left and right probe
                    plt.plot(base_x+2*x_offset,waveforms[cluster+1,:,loc]+y_offset,color=clust_color)
    
    if savedir !=None :
        fig1.savefig('{}/{}_Waveforms_changrp_{}.png'.format(savedir,name,chan_grp))
        with pd.ExcelWriter('{}/{}_waveforms_changrp_{}.xlsx'.format(savedir,name,chan_grp)) as writer:
            for cluster in np.unique(clust_id):
                if cluster==-1:
                    clust_WF = pd.DataFrame(waveforms[0,:,:])      
                    clust_WF.to_excel(writer,sheet_name='horizon')               
                    
                else:
                    clust_WF = pd.DataFrame(waveforms[cluster+1,:,:])      
                    clust_WF.to_excel(writer,sheet_name='cluster {}'.format(cluster))
                
    else : 
        print ('No savedir specified : nothing will be saved')

    
    if closefig==True:
        plt.close()
    
    #Spike Times extraction per cluster---------------------------------------- 
    fig2, ax =plt.subplots(2,1,figsize=(10,5))
    ax[0].set_title('{} All spike times (ch_group = {})'.format(name,chan_grp))
    ax[0].eventplot(spike_times)
    ax[1].set_xlabel('Time (s)')
    ax[1].set_ylabel('Cluster ID')
    
    for stim in stim_vector:
        ax[0].axvspan(stim,stim+stim_duration,color='skyblue',alpha=0.6)
        ax[1].axvspan(stim,stim+stim_duration,color='skyblue',alpha=0.6)
                    
    SPIKES = [] #To store all the spikes, one array per cluster 
    
    for cluster in np.unique(clust_id):
        if cluster==-1:
            continue
        
        clust_color = 'C{}'.format(cluster+1)
    
        temp_ = [] #To store spikes from each cluster
    
        for i,j in np.ndenumerate(clust_id):
            if j == cluster:
                temp_.append(spike_times[i])
                
        SPIKES.append(np.asarray(np.ravel(temp_)))
        
        ax[1].eventplot(np.ravel(temp_), lineoffsets=cluster+1, linelengths=0.5, color=clust_color)
       
    if closefig==True:
        plt.close()
        
    #SAVE THE SPIKE DATA (or not) ---------------------------------------------
    
    if savedir != None:
        sorted_spikes = pd.DataFrame(SPIKES)
        sorted_spikes.to_excel('{}/{}_Spike_times_changrp_{}.xlsx'.format(savedir,name,chan_grp))
        fig2.savefig('{}/{}_Spike_times_changrp_{}.png'.format(savedir,name,chan_grp))
        
    
        
    