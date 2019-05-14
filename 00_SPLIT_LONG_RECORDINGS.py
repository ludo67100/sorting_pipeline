# -*- coding: utf-8 -*-
"""
Created on Tue May 14 10:48:31 2019

Slices long recordings in the willed amount of epsisodes 

@author: ludovic.spaeth
"""

import os,sys
import numpy as np 

#------------------FILES info---------------------------------------------------

nb_channel = 16
sample_rate = 20000 #in Hz
t_start = 0 

split_in = 20 #Episodes 

#Folder for data 
folder  = '//equipe2-nas1/Public/Federica/Ephy/5101 (Baseline of 2s - Atlas - Male)/HDF5/rbf/P13'

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

sampling_period = 1./sample_rate

#All the files in folder
list_dir = os.listdir(folder)
#Only the .rbf files 
files = [file for file in list_dir if '.rbf' in file]

for file,idx in zip(files, range(len(files))):
    
    #To avoid last limit
    split_in = split_in+1
    
    #The file name
    new_path = '{}/{}'.format(folder, file)
    
    #Load signals in array
    sigs = np.fromfile(new_path, dtype='float64').reshape(nb_channel,-1)

    if sigs.shape[1]*sampling_period < 181: #Recording of 180s

        #Delimits the length and begining of each episodes 
        step = int(sigs.shape[1]/20)
        starts = np.arange(0,sigs.shape[1],step)
        
        #Remove the last step 
        starts = starts[:-1]
        
    else: 
        sys.exit('This file is over 180s, you should consider changing the parameters')
        
    for i in range(len(starts)):
              
        segment_start = int(starts[i])
        segment_stop = int(segment_start+step-1)
        
        print ('Saving segment from {}s to {}s'.format(segment_start*sampling_period,segment_stop*sampling_period))
        
        file_save = '{}/EP{}_{}'.format(folder,i+1,file)
        
        episode = sigs[:,segment_start:segment_stop]
        
        with open(file_save, mode='wb') as savefile: 
            
            episode.tofile(file_save,sep='')
    
