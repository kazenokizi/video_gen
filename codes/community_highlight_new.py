#!/usr/bin/env python
# -*- coding: utf-8 -*-
# community_highlight_new.py
#
# This code allows to color the PDB structure according to the 
# communities identified by Markov Stability. 
#
# The coloring is designed such that similar communities identified
# at different Markov times are assigned the same color.
# The color matching is not perfect and the results are thus not a 
# perfect representation of which communities are kept through
# different Markov times. It is mostly designed as a visual aid and
# the results should be interpreted accordingly.
#
# Author: Antoine Delmotte
# (Adapted from the original code by Yun William Yu)
#
#######################################################################
#
# Command description:
#
# From PyMol, or a pml script which can be called as 
# pymol -c pml_script.pml
# The code can be run via the following commands in PyMol: 
#
# run community_highlight_new.py
# part2png('Stability.mat','file_mask', start_MarkovTime, end_MarkovTime, rotating_angle, image_per_MarkovTime)
#
# Stability.mat is the output file from stability.m. It should contain 
# 4 variables: T, N, VI and C. Refer to the documentation of the Markov 
# Stability Matlab package for more information.
#
# Rotating angle is useful if the user wants to see the molecule rotating 
# by an angle rotating_angle between each image.
#
# image_per_MarkovTime allows, when combined with rotating_angle to rotate
# the molecule for different Markov times.
# 
# All the inputs are optional except for the first one.

import scipy
from scipy import *
from scipy import io
import os

def part2png(mat_file,filemask='output',start=0,end=5,angle_rot=1.5,img_per_MT=1,sig=1):
  
  # Create the directory where PNG files will be stored
  if (not os.path.isdir('./'+filemask+'_png')):
    os.mkdir('./'+filemask+'_png')
    
  # Extract data from .mat file
  data = scipy.io.loadmat(mat_file)
  Time=data['Time'];
  Communities=data['C'];
  NC=data['N'];
  
  # CD to PNG directory
  directory_old = os.getcwd()
  os.chdir('./'+filemask+'_png')
  
  # Initialisation
  zero_number='00000'
  Communities=Communities[:,nonzero((start<Time.T)&(Time.T<end))[0]]
  NC=NC[(Time>start)&(Time<end)] 
  maxi = NC.max();
  Time=Time[(Time>start)&(Time<end)] 
  comm_prec=zeros((len(Communities[:,0]),1)) #comm_prec=Communities[:,0]
  
  # Span all the partitions between times start and end    #for i in range(Communities.shape[1]):
  for i in arange(Communities.shape[1]-1,-1,-1):   
    comm_curr=renumber_comm(Communities[:,i],comm_prec,NC[i],maxi)
    community_highlight2(comm_curr, maxi)
    comm_prec=comm_curr
    for p in range(img_per_MT):
      prefix=zero_number[0:len(zero_number)-len(str(int(i*img_per_MT+p)))]+str(int(i*img_per_MT+p))
      comm_to_png(prefix+'_'+filemask+'_'+str(Time[i]),sig)
      cmd.rotate('y',angle_rot)
  os.chdir(directory_old)
      
       
      
def renumber_comm(comm_curr,comm_prec,NC,maxi):
  # Initialisation 
  comm_curr_new=zeros((comm_curr.shape[0],1)) 	# Output matrix
  labels_match=zeros((NC,3))			# Labels correspondence between original partition and relabeled partition	
  labels_match[:,0]=unique(comm_curr)
  labels_prec=unique(comm_prec);
  labels_taken=maxi;			# Labels that are already taken by a community
  count=0;
  
  # Find the best community label for each community
  for i in range(NC):
    # Extract the labels in the previous partition that correspond to the current community in the new partition
    assign=comm_prec[comm_curr==labels_match[i][0]]
    # For each community, store the best label to which it correspond in the previous partition
    for j in unique(assign):
      if ((assign==j).sum()>labels_match[i][2]):
	labels_match[i][1]=j
	labels_match[i][2]=(assign==j).sum()

  # Give the label only to the community that have the largest number of node in common with the community of the previous partition
  for i in range(NC):
    # If the label correspondance for community i is the best match, and if this label has not already been taken (to avoid equally good communities)
    if ((labels_match[i][2] == max(labels_match[labels_match[:,1]==labels_match[i][1],2])) and (not (labels_taken==labels_match[i][1]).any())):
      labels_taken=scipy.append(labels_taken,labels_match[i][1])  
    else:
      max_dist_label=0
      labels_taken_temp=scipy.append(labels_taken,unique(comm_prec))
      labels_taken_temp=scipy.append(labels_taken_temp,maxi)
      labels_taken_temp=scipy.append(labels_taken_temp,0)
      labels_taken_temp=sort(labels_taken_temp)
      old_label=labels_match[i][1]
      for k in arange(len(labels_taken_temp)-2,-1,-1):
	if (((labels_taken_temp[k+1]-labels_taken_temp[k])>max_dist_label) and comm_adj(comm_curr,labels_match,i,labels_taken_temp[k]+(labels_taken_temp[k+1]-labels_taken_temp[k])/2,maxi)):
	  max_dist_label=labels_taken_temp[k+1]-labels_taken_temp[k]
	  labels_match[i][1]=labels_taken_temp[k]+max_dist_label/2
      labels_taken=scipy.append(labels_match[i][1],labels_taken)
      
            
  # If it happens that two communities have the same color...
  for i in unique(labels_match[:,1]):
    labels_match_temp=labels_match[labels_match[:,1]==i,:];
    for j in range(len(labels_match_temp[:,0])-1):
      labels_taken_temp=unique(sorted(labels_taken))
      idx=random.randint(0,len(labels_taken_temp)-1) 
      labels_match[labels_match[:,0]==labels_match_temp[j,0]]=labels_taken_temp[idx]+(labels_taken_temp[idx+1]-labels_taken_temp[idx])/2
      labels_taken=scipy.append(labels_match[labels_match[:,0]==labels_match_temp[j,0]],labels_taken)
      
      
	    
  for i in unique(comm_curr):
    comm_curr_new[comm_curr==i]=labels_match[labels_match[:,0]==i,1]
    
  return comm_curr_new



def comm_to_png(filemask,sig):
#	cmd.set('ray_trace_mode', '3')
#  cmd.set('ray_trace_gain', '0.05')
#  cmd.hide('everything','all')
#  cmd.show('cartoon')
#  cmd.show('spheres','name MG')
  # bazinga
	# modified by Heng Zhang, Imperial College London, Oct 7 2015
	if sig == 2:
		'''
		these commands for producing all with cartoon and sticks
		'''
		cmd.show('sticks')
	elif sig == 3:
		'''
  		these commands for producing selected amino acids 
  		'''
		cmd.show('sticks')
		# choose residue group for showcase
		cmd.select('resigroup','resi 10+30+50+70+90+110+130+150+170+190+210')
		cmd.set('stick_transparency',0.7,'! resigroup')
		cmd.set_bond('stick_transparency',0.7,'! resigroup')
		cmd.set('cartoon_transparency',0.7)
	elif sig == 4:
  		'''
  		these commands for producing selected alpha helices
  		'''
		cmd.show('sticks')
		cmd.set('stick_transparency',0.7)
		cmd.set('cartoon_transparency',0.7)
		# choose residue group for showcase
		cmd.create('resigroup','resi 13-24')
		cmd.set('cartoon_transparency',0,'resigroup')
		# set view should be modified according to pml file  		
		cmd.set_view ('\
     		0.417308629,   -0.071697339,   -0.905933022,\
     0.437800139,    0.889436722,    0.131275997,\
     0.796357572,   -0.451400787,    0.402558148,\
    -0.000025660,    0.000025567, -198.414794922,\
    -2.271391630,   -3.156492949,  -16.358736038,\
   159.843780518,  236.987731934,  -20.000000000')
  
#  cmd.set('sphere_scale','2')
	cmd.set('antialias','1')
#  cmd.set('depth_cue','0')
	cmd.set('ray_trace_fog',"1")
	cmd.set("ray_trace_frames","1")
	cmd.set('cartoon_fancy_helices','1')
	cmd.bg_color(color='white')
#	cmd.ray(3300,3300)
	cmd.ray(600,600) # bazinga 
	cmd.png(filemask+'.png')
  

from pymol import cmd 
from colorsys import hsv_to_rgb

def community_highlight2(assignments, maxi):
  selections = cmd.get_names('selections')
  
  # Remove previous selections
  for i in range(len(selections)):
    cmd.delete(str( selections[i]))
  
  # Create the new selections
  for i in unique(assignments):
    cmd.select(i, 'none')

  # Put the correct nodes in the selection bins
  b = []
  prev_comm = assignments[0]
  length = 1
  for i in range(1,len(assignments)):
    if assignments[i] == prev_comm:
      length = length + 1
    else:
      b.append([i,length,prev_comm[0]])
      length = 1
    prev_comm = assignments[i]
  b.append([i,length,prev_comm[0]])
  
  for i in range(len(b)):
    print str(b[i][2])+'|id '+str(b[i][0]-b[i][1]+1)+'-'+str(b[i][0])
    cmd.select(b[i][2],str(b[i][2])+'|id '+str(b[i][0]-b[i][1]+1)+'-'+str(b[i][0]))
 # cmd.space('cmyk')
  for i in unique(assignments):
    cmd.set_color("ccc"+str(i),hsv_to_rgb(((i/maxi)+0.16)%1,0.8+((round(64.0*i/maxi)%2)+1)*0.1+((round(16.0*i/maxi)%2)+1)*0.1,1-0.35*(((i/maxi)+0.16)%1)-0.2*((round(32*i/maxi)%2)))) #+0.0278    +0.16    + 0.0972    + 0.16666666
    cmd.color("ccc"+str(i),i)
  
  
  
def comm_adj(curr_comm,label_match,comm,label,maxi):
  adj=zeros((curr_comm.max()+1,curr_comm.max()+1))
  ok=True
  for i in range(len(curr_comm)-1):
    if curr_comm[i] != curr_comm[i+1]:
      adj[curr_comm[i],curr_comm[i+1]]=1
      adj[curr_comm[i+1],curr_comm[i]]=1
  for i in range(len(label_match)):
    if ((adj[label_match[comm,0],label_match[i,0]]==1) and ((0.1*maxi+label_match[i,2]>label) and (label_match[i,2]-0.1*maxi<label))):
      ok=False
  return ok      
    
  
cmd.extend('part2png',part2png)
    
