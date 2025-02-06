#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np # Importing numpy for array useage
import astropy.units as u # Importing astropy for units


# In[9]:


def Read(filename):
    '''
    Function will read in a file and returns the time and
    total number of particles as variables and the particle
    type, mass, x, y, z, vx, vy, and vz columns as an array
    Input: file
    Output: As seen above 
    '''
    file = open(filename,"r") # Reads in the file
    
    line1 = file.readline() # Reads the first line
    label, value = line1.split() # Splits the labels and values
    time = float(value)*u.Myr # Changes the value or time to units
    # of millions of years
    
    line2 = file.readline() # Read the second line
    label, value = line2.split() # Splits the labels and values
    total_particles = int(value)# Take the integer of the total values to get the 
    # total number of particles
    
    file.close() # Close the file
    
    data = np.genfromtxt(filename,dtype=None,names=True,skip_header=3) #
    # Creates the arrays used to store the x,y,z,vx,vy,vz information
    
    return time, total_particles, data # Returns all desired values

