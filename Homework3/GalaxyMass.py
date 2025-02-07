#!/usr/bin/env python
# coding: utf-8


from ReadFile import Read # Import all required packages
import numpy as np


def ComponentMass(filename, particle_type):
    '''
    This function takes a given file and a particle type
    and gives a sum of rounded mass to the three decimals for
    the given particle type

    Input:
    filename(txt file): A txt file with data containing some mass colum

    particle_type(integer): A number containing the type of particle within
    the txt file
    
    '''
    _,_,data = Read(filename) # Gets data from the txt file
    index = np.where(data['type'] == particle_type)[0] # Creates index for the given
    # particle type
    mass = np.sum(data["m"][index]) # Gathers the sum of the mass for a give
    # particle type
    mass = mass/100 # Divides mass to get to units of 10**12 MSun 
    round_mass = np.round(mass,3) # Round mass to three decimals
    return round_mass # Returns the mass

