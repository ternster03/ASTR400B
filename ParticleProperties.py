#!/usr/bin/env python
# coding: utf-8

# In[6]:


from ReadFile import Read # Read in the Read file for useage
import numpy as np # Import numpy
from astropy import units as u # Import units from astropy


# In[7]:


def ParticleInfo(filename,particle_type,particle_number):
    '''
    This function takes a file, a particle type, particle number
    and gives the maginitude of distance, magnitude of velocity and 
    mass in units of solar mass
    
    Input: File, particle number, particle type
    
    Output: Magnitude of distnace, magnitude of velocity, and mass
    in units of Earth mass
    '''
    _, _, data = Read(filename) # Gathers the information needed contained
    # in data
    
    index = np.where(data['type'] == particle_type)[0]
    
    x_pos = data['x'][index] * u.kpc # x position
    y_pos = data['y'][index]  * u.kpc # y position
    z_pos = data['z'][index]  * u.kpc # z position
    
    x_new = x_pos[particle_number] # Chooses actually x,y,z pos based on particle
    y_new = y_pos[particle_number]
    z_new = z_pos[particle_number]
    
    mag_dis = np.sqrt(x_new**2 + y_new**2 + z_new**2) # Compute the magnitude
    # distance in kpc
    mag_dis = np.around(mag_dis,3) # Round to 3 decimal places
    
    vx = data['vx'][index]  * u.km / u.s # x velocity
    vy = data['vy'][index]  * u.km / u.s # y velocity
    vz = data['vz'][index]  * u.km / u.s # z velocity
    
    vx_new = vx[particle_number] # Chooses actually vx,vy,vz based on particle
    vy_new = vy[particle_number]
    vz_new = vz[particle_number]
    
    mag_vel = np.sqrt(vx_new**2 + vy_new**2 + vz_new**2) # Compute the magnitude
    # velocity in km/s
    mag_vel = np.around(mag_vel,3) # Round to 3 decimal places
    
    mass = data["m"][index] # Get mass from data from list of part_type
    new_mass = mass[particle_number] * 10**10 * u.M_sun # Get actual mass
    
    return mag_dis, mag_vel, new_mass # return our desired values
    


# In[ ]:





# In[ ]:





# In[ ]:




