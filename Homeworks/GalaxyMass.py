#!/usr/bin/env python
# coding: utf-8

# In[12]:


from ReadFile import Read
import numpy as np


# In[9]:


def ComponentMass(filename, particle_type):
    _,_,data = Read(filename)
    index = np.where(data['type'] == particle_type)[0]
    mass = np.sum(data["m"][index])
    mass = mass/100
    round_mass = np.round(mass,3)
    return round_mass


# In[10]:


M31_mass = ComponentMass("M31_000.txt",1)


# In[11]:


M31_mass


# In[ ]:


p

