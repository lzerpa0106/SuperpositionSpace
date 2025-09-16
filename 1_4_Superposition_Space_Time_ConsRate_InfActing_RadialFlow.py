#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 16 10:01:53 2025

Application of the Superposition principle to solutions of the diffusivity equation

Because the diffusivity equation is linear, the principle of superposition may 
be used to compute solutions for complex boundary conditions using only linear 
combinations of solutions for relatively simple boundary conditions. 

@author: luiszerpa
"""

# %% Superposition in space

''' Principle of superposition to calculate the pressure response at any point 
in an infinite-acting reservoir caused by production of multiple wells '''

import numpy as np
import scipy.special as sc
import matplotlib.pyplot as plt
import streamlit as st

def SPspace_InfActingRadialFlow(nwells,r,time,p_i,q,B,viscosity,totalCompressibility,porosity,permeability,thickness):
    '''This equation calculate the pressure response at any point in an 
    infinite-acting reservoir, for N wells producing at different rates (in 
    array q) at a given time. The distance to an specific point from each well 
    is given as the r array.'''
    
    p = p_i
    
    for j in range(0,nwells):
        
        p = p + 70.6* q[j] *B*viscosity/(permeability*thickness) * \
            sc.expi(-948*porosity*viscosity*totalCompressibility* r[j]**2/(permeability*time))
    
    return p


st.title("Application of the Superposition principle to solutions of the diffusivity equation in space")

# Input
nwells = 2
#r = np.array([100, 200]) # Distance in ft from wells to specific point within reservoir
r = np.zeros((2))
r[0] = st.slider("r_1", 0.0, 10000.0, 500.0)
r[1] = st.slider("r_2", 0.0, 10000.0, 500.0)

#q = np.array([1000, 2500]) # Production rate in STB/day for the wells
q = np.zeros((2))
q[0] = st.slider("q_1", 0.0, 10000.0, 500.0)
q[1] = st.slider("q_2", 0.0, 10000.0, 500.0)

p_i = 7500.0 # initial reservoir pressure, psia
B = 1.205 # formation volume factor, RB/STB
viscosity = 2.50 # oil viscosity, cp
totalCompressibility = 5.79e-6 # 1/psi
porosity = 0.22
permeability = 51.6 # md
thickness = 55.8 # reservoir thickness, ft

time = np.arange(0.1,10000.0,0.1)

p_response = SPspace_InfActingRadialFlow(nwells,r,time,p_i,q,B,viscosity,totalCompressibility,porosity,permeability,thickness)

fig, ax = plt.subplots(figsize=(8, 5))  # Create a figure containing a single axes.
ax.plot(time,p_response,linewidth=3.0) # Plot data on the axes.
ax.plot(time,p_i*np.ones(np.size(time)),'--r',linewidth=1.0)
ax.set_xlabel('Time (hr)',fontsize='12',fontweight='bold', fontname='Verdana')
ax.set_ylabel('Pressure (psia)',fontsize='12',fontweight='bold',fontname='Verdana')
ax.tick_params(direction='in',width=1.5,length=8.0, labelsize=12.0)
ax.set_xscale('log')
ax.set(xlim=(0.1, 10000.0), ylim=(4000.0,8000.0))
for tick in ax.get_xticklabels():
    tick.set_fontname('Verdana')
for tick in ax.get_yticklabels():
    tick.set_fontname('Verdana')
ax.grid(True)   
st.pyplot(fig)
