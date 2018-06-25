# -*- coding: utf-8 -*-
"""
Created on Sun Apr 1 12:54:36 2018
@author: Jochem van Rabenswaaij
"""

import numpy as np 
import matplotlib.pyplot as plt

#Initializations
α_init = 23.5 /180*np.pi
β_init = 19.5 /180*np.pi
AP = 0.85
a = 0.27

αv_init = 0
βv_init = 0
dt = 0.001
g = 9.81
t_init = 0
t_end = 10

n_steps = int(round((t_end-t_init)/dt)) # total number of timesteps
c3 = -AP - a

t_arr = np.zeros(n_steps + 1)
a_arr = np.zeros((2,n_steps+1))
v_arr = np.zeros((2,n_steps+1))
α_arr = np.zeros((2,n_steps+1))

t_arr[0] = t_init
v_arr[0,0] = αv_init
v_arr[1,0] = βv_init
α_arr[0,0] = α_init
α_arr[1,0] = β_init

#Euler's method
for i in range (1, n_steps + 1): 
    #Load everything
    t = t_arr[i-1]
    αv = v_arr[0,i-1]
    βv = v_arr[1,i-1]
    α = α_arr[0,i-1]
    β = α_arr[1,i-1]
    #Calculate new values
    p1 = AP * np.sin(β)
    p2 = AP * np.sin(α) * np.cos(β) + a * np.sin(α)
    p3 = -AP * np.cos(α) * np.cos(β) - a * np.cos(α)
    
    b1 = AP * np.sin(β)
    b3 = -AP * np.cos(β) - a
    
    αa_new = -g * abs(np.sin(α)) * np.sqrt((b3 * p2)**2 + (b3 * p1 - b1 * p3)**2 + (b1 * p2)**2) * abs(b3) / (np.sin(α) * (b1**2 + b3**2)**2 )
    βa_new = -g * np.sin(β) * np.cos(α) / AP
    
    αv_new = αv + αa_new * dt
    βv_new = βv + βa_new * dt
    
    α_new = α + αv_new *dt
    β_new = β + βv_new *dt
    #Store new values
    a_arr[0,i] = αa_new
    a_arr[1,i] = βa_new
    v_arr[0,i] = αv_new
    v_arr[1,i] = βv_new
    α_arr[0,i] = α_new
    α_arr[1,i] = β_new
    t_arr[i] = t + dt      
    
#Plot the results
fig = plt.figure() 
#ax = fig.gca(projection='3d')
plt.plot(v_arr[1,:], v_arr[0,:], linewidth = 1)
#plt.plot(t_arr[:], y[:], linewidth = 1, label="Phi vs time")

plt.title('ttt', fontsize = 20) # set title 
plt.xlabel('β', fontsize = 15) # name of horizontal axis 
plt.ylabel('α', fontsize = 15) # name of vertical axis
plt.xticks(fontsize = 15) # adjust the fontsize 
plt.yticks(fontsize = 15) # adjust the fontsize 
plt.axis([-1.5, 1.50, -1.50, 1.50]) # set the range of the axes
#plt.legend(fontsize=15) # show the legend 
plt.show() # necessary for some platforms

fig.savefig('Lissajous.jpg', dpi=fig.dpi, bbox_inches = "tight")
