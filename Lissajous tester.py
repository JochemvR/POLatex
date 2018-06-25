# -*- coding: utf-8 -*-
"""
Created on Sun Apr  1 12:54:36 2018

@author: Jochem van Rabenswaaij
"""

import numpy as np 
import matplotlib.pyplot as plt
import pandas as pd

#Initializations
α_init = 0.086635778
β_init = 0.1542188
AP = 0.97
a = 0.255

αv_init = 0
βv_init = 0
dt = 0.001
g = 9.81
t_init = 0
t_end = 9

n_steps = int(round((t_end-t_init)/dt)) # total number of timesteps
c3 = -AP - a

t_arr = np.zeros(n_steps + 1)
a_arr = np.zeros((2,n_steps+1))
v_arr = np.zeros((2,n_steps+1))
α_arr = np.zeros((2,n_steps+1))
x_arr = np.zeros((2,n_steps+1))

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
    
    αa_new = -g * abs(np.sin(α)) * np.sqrt((b3 * p2)**2 + (b3 * p1 - b1 * p3)**2
                      + (b1 * p2)**2) * abs(b3) / (np.sin(α) * (b1**2 + b3**2)**2 )
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
    x_arr[0,i-1] = p2
    x_arr[1,i-1] = p1
    t_arr[i] = t + dt  

x_arr[0,i] = p2
x_arr[1,i] = p1

csv_file = pd.read_csv('r6.csv')

#Plot the results
fig = plt.figure() 
plt.plot(x_arr[1,:], x_arr[0,:], linewidth = 2, label = "model")
plt.plot(csv_file['P1Y'], csv_file['P1X'], linewidth = 2, label = "metingen")
#plt.plot(csv_file['time'], csv_file['wy'], linewidth = 2, label = "wz")

#plt.title('Model', fontsize = 20) # set title 
plt.xlabel('x [m]', fontsize = 20) # name of horizontal axis 
plt.ylabel('y [m]', fontsize = 20) # name of vertical axis
plt.xticks(fontsize = 15) # adjust the fontsize 
plt.yticks(fontsize = 15) # adjust the fontsize 
plt.axis([-.2, .2, -.2, .2]) # set the range of the axes
plt.legend(fontsize=10) # show the legend 
plt.show() # necessary for some platforms

fig.savefig('Lissajous.jpg', dpi=fig.dpi, bbox_inches = "tight")


fig = plt.figure() 
plt.plot(t_arr[:], x_arr[0,:], linewidth = 2, label = "model")
plt.plot(csv_file['time']-0.4, csv_file['P1X'], linewidth = 2, label = "metingen")
#plt.plot(csv_file['time'], csv_file['wy'], linewidth = 2, label = "wz")

#plt.title('Model', fontsize = 20) # set title 
plt.xlabel('time [s]', fontsize = 20) # name of horizontal axis 
plt.ylabel('y [m]', fontsize = 20) # name of vertical axis
plt.xticks(fontsize = 15) # adjust the fontsize 
plt.yticks(fontsize = 15) # adjust the fontsize 
plt.axis([0, 5, -0.2, .2]) # set the range of the axes
plt.legend(fontsize=10) # show the legend 
plt.show() # necessary for some platforms

fig.savefig('Lissajous1.jpg', dpi=fig.dpi, bbox_inches = "tight")


fig = plt.figure() 
plt.plot(t_arr[:], x_arr[1,:], linewidth = 2, label = "model")
plt.plot(csv_file['time']-0.4, csv_file['P1Y'], linewidth = 2, label = "metingen")
#plt.plot(csv_file['time'], csv_file['wy'], linewidth = 2, label = "wz")

#plt.title('Model', fontsize = 20) # set title 
plt.xlabel('time [s]', fontsize = 20) # name of horizontal axis 
plt.ylabel('x [m]', fontsize = 20) # name of vertical axis
plt.xticks(fontsize = 15) # adjust the fontsize 
plt.yticks(fontsize = 15) # adjust the fontsize 
plt.axis([0, 5, -.2, .2]) # set the range of the axes
plt.legend(fontsize=10) # show the legend 
plt.show() # necessary for some platforms

fig.savefig('Lissajous2.jpg', dpi=fig.dpi, bbox_inches = "tight")