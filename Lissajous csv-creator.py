# -*- coding: utf-8 -*-
"""
Created on Sun Apr  1 12:54:36 2018

@author: Jochem van Rabenswaaij
"""

import numpy as np 
import matplotlib.pyplot as plt
import pandas as pd

#Initializations
α_init = 23.5/180*np.pi
β_init = 24/180*np.pi
R_init = 0.82
r = 0.27

vα_init = 0
vβ_init = 0
dt = 0.001
g = 9.81
t_init = 0
t_end = 6

n_steps = int(round((t_end-t_init)/dt)) # total number of timesteps
k = R_init - r

t_arr = np.zeros(n_steps + 1)
a_arr = np.zeros((2,n_steps+1))
v_arr = np.zeros((2,n_steps+1))
α_arr = np.zeros((2,n_steps+1))
x_arr = np.zeros((2,n_steps+1))

t_arr[0] = t_init
v_arr[0,0] = vα_init
v_arr[1,0] = vβ_init
α_arr[0,0] = α_init
α_arr[1,0] = β_init
x_init = r * np.sin(α_init)
y_init = R_init * np.sin(β_init)
x_arr[0,0] = x_init
x_arr[1,0] = y_init

R = np.sqrt( r**2 + k*( k + 2*r*np.cos(α_init) ) )

#Euler's method
for i in range (1, n_steps + 1): 
    #Load everything
    t = t_arr[i-1]
    vα = v_arr[0,i-1]
    vβ = v_arr[1,i-1]
    α = α_arr[0,i-1]
    β = α_arr[1,i-1]
    #Calculate new values
    aα_new = -g*np.cos(β)*np.sin(α) / r
    aβ_new = -g*( k + r*np.cos(α) )*np.sin(β) / R**2
    vα_new = vα + aα_new*dt
    vβ_new = vβ + aβ_new*dt
    α_new = α + vα_new *dt
    β_new = β + vβ_new *dt
    R = np.sqrt( r**2 + k*( k + 2*r*np.cos(α_new) ) )
    x_new = r*np.sin(α_new)
    y_new = R*np.sin(β_new)
    #Store new values
    a_arr[0,i] = aα_new
    a_arr[1,i] = aβ_new
    v_arr[0,i] = vα_new
    v_arr[1,i] = vβ_new
    α_arr[0,i] = α_new
    α_arr[1,i] = β_new
    x_arr[0,i] = x_new
    x_arr[1,i] = y_new
    t_arr[i] = t + dt     

#np.savetxt("metingen.csv", np.column_stack((t_arr[:], α_arr[0,:], α_arr[1,:])), delimiter=",", fmt='%s', header='time,x-hoek,y-hoek')

csv_file = pd.read_csv('r4.csv')


#Plot the results
fig = plt.figure() 
plt.plot(v_arr[1,:], v_arr[0,:], linewidth = 2, label = "mod")
plt.plot(csv_file['wy'], csv_file['wz'], linewidth = 2, label = "met")
#plt.plot(csv_file['time'], csv_file['wy'], linewidth = 2, label = "wz")

#plt.title('Model', fontsize = 20) # set title 
plt.xlabel('x', fontsize = 20) # name of horizontal axis 
plt.ylabel('y', fontsize = 20) # name of vertical axis
plt.xticks(fontsize = 15) # adjust the fontsize 
plt.yticks(fontsize = 15) # adjust the fontsize 
plt.axis([-1.5, 1.5, -1.5, 1.5]) # set the range of the axes
plt.legend(fontsize=10) # show the legend 
plt.show() # necessary for some platforms

fig.savefig('Lissajousoud.jpg', dpi=fig.dpi, bbox_inches = "tight")


fig = plt.figure() 
plt.plot(t_arr[:], v_arr[0,:], linewidth = 2, label = "mod")
plt.plot(csv_file['time'] - 0.1, csv_file['wz'], linewidth = 2, label = "α")
#plt.plot(csv_file['time'], csv_file['wy'], linewidth = 2, label = "wz")

#plt.title('Model', fontsize = 20) # set title 
plt.xlabel('x', fontsize = 20) # name of horizontal axis 
plt.ylabel('y', fontsize = 20) # name of vertical axis
plt.xticks(fontsize = 15) # adjust the fontsize 
plt.yticks(fontsize = 15) # adjust the fontsize 
plt.axis([0, 5, -2, 2]) # set the range of the axes
plt.legend(fontsize=10) # show the legend 
plt.show() # necessary for some platforms

fig.savefig('Lissajousoud1.jpg', dpi=fig.dpi, bbox_inches = "tight")


fig = plt.figure() 
plt.plot(t_arr[:], v_arr[1,:], linewidth = 2, label = "mod")
plt.plot(csv_file['time'] - 0.1, csv_file['wy'] * -1, linewidth = 2, label = "β")
#plt.plot(csv_file['time'], csv_file['wy'], linewidth = 2, label = "wz")

#plt.title('Model', fontsize = 20) # set title 
plt.xlabel('x', fontsize = 20) # name of horizontal axis 
plt.ylabel('y', fontsize = 20) # name of vertical axis
plt.xticks(fontsize = 15) # adjust the fontsize 
plt.yticks(fontsize = 15) # adjust the fontsize 
plt.axis([0, 5, -2, 2]) # set the range of the axes
plt.legend(fontsize=10) # show the legend 
plt.show() # necessary for some platforms

fig.savefig('Lissajousoud2.jpg', dpi=fig.dpi, bbox_inches = "tight")