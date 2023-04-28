#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 27 12:36:57 2023

@author: Samir Karam
"""

import numpy as np
import numpy.random as rd
from sklearn import linear_model
from scipy.integrate import quad, dblquad

import matplotlib.pyplot as plt

# Define the Fourier collocation basis
def fourier(nux,nuy,x,y):
    return np.exp(2*np.pi*(nux*x+nuy*y))
def fourierx(nux,nuy,x,y):
    return 
def fourierxx(nux,nuy,x,y):
    return 
def fouriery(nux,nuy,x,y):
    return 
def fourieryy(nux,nuy,x,y):
    return

# Define the exact solutions
def u_sparse(x,y):
    u = 0
    for k in range(10):
        u += d[k]*np.sin(2*np.pi*m[k]*x)*np.sin(2*np.pi*n[k]*y)
    return u
def u_nonsparse(x,y):
    def I():
        return dblquad(lambda y, x: np.exp(np.sin(2*np.pi*x)+
            np.sin(2*np.pi*y)),0 , 1, lambda x: 0, lambda x: 1)
    u = np.exp(np.sin(2*np.pi*x)+np.sin(2*np.pi*y)) - I()[0]
    return u

# Define the diffusion coefficients and their derivatives
def a_1(x,y):
    return 1
def a_2(x,y):
    return 1+0.25*np.sin(2*np.pi*x)*np.sin(2*np.pi*y)+0.25*np.sin(4*np.pi*x)
def a_2x(x,y):
    return 
def a_2y(x,y):
    return 
def a_3(x,y):
    return 1+0.2*np.exp(np.sin(2*np.pi*x)*np.sin(2*np.pi*y))
def a_3x(x,y):
    return 
def a_3y(x,y):
    return 

# Define and plot the hyperbolic truncation set
n = 39
lam = [] # lexicographic ordering of the indices
for x in range(-n+1, n):
    for y in range(-n+1,n):
        if (np.abs(x)+1)*(np.abs(y)+1) <= n and not (x == 0 and y == 0):
            lam.append(np.array([x,y]))

## PLOTTING

plt.xlabel('x')
plt.ylabel('y')
plt.scatter(*zip(*lam), s = 1, marker = '.')

#plt.savefig('hyperbolic_n=39', dpi = 300)
plt.show()

# PLOT EXACT SOLUTIONS
# Define problem parameters

np.random.seed(0)
d = rd.rand(10)
m = rd.randint(1, 6, size = 10)
n = rd.randint(1, 6, size = 10)

x = np.outer(np.linspace(0,1,1000), np.ones(1000))
y = x.copy().T
z_true = u_sparse(x,y)

# Plot z_true
fig = plt.figure()
ax = plt.axes(projection = '3d')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('u1')
ax.plot_surface(x,y,z_true,cmap='plasma')
ax.view_init(15,235)
#plt.savefig('Exact_u1', bbox_inches="tight", dpi=300)
plt.show()

## BEGIN COMPRESSIVE FOURIER COLLOCATION
# Build A



# Build b

