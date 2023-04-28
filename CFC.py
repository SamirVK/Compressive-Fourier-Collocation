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

# A function to compute the squared 2-norm of a vector nu
def euc2(nux,nuy):
    return (nux**2+nuy**2)

# Define the Fourier collocation basis and their derivatives
def fourier(nux,nuy,x,y):
    return np.exp(2j*np.pi*(nux*x+nuy*y))
def fourierx(nux,nuy,x,y):
    return 2j*np.pi*nux*np.exp(2j*np.pi*(nux*x+nuy*y))
def fourierxx(nux,nuy,x,y):
    return -4*np.pi**2*nux**2*np.exp(2j*np.pi*(nux*x+nuy*y))
def fouriery(nux,nuy,x,y):
    return 2j*np.pi*nuy*np.exp(2j*np.pi*(nux*x+nuy*y))
def fourieryy(nux,nuy,x,y):
    return -4*np.pi**2*nuy**2*np.exp(2j*np.pi*(nux*x+nuy*y))

# A function to compute the approximate solution with recovered coeffs
def u_approx(x, y, c_hat, lam, cardlam):
    u = 0
    for j in range(cardlam):
        nux = lam[j][0]
        nuy = lam[j][1]
        u += c_hat[0][j]*fourier(nux, nuy, x, y)/(4*np.pi**2*euc2(nux, nuy))
    return u

# Define the exact solutions and their derivative approximations
def u(x,y,sparsity):
    if sparsity == 'sparse':
        u = 0
        for k in range(10):
            u += d[k]*np.sin(2*np.pi*m[k]*x)*np.sin(2*np.pi*n[k]*y)
        return u
    else:
        def I():
            return dblquad(lambda y, x: np.exp(np.sin(2*np.pi*x)+
                np.sin(2*np.pi*y)),0 , 1, lambda x: 0, lambda x: 1)
        u = np.exp(np.sin(2*np.pi*x)+np.sin(2*np.pi*y)) - I()[0]
        return u
def ux(x,y,sparsity):
    if sparsity == 'sparse':
        u = 0
        for k in range(10):
            u += 2*np.pi*m[k]*d[k]*np.cos(2*np.pi*m[k]*x)*np.sin(2*np.pi*n[k]*y)
        return u
    else:
        u = 2*np.pi*np.cos(2*np.pi*x)*np.exp(np.sin(2*np.pi*x)+np.sin(2*np.pi*y))
        return u
def uy(x,y,sparsity):
    if sparsity == 'sparse':
        u = 0
        for k in range(10):
            u += 2*np.pi*n[k]*d[k]*np.sin(2*np.pi*m[k]*x)*np.cos(2*np.pi*n[k]*y)
        return u
    else:
        u = 2*np.pi*np.cos(2*np.pi*y)*np.exp(np.sin(2*np.pi*x)+np.sin(2*np.pi*y))
        return u
def uxx(x,y,sparsity):
    if sparsity == 'sparse':
        u = 0
        for k in range(10):
            u += -4*np.pi**2*m[k]**2*d[k]*np.sin(2*np.pi*m[k]*x)*np.sin(2*np.pi*n[k]*y)
        return u
    else:
        u = (4*np.pi**2*np.exp(np.sin(2*np.pi*x)+np.sin(2*np.pi*y))*
             (np.cos(2*np.pi*x)**2-np.sin(2*np.pi*x)))
        return u
def uyy(x,y,sparsity):
    if sparsity == 'sparse':
        u = 0
        for k in range(10):
            u += -4*np.pi**2*n[k]**2*d[k]*np.sin(2*np.pi*m[k]*x)*np.sin(2*np.pi*n[k]*y)
        return u
    else:
        u = (4*np.pi**2*np.exp(np.sin(2*np.pi*x)+np.sin(2*np.pi*y))*
             (np.cos(2*np.pi*y)**2-np.sin(2*np.pi*y)))
        return u 

# Define the diffusion coefficients and their derivatives
def a(x,y,num):
    if num == 1:
        return 1
    if num == 2:
        return 1+0.25*np.sin(2*np.pi*x)*np.sin(2*np.pi*y)+0.25*np.sin(4*np.pi*x)
    if num == 3:
        return 1+0.2*np.exp(np.sin(2*np.pi*x)*np.sin(2*np.pi*y)) 
    return
def ax(x,y,num):
    if num == 1:
        return 0
    if num == 2: 
        return np.pi*np.cos(4*np.pi*x)+np.pi/2*np.cos(2*np.pi*x)*np.sin(2*np.pi*y)
    if num == 3: 
        return (2/5*np.exp(np.sin(2*np.pi*x)*np.sin(2*np.pi*y))*np.pi
                *np.cos(2*np.pi*x)*np.sin(2*np.pi*y))
    return 
def ay(x,y,num):
    if num == 1:
        return 0
    if num == 2: 
        return np.pi/2*np.cos(2*np.pi*y)*np.sin(2*np.pi*x)
    if num == 3:
        return (2/5*np.exp(np.sin(2*np.pi*x)*np.sin(2*np.pi*y))*np.pi
                *np.cos(2*np.pi*y)*np.sin(2*np.pi*x))
    
## BEGIN COMPRESSIVE FOURIER COLLOCATION
## CFC takes parameters as input and outputs a vector of recovered coeffs
def CFC(colloc_points, cardinal_lambda, diff_num, exact_sparsity):
    num = diff_num
    sparsity = exact_sparsity
    M = colloc_points  
    N = cardinal_lambda
    colloc = rd.rand(M,2)
    
    # Build A in R^(M x N)
    A = np.zeros([M, N], dtype=complex)
    for i in range(M):
        x = colloc[i][0]
        y = colloc[i][1]
        for j in range(N):
            nux = lam[j][0]
            nuy = lam[j][1]
            
            c = -1/(np.sqrt(M)*4*np.pi**2*euc2(nux,nuy))
            s1 = ax(x,y,num)*fourierx(nux,nuy,x,y)
            s2 = a(x,y,num)*fourierxx(nux,nuy,x,y)
            s3 = ay(x,y,num)*fouriery(nux,nuy,x,y)
            s4 = a(x,y,num)*fourieryy(nux,nuy,x,y)
            A[i,j] = c*(s1+s2+s3+s4)
            
    # Build b
    b = [0]*(M)
    for i in range(M):
        x = colloc[i][0]
        y = colloc[i][1]
        s1 = -ax(x,y,num)*ux(x,y,sparsity)
        s2 = -a(x,y,num)*uxx(x,y,sparsity)
        s3 = -ay(x,y,num)*uy(x,y,sparsity)
        s4 = -a(x,y,num)*uyy(x,y,sparsity)
        b[i] = s1 + s2 + s3 + s4

    c_hat = np.linalg.lstsq(A, b, rcond=None)
    return c_hat

## PLOTTING
def exactPlots():
    plt.xlabel('x')
    plt.ylabel('y')
    plt.scatter(*zip(*lam), s = 1, marker = '.')
    
    #plt.savefig('hyperbolic_n=39', dpi = 300)
    plt.show()
    
    # PLOT EXACT SOLUTIONS
    x = np.outer(np.linspace(0,1,1000), np.ones(1000))
    y = x.copy().T
    z_true = u(x,y,'sparse')
    
    # Plot z_true
    fig = plt.figure()
    axi = plt.axes(projection = '3d')
    axi.set_xlabel('x')
    axi.set_ylabel('y')
    axi.set_zlabel('u1')
    axi.plot_surface(x,y,z_true,cmap='plasma')
    axi.view_init(15,235)
    #plt.savefig('Exact_u1', bbox_inches="tight", dpi=300)
    plt.show()

# Define the hyperbolic truncation set
n = 39
lam = [] # lexicographic ordering of the indices
for x in range(-n+1, n):
    for y in range(-n+1,n):
        if (np.abs(x)+1)*(np.abs(y)+1) <= n and not (x == 0 and y == 0):
            lam.append(np.array([x,y]))
cardlam = len(lam)

# PROBLEM PARAMETERS
np.random.seed(0)
d = rd.rand(10)
m = rd.randint(1, 6, size = 10)
n = rd.randint(1, 6, size = 10)

x = np.outer(np.linspace(0,1,1000), np.ones(1000))
y = x.copy().T

M = 2**9
c_hat = CFC(M, cardlam, 3, 'nonsparse')
z_approx = u_approx(x, y, c_hat, lam, cardlam)

# Plot z_true
fig = plt.figure()
axi = plt.axes(projection = '3d')
axi.set_xlabel('x')
axi.set_ylabel('y')
axi.set_zlabel('u1')
axi.plot_surface(x,y,z_approx,cmap='plasma')
axi.view_init(15,235)
#plt.savefig('Exact_u1', bbox_inches="tight", dpi=300)
plt.show()







