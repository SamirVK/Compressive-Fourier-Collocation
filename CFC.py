#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 27 12:36:57 2023

@author: Samir Karam
"""

import numpy as np
import numpy.random as rd
from scipy.integrate import dblquad
from scipy.stats import gmean

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
        u += c_hat[j]*fourier(nux, nuy, x, y)/(4*np.pi**2*euc2(nux, nuy))
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
def CFC(colloc_points, cardinal_lambda, diff_num, exact_sparsity, recovery):
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
        b[i] = (s1 + s2 + s3 + s4)/np.sqrt(M)
        
    if recovery == 'OLS':
        c_hat = np.linalg.lstsq(A, b, rcond=None)
        return c_hat[0]
    else:
        K = int(M/2)
        c_hat = OMP(A, b, K, N)
        return c_hat

## ORTHOGONAL MATCHING PURSUIT (OMP)
## OMP takes A, b, K, as input and returns a K-sparse coefficient vector
def OMP(A, b, K, N):
    # Normalize the columns of A
    A = A/np.sqrt(np.sum(np.square(np.abs(A))))
    # Initialize z with empty support
    S = []
    z = np.zeros(N, dtype=complex)
    for n in range(K):
        # Pick the index whose column of A is mostly correlated with the 
        # residual due to the current approximation z[n]
        r = b - A.dot(z) 
        j = np.argmax(np.abs((A.conj().T).dot(r)))
        S.append(j)
        z[S] = np.linalg.pinv(A[:,S]).dot(b)              
    c_hat = z    
    print(np.linalg.norm(c_hat,ord=0))
    return c_hat

## PLOT HYPERBOLIC SET AND EXACT SOLUTIONS
def exactPlots():
    plt.xlabel('x')
    plt.ylabel('y')
    plt.scatter(*zip(*lam), s = 1, marker = '.')
    
    #plt.savefig('hyperbolic_n=39', dpi = 300)
    plt.show()
    
    # PLOT EXACT SOLUTIONS
    x = np.outer(np.linspace(0,1,200), np.ones(200))
    y = x.copy().T
    z_true = u(x,y,'nonsparse')
    
    # Plot z_true
    axi = plt.axes(projection = '3d')
    axi.set_xlabel('x')
    axi.set_ylabel('y')
    axi.set_zlabel('u2')
    axi.plot_surface(x,y,z_true,cmap='plasma', cstride=1, rstride = 1,
                     alpha=None, linewidth = 0, antialiased=False)
    axi.view_init(15,235)
    #plt.savefig('Exact_u2', bbox_inches="tight", dpi=300)
    plt.show()
    
# PLOT APPROXIMATE SOLUTIONS
def approxPlots():      
    x = np.outer(np.linspace(0,1,100), np.ones(100))
    y = x.copy().T
    
    M = 2**9
    c_hat = CFC(M, cardlam, 1, 'sparse')
    z_approx = u_approx(x, y, c_hat, lam, cardlam)
    
    # Plot z_approx
    axi = plt.axes(projection = '3d')
    axi.set_xlabel('x')
    axi.set_ylabel('y')
    axi.set_zlabel('u1')
    axi.plot_surface(x,y,z_approx,cmap='plasma', cstride = 1, rstride = 1,
                     alpha=None, antialiased=True)
    axi.view_init(15,235)
    #plt.savefig('Approx_u1_a1', bbox_inches="tight", dpi=300)
    plt.show()

# Define the hyperbolic truncation set
n = 39
lam = [] # lexicographic ordering of the indices
for x in range(-n+1, n):
    for y in range(-n+1,n):
        if (np.abs(x)+1)*(np.abs(y)+1) <= n and not (x == 0 and y == 0):
            lam.append(np.array([x,y]))
cardlam = len(lam)

##############################################################
## Reproduce the results of the L2 loss plots! To do this, run 
## an outer loop over M = [8, 16, 32, 64, 128, 256, 512, 1024]
## and an inner loop that calculates the geometric mean of the 
## L2 loss over 25 runs for each value of M (as does Weiqi).

sparsity = 'sparse' ############## <<<<<<<<< Change sparsity
diff = 1 ###################### <<<<<<<<< Change diff coeff
r = 'OMP' ################## <<<<<<<<<< Change recovery method
L2_losses = []
collocation_points = [8,16,32,64,128,256,512]
for M in collocation_points:
    test_losses = []
    for test in range(1):
        
        # Provide random parameters (global variables)
        d = rd.rand(10)
        m = rd.randint(1, 6, size = 10)
        n = rd.randint(1, 6, size = 10)
        
        # Randomly sample points to evaluate the loss
        loss_points = rd.rand(2*cardlam,2)
        x = loss_points[:,0]
        y = loss_points[:,1]
        
        # Get the exact solution at the random points
        u_exact = u(x,y,sparsity)
        
        # Get the approx solution at the random points
        c_hat = CFC(M, cardlam, diff, sparsity, r)
        u_hat = u_approx(x, y, c_hat, lam, cardlam)
        
        # Compute the loss for this trial
        loss = (np.sqrt(sum(abs(u_exact-u_hat)**2)/np.sqrt(M))/
                np.sqrt(sum(abs(u_exact)**2)/np.sqrt(M)))
        test_losses.append(loss)

    # Compute the geometric mean of the test_losses
    L2_losses.append(gmean(test_losses))

plt.figure()
plt.plot(collocation_points, L2_losses, marker = '.', color='black', label='OLS')
plt.axvline(x = 444, linestyle='dotted', color='gray', label=r'$m=|\Lambda|$')
plt.yscale('log')
plt.legend()
plt.show

