#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 13 11:28:33 2019
@author: brent
"""
import numpy as np
import random as rand
#import matplotlib.pyplot as plt
import scipy.special as sp


rand.seed(8822731)

file = open('2T_18GHz.txt', 'a')

#-------------------Physical Constants-------------------#

pi = 3.14159265

c = 299792458

mu_0 = 12.566E-7 # N/A^2 or H/m

q = 1.602E-19 # electron charge in Coulombs

points = 20000 # Number of points to simulate


#--------------------Input Parameters--------------------#

omega = 2.0 * pi * 18.0E9 # Angular frequency

B_0 = 2.0  # Field in Tesla

a = 0.005625 # Guide radius in meters


#-------------------Derived Parameters-------------------#

gamma = (B_0 * q)/(omega * 9.109E-31) # Lorentz factor

beta = np.sqrt(1-(1/gamma**2)) # = v/c

r_cy = (c * beta)/omega # Relativistic cyclotron radius

k_cut = 1.841/a # Cutoff wavenumber in m^{-1}

zeta = np.sqrt((omega/c)**2-k_cut**2) # Wavenumber of TE11 mode at given omega, guide radius

Dt = 2*pi/(omega*360) #time step to ensure integration takes place over 1 orbital period

e_inside = 0
e_conf = 0

data=[]

for i in range(points):
    # choose random point in or near circle of radius a
    #x_c, y_c will be coordinates of the CENTER of the orbit    
    x_c = rand.uniform(-1.0*a, 1.0*a)
    y_c = rand.uniform(-1.0*a, 1.0*a)
    
    rho_c = np.sqrt(x_c**2+y_c**2) # radial coordinate of the orbit center
    
    if rho_c < a:
        e_inside = e_inside + 1
    
    if rho_c + r_cy < a:  # check combined radius of orbit center and cyclotron radius
        
        e_conf = e_conf + 1
        
        A = 0.0
        B = 0.0
        
        if x_c == 0.0: x_c = 0.00001
        if y_c == 0.0: y_c = 0.00001
        
        if ((x_c + r_cy > 0.0) & (y_c > 0.0)):    
            phi_0 = 0.0 * pi + np.arctan(y_c/(x_c + r_cy))
        if ((x_c + r_cy > 0.0) & (y_c > 0.0)):    
            phi_0 = 1.0 * pi + np.arctan(y_c/(x_c + r_cy))
        if ((x_c + r_cy > 0.0) & (y_c < 0.0)):    
            phi_0 = 1.0 * pi + np.arctan(y_c/(x_c + r_cy))
        if ((x_c + r_cy > 0.0) & (y_c < 0.0)):    
            phi_0 = 2.0 * pi + np.arctan(y_c/(x_c + r_cy))
        
        theta_0 = rand.uniform(0.0, 2.0 * pi) #random tangent vector for initial velocity

        for t in range(360): # calculate each orbit in 360 steps with time intervals of Dt
            
            x = x_c + r_cy * np.cos(omega*t*Dt + theta_0)
            y = y_c + r_cy * np.sin(omega*t*Dt + theta_0)
            
            rho = np.sqrt(x**2+y**2)
            
            if ((x > 0.0) & (y > 0.0)):    
                phi = 0.0 * pi + np.arctan(y/x)
            if ((x < 0.0) & (y > 0.0)):    
                phi = 1.0 * pi + np.arctan(y/x)
            if ((x < 0.0) & (y < 0.0)):    
                phi = 1.0 * pi + np.arctan(y/x)
            if ((x > 0.0) & (y < 0.0)):    
                phi = 2.0 * pi + np.arctan(y/x)
            
            j_x = q * -1.0 * r_cy * omega * np.sin(omega*t*Dt + theta_0)
            j_y = q *  1.0 * r_cy * omega * np.cos(omega*t*Dt + theta_0)
            
            j_rho = np.cos(phi) * j_x + np.sin(phi) * j_y
            j_phi = np.cos(phi) * j_y - np.sin(phi) * j_x
            
            E_rho_A = omega*mu_0/(k_cut**2*rho) *  1.0 * np.cos(phi) * sp.jv(1,k_cut*rho)  * np.sin(omega*t*Dt)
            E_rho_B = omega*mu_0/(k_cut**2*rho) * -1.0 * np.sin(phi) * sp.jv(1,k_cut*rho)  * np.sin(omega*t*Dt)
            
            E_phi_A = -1.0 * omega*mu_0/(2*k_cut) * np.sin(phi) * (sp.jv(0,k_cut*rho) - sp.jv(2,k_cut*rho))  * np.sin(omega*t*Dt)
            E_phi_B = -1.0 * omega*mu_0/(2*k_cut) * np.cos(phi) * (sp.jv(0,k_cut*rho) - sp.jv(2,k_cut*rho))  * np.sin(omega*t*Dt)
            
            A = A + Dt * (E_rho_A * j_rho + E_phi_A * j_phi) 
            B = B + Dt * (E_rho_B * j_rho + E_phi_B * j_phi)
                
        P_Nmw = (pi*omega*mu_0*zeta/(2*k_cut**4)) * (1.841**2-1.0) * sp.jv(1,k_cut*a)  *sp.jv(1,k_cut*a)
#        print("P_Nmw = ", P_Nmw)
        P_11 = ((pi*zeta)/(omega * mu_0*k_cut**2)) * (1.841**2 - 1) *sp.jv(1,1.841)
 #       print("P_11 = ", P_11)
        A = A * -1.0/(360 * Dt * P_Nmw)
        B = B * -1.0/(360 * Dt * P_Nmw)        
        P_tot = 2 * (A**2 + B**2) * P_Nmw
        data.append(P_tot)
        s = str(P_tot)
        file.write(s)
        file.write('\n')

print("N =",len(data))

#plt.figure(figsize=(20,10))    
#num_bins = 60
#n, bins, patches = plt.hist(data, num_bins, facecolor='blue', alpha=0.5)
#plt.xlabel('Power')
#plt.ylabel('Counts')
#plt.title(r'Histogram') 
# Tweak spacing to prevent clipping of ylabel
#plt.subplots_adjust(left=0.15)
#plt.show()

file.close()
