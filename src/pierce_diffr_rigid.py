# -*- coding: utf-8 -*-
"""
Created on Mon Dec 15 19:32:30 2014
calculate the rigid thick barrier with diffraction angle 1.5*pi
@author: Weigang Wei
"""
import numpy as np
from scipy import special
 
def f(x):
    x = np.abs(x)
    SC = special.fresnel(x)
    S = SC[0]  # fresnel sin
    C = SC[1]  # fresnel cos
    fx = (0.5-S)*np.cos(0.5*np.pi*x**2) - (0.5-C)*np.sin(0.5*np.pi*x**2)    
    return fx

def g(x):
    x = np.abs(x)
    SC = special.fresnel(x)
    S = SC[0]  # fresnel sin
    C = SC[1]  # fresnel cos
    gx = (0.5-C)*np.cos(0.5*np.pi*x**2) + (0.5-S)*np.sin(0.5*np.pi*x**2)
    return gx 

def calcualte_argin(rs, rr, thetas, thetar, w, waveLen):
    ''' the difault diffraction angles are 1.5pi'''
    L = rs + rr + w
    Ys = np.sqrt(3)*(np.cos(2.0/3.0*thetas)-0.5) * np.sqrt(2*rs*(w+rr)/(waveLen*L))
    Yr = np.sqrt(3)*(np.cos(2.0/3.0*thetar)-0.5) * np.sqrt(2*rr*(w+rs)/(waveLen*L))
    if Ys > Yr:
        Ygr = Ys   # greater Y
        Ysm = np.sqrt(3.0)*(np.cos(2.0/3.0*thetar)-0.5) * np.sqrt(2*rr*w/(waveLen*(w+rr)))   # smaller Y
    else:
        Ygr = Yr   # greater Y
        Ysm = np.sqrt(3.0)*(np.cos(2.0/3.0*thetas)-0.5) * np.sqrt(2*rs*w/(waveLen*(w+rs)))   # smaller Y
    return [Ys, Yr, Ygr, Ysm]            
        
       
def LdiffOverLff_theo(rs, rr, thetas, thetar, w, waveLen):
    [Ys, Yr, Ygr, Ysm] = calcualte_argin(rs, rr, thetas, thetar, w, waveLen)
    pSquare = (f(Ygr)**2 + g(Ygr)**2)*(f(Ysm)**2 + g(Ysm)**2)
    return 10.0*np.log10(pSquare)
    
def LdiffOverLff_simple(rs, rr, thetas, thetar, w, waveLen):
    [Ys, Yr, Ygr, Ysm] = calcualte_argin(rs, rr, thetas, thetar, w, waveLen)
    pSquare = (0.37/(0.37+Ygr))**2. * (0.37/(0.37+Ysm))**2.
    return 10.0*np.log10(pSquare)    
    

if __name__=='__main__':
    print LdiffOverLff_theo(10., 10., 0.25*np.pi, 0.25*np.pi, 10, 1.5)
