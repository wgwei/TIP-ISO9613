# -*- coding: utf-8 -*-
"""
Created on Tue Aug 19 10:43:35 2014
ISO 9613-2 Dz
@author: Weigang Wei
"""

import numpy as np

def Dz(C2, C3, wavelength, z, Kmet):
    return 10*np.log10(3+C2/wavelength*C3*z*Kmet)
    
def C2(groundIncluded):
    if groundIncluded==1:
        return 20.0
    else:
        return 40
        
def C3(barWidth, wavelength):
    if barWidth==0: #for single diffraction
        return 1.0
    else:
        return (1+(5*wavelength/barWidth)**2)/(1./3.+(5*wavelength/barWidth)**2)
        
def z(dss, dsr, a, d, barWidth):
    if barWidth==0: # for single diffraction
        return np.sqrt((dss+dsr)**2+a**2) - d
        
    else:
        return np.sqrt((dss+dsr+barWidth)**2+a**2) - d

def Kmet(dss, dsr, d, z):
    if z<=0:
        return 1.0
    else:
        return np.exp(-1./2000*np.sqrt(dss*dsr*d/(2*z)))

def calc_Dz(dss, dsr, barWidth, d, a, wavelength, groundIncluded):
    C2v = C2(groundIncluded)
    C3v = C3(barWidth, wavelength)
    zv = z(dss, dsr, a, d, barWidth)
    Kmetv = Kmet(dss, dsr, d, zv)
    Dzv =  Dz(C2v, C3v, wavelength, zv, Kmetv)
    if Dzv>20 and barWidth==0:
        Dzv = 20
    if Dzv>25 and barWidth>0:
        Dzv = 25
    return Dzv
    
def dist2D(p1, p2):
    return np.sqrt((p1[0]-p2[0])**2. + (p1[1]-p2[1])**2.)
    
def calc_Dz2(spos, rpos, barWidth, barpos, a, wavelength, groundIncluded):
    ''' barpos is the vertices of the barrier. 
        for single diffraction barpos =[x, y, z] or barpos = [x, y]; 
        for double diffraction, barpos = [[x1, y1, z1], [x2, y2, z2]]
        or barpos = [[x1, y1], [x2, y2]]
    '''    
    if barWidth==0:
        dss = dist2D(spos, barpos)
        dsr = dist2D(rpos, barpos)
    else:
        dss = dist2D(spos, barpos[0])
        dsr = dist2D(rpos, barpos[1])
    d = dist2D(spos, rpos)
    C2v = C2(groundIncluded)
    C3v = C3(barWidth, wavelength)
    zv = z(dss, dsr, a, d, barWidth)
    Kmetv = Kmet(dss, dsr, d, zv)
    
    Dzv = Dz(C2v, C3v, wavelength, zv, Kmetv)
    if Dzv>20 and barWidth==0:
        Dzv = 20
    if Dzv>25 and barWidth>0:
        Dzv = 25
    return Dzv
    
if __name__=='__main__':
    pass
    
    
    
    
