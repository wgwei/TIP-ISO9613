# -*- coding: utf-8 -*-
"""
Created on Thu Dec 11 21:51:58 2014
calculate diffraction by Pierce's theory
@author: Weigang Wei
"""

from scipy import special
from scipy.special import erfc
import numpy as np
import math
import matplotlib.pylab as plt

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
    
    
def  p2edge(fr, betai, thetaiS, thetaiR, beta2, theta2S, theta2R,\
    tauS, tauM, tauR, Rs, Rm, Rr, Zis, Z2r, \
    spos, rpos, diffEdgePosi, diffEdgePos2):
    '''' Zis is the surface impedance of the source side. it can be complex number
         Z2r is the surface impedance of the receiver side. It can be complex number
    '''
    tau = tauS+tauM+tauR
    el = Rs+Rm+Rr
    
    mainEdge = find_main_edge(spos, rpos, diffEdgePosi, diffEdgePos2)
    if mainEdge==1:  #first edge is primary
#        print 'mainEdge is 1'
        Dfi = Dwedge_nord2000(fr, betai, thetaiS, thetaiR, tau, tauS, tauM+tauR, el, Rs, Rm+Rr, Zis, Z2r)
        Df2 = Dwedge_nord2000(fr, beta2, theta2S, theta2R, tauM+tauR, tauM, tauR, Rm+Rr, Rm, Rr, Zis, Z2r)
    else:
#        print 'mainEdge is 2'
        Dfi = Dwedge_nord2000(fr, betai, thetaiS, thetaiR, tauS+tauM, tauS, tauM, Rs+Rm, Rs, Rm, Zis, Z2r)
        Df2 = Dwedge_nord2000(fr, beta2, theta2S, theta2R, tau, tauS+tauM, tauR, el, Rs+Rm, Rr, Zis, Z2r)
#    print 'Dfi, Df2', Dfi, Df2
    pf = 0.5*Dfi*Df2*np.exp(1j*2*np.pi*fr*tau)/el
    return pf

def p_diffr_ref_freefield(p2edgeValue, el, fr, tau):
    return p2edgeValue/(np.exp(1j*2*np.pi*fr*tau)/el)
    

def find_main_edge(spos, rpos, diffEdgePosi, diffEdgePos2):
    zmi = (rpos[1]-spos[1]) * (diffEdgePosi[0] - spos[0])/(rpos[0]-spos[0]) + spos[1]
    SPi = dist2D(spos, diffEdgePosi)
    RPi = dist2D(rpos, diffEdgePosi)
    SR = dist2D(spos, rpos)
    deltaL0 = H_switch(diffEdgePosi[1], zmi) * (SPi+RPi-SR)
    
    SP2 = dist2D(spos, diffEdgePos2)
    RP2 = dist2D(rpos, diffEdgePos2)
    deltaLi = H_switch(diffEdgePos2[1], zmi) * (SP2+RP2-SR)
    
    if deltaL0>deltaLi:
        mainEdge = 1 #'first' edge is primary
    else:
        mainEdge = 0 #'second' edge is primary
    return mainEdge
        

def H_switch(vi, v2):
    if vi>=v2:
        H = 1
    else:
        H = 0    
    return H


def dist2D(vtxi, vtx2):
    return np.sqrt((vtxi[0]-vtx2[0])**2 + (vtxi[1]-vtx2[1])**2)
    
    
def Dwedge_nord2000(fr, betav, thetaS, thetaR, tau, tauS, tauR, el, Rs, Rr, ZGs, ZGr):
    pdiffr = pdiffr_nord2000(fr, betav, thetaS, thetaR, tau, tauS, tauR, el, Rs, Rr, ZGs, ZGr)
    omiga = 2*np.pi*fr
    Df = pdiffr * el/np.exp(1j*omiga*tau)
    return Df


def pdiffr_nord2000(fr, betav, thetaS, thetaR, tau, tauS, tauR, el, Rs, Rr, ZGs, ZGr):
    ''' This is the diffraction function for a single barrier.
        fr is frequency
        betav is the outside angle of the barrier
        thetaS is the angle from the right face of the barrer to the conneting line from the source to the diffraction edge
        thetaR is the angle from the right face of the barrier to the conneing line from the receiver to the diffraction edge
        tauS is the travel time from teh source (S) to the top of the wedge (T)
        tauR is the travel time from T to the receiver (R)
        tau = tauS + tauR
        Rs is the distance from S to T
        Rr is the distance from T to R
        ZGs is the surface impedance of the source wedge face
        ZGr is the suface imepedance of the receiver wedge face
    '''
    omiga = 2*np.pi*fr;
    nv = np.pi/betav;
    thv = thetan(thetaS, thetaR, betav)

    prodQAE = 0
    for n in [1, 2, 3, 4]:
        if n==1:
            Q = 1
        elif n==2:
            Q = refl_coef_sphericalWave(fr, tauS+tauR, min(betav-thetaS, np.pi*0.5), ZGs)
        elif n==3:
            Q = refl_coef_sphericalWave(fr, tauS+tauR, min(thetaR, np.pi*0.5), ZGr)
        elif n==4:
            Qs = refl_coef_sphericalWave(fr, tauS+tauR, min(betav-thetaS, np.pi*0.5), ZGs)
            Qr = refl_coef_sphericalWave(fr, tauS+tauR, min(thetaR, np.pi*0.5), ZGr)
            Q = Qs*Qr
            
        else:
            print('error of calculating Q')        
                    
        AThetan = A_thetan(thv[n-1], betav, nv)
        Bvalue = B_value(thetaS, thetaR, tauS, tauR, tau, omiga, AThetan, nv);
        EnvHat = E_nv_hat(Bvalue, thv[n-1], betav, nv, tauS, tauR, tau)
        prodQAE = prodQAE + Q * AThetan * EnvHat
            
    pdiffr = -1/np.pi * prodQAE * np.exp(1j*omiga*tau)/el
    return pdiffr


def A_D_hat(B):
    AD = np.sign(B)*(f(B) - 1j*g(B))
    return AD


def E_nv_hat(B, thetanv, betav, nv, tauS, tauR, tau):
    AD = A_D_hat(B)
    ATheta = A_thetan(thetanv, betav, nv)
    EnvHat = np.pi/np.sqrt(2)*np.sin(abs(ATheta))/(abs(ATheta)) * \
            np.exp(1j*np.pi*0.25)*AD/np.sqrt(1+(2*tauS*tauR/tau**2+0.5) * \
            np.cos(abs(ATheta))**2/nv**2)    
    return EnvHat


def B_value(thetaS, thetaR, tauS, tauR, tau, omiga, AThetan, nv):
    return np.sqrt(4*omiga*tauS*tauR/(np.pi*tau)) * np.cos(abs(AThetan))/np.sqrt(nv**2+(2*tauS*tauR/tau**2 + 0.5)*np.cos(AThetan)**2)


def A_thetan(thetanv, betav, nv):
    if np.pi-thetanv>=0:
        H = 1
    else:
        H = 0
    ATheta = nv*0.5*(-betav-np.pi+thetanv)+np.pi*H
    return ATheta
    

def thetan (thetaS, thetaR, betav):
    thv = np.zeros(4)
    thv[0] = thetaS - thetaR
    thv[1] = thetaS + thetaR
    thv[2] = 2*betav - (thetaS + thetaR)
    thv[3] = 2*betav - (thetaS - thetaR)
    return thv


def refl_coef_sphericalWave(fr, tau2, grazingAngle, ZG):
    Rp = refl_coef_planeWave(fr, grazingAngle, ZG)
    omiga = 2*np.pi*fr
    rhoHat = rho_hat(omiga, grazingAngle, ZG, tau2)
    EHat = E_hat(rhoHat)
    Q = Rp+(1-Rp)*EHat
    return Q


def E_hat(rhoHat):
    return 1+1j*np.sqrt(np.pi)*rhoHat*np.exp(-rhoHat**2)* erfc(-1j*rhoHat);


def rho_hat(omiga, grazingAngle, ZG, tau2):
    return (1+1j)/2*np.sqrt(omiga*tau2)*(np.sin(grazingAngle) + 1/ZG);


def refl_coef_planeWave(fr, grazingAngle, ZG):
    try:
        Rp = (np.sin(grazingAngle) - 1/ZG)/(np.sin(grazingAngle) + 1.0/ZG)
    except:
        print('zero division encountered! return 1')
        Rp = 1
    return Rp
    

def Z_G(fr, flowResistivity):
    return 1+9.08*(1000*fr/flowResistivity)**(-0.75) + 11.9*(1000*fr/flowResistivity)**(-0.73)*1j


if __name__=='__main__':
    # an example
    speed = 340.0

    betai = 1.5*np.pi
    thetaiS = 1.5*np.pi - math.atan(4.8/10)
    thetaiR = 0
    beta2 = 1.5*np.pi
    theta2S = 1.5*np.pi
    theta2R = math.atan(4.5/5.7)
    
    grazingAngle = 1./6*np.pi
    flowResistivity = 20000 #// impedance class G in page 34. resignate normal alsphalt and concrete
    pref = 2.*10**(-5)
    
    spos = [-4.8,0]
    rpos = [14.5, 4.3]
    diffEdgePosi = [0, 10]
    diffEdgePos2 = [10, 10]
    
    Rs = np.sqrt(4.8**2+10**2.)
    Rm = diffEdgePos2[0] - diffEdgePosi[0]
    Rr = np.sqrt(4.5**2+5.7**2)
    el = Rs+Rr+Rm
    tauS = Rs/speed
    tauM = Rm/speed
    tauR = Rr/speed
    tau = tauS+tauM+tauR
    
    freqs = [63., 125., 250., 500., 1000., 2000., 4000., 8000.]
    pfs = []
    for n in range(len(freqs)):
        fr = freqs[n]
        print fr, ' Hz'
        ZG = Z_G(fr, flowResistivity)
        Zis = 200000.
        Z2r = 200000.
        temp = p2edge(fr, betai, thetaiS, thetaiR, beta2, theta2S, theta2R, \
            tauS, tauM, tauR, Rs, Rm, Rr, Zis, Z2r, \
            spos, rpos, diffEdgePosi, diffEdgePos2)
        pfs.append(abs(temp/(np.exp(1j*2.*np.pi*fr*tau)/el)))
     
    print pfs
    Lpreff = 10.*np.log10(abs(np.asarray(pfs))**2.)
    print(Lpreff)
    plt.semilogx(freqs, Lpreff)
    plt.xticks(freqs, freqs)
    plt.grid()
    plt.xlim([50,15000])
