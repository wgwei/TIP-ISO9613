# -*- coding: utf-8 -*-
"""
Created on Tue Jun  9 09:46:03 2015
    Ground absorption Agr in ISO9613-2
    Aim: 1) to figure out the attenuation in a soft ground
         2) to give suggestions when using Cadna/A
@author: Weigang Wei
"""

import numpy as np
import matplotlib.pylab as plt

class SRInfo():
    ''' Information of source and receiver '''
    def __init__(self, hs, hr, dp, Gs, Gr, Gm):
        ''' hs = height of source
            hr = height of receiver
            dp = distance from source to receiver
            G = ground factor from 0 to 1. Hard ground G=0, porose ground G=1
        '''
        self.hs = hs
        self.hr = hr
        self.dp = dp
        self.Gs = Gs
        self.Gr = Gr
        self.Gm = Gm
    
def a_h(h, dp):
    return 1.5 + 3.0*np.exp(-0.12*(h-5.)**2.) * (1 - np.exp(-dp/50.)) + 5.7*np.exp(-0.09*h**2.)*(1-np.exp(-2.8*(1e-6)*dp**2.)) 

def b_h(h, dp):
    return 1.5 + 8.6*np.exp(-0.09*h**2.)*(1-np.exp(-dp/50.))

def c_h(h, dp):
    return 1.5 + 14.0*np.exp(-0.46*h**2.)*(1-np.exp(-dp/50.))

def d_h(h, dp):
    return 1.5 + 5.0*np.exp(-0.9*h**2.)*(1-np.exp(-dp/50.))

def q_value(hs, hr, dp):
    if dp<=30.*(hs + hr):
        return 0.0
    else:
        return 1 - 30.*(hs + hr)/dp
        
def Agr_by_freq(fr, hs, hr, dp, Gs, Gm, Gr):
    if fr==63:
        As = -1.5
        Ar = -1.5
        Am = -3 * q_value(hs, hr, dp)
        Agr = As + Ar + Am
    elif  fr==125:
        As = -1.5 + Gs*a_h(hs, dp)
        Ar = -1.5 + Gr*a_h(hr, dp)
        Am = -3. * q_value(hs, hr, dp) * (1-Gm)
        Agr = As + Ar + Am
    elif fr==250:
        As = -1.5 + Gs*b_h(hs, dp)
        Ar = -1.5 + Gr*b_h(hr, dp)
        Am = -3. * q_value(hs, hr, dp) * (1-Gm)
        Agr = As + Ar + Am
    elif fr==500:
        As = -1.5 + Gs*c_h(hs, dp)
        Ar = -1.5 + Gr*c_h(hr, dp)
        Am = -3. * q_value(hs, hr, dp) * (1-Gm)
        Agr = As + Ar + Am
    elif fr==1000:
        As = -1.5 + Gs*d_h(hs, dp)
        Ar = -1.5 + Gr*d_h(hr, dp)
        Am = -3. * q_value(hs, hr, dp) * (1-Gm)
        Agr = As + Ar + Am
    elif fr==2000:
        As = -1.5*(1-Gs)
        Ar = -1.5*(1-Gr)
        Am = -3. * q_value(hs, hr, dp) * (1-Gm)
        Agr = As + Ar + Am
    elif fr==4000:
        As = -1.5*(1-Gs)
        Ar = -1.5*(1-Gr)
        Am = -3. * q_value(hs, hr, dp) * (1-Gm)
        Agr = As + Ar + Am
    elif fr==8000:
        As = -1.5*(1-Gs)
        Ar = -1.5*(1-Gr)
        Am = -3. * q_value(hs, hr, dp) * (1-Gm)
        Agr = As + Ar + Am
    else:
        print('Only octanve bands from 63 to 8000 Hz are avalaible in ISO9613')
    
    return Agr
    
def plot_abcd_h():
    dps = np.linspace(20, 2000, 2000)
    # plot a_h    
    plt.figure()
    for i, h in enumerate([1.5, 3.0, 6.0, 7.5, 10.]):
        ah = []
        for dp in dps:
            ah.append(a_h(h, dp))
        plt.semilogx(dps, ah)
        plt.xlim([20, 2000])
        plt.ylim([0, 10])
        plt.xlabel('Distance dp, m')
        plt.ylabel("a', dB")
    
    # plot b_h
    plt.figure()
    for i, h in enumerate([1.5, 2., 2.5, 3., 3.5, 4., 5.0, 10.]):
        bh = []
        for dp in dps:
            bh.append(b_h(h, dp))
        plt.semilogx(dps, bh)
        plt.xlim([20, 2000])
        plt.ylim([0, 10])
        plt.xlabel('Distance dp, m')
        plt.ylabel("b', dB")
    
    #plot c_h
    plt.figure()
    for i, h in enumerate([1.5, 1.75, 2.0, 2.5, 3.0]):
        ch = []
        for dp in dps:
            ch.append(c_h(h, dp))
        plt.semilogx(dps, ch)
        plt.xlim([20, 2000])
        plt.ylim([0, 10])
        plt.xlabel('Distance dp, m')
        plt.ylabel("c', dB")
    
    #plot d_h
    plt.figure()
    for i, h in enumerate([1.5, 3.0]):
        dh = []
        for dp in dps:
            dh.append(d_h(h, dp))
        plt.semilogx(dps, dh)
        plt.xlim([20, 2000])
        plt.ylim([0, 10])
        plt.xlabel('Distance dp, m')
        plt.ylabel("d', dB")          
                
    
def plot_As_by_hs_500Hz():
    ''' hss  = changes of source height
    '''
    hss = np.linspace(0, 5, 50)
    symbols = ['k-', '--', '-+', '-^', '-*', '-o', 'r-']
    Gs = 1    
    
    # plot c'(h)
    dps = [10, 50, 100, 200, 500, 1000]
    plt.figure()
    for i, dp in enumerate(dps):
        chList = []
        for h in hss:
            chList.append(c_h(h, dp))
        
        As = -1.5 + Gs * np.asarray(chList)
        plt.plot(hss, As, symbols[i])
        plt.ylabel("As [dB]")
        plt.xlabel("h [m]")
    plt.legend(['dp='+str(d)+'m' for d in dps])
    plt.text(np.max(hss/2.), np.max(As)/2., "G=1")
    
def plot_As_by_dp_500Hz():
    ''' dps  = changes of source receiver distance
    '''  
    dps = np.linspace(5, 500, 1000)
    symbols = ['k-', '--', 'm--', 'r--', 'c-', 'm-', 'r-']
    Gs = 1
    hss = [0, 0.2, 1, 1.5, 3]
    plt.figure()
    for i, h in enumerate(hss):
        chList2 = []
        for dp in dps:
            chList2.append(c_h(h, dp))        
        
        As = -1.5 + Gs * np.asarray(chList2)
        plt.plot(dps, As, symbols[i])
        plt.ylabel("As [dB]")
        plt.xlabel("dp [m]") 
    plt.legend(['h='+str(hh)+'m' for hh in hss])   
    plt.text(np.max(dps/2.), 5., "G=1")
    
def plot_Am_by_h_500Hz():
    hss = np.linspace(0, 5, 50)
    dps = [10, 50, 100, 200, 500, 1000]
    symbols = ['k-', '--', '-+', '-^', '-*', '-o', 'r-']
    hr = 1.5
    for Gm in [0, 0.2, 0.5, 0.8,1]:
        plt.figure()
        for i, dp in enumerate(dps):
            Am = []
            for hs in hss:
                Am.append(-3. * q_value(hs, hr, dp) * (1 - Gm))
            plt.plot(hss, Am, symbols[i])
            plt.ylabel("Am")
            plt.xlabel("hs [m]")
        plt.legend(['dp='+str(d)+'m' for d in dps], loc='best')
        plt.ylim([-3, 0])
        plt.title('Gm = ' + str(Gm))
    
def plot_A_along_propa_path_500Hz():
    plt.figure()
    Gs = 1.
    hs = 0.2
    Gm = Gs
    initGr = Gs
    hr = 1.5
    initDps = [10, 50, 100, 200, 500]
    d2edge = np.linspace(0, 30*hr+10, int(30*hr+10)+1)
    for dpi in initDps:
        print('dp -> ', dpi, '\n')
        AgrPlusAdivj = []
        #backward from the edge to source
        steps = np.linspace(1., dpi, int(dpi))
        for sp in steps:
            Gr = Gs
            As = -1.5 + Gs * c_h(hs, sp)
            Am = -3. * q_value(hs, hr, sp) * (1 - Gm)
            Ar = -1.5 + Gr * c_h(hr, sp)
            Adiv = 20.*np.log10(sp) + 11
            AgrPlusAdivj.append(As+Am+Ar+Adiv)
            
        # forwad from the edge to 30*hr
        AgrPlusAdiv = []
        for d2e in d2edge:
            if d2e<30*hr:
                Gr = ((30*hr-d2e)*initGr + d2e*0) / (30*hr) # d2e*0: suppose the extended part is rigid.
            else:
                Gr = 0.                
            dp = dpi + d2e
            
            As = -1.5 + Gs * c_h(hs, dp)
            Am = -3. * q_value(hs, hr, dp) * (1 - Gm)
            Ar = -1.5 + Gr * c_h(hr, dp)
            Adiv = 20.*np.log10(dp) + 11
            print('As: ', As, 'Ar: ', Ar, 'Adiv: ', Adiv)
            
            AgrPlusAdiv.append(As+Am+Ar+Adiv)            
        plt.plot(list(steps) + list(dpi+d2edge), AgrPlusAdivj + AgrPlusAdiv)
    plt.legend([str(d)+'m' for d in initDps], loc=4)
    plt.xlabel('Source receiver distance [m]')
    plt.ylabel('Agr + Adiv [dB]')
    plt.grid()
        
def plot_Agr_singleNum_eq10():
    heights =[0., 0.5, 1.5, 3., 10] 
    dps = np.linspace(10., 500., 50)
    symbols = ['k-', '--', '-+', '-^', '-*', '-o', 'r-']
    plt.figure()
    for i, hm in enumerate(heights):
        Agr = []
        for d in dps:
            agr = 4.8 - (2*hm/d) * (17 + 300./d)
            if agr>=0:
                Agr.append(agr)
            else:
                Agr.append(0.0)
        plt.plot(dps, Agr, symbols[i])
    plt.legend(['h='+str(v)+'m' for v in heights], loc=4)    
    plt.xlabel('d [m]')
    plt.ylabel('hm [m]')


def plot_Agr_spec():
    hes = np.asarray(range(10)) + 1
    dps = np.linspace(20, 2000, 2000)
#    hs = 0.05
    hr = 1.5
    Gs, Gm, Gr = 1., 1., 1.
    
    # plot every Agr,f
    for fr in [125, 250, 500, 1000, 2000, 4000]:
        plt.figure(fr)
        for i, h in enumerate(hes):
            Agrf = []
            for dp in dps:                
                Agrf.append(Agr_by_freq(fr, h, hr, dp, Gs, Gm, Gr))            
            plt.semilogx(dps, Agrf)
            plt.ylabel("As [dB]")
            plt.xlabel("dp [m]") 
            plt.xlim([min(dps), max(dps)])
            plt.ylim([0, 16])
        plt.legend(['h='+str(hh)+'m' for hh in hes], loc='best')   
        plt.text(np.max(dps/2.), np.max(Agrf)/2., "G=1")
        plt.grid()
        
    # for typical traffic source height
    plt.figure("Traffic")
    hs = 0.05
    hr = 1.5
    Gs, Gm, Gr = 1., 1., 1.
    trafficSpec = 100. + np.asarray([-14.5, -10.2, -7.2, -3.9, -6.4, -11.4]) #125, 250, 500, 1k, 2k, 4k
    LwTotal = 10*np.log10(np.sum([10**(0.1*Lv) for Lv in trafficSpec]))
    
    AgrSpec = []
    AgrSingle = []
    for dp in dps:
        Lvs = []
        for f, fr in enumerate([125, 250, 500, 1000, 2000, 4000]):   
            Agrf = Agr_by_freq(fr, hs, hr, dp, Gs, Gm, Gr)
            Lvs.append(trafficSpec[f] - Agrf)
        LwChange = 10*np.log10(np.sum([10**(0.1*Lv) for Lv in Lvs]))
        AgrSpec.append(LwTotal - LwChange)
        
        agr = 4.8 - (2*hs/dp) * (17 + 300./dp)
        if agr>=0:
            AgrSingle.append(agr)
        else:
            AgrSingle.append(0.0)
    
    plt.semilogx(dps, AgrSingle, '--')
    plt.semilogx(dps, AgrSpec, '-')
    
    plt.xlim([min(dps), max(dps)])
    plt.xlabel('dp [m]')
    plt.ylabel('Agr [dB]')
    plt.legend(['by 7.3.1', 'by 7.3.2'], loc=4)
    plt.grid()
    
    # for a typical plant noise
    plt.figure("Plant")
    hs = 2.
    hr = 1.5
    Gs, Gm, Gr = 1., 1., 1.
    
    steamTurbineHall = np.asarray([102., 90., 74., 64., 59., 59.]) + np.asarray([-16.1,  -8.6,  -3.2,   0,   1.2,   1.0]) #125, 250, 500, 1k, 2k, 4k
    LwTotal = 10*np.log10(np.sum([10**(0.1*Lv) for Lv in steamTurbineHall]))
    AgrSpec = []
    AgrSingle = []
    for dp in dps:
        Lvs = []
        for f, fr in enumerate([125, 250, 500, 1000, 2000, 4000]):
            Agrf = Agr_by_freq(fr, hs, hr, dp, Gs, Gm, Gr)
            Lvs.append(steamTurbineHall[f] - Agrf)
        LwChange = 10*np.log10(np.sum([10**(0.1*Lv) for Lv in Lvs]))
        AgrSpec.append(LwTotal - LwChange)
        
        agr = 4.8 - (2*hs/dp) * (17 + 300./dp)
        if agr>=0:
            AgrSingle.append(agr)
        else:
            AgrSingle.append(0.0)
    
    plt.semilogx(dps, AgrSingle, '--')        
    plt.plot(dps, AgrSpec, '-')
    
    plt.xlim([min(dps), max(dps)])
    plt.xlabel('dp [m]')
    plt.ylabel('Agr [dB]')
    plt.legend(['by 7.3.1', 'by 7.3.2'], loc=4)
    plt.grid()
    
if __name__=='__main__':
    # Test    
    plot_abcd_h()
    plot_As_by_hs_500Hz()
    plot_As_by_dp_500Hz()
    plot_Am_by_h_500Hz()
    plot_A_along_propa_path_500Hz()
    plot_Agr_singleNum_eq10()
    plot_Agr_spec()
    plt.show()
#
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
