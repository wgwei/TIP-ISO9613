# -*- coding: utf-8 -*-
"""
Created on Fri Dec 12 18:24:19 2014
main function to compare the ISO9613-2 and the pierce_diffr model
@author: Weigang Wei
"""
import numpy as np
import matplotlib.pylab as plt
from iso9613_2 import calc_Dz, calc_Dz2, dist2D
from line_x_poly import line_x_poly
from pierce_diffr_rigid import LdiffOverLff_theo, LdiffOverLff_simple
import math

class Building(object):
    def __init__(self, vertices, height=10.):
        ''' define the vertices of a building in a plane view
            vetices = [[x1, y1], [x2, y2], [x3, y3],..., [x1,y1]]'''
            
        self.vtx = np.asarray(vertices)[:,0]
        self.vty = np.asarray(vertices)[:,1]
        self.hi = height

class CNOSSOS_spectrum():
    def __init__(self):
        ''' spectrum condition: 70km/hour light vehicle
            traffic flow information: 40000 AAWT. average annual weekday traffic
            the Qm = 40000./24 = 1666
        '''
        Qm = 833.
        vm = 70.
        self.freq = [63., 125., 250., 500., 1000., 2000., 4000., 8000]
        self.spectrum = np.asarray([101.09, 97.17, 99.54, 99.76, 101.05, 96.57, 89.99, 84.75 ]) + 10*np.log10(Qm/1000./vm)
        self.Aweight =  np.array([-26.2, -16.1,-8.6, -3.2, 0, 1.2, 1, -1.1])   # {63:-26.2, 125:-16.1, 250:-8.6, 500:-3.2, 1000:0, 2000:1.2, 4000:1, 8000:-1.1}
        self.specA = self.spectrum+self.Aweight
    def plot_spectrum(self):
        plt.semilogx(self.freq, self.specA, '-', self.freq, self.spectrum, '--')
        plt.legend(['A-weighted', 'Linear'], loc=0)
        plt.xlabel('Frequency Hz')
        plt.ylabel('Sound power level dB or dB(A)')
        plt.xlim([50, 10000])
        plt.text(100, 73, '70km/hour, light vehicle, CNOSSOS emission model')
        plt.text(100, 70, '20000 average annual weekday traffic')
        plt.savefig('LW_spectrum.png')
      
        
def cmp_refl_iso_pierce_rigid_A():
    ''' calculate the multiple reflections in a 2D sections '''
    # case 1: in 2D section    
    cnos = CNOSSOS_spectrum()
    SPEED = 340.0        
    reflCoef = 0.8
    REFLNUM = 20
    SHEIGHT = 0.05
    COLOURS = ['r', 'b', 'c', 'm', 'k', 'g', 'y', 'r', 'b', 'c', 'm', 'k', 'g', 'y']
    HEIGHTS = [5., 6., 7., 8., 9., 10., 15]
    offsets = [5., 10., 20., 50.]
    for OFFSET in offsets:
        for BW in [10.]: #[5., 10.,15., 20.]: # building width
            for CW in [10.]: #[5., 10., 15., 20.]:  # canyon width
                print ' \n\n CW %0f' %CW      
                plt.figure()
                for h, height in enumerate(HEIGHTS):
                    reflN = []
                    totalLiso = []
                    totalLPrc = []
                    rposvp = [BW+OFFSET, 4.]
                    sposhp = [-0.5*CW, 100.]
                    rposhp = [BW+OFFSET, 100.]
                    bobj = Building([[0,0], [BW,0], [BW, 200], [0, 200], [0,0]], height) # building starts from 0, 0
                               
                    LAiso, LAPierce = 0.0, 0.0            
                    for n in range(REFLNUM):     
                        print '\n'
                        Lwf = 0.0
                        for c, fr in enumerate(cnos.freq):                    
                            Lwf = cnos.specA[c] + n*10*np.log10(reflCoef)+ 10. # 10 dB is used to convert the sound power level from dB/m to 10m segement
                            print 'refl %d    fr %0f Hz,    LW %0f dB' %(n, fr, Lwf)
                            sposvp = [-0.5*CW-n*CW, SHEIGHT]
                            interSect = line_x_poly(sposhp[0], sposhp[1], rposhp[0], rposhp[1],bobj.vtx, bobj.vty)
                            interSect = sorted(interSect)
                            barWidth = dist2D(interSect[0], interSect[1])
        #                    print 'interSect :', interSect
                            
                            # ISO9613-2
                            a = abs(interSect[0][1] - interSect[1][1])
                            wavelength = SPEED/fr
                            if n==0:
                                groundIncluded = 1
                            else:
                                groundIncluded = 0
                            barpos = [[interSect[0][0], bobj.hi], [interSect[1][0], bobj.hi]]
                            d = dist2D(sposvp, rposvp)
                            dz = calc_Dz2(sposvp, rposvp, barWidth, barpos, a, wavelength, groundIncluded)
                            LAiso = 10.*np.log10(10.**(0.1*LAiso)+10.**(0.1*(Lwf - dz - 20.*np.log10(d) - 11)))                
            
                            # pierce                
                            thetas = math.atan(abs(sposvp[0]/bobj.hi))
                            thetar = math.atan(abs((rposvp[0] - BW)/(bobj.hi-rposvp[1])))
                            rs = dist2D(sposvp, [interSect[0][0], bobj.hi])
                            rr = dist2D(rposvp, [interSect[1][0], bobj.hi])
                            print 'thetas: %0.3f thetar: %0.3f  rs: %0.2f  rr: %0.2f, BW: %0f waveL: %0.1f' %(thetas, thetar, rs, rr, BW, SPEED/fr)
                            
                            levelRefFF = LdiffOverLff_theo(rs, rr, thetas, thetar, BW, SPEED/fr)
                            print 'levelrefFF ', levelRefFF
                            if n==0:
                                # plus 3dB to include the ground reflection. the oter part corrects the path diffence
                                IL = - levelRefFF - 3. + 10*np.log10(dist2D(sposvp, rposvp)/(rs+rr+dist2D(interSect[0], interSect[1]))) 
                            else:
                                IL = - levelRefFF + 10*np.log10(dist2D(sposvp, rposvp)/(rs+rr+dist2D(interSect[0], interSect[1])))
                            print 'dz, IL: ', dz, IL
                            LAPierce = 10.*np.log10(10.**(0.1*LAPierce)+10.**(0.1*(Lwf - IL - 20.*np.log10(d) - 11)))
                        totalLiso.append(LAiso)
                        totalLPrc.append(LAPierce)                
                        reflN.append(n)
                    plt.plot(reflN, totalLiso, '-', color = COLOURS[h])
                    plt.plot(reflN, totalLPrc, '--', color = COLOURS[h])
                    
                plt.legend(['ISO h=5', 'Pierce h=5', 'ISO h=6', 'Pierce h=6', 'ISO h=7', 'Pierce h=7',\
                            'ISO h=8', 'Pierce h=8', 'ISO h=9', 'Pierce h=9', 'ISO h=10', 'Pierce h=10', \
                            'ISO h=15', 'Pierce h=15'], loc=4)
                plt.xticks(range(REFLNUM))
                plt.grid()
                plt.title('OFFSET %s CW %s BW %s' %(str(OFFSET), str(CW), str(BW)))
                plt.xlabel('Number of reflectionis')
                plt.ylabel('Sound pressure level [dB]')
                plt.savefig('OFFSET' + str(OFFSET) + '_CW'+str(CW)+'_BW'+str(BW)+'_A.png')
            

def cmp_refl_iso_pierce_rigid_A2():
    ''' calculate the multiple reflections in a 3D '''
    # case 1: in 2D section    
    cnos = CNOSSOS_spectrum()
    SPEED = 340.0        
    reflCoef = 0.8
    REFLNUM = 20
    COLOURS = ['r', 'b', 'c', 'm', 'k', 'g', 'y', 'r', 'b', 'c', 'm', 'k', 'g', 'y']
    HEIGHTS = [5., 6., 7., 8., 9., 10., 15]
    bobj = Building([[0,0], [10,0], [10, 200], [0, 200], [0,0]]) # building starts from 0, 0
    sposhpY = np.linspace(0, 200, 21) # position of sources/cars
    segementLen = sposhpY[1] - sposhpY[0]
    BW = 10. # building width
    OFFSETD = 10.
    RHEIGHT = 4.
    SHEIGHT = 0.05   
    
    for OFFSETD in [5., 10., 20., 50.]:
        rposhp = [BW+OFFSETD, 100.] # the receiver is supposed in the center of the building
        for CW in [10.]:#[5., 10., 15., 20.]:  # canyon width
            plt.figure()
            for h, height in enumerate(HEIGHTS):
                LAiso, LAPierce = 0.0, 0.0  
                print 'height -> ', height
                reflN = []
                totalLiso = []
                totalLPrc = []     
                for n in range(REFLNUM): 
                    # contribution of all sources in 200 m
                    for s,sv in enumerate(sposhpY):
                                   
                        bobj = Building([[0,0], [BW,0], [BW, 200], [0, 200], [0,0]], height) # building starts from 0, 0
                          
                        sposhp = [-0.5*CW - n*CW, sv]
                        a = abs(sv - rposhp[1])
                        interSect = line_x_poly(sposhp[0], sposhp[1], rposhp[0], rposhp[1], bobj.vtx, bobj.vty)
                        interSect = sorted(interSect)
        #                print 'interSect -> ', interSect
                        barWidth = dist2D(interSect[0], interSect[1])
                        sposvp = [-dist2D(rposhp, interSect[0]), SHEIGHT]            
                        rposvp = [barWidth+dist2D(rposhp, interSect[1]), RHEIGHT]
                        
                        Lwf = 0.0
                        # contribution of all frequencies
                        for c, fr in enumerate(cnos.freq):                    
                            Lwf = cnos.specA[c] + n*10*np.log10(reflCoef) + 10.*np.log10(segementLen) #  10.*np.log10(segementLen) is used to correct the sourec power level from dB/10m to dB/segment length
        #                    print 'refl %d    fr %0f Hz,    LW %0f dB' %(n, fr, Lwf)                
                            
                            # ISO9613-2
                            a = abs(sv - rposhp[1])
                            wavelength = SPEED/fr
                            if n==0:
                                groundIncluded = 1
                            else:
                                groundIncluded = 0
                            barpos = [[interSect[0][0], bobj.hi], [interSect[1][0], bobj.hi]]
                            
                            dss = np.sqrt(sposhp[0]**2.+(bobj.hi-SHEIGHT)**2.)
                            dsr = np.sqrt(OFFSETD**2.+(bobj.hi-RHEIGHT)**2.)
                            d = np.sqrt((abs(sposhp[0])+BW+OFFSETD)**2.+(abs(RHEIGHT-SHEIGHT))**2.)
                            
                            dz = calc_Dz(dss, dsr, BW, d, a, wavelength, groundIncluded)
                            print 'dz -> ', dz
                            LAiso = 10.*np.log10(10.**(0.1*LAiso)+10.**(0.1*(Lwf - dz - 20.*np.log10(d) - 11)))                
            
                            # pierce                
                            thetas = math.atan(abs(sposhp[0]/(bobj.hi-SHEIGHT)))
                            thetar = math.atan(abs((rposhp[0]-BW)/(bobj.hi-RHEIGHT)))
                            rsv = dist2D(sposvp, barpos[0])
                            rrv = dist2D(rposvp, barpos[1])
        #                    print 'thetas: %0.3f thetar: %0.3f  rs: %0.2f  rr: %0.2f, BW: %0f waveL: %0.1f' %(thetas, thetar, rs, rr, BW, SPEED/fr)
                            
                            levelRefFF = LdiffOverLff_theo(dss, dsr, thetas, thetar, BW, SPEED/fr)
                                      
                            if n==0:
                                # plus 3dB to include the ground reflection. the oter part corrects the path diffence
                                IL = -levelRefFF - 3. - 10*np.log10((dss+dsr+BW)/(rsv+rrv+barWidth)) 
                            else:
                                IL = -levelRefFF - 10*np.log10((dss+dsr+BW)/(rsv+rrv+barWidth)) 
                            print 'IL -> ', IL
                            LAPierce = 10.*np.log10(10.**(0.1*LAPierce)+10.**(0.1*(Lwf - IL - 20.*np.log10(d) - 11)))
                    totalLiso.append(LAiso)
                    totalLPrc.append(LAPierce)                
                    reflN.append(n)
                print 'totalLiso -> ', totalLiso
                print 'totalLPrc -> ', totalLPrc
                plt.plot(reflN, totalLiso, '-', color = COLOURS[h])
                plt.plot(reflN, totalLPrc, '--', color = COLOURS[h])
                
            plt.legend(['ISO h=5', 'Pierce h=5', 'ISO h=6', 'Pierce h=6', 'ISO h=7', 'Pierce h=7',\
                            'ISO h=8', 'Pierce h=8', 'ISO h=9', 'Pierce h=9', 'ISO h=10', 'Pierce h=10', \
                            'ISO h=15', 'Pierce h=15'], loc=4)
            plt.xticks(range(REFLNUM))
            plt.grid()
            plt.xlabel('Number of reflectionis')
            plt.ylabel('Sound pressure level [dB]')
            plt.title('OFFSET %s CW %s BW %s' %(str(OFFSETD), str(CW), str(BW)))
            plt.savefig('OFFSET' + str(OFFSETD) + '_CW'+str(CW)+'_BW'+str(BW)+'_A_line200m.png')
        

if __name__=='__main__':
#    plt.figure()
#    cnosObj = CNOSSOS_spectrum()
#    cnosObj.plot_spectrum()
#    cmp_refl_iso_pierce_rigid()
#    cmp_refl_iso_pierce_rigid_A()
    cmp_refl_iso_pierce_rigid_A2()
    plt.show()
            
        
