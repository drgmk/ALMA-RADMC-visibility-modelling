import os
import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as pyfits
import math as ma
from functions import dered


#os.system('radmc3d mctherm')
#os.system('python waves.py')
#os.system('radmc3d spectrum loadlambda incl 28.5 phi 0.0 posang 11.0  secondorder')

def BBody(lam,T):
    #lam in cm
    h_planck=6.62607e-27 # cgs 
    k=1.381e-16          #erg/K
    c=3.0e10             # cm/s
    
    return 2.0*h_planck*c/(lam**3.0*(ma.exp(h_planck*c/(lam*k*T))-1.0))

R_sun=6.955e+10 # cm
R_star=3.5*R_sun#3.0*R_sun
T_star=8250
pc=3.085677e18 # cm


#Av=0.22
#Rv=3.1
cum=2.99792458e14

# read sed
"""
arch = open('./HD-181327.txt','r')

nlines = len(arch.readlines())
arch.seek(0)
Np=0
waves=[]
flux=[]
eflux=[]

for l in xrange(nlines):
    line=arch.readline()
    dat=line.split()
    if dat[0]!='#':
        if dat[1]!='-' and dat[2]=='1' and dat[3]!='-':
            Np+=1
arch.seek(0)
SED=np.zeros((Np,3))
i=0
for l in xrange(nlines):
    line=arch.readline()
    dat=line.split()
    if dat[0]!='#':
        if dat[1]!='-' and dat[2]=='1' and dat[3]!='-':
            SED[i,0]=float(dat[1])
            SED[i,1]=float(dat[3])
            SED[i,2]=float(dat[4])
            i+=1
            
"""
SED=np.loadtxt('./HD-181327-2.txt')
# jansky to watts/mu
SED[:,1]=SED[:,1]*1.0e-26*cum/SED[:,0]**2.0
SED[:,2]=SED[:,2]*1.0e-26*cum/SED[:,0]**2.0

#SED=SEDred
#SED[:,1]=dered(SED[:,0],SEDred[:,1],Av,Rv)

c=2.99792458e18 # light speed in A/s
cms=2.99792458e8 # light speed in m/s 

#---------READ spectrum.out

arch=open('spectrum.out','r')
arch.readline()
Nw=int(arch.readline())
print Nw
SED1=np.zeros((Nw,2))
arch.readline()
for i in xrange(Nw):
    line=arch.readline()
    dat=line.split()
    SED1[i,0]=dat[0]
    SED1[i,1]=dat[1]

arch.close()

SED1[:,1]=SED1[:,1]*c/((SED1[:,0]*1.0e4)**2.0) # in erg/s A cm**2 
SED1[:,1]=SED1[:,1]*1.0e4 # in erg/s um cm**2.0

d=51.8 #pc
C=1.0e-3/(d**2.0) # erg/cm^2 to joule/m^2 and at 140 pc



R_BB=150.0*R_star
T_BB=75.0
Nwi=100
BB=np.zeros((Nwi,2))
wmax=5000.0
wmin=10.0


for i in xrange(len(BB[:,0])):
    wave=wmin*(wmax/wmin)**(i*1.0/(Nwi-1))
    BB[i,0]=wave
    BB[i,1]=BBody(wave*1.0e-4,T_BB)*ma.pi*(R_BB/(d*pc))**2.0 # Fnu in cgs

BB[:,1]=BB[:,1]*1.0e-3*cum/((BB[:,0])**2.0) # Flam in mks



"""
#------------- Herschel

Fh=np.loadtxt('./HD142527_herschelspectrum.dat') # GHz/Jansky

sedHer=Fh

sedHer[:,0]=cms*1.0e6/(sedHer[:,0]*1.0e9) # GHz to wavelength in microns
sedHer[:,1]=Fh[:,1]*1.0e-26*cms*1.0e6/(sedHer[:,0])**2.0         # Jansky to W/m2/um

#print sedHer
# unred

sedHer[:,1]=dered(sedHer[:,0],sedHer[:,1],Av,Rv)

sedHer[:,1]=sedHer[:,1]*sedHer[:,0] # flud density to sed erg/s/cm2 
"""
#-------------  PLOT

plt.figure(figsize=(12,6))


plt.plot(SED1[:,0],SED1[:,1]*C*SED1[:,0],color='red',linestyle='-', label='Dust Trap model')
#plt.plot(SED2[:,0],SED2[:,1]*C*SED2[:,0],color='green',marker='o',linestyle='--',lw=2.0, label='Dust Trap model 30')
#plt.plot(SED3[:,0],SED3[:,1]*C*SED3[:,0],color='orange',marker='o',linestyle='--',lw=2.0, label='Dust Trap 3, wo big grains')
#plt.plot(BB[:,0],BB[:,1]*BB[:,0],color='black',lw=2.0, label='Black body')

plt.errorbar(SED[:,0],SED[:,1]*SED[:,0],yerr=(SED[:,2]*SED[:,0]*1.0),color='blue',fmt='o',label='Deredened Observations')
#plt.plot(sedHer[:,0],sedHer[:,1],color='grey',label='Herschel spectrum')

plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'Wavelength [$\mu m$]')
plt.ylabel('$\lambda F_{\lambda}$ [ $W$ $m^{-2}$]')
plt.ylim(1.0e-18,1.0e-13)
plt.xlim(0.1,1.0e4)
plt.legend(numpoints=1)#,loc=3)
plt.savefig('sed.pdf')

plt.show()

