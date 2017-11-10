import numpy as np
import matplotlib.pyplot as plt
import math as ma
import cmath as cma
from scipy import interpolate

def Intextpol(x,y,xi):

    Nx=len(x)
    if xi<=x[0]: return y[0] # extrapol                                                                                                                                                                         
    elif xi<=x[Nx-1]: #interpol                                                                                                                                                                                                              
        for l in xrange(1,Nx):
            if xi<=x[l]:
                return y[l-1]+(xi-x[l-1])*(y[l]-y[l-1])/(x[l]-x[l-1])

    elif xi>x[Nx-1]:    #extrapol                                                                                                                                                                                                            
        alpha=ma.log(y[Nx-1]/y[Nx-2])/ma.log(x[Nx-1]/x[Nx-2])
        return y[Nx-1]*(xi/x[Nx-1])**alpha

def effnk(n1,k1,n2,k2,n3,k3,f2,f3): 

    # mixing rule Bruggeman http://en.wikipedia.org/wiki/Effective_medium_approximations
 
    np1=n1+k1*1j # matrix
    np2=n2+k2*1j # inclusion 1
    np3=n3+k3*1j # inclusion 2

    
    e1=np1**2.0  # n = sqrt(epsilon_r x mu_r) and mu_r is aprox 1
    e2=np2**2.0
    e3=np3**2.0


    # polynomial of third order
    
    p=np.zeros(4, dtype=complex)

    p[3]=e1*e2*e3*(1.0+f2+f3) # 0 order
    p[2]=-e1*e3*f2 -e1*e2*f3 - e2*e3 + 2*(e1*e2*f2 + e1*e2 +e1*e3*f3+e1*e3+e2*e3*f2+e2*e3*f3) # 1st order
    p[1]= -2.0*(e1*f2+e1*f3+e3*f2+e2*f3+e2+e3)+4.0*(e1+e2*f2+e3*f3)# 2nd order
    p[0]= -4.0*(1.0+f2+f3)

    roots=np.roots(p) 

    # check roots
    for i in xrange(len(roots)):
        effi=roots[i]
        if effi.real>0.0 and effi.imag>0.0:
            return cma.sqrt(effi)
    ### if nothins satisfy the above condition
        
    return -1.0

    #print "number of roots = ",len(roots)
    # A2=(e2-e1)/(e2+2.0*e1)
    # A3=(e3-e1)/(e3+2.0*e1)

    # A=f2*A2+f3*A3

    #print roots==e1
    # print roots[0]**3*p[0]+roots[0]**2*p[1]+roots[0]*p[2]+p[3]
    # print roots[1]**3*p[0]+roots[1]**2*p[1]+roots[1]*p[2]+p[3]
    # print roots[2]**3*p[0]+roots[2]**2*p[1]+roots[2]*p[2]+p[3]

    # print np.abs(roots-e1)/np.abs(e1)<0.01

    # eff=roots[0] #e1*(1+2.0*A)/(1-A)
    # print eff.imag
    
### there are 3 roots, we choose the first one as it is the most
    ### similar to the Maxwell Garnett solution (check
    ### http://lib.tkk.fi/Diss/2003/isbn9512260956/article3.pdf)

    return cma.sqrt(eff)

#print effnk(1.0,2.0,1.1,2.1,0.3)

path1='./astrosilicate_ext.lnk'
path2='./ac_opct.lnk'
path3='./ice_opct.lnk'

d1=4.0
d2=2.5
d3=1.0

m1=0.7 #0.5
m2=0.15 # 0.2
m3=0.15 #0.3

V1=m1/d1
V2=m2/d2
V3=m3/d3

df=(d1*V1+d2*V2+d3*V3)/(V1+V2+V3)

v1=V1/(V1+V2+V3) # volume fraction matrix
v2=V2/(V1+V2+V3) # volume fraction 1st inclusion
v3=V3/(V1+V2+V3) # volume fraction 2nd inclusion

f2=v2/v1
f3=v3/v1

print v2, v3
print "final density = ", df
pathout='./mix_ASiACIce.lnk'


O1=np.loadtxt(path1)
O2=np.loadtxt(path2)
O3=np.loadtxt(path3)


# logO1n=interpolate.interp1d(np.log10(O1[:,0]), np.log10(O1[:,1]))
# logO1k=interpolate.interp1d(np.log10(O1[:,0]), np.log10(O1[:,1]))

# logO2n=interpolate.interp1d(np.log10(O2[:,0]), np.log10(O2[:,1]))
# logO2k=interpolate.interp1d(np.log10(O2[:,0]), np.log10(O2[:,1]))

# logO3n=interpolate.interp1d(np.log10(O3[:,0]), np.log10(O3[:,1]))
# logO3k=interpolate.interp1d(np.log10(O3[:,0]), np.log10(O3[:,1]))


Nw=200

wmin=0.1
wmax=10000.0


Opct1=np.zeros((Nw,3))
Opct1[0,0]=wmin

i=0

# n1=10.0**logO1n(np.log10(Opct1[i,0]))
# n2=10.0**logO2n(np.log10(Opct1[i,0]))
# n3=10.0**logO3n(np.log10(Opct1[i,0]))

# k1=10.0**logO1k(np.log10(Opct1[i,0]))
# k2=10.0**logO2k(np.log10(Opct1[i,0]))
# k3=10.0**logO3k(np.log10(Opct1[i,0]))


n1=Intextpol(O1[:,0],O1[:,1],Opct1[i,0])
n2=Intextpol(O2[:,0],O2[:,1],Opct1[i,0])
n3=Intextpol(O3[:,0],O3[:,1],Opct1[i,0])

k1=Intextpol(O1[:,0],O1[:,2],Opct1[i,0])
k2=Intextpol(O2[:,0],O2[:,2],Opct1[i,0])
k3=Intextpol(O3[:,0],O3[:,2],Opct1[i,0])


eff=effnk(n1,k1,n2,k2,n3,k3,f2,f3)

Opct1[i,1]=eff.real
Opct1[i,2]=eff.imag

for i in xrange(1,Nw):
    #print i
    Opct1[i,0]=wmin*(wmax/wmin)**(i*1.0/(Nw-1))
    #print Opct[i,0]

    # n1=10.0**logO1n(np.log10(Opct1[i,0]))
    # n2=10.0**logO2n(np.log10(Opct1[i,0]))
    # n3=10.0**logO3n(np.log10(Opct1[i,0]))

    # k1=10.0**logO1k(np.log10(Opct1[i,0]))
    # k2=10.0**logO2k(np.log10(Opct1[i,0]))
    # k3=10.0**logO3k(np.log10(Opct1[i,0]))


    n1=Intextpol(O1[:,0],O1[:,1],Opct1[i,0])
    n2=Intextpol(O2[:,0],O2[:,1],Opct1[i,0])
    n3=Intextpol(O3[:,0],O3[:,1],Opct1[i,0])

    k1=Intextpol(O1[:,0],O1[:,2],Opct1[i,0])
    k2=Intextpol(O2[:,0],O2[:,2],Opct1[i,0])
    k3=Intextpol(O3[:,0],O3[:,2],Opct1[i,0])

    eff=effnk(n1,k1,n2,k2,n3,k3,f2,f3)

    #print eff
    Opct1[i,1]=eff.real
    Opct1[i,2]=eff.imag


# plt.plot(O1[:,0], O1[:,1], color='blue')
# plt.plot(Opct1[:,0], 10.0**logO1n(np.log10(Opct1[:,0])), color='red')
# plt.xscale('log')
# plt.yscale('log')
# plt.show()

# plt.plot(O2[:,0], O2[:,1], color='blue')
# plt.plot(Opct1[:,0], 10.0**logO2n(np.log10(Opct1[:,0])), color='red')
# plt.xscale('log')
# plt.yscale('log')
# plt.show()

# plt.plot(O3[:,0], O3[:,1], color='blue')
# plt.plot(Opct1[:,0], 10.0**logO3n(np.log10(Opct1[:,0])), color='red')
# plt.xscale('log')
# plt.yscale('log')
# plt.show()


np.savetxt(pathout,Opct1)

plt.subplot(1,2,1)
plt.plot(O1[:,0],O1[:,1], color='blue', label='Astrosilicates')
plt.plot(O2[:,0],O2[:,1], color='red', label='Amorphous carbon')
plt.plot(O3[:,0],O3[:,1], color='green', label='Ice')
plt.plot(Opct1[:,0],Opct1[:,1], color='black', label='Mix')

plt.xlabel(r'Wavelength [$\mu$m]')
plt.ylabel(r'n')
plt.xscale('log')
plt.yscale('log')
plt.legend()

plt.subplot(1,2,2)
plt.plot(O1[:,0],O1[:,2], color='blue', label='Astrosilicates')
plt.plot(O2[:,0],O2[:,2], color='red', label='Amorphous carbon')
plt.plot(O3[:,0],O3[:,2], color='green', label='Ice')
plt.plot(Opct1[:,0],Opct1[:,2], color='black', label='Mix', lw=3.0)

plt.xlabel(r'Wavelength [$\mu$m]')
plt.ylabel(r'k')
plt.xscale('log')
plt.yscale('log')
plt.legend()
plt.show()
