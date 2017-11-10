import numpy as np
import os, sys
import matplotlib.pyplot as plt 
#from matplotlib.colors import LogNorm
import astropy.io.fits as pyfits
import functions_radmc as fradmc


############### CONSTANTS
au=1.496e13           # [cm]
pc=3.08567758e18      # [cm]
M_sun=1.9891e33       # [g]
M_earth=5.97219e27    # [g]
R_sun=6.955e+10       # [cm]
G=6.67259e-8          # cgs
sig=5.6704e-5         # cgs


radmc_inputs=float(sys.argv[1]) # if 1 grid, wavelength grid,
                                # opacities and stellar spectrum are produced

##############################################
###################### STAR AND DISK PARAMETERS 
##############################################
R_star=1.0*R_sun    
T_star=5750.0  # [K] , multiple of 250.0 
M_star=1.0*M_sun      #  [g]
dir_stellar_template = '/data/sm2132/RADMC/Stellar_Spectra/'
 

M_dust = float(sys.argv[2])*M_earth  # total dust mass

grain_density=2.7    # g/cm3 (important as it is used to calculate opacities)
Sizedist_exp=-3.5
amin=0.0001 # [cm] minimum grain size
amax=1.0 # [cm] maximum grain size
Nspec=1 # number of dust species separated logarithmically in size
        # (each dust specie is sampled with 100 different sizes)
        # (relevant if fitting SED as small grains can be superheated)


###  DISK parameters

rmin=float(sys.argv[3])
rmax=float(sys.argv[4])
gamma_1=float(sys.argv[5])
rgap=float(sys.argv[6])
wgap=float(sys.argv[7]) # fwhm
dgap=float(sys.argv[8]) # peak 

flare_odisk=1.0 # 1.0 == no flaring
R_odisk = 100.0 # [AU]
H_odisk = float(sys.argv[9])*R_odisk

gamma_2=float(sys.argv[10])

sgap=wgap/(2.0*np.sqrt(2.0*np.log(2.0)))

##############################################
######################  DEFINING GRID of model
##############################################

Nr=100  # odd edges of the cells
Nth=20  # edges of one emisphere for the theta sampling inlcudding equator 
Nphi=2  # cells in phi (azimuthal angle)

Rmax=200.0     # [au] upper edge of the R sampling       
Rmin=20.0      # [au] 
 
h_c=0.1
Thmax=np.arctan(5.0*h_c) # maximum elevation angle in the grid [rad] (last edge)


#------Inner part of the disk

Redge=np.zeros(Nr) #from Rmin to Rmax
R=np.zeros(Nr-1)
Redge[0]=Rmin

######## linear sampling
Redge=np.linspace(Rmin,Rmax,Nr)
dR=Redge[1]-Redge[0]
R=Redge[1:]-dR/2.0

####### logarithmic (better for ppdisks when more than one order of
####### magnitude involve)

# Px_r=(Rmax/Rmin)**(1.0/(Nr-1))
# #print "Rin"
# for i in xrange(1,Nr):
#     Redge[i]=Rmin*Px_r**i   
#     R[i-1]=( Redge[i]+Redge[i-1])/2.0
    
#### logatithmic sampling in theta
# Thc=0.03 # minimum value of h
# Thmin=Thc/5.0  # second edge

# Pth=(Thmax/Thmin)**(1.0/(Nth-2)) 

# Thedge=np.zeros(Nth) #from Rmin to Rmax
# Th=np.zeros(Nth-1)
# Thedge[0]=0.0
# Thedge[1]=Thmin
# Th[0]=Thmin/2.0

# for i in xrange(1,Nth-1):
#     Thedge[i+1]=Thedge[i]*Pth
#     Th[i]=(Thedge[i+1]+Thedge[i])/2.0

### linear sampling in theta

Thedge=np.linspace(0.0,Thmax,Nth)
dth=Thedge[1]-Thedge[0]
Th=Thedge[:-1]+dth/2.0
Thmin=Th[0]  # second edge

#### linear sampling in phi
dphi=2*np.pi/Nphi
Phiedge=np.arange(0.0,2.0*np.pi,dphi) #
Phi=Phiedge+dphi/2.0


#### write grid in a file for RADMC

if radmc_inputs:

    path='amr_grid.inp' #'amr_grid.inp'

    arch=open(path,'w')

    arch.write('1 \n') # iformat: The format number, at present 1
    arch.write('0 \n') # Grid style (regular = 0)
    arch.write('101 \n') # coordsystem: If 100<=coordsystem<200 the
                         # coordinate system is spherical
    arch.write('0 \n') # gridinfo
    arch.write('1 \t 1 \t 1 \n') # incl x, incl y, incl z
    arch.write(str(Nr-1)+ '\t'+ str((Nth-1)*2)+'\t'+ str(Nphi)+'\n') 
    # arch.write('0 \t 0 \t 0 \n') # levelmax, nleafsmax, nbranchmax 

    for i in range(Nr):
        arch.write(str(Redge[i]*au)+'\t')
    arch.write('\n')

    for i in range(Nth):
        arch.write(str(np.pi/2.0-Thedge[Nth-1-i])+'\t')  # from northpole to equator
    for i in range(1,Nth):
        arch.write(str(np.pi/2.0+Thedge[i])+'\t')       # from 0 to -pi/2
    arch.write('\n')


    for i in range(Nphi):
        arch.write(str(Phiedge[i])+'\t')
    arch.write(str(2.0*np.pi)+'\t')
    arch.write('\n')

    arch.close()
    # ################################################
    # :::::::::::::::::::::::WAVELENGTHS IN THE MODEL
    # ###############################################
    print 'Setting Wavelengths...'
    wmin=0.09 #[um]
    wmax=10000.0 #[um]
    Nw=150
    Pw=(wmax/wmin)**(1.0/(Nw-1))

    waves=np.zeros(Nw)
    waves[0]=wmin
    for i in xrange(1,Nw):
        waves[i]=wmin*Pw**i

    # ----- write wavelength_micron.inp

    path='wavelength_micron.inp'
    arch=open(path,'w')
    arch.write(str(Nw)+'\n')
    for i in xrange(Nw):
        arch.write(str(waves[i])+'\n')
    arch.close()


# #########################################
#----- Density function of the dusty disk
# #########################################

def Sigma_p(r, rmin, rmax, gamma1, gamma2, rgap, wgap, dgap):
    
    # rmin
    # rmax
    # gamma
    # rgap
    # wgap
    # dgap
    
    if r>=rmin and r<=rmax:
    
        Sigma=(r/80.0)**gamma1
    
        # if abs(r-rgap)<=wgap/2.0:
        Sigma=Sigma*(1.0-dgap*np.exp(    - (r-rgap)**2.0/(2.0*sgap**2)      ))
        
        return Sigma
    elif r>rmax and r<Rmax:

        return ((rmax/80.0)**gamma1)*(r/rmax)**gamma2

    else: return 0.0

def rho_p(r,z, rmin, rmax, gamma1,gamm2, rgap, wgap, dgap):
    
    Hs=H_odisk*(r/R_odisk)**flare_odisk

    return Sigma_p(r, rmin, rmax, gamma1,gamm2, rgap, wgap, dgap)*np.exp(-z**2.0/(2.0*Hs**2.0))/(np.sqrt(2.0*np.pi)*Hs*au)

"""
### plot functions above to check
rs=np.linspace(20.0,200.0,100)
rhos=np.zeros(100)
for i in xrange(100):
    rhos[i]=Sigma_p(rs[i], rmin, rmax, gamma_1, gamma_2, rgap, wgap, dgap) #rho_p(rs[i],0.0,r0,sigr)
   
plt.plot(rs,rhos)
plt.show()
""" 

##################################
######## define grain size grid
##################################
As=np.zeros(Nspec) # mean grain size of a bin
Asedge=np.zeros(Nspec+1) # bin edges
Ms=np.zeros(Nspec) # dust mass in each bin


Asedge[0]=amin
Asedge[Nspec]=amax

Pa=(amax/amin)**(1.0/(Nspec)) # edges log spaced

for i in xrange(0,Nspec):
    Asedge[i+1]=amin*Pa**(i+1)
    As[i]=(Asedge[i]+Asedge[i+1])/2.0

for i in xrange(Nspec):
        Ms[i] = abs(Asedge[i+1]**(Sizedist_exp+3.0+1.0)-Asedge[i]**(Sizedist_exp+3.0+1.0)) 

if M_dust>0.0: Ms=Ms*M_dust/np.sum(Ms) # normalise to set total dust mass = M_dust

rho_d=np.zeros((Nspec,(Nth-1)*2,Nphi,Nr-1)) # density field

M_dust_temp= 0.0 #np.zeros(Nspec) 


# for ia in xrange(len(As)):
# print ia
# print "Dust species = ", ia
  
for k in xrange(Nth-1):
    theta=Th[Nth-2-k]
    for i in xrange(Nr-1):
        rho=R[i]*np.cos(theta)
        z=R[i]*np.sin(theta)
        # if rho>Rmin_out:

        # for j in xrange(Nphi):

        rho_d[:,k,:,i]=rho_p(rho,z, rmin, rmax, gamma_1,gamma_2, rgap, wgap, dgap)
        rho_d[:,2*(Nth-1)-1-k,:,i]=rho_d[:,k,:,i]

        M_dust_temp+=2.0*rho_d[0,k,0,i]*2.0*np.pi*rho*(Redge[i+1]-Redge[i])*(Thedge[Nth-2-k+1]-Thedge[Nth-2-k])*R[i]*au**3.0 

for ia in xrange(Nspec):
    rho_d[ia,:,:,:]=rho_d[ia,:,:,:]*Ms[ia]/M_dust_temp
        
        
#-------Save it in a file for RADMC

path='dust_density.inp'

dust_d=open(path,'w')

dust_d.write('1 \n') # iformat  
dust_d.write(str((Nr-1)*2*(Nth-1)*(Nphi))+' \n') # iformat n cells 
dust_d.write(str(Nspec)+' \n') # n species

for ai in xrange(Nspec):
    for j in range(Nphi):
        for k in range(2*(Nth-1)):
            for i in range(Nr-1):
                dust_d.write(str(rho_d[ai,k,j,i])+' \n')


dust_d.close()


#::::::::::::::::: DUST OPACITIES
if radmc_inputs:
    

    print "Opacities..."

    path='dustopac.inp'
    path_opac='./'


    arch=open(path,'w')
    arch.write("2               Format number of this file \n")
    arch.write(str(Nspec)+"              Nr of dust species \n")
    arch.write("============================================================================ \n")

    for i in xrange(Nspec):
        print i
        fradmc.opac_arg(Asedge[i]*1.0e4, Asedge[i+1]*1.0e4, grain_density, 100, 1,'opct_dustgrains', 'dust_'+str(i+1) )
        
        os.system('cp '+path_opac+'Tempkappa/* ./')

        # print i,Noliv-1
        arch.write("1               Way in which this dust species is read \n")
        arch.write("0               0=Thermal grain \n")
        arch.write("dust_"+str(i+1)+ " Extension of name of dustkappa_***.inp file \n")
        arch.write("---------------------------------------------------------------------------- \n")
    arch.close()
    

    # -------- STAR and wavelength


    def intflux(lami,flux,lam): # function to interpolate or extrapolate
        l=len(lam)
        if lami<lam[0]:
            m=0.0#(flux[1]-flux[0])/(lam[1]-lam[0])
            return m*(lami-lam[0])+flux[0]

        elif lami>lam[l-1]:
            m=0.0#(flux[l-1]-flux[l-2])/(lam[l-1]-lam[l-2])
            return m*(lami-lam[l-1])+flux[l-1]

        for j in xrange(l-1):
            if lami>=lam[j] and lam[i]<=lam[j+1]:
                m=(flux[j+1]-flux[j])/(lam[j+1]-lam[j])
                return m*(lami-lam[j])+flux[j]  
    
    print 'Setting Stellar flux...'

    def BBody(lam,T): # black body
        # lam in cm
        h_planck=6.62607e-27 # cgs 
        k=1.381e-16          #erg/K
        c=3.0e10             # cm/s
    
        return 2.0*h*c/(lam**3.0*np.exp(h*c/(lam*k*T))-1.0)

    # load template
    path=dir_stellar_template+'Kurucz'+str(int(T_star))+'-4.0.fits.gz'
    fit = pyfits.open(path) #abrir objeto cubo de datos
    data = fit[0].data # wave is in um, flux is in F_lam, units? #ergs cm**-2 s**-1 A**-1
    ltemp=len(data[0,:])

    tempflux=data[1,:]*(R_star/pc)**2.0 # computing flux at 1 pc
    tempwave=data[0,:]


    stellar_flux=np.zeros(Nw)

    I=0.0

    for i in range(Nw-1):
        wi=waves[i]
        if wi<tempwave[0]: stellar_flux[i]=0.0#tempflux[0]
    
        elif wi >=tempwave[0] and wi <tempwave[ltemp-6]:

            for j in range(ltemp-1):
                if wi>=tempwave[j] and wi<tempwave[j+1]:
                    m=(tempflux[j+1]-tempflux[j])/(tempwave[j+1]-tempwave[j])
                    stellar_flux[i]=m*(wi-tempwave[j])+tempflux[j]
                    break
        elif wi >=tempwave[ltemp-6]: #too sparse point to make an interpolation,  propto lam**-4
            stellar_flux[i]=tempflux[ltemp-6]*(wi/tempwave[ltemp-6])**-4.0
        dlam=(waves[i+1]-waves[i])*1.0e4 # from um to A
        I+=stellar_flux[i]*dlam

 
    I_star=sig*T_star**4.0*(R_star/pc)**2.0 
    
    stellar_flux[:]=stellar_flux[:]*I_star/I  # normalizing spectrum
                                              # so Lstar=sig*R**2*T**4

       
    # UV excess between 0.091 um and 0.25 um 

    fUV= 0.0 # fraction of uv excess
    slope=2.2  # slope for Fnu

    fluxuv=np.zeros(Nw)

    Iuv=0.0
    for i in range(Nw):
        if waves[i]>= 0.091 and waves[i]<=0.25:
        
            fluxuv[i]=(waves[i])**(slope-2.0)  # C * erg/ s cm**2 A
            dlam=(waves[i+1]-waves[i])*1.0e4 # from um to A
            Iuv+=fluxuv[i]*dlam

    fluxuv=fluxuv*I_star*fUV/Iuv # normalize uv excess


    # transform from F_lam to F_nu

    clight=2.99792458e18 # light speed in A/s
    Flux=np.zeros(Nw)
    Fluxwouv=np.zeros(Nw)
    for i in range(Nw):
        Flux[i]=(stellar_flux[i]+fluxuv[i])*(waves[i]*1.0e4)**2.0/clight
        Fluxwouv[i]=(stellar_flux[i])*(waves[i]*1.0e4)**2.0/clight

    path='stars.inp'
    arch_star=open(path,'w')
    arch_star.write('2 \n')
    arch_star.write('1  '+str(Nw)+'\n')
    arch_star.write(str(R_star)+'\t'+str(M_star)+'\t'+'0.0   0.0   0.0 \n')
    for i in xrange(Nw):
        arch_star.write(str(waves[i])+'\n')

    for i in xrange(Nw):
        arch_star.write(str(Flux[i])+'\n')

    arch_star.close()


    """
    # plot stellar spectrum
    plt.plot(waves,Flux,marker='o',color='blue')
    plt.plot(waves,Fluxwouv,marker='o',color='red')
    plt.xscale('log')
    plt.yscale('log')
    plt.show()
    """
