import sys, os, copy
import Chi2_simvis_c
import numpy as np
import functions_radmc as fradmc
from time import gmtime, strftime

print 'packages loaded'
Tin=strftime("%Y-%m-%d %H:%M:%S", gmtime())
print Tin
#----create density of disc model, opacities, star, etc.

dpc=27.5 # the intencisty is invariant of distance   
X0= 1.847762146765e+02  # [deg] RA star
Y0= 1.654758391174e+01  # [deg] DEC star
 
Mdust=0.38 # earth masses
rmin=32.4
rmax=135.8
gamma1=0.42 # sigma propto r**gamma
gamma2=-11.0 # sigma propto r**gamma

rgap=77.0
wgap=39.0
dgap=0.38

h=0.05#3      # vertical aspect ratio

PA=154.0
inc=19.0
offx=0.0 # RA offset, positive values to the left (east)
offy=0.0 # DEC offset

fsigma=2.3287


 
os.system('python setup_mcmc.py 0 %1.5f %1.5f %1.5f %1.5f %1.5f %1.5f %1.5f %1.5f %1.5f' %(Mdust,rmin,rmax, gamma1, rgap, wgap, dgap, h, gamma2))

# """
#---- define parameters for RADMC run 
arch=open('radmc3d.inp','w')
arch.write('nphot =       1000000 \n')
arch.write('nphot_scat=   100000 \n')
arch.write('nphot_spec=   10000 \n')
arch.write('nphot_mono=   100000 \n')
arch.write('scattering_mode_max = 0\n')
arch.write('modified_random_walk = 0\n')
arch.write('istar_sphere=0 \n')
arch.write('tgas_eq_tdust=1 \n')
arch.write('incl_lines = 0 \n')
arch.write('setthreads = 6')
arch.close()

##### compute temperature field
# os.system('radmc3d mctherm')

images=['ALMAB6']
wavelengths=[1300.0] # [um] wavelength
fields=[0]

Npix=512
dpix=0.05

canvas = 'canvas.b6.'+str(Npix)+'pixels'+str(dpix)+'arcsec.field'

print 'making images'
#### make image 

fradmc.Simimages_canvas_fields(dpc, 
                               X0, Y0, 
                               images, 
                               wavelengths,
                               fields, 
                               Npix, dpix, 
                               canvas, 
                               inc, PA, 
                               offx, offy, 
                               pb=1.0, 
                               tag='best_pb')

# fradmc.Simimages_canvas_fields(dpc, 
#                                X0, Y0, 
#                                images, 
#                                wavelengths,
#                                fields, 
#                                Npix, dpix, 
#                                canvas, 
#                                inc, PA, 
#                                offx, offy, 
#                                pb=0.0, 
#                                tag='best')

# os.system('python Simimages_canvas.py best2_pb '+canvas+' %1.5f %1.5f %1.5f %1.5f %1.5f' %(inc, PA-90.0, offx, offy, 1.0)) # the last 1.0 is there to multiply for the primary beam
# os.system('python Simimages_canvas.py best2_nopb '+canvas+' %1.5f %1.5f %1.5f %1.5f %1.5f' %(inc, PA-90.0, offx, offy, 0.0)) # the last 1.0 is there to multiply for the primary beam

##### COMPUTE CHI2 AND produce table with residuals

print 'loading table with visibilities'
table=np.load('tablevis.band6.field0.npy')
print 'loaded'

print 'computing simulated visibilities and chi2'
Chi2,  Nvis_noflags = Chi2_simvis_c.Model_simvisC_bilinear('images/image_ALMAB6_best_pb_field0.fits', table)
print 'done'

os.system('mv ./tablevismodel.npy ./Simulate_obs/')

print Nvis_noflags
print 'chired = ', Chi2/(Nvis_noflags*fsigma**2)

#### convolve image with a beam of 0.5"x0.5"

# os.system('python Convolve_ALMA.py test 1.39e-4 1.39e-4  0.0')

Tout=strftime("%Y-%m-%d %H:%M:%S", gmtime())
print Tin
print Tout

# """

