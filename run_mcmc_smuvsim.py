import os, sys
import numpy as np
import matplotlib.pyplot as plt
from random import shuffle
from time import gmtime, strftime
import uuid
import commands
import pickle
import emcee

import functions_radmc as fradmc
import Chi2_simvis_c


#import pickle


Tin=strftime("%Y-%m-%d %H:%M:%S", gmtime())

print Tin


ndim, nwalkers =14 , 70 # 8,100  # number of free parameters, walkers
Nit = 1500 # number of iterations
emceethreads=70 # number of threads (cores) used by emcee 

maindir=commands.getoutput('pwd')+'/'


############# STELLAR PARAMETERS

dpc=27.5 # the intencisty is invariant of distance   
X0= 1.847762146765e+02  # [deg] RA star
Y0= 1.654758391174e+01  # [deg] DEC star

########### IMAGE PARAMETERS

images=['ALMAB6']
wavelengths=[1300.0] # [um] wavelength
Npix=512
dpix=0.05
canvas = 'canvas.b6.512pixels0.05arcsec.field'


############ DATA PARAMETERS

fields=[0] # field indices in case of a mosaic (if only one leave as 0)
table_obs='tablevis.band6.field'  # 0.npy'
table_vis=[]
for i in xrange(len(fields)):
    table_vis.append(np.load(table_obs+str(fields[i])+'.npy'))



############## SEEDS 

Mdust0=0.38 # earth masses
rmin0=32.0
rmax0=135.0
gamma10=0.4 # sigma propto r**gamma
gamma20=-10.0 # sigma propto r**gamma
rgap0=77.0
wgap0=40.0
dgap0=0.38
h0=0.1      # vertical aspect ratio

PA0=155.0
inc0=19.0
offx0=0.0 # RA offset, positive values to the left (east)
offy0=0.0 # DEC offset

fsigma0=2.3287 # rescaling for new visibility errors

# deltas

sig_Mdust=0.02#0.33 # earth masses
sig_rmin=2.0
sig_rmax=2.0
sig_gamma1=0.05 # sigma propto r**gamma
sig_gamma2=1.0 # sigma propto r**gamma
sig_rgap=5.0
sig_wgap=5.0
sig_dgap=0.05
sig_h=0.03      # vertical aspect ratio
sig_PA=3.0
sig_inc=1.0
sig_offx=0.01 # RA offset, positive values to the left (east)
sig_offy=0.01 # DEC offset
sig_fsigma=0.002


#### INITIAL BALL
p0=np.zeros((nwalkers,ndim)) #array([4.3 , 10.0, 0.1])
for i in xrange(nwalkers):
    p0[i,0]=np.random.normal(Mdust0,sig_Mdust)
    p0[i,1]=np.random.normal(rmin0,sig_rmin)
    p0[i,2]=np.random.normal(rmax0,sig_rmax)
    p0[i,3]=np.random.normal(gamma10,sig_gamma1)
    p0[i,4]=np.random.normal(gamma20,sig_gamma2)
    p0[i,5]=np.random.normal(rgap0,sig_rgap)
    p0[i,6]=np.random.normal(wgap0,sig_wgap)
    p0[i,7]=np.random.normal(dgap0,sig_dgap)
    p0[i,8]=np.random.normal(h0,sig_h)
    p0[i,9]=np.random.normal(inc0,sig_inc)
    p0[i,10]=np.random.normal(PA0,sig_PA)
    p0[i,11]=np.random.normal(offx0,sig_offx)
    p0[i,12]=np.random.normal(offy0,sig_offy)
    p0[i,13]=np.random.normal(fsigma0,sig_fsigma)

def Simobs(dir_model, fsigma): 
    # ### This function computes the model visibilities given a model
    # ### image and calculates the chi**2 using the routine
    # ### Chi2_simvis developed by Sebastian Marino 
    
    chi2=0.0
    Nvis=0.0
    for fi in fields:

        image='images/image_'+images[0]+'_mcmc_field'+str(fi)+'.fits'
        table=table_vis[0]

        chi2i, Nvisi =  Chi2_simvis_c.Chi2_simvisC_bilinear(image, table)
        # chi2i, Nvisi =  Chi2_simvis(image, table)

        # N=vis2.size
        Nvis+=Nvisi
        chi2+=chi2i/(fsigma**2) #np.sum((np.absolute(vis2)**2.0)*wts0)*f_chi

    chi2_red=chi2/Nvis


    print "chi2: %1.8e N: %1.8e chi2_red: %1.8e fsigma: %1.2e" %(chi2, Nvis, chi2_red, fsigma)

    os.chdir(maindir)
    os.system('rm -rf '+dir_model)
    return -chi2/2.0-Nvis*np.log(fsigma) # return ln( EXP(-chi2/2) ) for emcee


def lnprob(par):
    # This function calculates the probability of having
    # parameters=par. To do so, it simulates a model using these
    # parameters in a directory with a random name that is deleted
    # afterwards

    Mdusti, rmini, rmaxi, gamma1i, gamma2i, rgapi, wgapi, dgapi, hi, inci, PAi, offxi, offyi, fsigmai  =par # curret values
    # of free parameters impose priors. In thi case, probability =0
    # (log(P)=-inf) outside domain 

    # These restrictions have to be consistent with the limits of the
    # model, e.g. limits of the spatial grid, or domain
    if wgapi<0.0 or hi > 0.3 or hi<0.03 or Mdusti<0.0 or rmini<10.0 or rmaxi>170.0 or rmini>rmaxi or inci<0.0 or inci>90.0 or PAi<0.0 or PAi>180.0 or gamma1i>5.0 or gamma1i<-5.0 or fsigmai<0.0 or gamma2i>0.0 or rgapi+wgapi/2>rmaxi or dgapi>1.0 or dgapi<0.0:
        print 'out of bounds'
        return -np.inf

    os.chdir(maindir)
    # create directory
    rand_name=str(uuid.uuid4())
    path_dir=maindir+'model_'+rand_name
    os.system('mkdir '+path_dir)
    os.system('mkdir '+path_dir+'/images')
    os.system('cp ./*.inp '+path_dir)
    os.system('cp setup_mcmc.py '+path_dir)
    os.system('cp dust_temperature.dat '+path_dir)
    os.system('cp '+canvas+'*.fits '+path_dir)
    os.system('cp functions_radmc.py '+path_dir)

    # os.system('cp Simimages_canvas.py '+path_dir)
    # os.system('cp functions.py '+path_dir)
    # os.system('cp '+uvfits_obs+ '*.uvfits '+path_dir)
    
    
    os.chdir(path_dir)
    print 'setup...'
    # setup 
    os.system('python setup_mcmc.py 0 %1.5f %1.5f %1.5f %1.5f %1.5f %1.5f %1.5f %1.5f %1.5f' %(Mdusti,rmini,rmaxi, gamma1i, rgapi, wgapi, dgapi, hi, gamma2i)) 
   
    # radmc3d.inp
    arch=open('radmc3d.inp','w')
    arch.write('nphot =       500000 \n')
    arch.write('nphot_scat=   100 \n')
    arch.write('nphot_spec=   10000 \n')
    arch.write('nphot_mono=   100000 \n')
    arch.write('scattering_mode_max = 0\n')
    arch.write('modified_random_walk = 0\n')
    arch.write('istar_sphere=0 \n')
    arch.write('tgas_eq_tdust=1 \n')
    arch.write('incl_lines = 0 \n')
    arch.write('setthreads = 1')
    arch.close()

    # print 'mctherm...'
    # # mctherm not necesary because we are not changing the stellar
    # # parameters nor the grid

    # os.system('radmc3d mctherm > mctherm.log')

    print 'simimage...'

    # ### make image 
    fradmc.Simimages_canvas_fields(dpc, 
                                   X0, Y0, 
                                   images, 
                                   wavelengths,
                                   fields, 
                                   Npix, dpix, 
                                   canvas, 
                                   inci, PAi, 
                                   offxi, offyi, 
                                   pb=1.0, 
                                   tag='mcmc')

    # # image
    # os.system('python Simimages_canvas.py mcmc '+canvas+' %1.5f %1.5f %1.5f %1.5f %1.1f > simimgaes.log' %(inci,PAi-90.0,offxi,offyi, 0) )

    return Simobs(path_dir, fsigmai)





######## CREATE SAMPLER

sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, threads=emceethreads)

####### RUN MCMC

pos, prob, state= sampler.run_mcmc(p0, Nit)


Tout=strftime("%Y-%m-%d %H:%M:%S", gmtime())
print Tin
print Tout

print "############################## FINISH"
samples = sampler.chain  # chain= array(nwalkers,nit,ndim)

lnprobs = sampler.lnprobability



######### SAVE SAMPLES
np.save('samples.dat',samples)
np.save('lnprobs.dat',lnprobs)
print("Mean acceptance fraction: {0:.3f} "  .format(np.mean(sampler.acceptance_fraction)))
f=open('Acceptance.dat', 'w')
f.write(str(Tin)+' \n')
f.write(str(Tout)+' \n')
f.write("Nit = "+str(Nit)+' \n')
f.write("nwalkers = "+str(nwalkers)+' \n')
f.write("ndim = "+str(ndim)+' \n')
f.write("Mean acceptance fraction: {0:.3f}"  .format(np.mean(sampler.acceptance_fraction)) +' \n')
f.close()

sampler.pool=[]
pickle.dump(sampler, open('./sampler_pickle', 'wb'))


try: 

    autocorr=sampler.get_autocorr_time(c=1, low=1)

    fautocorr=open('autocorrelation.dat', 'w')
    for item in autocorr:
        print>>fautocorr, item
    fautocorr.close()

    print autocorr

except:

    print 'no autocorrelation time'
