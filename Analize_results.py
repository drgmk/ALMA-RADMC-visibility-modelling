import numpy as np
import math as ma
import matplotlib.pyplot as plt
#import seaborn as sns
import emcee
import corner
from matplotlib import rc

font= {'family':'Times New Roman', 'size': 16}
rc('font', **font)
# rc('text', usetex=True)

Nruns=1
bad_runs=[]
Nbad=len(bad_runs)
#Nruns=Nruns-len(bad_runs)
Nit=1500
ndim, nwalkers=14, 70 #24
burn_in=300

par_labels=[r'$M_\mathrm{dust}$',
            r'$r_{\min}$',
            r'$r_{\max}$',
            r'$\gamma_1$',
            r'$\gamma_2$',
            r'$r_\mathrm{gap}$',
            r'$w_\mathrm{gap}$',
            r'$\delta_\mathrm{gap}$',
            r'$h$',
            r'$i$',
            r'PA',
            r'$\delta_\mathrm{RA}$',
            r'$\delta_\mathrm{Dec}$',
            r'$f_{\sigma}$']

chains=np.zeros(((Nruns-Nbad)*(Nit-burn_in)*nwalkers,ndim))
chains2=np.zeros((Nit-burn_in, (Nruns-Nbad)*nwalkers,ndim))

lnpchain=np.zeros(((Nruns-Nbad)*(Nit-burn_in)*nwalkers))
lnpchain2=np.zeros(((Nit-burn_in), (Nruns-Nbad)*nwalkers))

i2=0
for i in xrange(Nruns):
    print i
    if (i+1) not in bad_runs:
        
        path='./'#'./Runi_'+str(i+1)+'/'
        samplesi=np.load(path+'samples.dat.npy')                                                    
        lnprobsi=np.load(path+'lnprobs.dat.npy')

        #print np.shape(samplesi)

        #samples2=samples[:,burn_in:,:].
        #print np.shape(samplesi[:,burn_in:,:].reshape((nwalkers*(Nit-burn_in), ndim),order='C'))
        #print (Nit-burn_in)*nwalkers
        chains[(i-i2)*(Nit-burn_in)*nwalkers:(i-i2+1)*(Nit-burn_in)*nwalkers,:]=samplesi[:,burn_in:,:].reshape((nwalkers*(Nit-burn_in), ndim),order='C')
        lnpchain[(i-i2)*(Nit-burn_in)*nwalkers:((i-i2)+1)*(Nit-burn_in)*nwalkers]=lnprobsi[:,burn_in:].reshape((nwalkers*(Nit-burn_in)),order='C')

        for j in xrange(nwalkers):
            chains2[:,(i-i2)*nwalkers+j,:]=samplesi[j,burn_in:,:].reshape((Nit-burn_in, ndim),order='C')
            lnpchain2[:, (i-i2)*nwalkers+j]=lnprobsi[j,burn_in:].reshape(((Nit-burn_in)),order='C')

    else:
        i2+=1


# print "best chi_red = %1.5e" %(-2.0*np.max(lnpchain)/Nvis)        

fig=plt.figure(figsize=(10,8))

for i in xrange((Nruns-Nbad)*nwalkers):

    for ip in xrange(ndim):
        ax=fig.add_subplot(ndim+1,1,ip+1)
        ax.plot(chains2[:,i,ip])
        ax.set_ylabel(par_labels[ip])

    ax=fig.add_subplot(ndim+1,1,ndim+1)
    ax.plot(lnpchain2[:,i])
    ax.set_ylabel('ln(P)')
plt.show()

for i in xrange(ndim):
    print par_labels[i]+' mean and std = %1.4f , %1.4f' %(np.mean(chains[:,i]), np.std(chains[:,i]))


fig=corner.corner(chains, 
                  labels=par_labels,
                  quantiles=[0.16, 0.5,0.84], 
                  bins=20 ,
                  levels=[0.68, 0.95, 0.997],
                  show_titles=True,
                  # title_fmt=".2f",
                  title_kwards={"fontsize": 12})
                  # smooth=1.0,

plt.savefig('corner_all.pdf')


plt.show()

