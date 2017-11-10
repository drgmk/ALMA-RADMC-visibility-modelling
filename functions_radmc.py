import numpy as np
import sys
import os
from scipy.ndimage.interpolation import shift
from astropy.io import fits



####### OPACITIES

def opac_arg(amin, amax, density, N, Inte, lnk_file, Type):

    # amin      # microns
    # amax       # microns
    # density      # g/cm3
    # N             # bins
    # Inte # if to  average or not
    # lnk_file # file with optical constants (it assumes it ends with .lnk)
    # Type # output tag (dustkappa_'+Type+'.inp')


    exp=3.5

    path="./"
    pathout='Tempkappa/dustkappa_'+Type+'.inp'


    # ---------------------- MAIN 
    os.system('rm '+path+'Tempkappa/*')

    Pa=(amax/amin)**(1.0/(N-1.0))

    A=np.zeros(N)
    A[0]=amin


    for i in range(N):
        os.system('rm '+path+'param.inp')
        A[i]=amin*(Pa**(i))
        acm=A[i]*10.0**(-4.0)
        # print "a = %1.2e [um]"  %A[i]
        file_inp=open(path+'param.inp','w')
        file_inp.write(lnk_file+'\n')
        e=round(np.log10(acm))
        b=acm/(10.0**e)
        file_inp.write('%1.2fd%i \n' %(b,e))
        file_inp.write('%1.2f \n' %density) 
        file_inp.write('1')
        
        file_inp.close()

        os.system(path+'makeopac')
        os.system('mv '+path+'dustkappa_'+lnk_file+'.inp '+path+'Tempkappa/dustkappa_'+Type+'_'+str(i+1)+'.inp ')
    
   

    # --------- READ OPACITIES AND COMPUTE MEAN OPACITY

    if Inte==1:
        # read number of wavelengths
        opct=np.loadtxt(path+lnk_file+'.lnk')
        Nw=len(opct[:,0])


        Op=np.zeros((Nw,4)) # wl, kappa_abs, kappa_scat, g 
        Op[:,0]=opct[:,0]

        Ws_mass=np.zeros(N) # weigths by mass and abundances
        Ws_number=np.zeros(N) # wights by abundances

        for i in xrange(N):
            Ws_mass[i]=(A[i]**(-exp))*(A[i]**(3.0))*A[i]  # w(a) propto n(a)*m(a)*da and da propto a
            Ws_number[i]=A[i]**(-exp)*A[i]

        W_mass=Ws_mass/np.sum(Ws_mass)
        W_number=Ws_number/np.sum(Ws_number)

        for i in xrange(N):
            file_inp=open(path+'Tempkappa/dustkappa_'+Type+'_'+str(i+1)+'.inp','r')
            file_inp.readline()
            file_inp.readline()


            for j in xrange(Nw):
                line=file_inp.readline()
                dat=line.split()
                kabs=float(dat[1])
                kscat=float(dat[2])
                g=float(dat[3])

                Op[j,1]+=kabs*W_mass[i]
                Op[j,2]+=kscat*W_mass[i]
                Op[j,3]+=g*W_mass[i] 

            file_inp.close()
            os.system('rm '+path+'Tempkappa/dustkappa_'+Type+'_'+str(i+1)+'.inp')


        final=open(path+pathout,'w')
    
        final.write('3 \n')
        final.write(str(Nw)+'\n')
        for i in xrange(Nw):
            final.write('%f \t %f \t %f \t %f\n' %(Op[i,0],Op[i,1],Op[i,2],Op[i,3]))
        final.close()





def dered(lam,f,Av,Rv):  # unredenning according to cardelli's law


    wl=lam # microns

    # definir a y b

    npts = len(wl)
    a = np.zeros(npts)                
    b = np.zeros(npts)
    F_a=np.zeros(npts)
    F_b=np.zeros(npts)

    for i in range(npts):
        x=1.0/wl[i]
        
        if x<1.1:#x>=0.3 and x<1.1:
            a[i]=0.574*x**(1.61)
            b[i]=-0.527*x**(1.61)

        if x>=1.1 and x<3.3:
            y=x-1.82
            c1 = [ 1. , 0.104,   -0.609,    0.701,  1.137, -1.718,   -0.827,    1.647, -0.505 ]  
            c2 = [ 0.,  1.952,    2.908,   -3.989, -7.985, 11.102,    5.491,  -10.805,  3.347 ]
            for j in range(len(c1)): #polinomios
                a[i]=a[i]+c1[j]*y**j 
                b[i]=b[i]+c2[j]*y**j 

        if x>=3.3 and x<8:
            y=x
    
            if x>5.9: 
                y1=y-5.9
                F_a[i]=-0.04473 * y1**2 - 0.009779 *y1**3
                F_b[i]=0.2130 * y1**2 + 0.1207 * y1**3
    
            a[i]=1.752 - 0.316*y - (0.104 / ( (y-4.67)**2 + 0.341 )) + F_a[i]
            b[i]=-3.090 + 1.825*y + (1.206 / ( (y-4.62)**2 + 0.263 )) + F_b[i]

        if x>=8 and x<=11:

            y=x-8.0
            c1 = [ -1.073, -0.628,  0.137, -0.070 ]
            c2 = [ 13.670,  4.257, -0.420,  0.374 ]
            for j in range(len(c1)): #polinomios
                a[i]=a[i]+c1[j]*y**j 
                b[i]=b[i]+c2[j]*y**j 

    #print a
    #print b
    A=Av*(a+b/Rv)
    
    #plt.plot(wl,A,color='blue')
    #plt.plot(wl,a,color='red')
    #plt.plot(wl,b,color='green')
    #plt.show()
    Funred=np.zeros(npts)
    for i in xrange(npts):
        Funred[i]=f[i]*10.0**(0.4*A[i])
    

    return Funred


def convert_to_fits(path_image,path_fits,dpc, lam0, arcsec=False, mas=False, mx=0.0, my=0.0):

    cc  = 2.9979245800000e10      # Light speed             [cm/s]
    pc  = 3.08572e18              # Parsec                  [cm]

    lam0=float(lam0)

    reffreq=cc/(lam0*1.0e-4)

    f=open(path_image,'r')
  
    iformat=int(f.readline())

    if (iformat < 1) or (iformat > 4):
        print("ERROR: File format of {0:s} not recognized.".format(filename))
        return

    if (iformat == 1) or (iformat == 3):
        radian = (1 == 0)

    else:
        radian = (1 == 1)

    nx, ny = tuple(np.array(f.readline().split(),dtype=int))
    nf = int(f.readline()) # number of wavelengths
    sizepix_x, sizepix_y = tuple(np.array(f.readline().split(),dtype=float))

    lam = np.empty(nf)
    for i in range(nf):
        lam[i] = float(f.readline())
    
    f.readline()  

    image = np.zeros((nf,ny,nx), dtype=float)

    for k in range(nf):
        for j in range(ny):
            for i in range(nx):

                image[k,j,i] = float(f.readline())

                if (j == ny-1) and (i == nx-1):
                    f.readline()

    f.close()

    # Compute the flux in this image as seen at dpc (pc)

    
    pc = 3.0857200e18   # pc in cm

    pixdeg_x = 180.0*(sizepix_x/(dpc*pc))/np.pi
    pixdeg_y = 180.0*(sizepix_y/(dpc*pc))/np.pi

    print pixdeg_x
    # Compute the conversion factor from erg/cm^2/s/Hz/ster to erg/cm^2/s/Hz/ster at dpc
    #
    pixsurf_ster = pixdeg_x*pixdeg_y * (np.pi/180.)**2
    factor = 1e+23 * pixsurf_ster
    # And scale the image array accordingly:
    image_in_jypix = factor * image

    flux = np.sum(image_in_jypix)
   
    #### Shift image by mx and my in arcsec

    mvx_pix=(mx/(pixdeg_x*3600.0))
    mvy_pix=(my/(pixdeg_y*3600.0))
    shiftVector=(0.0, mvy_pix, -mvx_pix) # minus sign as left is positive 

    print np.shape(image_in_jypix)
    image_in_jypix_shifted=shift(image_in_jypix,shift=shiftVector, order=3)#,mode='wrap')
    #shiftx=np.roll(arrayToShift,3DshiftVector,axis=(0,1,2))

    # Make FITS header information:
    #
    header = fits.Header()
    
    #header['SIMPLE']='T'
    header['BITPIX']=-32
    # all the NAXIS are created automatically header['NAXIS']=2
    header['OBJECT']='HD109085'
    header['EPOCH']=2000.0
    header['LONPOLE']=180.0
    header['EQUINOX']=2000.0
    header['SPECSYS']='LSRK'
    header['RESTFREQ']=reffreq
    header['VELREF']=0.0
    header['CTYPE3']='FREQ'
    header['CRPIX3'] = 1.0
    header['CDELT3']  = 1.0
    header['CRVAL3']= reffreq


    header['FLUX']=flux
    print flux



    header['BTYPE'] = 'Intensity'
    header['BSCALE'] = 1
    header['BZERO'] = 0
    header['BUNIT'] = 'JY/PIXEL'#'erg/s/cm^2/Hz/ster'

    

    #header['EPOCH'] = 2000.
    #header['LONPOLE'] = 180.
    header['CTYPE1'] = 'RA---TAN'
    header['CRVAL1'] = 1.880158942527e2
    header['CTYPE2'] = 'DEC--TAN'
    header['CRVAL2'] = -1.619622628083e1
    

    # Figure out what units we are using per pixel:
    if arcsec:
        unit = 'arcsec'
        multiplier = 3600.
    elif mas:
        unit = 'mas'
        multiplier = 3.6e6
    else:
        unit = 'DEG'
        multiplier = 1
    print multiplier*pixdeg_x
    print unit
    header['CDELT1'] = -multiplier*pixdeg_x
    header['CUNIT1'] = unit
    #
    # ...Zero point of coordinate system
    #
    header['CRPIX1'] = 1.0*((nx+1)/2)

    
    header['CDELT2'] = multiplier*pixdeg_y
    header['CUNIT2'] = unit
    #
    # ...Zero point of coordinate system
    #
    header['CRPIX2'] = 1.0* ((ny+1)/2)




    if nf > 1:
        # multiple frequencies - set up the header keywords to define the
        #    third axis as frequency
        header['CTYPE3'] = 'VELOCITY'
        header['CUNIT3'] = 'km/s'
        header['CRPIX3'] = 1.0* ((nf+1)/2)
        header['CRVAL3'] = 0.0
        # Calculate the frequency step, assuming equal steps between all:
        delta_velocity = (lam[1] - lam[0])*cc*1e-5/lam0
        header['CDELT3'] = delta_velocity
        header['RESTWAVE']=lam0
    else:                # only one frequency
        header['RESTWAVE'] = lam[0]
        header['CUNIT3'] = 'Hz'
    # Make a FITS file!
    #

    image_in_jypix_float=image_in_jypix_shifted.astype(np.float32)
    print 'float'
    fits.writeto(path_fits, image_in_jypix_float, header, output_verify='fix')

def convert_to_fits_canvas(path_image,path_fits,path_canvas, dpc, lam0, arcsec=False, mas=False, mx=0.0, my=0.0 , pbm=False):
    cc  = 2.9979245800000e10      # Light speed             [cm/s]
    pc  = 3.08572e18              # Parsec                  [cm]

    lam0=float(lam0)

    reffreq=cc/(lam0*1.0e-4)

    f=open(path_image,'r')
  
    iformat=int(f.readline())

    if (iformat < 1) or (iformat > 4):
        print("ERROR: File format of {0:s} not recognized.".format(filename))
        return

    if (iformat == 1) or (iformat == 3):
        radian = (1 == 0)

    else:
        radian = (1 == 1)

    nx, ny = tuple(np.array(f.readline().split(),dtype=int))
    nf = int(f.readline()) # number of wavelengths
    sizepix_x, sizepix_y = tuple(np.array(f.readline().split(),dtype=float))

    lam = np.empty(nf)
    for i in range(nf):
        lam[i] = float(f.readline())
    
    f.readline()  

    image = np.zeros((1,nf,ny,nx), dtype=float)

    for k in range(nf):
        for j in range(ny):
            for i in range(nx):

                image[0,k,j,i] = float(f.readline())

                if (j == ny-1) and (i == nx-1):
                    f.readline()

    f.close()

    # Compute the flux in this image as seen at dpc (pc)

    
    pc = 3.0857200e18   # pc in cm

    pixdeg_x = 180.0*(sizepix_x/(dpc*pc))/np.pi
    pixdeg_y = 180.0*(sizepix_y/(dpc*pc))/np.pi

    print pixdeg_x
    # Compute the conversion factor from erg/cm^2/s/Hz/ster to erg/cm^2/s/Hz/ster at dpc
    #
    pixsurf_ster = pixdeg_x*pixdeg_y * (np.pi/180.)**2
    factor = 1e+23 * pixsurf_ster
    # And scale the image array accordingly:
    image_in_jypix = factor * image

    flux = np.sum(image_in_jypix)
   
    #### Shift image by mx and my in arcsec

    mvx_pix=(mx/(pixdeg_x*3600.0))
    mvy_pix=(my/(pixdeg_y*3600.0))
    shiftVector=(0.0, 0.0, mvy_pix, -mvx_pix) # minus sign as left is positive 


    # cp star and remove it
    istar=nx/2 
    Fstar=image_in_jypix[0,0,istar,istar]
    image_in_jypix[0,0,istar,istar]=0.0
    print "Fstar in Jy =",Fstar
    # shift
    image_in_jypix_shifted=shift(image_in_jypix,shift=shiftVector, order=3)#,mode='wrap')
    # add star in new position
    image_in_jypix_shifted[0,0,istar+int(mvy_pix),istar-int(mvx_pix)]=Fstar

    flux=np.sum(image_in_jypix_shifted)

    # Make FITS header information from canvas (this is necessary to
    # interact later with CASA). The commented lines were used when
    # not using the canvas
    

    canvas=fits.open(path_canvas)
    header = canvas[0].header # fits.Header()
    del header['Origin'] # necessary due to a comment that CASA adds automatically

    ###############################################
    ###### Commented, not used when using the canvas
    ################################################

    # #header['SIMPLE']='T'
    # header['BITPIX']=-32
    # # all the NAXIS are created automatically header['NAXIS']=2
    # header['OBJECT']='HD109085'
    # header['EPOCH']=2000.0
    # header['LONPOLE']=180.0
    # header['EQUINOX']=2000.0
    # header['SPECSYS']='LSRK'
    # header['RESTFREQ']=reffreq
    # header['VELREF']=0.0
    # header['CTYPE3']='FREQ'
    # header['CRPIX3'] = 1.0
    # header['CDELT3']  = 1.0
    # header['CRVAL3']= reffreq 
   
    # header['BTYPE'] = 'Intensity'
    # header['BSCALE'] = 1
    # header['BZERO'] = 0
    
    # #header['EPOCH'] = 2000.
    # #header['LONPOLE'] = 180.
    # header['CTYPE1'] = 'RA---TAN'
    # header['CRVAL1'] = 1.880158942527e2
    # header['CTYPE2'] = 'DEC--TAN'
    # header['CRVAL2'] = -1.619622628083e1
    

    # Figure out what units we are using per pixel:
    
    # if arcsec:
    #     unit = 'arcsec'
    #     multiplier = 3600.
    # elif mas:
    #     unit = 'mas'
    #     multiplier = 3.6e6
    # else:
    #     unit = 'DEG'
    #     multiplier = 1
    # print multiplier*pixdeg_x
    # print unit
    # header['CDELT1'] = -multiplier*pixdeg_x
    # header['CUNIT1'] = unit
    # #
    # # ...Zero point of coordinate system
    # #
    # header['CRPIX1'] = 1.0*((nx+1)/2)

    
    # header['CDELT2'] = multiplier*pixdeg_y
    # header['CUNIT2'] = unit
    # #
    # # ...Zero point of coordinate system
    # #
    # header['CRPIX2'] = 1.0* ((ny+1)/2)




    # if nf > 1:
    #     # multiple frequencies - set up the header keywords to define the
    #     #    third axis as frequency
    #     header['CTYPE3'] = 'VELOCITY'
    #     header['CUNIT3'] = 'km/s'
    #     header['CRPIX3'] = 1.0* ((nf+1)/2)
    #     header['CRVAL3'] = 0.0
    #     # Calculate the frequency step, assuming equal steps between all:
    #     delta_velocity = (lam[1] - lam[0])*cc*1e-5/lam0
    #     header['CDELT3'] = delta_velocity
    #     header['RESTWAVE']=lam0
    # else:                # only one frequency
    #     header['RESTWAVE'] = lam[0]
    #     header['CUNIT3'] = 'Hz'
    # # Make a FITS file!
    # #

    #########################################
    #########################################

    header['FLUX']=flux
    print "flux [Jy] = ", flux
   
    
    header['BUNIT'] = 'JY/PIXEL' 
    
    ##### multiply by primary beam or not. Only necessary when simulating visibilities in CASA with ft a posteriori

    if pbm:
        print "multiply by pb beam!"

        # multiply by primary beam
        pbfits=fits.open(path_canvas[:-4]+'pb.fits')
        pb=pbfits[0].data[0,0,:,:]
        image_in_jypix_shifted=image_in_jypix_shifted*pb

        inans= np.isnan(image_in_jypix_shifted)
        image_in_jypix_shifted[inans]=0.0
    else:
        print "don't multiply by pbm"

    image_in_jypix_float=image_in_jypix_shifted.astype(np.float32)
    print 'float'
    fits.writeto(path_fits, image_in_jypix_float, header, output_verify='fix')

def convert_to_fits_canvas_fields(path_image,path_fits,path_canvas, dpc, lam0, arcsec=False, mas=False, mx=0.0, my=0.0 , pbm=False, x0=0.0, y0=0.0):

    cc  = 2.9979245800000e10      # Light speed             [cm/s]
    pc  = 3.08572e18              # Parsec                  [cm]

    lam0=float(lam0)

    reffreq=cc/(lam0*1.0e-4)

    f=open(path_image,'r')
  
    iformat=int(f.readline())

    if (iformat < 1) or (iformat > 4):
        print("ERROR: File format of {0:s} not recognized.".format(filename))
        return

    if (iformat == 1) or (iformat == 3):
        radian = (1 == 0)

    else:
        radian = (1 == 1)

    nx, ny = tuple(np.array(f.readline().split(),dtype=int))
    nf = int(f.readline()) # number of wavelengths
    sizepix_x, sizepix_y = tuple(np.array(f.readline().split(),dtype=float))

    lam = np.empty(nf)
    for i in range(nf):
        lam[i] = float(f.readline())
    
    f.readline()  

    image = np.zeros((1,nf,ny,nx), dtype=float)

    for k in range(nf):
        for j in range(ny):
            for i in range(nx):

                image[0,k,j,i] = float(f.readline())

                if (j == ny-1) and (i == nx-1):
                    f.readline()

    f.close()

    # Compute the flux in this image as seen at dpc (pc)

    
    pc = 3.0857200e18   # pc in cm

    pixdeg_x = 180.0*(sizepix_x/(dpc*pc))/np.pi
    pixdeg_y = 180.0*(sizepix_y/(dpc*pc))/np.pi

    #print pixdeg_x
    # Compute the conversion factor from erg/cm^2/s/Hz/ster to erg/cm^2/s/Hz/ster at dpc
    #
    pixsurf_ster = pixdeg_x*pixdeg_y * (np.pi/180.)**2
    factor = 1e+23 * pixsurf_ster
    # And scale the image array accordingly:
    image_in_jypix = factor * image

   

    # Make FITS header information from canvas (this is necessary to
    # interact later with CASA). The commented lines were used when
    # not using the canvas
    

    canvas=fits.open(path_canvas)
    header = canvas[0].header # fits.Header()
    del header['Origin'] # necessary due to a comment that CASA adds automatically

    ###############################################
    ###### Commented, not used when using the canvas
    ################################################

    # #header['SIMPLE']='T'
    # header['BITPIX']=-32
    # # all the NAXIS are created automatically header['NAXIS']=2
    # header['OBJECT']='HD109085'
    # header['EPOCH']=2000.0
    # header['LONPOLE']=180.0
    # header['EQUINOX']=2000.0
    # header['SPECSYS']='LSRK'
    # header['RESTFREQ']=reffreq
    # header['VELREF']=0.0
    # header['CTYPE3']='FREQ'
    # header['CRPIX3'] = 1.0
    # header['CDELT3']  = 1.0
    # header['CRVAL3']= reffreq 
   
    # header['BTYPE'] = 'Intensity'
    # header['BSCALE'] = 1
    # header['BZERO'] = 0
    
    # #header['EPOCH'] = 2000.
    # #header['LONPOLE'] = 180.
    # header['CTYPE1'] = 'RA---TAN'
    # header['CRVAL1'] = 1.880158942527e2
    # header['CTYPE2'] = 'DEC--TAN'
    # header['CRVAL2'] = -1.619622628083e1
    

    # Figure out what units we are using per pixel:
    
    # if arcsec:
    #     unit = 'arcsec'
    #     multiplier = 3600.
    # elif mas:
    #     unit = 'mas'
    #     multiplier = 3.6e6
    # else:
    #     unit = 'DEG'
    #     multiplier = 1
    # print multiplier*pixdeg_x
    # print unit
    # header['CDELT1'] = -multiplier*pixdeg_x
    # header['CUNIT1'] = unit
    # #
    # # ...Zero point of coordinate system
    # #
    # header['CRPIX1'] = 1.0*((nx+1)/2)

    
    # header['CDELT2'] = multiplier*pixdeg_y
    # header['CUNIT2'] = unit
    # #
    # # ...Zero point of coordinate system
    # #
    # header['CRPIX2'] = 1.0* ((ny+1)/2)




    # if nf > 1:
    #     # multiple frequencies - set up the header keywords to define the
    #     #    third axis as frequency
    #     header['CTYPE3'] = 'VELOCITY'
    #     header['CUNIT3'] = 'km/s'
    #     header['CRPIX3'] = 1.0* ((nf+1)/2)
    #     header['CRVAL3'] = 0.0
    #     # Calculate the frequency step, assuming equal steps between all:
    #     delta_velocity = (lam[1] - lam[0])*cc*1e-5/lam0
    #     header['CDELT3'] = delta_velocity
    #     header['RESTWAVE']=lam0
    # else:                # only one frequency
    #     header['RESTWAVE'] = lam[0]
    #     header['CUNIT3'] = 'Hz'
    # # Make a FITS file!
    # #

    #########################################
    #########################################

   
    #### Shift image by mx and my in arcsec


    x1=header['CRVAL1']
    y1=header['CRVAL2']

    if x0==0.0 and y0==0.0:
        x0=x1
        y0=y1

    offx_arcsec=-(x1-x0)*np.cos(y0*np.pi/180.0)*3600.0 # offset projected in sky
    offy_arcsec=-(y1-y0)*3600.0

    mvx_pix=(mx+offx_arcsec)/(pixdeg_x*3600.0) 
    mvy_pix=(my+offy_arcsec)/(pixdeg_y*3600.0) 

    # print x0,x1,mx, offx_arcsec
    # print y0,y1,my, offy_arcsec

    shiftVector=(0.0, 0.0, mvy_pix, -mvx_pix) # minus sign as left is positive 


    # cp star and remove it
    istar=nx/2 
    Fstar=image_in_jypix[0,0,istar,istar]
    image_in_jypix[0,0,istar,istar]=0.0
    #print "Fstar in Jy =",Fstar
    
    # shift
    image_in_jypix_shifted=shift(image_in_jypix,shift=shiftVector, order=3)#,mode='wrap')
    # add star in new position
    image_in_jypix_shifted[0,0,istar+int(mvy_pix),istar-int(mvx_pix)]=Fstar

    flux=np.sum(image_in_jypix_shifted)

    
    header['BUNIT'] = 'JY/PIXEL' 
    
    ##### multiply by primary beam or not. Only necessary when simulating visibilities in CASA with ft a posteriori

    if pbm:
        print "multiply by pb beam!"
        # multiply by primary beam
        pbfits=fits.open(path_canvas[:-4]+'pb.fits')
        pb=pbfits[0].data[0,0,:,:]
        image_in_jypix_shifted=image_in_jypix_shifted*pb

        inans= np.isnan(image_in_jypix_shifted)
        image_in_jypix_shifted[inans]=0.0
    else:
        print "don't multiply by pbm"

    image_in_jypix_float=image_in_jypix_shifted.astype(np.float32)

    flux = np.sum(image_in_jypix_float)
    header['FLUX']=flux
    # print "flux [Jy] = ", flux
   

    fits.writeto(path_fits, image_in_jypix_float, header, output_verify='fix')






def Simimages_canvas_fields(dpc, X0, Y0, images, wavelengths, fields, Npix, dpix, canvas, inc, PA, offx, offy, pb=0.0, tag=''):

    # images: array of names for images produced at wavelengths
    # wavelgnths: wavelengths at which to produce images
    # fields: fields where to make images (=[0] unless observations are a mosaic)

    sau=Npix*dpix*dpc
    

    for im in xrange(len(images)):

        os.system('radmc3d image incl '+str(inc)+' phi 0.0 posang '+str(PA-90.0)+'  npix '+str(Npix)+' lambda '+str(wavelengths[im])+' sizeau '+str(sau)+' secondorder  > simimgaes.log')

        pathin ='image_'+images[im]+'_'+tag+'.out'
        pathout='image_'+images[im]+'_'+tag+'.fits'

        os.system('mv image.out '+pathin)

        for fi in fields:
            pathout='image_'+images[im]+'_'+tag+'_field'+str(fi)+'.fits'
            convert_to_fits_canvas_fields(pathin, pathout, canvas+str(fi)+'.fits' ,dpc, wavelengths[im], mx=offx, my=offy, pbm=pb, x0=X0, y0=Y0)
            os.system('mv '+pathout+' ./images')
            os.system('rm '+pathin)

    print "Done!"






