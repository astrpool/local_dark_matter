########################################################
##### FIASCO - Local (Torus) Velocity Distribution #####
########################################################

#-----Import the required libraries------------------------
import time
import numpy as np
import scipy as sci
from scipy.io.idl import readsav
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

#-----Set inital parameters for script---------------------
start_time = time.time()
plt.rcParams.update({'figure.max_open_warning': 0})
np.seterr(divide='ignore')  

#-----Create Gaussian function-----------------------------
def GaussianFit(x, *params):
    y = np.zeros_like(x)
    for i in range(0, len(params), 3):
        ctr = params[i]
        peak = params[i+1]
        sigma = params[i+2]
        y = y + peak * np.exp( -((x - ctr)/(2*sigma))**2)
    return y

#-----Read halo data for all halos-------------------------

n_halos = 12

for i in range (0, n_halos):

	d = readsav('dn_60p0_ds_1p0_halo_'+str(i)+'_local_torus.sav')

#-----Define the variables from .sav file------------------

#-----DM velocities-------------------------
	vx = d.dvx_torus
	vy = d.dvy_torus
	vz = d.dvz_torus
	v  = np.sqrt(vx**2+vy**2+vz**2)

	N  = len(vx)

#-----Plot velocity distributions--------------------------

#-----Plot v - single Gaussian-----------------
	plt.figure()

	bins = 50

	n = plt.hist(v, bins = bins,color = 'b', alpha=0.3)

	x = (n[1][:-1])  
	y = n[0]                         

	guess = [180, max(y), 50]                  								       #Guess the [centre, peak, sigma]

	popt, pcov = curve_fit(GaussianFit, x, y, p0=guess)

	fit_multi = GaussianFit(x, *popt)

#-----Test goodness of fit---
	chisq1g = np.sum(((y-GaussianFit(x, *popt))**2)/(GaussianFit(x, *popt)))	   #Calculated using Pearson’s Chi-square Test
	dof = (len(x)-len(popt))
	r_chisq1g = chisq1g/dof

	plt.plot(x, fit_multi , 'k--', label='Gaussian Fit', lw=2)

	plt.axis([0,(max(v)+50),0,(max(n[0])+50)])
	plt.title('Halo '+str(i)+' - Guassian Fit')
	plt.xlabel('$v_{total}\/[kms^{-1}]$')
	plt.ylabel('Count')
	plt.legend(frameon=False)
	plt.text((0.8*(max(v)+50)),(0.87*(max(n[0])+50)), '$\chi^2_\\nu = %4.3f$' %(r_chisq1g))
	#plt.savefig('Plots/Velocities/dn_60p0_ds_1p0_halo_'+str(i)+'_local_torus_vGfit.png', dpi=300)

#-----Plot v - Mixed-Gaussian-----------------   
	plt.figure()                    

	n = plt.hist(v, bins = bins,  color = 'b', alpha=0.3)

	nComponents = 2																   #Set number of Guassians to be fitted	
                								     
	for j in range(nComponents-1):                        						   #Loop over number of Gaussians to fit
   		guess += [120, 0.25*max(y), 40]               								       #2nd Gaussian guess

	popt, pcov = curve_fit(GaussianFit, x, y, p0=guess)

	fit_multi = GaussianFit(x, *popt)

#-----Test goodness of fit---
	chisq2g = np.sum(((y-GaussianFit(x, *popt))**2)/(GaussianFit(x, *popt)))	   #Calculated using Pearson’s Chi-square Test
	dof = (len(x)-len(popt))
	r_chisq2g = chisq2g/dof

	plt.plot(x, fit_multi , 'k--', label='Mixed-Gaussian Fit', lw=2)

	yfit1 = GaussianFit(x,popt[0],popt[1],popt[2])
	plt.plot(x,yfit1,'r:', lw =2, label='Components')

	yfit2 = GaussianFit(x,popt[3],popt[4],popt[5])
	plt.plot(x,yfit2,'r:', lw=2)

	plt.axis([0,(max(v)+50),0,(max(n[0])+50)])
	plt.title('Halo '+str(i)+' - Mixed Gaussian Fit')
	plt.xlabel('$v_{total}\/[kms^{-1}]$')
	plt.ylabel('Count')
	plt.legend(frameon=False)
	plt.text((0.72*(max(v)+50)),(0.81*(max(n[0])+50)),'$\chi^2_\\nu = %4.3f$' %(r_chisq2g))
	#plt.savefig('Plots/Velocities/dn_60p0_ds_1p0_halo_'+str(i)+'_local_torus_vMixGfit.png', dpi=300)

#plt.show()
#plt.close('all')

end_time = (time.time()-start_time)
print('Run time = %f seconds' %(end_time))

