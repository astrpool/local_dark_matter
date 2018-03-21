#################################################################
##### FIASCO - Averaged local (Torus) Velocity Distribution #####
#################################################################

#-----Import the required libraries------------------------
import numpy as np
import scipy as sci
import scipy.stats as stats
from scipy.io.idl import readsav
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from mpl_toolkits.mplot3d import Axes3D
from scipy.optimize import curve_fit


np.seterr(divide='ignore')  

#-----Select the FIASCO halo-------------------------------

#-----Read DM halo data---------------------
d0  = readsav('dn_60p0_ds_1p0_halo_0_local_torus.sav')
d1  = readsav('dn_60p0_ds_1p0_halo_1_local_torus.sav')
d2  = readsav('dn_60p0_ds_1p0_halo_2_local_torus.sav')
d3  = readsav('dn_60p0_ds_1p0_halo_3_local_torus.sav')
d4  = readsav('dn_60p0_ds_1p0_halo_4_local_torus.sav')
d5  = readsav('dn_60p0_ds_1p0_halo_5_local_torus.sav')
d6  = readsav('dn_60p0_ds_1p0_halo_6_local_torus.sav')
d7  = readsav('dn_60p0_ds_1p0_halo_7_local_torus.sav')
d8  = readsav('dn_60p0_ds_1p0_halo_8_local_torus.sav')
d9  = readsav('dn_60p0_ds_1p0_halo_9_local_torus.sav')
d10 = readsav('dn_60p0_ds_1p0_halo_10_local_torus.sav')
d11 = readsav('dn_60p0_ds_1p0_halo_11_local_torus.sav')

#-----Define the variables from .sav file------------------

#-----DM velocities and radial distance-------------------------
vx0 = d0.dvx_torus
vy0 = d0.dvy_torus
vz0 = d0.dvz_torus

vx1 = d1.dvx_torus
vy1 = d1.dvy_torus
vz1 = d1.dvz_torus

vx2 = d2.dvx_torus
vy2 = d2.dvy_torus
vz2 = d2.dvz_torus

vx3 = d3.dvx_torus
vy3 = d3.dvy_torus
vz3 = d3.dvz_torus

vx4 = d4.dvx_torus
vy4 = d4.dvy_torus
vz4 = d4.dvz_torus

vx5 = d5.dvx_torus
vy5 = d5.dvy_torus
vz5 = d5.dvz_torus

vx6 = d6.dvx_torus
vy6 = d6.dvy_torus
vz6 = d6.dvz_torus

vx7 = d7.dvx_torus
vy7 = d7.dvy_torus
vz7 = d7.dvz_torus

vx8 = d8.dvx_torus
vy8 = d8.dvy_torus
vz8 = d8.dvz_torus

vx9 = d9.dvx_torus
vy9 = d9.dvy_torus
vz9 = d9.dvz_torus

vx10 = d10.dvx_torus
vy10 = d10.dvy_torus
vz10 = d10.dvz_torus

vx11 = d11.dvx_torus
vy11 = d11.dvy_torus
vz11 = d11.dvz_torus

#-----Calculate the particle total velocity----------------
v0  = np.sqrt(vx0**2+vy0**2+vz0**2)
v1  = np.sqrt(vx1**2+vy1**2+vz1**2)
v2  = np.sqrt(vx2**2+vy2**2+vz2**2)
v3  = np.sqrt(vx3**2+vy3**2+vz3**2)
v4  = np.sqrt(vx4**2+vy4**2+vz4**2)
v5  = np.sqrt(vx5**2+vy5**2+vz5**2)
v6  = np.sqrt(vx6**2+vy6**2+vz6**2)
v7  = np.sqrt(vx7**2+vy7**2+vz7**2)
v8  = np.sqrt(vx8**2+vy8**2+vz8**2)
v9  = np.sqrt(vx9**2+vy9**2+vz9**2)
v10 = np.sqrt(vx10**2+vy10**2+vz10**2)
v11 = np.sqrt(vx11**2+vy11**2+vz11**2)

#-----Plot the velocity distributiomn for a halos----------
binsize = 50

n0 = plt.hist(v0, bins=binsize, histtype='step', lw=0.5, label='Halo0')
n1 = plt.hist(v1, bins=binsize, histtype='step', lw=0.5, label='Halo1')
n2 = plt.hist(v2, bins=binsize, histtype='step', lw=0.5, label='Halo2')
n3 = plt.hist(v3, bins=binsize, histtype='step', lw=0.5, label='Halo3')
n4 = plt.hist(v4, bins=binsize, histtype='step', lw=0.5, label='Halo4')
n5 = plt.hist(v5, bins=binsize, histtype='step', lw=0.5, label='Halo5')
n6 = plt.hist(v6, bins=binsize, histtype='step', lw=0.5, label='Halo6')
n7 = plt.hist(v7, bins=binsize, histtype='step', lw=0.5, label='Halo7')
n8 = plt.hist(v8, bins=binsize, histtype='step', lw=0.5, label='Halo8')
n9 = plt.hist(v9, bins=binsize, histtype='step', lw=0.5, label='Halo9')
n10 = plt.hist(v10, bins=binsize, histtype='step', lw=0.5, label='Halo10')
n11 = plt.hist(v11, bins=binsize, histtype='step', lw=0.5, label='Halo11')

y0 = n0[0]                  
x0 = (n0[1][:-1]) 
y1 = n1[0]                 
x1 = (n1[1][:-1])        
y2 = n2[0]
x2 = (n2[1][:-1]) 
y3 = n3[0]
x3 = (n3[1][:-1])
y4 = n4[0]
x4 = (n4[1][:-1])
y5 = n5[0]
x5 = (n5[1][:-1]) 
y6 = n6[0]
x6 = (n6[1][:-1])
y7 = n7[0]
x7 = (n7[1][:-1]) 
y8 = n8[0]
x8 = (n8[1][:-1]) 
y9 = n9[0]
x9 = (n9[1][:-1]) 
y10 = n10[0]
x10 = (n10[1][:-1])
y11 = n11[0]
x11 = (n11[1][:-1])

#-----Calculate the median of the distribution-------------
x_med = np.median([x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11],axis=0)
y_med = np.median([y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11],axis=0)

#-----Plot velocities--------
plt.figure(1)
plt.plot(x_med, y_med, 'b-', lw=2, label='Median')
plt.axis([0,500,0,900])
plt.title('Local Velocity Distribution for FIASCO Halos')
plt.xlabel('Velocity $[kms^{-1}]$')
plt.ylabel('Count')
plt.legend(loc='upper right', fontsize='small')
plt.savefig('Plots/Velocities/dn_60p0_ds_1p0_median_local_torus_v.png', dpi=300)

#-----Fit Guassians to the median distribution-------------
def GaussianFit(x, *params):
    y = np.zeros_like(x)
    for i in range(0, len(params), 3):
        ctr = params[i]
        amp = params[i+1]
        wid = params[i+2]
        y = y + amp * np.exp( -((x - ctr)/wid)**2)
    return y

#-----Fit single Gaussian to the distribution--------------
plt.figure(2)    
guess = [120, 650, 50]                  #Guess the [centre, amplitute, width]

popt, pcov = curve_fit(GaussianFit, x_med, y_med, p0=guess)

chisq2g = np.sum(((y_med-GaussianFit(x_med, *popt))**2)/(GaussianFit(x_med, *popt)))	   #Calculated using Pearson’s Chi-square Test
dof = (len(x_med)-len(popt))
r_chisq2g = chisq2g/dof

fit_multi = GaussianFit(x_med, *popt)

plt.plot(x_med, y_med, 'b-', lw=2, label='Median')
plt.plot(x_med, fit_multi , 'k--', label='Mixed-Gaussian Fit')

ci68 = np.array(np.percentile([y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11], (16, 84), axis=0))
ci95 = np.array(np.percentile([y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11], (2.5, 97.5), axis=0))

plt.fill_between(x_med, ci95[0], ci95[1], color='cyan', edgecolor='skyblue', label='95%')
plt.fill_between(x_med, ci68[0], ci68[1], color='deepskyblue', edgecolor='blue', label='68%')

plt.axis([0,500,0,900])
plt.title('Single Gaussian Fit')
plt.xlabel('Velocity $[kms^{-1}]$')
plt.ylabel('Count')
plt.text((0.75*(max(x_med))+50),(0.92*(max(y_med)+50)),'$\chi^2_\\nu = %4.3f$' %(r_chisq2g))
plt.legend(loc='upper right', fontsize='small',frameon=False)
plt.savefig('Plots/Velocities/dn_60p0_ds_1p0_median_local_torus_vGfit.png', dpi=300)

#-----Fitting mixed Gaussian----------------
plt.figure(3)
for i in range(1):                        #Loop over number of Gaussians to fit
    guess += [180, 200, 50]               #2nd Gaussian guess

popt, pcov = curve_fit(GaussianFit, x_med, y_med, p0=guess)

chisq2g = np.sum(((y_med-GaussianFit(x_med, *popt))**2)/(GaussianFit(x_med, *popt)))	   #Calculated using Pearson’s Chi-square Test
dof = (len(x_med)-len(popt))
r_chisq2g = chisq2g/dof

fit_multi = GaussianFit(x_med, *popt)

plt.plot(x_med, y_med, 'b-', lw=2, label='Median')
plt.plot(x_med, fit_multi , 'k--', label='Mixed-Gaussian Fit')

yfit1 = GaussianFit(x_med,popt[0],popt[1],popt[2])
plt.plot(x_med,yfit1,'r:', label='Components')

yfit2 = GaussianFit(x_med,popt[3],popt[4],popt[5])
plt.plot(x_med,yfit2,'r:')

ci68 = np.array(np.percentile([y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11], (16, 84), axis=0))
ci95 = np.array(np.percentile([y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11], (2.5, 97.5), axis=0))

plt.fill_between(x_med, ci95[0], ci95[1], color='cyan', edgecolor='skyblue', label='95%')
plt.fill_between(x_med, ci68[0], ci68[1], color='deepskyblue', edgecolor='blue', label='68%')

plt.axis([0,500,0,900])
plt.title('Mixed Gaussian Fit (2 components)')
plt.xlabel('Velocity $[kms^{-1}]$')
plt.ylabel('Count')
plt.text((0.75*(max(x_med))+50),(0.87*(max(y_med)+50)),'$\chi^2_\\nu = %4.3f$' %(r_chisq2g))
plt.legend(loc='upper right', fontsize='small', frameon=False)
plt.savefig('Plots/Velocities/dn_60p0_ds_1p0_median_local_torus_vMGfit.png', dpi=300)
plt.show(1)
plt.show(2)
plt.show(3)

