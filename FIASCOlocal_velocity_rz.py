#####################################################################
##### FIASCO - Local (Torus) Velocity Distribution (r,z) #####
#####################################################################

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
d = readsav('dn_60p0_ds_1p0_halo_0_local_torus.sav')
#d = readsav('dn_60p0_ds_1p0_halo_1_local_torus.sav')
#d = readsav('dn_60p0_ds_1p0_halo_2_local_torus.sav')
#d = readsav('dn_60p0_ds_1p0_halo_3_local_torus.sav')
#d = readsav('dn_60p0_ds_1p0_halo_4_local_torus.sav')
#d = readsav('dn_60p0_ds_1p0_halo_5_local_torus.sav')
#d = readsav('dn_60p0_ds_1p0_halo_6_local_torus.sav')
#d = readsav('dn_60p0_ds_1p0_halo_7_local_torus.sav')
#d = readsav('dn_60p0_ds_1p0_halo_8_local_torus.sav')
#d = readsav('dn_60p0_ds_1p0_halo_9_local_torus.sav')
#d = readsav('dn_60p0_ds_1p0_halo_10_local_torus.sav')
#d = readsav('dn_60p0_ds_1p0_halo_11_local_torus.sav')

#-----Define the variables from .sav file------------------

#-----DM velocities-------------------------
vx = d.dvx_torus
vy = d.dvy_torus
vz = d.dvz_torus
vr  = np.sqrt(vx**2+vy**2)

n  = len(vx)

#-----Create Gaussian function--------------
def GaussianFit(x, *params):
    y = np.zeros_like(x)
    for i in range(0, len(params), 3):
        ctr = params[i]
        peak = params[i+1]
        fwhm = params[i+2]
        y = y + peak * np.exp( -((x - ctr)/fwhm)**2)
    return y

#-----Plot velocity distributions--------------------------

#-----Plot vx, vy and vz--------------------
bins = 50

fig, axs = plt.subplots(1,2, figsize = (15,6))

vz_n = axs[1].hist(vz, bins = bins, histtype = 'step', label = 'FIASCO')
axs[1].set_aspect('equal')
axs[1].axis([-400,400,0,(max(vz_n[0])+50)])
axs[1].set(ylabel = 'Counts', xlabel = '$v_z [kms^{-1}]$')
axs[1].legend()

vr_n = axs[0].hist(vr, bins = bins, histtype = 'step', label = 'FIASCO')
axs[0].set_aspect('equal')
axs[0].axis([-400,400,0,(max(vz_n[0])+50)])
axs[0].set(ylabel = 'Counts', xlabel = '$v_r [kms^{-1}]$')

plt.tight_layout()
plt.close()

#-----Define x and y for Gaussian fit-------
vr_x = vr_n[1][:-1]
vr_y = vr_n[0]

vz_x = vz_n[1][:-1]
vz_y = vz_n[0]

#-----Guess Gaussian fit parameters---------
vr_guess = [1, max(vr_n[0]), 250]
vz_guess = [1, max(vz_n[0]), 250]

#-----Obtain best fit parameters------------
vr_popt, vr_pcov = curve_fit(GaussianFit, vr_x, vr_y, p0 = vr_guess)
vz_popt, vz_pcov = curve_fit(GaussianFit, vz_x, vz_y, p0 = vz_guess)

#-----Obtain Gaussian y coordinates----------
vr_fit_multi = GaussianFit(vr_x, *vr_popt)
vz_fit_multi = GaussianFit(vz_x, *vz_popt)

#-----Plot velocity distributions with Gaussian fits-------

#-----Plot vx, vy and vz--------------------
fig, axs = plt.subplots(1,2, figsize = (15,6))

axs[0].hist(vr, bins = bins, label = 'FIASCO', alpha=0.3,  color = 'b')
axs[0].plot(vr_x, vr_fit_multi, 'k--', label = 'Gaussian Fit')
axs[0].set_aspect('equal')
axs[0].axis([0,400,0,(max(vz_n[0])+50)])
axs[0].set(ylabel = 'Counts', xlabel = '$v_r \/[kms^{-1}]$')
axs[0].text((180),(max(vz_n[0])), '$N_{DM} = %d$' %(n))

axs[1].hist(vz, bins = bins, label = 'FIASCO', alpha=0.3,  color = 'b')
axs[1].plot(vz_x, vz_fit_multi, 'k--', label = 'Gaussian Fit')
axs[1].set_aspect('equal')
axs[1].axis([-400,400,0,(max(vz_n[0])+50)])
axs[1].set(ylabel = 'Counts', xlabel = '$v_z \/[kms^{-1}]$')
axs[1].legend()

plt.tight_layout()
#plt.savefig('Plots/Velocities/dn_60p0_ds_1p0_halo_0_local_torus_vxyzGfit.png', dpi=300)
plt.show()

