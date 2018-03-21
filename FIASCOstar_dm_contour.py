###################################################
##### FIASCO - Stellar disk vs. DM halo ###########
###################################################

#-----Import the required libraries------------------------
import numpy as np
from numpy import linspace
import scipy as sci
import scipy.stats as stats
from scipy.io.idl import readsav
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from mpl_toolkits.mplot3d import Axes3D
from scipy.optimize import curve_fit
from matplotlib.colors import LogNorm

np.seterr(divide='ignore')            

#-----Select the FIASCO halo-------------------------------

#-----Read DM halo data---------------------
#d = readsav('dn_60p0_ds_1p0_halo_0_dm_re.sav')               #LARGE DISK - MW-like
#d = readsav('dn_60p0_ds_1p0_halo_1_dm_re.sav')               #MERGER????
#d = readsav('dn_60p0_ds_1p0_halo_2_dm_re.sav')               #GOOD - MW-like
#d = readsav('dn_60p0_ds_1p0_halo_3_dm_re.sav')               #VERY FLAT DISK
#d = readsav('dn_60p0_ds_1p0_halo_4_dm_re.sav')               #GOOD - MW-like, POSSIBLE LATE POST-MERGER??
#d = readsav('dn_60p0_ds_1p0_halo_5_dm_re.sav')               #POST-MERGING EVENT????
#d = readsav('dn_60p0_ds_1p0_halo_6_dm_re.sav')               #POSSIBLE LATE SPIRAL - BARRED??
#d = readsav('dn_60p0_ds_1p0_halo_7_dm_re.sav')               #SMALL DISK - OFF-HORIZONAL
#d = readsav('dn_60p0_ds_1p0_halo_8_dm_re.sav')               #LARGE DISK - MW-like
#d = readsav('dn_60p0_ds_1p0_halo_9_dm_re.sav')               #POST-MERGER
#d = readsav('dn_60p0_ds_1p0_halo_10_dm_re.sav')              #GOOD - MW-like
d = readsav('dn_60p0_ds_1p0_halo_11_dm_re.sav')               #SMALL DISK - MW-like

#--------Read stellar disk data-------------
#s = readsav('dn_60p0_ds_1p0_halo_0_disk_re.sav')
#s = readsav('dn_60p0_ds_1p0_halo_1_disk_re.sav')
#s = readsav('dn_60p0_ds_1p0_halo_2_disk_re.sav')
#s = readsav('dn_60p0_ds_1p0_halo_3_disk_re.sav')
#s = readsav('dn_60p0_ds_1p0_halo_4_disk_re.sav')
#s = readsav('dn_60p0_ds_1p0_halo_5_disk_re.sav')
#s = readsav('dn_60p0_ds_1p0_halo_6_disk_re.sav')
#s = readsav('dn_60p0_ds_1p0_halo_7_disk_re.sav')
#s = readsav('dn_60p0_ds_1p0_halo_8_disk_re.sav')
#s = readsav('dn_60p0_ds_1p0_halo_9_disk_re.sav')
#s = readsav('dn_60p0_ds_1p0_halo_10_disk_re.sav')
s = readsav('dn_60p0_ds_1p0_halo_11_disk_re.sav')

#-----Define the variables from .sav file------------------

#-----DM  coordinates-----------------------
dx = d.dx_dm_new_re
dy = d.dy_dm_new_re
dz = d.dz_dm_new_re

#-----Stellar/disk coordinates--------------

dx_disk = s.x_disk_new
dy_disk = s.y_disk_new
dz_disk = s.z_disk_new

#-----Plot 2D histograms of the DM, gas and stars----------

#-----Define the histogram bins-------------
xmin   = min(dx)
xmax   = max(dx)
ymin   = min(dy)
ymax   = max(dy)
zmin   = min(dz)
zmax   = max(dz)
nxbins = 600
nybins = 600
nzbins = 600
xbins  = linspace(start = xmin, stop = xmax, num = nxbins)
ybins  = linspace(start = ymin, stop = ymax, num = nybins)
zbins  = linspace(start = zmin, stop = zmax, num = nzbins)

xmin_disk   = min(dx_disk)
xmax_disk   = max(dx_disk)
ymin_disk   = min(dy_disk)
ymax_disk   = max(dy_disk)
zmin_disk   = min(dz_disk)
zmax_disk   = max(dz_disk)
nxbins_disk = 600
nybins_disk = 600
nzbins_disk = 600
xbins_disk  = linspace(start = xmin_disk, stop = xmax_disk, num = nxbins_disk)
ybins_disk  = linspace(start = ymin_disk, stop = ymax_disk, num = nybins_disk)
zbins_disk  = linspace(start = zmin_disk, stop = zmax_disk, num = nzbins_disk)

#-----Plot dark matter component--------------------

Dxy, xedges, yedges = np.histogram2d(dx, dy, bins=(xbins,ybins))
Dxz, xedges, zedges = np.histogram2d(dx, dz, bins=(xbins,zbins))
Dyz, yedges, zedges = np.histogram2d(dy, dz, bins=(ybins,zbins))

#-----Plot stellar component--------------------

Sxy, sxedges, syedges = np.histogram2d(dx_disk, dy_disk, bins=(xbins_disk,ybins_disk))
Sxz, sxedges, szedges = np.histogram2d(dx_disk, dz_disk, bins=(xbins_disk,zbins_disk))
Syz, syedges, szedges = np.histogram2d(dy_disk, dz_disk, bins=(ybins_disk,zbins_disk))

levels = [0,0.5,1,1.5,2,2.5,3,3.5,4,4.5]

plt.figure()
plt.imshow(np.log10(Sxz.transpose()), extent=[xmin_disk,xmax_disk,zmin_disk,zmax_disk], interpolation='nearest', cmap=plt.cm.viridis, origin="lower")
plt.xlabel('x [kpc]')
plt.ylabel('z [kpc]')
plt.clim(0,5)
cb = plt.colorbar()
cb.set_ticks([0,1,2,3,4,5])
cb.set_ticklabels([r'$10^{0}$', r'$10^{1}$', r'$10^{2}$', r'$10^{3}$', r'$10^{4}$', r'$10^{5}$'])
cb.set_label('Counts')
cs = plt.contour(np.log10(Dxz.transpose()) , extent=[xmin,xmax,zmin,zmax], colors='k', linestyles = 'solid', levels=levels)
plt.clabel(cs, cs.levels, inline=True)
plt.axis('scaled')
plt.axis([-50, 50, -25, 25])
#plt.savefig('Plots/Contours/dn_60p0_ds_1p0_halo_11_contour.png', dpi=300)
plt.show()


