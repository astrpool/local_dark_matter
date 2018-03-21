############################################
##### FIASCO Halos - Gas, DM and Stars #####
############################################

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
d = readsav('dn_60p0_ds_1p0_halo_11_dm_re.sav')              #SMALL DISK - MW-like

#-------Read gas halo data------------------
#g = readsav('dn_60p0_ds_1p0_halo_0_gas_re.sav')
#g = readsav('dn_60p0_ds_1p0_halo_1_gas_re.sav')
#g = readsav('dn_60p0_ds_1p0_halo_2_gas_re.sav')
#g = readsav('dn_60p0_ds_1p0_halo_3_gas_re.sav')
#g = readsav('dn_60p0_ds_1p0_halo_4_gas_re.sav')
#g = readsav('dn_60p0_ds_1p0_halo_5_gas_re.sav')
#g = readsav('dn_60p0_ds_1p0_halo_6_gas_re.sav')
#g = readsav('dn_60p0_ds_1p0_halo_7_gas_re.sav')
#g = readsav('dn_60p0_ds_1p0_halo_8_gas_re.sav')
#g = readsav('dn_60p0_ds_1p0_halo_9_gas_re.sav')
#g = readsav('dn_60p0_ds_1p0_halo_10_gas_re.sav')
g = readsav('dn_60p0_ds_1p0_halo_11_gas_re.sav')

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

#-----Gas coordinates-----------------------

dx_gas = g.dx_gas_new_re
dy_gas = g.dy_gas_new_re
dz_gas = g.dz_gas_new_re

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

xmin_gas   = min(dx_gas)
xmax_gas   = max(dx_gas)
ymin_gas   = min(dy_gas)
ymax_gas   = max(dy_gas)
zmin_gas   = min(dz_gas)
zmax_gas   = max(dz_gas)
nxbins_gas = 600
nybins_gas = 600
nzbins_gas = 600
xbins_gas  = linspace(start = xmin_gas, stop = xmax_gas, num = nxbins_gas)
ybins_gas  = linspace(start = ymin_gas, stop = ymax_gas, num = nybins_gas)
zbins_gas  = linspace(start = zmin_gas, stop = zmax_gas, num = nzbins_gas)

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

#-----Plot gas component--------------------

Gxy, gxedges, gyedges = np.histogram2d(dx_gas, dy_gas, bins=(xbins_gas,ybins_gas), normed=True)
Gxz, gxedges, gzedges = np.histogram2d(dx_gas, dz_gas, bins=(xbins_gas,zbins_gas), normed=True)
Gyz, gyedges, gzedges = np.histogram2d(dy_gas, dz_gas, bins=(ybins_gas,zbins_gas), normed=True)

fig, axs = plt.subplots(3,3)

#-----x-y plane--------------
axs[0,0].imshow(np.log10(Gxy.transpose())+0.00000001, extent=[xmin_gas,xmax_gas,ymin_gas,ymax_gas], interpolation='nearest', cmap=plt.cm.gist_heat, origin="lower")
axs[0,0].set_aspect('equal')
axs[0,0].axis([-300,300,-300,300])
axs[0,0].set(xlabel = 'x [kpc]', ylabel = 'y [kpc]')

axs[0,1].imshow(np.log10(Gxy.transpose())+0.00000001, extent=[xmin_gas,xmax_gas,ymin_gas,ymax_gas], interpolation='nearest', cmap=plt.cm.gist_heat, origin="lower")
axs[0,1].set_aspect('equal')
axs[0,1].axis([-150,150,-150,150])
axs[0,1].set(xlabel = 'x [kpc]', ylabel = 'y [kpc]')

axs[0,2].imshow(np.log10(Gxy.transpose())+0.00000001, extent=[xmin_gas,xmax_gas,ymin_gas,ymax_gas], interpolation='nearest', cmap=plt.cm.gist_heat, origin="lower")
axs[0,2].set_aspect('equal')
axs[0,2].axis([-30,30,-30,30])
axs[0,2].set(xlabel = 'x [kpc]', ylabel = 'y [kpc]')

#-----x-z plane--------------
axs[1,0].imshow(np.log10(Gxz.transpose())+0.00000001, extent=[xmin_gas,xmax_gas,zmin_gas,zmax_gas], interpolation='nearest', cmap=plt.cm.gist_heat, origin="lower")
axs[1,0].set_aspect('equal')
axs[1,0].axis([-300,300,-300,300])
axs[1,0].set(xlabel = 'x [kpc]', ylabel = 'z [kpc]')

axs[1,1].imshow(np.log10(Gxz.transpose())+0.00000001, extent=[xmin_gas,xmax_gas,zmin_gas,zmax_gas], interpolation='nearest', cmap=plt.cm.gist_heat, origin="lower")
axs[1,1].set_aspect('equal')
axs[1,1].axis([-150,150,-150,150])
axs[1,1].set(xlabel = 'x [kpc]', ylabel = 'z [kpc]')

axs[1,2].imshow(np.log10(Gxz.transpose())+0.00000001, extent=[xmin_gas,xmax_gas,zmin_gas,zmax_gas], interpolation='nearest', cmap=plt.cm.gist_heat, origin="lower")
axs[1,2].set_aspect('equal')
axs[1,2].axis([-30,30,-30,30])
axs[1,2].set(xlabel = 'x [kpc]', ylabel = 'z [kpc]')

#-----y-z plane--------------
axs[2,0].imshow(np.log10(Gyz.transpose())+0.00000001, extent=[ymin_gas,ymax_gas,zmin_gas,zmax_gas], interpolation='nearest', cmap=plt.cm.gist_heat, origin="lower")
axs[2,0].set_aspect('equal')
axs[2,0].axis([-300,300,-300,300])
axs[2,0].set(xlabel = 'y [kpc]', ylabel = 'z [kpc]')

axs[2,1].imshow(np.log10(Gyz.transpose())+0.00000001, extent=[ymin_gas,ymax_gas,zmin_gas,zmax_gas], interpolation='nearest', cmap=plt.cm.gist_heat, origin="lower")
axs[2,1].set_aspect('equal')
axs[2,1].axis([-150,150,-150,150])
axs[2,1].set(xlabel = 'y [kpc]', ylabel = 'z [kpc]')

axs[2,2].imshow(np.log10(Gyz.transpose())+0.00000001, extent=[ymin_gas,ymax_gas,zmin_gas,zmax_gas], interpolation='nearest', cmap=plt.cm.gist_heat, origin="lower")
axs[2,2].set_aspect('equal')
axs[2,2].axis([-30,30,-30,30])
axs[2,2].set(xlabel = 'y [kpc]', ylabel = 'z [kpc]')

plt.tight_layout()
#plt.savefig('Plots/dn_60p0_ds_1p0_halo_11_gas.png', dpi=300)
plt.show()

#-----Plot dark matter component--------------------

Dxy, xedges, yedges = np.histogram2d(dx, dy, bins=(xbins,ybins), normed=True)
Dxz, xedges, zedges = np.histogram2d(dx, dz, bins=(xbins,zbins), normed=True)
Dyz, yedges, zedges = np.histogram2d(dy, dz, bins=(ybins,zbins), normed=True)

fig, axs = plt.subplots(3,3)

#-----x-y plane--------------
axs[0,0].imshow(np.log10(Dxy.transpose())+0.00000001, extent=[xmin,xmax,ymin,ymax], interpolation='nearest', cmap=plt.cm.bone, origin="lower")
axs[0,0].set_aspect('equal')
axs[0,0].axis([-300,300,-300,300])
axs[0,0].set(xlabel = 'x [kpc]', ylabel = 'y [kpc]')

axs[0,1].imshow(np.log10(Dxy.transpose())+0.00000001, extent=[xmin,xmax,ymin,ymax], interpolation='nearest', cmap=plt.cm.bone, origin="lower")
axs[0,1].set_aspect('equal')
axs[0,1].axis([-150,150,-150,150])
axs[0,1].set(xlabel = 'x [kpc]', ylabel = 'y [kpc]')

axs[0,2].imshow(np.log10(Dxy.transpose())+0.00000001, extent=[xmin,xmax,ymin,ymax], interpolation='nearest', cmap=plt.cm.bone, origin="lower")
axs[0,2].set_aspect('equal')
axs[0,2].axis([-30,30,-30,30])
axs[0,2].set(xlabel = 'x [kpc]', ylabel = 'y [kpc]')

#-----x-z plane--------------
axs[1,0].imshow(np.log10(Dxz.transpose())+0.00000001, extent=[xmin,xmax,zmin,zmax], interpolation='nearest', cmap=plt.cm.bone, origin="lower")
axs[1,0].set_aspect('equal')
axs[1,0].axis([-300,300,-300,300])
axs[1,0].set(xlabel = 'x [kpc]', ylabel = 'z [kpc]')

axs[1,1].imshow(np.log10(Dxz.transpose())+0.00000001, extent=[xmin,xmax,zmin,zmax], interpolation='nearest', cmap=plt.cm.bone, origin="lower")
axs[1,1].set_aspect('equal')
axs[1,1].axis([-150,150,-150,150])
axs[1,1].set(xlabel = 'x [kpc]', ylabel = 'z [kpc]')

axs[1,2].imshow(np.log10(Dxz.transpose())+0.00000001, extent=[xmin,xmax,zmin,zmax], interpolation='nearest', cmap=plt.cm.bone, origin="lower")
axs[1,2].set_aspect('equal')
axs[1,2].axis([-30,30,-30,30])
axs[1,2].set(xlabel = 'x [kpc]', ylabel = 'z [kpc]')

#-----y-z plane--------------
axs[2,0].imshow(np.log10(Dyz.transpose())+0.00000001, extent=[ymin,ymax,zmin,zmax], interpolation='nearest', cmap=plt.cm.bone, origin="lower")
axs[2,0].set_aspect('equal')
axs[2,0].axis([-300,300,-300,300])
axs[2,0].set(xlabel = 'y [kpc]', ylabel = 'z [kpc]')

axs[2,1].imshow(np.log10(Dyz.transpose())+0.00000001, extent=[ymin,ymax,zmin,zmax], interpolation='nearest', cmap=plt.cm.bone, origin="lower")
axs[2,1].set_aspect('equal')
axs[2,1].axis([-150,150,-150,150])
axs[2,1].set(xlabel = 'y [kpc]', ylabel = 'z [kpc]')

axs[2,2].imshow(np.log10(Dyz.transpose())+0.00000001, extent=[ymin,ymax,zmin,zmax], interpolation='nearest', cmap=plt.cm.bone, origin="lower")
axs[2,2].set_aspect('equal')
axs[2,2].axis([-30,30,-30,30])
axs[2,2].set(xlabel = 'y [kpc]', ylabel = 'z [kpc]')

plt.tight_layout()
#plt.savefig('Plots/dn_60p0_ds_1p0_halo_11_dm.png', dpi=300)
plt.show()

#-----Plot stellar component--------------------

Sxy, sxedges, syedges = np.histogram2d(dx_disk, dy_disk, bins=(xbins_disk,ybins_disk), normed=True)
Sxz, sxedges, szedges = np.histogram2d(dx_disk, dz_disk, bins=(xbins_disk,zbins_disk), normed=True)
Syz, syedges, szedges = np.histogram2d(dy_disk, dz_disk, bins=(ybins_disk,zbins_disk), normed=True)

fig, axs = plt.subplots(3,3)

#-----x-y plane--------------
axs[0,0].imshow(np.log10(Sxy.transpose())+0.00000001, extent=[xmin_disk,xmax_disk,ymin_disk,ymax_disk], interpolation='nearest', cmap=plt.cm.inferno, origin="lower")
axs[0,0].set_aspect('equal')
axs[0,0].axis([-300,300,-300,300])
axs[0,0].set(xlabel = 'x [kpc]', ylabel = 'y [kpc]')

axs[0,1].imshow(np.log10(Sxy.transpose())+0.00000001, extent=[xmin_disk,xmax_disk,ymin_disk,ymax_disk], interpolation='nearest', cmap=plt.cm.inferno, origin="lower")
axs[0,1].set_aspect('equal')
axs[0,1].axis([-150,150,-150,150])
axs[0,1].set(xlabel = 'x [kpc]', ylabel = 'y [kpc]')

axs[0,2].imshow(np.log10(Sxy.transpose())+0.00000001, extent=[xmin_disk,xmax_disk,ymin_disk,ymax_disk], interpolation='nearest', cmap=plt.cm.inferno, origin="lower")
axs[0,2].set_aspect('equal')
axs[0,2].axis([-30,30,-30,30])
axs[0,2].set(xlabel = 'x [kpc]', ylabel = 'y [kpc]')

#-----x-z plane--------------
axs[1,0].imshow(np.log10(Sxz.transpose())+0.00000001, extent=[xmin_disk,xmax_disk,zmin_disk,zmax_disk], interpolation='nearest', cmap=plt.cm.inferno, origin="lower")
axs[1,0].set_aspect('equal')
axs[1,0].axis([-300,300,-300,300])
axs[1,0].set(xlabel = 'x [kpc]', ylabel = 'z [kpc]')

axs[1,1].imshow(np.log10(Sxz.transpose())+0.00000001, extent=[xmin_disk,xmax_disk,zmin_disk,zmax_disk], interpolation='nearest', cmap=plt.cm.inferno, origin="lower")
axs[1,1].set_aspect('equal')
axs[1,1].axis([-150,150,-150,150])
axs[1,1].set(xlabel = 'x [kpc]', ylabel = 'z [kpc]')

axs[1,2].imshow(np.log10(Sxz.transpose())+0.00000001, extent=[xmin_disk,xmax_disk,zmin_disk,zmax_disk], interpolation='nearest', cmap=plt.cm.inferno, origin="lower")
axs[1,2].set_aspect('equal')
axs[1,2].axis([-30,30,-30,30])
axs[1,2].set(xlabel = 'x [kpc]', ylabel = 'z [kpc]')

#-----y-z plane--------------
axs[2,0].imshow(np.log10(Syz.transpose())+0.00000001, extent=[ymin_disk,ymax_disk,zmin_disk,zmax_disk], interpolation='nearest', cmap=plt.cm.inferno, origin="lower")
axs[2,0].set_aspect('equal')
axs[2,0].axis([-300,300,-300,300])
axs[2,0].set(xlabel = 'y [kpc]', ylabel = 'z [kpc]')

axs[2,1].imshow(np.log10(Syz.transpose())+0.00000001, extent=[ymin_disk,ymax_disk,zmin_disk,zmax_disk], interpolation='nearest', cmap=plt.cm.inferno, origin="lower")
axs[2,1].set_aspect('equal')
axs[2,1].axis([-150,150,-150,150])
axs[2,1].set(xlabel = 'y [kpc]', ylabel = 'z [kpc]')

axs[2,2].imshow(np.log10(Syz.transpose())+0.00000001, extent=[ymin_disk,ymax_disk,zmin_disk,zmax_disk], interpolation='nearest', cmap=plt.cm.inferno, origin="lower")
axs[2,2].set_aspect('equal')
axs[2,2].axis([-30,30,-30,30])
axs[2,2].set(xlabel = 'y [kpc]', ylabel = 'z [kpc]')

plt.tight_layout()
#plt.savefig('Plots/dn_60p0_ds_1p0_halo_11_star.png', dpi=300)
plt.show()
