#/usr/lib/python3.5
#==================================================
import pMonteCarlo as pmc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
    
#========== Set up parameters to pass to MCML ==========
# This is the example from Wang1995 from Tab. 1 and Fig. 3
writeDetData = 1 # 0 for no output, 1 for binary, 2 for ascii
detDataFilename = 'detData.pmc'
numPhotons = 10000
globalSeed = 1
backgroundIndex = 1.0
layer_leftZ_mm = [0.0]
layer_rightZ_mm = [5.0]
layer_index = [1.6]
layer_g = [0.6]
layer_us_permm = [2083.0]
layer_ua_permm = [0.001]
log_abs_profile = 0 # 0 for no logging, 1 for binary, 2 for ascii
absDataFilename = 'absData.pmc'
grid_x_min_mm = -1.0
grid_x_max_mm = 1.0
grid_x_n = 300
grid_y_min_mm = -1.0
grid_y_max_mm = 1.0
grid_y_n = 400
grid_z_min_mm = 0.0
grid_z_max_mm = 0.2
grid_z_n = 50

paramList=[numPhotons, globalSeed, backgroundIndex, layer_leftZ_mm,\
           layer_rightZ_mm, layer_index, layer_g, layer_us_permm, \
           layer_ua_permm, writeDetData, detDataFilename, \
           log_abs_profile, absDataFilename, grid_x_min_mm, \
           grid_x_max_mm, grid_x_n, grid_y_min_mm, grid_y_max_mm, \
           grid_y_n, grid_z_min_mm, grid_z_max_mm, grid_z_n]
paramString=','.join(str(x) for x in paramList) # see Notes

#========== Run MCML Code ==========
pmc.mcml(paramString)

#========== Read in detData ==========
if writeDetData != 0:
    pdata=np.nan
    if writeDetData == 1:
        pdata=pmc.input_det_bin(detDataFilename)
    elif writeDetData == 2:
        pdata=pmc.input_det_ascii(detDataFilename)

#========== Read in absData ==========
X={'z':np.linspace(grid_z_min_mm,grid_z_max_mm,grid_z_n),\
   'y':np.linspace(grid_y_min_mm,grid_y_max_mm,grid_y_n),\
   'x':np.linspace(grid_x_min_mm,grid_x_max_mm,grid_x_n)}
gdata=np.nan
if log_abs_profile != 0:
    grid_shape=(grid_z_n,grid_y_n,grid_x_n)
    if log_abs_profile == 1:
        gdata=pmc.input_abs_bin(absDataFilename,grid_shape)
    elif log_abs_profile == 2:
        gdata=pmc.input_abs_ascii(absDataFilename,grid_shape)

#========== Seperate out reflection from transmission ==========
dfr=pdata.drop(pdata[pdata['det flag']!=1].index) # Ref photons
dft=pdata.drop(pdata[pdata['det flag']!=2].index) # Trans photons
dfs=pdata.drop(pdata[pdata['det flag']!=3].index) # Specular photons

#========== Compute total reflection and transmission ==========
diffRef=sum(dfr['weight'].values)/numPhotons
specRef=sum(dfs['weight'].values)/numPhotons
totalRef=diffRef+specRef
print('Total reflection coefficient:', totalRef)
totalTrans=sum(dft['weight'].values)/numPhotons
print('Total transmission coefficient:', totalTrans)

#========== Plot angular diffuse reflectance ==========
refAngleHist_psr, refAngleBins_rad = pmc.angular_diff(dfr, 30, numPhotons) 
plt.figure(1)
plt.plot(refAngleBins_rad, refAngleHist_psr, 'ro')
plt.xlabel(r"$\pi$ radians")
plt.ylabel(r"$R ( \theta )$ sr$^{-1}$")

#========== Plot angular diffuse transmission ==========
transAngleHist_psr, transAngleBins_rad = pmc.angular_diff(dft, 30, numPhotons) 
plt.figure(2)
plt.plot(transAngleBins_rad, transAngleHist_psr, 'ro')
plt.xlabel(r"$\pi$ radians")
plt.ylabel(r"$T ( \theta )$ sr$^{-1}$")

#========== Plot temporal reflectance ==========
refTimeHist_perps, refTimeBins_ps = pmc.temporal_diff(dfr, 0.0, 1000.0, 30, numPhotons)
plt.figure(3)
plt.semilogy(refTimeBins_ps, refTimeHist_perps, 'ro')
plt.xlabel(r"$t$ (ps)")
plt.ylabel(r"$R ( t )$ ps$^{-1}$")

plt.show()

    #pdata[pdata['det flag']==1]['time'].hist(bins=100, weights=pdata[pdata['det flag']==1]['weight'], range=(0.0, 10.0))#Hist of time
    
    #pdata[pdata['det flag']==1].plot(x='x-pos',y='y-pos',kind='scatter')
    
    #========== Examples of plotting gdata ==========
    # pcolormesh(X['x'],X['y'],gdata[n]) # slice at nth plane in z-dir 
    # pcolormesh(X['z'],X['x'],gdata[:,X['y'].size/2]) # Slice at middle y-plane
    # plot(X['z'],sum(gdata.values,axis=(1,2))) # Sum x&y, plot as f(z)



#==================== Notes ====================
# Note: det flag is set as 
#       0=not_detected, 1=reflection, 2=transmission, 
#       3=specular reflection, 4=killed by absorption

# Note on slicing: pdata.iloc[irow,icol] or pdata['colname']
# If you just want a matrix with all of the data: Matrix=pdata.values

