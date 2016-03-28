#/usr/lib/python3.5
#==================================================
import pmontecarlo as pmc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
    
#========== Run index matched example ==========
# This is the example from Wang1995 from Tab. 1 and Fig. 3
# Note that in Wang1995 he only plots the scattered photons in the 
# transmission plot, we make no distinction between scattered and 
# unscattered in transmission. This causes only the first point to 
# differ from his example.
print('Running index matched example')
writeDetData = 1 # 0 for no output, 1 for binary, 2 for ascii
detDataFilename = 'indexMatched.pmc'
numPhotons = 1000000
globalSeed = 1
backgroundIndex = 1.0
layer_leftZ_mm = [0.0]
layer_rightZ_mm = [0.2]
layer_index = [1.0]
layer_g = [0.75]
layer_us_permm = [9.0]
layer_ua_permm = [1.0]
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
plt.ylim([0,0.025])

#========== Plot angular diffuse transmission ==========
transAngleHist_psr, transAngleBins_rad = pmc.angular_diff(dft, 30, numPhotons) 
plt.figure(2)
plt.plot(transAngleBins_rad, transAngleHist_psr, 'ro')
plt.xlabel(r"$\pi$ radians")
plt.ylabel(r"$T ( \theta )$ sr$^{-1}$")
plt.ylim([0,0.8])

#========== Run index mismatched example ==========
# This is the example from Wang1995 from Tab. 2
print('Running index mismatched example')
writeDetData = 1 # 0 for no output, 1 for binary, 2 for ascii
detDataFilename = 'indexMismatched.pmc'
numPhotons = 250000
globalSeed = 1
backgroundIndex = 1.0
layer_leftZ_mm = [0.0]
layer_rightZ_mm = [100.0]
layer_index = [1.5]
layer_g = [0.0]
layer_us_permm = [9.0]
layer_ua_permm = [1.0]
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

plt.show()
