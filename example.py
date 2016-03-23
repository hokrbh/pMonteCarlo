#/usr/lib/python3.5
#==================================================
#import future # Uncomment this line for python2 (untested)
import sys
import struct
import pMonteCarlo
import pandas as pd
import numpy as np

#==================== extra Functions ====================
def input_det_bin(detDataFilename): 
    dt=np.dtype([('x-pos','f8'),('y-pos','f8'),('x-vel','f8'),\
                     ('y-vel','f8'),('z-vel','f8'),('weight','f8'),\
                     ('time','f8'),('det flag','i'),('seed1','uint32'),\
                     ('seed2','uint32'),('seed3','uint32'),\
                     ('seed4','uint32'),('pad','4V')])
    records=np.fromfile(detDataFilename,dt)
    return pd.DataFrame(records).drop('pad',axis=1)

def input_abs_bin(absDataFilename,grid_shape): # faster than ascii
    records=np.fromfile(absDataFilename)
    return pd.Panel(records.reshape(grid_shape))

def input_det_ascii(detDataFilename):
    '''Import detData (return as a pd.DataFrame)'''
    pheader=['x-pos','y-pos','x-vel','y-vel','z-vel','weight',\
             'time','det flag','seed1','seed2','seed3','seed4']
    return pd.read_csv(detDataFilename,sep='\t',header=None,\
                       names=pheader)

def input_abs_ascii(absDataFilename,grid_shape): # slower than bin
    '''Import absData (return as a pd.Panel)'''
    idx=pd.MultiIndex.from_product([np.arange(grid_shape[0]),\
                                   np.arange(grid_shape[1])])
    df2=pd.read_csv(absDataFilename,sep='\t',header=None).dropna(axis=1) 
    df2.set_index(idx,inplace=True)                                      
    return pd.DataFrame.to_panel(df2).transpose(1,2,0) 

#========================= MAIN =========================
def main():
    writeDetData = 2 # 0 for no output, 1 for binary, 2 for ascii
    detDataFilename = 'detData.pmc'
    numPhotons = 10000
    globalSeed = 1
    backgroundIndex = 1.0
    layer_leftZ_mm = [0.0]
    layer_rightZ_mm = [5.0]
    layer_index = [1.6]
    layer_g = [0.0]
    layer_us_permm = [2083.0]#For BaSO4
    layer_ua_permm = [0.01]
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

    print(pMonteCarlo.mcml(paramString)) # Run MCML Code

    #========== Read in detData ==========
    if writeDetData != 0:
        pdata=np.nan
        if writeDetData == 1:
            pdata=input_det_bin(detDataFilename)
        elif writeDetData == 2:
            pdata=input_det_ascii(detDataFilename)

    #========== Read in absData ==========
    X={'z':np.linspace(grid_z_min_mm,grid_z_max_mm,grid_z_n),\
       'y':np.linspace(grid_y_min_mm,grid_y_max_mm,grid_y_n),\
       'x':np.linspace(grid_x_min_mm,grid_x_max_mm,grid_x_n)}
    gdata=np.nan
    if log_abs_profile != 0:
        grid_shape=(grid_z_n,grid_y_n,grid_x_n)
        if log_abs_profile == 1:
            gdata=input_abs_bin(absDataFilename,grid_shape)
        elif log_abs_profile == 2:
            gdata=input_abs_ascii(absDataFilename,grid_shape)


    return pdata,gdata,X

#==================================================
if __name__ == '__main__':
    pdata,gdata,X=main()
    #pMatrix=pdata.values# Use this to convert DataFrame to numpy Array

    dfr=pdata.drop(pdata[pdata['det flag']!=1].index) # Ref photons
    dft=pdata.drop(pdata[pdata['det flag']!=2].index) # Trans photons
    dfa=pdata.drop(pdata[pdata['det flag']!=4].index) # absorb photons
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

