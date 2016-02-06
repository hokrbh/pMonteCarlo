#/usr/lib/python3.5
#==================================================
#import future # Uncomment this line for python2 (untested)
import sys
import struct
import pMonteCarlo
import pandas as pd
import numpy as np

#==================== extra Functions ====================
def input_bin(detDataFilename): 
    dt=np.dtype([('x-pos','f8'),('y-pos','f8'),('x-vel','f8'),\
                     ('y-vel','f8'),('z-vel','f8'),('weight','f8'),\
                     ('time','f8'),('det flag','i'),('seed1','uint32'),\
                     ('seed2','uint32'),('seed3','uint32'),\
                     ('seed4','uint32'),('pad','4V')])
    records=np.fromfile(detDataFilename,dt)
    df=pd.DataFrame(records).drop('pad',axis=1)
    return df

def input_ascii(detDataFilename):
    pheader=['x-pos','y-pos','x-vel','y-vel','z-vel','weight',\
             'time','det flag','seed1','seed2','seed3','seed4']
    df=pd.read_csv(detDataFilename,sep='\t',header=None,\
                   names=pheader)
    return df

#========================= Main =========================
def main():
    writeDetData = 1 # 0 for no output, 1 for binary, 2 for ascii
    detDataFilename = 'detData.pmc'
    numPhotons = 100000
    globalSeed = 1
    backgroundIndex = 1.0
    layer_leftZ_mm = [0.0]
    layer_rightZ_mm = [0.2]
    layer_index = [1.0]#5]
    layer_g = [0.750]
    layer_us_permm = [9.0]
    layer_ua_permm = [1.0]
    log_abs_profile = 2 # 0 for no logging, 1 for binary, 2 for ascii
    absDataFilename = 'absData.pmc'
    grid_x_min_mm = -1.0
    grid_x_max_mm = 1.0
    grid_x_n = 10
    grid_y_min_mm = -1.0
    grid_y_max_mm = 1.0
    grid_y_n = 10
    grid_z_min_mm = 0.0
    grid_z_max_mm = 0.2
    grid_z_n = 10

    paramList=[numPhotons, globalSeed, backgroundIndex, layer_leftZ_mm,\
               layer_rightZ_mm, layer_index, layer_g, layer_us_permm, \
               layer_ua_permm, writeDetData, detDataFilename, \
               log_abs_profile, absDataFilename, grid_x_min_mm, \
               grid_x_max_mm, grid_x_n, grid_y_min_mm, grid_y_max_mm, \
               grid_y_n, grid_z_min_mm, grid_z_max_mm, grid_z_n]
    paramString=','.join(str(x) for x in paramList) # see Notes

    print(pMonteCarlo.mcml(paramString)) # Run MCML Code

    #========== Read in Data ==========
    if writeDetData != 0:
        if writeDetData == 1:
            pdata=input_bin(detDataFilename)
        elif writeDetData == 2:
            pdata=input_ascii(detDataFilename)

    return pdata

#==================================================
if __name__ == '__main__':
    pdata=main()

    #pMatrix=pdata.values# Use this to convert DataFrame to numpy Matrix

    dfr=pdata.drop(pdata[pdata['det flag']!=1].index) # Ref photons
    dft=pdata.drop(pdata[pdata['det flag']!=2].index) # Trans photons
    dfa=pdata.drop(pdata[pdata['det flag']!=4].index) # absorb photons
    #pdata[pdata['det flag']==1]['time'].hist(bins=100)#Hist of time
    #pdata[pdata['det flag']==1].plot(x='x-pos',y='y-pos',kind='scatter')



#==================== Notes ====================
# Note: det flag is set as 
#       0=not_detected, 1=reflection, 2=transmission, 
#       3=specular reflection, 4=killed by absorption

# Note on slicing: pdata.iloc[irow,icol] or pdata['colname']
# If you just want a matrix with all of the data: Matrix=pdata.values

#========== Alternate String Construction ==========
# Old String method:
# paramString = '{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10}'.format(numPhotons, globalSeed, backgroundIndex, layer_leftZ_mm, layer_rightZ_mm, layer_index, layer_g, layer_us_permm, layer_ua_permm, writeDetData, detDataFilename)

# If we don't mind the string being sorted alphabetically by key
# we can use a dict:  (Would need to modify the parsing on c-side)
# myparams={'filename':'detData.pmc',\
#         'numPhotons':10000,\
#         'globalSeed':1,\
#         'backgroundIndex': 1.0,\
#         'layer_leftZ_mm':[0.0]}
# myl=list(myparams.values())     # convert dict to list
# paramString1=','.join(str(x) for x in myl) # convert list to string
# Or more simply, just build a list from variables like I did in program

#==================================================
# def input_bin2(defFilename): # Better for short files(low numPhoton)
#     with open( detDataFilename , 'rb' ) as detF:
#         detFContent = detF.read()
#         lineLength = 80 # Number of bytes per photon in data
#         photonData=[]
#         for i in range(0, numPhotons):	
#             photonData.append(struct.unpack("dddddddiIIIIxxxx", detFContent[lineLength*i:lineLength*(i+1)]))
#         detF.closed
#         pheader=['x-pos','y-pos','x-vel','y-vel','z-vel','weight',\
#                 'time','det flag','seed1','seed2','seed3','seed4']
#         df=pd.DataFrame(photonData)
#         df.columns=pheader
#     return df
