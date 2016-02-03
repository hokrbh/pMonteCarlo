#!/usr/local/lib/python2.7/

import sys
import struct
import pMonteCarlo
import struct

def main():
    writeDetData = 1 # 0 for no output, 1 for binary, 2 for ascii
    detDataFilename = 'detData.pmc'
    numPhotons = 2
    globalSeed = 1
    backgroundIndex = 1.0
    layer_leftZ_mm = [0.0]
    layer_rightZ_mm = [100.0]
    layer_index = [1.5]
    layer_g = [0.0]
    layer_us_permm = [9.0]
    layer_ua_permm = [1.0]
    log_abs_profile = 1 # 0 for no logging, 1 for binary, 2 for ascii
    absDataFilename = 'absData.pmc'
    grid_x_min_mm = -1.0
    grid_x_max_mm = 1.0
    grid_x_n = 10
    grid_y_min_mm = -1.0
    grid_y_max_mm = 1.0
    grid_y_n = 10
    grid_z_min_mm = 0.0
    grid_z_max_mm = 1.0
    grid_z_n = 10

    paramString = '{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14},{15},{16},{17},{18},{19},{20},{21}'.format(numPhotons, globalSeed, backgroundIndex, layer_leftZ_mm, layer_rightZ_mm, layer_index, layer_g, layer_us_permm, layer_ua_permm, writeDetData, detDataFilename, log_abs_profile, absDataFilename, grid_x_min_mm, grid_x_max_mm, grid_x_n, grid_y_min_mm, grid_y_max_mm, grid_y_n, grid_z_min_mm, grid_z_max_mm, grid_z_n)
    
    #print(paramString)
    
    print(pMonteCarlo.mcml(paramString))
    
    if writeDetData != 0:
        with open( detDataFilename , 'rb' ) as detF:
            if writeDetData == 1:
                detFContent = detF.read()
                lineLength = 80 # Number of bytes per photon in data
                for i in range(0, numPhotons):	
                    photonData = struct.unpack("dddddddiIIIIxxxx", detFContent[lineLength*i:lineLength*(i+1)])
                    print(photonData)
            elif writeDetData == 2:
                for i in range(0, numPhotons):
                    photonData = detF.readline()
                    print(photonData)
            detF.closed
    
    

if __name__ == '__main__':
    main()
