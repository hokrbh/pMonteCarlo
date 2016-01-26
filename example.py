#!/usr/local/lib/python2.7/

import sys
import struct
import pMonteCarlo
from array import *

def main():
    writeDetData = 1
    detDataFilename = 'detData.pmc'
    numPhotons = 1000000
    globalSeed = 1
    backgroundIndex = 1.0
    layer_leftZ_mm = [0.0]
    layer_rightZ_mm = [100.0]
    layer_index = [1.5]
    layer_g = [0.0]
    layer_us_permm = [9.0]
    layer_ua_permm = [1.0]
    
    
    paramString = '{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10}'.format(numPhotons, globalSeed, backgroundIndex, layer_leftZ_mm, layer_rightZ_mm, layer_index, layer_g, layer_us_permm, layer_ua_permm, writeDetData, detDataFilename)
    
    #print(paramString)
    
    print(pMonteCarlo.mcml(paramString))
    
    with open( detDataFilename , 'r' ) as detF:
        
    detF.closed
    
    

if __name__ == '__main__':
    main()
