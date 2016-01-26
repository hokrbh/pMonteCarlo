import sys
import struct
import pMonteCarlo

def main():
    numPhotons = 100000
    globalSeed = 1
    backgroundIndex = 1.0
    layer_leftZ_mm = [0.0]
    layer_rightZ_mm = [1000.0]
    layer_index = [1.5]
    layer_g = [0.0]
    layer_us_permm = [9.0]
    layer_ua_permm = [1.0]
    
    
    paramString = '{0},{1},{2},{3},{4},{5},{6},{7},{8}'.format(numPhotons, globalSeed, backgroundIndex, layer_leftZ_mm, layer_rightZ_mm, layer_index, layer_g, layer_us_permm, layer_ua_permm)
    
    #print(paramString)
    
    print(pMonteCarlo.mcml(paramString))
