#/usr/lib/python3.5
#==================================================

from pMonteCarloC import *
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

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

def angular_diff(data, numAngleBins, numPhotons): 
    # Build a histogram of the angles of reflection photons
    # Bins are assumed to go from 0 to pi/2
    dalpha = np.pi/(2.0*numAngleBins) # Bin width of the angles
    dOmega = np.empty(numAngleBins)
    for i in range(0,numAngleBins):
        dOmega[i] = 4.0*np.pi*np.sin((i+0.5)*dalpha)*np.sin(dalpha*0.5)
    angles = np.arccos(np.fabs(data['z-vel'].values)) # Angles of Ref photons
    weights = data['weight'].values
    angleHist, angleBins = np.histogram(angles, bins=numAngleBins,\
        weights=weights, range=(0.0,0.5*np.pi))
    angleHist_psr = np.divide(angleHist, dOmega)/numPhotons
    angleBins_centers = (angleBins[0:numAngleBins]+0.5*dalpha)/np.pi
    return angleHist_psr, angleBins_centers
    
def temporal_diff(data, time_min, time_max, numTimeBins, numPhotons):
    # Build a histogram of the emission time of reflection photons
    dt = (time_max-time_min)/numTimeBins # Bin width
    times = data['time'].values
    weights = data['weight'].values
    timeHist, timeBins = np.histogram(times, bins=numTimeBins,\
        weights=weights, range=(time_min,time_max))
    timeHist_ps = timeHist/(dt*numPhotons)
    timeBins_centers = timeBins[0:numTimeBins]+0.5*dt
    return timeHist_ps, timeBins_centers
    
    
