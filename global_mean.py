from __future__ import division
import numpy as np
import pdb
from scipy import integrate

__author__ = 'Penelope Maher'

def area_weight_avg(data, lat, lat_axis):
    '''Only use this for testing or plotting. This is a rough test. 
       Use calc_global_mean instead'''

    weights = np.cos(np.radians(lat))

    return np.average(data, weights=weights,
                          axis=lat_axis)

def area_weight_data(data, lat):
    '''Used for plotting '''

    weights = np.cos(np.radians(lat))

    if len(data.shape) == 1:        
        data_weighted = data * weights
    elif len(data.shape) ==2:
        data_weighted = data * weights[:,None]
    else:
        print('Check dimension of data')
        pdb.set_trace()

    return data_weighted

def calc_global_mean(data, lat):
    '''Why integrate to find an average?
       The average is an integral. It is more accurate to
       take the integral than to 'brute force' it with an average.
       The avergae will be smaller unless dlat is infinitly small.'''
        
    lat_rad         = np.deg2rad(lat)

    # area weight the latitude to account for differences in latitude 
    weights         = np.cos(lat_rad)

    # find the weights and then integrate
    area_weight     = integrate.trapz(weights,lat_rad)
    global_integral = integrate.trapz(data * weights, lat_rad)
    
    return global_integral / area_weight

