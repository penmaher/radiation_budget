from __future__ import division
import numpy as np
import pdb


# Assumptions about the input data.
#  - All data is the time mean and zonal mean. 

class EnergyBudget():

    def __init__(self):
        """ Initialize the labels for the budget dictionary.

        List of acronyms
        -----------------
            lw: longwave
            sw: shortwave
           u,d: up and down (direction of fluxes to make sense of sign)
           t,s: top of atmosphere and surface
           swa: shortwave absorption
           lwc: longwave cooling
            pl: precip * L (budget output not input data)
            sh: sensible heat
           cre: cloud radiative effect
           crf: cloud radiative forcing
    
        List of papers
        --------------
            Flux equations:
               DeAngelis et al. 2015: An observational radiative constraint 
                                      on hydrologic cycle intensification
            CRE equation:
              Miller et al 2011: The Radiation Budget of the West African Sahel
                                 and Its Controls: A Perspective from Observations
                                 and Global Climate Models
        """

        self.flux_names     = ['lwds', 'lwut', 'lwus',
                               'swus', 'swut', 'swds', 'swdt']

        self.budget_terms   = ['swa', 'pl', 'lwc', 'sh'] 

    def area_weight_avg(self, data, lat, lat_axis):
        #this is only needed as a test and not applied to any data
        return np.average(data, weights=np.cos(np.radians(lat)),
                          axis=lat_axis).mean()

    def compute_radiation_budget(self, data):
        """
        Compute the radiatve budget

        Parameters
        ----------
        data : dictionary of data with keys:  [self.flux_names, lat]
               Each dictionary element is the time mean, zonal mean and is a 
               1D array of shate latitude 
        Returns
        -------
        output : a dictionary with keys: lwc, swa, pl

        """

        # long wave cooling : positive -> cooling.
        # net LW radiation lost to space
        lwc = data['lwut'] + (data['lwds'] - data['lwus'])      
    
        # short wave absorption : positive -> heating.
        # net sw absorbed by atmosphere
        swa = (data['swdt'] - data['swut']) - (data['swds'] - data['swus']) 

        # lwc is the biggest term by far
        pl = ( lwc - swa - data['sh'] )  #units are #W/m2

        #testing magnitude
        rad_loss_space = ( - self.area_weight_avg(data['lwut'],data['lat'], lat_axis=0)
                           + self.area_weight_avg(data['swdt'],data['lat'], lat_axis=0)
                           - self.area_weight_avg(data['swut'],data['lat'], lat_axis=0) )
    
        output = {'lwc':lwc, 'swa':swa, 'pl':pl}  

        return output 

    def compute_cre(self, data, data_clear_sky):
        """
        Compute the cloud radiative effect

        Parameters
        ----------
        data : dictionary of data with keys:  [self.flux_names, lat]
               Each dictionary element is the time mean, zonal mean and is a 
               1D array of shate latitude.
               Fluxes are all-sky (all sky = clear-sky + cloudy-sky) 
        data : as in data, but for clear-sky fluxes

        Returns
        -------
        output : a dictionary with keys: crf_s, crf_t, cre

        """

        # compute the cloud radiative forcing:
        #    - at the surface
        crf_surf = data['lwds'] - data_clear_sky['lwds'] + data['swds'] - data_clear_sky['swds']

        #    - at the top of atmosphere
        crf_toa  = data['lwut'] - data_clear_sky['lwut'] + data['swut'] - data_clear_sky['swut']

        cre = crf_toa - crf_surf + data['swus'] - data_clear_sky['swus']

        output = {'crf_s':crf_surf, 'crf_t':crf_toa, 'cre':cre}

        return output
