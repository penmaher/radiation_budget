from __future__ import division
import numpy as np
import pdb

__author__ = 'Penelope Maher'

def area_weight_avg(data, lat, lat_axis):
    #this is only needed as a test and not applied to any data
    return np.average(data, weights=np.cos(np.radians(lat)),
                          axis=lat_axis)

def area_weight_data(data, lat):

    weights = np.cos(np.radians(lat))
    if len(data.shape) == 1:        
        data_weighted = data * weights
    elif len(data.shape) ==2:
        data_weighted = data * weights[:,None]
    else:
        print('Check dimension of data')
        pdb.set_trace()
    #avg_data = np.sum(data * weights) / np.sum(weights)
    return data_weighted

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
    
        List of papers
        --------------
            Flux equations:
                DeAngelis et al. 2015: An observational radiative constraint 
                                       on hydrologic cycle intensification
            CRE equation:
                Allen 2011 "Combining satellite data and models to estimate cloud
                            radiative effect at the surface and in the atmosphere"  

        """

        self.flux_names     = ['lwds', 'lwut', 'lwus',
                               'swus', 'swut', 'swds', 'swdt']

        self.budget_terms   = ['swa', 'pl', 'lwc', 'sh'] 


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
        rad_loss_space = ( - area_weight_avg(data['lwut'],data['lat'], lat_axis=0)
                           + area_weight_avg(data['swdt'],data['lat'], lat_axis=0)
                           - area_weight_avg(data['swut'],data['lat'], lat_axis=0) )
    
        output = {'lwc':lwc, 'swa':swa, 'pl':pl, 'sh':data['sh']}  

        return output 

    def compute_cre(self, data_all_sky, data_clear_sky):
        """
        Compute the cloud radiative effect

        Parameters
        ----------
            data_all_sky : dict keys are self.flux_names and lat
                           each dict element is time mean and zonal mean
                           each dict element is a 1D array of shape latitude.
            data_clear_sky : as above but for clear-sky fluxes

        Returns
        -------
            output : a dictionary with each of the computed componets


        Balance values:
        -------

            IPCC AR5 Chapter 7 p 580 global average (TOA)
                SWCRE = -50 W/m2
                LWCRE = +30 W/m2
                CRE   = -20 W/m2 (ie the net effect of clouds is a cooling)

        Notes:
        -------
          - we have for a long time computed CRE at TOA. Improvements in 
            remote sensing has made it possible to estimate the surface CRE
            and hence within the atmosphere.
          - cloud forcing  F = LWF - SWF (the sign depends on how you define the fluxes)
          - all sky forcing = cloudy sky forcing + clear sky forcing
          - cloud radiatve effect (CRE) and cloud radiative forcing (CRF) are used interchangeably

        """

        #TOA CRE 
        lwcre = data_clear_sky['lwut'] - data_all_sky['lwut']
        swcre = data_clear_sky['swut'] - data_all_sky['swut']
        cre  = lwcre + swcre 

        #surf CRE
        lwcre_surf = (data_all_sky['lwds'] - data_clear_sky['lwds'] -
                      data_all_sky['lwus'] + data_clear_sky['lwus']) 
        swcre_surf = (data_all_sky['swds'] - data_clear_sky['swds'] - 
                      data_all_sky['swus'] + data_clear_sky['swus'])
        cre_surf   = lwcre_surf + swcre_surf

        #atmospheric CRE 
        alwcre =  lwcre - lwcre_surf
        aswcre =  swcre - swcre_surf
        acre   = alwcre + aswcre

        output = {'lwcre':lwcre,'swcre':swcre,'cre':cre, #TOA
                  'lwcre_surf':lwcre_surf,'swcre_surf':swcre_surf,'cre_surf':cre_surf, #surf
                   'alwcre':alwcre,'aswcre':aswcre,'acre':acre} #atmosphere

        return output
    
    def global_avg_cre(self, cre, lat):

        print_values = True
        global_cre = {}
        for var in cre.keys():
            weighted_lat = area_weight_data(cre[var], lat)
            global_cre[var] = weighted_lat.mean()
            if print_values:
                print('Global CRE for {} is: {:8.2f}'.format(var,global_cre[var]))
        
        return global_cre

    def global_atmos_budget(self, budget, lat):
        budget_weighted = {}
        for var in ['pl', 'lwc', 'swa', 'sh']:
            weighted_lat = area_weight_data(budget[var], lat) 
            budget_weighted[var] = weighted_lat.mean()
            print('Global energy budget {} is: {:8.2f}'.format(var,budget_weighted[var]))

        return budget_weighted
