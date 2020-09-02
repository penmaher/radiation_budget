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

    #avg_data = np.sum(data * weights) / np.sum(weights)
    return data_weighted

def calc_global_mean(data, lat):
    '''Why integrate to find an average?
       The average is a course integral. It is more accurate to
       take the integral than to 'brute force' it with an average.
       The avergae will be smaller unless dlat is infinitly small.'''
        
    lat_rad         = np.deg2rad(lat)

    # area weight the latitude to account for differences in latitude 
    weights         = np.cos(lat_rad)

    # find the weights and then integrate
    area_weight     = integrate.trapz(weights,lat_rad)
    global_integral = integrate.trapz(data * weights, lat_rad)
    
    return global_integral / area_weight

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
            pl: precip * L
            sh: sensible heat
           cre: cloud radiative effect
        Units
        --------------
        All input data should be fluxes (including precip). All units are W/m2
    
        List of papers
        --------------
            Energy budget equations:
                DeAngelis et al. 2015: An observational radiative constraint 
                                       on hydrologic cycle intensification
            CRE equation:
                Allan 2011 "Combining satellite data and models to estimate cloud
                            radiative effect at the surface and in the atmosphere"  

        """

        self.flux_names     = ['lwds', 'lwut', 'lwus',
                               'swus', 'swut', 'swds', 'swdt']

        self.other_names     = ['sh','precip']

        self.budget_terms   = ['swa', 'pl', 'lwc', 'sh','net'] 

    def compute_energy_budget(self, data):
        """
        Compute the Atmospheric energy budget

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

        #lwc-swa = p*l (ie lh) + sh
        #net =  lwc - swa  - data['p'] - data['sh']  The form in Deangelis 2015 is 
        # not consistent with the totoal atmospheric forcing code below (ref pending)

        #testing magnitude - not nedded otherwise
        rad_loss_space = ( - area_weight_avg(data['lwut'],data['lat'], lat_axis=0)
                           + area_weight_avg(data['swdt'],data['lat'], lat_axis=0)
                           - area_weight_avg(data['swut'],data['lat'], lat_axis=0) )

        output = {'lwc':lwc, 'swa':swa, 'pl':data['p'], 'sh':data['sh']}  

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
          - cre, cre_suf, and acre can be validated against Allan 2011 Fig 5. (see ref above)

        """

        #TOA CRE all sky
        lwcre = data_clear_sky['lwut'] - data_all_sky['lwut']
        swcre = data_clear_sky['swut'] - data_all_sky['swut']
        cre  = lwcre + swcre 

        #surf CRE all sky
        lwcre_surf = (data_all_sky['lwds'] - data_clear_sky['lwds'] 
                     -data_all_sky['lwus'] + data_clear_sky['lwus']) 
        swcre_surf = (data_all_sky['swds'] - data_clear_sky['swds'] 
                     -data_all_sky['swus'] + data_clear_sky['swus'])
        cre_surf   = lwcre_surf + swcre_surf

        #atmospheric CRE all sky 
        alwcre =  lwcre - lwcre_surf
        aswcre =  swcre - swcre_surf
        acre   = alwcre + aswcre

        output = {'lwcre':lwcre,'swcre':swcre,'cre':cre, #TOA
                  'lwcre_surf':lwcre_surf,'swcre_surf':swcre_surf,'cre_surf':cre_surf, #surf
                   'alwcre':alwcre,'aswcre':aswcre,'acre':acre} #atmosphere

        return output

    def total_atmos_forcing(self, data_all_sky, data_clear_sky, data_budget):    

        # get the cloud forcing
        cre_output = self.compute_cre(data_all_sky, data_clear_sky)

        #get the clear sky forcing
        cs_forcing =  atm_cs_forcing(self, data_all_sky, data_clear_sky)

        #atmos forcing = toa forcing - surface forcing
        atm_sw_crf_cld = output['swcre'] - output['swcre_surf']
        atm_sw_crf_cs  = cs_forcing['swcre'] - cs_forcing['swcre_surf']
        #atmospheric all sky sw forcing
        atm_sw_crf     = atm_sw_crf_cld + atm_sw_crf_cs

        atm_lw_crf_cld = output['lwcre'] - output['lwcre_surf']
        atm_lw_crf_cs  = cs_forcing['lwcre'] - cs_forcing['lwcre_surf']
        #atmospheric all sky lw forcing
        atm_lw_crf     = atm_lw_crf_cld + atm_lw_crf_cs 

        #total = cre - cre_surf + PL + SH
        total_forcing = cre_output['cre'] - cre_output['cre_surf'] + data['p'] + data['sh']

        forcing = {'sw_crf_cld':atm_sw_crf_cld, 'sw_crf_cs':atm_sw_crf_cs, 'sw_crf_all':atm_sw_crf,
                   'lw_crf_cld':atm_lw_crf_cld, 'lw_crf_cs':atm_lw_crf_cs, 'lw_crf_all':atm_lw_crf,
                   'total': total_forcing}

        return forcing

    def atm_cs_forcing(self, data_all_sky, data_clear_sky):
        # the clear sky fluxes are known (unlike the cloud fluxes above)
        # so the clear sky forcing is the net flux (i.e. down - up) 

        sw_crf_toa_cs  =     data_all_sky['swdt'] - data_clear_sky['swut'] #SWDT_all = SWDT_clear 
        lw_crf_surf_cs =   data_clear_sky['lwds'] - data_clear_sky['lwus']
        sw_crf_surf_cs =   data_clear_sky['swds'] - data_clear_sky['swus']
        lw_crf_toa_cs  = - data_clear_sky['lwut']  #LWDT_cs = 0 (as does all-sky)

        cs_forcing = {'sw_toa':sw_crf_toa_cs,   'lw_toa':lw_crf_toa_cs, 
                     'sw_surf':sw_crf_surf_cs, 'lw_surf':lw_crf_surf_cs}

        return cs_forcing

    def global_avg_cre(self, cre, lat):

        global_cre = {}
        for var in cre.keys():
            global_cre[var] = calc_global_mean(cre[var], lat)    
            print('Global CRE for {} is: {:8.2f}'.format(var,global_cre[var]))
        
        return global_cre

    def global_avg_cre_comp(self, data_all_sky, data_clear_sky):
        '''Compute the global average of the CRE components.'''

        cre_output  = self.compute_cre(data_all_sky, data_clear_sky)
        glb_cre_comp = self.global_avg_cre(cre_output, data_all_sky['lat'])

        return glb_cre_comp

    def global_atmos_budget(self, budget, lat):

        global_budget = {}
        for var in ['pl', 'lwc', 'swa', 'sh', 'net']:
            global_budget[var] = calc_global_mean(budget[var], lat)
            print('Global energy budget {} is: {:8.2f}'.format(var,global_budget[var]))

        return global_budget

    def global_avg_flux_comp(self, data, lat):

        global_flux_comp = {}
        for var in ['lwut', 'lwus', 'lwds', 'swut', 'swdt', 'swus', 'swds']:
            global_flux_comp[var] = calc_global_mean(data[var], lat)

        return global_flux_comp

