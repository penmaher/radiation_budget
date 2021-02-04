from __future__ import division
import numpy as np
import pdb

from .global_mean import calc_global_mean

__author__ = 'Penelope Maher'


class AtmosEnergyBudget():

    def __init__(self, data, lh, sh):
        """ Initialize the labels for the budget dictionary.

        List of acronyms
        -----------------
            lw: longwave
            sw: shortwave
           u,d: up and down (direction of fluxes to make sense of sign)
           t,s: top of atmosphere and surface
            lh: latent heat
            sh: sensible heat
        Units
        --------------
            All input data should be fluxes (including precip). 
            All units are W/m2.
    
        List of papers
        --------------
            Energy budget equation:
                DeAngelis et al. 2015: An observational radiative constraint 
                                       on hydrologic cycle intensification
                Pendergrass and Hartmann 2013: The Atmospheric Energy Constraint
                                       on Global-Mean Precipitation Change
        """

        self.data     = data # a dictionary of all-sky flux arrays
        self.lh       = lh           # lh = (L_f x snow + L_c x precip) in W/m2
        self.sh       = sh

        #budget terms
        self.lwc = None
        self.swa = None
        self.net = None

    def compute_lwc(self):
        # long wave cooling : positive -> cooling.
        # net LW radiation lost to space
        self.lwc = self.data['lwut'] + (self.data['lwds'] - self.data['lwus'])      
    

    def compute_swa(self):
        # short wave absorption : positive -> heating.
        # net sw absorbed by atmosphere
        self.swa = (self.data['swdt'] - self.data['swut']) - (self.data['swds'] - self.data['swus']) 


    def compute_atmos_budget(self):

        self.compute_lwc()
        self.compute_swa()

        self.net = - self.lwc + self.swa + self.lh + self.sh

    def atm_cs_forcing(self, data_all_sky, data_clear_sky):
        # the clear sky fluxes are known (unlike the cloud fluxes above)
        # so the clear sky forcing is the net flux (i.e. down - up) 

        sw_crf_toa_cs  =   data_clear_sky['swdt'] - data_clear_sky['swut'] 
        lw_crf_surf_cs =   data_clear_sky['lwds'] - data_clear_sky['lwus']
        sw_crf_surf_cs =   data_clear_sky['swds'] - data_clear_sky['swus']
        lw_crf_toa_cs  = - data_clear_sky['lwut']  #LWDT_cs = 0 (as does all-sky)

        cs_forcing = {'sw_toa':sw_crf_toa_cs,   'lw_toa':lw_crf_toa_cs, 
                     'sw_surf':sw_crf_surf_cs, 'lw_surf':lw_crf_surf_cs}

        return cs_forcing

    def total_atmos_forcing(self, data_all_sky, data_clear_sky, data_budget):    

        # get the cloud forcing
        cre_output = self.compute_cre(data_all_sky, data_clear_sky)

        #get the clear sky forcing
        cs_forcing =  self.atm_cs_forcing(data_all_sky, data_clear_sky)

        #atmos sw forcing for clouds, clear sky and all sky
        atm_sw_crf_cld = cre_output['swcre'] - cre_output['swcre_surf']
        atm_sw_crf_cs  = cs_forcing['sw_toa'] - cs_forcing['sw_surf']
        atm_sw_crf     = atm_sw_crf_cld + atm_sw_crf_cs

        #atmos lw forcing for clouds, clear sky and all sky
        atm_lw_crf_cld = cre_output['lwcre'] - cre_output['lwcre_surf']
        atm_lw_crf_cs  = cs_forcing['lw_toa'] - cs_forcing['lw_surf']
        atm_lw_crf     = atm_lw_crf_cld + atm_lw_crf_cs 

        total_forcing = ( atm_lw_crf +  atm_sw_crf +
                          data_budget['p'] + data_budget['sh'])

        forcing = {'sw_crf_cld':atm_sw_crf_cld, 'sw_crf_cs':atm_sw_crf_cs, 
                   'sw_crf_all':atm_sw_crf, 'lw_crf_cld':atm_lw_crf_cld, 
                   'lw_crf_cs':atm_lw_crf_cs, 'lw_crf_all':atm_lw_crf,
                   'total': total_forcing}

        return forcing


    def global_avg_flux_comp(self, data, lat):

        global_flux_comp = {}
        for var in ['lwut', 'lwus', 'lwds', 'swut', 'swdt', 'swus', 'swds']:
            global_flux_comp[var] = calc_global_mean(data[var], lat)

        return global_flux_comp

