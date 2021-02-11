from __future__ import division
import numpy as np
import pdb

from .global_mean import calc_global_mean
from .cre import ComputeCloudRadiativeEffect

__author__ = 'Penelope Maher'


#constants
Lc =  2.50e6 #latent heat of condensation [J/kg]
Lf =  3.34e5#latent heat of fusion [J/kg]
sec_in_day = 60.*60.*24

class AtmosEnergyBudget(ComputeCloudRadiativeEffect):

    def __init__(self, data):
        """ Initialize the labels for the budget dictionary.

        Input
        -----------------
            data is a dictionary with the following keys:
                - all-sky fluxes: swut, swdt, swus, swds, lwut, lwus, lwds
                - sh
                - lh 
                - if lh is not provided: precipitation from rain P_r
                                         and snow P_s is needed. Both in mm/day

        List of acronyms
        -----------------
            lw: longwave
            sw: shortwave
           u,d: up and down (direction of fluxes to make sense of sign)
           t,s: top of atmosphere and surface
            cs: clear-sky fluxes (if not cs then fluxes are all-sky)
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

        super(AtmosEnergyBudget, self).__init__(data)

        self.data     = data # a dictionary of all-sky flux arrays
        #self.lh       = lh           # lh = (L_f x snow + L_c x precip) in W/m2

        #budget terms that are defined in this class
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

        #if 'lh' in self.data.keys():
        self.lh = self.data['lh']

        #elif ('rain' in self.data.keys()) and ('snow' in self.data.keys()):
        #rain and snow are in mm/day.
        #convert to kg/m2/3 then multiply by LH for W/m2

        self.lh_p = ((self.data['rain']/sec_in_day)*Lc +
                     (self.data['snow']/sec_in_day)*Lf)

        self.net = - self.lwc + self.swa +self.data['sh'] + self.lh
        self.net_p = - self.lwc + self.swa + self.data['sh'] + self.lh_p

    def atm_cs_forcing(self, data):
        # the clear sky fluxes are known (unlike the cloud fluxes above)
        # so the clear sky forcing is the net flux (i.e. down - up) 

        self.sw_crf_toa_cs  =   data['swdt_cs'] - data['swut_cs'] 
        self.lw_crf_surf_cs =   data['lwds_cs'] - data['lwus_cs']
        self.sw_crf_surf_cs =   data['swds_cs'] - data['swus_cs']
        self.lw_crf_toa_cs  = - data['lwut_cs']

        #cs_forcing = {'sw_toa':sw_crf_toa_cs,   'lw_toa':lw_crf_toa_cs, 
        #             'sw_surf':sw_crf_surf_cs, 'lw_surf':lw_crf_surf_cs}

        #return cs_forcing

    def total_atmos_forcing(self, data):    
        # there are all class attributes rather than variables passes out
        # as it is common to want to plot the components and this is a cleaner
        # way of passing out all the variables.

        # get the cloud forcing
        self.compute_cre()

        #get the clear sky forcing
        self.atm_cs_forcing(data)
        
        #atmos sw forcing for clouds, clear sky and all sky
        self.atm_sw_crf_cld = self.swcre - self.swcre_surf
        self.atm_sw_crf_cs  = self.sw_crf_toa_cs - self.sw_crf_surf_cs
        self.atm_sw_crf     = self.atm_sw_crf_cld + self.atm_sw_crf_cs

        #atmos lw forcing for clouds, clear sky and all sky
        self.atm_lw_crf_cld = self.lwcre - self.lwcre_surf
        self.atm_lw_crf_cs  = self.lw_crf_toa_cs - self.lw_crf_surf_cs
        self.atm_lw_crf     = self.atm_lw_crf_cld + self.atm_lw_crf_cs 

        self.total_forcing = ( self.atm_lw_crf +  self.atm_sw_crf +
                          (data['precip']/sec_in_day)*Lc + data['sh'])

        
        #forcing = {'sw_crf_cld':atm_sw_crf_cld, 'sw_crf_cs':atm_sw_crf_cs, 
        #           'sw_crf_all':atm_sw_crf, 'lw_crf_cld':atm_lw_crf_cld, 
        #           'lw_crf_cs':atm_lw_crf_cs, 'lw_crf_all':atm_lw_crf,
        #           'total': total_forcing}

        #return forcing


    def global_avg_flux_comp(self, data, lat):

        global_flux_comp = {}
        for var in ['lwut', 'lwus', 'lwds', 'swut', 'swdt', 'swus', 'swds']:
            global_flux_comp[var] = calc_global_mean(data[var], lat)

        return global_flux_comp

