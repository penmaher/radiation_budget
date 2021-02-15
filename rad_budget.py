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

        Code purpose
        --------------
            - self.compute_lwc() and compute_swa calculate the LWC and SWA
            - there are two methods for calculating the atmospheric energy budget
              Both methods are equivalent but use different was to express the
              LW and SW all-ky atmospheric forcing.  

        List of acronyms
        -----------------
            lw: longwave
            sw: shortwave
           u,d: up and down (direction of fluxes to make sense of sign)
           t,s: top of atmosphere and surface
            cs: clear-sky fluxes (if not cs then fluxes are all-sky)
            lh: latent heat
            sh: sensible heat
           lwc: longwave cooling (positive -> cooling)
           swa: shortwave absorption (positive -> heating)

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
        self.lwc   = None
        self.swa   = None
        self.lh    = None
        self.lh_p  = None  #LH computed using precip
        self.net   = None  #budget residual
        self.net_p = None  #budget residual computed using self.lh_p
         
        # Note that because this class inherits ComputeCloudRadiativeEffect, 
        # all of the CRE variables are also available in this class.

        #clear sky forcing
        self.sw_crf_toa_cs  = None
        self.lw_crf_surf_cs = None
        self.sw_crf_surf_cs = None
        self.lw_crf_toa_cs  = None

        #atmos sw forcing
        self.atm_sw_crf_cld = None
        self.atm_sw_crf_cs  = None
        self.atm_sw_crf     = None #swacre

        # atmos lw forcing
        self.atm_lw_crf_cld = None 
        self.atm_lw_crf_cs  = None
        self.atm_lw_crf     = None  

        self.total_forcing  = None # budget residual 
                         
        #to keep track of what variables there are
        self.var_list = ['lwc','swa','lh','lh_p','net','net_p',
                         'sw_crf_toa_cs','lw_crf_surf_cs','sw_crf_surf_cs','lw_crf_toa_cs',
                         'atm_sw_crf_cld','atm_sw_crf_cs','atm_sw_crf',
                         'atm_lw_crf_cld','atm_lw_crf_cs','atm_lw_crf',
                         'total_forcing']

        # The data are managed as class attributes rather than variables passed
        # out as it is common to plot the components and this is a cleaner
        # way of passing out all the variables.

    def compute_lwc(self):
        # net LW radiation lost to space
        self.lwc = self.data['lwut'] + (self.data['lwds'] - self.data['lwus'])      

    def compute_swa(self):
        # net sw absorbed by atmosphere
        self.swa = (self.data['swdt'] - self.data['swut']) - (self.data['swds'] - self.data['swus']) 

    def atmos_budget_lwc_swa(self):
        # Method 1 for computing the atmospheric energy budget: using LWC and SWA

        self.compute_lwc()
        self.compute_swa()

        #if 'lh' in self.data.keys():
        self.lh = self.data['lh']
        self.sh = self.data['sh']


        #elif ('rain' in self.data.keys()) and ('snow' in self.data.keys()):
        #rain and snow are in mm/day.
        #convert to kg/m2/3 then multiply by LH for W/m2

        self.lh_p = ((self.data['rain']/sec_in_day)*Lc +
                     (self.data['snow']/sec_in_day)*Lf)

        self.lh_p_only = (self.data['precip']/sec_in_day)*Lc


        self.net = - self.lwc + self.swa +self.sh + self.lh
        self.net_p = - self.lwc + self.swa + self.sh + self.lh_p
        self.net_p_only = - self.lwc + self.swa + self.sh + self.lh_p_only

    def atm_cs_forcing(self, data):
        # the clear sky fluxes are known (unlike the cloud fluxes above)
        # so the clear sky forcing is the net flux (i.e. down - up) 

        self.sw_crf_toa_cs  =   data['swdt_cs'] - data['swut_cs'] 
        self.lw_crf_surf_cs =   data['lwds_cs'] - data['lwus_cs']
        self.sw_crf_surf_cs =   data['swds_cs'] - data['swus_cs']
        self.lw_crf_toa_cs  = - data['lwut_cs']

    def total_atmos_forcing(self, data):    
        # Method 1 for computing the atmospheric energy budget: using 
        # the cloudy and clear-sky forcing.

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
   
        self.turb_flux      = self.lh_p + self.sh

        self.total_forcing =  self.atm_lw_crf + self.atm_sw_crf + self.turb_flux
        self.total_forcing_p_only =  self.atm_lw_crf + self.atm_sw_crf + self.lh_p_only + self.sh
        self.total_forcing_p =  self.atm_lw_crf + self.atm_sw_crf + self.lh_p + self.sh

    def global_avg_flux_comp(self, data, lat):

        global_flux_comp = {}
        for var in ['lwut', 'lwus', 'lwds', 'swut', 'swdt', 'swus', 'swds']:
            global_flux_comp[var] = calc_global_mean(data[var], lat)

        return global_flux_comp

