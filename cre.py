from .global_mean import calc_global_mean

__author__ = 'Penelope Maher'

class ComputeCloudRadiativeEffect():

    def __init__(self, data):
        """ Initialize the labels for calculation.

        Input
        ----------------
            data: is a dictionary that contain all-sky and clear-sky fluxes
                  the keys names are assumed to be 
                  ['swut', 'swdt', 'swus', 'swds', 'lwut','lwus','lwds',
                   'swut_cs', 'swdt_cs', 'swus_cs', 'swds_cs', 'lwut_cs',
                   'lwus_cs','lwds_cs']

        List of acronyms
        -----------------
            lw: longwave
            sw: shortwave
           u,d: up and down (direction of fluxes to make sense of sign)
           t,s: top of atmosphere and surface
           cre: cloud radiative effect (at top of atmosphere)
      cre_surf: cloud radiative effect (at the surface)
          acre: cloud radiative effect (of the atmospheric column)
            cs: clear sky fluxes
           all: all-sky fluxes (where all-sky = cloudy-sky + clear-sky)

        Units
        --------------
        All input data should be fluxes. All units are W/m2.
    
        List of papers
        --------------

            CRE equation:
                Allan 2011 "Combining satellite data and models to estimate cloud
                            radiative effect at the surface and in the atmosphere"  

        """

        self.data = data

        #cre terms that are calculated in this class
        self.lwcre = None
        self.swcre = None
        self.cre   = None

        self.lwcre_surf = None
        self.swcre_surf = None
        self.cre_surf   = None

        self.alwcre = None
        self.aswcre = None
        self.acre   = None

        #to keep track of what variables there are
        self.var_list = ['cre', 'lwcre', 'swcre',
                         'cre_surf', 'lwcre_surf', 'swcre_surf',
                         'acre', 'alwcre', 'aswcre']

    def compute_cre(self):
        """
        Compute the cloud radiative effect

        Parameters needed
        ----------
            self.data : Dict of all-sky and clear-sky fluxes. 
                        Fluxes are assumed to be 1D and have the shape 
                        self.data['flux_name'][lat], i.e. time-mean zonal-mean.
                        But the code will work on other shapped data.
 
        Returns
        -------
                          : no returns

        Balance values:
        -------

            IPCC AR5 Chapter 7 p 580 global average (TOA)
                SWCRE = -50 W/m2
                LWCRE = +30 W/m2
                CRE   = -20 W/m2 (ie the net effect of clouds is a cooling)

        Notes:
        -------
          - cre, cre_surf, and acre can be validated against Allan 2011 Fig 5.

        """

        #TOA CRE all sky
        self.lwcre = self.data['lwut_cs'] - self.data['lwut']
        self.swcre = self.data['swut_cs'] - self.data['swut']
        self.cre  = self.lwcre + self.swcre 

        #surf CRE all sky
        self.lwcre_surf = (self.data['lwds'] - self.data['lwds_cs'] 
                          -self.data['lwus'] + self.data['lwus_cs']) 
        self.swcre_surf = (self.data['swds'] - self.data['swds_cs'] 
                          -self.data['swus'] + self.data['swus_cs'])
        self.cre_surf   = self.lwcre_surf + self.swcre_surf

        #atmospheric CRE all sky 
        self.alwcre = self.lwcre - self.lwcre_surf
        self.aswcre = self.swcre - self.swcre_surf
        self.acre   = self.alwcre + self.aswcre

    def global_avg_cre(self, lat):
        # returns the global mean 
        global_cre = {}
        for var in self.cre.keys():
            global_cre[var] = calc_global_mean(self.cre[var], lat)    
            print('Global CRE for {} is: {:8.2f}'.format(var,global_cre[var]))
        
        return global_cre


