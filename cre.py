from .global_mean import calc_global_mean

__author__ = 'Penelope Maher'

class ComputeCloudRadiativeEffect():

    def __init__(self, data_all_sky, data_cs):
        """ Initialize the labels for calculation.

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

        self.data_all = data_all_sky # a dictionary of all-sky flux arrays
        self.data_cs  = data_cs      # a dictionary of clear-sky flux arrays

        #cre terms
        self.lwcre = None
        self.swcre = None
        self.cre   = None

        self.lwcre_surf = None
        self.swcre_surf = None
        self.cre_surf   = None

        self.alwcre = None
        self.aswcre = None
        self.acre   = None

    def compute_cre(self):
        """
        Compute the cloud radiative effect

        Parameters needed
        ----------
            self.data_all : dict of all-sky fluxes. Fluxes are assumed to be 1D
                            and have the shape self.data['flux_name'][lat],
                            i.e. they are time-mean zonal-mean flux.
                            But the code will work on other shapped data.
 
            self.data_clear_sky : as above but for clear-sky fluxes

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
          - cre, cre_surf, and acre can be validated against Allan 2011 Fig 5. (see ref above)

        """

        #TOA CRE all sky
        self.lwcre = self.data_cs['lwut'] - self.data['lwut']
        self.swcre = self.data_cs['swut'] - self.data['swut']
        self.cre  = self.lwcre + self.swcre 

        #surf CRE all sky
        self.lwcre_surf = (self.data['lwds'] - self.data_cs['lwds'] 
                          -self.data['lwus'] + self.data_cs['lwus']) 
        self.swcre_surf = (self.data['swds'] - self.data_cs['swds'] 
                          -self.data['swus'] + self.data_cs['swus'])
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


