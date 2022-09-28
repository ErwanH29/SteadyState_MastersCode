from evol_func import *
from physics_func import *
from amuse.ext.LagrangianRadii import LagrangianRadii
import pandas as pd

class data_initialiser(object):

    def IMBH_tracker(self, pset, time, init_pop):
        """
        Data set which holds information on each individual particles
        
        Inputs:
        pset:     The particle set
        time:     The initial time of the simulation
        init_pop: The initial population of the simulation
        """

        tdyn_val = tdyn_calc(pset)

        IMBH_array = pd.DataFrame()
        df_IMBH    = pd.DataFrame()
        for i in range(init_pop):
            parti_KE = 0.5*pset[i].mass*pset[i].velocity.length()**2
            temp_PE = []
            temp_PE = indiv_PE_BH(pset[i], pset, temp_PE)
            parti_PE = max(temp_PE)
            df_IMBH_vals = pd.Series({'key_tracker': pset[i].key_tracker, 
                                      '{}'.format(time): [pset[i].mass, pset[i].position, 
                                      parti_KE, parti_PE, tdyn_val[i] | units.yr]})
            df_IMBH = df_IMBH.append(df_IMBH_vals, ignore_index=True)
        IMBH_array = IMBH_array.append(df_IMBH, ignore_index=True)

        return IMBH_array

    def com_tracker(self, gravity):
        """
        Data set which holds information on the cluster c.o.m

        Inputs:
        gravity: The integrator used in the simulation
        """

        com_tracker = pd.DataFrame()
        df_com_tracker = pd.Series({'x': SMBH_filter(gravity.particles).center_of_mass()[0].in_(units.parsec),
                                    'y': SMBH_filter(gravity.particles).center_of_mass()[1].in_(units.parsec),
                                    'z': SMBH_filter(gravity.particles).center_of_mass()[2].in_(units.parsec)})
        com_tracker = com_tracker.append(df_com_tracker, ignore_index=True)

        return com_tracker

    def LG_tracker(self, clust_rad, clust_mass, no_stars, pset, time, gravity):
        """
        Data set which tracks the Lagrangian radius and tidal radius of the cluster.
        
        Inputs:
        clust_rad:  The cluster radius
        clust_mass: The cluster mass
        no_stars:   The number of stars in the cluster
        pset:       The particle set
        time:       The initial time of the simulation
        gravity:    The integrator used for the simulation
        """

        LG_array = pd.DataFrame()
        df_LG_tracker = pd.Series({'Time': time.in_(units.kyr),
                                   'LG25': LagrangianRadii(gravity.particles[1:])[5].in_(units.parsec),
                                   'LG50': LagrangianRadii(gravity.particles[1:])[6].in_(units.parsec),
                                   'LG75': LagrangianRadii(gravity.particles[1:])[7].in_(units.parsec),
                                   'Tidal Radius': tidal_radius(pset).in_(units.parsec),
                                   'Relaxation Time': relax_timescale(clust_rad, clust_mass, no_stars).in_(units.yr)})
        LG_array = LG_array.append(df_LG_tracker , ignore_index=True)

        return LG_array

    
    def energy_tracker(self, E0, time, app_time):
        energy_tracker = pd.DataFrame()
        df_energy_tracker = pd.Series({'Time': time.in_(units.kyr), 'Et': E0 , 'dE': 0, 'dEs': 0, 'Appearance': app_time, 
                                    'Collision Time': 0 | units.s, 'Collision Mass': 0 | units.MSun })
        energy_tracker = energy_tracker.append(df_energy_tracker, ignore_index=True)
        return energy_tracker

    def parti_energy_tracker(self, BE, KE, E0, time):
        parti_energy_tracker = pd.DataFrame()
        df_parti_energy_tracker = pd.Series({'Time': time.in_(units.kyr), "BE": BE.in_(units.J), 
                                            "KE": KE.in_(units.J), "Total E": E0})
        parti_energy_tracker = parti_energy_tracker.append(df_parti_energy_tracker, ignore_index=True)
        return parti_energy_tracker

    def coll_tracker(self):
        coll_tracker = pd.DataFrame()
        df_coll_tracker = pd.Series({'Collision Time': 0 | units.s, 'Collided Particles': [0, 0], 'Initial Mass': [0, 0] | units.MSun,
                                    'Emergent Particle': 0, 'Collision Mass': 0 | units.MSun})
        coll_tracker = coll_tracker.append(df_coll_tracker, ignore_index = True)
        return coll_tracker

    def event_tracker(self, pset):
        eventstab_tracker = pd.DataFrame()
        df_eventstab_tracker = pd.Series({'Initial No. Particles': len(pset)-1, 'Merger Event': 0, 
                                        'Collision Time': 0 | units.s, 'Injected Event': 0, 'Injected Time': 0 |units.s})
        eventstab_tracker = eventstab_tracker.append(df_eventstab_tracker, ignore_index = True)
        return eventstab_tracker

