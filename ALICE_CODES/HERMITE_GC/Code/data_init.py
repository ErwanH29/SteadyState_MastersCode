from evol_func import *
from physics_func import *
from file_logistics import file_counter
from amuse.ext.LagrangianRadii import LagrangianRadii
import pandas as pd
import os

class data_initialiser(object):
    def chaotic_sim_tracker(self, pset, init_pop, Nmerge, cum_mergermass, time, ejected_key, 
                            stab_time, added_mass, ejected_mass, comp_time, int_string):

        """
        Data set which tracks statistics on the simulations resulting in ejections
        
        Inputs:
        pset:           The particle set
        init_pop:       The initial population of the simulation
        Nmerge:         The number of mergers encountered since last particle added
        cum_mergermass: The cumulative mass merged
        time:           The simulation end time
        ejected_key:    The key of the ejected particle
        stab_time:      The length the system of N particle was stable for
        added_mass:     The mass of the most recently added particle
        ejected_mass:   The mass of the ejected particle
        comp_time:      Total time simulation lasted
        int_string:     String describing integrator used
        """

        count = file_counter(int_string)
        path = '/home/s2009269/data1/HERMITE_GC_data/no_addition/chaotic_simulation'
        gc_code = globular_cluster()
        stab_tracker = pd.DataFrame()
        df_stabtime = pd.Series({'Initial Particles': (len(init_pop)-2), 'Final Particles': (len(pset)-2), 
                                 'Number of Mergers': Nmerge, 'Cumulative Merger Mass': cum_mergermass.in_(units.MSun),
                                 'Simulated Till': time.in_(units.yr),
                                 'Ejected Particle': ejected_key, 'Stability Time': stab_time, 
                                 'Initial Distance': gc_code.gc_dist.in_(units.parsec), 
                                 'Cluster Radius': gc_code.gc_rad.in_(units.pc),
                                 'Initial Particle Mass': init_pop[2:].mass.in_(units.MSun),
                                 'Added Particle Mass': added_mass.in_(units.MSun),
                                 'Ejected Mass': ejected_mass.in_(units.MSun),
                                 'Computation Time': str(comp_time),
                                 'Relaxation Time': relax_timescale(gc_code.gc_rad, gc_code.gc_mass, 10**5).in_(units.yr)})
        stab_tracker = stab_tracker.append(df_stabtime, ignore_index = True)
        stab_tracker.to_pickle(os.path.join(path, 'IMBH_'+str(int_string)+'_Local_'+str(len(init_pop)-2)+'_sim'+str(count)+'_init_dist'
                               +str('{:.3f}'.format(gc_code.gc_dist.value_in(units.parsec)))+'_equal_mass_'+str('{:.3f}'.format(init_pop[2].mass.value_in(units.MSun)))+'.pkl'))
             
    def coll_tracker(self):
        coll_tracker = pd.DataFrame()
        df_coll_tracker = pd.Series({'Collision Time': 0 | units.s, 'Collided Particles': [0, 0], 
                                     'Initial Mass': [0, 0] | units.MSun, 'Emergent Particle': 0, 
                                     'Collision Mass': 0 | units.MSun})
        coll_tracker = coll_tracker.append(df_coll_tracker, ignore_index = True)
        return coll_tracker

    def energy_tracker(self, E0, time, app_time):
        """
        Data set to track the energy evolution of the system
        
        Inputs:
        E0:       The initial energy
        time:     The initial time
        app_time: The time a new particle appears
        """

        energy_tracker = pd.DataFrame()
        df_energy_tracker = pd.Series({'Time': time.in_(units.kyr), 'Et': E0 , 'dE': 0, 'dEs': 0, 
                                       'Appearance': app_time, 'Collision Time': 0 | units.s, 
                                       'Collision Mass': 0 | units.MSun })
        energy_tracker = energy_tracker.append(df_energy_tracker, ignore_index=True)

        return energy_tracker

    def event_tracker(self, pset):
        eventstab_tracker = pd.DataFrame()
        df_eventstab_tracker = pd.Series({'Initial No. Particles': len(pset)-2, 'Merger Event': 0, 
                                          'Collision Time': 0 | units.s, 'Injected Event': 0, 
                                          'Injected Time': 0 | units.s})
        eventstab_tracker = eventstab_tracker.append(df_eventstab_tracker, ignore_index = True)
        return eventstab_tracker

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
            if i == 0 :
                df_IMBH_vals = pd.Series({#'key_tracker': pset[i].key_tracker, 
                                          '{}'.format(time): [pset[i].key_tracker, pset[i].mass, pset[i].position, pset[i].velocity, 
                                          0 | units.J, 0 | units.J, tdyn_val[i] | units.yr]})
                df_IMBH = df_IMBH.append(df_IMBH_vals, ignore_index=True)
            else:
                gcframe_vel = pset[i].velocity - pset[1].velocity
                parti_KE = 0.5*pset[i].mass*pset[i].velocity.length()**2
                temp_PE = []
                temp_PE = indiv_PE_all(pset[i], pset, temp_PE)
                parti_PE = max(temp_PE)
                df_IMBH_vals = pd.Series({'{}'.format(time): [pset[i].key_tracker, pset[i].mass, pset[i].position, gcframe_vel,
                                        parti_KE, parti_PE, tdyn_val[i] | units.yr]})
                df_IMBH = df_IMBH.append(df_IMBH_vals, ignore_index=True)
        IMBH_array = IMBH_array.append(df_IMBH, ignore_index=True)

        return IMBH_array

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

    def parti_energy_tracker(self, BE, KE, E0, time):
        """
        Data set which tracks the inidivual particle energy
        
        Inputs:
        BE:     The particles maximum binding energy
        KE:     The particles kinetic energy
        E0:     The particles total energy
        time:   The time of the
        """
        
        parti_energy_tracker = pd.DataFrame()
        df_parti_energy_tracker = pd.Series({'Time': time.in_(units.kyr), 'BE': BE.in_(units.J), 
                                             'KE': KE.in_(units.J), 'Total E': E0})
        parti_energy_tracker = parti_energy_tracker.append(df_parti_energy_tracker, ignore_index=True)

        return parti_energy_tracker

    def stable_sim_tracker(self, pset, Ninj, Nmerge, merger_mass, time, int_string):

        """
        Function which tracks information on the simulations which ended in stable state.
        This occurs if time = tend, new particle is injected or a merging event occurs.
        
        Inputs:
        pset:        The particle set
        Ninj:        Binary stating if simulation ended due to injection (1 = injection, 0 = none)
        Nmerge:      Binary stating if simulation ended due to merging event (1 = merge, 0 = none)
        merger_mass: The merger event
        time:        Duration of the stable system
        int_string:  String describing integrator used
        """

        gc_code = globular_cluster()
        count = file_counter(int_string)

        stable_sim_tracker = pd.DataFrame()
        df_stablesim = pd.Series({'Initial Particles': (len(pset)-2), 'Injected Event': Ninj,
                                  'Merged Event': Nmerge, 'Merger Mass': merger_mass,
                                  'Simulated Till/Stab Time': time.in_(units.yr),
                                  'Initial Distance': gc_code.gc_dist.in_(units.parsec), 
                                  'Cluster Radius': gc_code.gc_rad.in_(units.pc),
                                  'Initial Particle Mass': pset[2:].mass.in_(units.MSun),
                                  'Relaxation Time': relax_timescale(gc_code.gc_rad, gc_code.gc_mass, 10**5).in_(units.yr)})
        stable_sim_tracker = stable_sim_tracker.append(df_stablesim, ignore_index = True)
        stable_sim_tracker.to_pickle('data/'+str(int_string)+'/stable_simulation/IMBH_'+str(int_string)+'_Ni'+str(len(pset)-2)+'_sim'+str(count)
                                     +'_init_dist' +str('{:.3f}'.format(gc_code.gc_dist.value_in(units.parsec)))
                                     +'_equal_mass_' +str('{:.3f}'.format(pset[2].mass.value_in(units.MSun)))
                                     +'_inj_'+str(Ninj)+'_merge_'+str(Nmerge)+'.pkl')