from evol_func import *
from physics_func import *
from file_logistics import file_counter
from amuse.ext.LagrangianRadii import LagrangianRadii
import pandas as pd
import os

class data_initialiser(object):
    def chaotic_sim_tracker(self, pset, init_pop, Nmerge, cum_mergermass, time, ejected_key, 
                            stab_time, added_mass, ejected_mass, comp_time, ejected, int_string):

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
        ejected:        0 for not ejected (sim. ended with merger) or 1 for ejection away from SMBH
        int_string:     String describing integrator used
        """

        SMBH_code = MW_SMBH()
        count = file_counter(int_string)
        path = '/home/s2009269/data1/HERMITE_Orbit_data/'
        stab_tracker = pd.DataFrame()
        df_stabtime = pd.Series({'Initial Particles': (len(init_pop)-1), 'Final Particles': (len(pset)-1), 
                                 'Number of Mergers': Nmerge, 'Cumulative Merger Mass': cum_mergermass.in_(units.MSun),
                                 'Ejection': ejected, 'Simulated Till': time.in_(units.yr),
                                 'Ejected Particle': ejected_key, 'Stability Time': stab_time, 
                                 'Initial Distance': SMBH_code.distance.in_(units.parsec),
                                 'Initial Particle Mass': init_pop[2:].mass.in_(units.MSun),
                                 'Added Particle Mass': added_mass.in_(units.MSun),
                                 'Ejected Mass': ejected_mass.in_(units.MSun),
                                 'Computation Time': str(comp_time)})
        stab_tracker = stab_tracker.append(df_stabtime, ignore_index = True)
        stab_tracker.to_pickle(os.path.join(path+str('/no_addition/chaotic_simulation'), 'IMBH_'+str(int_string)+'_Local_'+str(init_pop)
                                            +'_sim'+str(count)+'_init_dist'+str('{:.3f}'.format(SMBH_code.distance.value_in(units.parsec)))
                                            +'_equal_mass_'+str('{:.3f}'.format(pset[2].mass.value_in(units.MSun)))+'.pkl'))

    def energy_tracker(self, E0, Ek, Ep, time, app_time):
        """
        Data set to track the energy evolution of the system
        
        Inputs:
        E0:       The total initial energy of the particles
        Ek/Ep:    The total kinetic/pot. energy of the particles
        time:     The initial time
        app_time: The time a new particle appears
        """

        energy_tracker = pd.DataFrame()
        df_energy_tracker = pd.Series({'Time': time.in_(units.kyr), 'Et': E0 , 'dE': 0, 'dEs': 0,
                                       'Kinetic E': Ek.in_(units.J), 'Pot. E': Ep.in_(units.J), 
                                       'Appearance': app_time, 'Collision Time': 0 | units.s, 
                                       'Collision Mass': 0 | units.MSun })
        energy_tracker = energy_tracker.append(df_energy_tracker, ignore_index=True)

        return energy_tracker

    def IMBH_tracker(self, pset, time, init_pop):
        """
        Data set which holds information on each individual particles
        
        Inputs:
        pset:        The particle set
        time:        The initial time of the simulation
        init_pop:    The initial population of the simulation
        """

        IMBH_array = pd.DataFrame()
        df_IMBH    = pd.DataFrame()
        for i in range(init_pop):
            semimajor = []
            eccentric = []
            inclinate = []
            arg_peri  = []
            asc_node  = []
            true_anom = []

            if i == 0 :
                df_IMBH_vals = pd.Series({#'key_tracker': pset[i].key_tracker, 
                                          '{}'.format(time): [pset[i].key_tracker, pset[i].mass, pset[i].position, pset[i].velocity, 
                                                              0 | units.J, 0 | units.J, [0,0], [0,0], [0,0], [0,0], 
                                                             [0,0], [0,0]]})
                df_IMBH = df_IMBH.append(df_IMBH_vals, ignore_index=True)

            else:
                bin_sys = Particles()
                bin_sys.add_particle(pset[i])
                bin_sys.add_particle(pset[0])
                kepler_elements = orbital_elements_from_binary(bin_sys, G=constants.G)
                semimajor.append(kepler_elements[2].value_in(units.parsec))
                eccentric.append(kepler_elements[3])
                inclinate.append(kepler_elements[4])
                arg_peri.append(kepler_elements[5])
                asc_node.append(kepler_elements[6])
                true_anom.append(kepler_elements[7])
                
                neighbour_dist, nearest_parti = nearest_neighbour(pset[i], pset)
                bin_sys = Particles()
                bin_sys.add_particle(pset[i])
                bin_sys.add_particle(nearest_parti)
                kepler_elements = orbital_elements_from_binary(bin_sys, G=constants.G)
                semimajor.append(kepler_elements[2].value_in(units.parsec))
                eccentric.append(kepler_elements[3])
                inclinate.append(kepler_elements[4])
                arg_peri.append(kepler_elements[5])
                asc_node.append(kepler_elements[6])
                true_anom.append(kepler_elements[7])

                parti_KE = 0.5*pset[i].mass*pset[i].velocity.length()**2
                temp_PE = indiv_PE_all(pset[i], pset)
                parti_PE = np.sum(temp_PE)
                df_IMBH_vals = pd.Series({'{}'.format(time): [pset[i].key_tracker, pset[i].mass, pset[i].position, 
                                                              pset[i].velocity, parti_KE, parti_PE, semimajor * 1 | units.parsec, 
                                                              eccentric, inclinate, arg_peri, asc_node, true_anom, neighbour_dist]})
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

    def stable_sim_tracker(self, pset, Ninj, Nmerge, merger_mass, time, int_string, deltaE):

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
        deltaE:      Change in energy
        """

        count = file_counter(int_string)
        SMBH_code = MW_SMBH()

        stable_sim_tracker = pd.DataFrame()
        df_stablesim = pd.Series({'Initial Particles': (len(pset)-1), 'Injected Event': Ninj,
                                  'Merged Event': Nmerge, 'Merger Mass': merger_mass,
                                  'Simulated Till/Stab Time': time.in_(units.yr),
                                  'Initial Distance': SMBH_code.distance.in_(units.parsec),
                                  'Initial Particle Mass': pset[2:].mass.in_(units.MSun),
                                  'Change in Energy': deltaE})
        stable_sim_tracker = stable_sim_tracker.append(df_stablesim, ignore_index = True)
        stable_sim_tracker.to_pickle('data/'+str(int_string)+'/stable_simulation/IMBH_'+str(int_string)+'_Ni'+str(len(pset)-1)+'_sim'+str(count)
                                     +'_init_dist' +str('{:.3f}'.format(SMBH_code.distance.value_in(units.parsec)))
                                     +'_equal_mass_' +str('{:.3f}'.format(pset[2].mass.value_in(units.MSun)))
                                     +'_inj_'+str(Ninj)+'_merge_'+str(Nmerge)+'.pkl')