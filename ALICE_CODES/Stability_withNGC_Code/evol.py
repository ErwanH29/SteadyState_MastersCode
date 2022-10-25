from amuse.lab import *
from amuse.units import units
from parti_initialiser import *
from data_init import *
from physics_func import *
from evol_func import *
from amuse.ext.LagrangianRadii import LagrangianRadii
from amuse.ext.orbital_elements import orbital_elements_from_binary
from amuse.community.hermite import Hermite
#from amuse.community.hermite_grx.interface import *
import numpy as np
import pandas as pd
import time as cpu_time

# FOR ALICE SIMULATIONS:
#   - REMOVE ALL PRINT STATEMENTS + EXCESS FILE (ONLY KEEP CHAOTIC + STABLE + SIMULATION FILES)
#   - CHANGE number_of_workers
#   - ENSURE CORRECT FILE SAVING FORMAT (NAME)

def evolve_system(parti, tend, eta, init_dist, converter, int_string):
    """
    Bulk of the simulation. Uses Hermite integrator to evolve the system while bridging it.
    Keeps track of the particles' position and stores it in a pickle file.
    
    Inputs:
    parti:       The particle set needed to simulate
    tend:        The end time of the simulation
    eta:         The step size
    init_dist:   The initial distance between IMBH and SMBH
    converter:   Variable used to convert between nbody units and SI
    int_string:  String to dictate whether using Hermite or Hermite GRX
    """
    
    set_printing_strategy("custom", preferred_units = [units.MSun, units.pc, units.kyr, units.AU/units.yr, units.J],
                                                       precision = 20, prefix = "", separator = " ", suffix = "")

    comp_start = cpu_time.time()
    #SMBH_code = MW_SMBH()
    #IMBH_adder = IMBH_init()

    initial_set = parti.copy()

    code = Hermite(converter, number_of_workers = 6)
    #code.parameters.light_speed = constants.c
    code.particles.add_particles(parti)
    code.commit_particles()
    stopping_condition = code.stopping_conditions.collision_detection
    stopping_condition.enable()
    
    channel_IMBH = {"from_gravity": 
                    code.particles.new_channel_to(parti,
                    attributes=["x", "y", "z", "vx", "vy", "vz", "mass"],
                    target_names=["x", "y", "z", "vx", "vy", "vz", "mass"]),
                    "to_gravity": 
                    parti.new_channel_to(code.particles,
                    attributes=["mass", "collision_radius"],
                    target_names=["mass", "radius"])} 
      
    time = 0 | units.yr
    iter = 0
    add_iter = 0
    Nenc = 0
    app_time2 = time #To activate when allowing ejections to occur in a continuous simulation
    cum_merger_mass = 0 | units.MSun
    stab_timescale = time

    N_parti = len(parti)
    init_IMBH = N_parti-1
    ejected = False
    extra_note = ''

    parti_KE = code.kinetic_energy
    parti_BE = code.potential_energy
    E0p = parti_KE + parti_BE 
    E0 = E0p

    data_trackers = data_initialiser()
    energy_tracker = data_trackers.energy_tracker(E0p, parti_KE, parti_BE, time, 0 | units.s)
    IMBH_tracker = data_trackers.IMBH_tracker(parti, time, N_parti)
    #LG_tracker = data_trackers.LG_tracker(cluster_radi, IMBH_adder.mass, gc_code.gc_pop, parti, time, code)
    bounded_SMBH = 0
    
    while time < tend:
        eject  = 0
        merger = 0
        iter += 1
        app_time = 0 | units.s

        merger_mass = 0 | units.MSun
        added_mass = 0 | units.MSun
        ejected_mass = 0 | units.MSun
        tcoll = 0 | units.yr

        if iter %100 == 0:
            print('Iteration', iter)

        time += eta*tend
        channel_IMBH["to_gravity"].copy()
        code.evolve_model(time)

        for particle in SMBH_filter(parti):
            rel_vel = particle.velocity - parti[0].velocity
            dist_core = particle.position - parti[0].position

            dist_vect = np.sqrt(np.dot(dist_core, dist_core))
            vel_vect  = np.sqrt(np.dot(rel_vel, rel_vel))
            curr_traj = (np.dot(dist_core, rel_vel))/(dist_vect * vel_vect) #Movement towards SMBH

            parti_KE = 0.5*particle.mass*(particle.velocity.length())**2
            temp_PE = indiv_PE_all(particle, parti)
            parti_BE = np.sum(temp_PE)

            if parti_KE > abs(parti_BE) and particle.position.length() > 0.4 | units.parsec and curr_traj > 0:
                eject = 1
                ejected_key_track = particle.key_tracker
                ejected_mass = particle.mass
                extra_note = 'Stopped due to particle ejection'
                break

        injbin = 0
        if ejected == 0:
            if stopping_condition.is_set():
                print("........Encounter Detected........")
                print('Collision at step: ', iter)
                print('Simulation will now stop')
                energy_before = code.potential_energy + code.kinetic_energy
                merge = 1
                for ci in range(len(stopping_condition.particles(0))):
                    Nenc += 1

                    enc_particles_set = Particles(particles=[stopping_condition.particles(0)[ci],
                                                         stopping_condition.particles(1)[ci]])
                    enc_particles = enc_particles_set.get_intersecting_subset_in(parti)
                    merged_parti = merge_IMBH(parti, enc_particles, code.model_time)
                    merger_mass = merged_parti.mass.sum()
                    cum_merger_mass += merger_mass
                    parti.synchronize_to(code.particles)

                    tcoll = time.in_(units.s) - eta*tend
                    energy_after = code.potential_energy + code.kinetic_energy
                    deltaE = abs(energy_after-energy_before)/abs(energy_before)
                    data_trackers.stable_sim_tracker(parti, injbin, merge, merger_mass, stab_timescale, int_string, deltaE) 

                    df_eventstab_tracker = pd.Series({'Initial No. Particles': len(parti)-1, 'Merger Event': 1, 
                                                      'Collision Time': tcoll.in_(units.kyr), 'Injected Event': 0, 'Injected Time': 0 |units.s, 'Bounded to SMBH': bounded_SMBH})
                    eventstab_tracker = eventstab_tracker.append(df_eventstab_tracker, ignore_index = True)

                    ejected_key_track = parti[-1].key_tracker
                    df_coll_tracker = pd.Series({'Collision Time': tcoll.in_(units.kyr), 
                                                 'Collided Particles': [enc_particles_set[0].key, enc_particles_set[1].key], 
                                                 'Initial Mass': [enc_particles_set[0].mass, enc_particles_set[1].mass] | units.MSun, 
                                                 'Emergent Particle': ejected_key_track, 'Final Mass': merger_mass | units.MSun})                    
                    coll_tracker = coll_tracker.append(df_coll_tracker, ignore_index = True) 
                    extra_note = 'Stopped due to merger'

        channel_IMBH["from_gravity"].copy()     
        rows = (len(parti))

        df_IMBH = pd.DataFrame()
        for i in range(len(parti)):
            for j in range(len(parti)):
                if IMBH_tracker.iloc[i][0][0] == parti[j].key_tracker:
                    neighbour_dist, nearest_parti = nearest_neighbour(parti[j], parti)
                    semimajor = []
                    eccentric = []
                    inclinate = []
                    arg_peri  = []
                    asc_node  = []
                    true_anom = []

                    if j == 0:
                        semimajor = [0, 0]
                        eccentric = [0, 0]
                        inclinate = [0, 0]
                        arg_peri  = [0, 0]
                        asc_node  = [0, 0]
                        true_anom = [0, 0]
                    else:
                        bin_sys = Particles()                # First elements of the orbital arrays will be parti + SMBH
                        bin_sys.add_particle(parti[j])
                        bin_sys.add_particle(parti[0])
                        kepler_elements = orbital_elements_from_binary(bin_sys, G=constants.G)
                        semimajor.append(kepler_elements[2].value_in(units.parsec))
                        eccentric.append(kepler_elements[3])
                        inclinate.append(kepler_elements[4])
                        arg_peri.append(kepler_elements[5])
                        asc_node.append(kepler_elements[6])
                        true_anom.append(kepler_elements[7])

                        bin_sys = Particles()                 # Second elements correspond to IMBH + IMBH
                        bin_sys.add_particle(parti[j])
                        bin_sys.add_particle(nearest_parti)
                        kepler_elements = orbital_elements_from_binary(bin_sys, G=constants.G)
                        semimajor.append(kepler_elements[2].value_in(units.parsec))
                        eccentric.append(kepler_elements[3])
                        inclinate.append(kepler_elements[4])
                        arg_peri.append(kepler_elements[5])
                        asc_node.append(kepler_elements[6])
                        true_anom.append(kepler_elements[7])
                        
                    parti_KE = 0.5*parti[j].mass*parti[j].velocity.length()**2
                    temp_PE = indiv_PE_all(parti[j], parti)
                    parti_PE = np.sum(temp_PE)
                    df_IMBH_vals = pd.Series({'{}'.format(time): [parti[j].key_tracker, parti[j].mass, parti[j].position,
                                                                  parti[j].velocity, parti_KE, parti_PE, semimajor * 1 | units.parsec, 
                                                                  eccentric, inclinate, arg_peri, asc_node, true_anom, neighbour_dist]})
                    break
                else:
                    df_IMBH_vals = pd.Series({'{}'.format(time): [np.NaN, np.NaN | units.MSun,
                                                                 [np.NaN | units.parsec, np.NaN | units.parsec, np.NaN | units.parsec],
                                                                 [np.NaN | units.kms, np.NaN | units.kms, np.NaN | units.kms],
                                                                  np.NaN | units.J, np.NaN | units.J, np.NaN | units.m, np.NaN,
                                                                  np.NaN, np.NaN, np.NaN , np.NaN, np.NaN]})
            df_IMBH = df_IMBH.append(df_IMBH_vals, ignore_index=True)

        IMBH_tracker = IMBH_tracker.append(df_IMBH, ignore_index=True)
        IMBH_tracker['{}'.format(time)] = IMBH_tracker['{}'.format(time)].shift(-rows)
        IMBH_tracker = IMBH_tracker.dropna(axis='rows')

        parti_KE = code.kinetic_energy
        parti_BE = code.potential_energy
        Etp  = parti_KE + parti_BE
        Et = Etp
        de = abs(Et-E0)/abs(E0)
        if 20 < iter:
            dEs = abs(Et-energy_tracker.iloc[19][1])/abs(energy_tracker.iloc[19][1])
            df_energy_tracker = pd.Series({'Time': time.in_(units.kyr), 'Et': Et, 'dE': de, 'dEs': dEs, 
                                           'Kinetic E': parti_KE.in_(units.J), 'Pot.E': parti_BE.in_(units.J),
                                           'Appearance': app_time, 'Collision Time': tcoll, 
                                           'Collision Mass': merger_mass})
        else:
            df_energy_tracker = pd.Series({'Time': time.in_(units.kyr), 'Et': Et, 'dE': de, 'dEs': 0,
                                           'Kinetic E': parti_KE.in_(units.J), 'Pot.E': parti_BE.in_(units.J),
                                           'Appearance': app_time, 'Collision Time': tcoll, 
                                           'Collision Mass': merger_mass })
        energy_tracker = energy_tracker.append(df_energy_tracker, ignore_index=True)

        time1 = time
        if ejected > 0 or Nenc > 0:
            time = tend   
            
    code.stop()
    comp_end = cpu_time.time()
    comp_time = comp_end-comp_start 
    
    if iter > 5 + add_iter and iter > 10:
        print('Saving Data')
        no_plot = False
        #chaos_stab_timescale = time1 - app_time2
        chaos_stab_timescale = time1
        count = file_counter(int_string)

        if time1 == tend:
            ejected_key_track = parti[1].key_tracker
       
        path = '/home/s2009269/data1/HERMITE_Orbit_data/'
        file_names = 'IMBH_'+str(int_string)+'_Local_'+str(init_IMBH)+'_sim'+str(count)+'_init_dist'+str('{:.3f}'.format(init_dist.value_in(units.parsec)))+'_equal_mass_'+str('{:.3f}'.format(parti[2].mass.value_in(units.MSun)))+'.pkl'
        IMBH_tracker.to_pickle(os.path.join(path+str('particle_trajectory'), file_names))
        energy_tracker.to_pickle(os.path.join(path+str('energy'), file_names))
        coll_tracker.to_pickle(os.path.join(path+str('collision_events'), file_names))
        eventstab_tracker.to_pickle(os.path.join(path+str('event_tracker'), file_names))   #For different dependent variables [clust_dist, clust_rad, ALICE/local...], change output file name
        data_trackers.chaotic_sim_tracker(parti, initial_set, Nenc, cum_merger_mass, time1, ejected_key_track, 
                                          chaos_stab_timescale, added_mass, ejected_mass, comp_time, eject, int_string)                   #For different dependent variables [clust_dist, clust_rad, ALICE/local...], change output file name

        lines = ['Simulation: ', "Total CPU Time: "+str(comp_end-comp_start)+' seconds', 
                 'Timestep: '+str(eta),
                 'Simulated until: '+str(time1.value_in(units.yr))+str(' years'), 
                 'Cluster Distance: '+str(init_dist.value_in(units.parsec))+' parsecs', 
                 'Masses of IMBH: '+str(parti.mass.value_in(units.MSun))+' MSun',
                 "No. of initial IMBH: "+str(init_IMBH), 
                 'Number of new particles: '+str(len(parti)-1-init_IMBH),
                 'Total Number of (Final) IMBH: '+str(len(parti)-2), 
                # 'IMBH Appearance Rate: '+str(IMBHapp.value_in(units.yr))+' years',    #To add when continuous ejecting simulations occur
                 'Number of mergers: '+str(Nenc), 'End Time: '+str(tend.value_in(units.yr))+' years', 
                 'Integrator: Hermite (NO PN)',
                 'Extra Notes: ', extra_note]

        with open(os.path.join(path+str('simulation_stats'), 'simulation'+str(count)+'.txt'), 'w') as f:
            for line in lines:
                f.write(line)
                f.write('\n')
    else:
        no_plot = True
        print('...No stability timescale - simulation ended too quick...')

    return no_plot