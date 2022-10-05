from amuse.lab import *
from amuse.units import units
from parti_initialiser import *
from data_init import *
from physics_func import *
from evol_func import *
from amuse.ext.LagrangianRadii import LagrangianRadii
from amuse.community.hermite import Hermite
import numpy as np
import pandas as pd
import time as cpu_time

def evolve_system(parti, tend, eta, cluster_distance, cluster_radi, cluster_mass, cluster_rhmass, converter):
    """
    Bulk of the simulation. Uses Hermite integrator to evolve the system while bridging it.
    Keeps track of the particles' position and stores it in a pickle file.
    
    Inputs:
    parti:            The particle set needed to simulate
    tend:             The end time of the simulation
    eta:              The step size
    cluster_distance: Initialised distance between cluster and SMBH
    cluster_radi:     Initialised cluster radius.
    converter:        Variable used to convert between nbody units and SI
    """

    set_printing_strategy("custom", preferred_units = [units.MSun, units.pc, units.kyr, units.AU/units.yr, units.J],
                                    precision = 20, prefix = "", separator = " ", suffix = "")

    comp_start = cpu_time.time()

    SMBH_code = MW_SMBH()
    gc_code = globular_cluster()
    IMBH_adder = IMBH_init()

    initial_set = parti.copy()

    code = Hermite(converter, number_of_workers = 20)
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
    Nenc = 0
    app_time2 = time
    cum_merger_mass = 0 | units.MSun
    stab_time_init = 0 | units.s

    init_IMBH_pop = len(parti)
    add_iter = 0
    N_parti = init_IMBH_pop
    N_parti_init = N_parti
    ejected = False
    extra_note = ''

    parti_KE = code.kinetic_energy
    parti_BE = code.potential_energy
    E0p = parti_KE + parti_BE 
    E0 = E0p

    data_trackers = data_initialiser()
    IMBH_tracker = data_trackers.IMBH_tracker(parti, time, init_IMBH_pop)
    com_tracker = data_trackers.com_tracker(code)
    LG_tracker = data_trackers.LG_tracker(cluster_radi, IMBH_adder.mass, gc_code.gc_pop, parti, time, code)
    energy_tracker = data_trackers.energy_tracker(E0p, time, 0 | units.s)
    parti_energy_tracker = data_trackers.parti_energy_tracker(parti_BE, parti_KE, E0p, time)
    coll_tracker = data_trackers.coll_tracker()
    eventstab_tracker = data_trackers.event_tracker(parti)

    """
    IMBHapp = df_timescale(cluster_rhmass)
    decision_scale = (eta*tend)/IMBHapp #
    print('\nNew particle every: '+str('{:.6}'.format(IMBHapp.value_in(units.yr))+' years'))
    """
    print('One timestep:       '+str('{:.4}'.format(eta*tend.value_in(units.yr))+' years'))

    while time < tend:
        iter += 1
        app_time = 0 | units.s
        if iter %100 == 0:
            print('Iteration: ', iter)

        merger_mass = 0 | units.MSun
        added_mass = 0 | units.MSun
        ejected_mass = 0 | units.MSun
        tcoll = 0 | units.yr

        time += eta*tend
        channel_IMBH["to_gravity"].copy()
        code.evolve_model(time)

        gc_code.d_update(parti[1].x, parti[1].y, parti[1].z)

        rtide = tidal_radius(parti)
        for particle in SMBH_filter(parti):
            rel_vel = ((particle.vx-parti[1].vx)**2+(particle.vy-parti[1].vy)**2+(particle.vz-parti[1].vz)**2).sqrt()
            parti_KE = 0.5*particle.mass*(rel_vel)**2
            temp_PE = [ ]
            temp_PE = indiv_PE_closest(particle, parti, temp_PE)
            parti_BE = max(temp_PE)
            dist_core = particle.position - parti[1].position
            dist_SMBH = particle.position - parti[0].position

            if parti_KE > parti_BE and dist_core.length() < rtide:
                print('KE > BE but inside tidal')
            if parti_KE < parti_BE and dist_core.length() > rtide:
                print('KE < BE but outside tidal')

            if dist_SMBH.length() < 0.3*gc_code.gc_dist:
                ejected = True
                ejected_key_track = particle.key_tracker
                ejected_mass = particle.mass
                print('...Leaving Simulation - Particle bound to SMBH...')
                print('Ejected time:         ', time)
                print('Ejected iter:         ', iter)
                extra_note = 'Stopped due to particle bound with SMBH'
                break
            
            if parti_KE > parti_BE and dist_core.length() > rtide:
                ejected = True
                ejected_key_track = particle.key_tracker
                ejected_mass = particle.mass
                print('...Leaving Simulation - Ejection Occured...')
                print('Ejected time:         ', time)
                print('Ejected iter:         ', iter)
                extra_note = 'Stopped due to ejection'
                break

            else:
                pass

        mergebin = 0
        injbin = 0
        if ejected == False:
            if stopping_condition.is_set():
                print("........Encounter Detected........")
                print('Collision at step: ', iter)
                print('Simulation will now stop')
                energy_before = code.potential_energy + code.kinetic_energy
                for ci in range(len(stopping_condition.particles(0))):
                    mergebin = 1
                    Nenc += 1
                    app_time = time
                    app_time2 = time

                    enc_particles = Particles(particles=[stopping_condition.particles(0)[ci],
                                                         stopping_condition.particles(1)[ci]])
                    enc_particles = enc_particles.get_intersecting_subset_in(parti)
                    ejected_key_track = parti[-1].key_tracker
                    merged_parti  = merge_IMBH(parti, enc_particles, code.model_time)
                    merger_mass = merged_parti.mass.sum()
                    cum_merger_mass += merger_mass
                    df_coll_tracker = pd.Series({'Collision Time': tcoll.in_(units.kyr), 
                                                 'Collided Particles': [enc_particles[0].key_tracker, enc_particles[1].key_tracker], 
                                                 'Initial Mass': [enc_particles[0].mass, enc_particles[1].mass] | units.MSun, 
                                                 'Emergent Particle': ejected_key_track, 'Collision Mass': merger_mass | units.MSun})                    
                    coll_tracker = coll_tracker.append(df_coll_tracker, ignore_index = True)   

                    

                    stab_timescale = time - stab_time_init
                    data_trackers.stable_sim_tracker(parti, injbin, mergebin, merger_mass, stab_timescale)
                    stab_time_init = time

                    parti.synchronize_to(code.particles)
                    energy_after = code.potential_energy + code.kinetic_energy
                    print('Energy change:     ', energy_after/energy_before - 1)
                    tcoll = time.in_(units.s) - eta*tend

                    df_eventstab_tracker = pd.Series({'Initial No. Particles': len(parti)-1, 'Merger Event': 1, 
                                                      'Collision Time': tcoll.in_(units.kyr), 'Injected Event': 0, 'Injected Time': 0 |units.s})
                    eventstab_tracker = eventstab_tracker.append(df_eventstab_tracker, ignore_index = True)

            """df_IMBH = pd.DataFrame()
            #if IMBH_adder.decision(time, decision_scale)==True:
            if time % IMBHapp < (eta*tend) and iter > 1:
                add_iter = iter
                print('.........New particle added.........')
                print('Added on step: ', add_iter)
                print('Added at time: ', str('{:.3f}'.format(time.value_in(units.kyr))), ' kyr')
                injbin = 1
                N_parti += 1
                Nenc = 0

                app_time = time
                app_time2 = time
                temp_E1 = code.kinetic_energy + code.potential_energy

                add_IMBH = IMBH_adder.add_IMBH(parti[1])
                added_mass = add_IMBH.mass
                parti.add_particle(add_IMBH)
                code.particles.add_particles(add_IMBH)
                temp_E2 = code.kinetic_energy + code.potential_energy

                tinj = time.in_(units.s) - eta*tend
                df_IMBH_vals = pd.Series({'key_tracker': add_IMBH.key_tracker[0]})
                df_IMBH = df_IMBH.append(df_IMBH_vals, ignore_index=True)

                stab_timescale = time - stab_time_init
                data_trackers.stable_sim_tracker(parti, injbin, mergebin, merger_mass, stab_timescale)
                stab_time_init = time

                E0 += (temp_E2-temp_E1) 

                df_eventstab_tracker = pd.Series({'Initial No. Particles': len(parti)-1, 'Merger Event': 0, 
                                                  'Collision Time': 0 | units.s, 'Injected Event': 1, 
                                                  'Injected Time': tinj.in_(units.kyr)})
                eventstab_tracker = eventstab_tracker.append(df_eventstab_tracker, ignore_index = True)

            IMBH_tracker = IMBH_tracker.append(df_IMBH, ignore_index=True)"""

        channel_IMBH["from_gravity"].copy()     
        rows = (N_parti)
        tdyn_val = tdyn_calc(parti) | units.yr

        df_IMBH = pd.DataFrame()
        for i in range(N_parti):
            for j in range(len(parti)):
                if IMBH_tracker.iloc[[i][0]][0] == parti[j].key_tracker:
                    parti_KE = 0.5*parti[j].mass*parti[j].velocity.length()**2
                    temp_PE = []
                    temp_PE = indiv_PE_all(parti[j], parti, temp_PE)
                    parti_PE = max(temp_PE)
                    df_IMBH_vals = pd.Series({'{}'.format(time): [parti[j].mass, parti[j].position, 
                                              parti_KE, parti_PE, tdyn_val[j]]})
                    break
                else:
                    df_IMBH_vals = pd.Series({'{}'.format(time): [np.NaN | units.MSun,
                                                                 [np.NaN | units.parsec, 
                                                                  np.NaN | units.parsec, 
                                                                  np.NaN | units.parsec],
                                                                  np.NaN | units.J,
                                                                  np.NaN | units.J,
                                                                  np.NaN | units.yr]})
            df_IMBH = df_IMBH.append(df_IMBH_vals, ignore_index=True)

        IMBH_tracker = IMBH_tracker.append(df_IMBH, ignore_index=True)
        IMBH_tracker['{}'.format(time)] = IMBH_tracker['{}'.format(time)].shift(-rows)
        IMBH_tracker = IMBH_tracker[pd.notnull(IMBH_tracker["key_tracker"])]

        df_com_tracker = pd.Series({'x': SMBH_filter(code.particles).center_of_mass()[0].in_(units.parsec),
                                    'y': SMBH_filter(code.particles).center_of_mass()[1].in_(units.parsec),
                                    'z': SMBH_filter(code.particles).center_of_mass()[2].in_(units.parsec)})
        com_tracker = com_tracker.append(df_com_tracker, ignore_index=True)

        df_LG_tracker = pd.Series({'Time': time.in_(units.kyr),
                                   'LG25': LagrangianRadii(code.particles[2:])[5].in_(units.parsec),
                                   'LG50': LagrangianRadii(code.particles[2:])[6].in_(units.parsec),
                                   'LG75': LagrangianRadii(code.particles[2:])[7].in_(units.parsec),
                                   'Tidal Radius': rtide.in_(units.parsec),
                                   'Relaxation Time': relax_timescale(cluster_radi, IMBH_adder.mass, gc_code.gc_pop).in_(units.yr)})
        LG_tracker = LG_tracker.append(df_LG_tracker , ignore_index=True)

        parti_KE = code.kinetic_energy
        parti_BE = code.potential_energy
        Etp  = parti_KE + parti_BE
        Et = Etp
        de = abs(Et-E0)/abs(E0)
        if 20 < iter:
            dEs = abs(Et-energy_tracker.iloc[19][1])/abs(energy_tracker.iloc[19][1])
            df_energy_tracker = pd.Series({'Time': time.in_(units.kyr), 'Et': Et, 'dE': de, 'dEs': dEs, 
                                           'Appearance': app_time, 'Collision Time': tcoll, 
                                           'Collision Mass': merger_mass})
        else:
            df_energy_tracker = pd.Series({'Time': time.in_(units.kyr), 'Et': Et, 'dE': de, 'dEs': 0, 
                                           'Appearance': app_time, 'Collision Time': tcoll, 
                                           'Collision Mass': merger_mass })
        energy_tracker = energy_tracker.append(df_energy_tracker, ignore_index=True)

        df_parti_energy_tracker = pd.Series({'Iteration': 0, "BE": parti_BE.in_(units.J), "KE": parti_KE.in_(units.J), "Total E": Etp})
        parti_energy_tracker = parti_energy_tracker.append(df_parti_energy_tracker, ignore_index=True)

        time1 = time
        if (ejected) or Nenc > 0:
            time = tend
    
    code.stop()
    comp_end = cpu_time.time()
    comp_time = comp_end-comp_start

    if iter > 5 + add_iter and iter > 10:
        no_plot = False
        #chaos_stab_timescale = time1 - app_time2
        chaos_stab_timescale = time1
        print("Total Merging Events: ", Nenc)
        print('Total integration time: ', comp_end-comp_start)
        print('...Dumping Files...')
        count = file_counter()

        if time1 == tend:
            ejected_key_track = parti[1].key_tracker
       
       # com_tracker.to_pickle('data/center_of_mass/IMBH_com_parsecs_'+str(N_parti_init)+str(count)+'_equal_mass.pkl')
       # IMBH_tracker.to_pickle('data/positions_IMBH/IMBH_positions_'+str(N_parti_init)+str(count)+'_equal_mass.pkl')
       # energy_tracker.to_pickle('data/energy/IMBH_energy_'+str(N_parti_init)+str(count)+'_equal_mass.pkl')
       # parti_energy_tracker.to_pickle('data/particle_energies/particle_energies_'+str(N_parti_init)+str(count)+'_equal_mass.pkl')
       # LG_tracker.to_pickle('data/lagrangians/IMBH_Lagrangian_'+str(N_parti_init)+str(count)+'_equal_mass.pkl')
       # coll_tracker.to_pickle('data/collision_events/IMBH_merge_events'+str(N_parti_init)+str(count)+'_equal_mass.pkl')
       # eventstab_tracker.to_pickle('data/event_tracker/IMBH_events'+str(N_parti_init)+str(count)+'_equal_mass.pkl')
        data_trackers.chaotic_sim_tracker(parti, initial_set, Nenc, cum_merger_mass, time1, ejected_key_track, 
                                          chaos_stab_timescale, added_mass, ejected_mass, comp_time)

        lines = ['Simulation: ', "Total CPU Time: "+str(comp_end-comp_start)+' seconds', 
                 'Timestep: '+str(eta),
                 'Simulated until: '+str(time1.value_in(units.yr))+str(' years'), 
                 'Cluster Radius: '+str(gc_code.gc_rad.value_in(units.parsec))+' parsecs', 
                 'Cluster Distance: '+str(cluster_distance.value_in(units.parsec))+' parsecs', 
                 'Masses of IMBH: '+str(parti.mass.value_in(units.MSun))+' MSun',
                 "No. of initial IMBH: "+str(init_IMBH_pop-2), 
                 'Number of new particles: '+str(N_parti-N_parti_init),
                 'Total Number of (Final) IMBH: '+str(len(parti)-2), 
                # 'IMBH Appearance Rate: '+str(IMBHapp.value_in(units.yr))+' years',
                 'Number of mergers: '+str(Nenc), 'End Time: '+str(tend.value_in(units.yr))+' years', 
                 'Integrator: Hermite (NO PN)',
                 'Extra Notes: ', extra_note]

        with open('data/simulation_stats/simulation'+str(count)+'.txt', 'w') as f:
            for line in lines:
                f.write(line)
                f.write('\n')        
    
    else:
        no_plot = True
        print('...No stability timescale - simulation ended too quick...')

    return no_plot