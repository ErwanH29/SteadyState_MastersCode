from amuse.lab import *
from amuse.units import units
from parti_initialiser import *
from data_init import *
from physics_func import *
from evol_func import *
from amuse.ext.LagrangianRadii import LagrangianRadii
from amuse.community.hermite import Hermite
#from amuse.community.hermite_grx.interface import *
import numpy as np
import pandas as pd
import time as cpu_time

# FOR ALICE SIMULATIONS:
#   - REMOVE ALL PRINT STATEMENTS + EXCESS FILE (ONLY KEEP CHAOTIC + STABLE + SIMULATION FILES)
#   - CHANGE number_of_workers
#   - ENSURE CORRECT FILE SAVING FORMAT (NAME)

def evolve_system(parti, tend, eta, cluster_distance, cluster_radi, converter, int_string):
    """
    Bulk of the simulation. Uses Hermite integrator to evolve the system while bridging it.
    Keeps track of the particles' position and stores it in a pickle file.
    
    Inputs:
    parti:            The particle set needed to simulate
    tend:             The end time of the simulation
    eta:              The step size
    converter:        Variable used to convert between nbody units and SI
    """
    
    set_printing_strategy("custom", preferred_units = [units.MSun, units.pc, units.kyr, units.AU/units.yr, units.J],
                                                       precision = 20, prefix = "", separator = " ", suffix = "")

    comp_start = cpu_time.time()

    SMBH_code = MW_SMBH()
    gc_code = globular_cluster()
    IMBH_adder = IMBH_init()

    initial_set = parti.copy()

    converter = nbody_system.nbody_to_si(parti.mass.sum(), gc_code.gc_dist)
    code = Hermite(converter, number_of_workers = 6) #To change when going into ALICE
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
    app_time2 = time #To activate when allowing ejections to occur in a continuous simulation
    cum_merger_mass = 0 | units.MSun
    stab_timescale = time

    init_IMBH_pop = len(parti)
    add_iter = 0
    N_parti = init_IMBH_pop
    init_pop = N_parti-2
    ejected = False
    extra_note = ''

    parti_KE = code.kinetic_energy
    parti_BE = code.potential_energy
    E0p = parti_KE + parti_BE 
    E0 = E0p

    data_trackers = data_initialiser()
    IMBH_tracker = data_trackers.IMBH_tracker(parti, time, init_IMBH_pop)
    LG_tracker = data_trackers.LG_tracker(cluster_radi, IMBH_adder.mass, gc_code.gc_pop, parti, time, code)
    #energy_tracker = data_trackers.energy_tracker(E0p, time, 0 | units.s)
    #parti_energy_tracker = data_trackers.parti_energy_tracker(parti_BE, parti_KE, E0p, time)
    coll_tracker = data_trackers.coll_tracker()
    eventstab_tracker = data_trackers.event_tracker(parti)

    while time < tend:
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
            rel_vx = particle.vx#-SMBH_filter(parti).center_of_mass_velocity()[0]
            rel_vy = particle.vy#-SMBH_filter(parti).center_of_mass_velocity()[1]
            rel_vz = particle.vz#-SMBH_filter(parti).center_of_mass_velocity()[2]
            rel_vvect = [rel_vx, rel_vy, rel_vz]
            rel_vel = (rel_vx**2+rel_vy**2+rel_vz**2).sqrt()

            dist_core = particle.position - parti[0].position

            dist_vect = np.sqrt(np.dot(dist_core, dist_core))
            vel_vect  = np.sqrt(np.dot(rel_vvect, rel_vvect))
            curr_traj = abs(np.dot(dist_core, rel_vvect))/(dist_vect * vel_vect)

            parti_KE = 0.5*particle.mass*(rel_vel)**2
            temp_PE = [ ]
            temp_PE = indiv_PE_closest(particle, parti, temp_PE)
            parti_BE = np.sum(temp_PE)
            SMBH_potential = particle.mass*SMBH_code.get_potential_at_point(0, particle.x, particle.y, particle.z)
            tot_PE = abs(SMBH_potential) + abs(parti_BE)

            dist_cluster_core = particle.position - SMBH_filter(parti).center_of_mass()

            if parti_KE > tot_PE and dist_cluster_core.length() > 4*gc_code.gc_rad:
                print('High enough energy and far enough, but towards cluster', curr_traj)
            if dist_cluster_core.length() > 4*gc_code.gc_rad and curr_traj > 0:
                print('not enough E', parti_KE, tot_PE)
            if parti_KE > tot_PE and curr_traj > 0:
                print('Not far enough')

            if parti_KE > tot_PE and dist_cluster_core.length() > 4*gc_code.gc_rad and curr_traj > 0:
                ejected = True
                ejected_key_track = particle.key_tracker
                ejected_mass = particle.mass
                extra_note = 'Stopped due to particle bound with SMBH'
                break

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

                    enc_particles_set = Particles(particles=[stopping_condition.particles(0)[ci],
                                                         stopping_condition.particles(1)[ci]])
                    enc_particles = enc_particles_set.get_intersecting_subset_in(parti)
                    merged_parti = merge_IMBH(parti, enc_particles, code.model_time)
                    merger_mass = merged_parti.mass.sum()
                    cum_merger_mass += merger_mass

                    stab_timescale = time - stab_timescale   #To record the time of the last ejection/merging event
                    data_trackers.stable_sim_tracker(parti, injbin, mergebin, merger_mass, stab_timescale) 

                    parti.synchronize_to(code.particles)
                    energy_after = code.potential_energy + code.kinetic_energy
                    tcoll = time.in_(units.s) - eta*tend

                    df_eventstab_tracker = pd.Series({'Initial No. Particles': len(parti)-1, 'Merger Event': 1, 
                                                      'Collision Time': tcoll.in_(units.kyr), 'Injected Event': 0, 'Injected Time': 0 |units.s})
                    eventstab_tracker = eventstab_tracker.append(df_eventstab_tracker, ignore_index = True)

                    ejected_key_track = parti[-1].key_tracker
                    df_coll_tracker = pd.Series({'Collision Time': tcoll.in_(units.kyr), 
                                                 'Collided Particles': [enc_particles_set[0].key, enc_particles_set[1].key], 
                                                 'Initial Mass': [enc_particles_set[0].mass, enc_particles_set[1].mass] | units.MSun, 
                                                 'Emergent Particle': ejected_key_track, 'Collision Mass': merger_mass | units.MSun})                    
                    coll_tracker = coll_tracker.append(df_coll_tracker, ignore_index = True) 

            mergebin = 0

        channel_IMBH["from_gravity"].copy()     
        rows = (N_parti)
        tdyn_val = tdyn_calc(parti) | units.yr
        df_IMBH = pd.DataFrame()
        for i in range(N_parti):
            for j in range(len(parti)):
                if IMBH_tracker.iloc[i][0][0] == parti[j].key_tracker:
                    parti_KE = 0.5*parti[j].mass*parti[j].velocity.length()**2
                    temp_PE = []
                    temp_PE = indiv_PE_all(parti[j], parti, temp_PE)
                    parti_PE = max(temp_PE)
                    gcframe_vel = parti[j].velocity - parti[1].velocity
                    df_IMBH_vals = pd.Series({'{}'.format(time): [parti[j].key_tracker, parti[j].mass, parti[j].position,
                                              gcframe_vel, parti_KE, parti_PE, tdyn_val[j]]})
                    break
                else:
                    df_IMBH_vals = pd.Series({'{}'.format(time): [np.NaN, np.NaN | units.MSun,
                                                                 [np.NaN | units.parsec, np.NaN | units.parsec, np.NaN | units.parsec],
                                                                 [np.NaN | units.kms, np.NaN | units.kms, np.NaN | units.kms],
                                                                  np.NaN | units.J,
                                                                  np.NaN | units.J,
                                                                  np.NaN | units.yr]})
            df_IMBH = df_IMBH.append(df_IMBH_vals, ignore_index=True)
            
        IMBH_tracker = IMBH_tracker.append(df_IMBH, ignore_index=True)
        IMBH_tracker['{}'.format(time)] = IMBH_tracker['{}'.format(time)].shift(-rows)
        IMBH_tracker = IMBH_tracker.dropna(axis='rows')

        df_LG_tracker = pd.Series({'Time': time.in_(units.kyr),
                                   'LG25': LagrangianRadii(code.particles[2:])[5].in_(units.parsec),
                                   'LG50': LagrangianRadii(code.particles[2:])[6].in_(units.parsec),
                                   'LG75': LagrangianRadii(code.particles[2:])[7].in_(units.parsec),
                                   'Tidal Radius': 4*gc_code.gc_rad.in_(units.parsec),
                                   'Relaxation Time': relax_timescale(cluster_radi, IMBH_adder.mass, gc_code.gc_pop).in_(units.yr)})
        LG_tracker = LG_tracker.append(df_LG_tracker , ignore_index=True)

        """parti_KE = code.kinetic_energy
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
        parti_energy_tracker = parti_energy_tracker.append(df_parti_energy_tracker, ignore_index=True)"""

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
        count = file_counter(int_string)

        if time1 == tend:
            ejected_key_track = parti[1].key_tracker
       
        path = '/home/s2009269/data1/HERMITE_NGC_data/'
        IMBH_tracker.to_pickle(os.path.join(path+str('particle_trajectory'), 'IMBH_'+str(int_string)+'_Local_'+str(init_pop)+'_sim'+str(count)+'_init_dist'+str('{:.3f}'.format(gc_code.gc_dist.value_in(units.parsec)))
                               +'_equal_mass_'+str('{:.3f}'.format(parti[2].mass.value_in(units.MSun)))+'.pkl'))
       # energy_tracker.to_pickle('data/'+str(int_string)+'/energy/IMBH_'+str(int_string)+'energy_'+str(init_pop)+str(count)+'_equal_mass.pkl')
       # parti_energy_tracker.to_pickle('data/'+str(int_string)+'/particle_energies/particle_energies_'+str(int_string)+str(init_pop)+str(count)+'_equal_mass.pkl')
        LG_tracker.to_pickle(os.path.join(path+str('lagrangians'), 'IMBH_'+str(int_string)+'Lagrangian_'+str(init_pop)+str(count)+'_equal_mass.pkl'))
        coll_tracker.to_pickle(os.path.join(path+str('collision_events'), 'IMBH_'+str(int_string)+'merge_events'+str(init_pop)+str(count)+'_equal_mass.pkl'))
        eventstab_tracker.to_pickle(os.path.join(path+str('event_tracker'), 'IMBH_'+str(int_string)+'events'+str(init_pop)+str(count)+'_equal_mass.pkl'))   #For different dependent variables [clust_dist, clust_rad, ALICE/local...], change output file name"""
        data_trackers.chaotic_sim_tracker(parti, initial_set, Nenc, cum_merger_mass, time1, ejected_key_track, 
                                          chaos_stab_timescale, added_mass, ejected_mass, comp_time, int_string)                   #For different dependent variables [clust_dist, clust_rad, ALICE/local...], change output file name

        lines = ['Simulation: ', "Total CPU Time: "+str(comp_end-comp_start)+' seconds', 
                 'Timestep: '+str(eta),
                 'Simulated until: '+str(time1.value_in(units.yr))+str(' years'), 
                 'Cluster Radius: '+str(gc_code.gc_rad.value_in(units.parsec))+' parsecs', 
                 'Cluster Distance: '+str(cluster_distance.value_in(units.parsec))+' parsecs', 
                 'Masses of IMBH: '+str(parti.mass.value_in(units.MSun))+' MSun',
                 "No. of initial IMBH: "+str(init_IMBH_pop-2), 
                 'Number of new particles: '+str(N_parti-init_pop),
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