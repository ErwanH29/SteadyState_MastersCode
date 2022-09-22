from cmath import nan
from amuse.lab import *
from amuse.units import units
from parti_initialiser import *
from physics_func import *
from evol_func import *
from amuse.ext.LagrangianRadii import LagrangianRadii
from amuse.community.hermite import Hermite
import numpy as np
import pandas as pd
import time as cpu_time

def evolve_system(parti, tend, eta, cluster_distance, cluster_radi, distr_string, converter):
    """
    Bulk of the simulation. Uses Mikkola integrator to evolve the system while bridging it.
    Keeps track of the particles' position and stores it in a pickle file.
    
    Inputs:
    parti:          The particle set needed to simulate
    tend:           The end time of the simulation
    eta:            The step size
    cluster_distance: Initialised cluster radius.
    converter:      Variable used to convert between nbody units and SI
    output:         The evolved simulation
    """

    set_printing_strategy("custom", preferred_units = [units.MSun, units.AU, units.yr, units.AU/units.yr, units.J],
                                    precision = 20, prefix = "", separator = " ", suffix = "")

    comp_start = cpu_time.time()

    code = Hermite(converter, number_of_workers = 3)
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
    init_IMBH_pop = len(parti)

    IMBH_array = pd.DataFrame()
    df_IMBH    = pd.DataFrame()

    for i in range(init_IMBH_pop):
        parti_KE = 0.5*parti[i].velocity.length()**2
        temp_PE = []
        temp_PE = indiv_PE_BH(parti[i], parti, temp_PE)
        parti_PE = max(temp_PE)
        df_IMBH_vals = pd.Series({'key_tracker': parti[i].key_tracker, '{}'.format(time): [parti[i].mass, parti[i].position, parti_KE, parti_PE]})
        df_IMBH      = df_IMBH.append(df_IMBH_vals, ignore_index=True)
    IMBH_array = IMBH_array.append(df_IMBH, ignore_index=True)

    com_tracker = pd.DataFrame()
    df_com_tracker = pd.Series({'x': SMBH_filter(code.particles).center_of_mass()[0].in_(units.parsec),
                                'y': SMBH_filter(code.particles).center_of_mass()[1].in_(units.parsec),
                                'z': SMBH_filter(code.particles).center_of_mass()[2].in_(units.parsec)})
    com_tracker = com_tracker.append(df_com_tracker, ignore_index=True)

    LG_array = pd.DataFrame()
    df_LG_tracker = pd.Series({'Time': time.in_(units.kyr),
                               'LG25': LagrangianRadii(code.particles[1:])[5].in_(units.parsec),
                               'LG50': LagrangianRadii(code.particles[1:])[6].in_(units.parsec),
                               'LG75': LagrangianRadii(code.particles[1:])[7].in_(units.parsec)})
    LG_array = LG_array.append(df_LG_tracker , ignore_index=True)

    energy_tracker = pd.DataFrame()

    SMBH = MW_SMBH()
    parti_KE = code.kinetic_energy
    parti_BE = code.potential_energy
    E0p  = parti_KE + parti_BE 
    E0 = E0p

    app_time = np.NaN | units.s
    df_energy_tracker = pd.Series({'Time': time.in_(units.kyr), 'Et': E0 , 'dE': 0, 'dEs': 0, 'Appearance': app_time, 
                                   'Collision Time': 0 | units.s, 'Collision Mass': 0 | units.MSun })
    energy_tracker = energy_tracker.append(df_energy_tracker, ignore_index=True)
    parti_energy_tracker = pd.DataFrame()
    df_parti_energy_tracker = pd.Series({'Time': time.in_(units.kyr), "BE": parti_BE.in_(units.J), 
                                         "KE": parti_KE.in_(units.J), "Total E": E0p})
    parti_energy_tracker = parti_energy_tracker.append(df_parti_energy_tracker, ignore_index=True)

    tdyn_tracker = pd.DataFrame()
    df_tdyn = pd.DataFrame()
    tdyn_val = tdyn_calc(parti) | units.s

    for i in range(init_IMBH_pop):
        if i == 0:
            df_tdyn_vals = pd.Series({'key_tracker': parti[i].key_tracker, '{}'.format(time): [0 | units.s]})
            df_tdyn = df_tdyn.append(df_tdyn_vals, ignore_index = True)
        else:
            df_tdyn_vals = pd.Series({'key_tracker': parti[i].key_tracker, '{}'.format(time): [tdyn_val[i]]})
            df_tdyn = df_tdyn.append(df_tdyn_vals, ignore_index = True)
    tdyn_tracker = tdyn_tracker.append(df_tdyn, ignore_index = True)
    com_tracker = com_tracker.append(df_com_tracker, ignore_index=True)

    IMBHapp = df_timescale(parti[1], cluster_distance, LagrangianRadii(code.particles[1:])[6].in_(units.parsec))
    decision_scale = (eta*tend)/IMBHapp #
    IMBH_adder = IMBH_init()
    N_parti = init_IMBH_pop
    N_parti_init = N_parti
    extra_note = ''

    ejected = False
    new_parti = False

    while time < tend:
        iter += 1
        if iter %100 == 0:
            print('Iteration: ', iter)
        time += eta*tend
        channel_IMBH["to_gravity"].copy()
        code.evolve_model(time)

        avgx_cluster = SMBH_filter(parti).x.mean().value_in(units.AU)
        avgy_cluster = SMBH_filter(parti).y.mean().value_in(units.AU)
        avgz_cluster = SMBH_filter(parti).z.mean().value_in(units.AU)
        avg_clustpos = [avgx_cluster, avgy_cluster, avgz_cluster] * 1 | units.AU


        for particle in SMBH_filter(parti):
            parti_KE  = 0.5*particle.mass*(particle.velocity.length())**2
            temp_PE   = [ ]
            temp_PE   = indiv_PE(particle, parti, temp_PE)
            parti_BE  = max(temp_PE)
            dist_core = particle.position - avg_clustpos

            if parti_KE > parti_BE and dist_core.length() > 10*cluster_radi:
                ejected = True
                ejected_key_track = particle.key_tracker
                print('...Leaving Simulation - Ejection Occured...')
                print('Final set:        ', parti)
                print('Ejected particle: ', ejected_key_track)
                print('Ejected time:     ', time)
                print('Ejected iter:     ', iter)
                extra_note = 'Stopped due to ejection'
                break
            else:
                pass

        if ejected == False:
            merger_mass = 0 | units.MSun
            tcoll = 0 | units.yr
            if stopping_condition.is_set():
                print("........Encounter Detected........")
                print('Collision at step: ', iter)
                energy_before = code.potential_energy + code.kinetic_energy
                print('Energy before:     ', energy_before)
                for ci in range(len(stopping_condition.particles(0))):
                    Nenc += 1
                    enc_particles = Particles(particles=[stopping_condition.particles(0)[ci],
                                                        stopping_condition.particles(1)[ci]])
                    enc_particles = enc_particles.get_intersecting_subset_in(parti)
                    print('Particles in enc:  ', enc_particles)
                    merged_parti  = merge_IMBH(parti, enc_particles, code.model_time)
                    parti.synchronize_to(code.particles)
                    merger_mass = merged_parti.mass.sum()
                    energy_after = code.potential_energy + code.kinetic_energy
                    print('Energy after:      ', energy_after)
                    print('Energy change:     ', energy_after/energy_before - 1)
                    tcoll = time.in_(units.s) - eta*tend

        parti.move_to_center()
        channel_IMBH["from_gravity"].copy()      

        if ejected == False:
            df_IMBH = pd.DataFrame()
            df_tdyn = pd.DataFrame()
            app_time = np.NaN | units.s
            if IMBH_adder.decision(time, decision_scale)==True:
                print('.........New particle added.........')
                print('Added on step: ', iter)
                N_parti += 1
                app_time = time
                temp_E1 = code.kinetic_energy + code.potential_energy
                temp_pos = SMBH_filter(code.particles).center_of_mass()
                add_IMBH = IMBH_adder.add_IMBH(temp_pos, distr_string, converter)
                add_IMBH.velocity += SMBH_filter(code.particles).center_of_mass_velocity()
                parti.add_particle(add_IMBH)
                code.particles.add_particles(add_IMBH)
                temp_E2 = code.kinetic_energy + code.potential_energy

                df_IMBH_vals = pd.Series({'key_tracker': add_IMBH.key_tracker[0]})
                df_IMBH      = df_IMBH.append(df_IMBH_vals, ignore_index=True)
                df_tdyn_vals = pd.Series({'key_tracker': add_IMBH.key_tracker})
                df_tdyn = df_tdyn.append(df_tdyn_vals, ignore_index = True)
                E0 += (temp_E2-temp_E1)  

            tdyn_tracker = tdyn_tracker.append(df_tdyn, ignore_index = True)
            IMBH_array = IMBH_array.append(df_IMBH, ignore_index=True)
            rows = (N_parti)

        df_IMBH = pd.DataFrame()
        df_tdyn = pd.DataFrame()
        tdyn_val = tdyn_calc(parti) | units.yr
        for i in range(N_parti):
            for j in range(len(parti)):
                if IMBH_array.iloc[[i][0]][0] == parti[j].key_tracker:
                    parti_KE = 0.5*parti[j].velocity.length()**2
                    temp_PE = []
                    temp_PE = indiv_PE_BH(parti[j], parti, temp_PE)
                    parti_PE = max(temp_PE)
                    df_IMBH_vals = pd.Series({'{}'.format(time): [parti[j].mass, parti[j].position, parti_KE, parti_PE]})
                    df_tdyn_vals = pd.Series({'key_tracker': parti[j].key_tracker, '{}'.format(time): [tdyn_val[j]]})
                    break
                else:
                    df_IMBH_vals = pd.Series({'{}'.format(time): [np.NaN | units.MSun,
                                                                 [np.NaN | units.parsec, 
                                                                  np.NaN | units.parsec, 
                                                                  np.NaN | units.parsec]]})
                    df_tdyn_vals = pd.Series({'key_tracker': nan, '{}'.format(time): [np.NaN | units.s]})
            df_IMBH = df_IMBH.append(df_IMBH_vals, ignore_index=True)
            df_tdyn = df_tdyn.append(df_tdyn_vals, ignore_index = True)

        IMBH_array = IMBH_array.append(df_IMBH, ignore_index=True)
        IMBH_array['{}'.format(time)] = IMBH_array['{}'.format(time)].shift(-rows)
        IMBH_array = IMBH_array[pd.notnull(IMBH_array["key_tracker"])]

        tdyn_tracker = tdyn_tracker.append(df_tdyn, ignore_index = True)
        tdyn_tracker['{}'.format(time)] = tdyn_tracker['{}'.format(time)].shift(-rows)
        tdyn_tracker = tdyn_tracker[pd.notnull(tdyn_tracker['{}'.format(time)])]

        df_com_tracker = pd.Series({'x': SMBH_filter(code.particles).center_of_mass()[0].in_(units.parsec),
                                    'y': SMBH_filter(code.particles).center_of_mass()[1].in_(units.parsec),
                                    'z': SMBH_filter(code.particles).center_of_mass()[2].in_(units.parsec)})
        com_tracker = com_tracker.append(df_com_tracker, ignore_index=True)

        df_LG_tracker = pd.Series({'Time': time.in_(units.kyr),
                                   'LG25': LagrangianRadii(SMBH_filter(code.particles))[5].in_(units.parsec),
                                   'LG50': LagrangianRadii(SMBH_filter(code.particles))[6].in_(units.parsec),
                                   'LG75': LagrangianRadii(SMBH_filter(code.particles))[7].in_(units.parsec)})
        LG_array = LG_array.append(df_LG_tracker , ignore_index=True)

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
        if (ejected):
            time = tend
        if (new_parti):
            time = tend
        
    
    code.stop()
    comp_end = cpu_time.time()

    if len(parti) < 4:
        print('Ejected at two particles - ignored results')
        no_plot = True
        return no_plot
    else:
        stab_tracker = pd.DataFrame()
        df_stabtime = pd.Series({'Initial Particles': init_IMBH_pop, 'Final Particles': (len(parti)-1), 
                                'Number of Mergers': Nenc, 'Simulated Till': time.in_(units.yr),
                                'Ejected Particle': ejected_key_track})
        stab_tracker = stab_tracker.append(df_stabtime, ignore_index = True)

        print("Total Merging Events: ", Nenc)
        print('Total integration time: ', comp_end-comp_start)
        print('...Dumping Files...')
        count = file_counter()

        stab_tracker.to_pickle('data/stability_time/IMBH_Hermite_'+str(count)+'.pkl')
        com_tracker.to_pickle('data/center_of_mass/IMBH_com_parsecs_'+str(count)+'.pkl')
        IMBH_array.to_pickle('data/positions_IMBH/IMBH_positions_'+str(count)+'.pkl')
        energy_tracker.to_pickle('data/energy/IMBH_energy_'+str(count)+'.pkl')
        parti_energy_tracker.to_pickle('data/particle_energies/particle_energies_'+str(count)+'.pkl')
        LG_array.to_pickle('data/lagrangians/IMBH_Lagrangian_'+str(count)+'.pkl')
        tdyn_tracker.to_pickle('data/dynamical_time/IMBH_Dynamical_Time'+str(count)+'.pkl')
        
        lines = ['Simulation: ', "Total CPU Time: "+str(comp_end-comp_start)+' seconds', 'Timestep: '+str(eta),
                'Simulated until: '+str(time1.value_in(units.yr))+str('year'), 
                'Cluster Radius: '+str(cluster_distance.value_in(units.parsec))+' parsecs', 
                'Masses of IMBH: '+str(parti.mass.value_in(units.MSun))+' MSun',
                "No. of initial IMBH: "+str(init_IMBH_pop-1), 'Number of new particles: '+str(N_parti-N_parti_init),
                'Total Number of (Final) IMBH: '+str(len(parti)), 'IMBH Appearance Rate: '+str(IMBHapp.value_in(units.yr))+' years',
                'Number of mergers: '+str(Nenc), 'End Time: '+str(tend.value_in(units.yr))+' years', 'Integrator: Hermite (NO PN)',
                'Extra Notes: ', extra_note]

        with open('data/simulation_stats/simulation'+str(count)+'.txt', 'w') as f:
            for line in lines:
                f.write(line)
                f.write('\n')