from amuse.lab import *
from amuse.units import units
from amuse.community.mikkola.interface import Mikkola
from amuse.community.ph4.interface import ph4
from amuse.community.brutuspn.interface import Brutus
from amuse.couple import bridge
from initialiser import *
from evol_func import *
from amuse.ext.galactic_potentials import MWpotentialBovy2015
from datetime import datetime
import numpy as np
import pickle as pkl
import pandas as pd
import math


def evolve_system(parti, tend, eta, grav_solver, converter):
    """
    Bulk of the simulation. Uses Mikkola integrator to evolve the system while bridging it.
    Keeps track of the particles' position and stores it in a pickle file.
    
    Inputs:
    parti:   The particle set needed to simulate
    tend:    The end time of the simulation
    eta:     The step size
    output:  The evolved simulation
    """
    #MWG            = MWG_parameters()
    #print(MWG_code.parameters)
    #MWG_code.kinetic_energy = quantities.zero
    #MWG_code.potential_energy = quantities.zero
    #MWG_code.get_potential_at_point

    SMBH_code      = MW_SMBH()
    MWG_code       = MWpotentialBovy2015()
    GC_code        = GC_pot()
    GC_parti_track = GC_init()

    GC_tracker = GC_parti_track.GC_tracer(parti)                        #Initial position/vel. of the GC
    conv_gc = nbody_system.nbody_to_si(GC_tracker.mass.sum(),
                                    GC_tracker.position.sum())

    gravity_code_gc = drift_without_gravity(conv_gc)
    gravity_code_gc.particles.add_particles(GC_tracker)                 #Initiates how GC will evolve in time
    gc_track_array = GC_parti_track.gal_path_init()
    channel_gc = gravity_code_gc.particles.new_channel_to(GC_tracker)
    parti.collision_radius = parti.radius * 3000

    code = grav_solver(converter, number_of_workers = 3)
    code.particles.add_particles(parti)
    code.commit_particles()
    epsilon = 1e-16
    Lw = 4 * abs(math.log10(epsilon)) + 32
    code.set_PN_terms(1,0,0,0)
    code.set_bs_tolerance(1e-16)
    code.calculate_word_length()

    #brd = bridge.Bridge(timestep=1e-4 | units.yr)
    #brd.add_system(gravity_code_gc, (SMBH_code, MWG_code))
    #brd.add_system(code, (SMBH_code, MWG_code, GC_code))

    channel_IMBH = {"from_gravity": 
                code.particles.new_channel_to(parti,
                attributes=["x", "y", "z", "vx", "vy", "vz", "mass"],
                target_names=["x", "y", "z", "vx", "vy", "vz", "mass"]),
                "to_gravity": 
                parti.new_channel_to(code.particles,
                attributes=["mass", "collision_radius"],
                target_names=["mass", "radius"])} 

    #stopping_condition = code.stopping_conditions.collision_detection
    #stopping_condition.enable()

    adaptive_time = False
    time = 0 | units.yr
    iter = 0
    Nenc = 0

    IMBH_array = pd.DataFrame()
    df_IMBH    = pd.DataFrame()
    for i in range(len(parti)):
        df_IMBH_vals = pd.Series({'key_tracker': parti[i].key_tracker, '{}'.format(time): [parti[i].position]})
        df_IMBH      = df_IMBH.append(df_IMBH_vals, ignore_index=True)
        IMBH_array = IMBH_array.append(df_IMBH, ignore_index=True)


    GC_array = pd.DataFrame()
    df_GC_tracker = pd.Series({'x': GC_tracker.position.x.in_(units.parsec),
                            'y': GC_tracker.position.y.in_(units.parsec),
                            'z': GC_tracker.position.z.in_(units.parsec)})
    GC_array = GC_array.append(df_GC_tracker, ignore_index=True)

    com_tracker = pd.DataFrame()
    df_com_tracker = pd.Series({'x': parti.center_of_mass()[0].in_(units.parsec),
                            'y': parti.center_of_mass()[1].in_(units.parsec),
                            'z': parti.center_of_mass()[2].in_(units.parsec)})
    com_tracker = com_tracker.append(df_com_tracker, ignore_index=True)

    #E0 = brd.kinetic_energy
    E0 = (parti[:].mass * (SMBH_code.get_potential_at_point(0, parti[:].x, parti[:].y, parti[:].z)
                            +  MWG_code.get_potential_at_point(0 | units.kpc, parti[:].x, parti[:].y, parti[:].z)
                            +  GC_code.get_potential_at_point(0, parti[:].x, parti[:].y, parti[:].z))).sum()
    E0 += GC_code.gc_mass * (SMBH_code.get_potential_at_point(0, GC_tracker.x, GC_tracker.y, GC_tracker.z)
                            + MWG_code.get_potential_at_point(0 | units.kpc, GC_tracker.x, GC_tracker.y, GC_tracker.z))

    energy_tracker = pd.DataFrame()
    df_energy_tracker = pd.Series({'Iteration': 0, 'Et': E0 , 'dE': 0, 'dEs': 0 })
    energy_tracker = energy_tracker.append(df_energy_tracker, ignore_index=True)

    encounter_particles = []

    tend = 0.1 | units.yr
    eta = 10**-3

    while time < tend:
        iter += 1
        rows = (Nenc+len(parti))

        if iter%100 == 0:
            print('Iteration: ', iter)
        if (adaptive_time):
            adaptive_eta = adaptive_dt(eta, tend, parti).in_(units.yr)
            time += adaptive_eta
        else:
            time += eta*tend

        """
        if new_particle:
            key_identifier.append(new_particle.key)
        """

        GC_code.d_update(GC_tracker.x, GC_tracker.y, GC_tracker.z)

        channel_IMBH["to_gravity"].copy()
        code.evolve_model(time)

        for i in range(len(parti)):
            for j in range(len(parti)):
                if i == j or i > j:
                    pass
                else:
                    print(i,j)
                    distance = np.sqrt((parti[i].position.x-parti[j].position.x)**2
                                        +(parti[i].position.y-parti[j].position.y)**2
                                        +(parti[i].position.z-parti[j].position.z)**2)
                    if 0.5*(parti[i].collision_radius + parti[j].collision_radius) > distance:
                        print("...Encounter Detected")
                        Nenc += 1
                        for ci in range(2):
                            enc_particles = Particles(particles=[parti[i], parti[j]])
                            print('Complete particle set:   ', parti)
                            merged_parti  = merge_IMBH(parti, enc_particles, code.model_time)
                            parti.synchronize_to(code.particles)
                            print('Updated particle set:    ', parti)
                            print('Merger mass:  ', merged_parti.mass.sum())
        """
        if stopping_condition.is_set():
            print("........Encounter Detected........")
            print('Collision at step: ', iter)
            
            for ci in range(len(stopping_condition.particles(0))):
                Nenc += 1
                enc_particles = Particles(particles=[stopping_condition.particles(0)[ci],
                                            stopping_condition.particles(1)[ci]])
                enc_particles = enc_particles.get_intersecting_subset_in(parti)
                print('Complete particle set:   ', parti)

                merged_parti  = merge_IMBH(parti, enc_particles, code.model_time)
                parti.synchronize_to(code.particles)
                print('Updated particle set:    ', parti)
                print('Merger mass:  ', merged_parti.mass.sum())
                encounter_particles.append(merged_parti)
                adaptive_eta = 10**-3 | units.yr
        """
                
        channel_gc.copy()
        channel_IMBH["from_gravity"].copy()  

        gc_track_array[0].append(GC_tracker.x)
        gc_track_array[1].append(GC_tracker.y)
        gc_track_array[2].append(GC_tracker.z)

        df_IMBH = pd.DataFrame()
        for i in range(len(parti)+Nenc):
            for j in range(len(parti)):
                if IMBH_array.iloc[[i][0]][0] == parti[j].key_tracker:
                    df_IMBH_vals = pd.Series({'{}'.format(time): [parti[j].position]})
                    break
                else:
                    df_IMBH_vals = pd.Series({'{}'.format(time): [[np.NaN | units.parsec, 
                                                                    np.NaN | units.parsec, 
                                                                    np.NaN | units.parsec]]})
            
            df_IMBH = df_IMBH.append(df_IMBH_vals, ignore_index=True)

        IMBH_array = IMBH_array.append(df_IMBH, ignore_index=True)
        IMBH_array['{}'.format(time)] = IMBH_array['{}'.format(time)].shift(-rows)
        IMBH_array = IMBH_array[pd.notnull(IMBH_array["key_tracker"])]

        df_GC_tracker = pd.Series({'x': GC_tracker.position.x.in_(units.parsec),
                                'y': GC_tracker.position.y.in_(units.parsec),
                                'z': GC_tracker.position.z.in_(units.parsec)})
        GC_array = GC_array.append(df_GC_tracker, ignore_index=True)

        df_com_tracker = pd.Series({'x': parti.center_of_mass()[0].in_(units.parsec),
                                    'y': parti.center_of_mass()[1].in_(units.parsec),
                                    'z': parti.center_of_mass()[2].in_(units.parsec)})
        com_tracker = com_tracker.append(df_com_tracker, ignore_index=True)

        Et = brd.kinetic_energy 
        Et += (parti[:].mass * (SMBH_code.get_potential_at_point(0, parti[:].x, parti[:].y, parti[:].z)
                                + MWG_code.get_potential_at_point(0 | units.kpc, parti[:].x, parti[:].y, parti[:].z)
                                + GC_code.get_potential_at_point(0, parti[:].x, parti[:].y, parti[:].z))).sum()
        Et += GC_code.gc_mass * (SMBH_code.get_potential_at_point(0, GC_tracker.x, GC_tracker.y, GC_tracker.z)
                                + MWG_code.get_potential_at_point(0 | units.kpc, GC_tracker.x, GC_tracker.y, GC_tracker.z))
        de = abs(Et-E0)/abs(E0)
        print(de)
        if 16 < iter:
            dEs = abs(Et-energy_tracker.iloc[15][1])/abs(energy_tracker.iloc[15][1])
            df_energy_tracker = pd.Series({'Iteration': iter, 'Et': Et, 'dE': de, 'dEs': dEs})
        else:
            df_energy_tracker = pd.Series({'Iteration': iter, 'Et': Et, 'dE': de, 'dEs': 0 })
        energy_tracker = energy_tracker.append(df_energy_tracker, ignore_index=True)
    brd.stop()

    print('...Dumping Files...')
    com_tracker.to_pickle('data/center_of_mass/IMBH_com_parsecs_'+str(datetime.now())+'.pkl')
    IMBH_array.to_pickle('data/positions_IMBH/IMBH_positions_'+str(datetime.now())+'.pkl')
    GC_array.to_pickle('data/positions_GC/GC_positions_'+str(datetime.now())+'.pkl')
    energy_tracker.to_pickle('data/energy/IMBH_energy_'+str(datetime.now())+'.pkl')
    