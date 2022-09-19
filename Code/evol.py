from cmath import nan
from amuse.lab import *
from amuse.units import units
from initialiser import *
from physics_func import *
from evol_func import *
from amuse.ext.LagrangianRadii import LagrangianRadii
import numpy as np
import pandas as pd
import math
import time as cpu_time
import glob
import os
import pickle as pkl

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
    set_printing_strategy("custom", preferred_units = [units.MSun, units.AU, units.yr, units.AU/units.yr],
                                    precision = 4, prefix = "", separator = "[", suffix = "]")

    comp_start = cpu_time.time()

    """
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
    """

    code = grav_solver(converter, number_of_workers = 3)
    code.particles.add_particles(parti)
    code.commit_particles()
    epsilon = 1e-16
    Lw = 4 * abs(math.log10(epsilon)) + 32
    code.set_PN_terms(1,1,1,1)
    code.set_bs_tolerance(1e-16)
    code.calculate_word_length()
    
    """
    brd = bridge.Bridge(timestep=1e-4 | units.yr)
    brd.add_system(gravity_code_gc, (SMBH_code, MWG_code))
    brd.add_system(code, (SMBH_code, MWG_code, GC_code))
    """
    
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
    filename = glob.glob('data/preliminary_calcs/*')
    with open(os.path.join(max(filename, key=os.path.getctime)), 'rb') as input_file:
        temp_data = pkl.load(input_file)
    pos_values = temp_data.iloc[0][0]
    mas_values = temp_data.iloc[0][1]
    vel_values = temp_data.iloc[0][2]
    for particle in parti:
        temp_vel = dynamical_fric(pos_values, vel_values, mas_values, particle, 10**-3, tend)
        particle.velocity = particle.velocity + temp_vel

    IMBH_array = pd.DataFrame()
    df_IMBH    = pd.DataFrame()
    for i in range(init_IMBH_pop):
        parti_KE = 0.5*parti[i].velocity.length()**2
        temp_PE = []
        temp_PE = indiv_PE(parti[i], parti, temp_PE)
        parti_PE = max(temp_PE)
        df_IMBH_vals = pd.Series({'key_tracker': parti[i].key_tracker, '{}'.format(time): [parti[i].mass, parti[i].position, parti_KE, parti_PE]})
        df_IMBH      = df_IMBH.append(df_IMBH_vals, ignore_index=True)
    IMBH_array = IMBH_array.append(df_IMBH, ignore_index=True)

    """
    GC_array = pd.DataFrame()
    df_GC_tracker = pd.Series({'x': GC_tracker.position.x.in_(units.parsec),
                               'y': GC_tracker.position.y.in_(units.parsec),
                               'z': GC_tracker.position.z.in_(units.parsec)})
    GC_array = GC_array.append(df_GC_tracker, ignore_index=True)
    """

    com_tracker = pd.DataFrame()
    df_com_tracker = pd.Series({'x': parti.center_of_mass()[0].in_(units.parsec),
                                'y': parti.center_of_mass()[1].in_(units.parsec),
                                'z': parti.center_of_mass()[2].in_(units.parsec)})
    com_tracker = com_tracker.append(df_com_tracker, ignore_index=True)

    LG_array = pd.DataFrame()
    df_LG_tracker = pd.Series({'Time': time.in_(units.kyr),
                               'LG25': LagrangianRadii(code.particles)[5].in_(units.parsec),
                               'LG75': LagrangianRadii(code.particles)[7].in_(units.parsec)})
    LG_array = LG_array.append(df_LG_tracker , ignore_index=True)

    energy_tracker = pd.DataFrame()

    """
    E0 = brd.kinetic_energy #+ brd.potential_energy
    E0 += (parti[:].mass * (SMBH_code.get_potential_at_point(0, parti[:].x, parti[:].y, parti[:].z)
                             +  MWG_code.get_potential_at_point(0 | units.kpc, parti[:].x, parti[:].y, parti[:].z)
                             +  GC_code.get_potential_at_point(0, parti[:].x, parti[:].y, parti[:].z))).sum()
    E0 += GC_code.gc_mass * (SMBH_code.get_potential_at_point(0, GC_tracker.x, GC_tracker.y, GC_tracker.z)
                             + MWG_code.get_potential_at_point(0 | units.kpc, GC_tracker.x, GC_tracker.y, GC_tracker.z))
    """

    parti_KE = code.kinetic_energy
    parti_BE = code.potential_energy
    E0p = parti_KE + parti_BE
    E0 = E0p

    df_energy_tracker = pd.Series({'Time': time.in_(units.kyr), 'Et': E0 , 'dE': 0, 'dEs': 0, 'Collision Time': 0 | units.s, 'Collision Mass': 0 | units.MSun })
    energy_tracker = energy_tracker.append(df_energy_tracker, ignore_index=True)
    parti_energy_tracker = pd.DataFrame()
    df_parti_energy_tracker = pd.Series({'Time': time.in_(units.kyr), "BE": parti_BE.in_(units.J), 
                                         "KE": parti_KE.in_(units.J), "Total E": E0p})
    parti_energy_tracker = parti_energy_tracker.append(df_parti_energy_tracker, ignore_index=True)

    binary_array = pd.DataFrame()

    tdyn_tracker = pd.DataFrame()
    df_tdyn = pd.DataFrame()
    tdyn_val = tdyn_calc(parti) | units.s
    for i in range(init_IMBH_pop):
        df_tdyn_vals = pd.Series({'key_tracker': parti[i].key_tracker, '{}'.format(time): [tdyn_val[i]]})
        df_tdyn = df_tdyn.append(df_tdyn_vals, ignore_index = True)
    tdyn_tracker = tdyn_tracker.append(df_tdyn, ignore_index = True)
    com_tracker = com_tracker.append(df_com_tracker, ignore_index=True)

    while time < tend:
        iter += 1
        rows = (Nenc+len(parti))
        print('Iteration:    ', iter)
        time += eta*tend

        """
        if new_particle:
            key_identifier.append(new_particle.key)

        GC_code.d_update(GC_tracker.x, GC_tracker.y, GC_tracker.z)
        """

        channel_IMBH["to_gravity"].copy()
        #brd.evolve_model(time)
        code.evolve_model(time)

        for particle in parti:
            temp_vel = dynamical_fric(pos_values, vel_values, mas_values, particle, 10**-3, tend)
            particle.velocity = particle.velocity + temp_vel

        """
        #Rudimentary binary detector
        binaries = parti.get_binaries()
        df_bin   = pd.DataFrame()
        if (len(binaries)>0):
            print("...Binaries Detected...")
            for i in range(2):
                df_bin_vals = pd.Series({'key_tracker': binaries[0][i].key_tracker, '{}'.format(time): [binaries[0][i].position]})
            df_bin = df_bin.append(df_bin_vals, ignore_index=True)
        binary_array = binary_array.append(df_bin, ignore_index = True)
        """

        merger_mass = 0 | units.MSun
        enc = False
        tcoll = 0 | units.yr

        #Rudimentary collision detector
        for i in range(len(parti)):
            for j in range(len(parti)):
                if i == j or i > j:
                    pass
                else:
                    distance = abs(parti[i].position.length() - parti[j].position.length())
                    if (parti[i].collision_radius + parti[j].collision_radius) > distance:
                        energy1 = code.kinetic_energy + code.potential_energy
                        enc = True
                        Nenc += 1
                        print("...Encounter Detected...")
                        print("Merging Events: ", Nenc)
                        print('Particles in encounter:  ', enc_particles)
                        enc_particles = Particles(particles=[parti[i], parti[j]])
                        merged_parti  = merge_IMBH(parti, enc_particles, code.model_time)
                        print('Merged Particle: ', merged_parti)
                        print('Merger mass:    ', merged_parti.mass.sum())
                        tcoll = time.in_(units.yr)
                        break
                        
            if (enc):
                code.stop()
                code, channel_IMBH = reset_grav(parti, grav_solver, converter)
                energy2 = code.kinetic_energy + code.potential_energy
                print('Coll. DeltaE: ', abs(energy1 - energy2)/abs(energy1))
                print('----------------------------------------------------------- \\'
                      '-----------------------------------------------------------')
                break

        parti.move_to_center()
        channel_IMBH["from_gravity"].copy()      
        """
        channel_gc.copy()
        gc_track_array[0].append(GC_tracker.x)
        gc_track_array[1].append(GC_tracker.y)
        gc_track_array[2].append(GC_tracker.z)
        """
        
        df_IMBH = pd.DataFrame()
        df_tdyn = pd.DataFrame()
        tdyn_val = tdyn_calc(parti) | units.yr
        for i in range(len(parti)+Nenc):
            for j in range(len(parti)):
                if IMBH_array.iloc[[i][0]][0] == parti[j].key_tracker:
                    parti_KE = 0.5*parti[j].velocity.length()**2
                    temp_PE = []
                    temp_PE = indiv_PE(parti[j], parti, temp_PE)
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

        """
        df_GC_tracker = pd.Series({'x': GC_tracker.position.x.in_(units.parsec),
                                'y': GC_tracker.position.y.in_(units.parsec),
                                'z': GC_tracker.position.z.in_(units.parsec)})
        GC_array = GC_array.append(df_GC_tracker, ignore_index=True)
        """

        df_com_tracker = pd.Series({'x': parti.center_of_mass()[0].in_(units.parsec),
                                    'y': parti.center_of_mass()[1].in_(units.parsec),
                                    'z': parti.center_of_mass()[2].in_(units.parsec)})
        com_tracker = com_tracker.append(df_com_tracker, ignore_index=True)

        df_LG_tracker = pd.Series({'Time': time.in_(units.kyr),
                                   'LG25': LagrangianRadii(code.particles)[5].in_(units.parsec),
                                   'LG75': LagrangianRadii(code.particles)[7].in_(units.parsec)})
        LG_array = LG_array.append(df_LG_tracker , ignore_index=True)
        """
        Et = brd.kinetic_energy #+ brd.potential_energy
        Et += (parti[:].mass * (SMBH_code.get_potential_at_point(0, parti[:].x, parti[:].y, parti[:].z)
                                + MWG_code.get_potential_at_point(0 | units.kpc, parti[:].x, parti[:].y, parti[:].z)
                                + GC_code.get_potential_at_point(0, parti[:].x, parti[:].y, parti[:].z))).sum()
        Et += GC_code.gc_mass * (SMBH_code.get_potential_at_point(0, GC_tracker.x, GC_tracker.y, GC_tracker.z)
                                + MWG_code.get_potential_at_point(0 | units.kpc, GC_tracker.x, GC_tracker.y, GC_tracker.z))
        """

        parti_KE = code.kinetic_energy
        parti_BE = code.potential_energy
        Etp = parti_KE + parti_BE
        Et = Etp

        de = abs(Et-E0)/abs(E0)
        print(de)
        if 16 < iter:
            dEs = abs(Et-energy_tracker.iloc[15][1])/abs(energy_tracker.iloc[15][1])
            df_energy_tracker = pd.Series({'Time': time.in_(units.kyr), 'Et': Et, 'dE': de, 'dEs': dEs, 
                                           'Collision Time': tcoll, 'Collision Mass': merger_mass})
        else:
            df_energy_tracker = pd.Series({'Time': time.in_(units.kyr), 'Et': Et, 'dE': de, 'dEs': 0, 
                                           'Collision Time': tcoll, 'Collision Mass': merger_mass })
        energy_tracker = energy_tracker.append(df_energy_tracker, ignore_index=True)

        df_parti_energy_tracker = pd.Series({'Iteration': 0, "BE": parti_BE.in_(units.J), "KE": parti_KE.in_(units.J), "Total E": Etp})
        parti_energy_tracker = parti_energy_tracker.append(df_parti_energy_tracker, ignore_index=True)

    code.stop()
    comp_end = cpu_time.time()

    print("Total Merging Events: ", Nenc)
    print('Total integration time: ', comp_end-comp_start)
    print('...Dumping Files...')
    count = file_counter()

    com_tracker.to_pickle('data/center_of_mass/IMBH_com_parsecs_'+str(count)+'.pkl')
    IMBH_array.to_pickle('data/positions_IMBH/IMBH_positions_'+str(count)+'.pkl')
    """
    GC_array.to_pickle('data/positions_GC/GC_positions_'+str(count)+'.pkl')
    """
    energy_tracker.to_pickle('data/energy/IMBH_energy_'+str(count)+'.pkl')
    parti_energy_tracker.to_pickle('data/particle_energies/particle_energies_'+str(count)+'.pkl')
    LG_array.to_pickle('data/lagrangians/IMBH_Lagrangian_'+str(count)+'.pkl')
    tdyn_tracker.to_pickle('data/dynamical_time/IMBH_Dynamical_Time'+str(count)+'.pkl')
    
    lines = ['Simulation: ', "Total CPU Time: "+str(comp_end-comp_start), 'Timestep: '+str(eta),
             'End Time: '+str(tend.value_in(units.yr))+' years', 'Integrator: '+str(grav_solver), 
             'Epsilon: '+str(epsilon), 'PN Terms: ' + str([1,1,1,1]),
             "No. of initial IMBH: "+str(init_IMBH_pop), 'Number of mergers: '+str(Nenc)]

    with open('data/simulation_stats/simulation'+str(count)+'.txt', 'w') as f:
        for line in lines:
            f.write(line)
            f.write('\n')