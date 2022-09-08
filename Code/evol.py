from amuse.lab import *
from amuse.units import units
from amuse.community.mikkola.interface import Mikkola
from amuse.couple import bridge
from initialiser import *
from evol_func import *
from amuse.ext.galactic_potentials import MWpotentialBovy2015
from datetime import datetime
import numpy as np
import pickle as pkl

def evolve_system(parti, tend, eta, converter):
    """
    Bulk of the simulation. Uses Mikkola integrator to evolve the system while bridging it.
    Keeps track of the particles' position and stores it in a pickle file.
    
    Inputs:
    parti:   The particle set needed to simulate
    tend:    The end time of the simulation
    eta:     The step size
    output:  The evolved simulation
    """

    SMBH_code      = MW_SMBH()
    MWG_code       = MWpotentialBovy2015()
    GC_code        = GC_pot()
    GC_parti_track = GC_init()
    conv = converter
    
    GC_tracker = GC_parti_track.GC_tracer(parti)                        #Initial position/vel. of the GC
    conv_gc = nbody_system.nbody_to_si(GC_tracker.mass.sum(),
                                       GC_tracker.position.sum())
    gravity_code_gc = drift_without_gravity(conv_gc)
    gravity_code_gc.particles.add_particles(GC_tracker)                 #Initiates how GC will evolve in time
    gc_track_array = GC_parti_track.gal_path_init()
    channel_gc = gravity_code_gc.particles.new_channel_to(GC_tracker)

    code = Mikkola(conv, number_of_workers = 3)
    code.particles.add_particles(parti)

    channel_IMBH = {"from_gravity": 
                    code.particles.new_channel_to(parti,
                    attributes=["x", "y", "z", "vx", "vy", "vz", "mass"],
                    target_names=["x", "y", "z", "vx", "vy", "vz", "mass"]),
                    "to_gravity": 
                    parti.new_channel_to(code.particles,
                    attributes=["mass", "collision_radius"],
                    target_names=["mass", "radius"])}

    brd = bridge.Bridge(timestep=1e-4 | units.yr)
    brd.add_system(gravity_code_gc, (SMBH_code, MWG_code))
    brd.add_system(code, (SMBH_code, MWG_code, GC_code))
   
    time = 0 | units.yr
    iter = 0
    Nenc = 0

    stopping_condition = code.stopping_conditions.collision_detection
    stopping_condition.enable()

    key_identify = [parti[i].key for i in range(len(parti))]
    IMBH_tracker = np.zeros((len(parti)+1, 20002, 3))
    for i in range(len(parti)):
        IMBH_tracker[i][0,:] = parti[i].position.in_(units.parsec).number
    IMBH_tracker[-1][0,:] = GC_tracker.position.in_(units.parsec).number

    com_tracker = np.zeros((1, 20002, 3))
    com_tracker[0][0] = parti.center_of_mass().in_(units.parsec).number

    E0 = brd.kinetic_energy + code.get_radiated_gravitational_energy()
    E0 += (parti[:].mass * (SMBH_code.get_potential_at_point(0, parti[:].x, parti[:].y, parti[:].z)
                             +  MWG_code.get_potential_at_point(0, parti[:].x, parti[:].y, parti[:].z)
                             +  GC_code.get_potential_at_point(0, parti[:].x, parti[:].y, parti[:].z))).sum()
    E0 += GC_code.gc_mass * (SMBH_code.get_potential_at_point(0, GC_tracker.x, GC_tracker.y, GC_tracker.z)
                             + MWG_code.get_potential_at_point(0, GC_tracker.x, GC_tracker.y, GC_tracker.z))

    arr_Et        = [ ]
    arr_de_stab   = [ ]
    arr_de        = [ ]
    arr_time      = [ ]
    parti.collision_radius = parti.radius * 3000
    encounter_particles = []
    
    adaptive_time = False

    while time < tend:
        iter += 1
        if iter%100 == 0:
            print('Iteration: ', iter)

        if (adaptive_time):
            adaptive_eta = adaptive_dt(eta, tend, parti).in_(units.yr)
            time += adaptive_eta

        else:
            time += eta*tend

        GC_code.d_update(GC_tracker.x, GC_tracker.y, GC_tracker.z)

        channel_IMBH["to_gravity"].copy()
        brd.evolve_model(time)
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
                
        channel_gc.copy()
        channel_IMBH["from_gravity"].copy()  

        gc_track_array[0].append(GC_tracker.x)
        gc_track_array[1].append(GC_tracker.y)
        gc_track_array[2].append(GC_tracker.z)
        
        for i in range(len(parti)+Nenc):
            for j in range(len(parti)):
                if parti[j].key_tracker == key_identify[i]:
                    IMBH_tracker[i][iter,:] = parti[j].position.in_(units.parsec).number# 
        IMBH_tracker[-1][iter,:] = GC_tracker.position.in_(units.parsec).number
        com_tracker[0][iter] = parti.center_of_mass().in_(units.parsec).number
        
        Et = brd.kinetic_energy + code.get_radiated_gravitational_energy()
        Et += (parti[:].mass * (SMBH_code.get_potential_at_point(0, parti[:].x, parti[:].y, parti[:].z)
                                + MWG_code.get_potential_at_point(0, parti[:].x, parti[:].y, parti[:].z)
                                + GC_code.get_potential_at_point(0, parti[:].x, parti[:].y, parti[:].z))).sum()
        Et += GC_code.gc_mass * (SMBH_code.get_potential_at_point(0, GC_tracker.x, GC_tracker.y, GC_tracker.z)
                                + MWG_code.get_potential_at_point(0, GC_tracker.x, GC_tracker.y, GC_tracker.z))
        de = abs(Et-E0)/abs(E0)
        
        arr_Et.append(Et)
        arr_de.append(de)
        arr_time.append(iter)
        
        if 15 < iter:
            arr_de_stab.append(abs(Et-arr_Et[13])/abs(arr_Et[13]))

    IMBH_tracker = IMBH_tracker | units.parsec
    com_tracker  = com_tracker  | units.parsec
    energy_tracker = [[arr_time], [arr_de], [arr_de_stab]]

    brd.stop()

    print('...Dumping Files...')

    with open('data/center_of_mass/IMBH_com_parsecs_'+str(datetime.now())+'.pkl', 'wb') as file:
        pkl.dump(com_tracker, file)

    with open('data/positions/IMBH_positions_parsecs_'+str(datetime.now())+'.pkl', 'wb') as file:
        pkl.dump(IMBH_tracker, file)

    with open('data/energy/IMBH_energy_'+str(datetime.now())+'.pkl', 'wb') as file:
        pkl.dump(energy_tracker, file)
