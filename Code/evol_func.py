from amuse.lab import *
from amuse.units import units
from parti_initialiser import *
from evol_func import *
import numpy as np
import math
import fnmatch
import os
import matplotlib.pyplot as plt

class drift_without_gravity(object):
    """
    Class which induces a kick-drift-kick mechanism so IMBH feel 
    the background potential
    """

    def __init__(self, convert_nbody, time= 0 |units.Myr):
        self.model_time = time
        self.convert_nbody = convert_nbody
        self.particles = Particles()

    def evolve_model(self, t_end):
        """
        Evolves the GC in time gravitationally. 
        Note the particles here are NOT IMBH, but sample the GC
        defined particle (see initialiser.py)
        """
        
        dt = t_end - self.model_time
        self.particles.position += self.particles.velocity*dt
        self.model_time = t_end

    @property
    def potential_energy(self):
        return quantities.zero

    @property 
    def kinetic_energy(self):
        return (0.5*self.particles.mass \
                *self.particles.velocity.lengths()**2).sum()

    def stop(self):
        pass

def calc_momentum(parti):
    """
    Function which calculates the momentum of the particles in collision

    Inputs:
    parti:    The colliding particles
    output:   The momentum value for each particle
    """

    mom_value = (parti.mass * parti.velocity).sum()
    return mom_value

def df_timescale(particle, clust_rad, halfmass):
    """
    Dynamical friction timescale to compute the time needed before IMBH enters cluster.
    Uses eqn. 8.12 of Tremaine and Binney (2007).
    
    Inputs:
    particle:   Single particle of the initialised set to extract approx. value from.
    halfmass:   The half-mass radius of the complete particle set (IMBH)
    dist_const: The radius constant used in TB07
    vel_const:  The vel. constant used in TB07
    mass_const: The mass constant used in TB07
    sigmav:     Assumed to be ~ (circular_velocity)/(sqrt(2)) as per TB07
    clust_rad:  The initial radius of the cluster particles orbit at
    """

    dist_const = 5 | units.kpc
    vel_const  = 200 | units.kms
    mass_const = 10**8 | units.MSun
    sigmav = (np.sqrt(2))**-1 * particle.velocity.length()
    rh = halfmass
    rtyp = (constants.G*particle.mass)/(particle.velocity.length()**2)
    Lambda = (particle.position.length()/max(rh, rtyp))
    return 19/(np.log(Lambda))*(clust_rad/dist_const)**2 * (sigmav / vel_const) * (mass_const/particle.mass) * 1 | units.Gyr

def file_counter():
    """
    Function which counts the number of files in a directory.
    """

    dir_path = r'data/simulation_stats/'
    count = len(fnmatch.filter(os.listdir(dir_path), '*.*'))
    return count

def find_nearest(array, value):
    """
    Function to find the nearest value in an array for a particular element.
    outputs: The index where the nearest array-value is found.
    """

    array = np.asarray(array)
    index = (np.abs(array - value)).argmin()
    return index

def merge_IMBH(parti, particles_in_encounter, tcoll):
    """
    Function which merges two particles if the collision stopping condition has been met
    
    Inputs:
    parti:                    The complete particle set being simulated
    particles_in_encounter:   The particles in the collision
    tcoll:                    The time-stamp for which the particles collide at
    outputs:                  Removal of two colliding particles, while addition of the merging product
    """

    com_pos = particles_in_encounter.center_of_mass()
    com_vel = particles_in_encounter.center_of_mass_velocity()

    new_particle  = Particles(1)

    if calc_momentum(particles_in_encounter[0]) > calc_momentum(particles_in_encounter[1]):
        new_particle.key_tracker = particles_in_encounter[0].key
    else: 
        new_particle.key_tracker = particles_in_encounter[1].key

    new_particle.mass = particles_in_encounter.total_mass()
    new_particle.collision_time = tcoll
    new_particle.position = com_pos
    new_particle.velocity = com_vel
    new_particle.radius = (2*constants.G*new_particle.mass)/(constants.c**2)
    new_particle.collision_radius = new_particle.radius * 10
    
    parti.add_particles(new_particle)
    parti.remove_particles(particles_in_encounter)
    return new_particle