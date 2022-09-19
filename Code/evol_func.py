from amuse.lab import *
from amuse.units import units
from initialiser import *
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

def file_counter():
    """
    Function which counts the number of files in a directory.
    """

    dir_path = r'data/center_of_mass/'
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

def dynamical_fric(pos_distr, vel_distr, mass_distr, particle, eta, tend):
    """
    Function to compute the dynamical friction as per Petts et al. 2012 eqn 9.
    
    Inputs:
    pos_distr:  A theoretical Plummer model positional distribution sampled over 1e7 stars
    vel_distr:  A theoretical Plummer model velocity distribution
    mass_distr: A theoretical Plummer model mass distribution
    particle:   The individual particle who is having its dynamical friction computed
    eta:        Simulation time-step
    tend:       Simulation final time
    output:     Value of the dynamical friction
    """

    systv = (4/3)*np.pi*(10**-3 | units.parsec)**3       # The volume for which the BHs are in (line 44 interface.py)
    systm = 10**7 | units.MSun                           # The mass of the GC which the BHs are in (line 44 interface.py)

    index = find_nearest(pos_distr, particle.position.length().value_in(units.AU))
    enc_mass = mass_distr[index] | units.MSun
    index2 = find_nearest(vel_distr, particle.velocity.length().value_in(units.kms))
    frac_vel = index2/len(vel_distr)

    value = -2*np.pi*(constants.G)**2*particle.mass*(systm/systv)*np.log((particle.mass/enc_mass)**2+1) \
            * frac_vel*(particle.velocity.length())**-3 * eta * tend * particle.velocity
    return value
    
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

def reset_grav(particle, integrator, conv):
        code = integrator(conv, number_of_workers = 3)
        code.particles.add_particles(particle)
        code.commit_particles()
        epsilon = 1e-16
        Lw = 4 * abs(math.log10(epsilon)) + 32
        code.set_PN_terms(1,1,1,0)
        code.set_bs_tolerance(1e-16)

        channel_IMBH = {"from_gravity": 
                        code.particles.new_channel_to(particle,
                        attributes=["x", "y", "z", "vx", "vy", "vz", "mass"],
                        target_names=["x", "y", "z", "vx", "vy", "vz", "mass"]),
                        "to_gravity": 
                        particle.new_channel_to(code.particles,
                        attributes=["mass", "collision_radius"],
                        target_names=["mass", "radius"])} 

        return code, channel_IMBH