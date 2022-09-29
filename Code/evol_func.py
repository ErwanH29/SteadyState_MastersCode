from amuse.lab import *
from amuse.units import units
from amuse.ext.orbital_elements import orbital_elements_from_binary
from parti_initialiser import *
from evol_func import *
import numpy as np
import fnmatch
import os

def bin_global(parti1, parti2):
    """
    Function which computes the Kepler elements for a specific binary.
    
    Inputs:
    parti1: The first particle in the binary
    parti2: The second particle in the binary
    """

    bin_sys =  Particles()
    bin_sys.add_particle(parti1)
    bin_sys.add_particle(parti2)

    kepler_elements = orbital_elements_from_binary(bin_sys, G=constants.G)    
    mass1 = kepler_elements[0]
    mass2 = kepler_elements[1]
    semimajor = kepler_elements[2]
    eccentric = kepler_elements[3]
    inclinate = kepler_elements[4]
    arg_peri  = kepler_elements[5]
    asc_node  = kepler_elements[6]
    true_anom = kepler_elements[7]

    return  mass1, mass2, semimajor, eccentric, inclinate, arg_peri, asc_node, true_anom

def calc_momentum(indivp):
    """
    Function which calculates the momentum of the particles in collision

    Inputs:
    indivp: The colliding particles
    """

    mom_value = (indivp.mass * indivp.velocity).sum()
    return mom_value

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

    Inputs:
    array:  The array for which has all the elements
    value:  The value to compare the elements with
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
    if particles_in_encounter.total_mass() > 10**6 | units.MSun:
        new_particle.key_tracker = SMBH_filter(particles_in_encounter).key_tracker
    else:
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

def nearest_neighbour(indivp, pset):
    """
    Function to find the nearest particle to some individual.
    
    Inputs:
    indivp: The individual particle
    pset:   The complete particle set
    """

    min_dist = [ ]
    for i in range(len(pset)):
        if indivp == pset[i]:
            pass
        else:
            vec_x = indivp.x - pset[i].x
            vec_y = indivp.y - pset[i].y
            vec_z = indivp.z - pset[i].z
            dist = (vec_x**2+vec_y**2+vec_z**2).sqrt()
            min_dist.append(dist)
    return min(min_dist)

def SMBH_filter(pset):
    return pset[pset.mass < 10**6 | units.MSun]


def tidal_radius1(parti):
    SMBH = MW_SMBH()
    return (2/3)*(1/3 * SMBH_filter(parti).mass.sum()/SMBH.mass * (SMBH_filter(parti).center_of_mass()).length()**3)**(1/3)

def tidal_radius(pset):
    """
    Function to outline the tidal radius. Uses equation 5.8 and 5.10 of Spitzer 1969.
    
    Inputs:
    pset:  The complete particle set
    """

    SMBH = MW_SMBH()

    new_parti = Particle()
    new_parti.mass = SMBH_filter(pset).mass.sum()
    new_parti.position = SMBH_filter(pset).center_of_mass()
    new_parti.velocity = SMBH_filter(pset).center_of_mass_velocity()

    m1, m2, semimajor, ecc, inc, argp, ascn, tanom = bin_global(pset[0], new_parti)
    perigal = semimajor*(1-ecc)
    xe = ((3+ecc)**-1 * (new_parti.mass)/SMBH.mass * (perigal)**3)**(1/3)

    return (2/3)*(xe)
