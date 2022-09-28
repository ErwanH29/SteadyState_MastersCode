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

def merge_IMBH(parti, enc_part, tcoll):
    """
    Function which merges two particles if the collision stopping condition has been met
    
    Inputs:
    parti:     The complete particle set being simulated
    enc_part:  The particles in the collision
    tcoll:     The time-stamp for which the particles collide at
    """

    com_pos = enc_part.center_of_mass()
    com_vel = enc_part.center_of_mass_velocity()

    new_particle  = Particles(1)
    if enc_part.total_mass() > 10**6 | units.MSun:
        new_particle.key_tracker = SMBH_filter(enc_part).key_tracker
    else:
        if calc_momentum(enc_part[0]) > calc_momentum(enc_part[1]):
            new_particle.key_tracker = enc_part[0].key
        else: 
            new_particle.key_tracker = enc_part[1].key

    new_particle.mass = enc_part.total_mass()
    new_particle.collision_time = tcoll
    new_particle.position = com_pos
    new_particle.velocity = com_vel
    new_particle.radius = (2*constants.G*new_particle.mass)/(constants.c**2)
    new_particle.collision_radius = new_particle.radius * 10
    new_particle.coll_events = 1 + enc_part.coll_events.sum()
    
    parti.add_particles(new_particle)
    parti.remove_particles(enc_part)
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
