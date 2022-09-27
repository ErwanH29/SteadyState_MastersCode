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
    clust_rad:  The initial radius of the cluster particles orbit at
    halfmass:   The half-mass radius of the complete particle set (IMBH)
    dist_const: The radius constant used in TB07
    vel_const:  The vel. constant used in TB07
    mass_const: The mass constant used in TB07
    sigmav:     Assumed to be ~ (circular_velocity)/(sqrt(2)) as per TB07
    """

    SMBH_parti = MW_SMBH()

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

def nearest_neighbour(indiv, set):
    """
    Function to find the nearest particle to some individual.
    
    Inputs:
    indiv:   The individual particle
    set:     The complete particle set
    """

    min_dist = [ ]
    for i in range(len(set)):
        if indiv == set[i]:
            pass
        else:
            vec_x = indiv.x - set[i].x
            vec_y = indiv.y - set[i].y
            vec_z = indiv.z - set[i].z
            dist = (vec_x**2+vec_y**2+vec_z**2).sqrt()
            min_dist.append(dist)
    return min(min_dist)


def SMBH_filter(parti):
    """
    Function to remove the SMBH from the particle set
    
    Inputs:
    parti:  The particle set
    """
    return parti[parti.mass < 10**6 | units.MSun]


def tidal_radius1(parti):
    SMBH = MW_SMBH()
    return (2/3)*(1/3 * SMBH_filter(parti).mass.sum()/SMBH.mass * (SMBH_filter(parti).center_of_mass()).length()**3)**(1/3)

def tidal_radius(parti):
    """
    Function to outline the tidal radius. Uses equation 5.8 and 5.10 of Spitzer 1969.
    
    Inputs:
    parti:  The complete particle set
    """

    SMBH = MW_SMBH()

    new_parti = Particle()
    new_parti.mass = SMBH_filter(parti).mass.sum()
    new_parti.position = SMBH_filter(parti).center_of_mass()
    new_parti.velocity = SMBH_filter(parti).center_of_mass_velocity()

    m1, m2, semimajor, ecc, inc, argp, ascn, tanom = bin_global(parti[0], new_parti)
    perigal = semimajor*(1-ecc)
    xe = ((3+ecc)**-1 * (new_parti.mass)/SMBH.mass * (perigal)**3)**(1/3)

    return (2/3)*((3+ecc)**-1 * (new_parti.mass)/SMBH.mass * (perigal)**3)**(1/3)
