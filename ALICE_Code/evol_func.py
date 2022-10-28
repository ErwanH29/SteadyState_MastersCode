from amuse.lab import *
from amuse.units import units
from amuse.ext.orbital_elements import orbital_elements_from_binary
from parti_initialiser import *
from evol_func import *
import numpy as np

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

def merge_IMBH(parti, parti_in_enc, tcoll, int_string, code):
    """
    Function which merges two particles if the collision stopping condition has been met
    
    Inputs:
    parti:          The complete particle set being simulated
    parti_in_enc:   The particles in the collision
    tcoll:          The time-stamp for which the particles collide at
    int_string:     String telling whether simulating Hermite or GRX (for the PN cross-terms)
    code:           The integrator used
    """

    com_pos = parti_in_enc.center_of_mass()
    com_vel = parti_in_enc.center_of_mass_velocity()

    new_particle  = Particles(1)
    if calc_momentum(parti_in_enc[0]) > calc_momentum(parti_in_enc[1]):
        new_particle.key_tracker = parti_in_enc[0].key
    else: 
        new_particle.key_tracker = parti_in_enc[1].key
    
    new_particle.mass = parti_in_enc.total_mass()
    new_particle.collision_time = tcoll
    new_particle.position = com_pos
    new_particle.velocity = com_vel
    new_particle.radius = (2*constants.G*new_particle.mass)/(constants.c**2)
    new_particle.collision_radius = 3 * new_particle.radius
    parti.add_particles(new_particle)
    parti.remove_particles(parti_in_enc)

    if int_string != 'Hermite':
        code.small_particles.remove_particle(parti_in_enc[1])
        if new_particle.mass > 1e6 | units.MSun:
            code.large_particles.remove_particles(parti_in_enc[0])
            code.large_particles.add_particles(new_particle)
        else:
            code.small_particles.remove_particles(parti_in_enc[0])
            code.small_particles.add_particles(new_particle)

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
            rel_pos = indivp.position - pset[i].position
            min_dist.append(rel_pos.length().value_in(units.parsec))
    temp = np.sort(min_dist)
    index = np.where(min_dist == temp[0])[0]
    index2 = np.where(min_dist == temp[1])[0]

    return min(min_dist), pset[index], pset[index2]

def SMBH_filter(pset):
    return pset[pset.mass < 5*10**4 | units.MSun]

def tidal_radius(pset):
    """
    Function to outline the tidal radius. Uses equation 5.8 and 5.10 of Spitzer 1969.
    
    Inputs:
    pset:  The complete particle set
    """

    SMBH = MW_SMBH()
    gc_code = globular_cluster()
    com_particle = Particles(1)
    com_particle.mass = SMBH_filter(pset).total_mass()
    com_particle.radius = 0 | units.RSun
    com_particle.velocity = SMBH_filter(pset).center_of_mass_velocity()
    com_particle.position = SMBH_filter(pset).center_of_mass() - pset[0].position

    m1, m2, semimajor, ecc, inc, argp, ascn, tanom = bin_global(pset[0], com_particle)
    perigal = semimajor*(1-ecc)
    xe = ((3+ecc)**-1 * (com_particle.mass)/SMBH.mass * (perigal)**3)**(1/3)
#    xe = ((pset[1].mass)/SMBH.mass * gc_code.gc_dist**3)**(1/3)

    return xe

