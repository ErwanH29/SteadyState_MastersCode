from amuse.lab import *
from amuse.units import units
from parti_initialiser import *
from evol_func import *
import numpy as np


def df_timescale(indivp, clust_rad, halfmass):
    """
    Dynamical friction timescale to compute the time needed before IMBH enters cluster.
    Uses eqn. 8.12 of Tremaine and Binney (2007).
    
    Inputs:
    indivp:     Single particle of the initialised set to extract approx. value from.
    clust_rad:  The initial radius of the cluster particles orbit at
    halfmass:   The half-mass radius of the complete particle set (IMBH)
    dist_const: The radius constant used in TB07
    vel_const:  The vel. constant used in TB07
    mass_const: The mass constant used in TB07
    sigmav:     Assumed to be ~ (circular_velocity)/(sqrt(2)) as per TB07
    """

    dist_const = 5 | units.kpc
    vel_const  = 200 | units.kms
    mass_const = 10**8 | units.MSun
    sigmav = (np.sqrt(2))**-1 * indivp.velocity.length()
    rh = halfmass
    rtyp = (constants.G*indivp.mass)/(indivp.velocity.length()**2)
    Lambda = (indivp.position.length()/max(rh, rtyp))
    return 19/(np.log(Lambda))*(clust_rad/dist_const)**2 * (sigmav / vel_const) * (mass_const/indivp.mass) * 1 | units.Gyr

def df_velchange(indivp, vel_disp, no_stars):
    """
    Function which computes the speed at which particles sink inwards.
    Use of BT 2008 eqn 8.11.
    
    Inputs:
    indivp:   The individual particle who will have a change in velocity
    vel_disp: The velocity dispersion of the system
    no_stars: The number of stars in the cluster
    """

    avgx_cluster = SMBH_filter(indivp).x.mean().value_in(units.AU)
    avgy_cluster = SMBH_filter(indivp).y.mean().value_in(units.AU)
    avgz_cluster = SMBH_filter(indivp).z.mean().value_in(units.AU)
    avg_clustpos = [avgx_cluster, avgy_cluster, avgz_cluster] * 1 | units.AU

    return (-0.428*constants.G*indivp.mass)/(vel_disp*(indivp.position.length()-avg_clustpos.length()))*np.log(0.1*no_stars)

def indiv_max_PE(indivp, set):
    """
    Finding a particles' maximum PE

    Input:
    indivp:  The individual particle computing BE for
    set:     The complete particle set
    """

    SMBH = MW_SMBH()

    temp_PE_array = []
    for j in range(len(set)):
        if indivp == set[j]:
            pass
        else:
            distance = (indivp.position.length()-set[j].position.length())
            temp_PE  = ((constants.G*indivp.mass*set[j].mass)/abs(distance) \
                        + indivp.mass * SMBH.get_potential_at_point(0, indivp.x, indivp.y, indivp.z))
            temp_PE_array.append(temp_PE)

    return max(temp_PE_array)

def indiv_PE_BH(indivp, set, array):
    """
    Finding a particles' individual PE based on its closest binary

    Input:
    indivp:  The individual particle computing BE for
    set:     The complete particle set
    """

    SMBH = MW_SMBH()
    for comp_ in set:
        if indivp == comp_:
            pass
        else:
            distance = (indivp.position.length()-comp_.position.length())
            temp_PE  = abs(((constants.G*indivp.mass*comp_.mass)/abs(distance) \
                            + indivp.mass * SMBH.get_potential_at_point(0, indivp.x, indivp.y, indivp.z)))
            array.append(temp_PE)

    return array

def indiv_PE(indivp, set, array):
    """
    Finding a particles' individual PE based on its closest binary

    Input:
    indivp:  The individual particle computing BE for
    set:     The complete particle set
    """

    SMBH = MW_SMBH()

    for comp_ in set:
        if indivp == comp_ or comp_.mass > 10**6 | units.MSun:
            pass
        else:
            x_pos = indivp.x - comp_.position.x
            y_pos = indivp.y - comp_.position.y
            z_pos = indivp.z - comp_.position.z
            distance = (x_pos*x_pos + y_pos*y_pos + z_pos*z_pos).sqrt()
            
            temp_PE  = abs(((constants.G*indivp.mass*comp_.mass)/abs(distance)))
            array.append(temp_PE)

    return array

def relax_timescale(clust_rad, clust_mass, no_stars):
    """
    Function to record the relaxation timescale of the cluster.
    The equation is taken from Spitzer 1987 (eqn. 2.63) and is for
    isotehrmal systems.
    
    Inputs:
    clust_rad:  The cluster half-mass radius
    clust_mass: The cluster mass
    n_stars:    The number of stars present in the cluster
    """

    return (np.sqrt(clust_rad**3/(constants.G*clust_mass))*no_stars/(8*np.log(0.1*no_stars)))

def tdyn_calc(set):
    """
    Function to compute the dynamical timescale for each particle relative to one another.
    In all its permutations, it only keeps the minimal value.
    
    Inputs:
    set:    The particle set currently being evolved
    """
    
    tdyn_array = [ ]
    for i in range(len(set)):
        tdyn_temp_arr = [ ]
        for j in range(len(set)):
            if i == j:
                pass
            else:
                dist_temp = abs(set[i].position.length()-set[j].position.length())
                value = ((4*np.pi*dist_temp**3)/(3*constants.G*(set[i].mass+set[j].mass))).sqrt()
                value = value.value_in(units.yr)
                tdyn_temp_arr.append(value)
        tdyn_array.append(min(tdyn_temp_arr))

    return tdyn_array