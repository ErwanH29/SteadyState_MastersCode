from amuse.lab import *
from amuse.units import units
from parti_initialiser import *
from evol_func import *
import numpy as np

def indiv_PE(indiv, set, array):
    for j in range(len(set)):
        if indiv == set[j]:
            pass
        else:
            distance = (indiv.position.length()-set[j].position.length())
            temp_PE  = ((constants.G*indiv.mass*set[j].mass)/abs(distance))
            array.append(temp_PE)

    return array

def tdyn_calc(particle):
    """
    Function to compute the dynamical timescale for each particle relative to one another.
    In all its permutations, it only keeps the minimal value.
    
    Inputs:
    particle: The particle set currently being evolved
    outputs:  The dynamical timescale
    """
    
    tdyn_array = [ ]
    for i in range(len(particle)):
        tdyn_temp_arr = [ ]
        for j in range(len(particle)):
            if i == j:
                pass
            else:
                dist_temp = abs(particle[i].position.length()-particle[j].position.length())
                value = ((4*np.pi*dist_temp**3)/(3*constants.G*(particle[i].mass+particle[j].mass))).sqrt()
                value = value.value_in(units.yr)
                tdyn_temp_arr.append(value)
        tdyn_array.append(min(tdyn_temp_arr))

    return tdyn_array