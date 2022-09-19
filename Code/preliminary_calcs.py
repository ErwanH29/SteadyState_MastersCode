from amuse.lab import *
from amuse.units import units
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from amuse.ext.galactic_potentials import MWpotentialBovy2015
from initialiser import *
import numpy as np
import matplotlib.pyplot as plt

def prelim_dynamical_friction():
    """
    Function to form a massive data set of a hypothetical GC.
    Used in dynamical friction calculations for the enclosed mass
    and velocity fraction population.
    """
    
    conv = nbody_system.nbody_to_si(10**7 | units.MSun, 3 | units.parsec)
    m_stars = new_powerlaw_mass_distribution(10**7, 0.01 | units.MSun, 500 |units.MSun, 
                                             alpha = - 2.50)
    dynamicf_calc = new_plummer_model(10**7, convert_nbody = conv)
    dynamicf_calc.mass = m_stars
    dynamicf_mass = []
    dynamicf_pos  = []
    dynamicf_vel  = []
    for i in range(len(dynamicf_calc)):
        dynamicf_mass.append(dynamicf_calc[i].mass.value_in(units.MSun))
        dynamicf_pos.append(dynamicf_calc[i].position.length().value_in(units.AU))
        dynamicf_vel.append(dynamicf_calc[i].velocity.length().value_in(units.kms))
    dynamicf_mass = np.array(dynamicf_mass)
    dynamicf_pos = np.array(dynamicf_pos)
    dynamicf_vel = np.array(dynamicf_vel)
    indices = dynamicf_pos.argsort()

    sortedpos_array  = dynamicf_pos[indices]
    sortedmass_array = dynamicf_mass[dynamicf_pos.argsort()]
    sortedmass_array = np.cumsum(dynamicf_mass)

    massdist_Array = pd.DataFrame()
    df_massdist = pd.Series({'Position': sortedpos_array, 
                             'Mass': sortedmass_array, 
                             'Velocity': dynamicf_vel})
    massdist_Array = massdist_Array.append(df_massdist, ignore_index=True)
    massdist_Array.to_pickle('data/preliminary_calcs/Mass_Dist_Distr.pkl')

prelim_dynamical_friction()
