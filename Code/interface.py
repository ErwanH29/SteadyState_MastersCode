from parti_initialiser import *
from plotters import *
from evol import *
from amuse.community.brutuspn.interface import Brutus

"""
pop_choice   = int(input('How many IMBH do you wish to simulate?  '))
birth_param  = float(input('How often should IMBH appear \'near\' the cluster [years]?  ')) | units.yr
tend = pop_choice * birth_param + 10 * birth_param

print('Simulation will evolve until: ', tend.value_in(units.yr), 'yr.')

mass_choice = input('Which mass function to sample? \n'
                    'Constants (C) \n'
                    'Salpeter  (S) \n'
                    'Scalo     (Scalo) \n'
                    'Kroupa    (K) \n')
if mass_choice == 'S' or mass_choice == 's':
    alpha = -abs(float(input('What power law do you wish to simulate the distribution with? \n'
                             'Standard (Salpeter):     -2.35 \n'
                             'Custom choice:        (-1.0,-5.0]\n')))
if mass_choice =='C' or mass_choice == 'c':
    alpha = abs(float(input('Choose a value for which the IMBH masses are [in MSun]?)))'

else:
    alpha = None

distr_choice = input('Which positional distribution function to sample? \n'
                  'Plummer (P) \n'
                  'King    (K) \n')

if distr_choice == 'K' or distr_choice == 'k':
    beta = abs(float(input('Pick a potential depth (W0):   ')))

init_dist = float(input('Where do you wish to simulate the cluster (distance from SMBH in parsecs)? '))

"""

tend = 10**6 | units.yr
eta = 10**-5
cluster_mass = 10**7  | units.MSun
cluster_radi = 10**-3 | units.parsec
cluster_dist = 10**-2 | units.parsec
conv = nbody_system.nbody_to_si(cluster_mass, cluster_radi)

IMBH_code = IMBH_init()
IMBH_parti = IMBH_code.IMBH_first(mass_string = 'S', dist_string = 'P', alpha = -2.35,
                                  init_dist = cluster_dist, converter = conv)

evolve_system(IMBH_parti, tend, eta, cluster_dist, conv)
print('...Plotting Figures...')
spatial_plotter(cluster_dist)
energy_plotter()

anim = True
if (anim):
    animator(tend, eta)