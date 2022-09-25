from parti_initialiser import *
from plotters import *
from file_logistics import *
from evol import *

remove_files = True
if (remove_files):
    file_reset('data/center_of_mass')
    file_reset('data/dynamical_time')
    file_reset('data/energy')
    file_reset('data/lagrangians')
    file_reset('data/particle_energies')
    file_reset('data/positions_IMBH')
    file_reset('data/simulation_stats')
    #file_reset('data/stability_time')
    file_reset('figures')

tend = 10**7 | units.yr
eta  = 10**-6
cluster_mass = 10**7  | units.MSun
cluster_radi = 3e-3 | units.parsec
cluster_dist = 0.1 | units.parsec
conv = nbody_system.nbody_to_si(cluster_mass, cluster_radi)
no_sim = 4000

for j in [3, 4, 5, 6, 7, 8, 9, 10]:
    initial_pop = j
    file_reset('data/center_of_mass')
    file_reset('data/dynamical_time')
    file_reset('data/energy')
    file_reset('data/lagrangians')
    file_reset('data/particle_energies')
    file_reset('data/positions_IMBH')
    file_reset('data/simulation_stats')
    file_reset('data/stability_time')
    file_reset('figures')

    if j < 5:
        eta = 1e-6
    else:
        eta = 1e-7

    for i in range(no_sim):
        print('=========== Simulation '+str(i+1)+'/'+str(no_sim)+' Running ===========')
        IMBH_code = IMBH_init()
        IMBH_parti = IMBH_code.IMBH_first(cluster_dist, initial_pop, conv)
        failed_simul = evolve_system(IMBH_parti, tend, eta, cluster_dist, cluster_radi , conv)
        if (failed_simul):
            pass

        else:
           # print('...Plotting Figures...')
            #spatial_plotter(1.25*cluster_dist)
            #energy_plotter()

            anim = False
            if (anim):
                animator(2*cluster_dist)

steadytime_plotter()