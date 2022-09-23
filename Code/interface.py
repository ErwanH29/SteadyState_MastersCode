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
    file_reset('data/stability_time')
    file_reset('figures')

tend = 10**5 | units.yr
eta  = 10**-4
cluster_mass = 10**7  | units.MSun
cluster_radi = 10**-3 | units.parsec
cluster_dist = 0.03 | units.parsec
conv = nbody_system.nbody_to_si(cluster_mass, cluster_radi)
no_sim = 1000

for i in range(no_sim):
    print('=========== Simulation '+str(i+1)+'/'+str(no_sim)+' Running ===========')
    IMBH_code = IMBH_init()
    IMBH_parti = IMBH_code.IMBH_first(init_dist = cluster_dist, converter = conv)
    failed_simul = evolve_system(IMBH_parti, tend, eta, cluster_dist, cluster_radi , conv)
    if (failed_simul):
        pass

    else:
        print('...Plotting Figures...')
        spatial_plotter(1.25*cluster_dist)
        energy_plotter()

        anim = False
        if (anim):
            animator(2*cluster_dist)

steadytime_plotter()