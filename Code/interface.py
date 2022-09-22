from parti_initialiser import *
from plotters import *
from evol import *
from amuse.community.brutuspn.interface import Brutus

tend = 10**5 | units.yr
eta = 10**-4
cluster_mass = 10**7  | units.MSun
cluster_radi = 10**-3 | units.parsec
cluster_dist = 0.03 | units.parsec
conv = nbody_system.nbody_to_si(cluster_mass, cluster_radi)

for i in range(3):
    IMBH_code = IMBH_init()
    IMBH_parti = IMBH_code.IMBH_first(init_dist = cluster_dist, converter = conv)

    no_plot = evolve_system(IMBH_parti, tend, eta, cluster_dist, cluster_radi , 'P', conv)
    if (no_plot):
        pass
    else:
        print('...Plotting Figures...')
        spatial_plotter(1.25*cluster_dist)
        energy_plotter()

        anim = True
        if (anim):
            animator(2*cluster_dist)