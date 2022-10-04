from parti_initialiser import *
from spatial_plotters import *
from steady_plotters import *
from file_logistics import *
from evol import *

eta  = 1e-4
tend = 10e6 | units.yr

SMBH_code = MW_SMBH()
gc_code = globular_cluster()
code_conv = nbody_system.nbody_to_si((gc_code.gc_mass+SMBH_code.mass), gc_code.gc_dist)

"""prompt = input(('WARNING: About to delete all files. Are you sure (y|n)?'))
if prompt == 'y':
    file_reset('data/stable_simulation')
    file_reset('data/chaotic_simulation')
    file_reset('data/stability_time')"""

no_sim = 2000

initial_pop = 3
remove_files = True

if (remove_files):
    file_reset('data/center_of_mass')
    #file_reset('data/chaotic_simulation')
    file_reset('data/collision_events')
    file_reset('data/dynamical_time')
    file_reset('data/energy')
    file_reset('data/event_tracker')
    file_reset('data/lagrangians')
    file_reset('data/particle_energies')
    file_reset('data/positions_IMBH')
    file_reset('data/collision_events')
    #file_reset('data/simulation_stats')
    #file_reset('data/stable_simulation')
    file_reset('figures')

for j in [3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]:
    initial_pop = j
    for i in range(no_sim):
        print('=========== Simulation '+str(i+1)+'/'+str(no_sim)+' Running ===========')
        IMBH_code = IMBH_init()
        IMBH_parti, rhmass = IMBH_code.IMBH_first(initial_pop)
        failed_simul = evolve_system(IMBH_parti, tend, eta, gc_code.gc_dist, gc_code.gc_rad, 
                                        gc_code.gc_mass, rhmass, code_conv)
        if (failed_simul):
            pass

        else:
            #print('...Plotting Figures...')
            #spatial_plotter(1.25*gc_code.gc_dist)
            #energy_plotter()

            anim = False
            if (anim):
                animator(2*gc_code.gc_dist)