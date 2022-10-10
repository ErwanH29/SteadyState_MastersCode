from parti_initialiser import *
#from spatial_plotters import *   #Comment out for local simulations in case you wish to plot
#from steady_plotters import *
from file_logistics import *
from evol import *

eta  = 5e-5
tend = 3e7 | units.yr

SMBH_code = MW_SMBH()
gc_code = globular_cluster()
code_conv = nbody_system.nbody_to_si((gc_code.gc_mass+SMBH_code.mass), gc_code.gc_dist)
no_sim = 100

initial_pop = 3
remove_files = True
int_string = 'GRX'

"""prompt = input(('WARNING: About to delete all files. Are you sure (y|n)?'))
if prompt == 'y':
    file_reset('data/'+str(int_string)+'/stable_simulation')
    file_reset('data/'+str(int_string)+'/chaotic_simulation')
    file_reset('data/'+str(int_string)+'/stability_time')
    file_reset('data/'+str(int_string)+'/simulation_stats')"""

if (remove_files):
    for string_ in ['Hermite', 'GRX']:
        file_reset('data/'+str(string_)+'/center_of_mass')
        file_reset('data/'+str(string_)+'/collision_events')
        file_reset('data/'+str(string_)+'/dynamical_time')
        file_reset('data/'+str(string_)+'/energy')
        file_reset('data/'+str(string_)+'/event_tracker')
        file_reset('data/'+str(string_)+'/lagrangians')
        file_reset('data/'+str(string_)+'/particle_energies')
        #file_reset('data/'+str(int_string)+'/particle_trajectory')
        file_reset('figures')

for j in [7, 6, 5, 4, 3, 2, 1]:
    initial_pop = j

    #if j > 8:
        #tend = 1e7 | units.yr

    for i in range(no_sim):
        print('=========== Simulation '+str(i+1)+'/'+str(no_sim)+' Running ===========')
        IMBH_code = IMBH_init()
        IMBH_parti, rhmass = IMBH_code.IMBH_first(initial_pop)
        failed_simul = evolve_system(IMBH_parti, tend, eta, gc_code.gc_dist, gc_code.gc_rad, 
                                     code_conv, int_string)
        if (failed_simul):
            pass

        else:
            #print('...Plotting Figures...')
            #spatial_plotter(1.25*gc_code.gc_dist)
            #energy_plotter()

            anim = False
            #if (anim):
            #    animator(2*gc_code.gc_dist)