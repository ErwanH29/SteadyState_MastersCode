from parti_initialiser import *
from file_logistics import *
from evol import *

eta  = 1e-6
tend = 1 | units.Gyr

SMBH_code = MW_SMBH()
int_string = 'Hermite'

pops = [10, 20, 30, 40, 50, 60, 70, 80, 90]
sims = [20, 20, 20, 20, 20, 20, 20, 20, 20]

for k in range(3):
    for ipop_, sims_ in zip(pops, sims):
        initial_pop = ipop_
        for k in range(sims_):
            print('=========== Simulation '+str(k)+'/'+str(sims_)+' Running ===========')
            IMBH_code = IMBH_init()
            code_conv = nbody_system.nbody_to_si((ipop_*IMBH_code.mass + SMBH_code.mass), SMBH_code.distance)
            IMBH_parti, rhmass = IMBH_code.IMBH_first(initial_pop)
            failed_simul = evolve_system(IMBH_parti, tend, eta, SMBH_code.distance, code_conv, int_string)

from amuse.ext.orbital_elements import orbital_elements_from_binary
print(dir(orbital_elements_from_binary))