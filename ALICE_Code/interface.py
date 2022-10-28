from parti_initialiser import *
from file_logistics import *
from evol import *

eta  = 1e-4
tend = 100 | units.Myr

SMBH_code = MW_SMBH()
int_string = 'Hermite'

pops = [10, 20]
sims = [20, 20]
seeds = [2929, 8888, 3600, 2706, 1514, 3781, 5482, 9399, 9999, 3212, 3583, 5911, 5767, 2312, 54314, 30125, 65863, 5812, 32131, 882988]

for k in range(3):
    for ipop_, sims_, seed_ in zip(pops, sims, seeds):
        initial_pop = ipop_
        for k in range(sims_):
            print('=========== Simulation '+str(k)+'/'+str(sims_)+' Running ===========')
            IMBH_code = IMBH_init()
            code_conv = nbody_system.nbody_to_si((ipop_*IMBH_code.mass + SMBH_code.mass), SMBH_code.distance)
            IMBH_parti, rhmass = IMBH_code.IMBH_first(initial_pop, seed_)
            failed_simul = evolve_system(IMBH_parti, tend, eta, SMBH_code.distance, code_conv, int_string)

#from amuse.ext.orbital_elements import orbital_elements_from_binary
#print(dir(orbital_elements_from_binary))
