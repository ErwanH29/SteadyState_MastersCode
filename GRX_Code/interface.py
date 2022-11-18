from parti_initialiser import *
from file_logistics import *
from evol import *

"""
9999, 3212, 3213213, 594210, 567810, 213210, 
943218, 8320021, 1010, 4206988666,  3583, 
5911, 5767, 2312, 54314, 30125, 65863, 5812, 
32131, 882988
"""

SMBH_code = MW_SMBH()

eta  = 1e-5
tend = 100 | units.Myr
int_string = 'GRX'

pops = [10]
seeds = [567810, 213210, 
         943218, 8320021, 1010, 4206988666,  3583, 
         5911, 5767, 2312, 54314, 30125, 65863, 5812, 
         32131, 882988]

for ipop_ in pops:
    iter = 0
    for seed_ in seeds:
        iter += 1
        print('=========== Simulation '+str(iter)+'/'+str(len(seeds))+' Running ===========')
        IMBH_code = IMBH_init()
        code_conv = nbody_system.nbody_to_si((ipop_*IMBH_code.mass + SMBH_code.mass), SMBH_code.distance)
        IMBH_parti, rhmass = IMBH_code.IMBH_first(ipop_, seed_)
        
        if int_string == 'GRX':
            SMBH = Particles(1)
            SMBH.mass = IMBH_parti[0].mass
            SMBH.velocity = IMBH_parti[0].velocity
            SMBH.position = IMBH_parti[0].position
            SMBH.key_tracker = IMBH_parti[0].key_tracker
            SMBH.collision_radius = 2*IMBH_parti[0].collision_radius
            SMBH.radius = IMBH_parti[0].radius
            SMBH.ejection = 0
            SMBH.collision_events = 0
            failed_simul = evolve_system(IMBH_parti, tend, eta, SMBH_code.distance, 
                                         code_conv, int_string, SMBH)

        else:
            failed_simul = evolve_system(IMBH_parti, tend, eta, SMBH_code.distance, 
                                         code_conv, int_string, None)