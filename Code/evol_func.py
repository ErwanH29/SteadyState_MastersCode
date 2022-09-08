from amuse.lab import *
from amuse.units import units
from initialiser import *
from evol_func import *
import numpy as np


class drift_without_gravity(object):
    """
    Class which induces a kick-drift-kick mechanism so IMBH feel 
    the background potential
    """

    def __init__(self, convert_nbody, time= 0 |units.Myr):
        self.model_time = time
        self.convert_nbody = convert_nbody
        self.particles = Particles()

    def evolve_model(self, t_end):
        """
        Evolves the GC in time gravitationally. 
        Note the particles here are NOT IMBH, but sample the GC
        defined particle (see initialiser.py)
        """
        
        dt = t_end - self.model_time
        self.particles.position += self.particles.velocity*dt
        self.model_time = t_end

    @property
    def potential_energy(self):
        return quantities.zero

    @property 
    def kinetic_energy(self):
        return (0.5*self.particles.mass \
                   *self.particles.velocity.lengths()**2).sum()

    def stop(self):
        pass
                   

def adaptive_dt(eta, tend, parti):
    """
    In case one wants, the time-scale can change depending on the minimum 
    interparticle collision time. This is calculated from the unaccelerated 
    linear motion and free-fall time using Eqn (A.6) of Zwart et al. 2021 
    
    Input:
    eta:     Rigid time-step used
    tend:    The final time simulation evolves until
    parti:   The IMBH particle set
    output:  A new time-step once particles come closer than a certain distance.
    """
    
    comp_val = eta * tend
    dr_ratio = [ ]
    for i in range(len(parti)):
        for j in range(len(parti)):
            if j != i:
                dist  = (abs(parti[i].position.length()-parti[j].position.length()))
                
                if dist < 1000 * (parti[i].collision_radius + parti[j].collision_radius)/2:
                    vel   = (abs(parti[i].velocity.length()-parti[j].velocity.length()))
                    value = dist/vel

                    if value > 10**-7 | units.yr:
                        value * eta
                        dr_ratio.append(value)
                    else:
                        dr_ratio.append(comp_val)

            else:
                j += 1

    return min(comp_val, np.min(dr_ratio))

def calc_momentum(parti):
    """
    Function which calculates the momentum of the particles in collision

    Inputs:
    parti:    The colliding particles
    output:   The momentum value for each particle
    """

    mom_value = (parti.mass * parti.velocity).sum()
    return mom_value

def merge_IMBH(parti, particles_in_encounter, tcoll):
    """
    Function which merges two particles if the collision stopping condition has been met
    
    Inputs:
    parti:                    The complete particle set being simulated
    particles_in_encounter:   The particles in the collision
    tcoll:                    The time-stamp for which the particles collide at
    outputs:                  Removal of two colliding particles, while addition of the merging product
    """
    com_pos = particles_in_encounter.center_of_mass()
    com_vel = particles_in_encounter.center_of_mass_velocity()

    new_particle = Particles(1)
    if calc_momentum(particles_in_encounter[0]) > calc_momentum(particles_in_encounter[1]):
        new_particle.key_tracker = particles_in_encounter[0].key
    else: 
        new_particle.key_tracker = particles_in_encounter[1].key

    new_particle.mass = particles_in_encounter.total_mass()
    new_particle.collision_time = tcoll
    new_particle.position = com_pos
    new_particle.velocity = com_vel
    new_particle.radius = (2*constants.G*new_particle.mass)/(constants.c**2)
    new_particle.collision_radius = new_particle.radius * 3000

    
    parti.add_particles(new_particle)
    parti.remove_particles(particles_in_encounter)
    return new_particle

