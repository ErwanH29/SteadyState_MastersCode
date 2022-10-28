from amuse.lab import *
from amuse.units import (units, constants)
from amuse.ic.plummer import new_plummer_model
from amuse.ic.kingmodel import new_king_model
from amuse.ext.galactic_potentials import Plummer_profile
from amuse.ic.scalo import new_scalo_mass_distribution
from random import random, randint, choices
import numpy as np
from amuse.ext.LagrangianRadii import LagrangianRadii

class MW_SMBH(object):
    """
    Class which defines the central SMBH
    """
    def __init__(self, mass = 4.e6 | units.MSun,
                 position = [0, 0, 0] | units.parsec,
                 velocity = [0, 0, 0] | (units.AU/units.yr),
                 distance = 0.2 | units.parsec):

        self.distance = distance
        self.mass = mass
        self.position = position
        self.velocity = velocity
        self.bh_rad = (2*constants.G*mass)/(constants.c**2)

    def get_potential_at_point(self, eps, x, y, z):

        dx = x - self.position.x
        dy = y - self.position.y
        dz = z - self.position.z
        radius = (dx*dx + dy*dy + dz*dz).sqrt()
        phi = -constants.G*self.mass/radius

        return phi

class IMBH_init(object):
    def __init__(self):
        self.N = 0
        self.mass = 1000 | units.MSun                                    #Change this for different mass simulations

    def N_count(self):
        """
        Function which counts the number of particles in the simulation
        """
        return int(self.N)
    
    def IMBH_radius(self, mass):
        """
        Function which sets the IMBH radius based on the Schwarzschild radius
        
        Inputs:
        mass:   The mass of the input particle
        """
        return (2*constants.G*mass)/(constants.c**2)

    def coll_radius(self, radius):
        return 6*radius

    def decision(self, time, app_rate):
        """
        Function which chooses through a Gaussian probability whether an IMBH
        appears in the simulation.

        Inputs:
        time:     The time the simulation is evolving to
        app_rate: The computed birth rate (step*tend)/(df)
        """
        
        if time == 0.0 | units.Myr: # because at 0.0 particles already appear
            c = False
        else:
            c = random() < (app_rate) # from  1/(df_timescale)

        return c

    def ProbFunc(self, vel):
        """
        Function which initialises the velocity distribution [Maxwell distribution]
        
        Inputs:
        vel:    The velocity range for which to sample the weights from
        """

        sigmaV = 15 # in kms (if changing, change line 26 of physics_func.py)

        return np.sqrt(2/np.pi)*(vel**2/sigmaV**3)*np.exp(-vel**2/(2*sigmaV**2))

    def velocityList(self):
        """
        Function to give a velocity for an initialised particle
        """

        vrange = np.linspace(0, 50) # in kms
        r=[-1,1]
        w = self.ProbFunc(vrange)
        scalex = [np.random.choice(r)]
        scaley = [np.random.choice(r)]
        scalez = [np.random.choice(r)]
        vx = np.array(choices(vrange, weights=w, k = 1))*scalex
        vy = np.array(choices(vrange, weights=w, k = 1))*scaley
        vz = np.array(choices(vrange, weights=w, k = 1))*scalez
 
        return np.concatenate((vx,vy,vz))

    def kroupa_mass(self):
        return new_kroupa_mass_distribution(1, 50 | units.MSun, 10**5 | units.MSun)

    def salpeter_mass(self):
        alpha = -2.35
        return new_powerlaw_mass_distribution(1, 50 | units.MSun, 10**5 | units.MSun, alpha)

    def scalo_mass(self):
        return new_scalo_mass_distribution(1, 50 | units.MSun, 10**5 | units.MSun)

    def custom_mass(self, constant):
        return constant

    def plummer_distr(self, N, seed):
        np.random.seed(seed)
        distr = new_plummer_model(N, convert_nbody = self.code_conv)
        rhmass = LagrangianRadii(distr)[6].in_(units.parsec)
        return distr, rhmass

    def king_distr(self, converter):
        N = 20
        beta = -9
        return new_king_model(N, W0 = beta, convert_nbody = converter)

    def IMBH_first(self, init_parti, seed):
        """
        Function to initialise the first IMBH population.
        The first particle forms the center of the cluster

        Inputs:
        init_parti: The (2+N) number of IMBH particles you wish to simulate
        seed:       Seed for which initialises the system
        """
        
        SMBH_parti = MW_SMBH()
        self.N += init_parti+1
        sys_mass = SMBH_parti.mass + self.N * self.mass
        self.code_conv = nbody_system.nbody_to_si(self.N * self.mass, SMBH_parti.distance)
        crazy_ecc = True

        particles, rhmass = self.plummer_distr(self.N, seed)
        particles.velocity = self.velocityList() * (1 | units.AU/units.yr)
        particles.ejection = 0
        particles.coll_events = 0
        particles.z *= 0.1
        particles.key_tracker = particles.key
        particles[1:].mass = self.mass
        #particles[0].mass = SMBH_parti.mass
        particles.radius = self.IMBH_radius(particles.mass)
        particles.collision_radius = self.coll_radius(particles.radius)

        min_dist = 0.17 | units.parsec
        max_dist = 0.23 | units.parsec
        
        particles[0].position = [0, 0, 0] | units.parsec
        particles[0].velocity = [0, 0, 0] | units.ms
        for parti_ in particles[1:]:
            if parti_.position.length() < min_dist:
                parti_.position *= min_dist/parti_.position.length()
            if parti_.position.length() > max_dist:
                parti_.position *= max_dist/parti_.position.length()

        particles.scale_to_standard(convert_nbody=self.code_conv)

        if (crazy_ecc):
            for parti_ in particles[1:]:
                parti_.vx += (constants.G*SMBH_parti.mass * (abs(parti_.y)/parti_.position.length()**2)).sqrt()
                parti_.vy += (constants.G*SMBH_parti.mass * (abs(parti_.x)/parti_.position.length()**2)).sqrt()
        else:
            for parti_ in particles[1:]:
                parti_.vx = (constants.G*SMBH_parti.mass * (abs(parti_.y)/parti_.position.length()**2)).sqrt()
                parti_.vy = (constants.G*SMBH_parti.mass * (abs(parti_.x)/parti_.position.length()**2)).sqrt()

        particles[1:].mass = self.mass
        particles[0].mass = SMBH_parti.mass
        particles.radius = self.IMBH_radius(particles.mass)
        particles.collision_radius = self.coll_radius(particles.radius)
        
        return particles, rhmass

    def add_IMBH(self, globular):
        """
        Function which adds an IMBH particle to the cluster
        
        Inputs:
        pos:       The current c.o.m position of the cluster
        """

        SMBH_parti = MW_SMBH()
        self.N += 1
        add_IMBH = Particles(1)
        add_IMBH.mass = self.mass
        add_IMBH.position = [1, 1, 0.1] * SMBH_parti.distance * [np.random.uniform(0.85,1)*np.random.choice([-1,1]), 
                                                                 np.random.uniform(0.85,1)*np.random.choice([-1,1]), 
                                                                 np.random.uniform(0.85,1)*np.random.choice([-1,1])]

        add_IMBH.velocity = self.velocityList() * (1 | units.AU/units.yr)
        add_IMBH.velocity = np.sqrt(abs(gc_code.get_potential_at_point(0, add_IMBH.x, add_IMBH.y, add_IMBH.z)) )* (add_IMBH.velocity)/(add_IMBH.velocity.length())
        add_IMBH.velocity += globular.velocity
        add_IMBH.key_tracker = add_IMBH.key
        add_IMBH.radius = self.IMBH_radius(add_IMBH.mass)
        add_IMBH.collision_radius = self.coll_radius(add_IMBH.radius)
        add_IMBH.ejection = 0
        add_IMBH.coll_events = 0

        return add_IMBH
