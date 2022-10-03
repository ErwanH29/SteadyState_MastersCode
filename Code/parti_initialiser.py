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
                 velocity = [0, 0, 0] | (units.AU/units.yr)):
    
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

class globular_cluster(object):
    """
    Class which defines the GC you want to simulate
    """

    def __init__(self, no_stars = 1e5,
                 mass = 1.e6 | units.MSun,
                 cluster_radi = 1e-1 | units.parsec,
                 cluster_dist = 1.4 | units.parsec,
                 position = [0, 0, 0] | units.parsec,
                 velocity = [0, 0, 0] | (units.AU/units.yr)):

        self.gc_pop = no_stars
        self.gc_mass = mass
        self.gc_rad  = cluster_radi
        self.gc_dist = cluster_dist
        self.gc_pos = position
        self.gc_velocity = velocity
        self.gc_conv = nbody_system.nbody_to_si(mass, cluster_radi)

        self.plum = Plummer_profile(self.gc_mass, self.gc_rad)

    def d_update(self, x_c, y_c, z_c):
        """
        Function which defines and updates the cluster center
        """

        self.gc_pos[0] = x_c
        self.gc_pos[1] = y_c
        self.gc_pos[2] = z_c
        
    def get_potential_at_point(self, eps, x, y, z):
        return self.plum.get_potential_at_point(eps, 
                                        x-self.gc_pos[0], 
                                        y-self.gc_pos[1], 
                                        z-self.gc_pos[2])

class IMBH_init(object):
    def __init__(self):
        self.N = 0
        self.mass = 1000 | units.MSun

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
        return 10*radius

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

        sigmaV = 15 # in kms

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

    def plummer_distr(self, N):
        gc_code = globular_cluster()
        distr = new_plummer_model(N, convert_nbody = gc_code.gc_conv)
        rhmass = LagrangianRadii(distr)[6].in_(units.parsec)
        return distr, rhmass

    def king_distr(self, converter):
        N = 20
        beta = -9
        return new_king_model(N, W0 = beta, convert_nbody = converter)

    def IMBH_first(self, init_parti):
        """
        Function to initialise the first IMBH population.
        The first particle forms the center of the cluster

        Inputs:
        init_parti: The (2+N) number of IMBH particles you wish to simulate
        converter:  Converter used to translate nbody_system to SI units
        """
        
        SMBH_parti = MW_SMBH()
        gc_code = globular_cluster()
        self.N += init_parti+2

        gc_particles, rhmass = self.plummer_distr(10**4)
        gc_particles.velocity = self.velocityList() * (1 | units.AU/units.yr)
        gc_particles.ejection = 0
        gc_particles.coll_events = 0

        gc_particles[:self.N].radius = self.IMBH_radius(gc_particles[:self.N].mass)
        gc_particles[:self.N].collision_radius = self.coll_radius(gc_particles[:self.N].radius)
        gc_particles[:self.N].key_tracker = gc_particles[:self.N].key
        gc_particles[2:self.N].mass = self.mass
        gc_particles[0:2].velocity  = [0, 0, 0] | units.AU/units.yr
        for part_ in gc_particles[2:self.N]:
            if part_.x > gc_code.gc_rad:
                part_.x = gc_code.gc_rad
            if part_.y > gc_code.gc_rad:
                part_.y = gc_code.gc_rad
            if part_.z > gc_code.gc_rad:
                part_.z = gc_code.gc_rad
        gc_particles[1].position = gc_particles.center_of_mass()
        gc_particles.scale_to_standard(convert_nbody=gc_code.gc_conv)

        gc_particles[2:self.N].mass = self.mass
        gc_particles[0].mass = SMBH_parti.mass
        gc_particles[0].radius = self.IMBH_radius(gc_particles[0].mass)
        gc_particles[0].collision_radius = self.coll_radius(gc_particles[0].radius)        
        gc_particles[0].position = [0, 0, 0] | units.AU
        gc_particles[0].velocity = [0, 0, 0] | units.AU/units.yr

        gc_particles[1].mass = gc_code.gc_mass
        gc_particles[1].radius = 0 | units.m
        gc_particles[1].collision_radius = 0 | units.m

        gc_particles[1:].x += gc_code.gc_dist
        for part_ in gc_particles[1:]:
            part_.vy += (constants.G*SMBH_parti.mass/gc_code.gc_dist).sqrt()
        sim_particles = gc_particles[:self.N]

        return sim_particles, rhmass

    def add_IMBH(self, globular):
        """
        Function which adds an IMBH particle to the cluster
        
        Inputs:
        pos:       The current c.o.m position of the cluster
        """

        gc_code = globular_cluster()

        self.N += 1
        add_IMBH = Particles(1)
        add_IMBH.mass = self.mass
        add_IMBH.position = gc_code.gc_rad * [np.random.uniform(0.85,1)*np.random.choice([-1,1]), 
                                              np.random.uniform(0.85,1)*np.random.choice([-1,1]), 
                                              np.random.uniform(0.85,1)*np.random.choice([-1,1])]
        add_IMBH.position += globular.position
        add_IMBH.velocity = self.velocityList() * (1 | units.AU/units.yr)
        add_IMBH.velocity = np.sqrt(abs(gc_code.get_potential_at_point(0, add_IMBH.x, add_IMBH.y, add_IMBH.z)) )* (add_IMBH.velocity)/(add_IMBH.velocity.length())
        add_IMBH.velocity += globular.velocity
        add_IMBH.key_tracker = add_IMBH.key
        add_IMBH.radius = self.IMBH_radius(add_IMBH.mass)
        add_IMBH.collision_radius = self.coll_radius(add_IMBH.radius)
        add_IMBH.ejection = 0
        add_IMBH.coll_events = 0

        return add_IMBH