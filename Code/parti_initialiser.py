from amuse.lab import *
from amuse.units import (units, constants)
from amuse.ic.plummer import new_plummer_model
from amuse.ic.kingmodel import new_king_model
from amuse.ic.scalo import new_scalo_mass_distribution
from random import random, randint, choices
import numpy as np

class MW_SMBH(object):
    """
    Class which describes the central (MW) SMBH.
    """

    def __init__(self, mass=4.e6 | units.MSun,
                 position=[0, 0, 0] | units.parsec,
                 velocity=[0, 0, 0] | (units.AU/units.yr)):
        """
        Initialising function for the SMBH class.
        """
    
        self.mass = mass
        self.position = position
        self.velocity = velocity
        self.bh_rad = (2*constants.G*mass)/(constants.c**2)

    def get_gravity_at_point(self, eps, x, y, z):
        """
        Function which gathers the gravitational acceleration induced
        at any given point.
        Inputs:
        eps:     Integrator softening parameter
        x, y, z: Cartesian coordinates of the object wanting to compute
                 potential energy of.
        """     

        dx = x - self.position.x
        dy = y - self.position.y
        dz = z - self.position.z
        radius = (dx*dx + dy*dy + dz*dz).sqrt()
        radius3 = radius*radius*radius
        fr = -constants.G*self.mass/radius3
        ax = fr*dx
        ay = fr*dy
        az = fr*dz

        return ax, ay, az

    def get_potential_at_point(self, eps, x, y, z):
        """
        Function calculating the potential of the SMBH at any point
        Inputs:
        eps:     Integrator softening parameter
        x, y, z: Cartesian coordinates of the object wanting to compute
                 potential energy of.
        """

        dx = x - self.position.x
        dy = y - self.position.y
        dz = z - self.position.z
        radius = (dx*dx + dy*dy + dz*dz).sqrt()
        phi = -constants.G*self.mass/radius

        return phi

class IMBH_init(object):
    """
    Class which initialises the IMBH population.
    """

    def __init__(self):
        self.N = 0
        self.mass = 1000 | units.MSun

        return

    def N_count(self):
        """
        Function which counts the number of particles in the simulation
        """
        
        N = self.N

        return int(N)
    
    def IMBH_radius(self, mass):
        """
        Function which sets the IMBH radius based on the Schwarzschild radius
        
        Inputs:
        mass:   The mass of the input particle
        """

        radius = (2*constants.G*mass)/(constants.c**2)

        return radius

    def coll_radius(self, radius):
        """
        Function which sets the IMBH collision radius based on Zwart et al. 2021
        """

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

        sigmaV = 6 # in kms
        return np.sqrt(2/np.pi)*(vel**2/sigmaV**3)*np.exp(-vel**2/(2*sigmaV**2))

    def velocityList(self):
        """
        Function to give a velocity for an initialised particle
        """

        vrange = np.linspace(0, 5*np.sqrt(6)) # in kms
        r=[-1,1]
        w = self.ProbFunc(vrange)
        scalex = [np.random.choice(r)]
        scaley = [np.random.choice(r)]
        scalez = [np.random.choice(r)]
        vx = np.array(choices(vrange, weights=w, k = 1))*scalex
        vy = np.array(choices(vrange, weights=w, k = 1))*scaley
        vz = np.array(choices(vrange, weights=w, k = 1))*scalez
        velocity = np.concatenate((vx,vy,vz))
        return velocity

    def kroupa_mass(self):
        return new_kroupa_mass_distribution(1, 50 | units.MSun, 10**5 | units.MSun)

    def salpeter_mass(self):
        alpha = -2.35
        return new_powerlaw_mass_distribution(1, 50 | units.MSun, 10**5 | units.MSun, alpha)

    def scalo_mass(self):
        return new_scalo_mass_distribution(1, 50 | units.MSun, 10**5 | units.MSun)

    def custom_mass(self, constant):
        return constant

    def plummer_distr(self, converter):
        N = 100
        return new_plummer_model(N, radius_cutoff = 50, convert_nbody = converter)

    def king_distr(self, converter):
        N = 100
        beta = -9
        return new_king_model(N, W0 = beta, convert_nbody = converter)

    def IMBH_first(self, init_dist, converter):
        """
        Function to initialise the first IMBH population.
        The first particle forms the center of the cluster

        Inputs:
        init_dist:  The distance from central SMBH the cluster will be
        converter:  Converter used to translate nbody_system to SI units
        """
        
        SMBH_parti = MW_SMBH()
        IMBH = Particles(3)
        self.N += len(IMBH)

        IMBH[0].position = SMBH_parti.position
        IMBH[0].velocity = SMBH_parti.velocity
        IMBH[0].mass     = SMBH_parti.mass
        
        for i in range(self.N):
            #IMBH[i].mass = self.mass_func()
            IMBH[i].position = 0.85*self.plummer_distr(converter)[randint(0,self.N)].position
            IMBH[i].velocity = self.velocityList() * (1 | units.AU/units.yr)

        IMBH[1].position  = [0, 0, 0] | units.AU
        IMBH[1].velocity  = [0, 0, 0] | units.AU/units.yr

        for i in range(self.N-2):
            velx_vect = float((IMBH[i+2].x - IMBH[1].x).value_in(units.AU))
            vely_vect = float((IMBH[i+2].y - IMBH[1].y).value_in(units.AU))
            velz_vect = float((IMBH[i+2].z - IMBH[1].z).value_in(units.AU))
        vel_vect  = [velx_vect, vely_vect, velz_vect]
        veldist   = np.sqrt((velx_vect**2+vely_vect**2+velz_vect**2))
        
        IMBH[2:].velocity = -1 * IMBH[2:].velocity * (vel_vect)/(veldist)
        IMBH[1:].mass     = self.mass
        IMBH[1:].x += init_dist
        IMBH[1:].vy += 1.15*(constants.G*SMBH_parti.mass/IMBH.position.length()).sqrt()
        IMBH.radius = self.IMBH_radius(IMBH.mass)
        IMBH.collision_radius = self.coll_radius(IMBH.radius)
        IMBH.key_tracker = IMBH.key
        IMBH.ejection = 0
        IMBH.move_to_center()
        
        return IMBH

    def add_IMBH(self, pos, converter):
        """
        Function which adds an IMBH particle to the cluster
        
        Inputs:
        pos:       The current c.o.m position of the cluster
        converter: The converter to go between SI and Nbody units
        """

        self.N += 1
        add_IMBH = Particles(1)
        add_IMBH.mass = self.mass
        add_IMBH.position = 2 * self.plummer_distr(converter)[randint(0,self.N)].position
        add_IMBH.position += pos
        add_IMBH.velocity = self.velocityList() * (1 | units.AU/units.yr)
        add_IMBH.velocity *= (1 * add_IMBH.position)/(add_IMBH.position.length()) 
        add_IMBH.key_tracker = add_IMBH.key
        add_IMBH.radius = self.IMBH_radius(add_IMBH.mass)
        add_IMBH.collision_radius = self.coll_radius(add_IMBH.radius)
        add_IMBH.ejection = 0

        return add_IMBH