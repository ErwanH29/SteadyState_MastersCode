from amuse.lab import *
from amuse.units import (units, constants, quantities)
from amuse.ic.plummer import new_plummer_model
from amuse.ic.kingmodel import new_king_model
from amuse.ic.scalo import new_scalo_mass_distribution
from amuse.ext.galactic_potentials import Plummer_profile
from random import random, randint, choices
from amuse.community.galaxia.interface import BarAndSpirals3D
import numpy as np

class MW_SMBH(object):
    """
    Class which describes the central (MW) SMBH.
    """

    def __init__(self,
                 mass=4.e6 | units.MSun,
                 position=[0, 0, 0] | units.parsec,
                 velocity=[0, 0, 0] | (units.AU/units.yr)):
        """
        Initialising function.
        outputs: The SMBH mass, radius and position. 
        """
    
        self.mass = mass
        self.position = position
        self.velocity = velocity
        self.bh_rad = (2*constants.G*mass)/(constants.c**2)

    def get_gravity_at_point(self, eps, x, y, z):
        """
        Function which gathers the gravitational acceleration induced
        at any given point
        output: The gravitational acceleration induced onto the particle
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
        output: The potential at any given point
        """

        dx = x - self.position.x
        dy = y - self.position.y
        dz = z - self.position.z
        radius = (dx*dx + dy*dy + dz*dz).sqrt()
        phi = -constants.G*self.mass/radius
        return phi

class IMBH_init(object):
    """
    Class which initialises the IMBH population (starts with N = 2 for Mikkola)
    """

    def __init__(self):
        self.N = 0
        self.mass = 1000 | units.MSun
        return

    def N_count(self):
        N = self.N
        return int(N)
    
    def IMBH_radius(self, mass):
        """
        Function which sets the IMBH radius based on the Schwarzschild radius
        
        Inputs:
        parti:   The particle set being simulated
        output:  The particle set with a refined radius
        """

        radius = (2*constants.G*mass)/(constants.c**2)
        return radius

    def coll_radius(self, radius):
        """
        Function which sets the IMBH collision radius based on Zwart et al. 2021
        """

        return 10*radius

    def decision(self, time, app_rate):
        if time == 0.0 | units.Myr: # because at 0.0 particles already appear
            c = False
        else:
            c = random() < (app_rate) # from  1/(df_timescale)

        return c

    def ProbFunc(self, vel):
        sigmaV = 6 # in kms
        return np.sqrt(2/np.pi)*(vel**2/sigmaV**3)*np.exp(-vel**2/(2*sigmaV**2))

    def velocityList(self):
        vrange = np.linspace(0, 500/3) # in kms
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

    def mass_func(self, mass_string, alpha):
        """
        Function to define the mass of BHs present in simulation
        
        Inputs:
        mass_string: The mass function wanting to simulate.
        alpha:       The power-law/constant value to model the distribution from
        output:      Defines the mass of the particle
        """

        if mass_string == 'C':              # Can change these values based on simulation trial (Range [1e2, 1e5])
            mass_bh = alpha    | units.MSun
            return mass_bh

        if mass_string == 'S' or mass_string == 's':
            mass_distr = new_powerlaw_mass_distribution(1, 50 | units.MSun, 
                                                        10**5 | units.MSun, alpha)
            return mass_distr

        if mass_string == 'Scalo' or mass_string == 'scalo':
            mass_distr = new_scalo_mass_distribution(1, 50 | units.MSun, 
                                                     10**5 | units.MSun)
            return mass_distr

        if mass_string == 'K' or mass_string == 'k':
            mass_distr = new_kroupa_mass_distribution(1, 50 | units.MSun, 
                                                      10**5 | units.MSun)
            return mass_distr

        else:
            return print('Error when choosing mass function')

    def IMBH_posinit(self, N, dist_string, beta, converter):
        """
        Function which sets the distribution of the IMBH particles when initialised
        
        Inputs:
        N:          Number of particles. Hard-coded to a large value for greater variety
        converter:  Converter used to translate nbody_system to SI units
        beta:       Parameter needed for the dimensionless
        output:     The initial position of the particle
        """

        N = 100
        
        if dist_string == 'P' or dist_string == 'p':
            distributer = new_plummer_model(N, convert_nbody = converter)
            return distributer

        if dist_string == 'K' or dist_string == 'k':
            distributer = new_king_model(N, W0 = beta, convert_nbody = converter)
            return distributer

        else:
            return print('Error when choosing a distribution function.')

    def IMBH_first(self, mass_string, dist_string, alpha, init_dist, converter):
        """
        Function to initialise the first IMBH population
        The first particle forms the center of the cluster

        Inputs:
        mass_string:  The mass distribution wanting to be modelled
        dist_string:  The spatial distribution wanting to be modelled
        alpha:        The power law, influencing the mass-distribtion
        init_dist:    The distance from central SMBH the cluster will be
        converter:    Converter used to translate nbody_system to SI units
        output:       The initial particle set for the simulation (as of now, N0 = 3)
        """
        
        SMBH_parti = MW_SMBH()
        IMBH = Particles(4)
        self.N += len(IMBH)

        IMBH[0].position = SMBH_parti.position
        IMBH[0].velocity = SMBH_parti.velocity
        IMBH[0].mass     = SMBH_parti.mass

        for i in range(self.N):
            #IMBH[i].mass = self.mass_func(mass_string,alpha)
            IMBH[i].position = self.IMBH_posinit(self.N, dist_string, 5, converter)[randint(0,self.N)].position
            IMBH[i].velocity = self.velocityList() * (1 | units.AU/units.yr)

        IMBH[1].position = [0, 0, 0] | units.AU
        IMBH[1].velocity = [0, 0, 0] | units.AU/units.yr
        IMBH[1:].mass    = self.mass
        IMBH[1:].x += init_dist
        IMBH[1:].vy += (constants.G*SMBH_parti.mass/init_dist).sqrt()
        IMBH.radius = self.IMBH_radius(IMBH.mass)
        IMBH.collision_radius = self.coll_radius(IMBH.radius)
        IMBH.key_tracker = IMBH.key
        IMBH.move_to_center()
        
        print(IMBH)
        return IMBH

    def add_IMBH(self, pos, distance):
                
        SMBH_parti = MW_SMBH()
        self.N += 1
        add_IMBH = Particles(1)
        add_IMBH.mass = self.mass
        add_IMBH.position = pos
        add_IMBH.velocity = self.velocityList() * (1 | units.AU/units.yr)
        add_IMBH.vy += (constants.G*SMBH_parti.mass/distance).sqrt()
        add_IMBH.key_tracker = add_IMBH.key
        add_IMBH.radius = self.IMBH_radius(add_IMBH.mass)
        add_IMBH.collision_radius = self.coll_radius(add_IMBH.radius)
        
        return add_IMBH