from amuse.lab import *
from amuse.units import (units, constants, quantities)
from amuse.ic.plummer import new_plummer_model
from amuse.ic.kingmodel import new_king_model
from amuse.ic.scalo import new_scalo_mass_distribution
from amuse.ext.galactic_potentials import Plummer_profile
from random import random, randint

class GC_pot(object):
    """
    Class which defines the GC you want to simulate
    """

    def __init__(self, gc_rad  = 50 | units.AU):
        """
        Initialising function.
        outputs: The GC mass, potential profile and positions. 
        """

        self.gc_mass = 2e5 | units.MSun
        self.plum = Plummer_profile(self.gc_mass, gc_rad)
        self.gc_pos  = [0, 0, 0] | units.parsec

    def d_update(self, x_c, y_c, z_c):
        """
        Function which defines the cluster center
        outputs: Array with x, y, z position of the cluster center
        """

        self.gc_pos[0] = x_c
        self.gc_pos[1] = y_c
        self.gc_pos[2] = z_c
        
    def get_potential_at_point(self, eps, x, y, z):
        """
        Function which gathers the potential of the GC at any given point
        output: The potential at any given point
        """

        return self.plum.get_potential_at_point(eps, 
                                        x-self.gc_pos[0], 
                                        y-self.gc_pos[1], 
                                        z-self.gc_pos[2])
    
    def get_gravity_at_point(self, eps, x, y, z):
        """
        Function which gathers the gravitational acceleration induced
        at any given point
        output: The gravitational acceleration induced onto the particle
        """

        ax_p, ay_p, az_p = self.plum.get_gravity_at_point(eps, 
                                        x-self.gc_pos[0], 
                                        y-self.gc_pos[1], 
                                        z-self.gc_pos[2])
        return ax_p, ay_p, az_p
        
class MW_SMBH(object):
    """
    Class which describes the central (MW) SMBH.
    """

    def __init__(self,
                 bh_mass=4.31e6 | units.MSun,
                 position=[0, 0, 0] | units.parsec):
        """
        Initialising function.
        outputs: The SMBH mass, radius and position. 
        """
    
        self.bh_mass = bh_mass
        self.position = position
        self.bh_rad = (2*constants.G*bh_mass)/(constants.c**2)

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
        fr = -constants.G*self.bh_mass/radius3
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
        phi = -constants.G*self.bh_mass/radius
        return phi

class IMBH_init(object):
    """
    Class which initialises the IMBH population (starts with N = 2 for Mikkola)
    """

    def __init__(self):
        self.N = 3

    def N_count(self):
        """
        Function to keep track of the number of IMBH
        """

        N = self.N
        return int(N)
    
    def IMBH_radius(self, parti):
        """
        Function which sets the IMBH radius based on the Schwarzschild radius
        
        Inputs:
        parti:   The particle set being simulated
        output:  The particle set with a refined radius
        """

        parti.radius = (2*constants.G*parti.mass)/(constants.c**2)
        return parti

    def decision(self, time):
        """
        Function to define if a IMBH is kicked into the cluster or not

        Inputs:
        time:    The current time of the simulation
        outputs: Decision whether to throw a IMBH into the mix
        """

        step_size = 10**-3
        birth_param = 10
        if time == 0.0 | units.Myr:
            c = False
        else:
            c = random.random() < step_size/(birth_param)

    def velocity_func(self):
        """
        Function to initialise the velocity of the IMBH sinking towards center
        output: The velocity of the particle
        """

        return [6.227667635500401, 12.998565566927823, 0] * (1 | units.AU/units.yr)

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
        
        N = 1000         # Hard-coded value to show how many IMBH in simulation

        if dist_string == 'P' or dist_string == 'p':
            distributer = new_plummer_model(N, convert_nbody = converter, 
                                            radius_cutoff = 15)
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

        SMBH_code = MW_SMBH()
        r = init_dist | units.parsec
        self.N = 2

        IMBH = Particles(3)
        
        for i in range(len(IMBH)):
            IMBH[i].mass = self.mass_func(mass_string,alpha)
            IMBH[i].position = self.IMBH_posinit(100, dist_string, 5, converter)[randint(0,100)].position
        IMBH[0].position  = self.IMBH_posinit(100, dist_string, 5, converter)[0].position  
          
        IMBH.position += (1, 0, 0) * r

        IMBH.velocity  = [[23.518569826682455, 10.52227851919019, 0],    #Still needs development
                          [-20.74623746218285, -22.52084408611801, 0],#] | units.AU/units.yr
                          [-10.74623746218285, 15.52084408611801, 0], ] | units.AU/units.yr

        IMBH.velocity += (0, 1, 0) * (constants.G*SMBH_code.bh_mass/r).sqrt()
        IMBH.key_tracker = IMBH.key

        return IMBH


    def add_part(self, mass_func, init_dist):
        """
        Function to initialise parameters of IMBH who come in time to the cluster

        Inputs:
        mass_func:   The mass function used to sample the particles from
        init_dist:   Used to ensure the added particle is near the cluster
        output:      A new particle to add to the simulation
        """

        rand_place = 1000 * random.random() | units.parsec
        self.N += 1
        add_IMBH = Particles(1)
        add_IMBH.mass      = mass_func
        add_IMBH.velocity  = self.velocity_func() * (1 | units.AU / units.yr)
        add_IMBH.position  = init_dist * rand_place

        return add_IMBH

class GC_init(object):
    """
    Class which initiates the GC position
    """

    def GC_tracer(self, parti):
        """
        Function which tracks the position of the GC
        
        Inputs:
        parti:    The particle set for which the first particle describes
                  GC center
        outputs:  Initial kinematic property of the GC
        """

        gc_parti  = Particles(1)
        gc_parti.mass = 1.7e5 | units.MSun #Choose total masses of the cluster (Also line 11)
        gc_parti.position = parti[0].position
        gc_parti.velocity = parti[0].velocity

        return gc_parti
    
    def gal_path_init(self):
        """
        Function which tracks the cartesian coordinates of the GC
        outputs: Cartesian coordinates of the GC
        """

        x_gc = [] | units.kpc
        y_gc = [] | units.kpc
        z_gc = [] | units.kpc

        return x_gc, y_gc, z_gc