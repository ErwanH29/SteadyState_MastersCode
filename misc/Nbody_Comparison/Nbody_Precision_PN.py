from amuse.lab import new_plummer_model
from amuse.units import (nbody_system, units, constants, quantities)
from amuse.community.mikkola.interface import Mikkola
from amuse.community.brutuspn.interface import Brutus
from amuse.datamodel import Particles
from amuse.datamodel.particles import Channels
from matplotlib import pyplot
import numpy as np

def particle_initialise(distr, N, m):
    """
    Function to initialise the particle set
    
    Inputs:
    distr: The type of distribution needing simulation
    N:     The number of particles to be simulated
    m:     The particle mass
    """

    r   = 0.001 | units.parsec

    if distr == 'Triple':
        sys = Particles(N, mass = m)
        sys.position = [[0.3070640329532462, 0.2392405139186416, 0],                    
                        [-0.3865927585370948, -0.10634061110834962, 0],                   
                        [0.27952872558384859, -0.4328999028102919, 0]] | units.AU    

        sys.velocity = [[23.518569826682455, 10.52227851919019, 0],                       
                        [-20.74623746218285, -22.52084408611801, 0],                      
                        [6.227667635500401, 12.998565566927823, 0]] | units.AU/units.yr

       # sys.position += (1, 0, 0) * r
        sys.velocity += (0, 1, 0) * (constants.G*bh_mass/r).sqrt()
        sys.particles = 1 | units.RSun
    
    if distr == 'Plummer':
        np.random.seed(31415)
        mass_stars = m
        converter = nbody_system.nbody_to_si(mass_stars * 100, 0.5 | units.AU)
        sys = new_plummer_model(15, convert_nbody = converter)
        sys.mass = mass_stars
        sys.velocity += (0, 1, 0) * (constants.G*bh_mass/r).sqrt()
        sys.particles = 1 | units.RSun

    return sys

def integrator_initialise(string, grav_int, partic, pri):
    """
    Function to initialise the gravitational integration method.
    
    Inputs:
    string:   String defining the integrator being used.
    grav_int: The class of the gravitational integrator being used
    partic:   The particle set to evolve.
    pri:      The time-step/precision of the simulation
    """
    
    converter = nbody_system.nbody_to_si(1 | units.MSun, 1 | units.AU)
    code = grav_int(converter)
#    stopping_condition = code.stopping_conditions.collision_detection

    if string == 'Brutus 2.5':
        print('Simulating Brutus (PN2.5)')
        code.set_PN_terms(1,1,1,1)
        code.particles.add_particles(partic)
        code.set_eta(dt_param = pri.number)

    elif string == 'Brutus 0.0':
        print('Simulating Brutus (PN0.0)')
        code.set_PN_terms(0,0,0,0)
        code.particles.add_particles(partic)
        code.set_eta(dt_param = pri.number)

    elif string == 'Mikkola':
        print('Simulating Mikkola')
        code.particles.add_particles(partic)
        code.set_time_step(pri)
    
    channel = Channels()
    channel.add_channel(code.particles.new_channel_to(partic))

    return code, channel

def evolve_model(string, brd, pri, part, t_end = 0.01 | units.yr):
    """
    Function to evolve the dynamics of a given distribution.
    
    Inputs:
    string:   String defining the integrator being used.
    brd:      Gravitational code being used
    pri:      The precisionm/time-step parameter
    part:     Particle array.
    t_end:    How long the simulation will evolve till.
    """
    channel = Channels()
    channel.add_channel(brd.particles.new_channel_to(part))

    time = 0 | units.yr
    E0 = brd.particles.potential_energy(G=constants.G)  + brd.particles.kinetic_energy()
    print('Initial State:   ', brd.particles,
          '\n E0:           ', E0)
    while time < t_end:
        time += pri
        brd.evolve_model(time)
        channel.copy()
    
    if string == 'Brutus 2.5':
        Et = brd.particles.potential_energy(G=constants.G)  + brd.particles.kinetic_energy()
        de = abs(Et-E0)/abs(E0)

    elif string == 'Brutus 0.0':
        Et = brd.particles.potential_energy(G=constants.G)  + brd.particles.kinetic_energy()
        de = abs(Et-E0)/abs(E0) 

    elif string == 'Mikkola':
        Et = brd.particles.potential_energy(G=constants.G)  + brd.particles.kinetic_energy()
        Et += brd.get_radiated_gravitational_energy()
        de = abs(Et-E0)/abs(E0)
    brd.stop()

    print('Final State:    ', brd.particles,
        '\nFinal Energy:   ', Et)

    return de

for i in range(3):
    print(i)

############################# GLOBALS
bh_mass = 4.31e6    | units.MSun
bh_pos  = [0, 0, 0] | units.parsec  

precision = 0.01*10.**np.linspace(-5.0, 0, 6) | units.yr
prec_plot = 10.**np.linspace(-5., 0, 6)
N = 3

string_array = ['Mikkola', 'Brutus 2.5', 'Brutus 0.0']
integ_array  = [Mikkola, Brutus, Brutus]

############################# Evolving Various Models


distr_str = 'Plummer'
mass = 1 | units.MSun
dE_Pl = [[ ], [ ], [ ]]
print('Simulating with Plummer distribution')

for i in range(3):
    for j in range(len(precision)):
        particles = particle_initialise(distr_str, N, mass)                                   
        brd, channel = integrator_initialise(string_array[i], integ_array[i], particles, precision[j])
        dE_Pl[i] = evolve_model(string_array[i], brd, precision[j], particles)
    print(dE_Pl)
    pyplot.scatter(precision, dE_Pl[i], label = string_array[i])

precision = 10.**np.linspace(-5., 0, 6)
pyplot.title('Energy Conservation Plot for Plummer Distribution')
pyplot.xlabel(r'Time Step [$\delta t$]')
pyplot.ylabel(r'$\frac{|E_0-E_f|}{|E_0|}$')
pyplot.legend()
pyplot.loglog()
pyplot.savefig('Nbody_Precision_Plummer.pdf', dpi = 300)
pyplot.clf()


distr_str = 'Triple'
mass = 1 | units.MSun
dE_low_array = [ [ ], [ ], [ ] ]
print('Simulating for low mass (Triple)')

for i in range(3):
    for prec in precision:
        print('Simulating for: ', prec.number/0.005)
        particles = particle_initialise(distr_str, N, mass)                                   
        brd, channel = integrator_initialise(string_array[i], integ_array[i], particles, prec)
        dE_low_array[i].append(evolve_model(string_array[i], brd, prec, particles))
    print(dE_low_array)
    print(dE_low_array[i])
    pyplot.scatter(prec_plot, dE_low_array[i], label = string_array[i])

pyplot.title('Energy Conservation Plot for Triple System')
pyplot.xlabel(r'Time Step [$\eta$]')
pyplot.ylabel(r'$\frac{|E_0-E_f|}{|E_0|}$')
pyplot.legend()
pyplot.loglog()
pyplot.savefig('Nbody_Precision_Low.pdf', dpi = 300)
pyplot.clf()

print('Simulating high mass stars')
mass = 100 | units.MSun
dE_array = [ [ ], [ ], [ ] ]

for i in range(3):
    for prec in precision:
        print('Simulating for: ', prec.number/0.005)
        particles = particle_initialise(distr_str, N, mass)                                   
        brd, channel = integrator_initialise(string_array[i], integ_array[i], particles, prec)
        dE_array[i].append(evolve_model(string_array[i], brd, prec, particles))
    print(dE_array)
    print(dE_array[i])
    pyplot.scatter(prec_plot, dE_array[i], label = string_array[i])

pyplot.title('Energy Conservation Plot for Triple System')
pyplot.xlabel(r'Time Step [$\eta$]')
pyplot.ylabel(r'$\frac{|E_0-E_f|}{|E_0|}$')
pyplot.legend()
pyplot.loglog()
pyplot.savefig('Nbody_Precision.pdf', dpi = 300)
pyplot.clf()
