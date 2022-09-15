import math, sys

from amuse.units import nbody_system, units, constants
from amuse.datamodel import Particles, Particle

from amuse.ext.orbital_elements import new_binary_from_orbital_elements, orbital_elements_from_binary

from amuse.community.brutuspn.interface import Brutus

import matplotlib.pyplot as plt  

from mpmath import mp

#-----------------------------------------------------------------------

def get_elliptic_orbit(m1, m2, a, e, i, O, w, TA):
    p = new_binary_from_orbital_elements(m1, m2, a, e, true_anomaly=TA, inclination=i, longitude_of_the_ascending_node=O, argument_of_periapsis=w, G=constants.G)
    return p

def get_elements(p):
    mass1, mass2, semimajor_axis, eccentricity, true_anomaly, inclination, long_asc_node, arg_per = orbital_elements_from_binary(p, G=constants.G)
    return mass1, mass2, semimajor_axis, eccentricity, inclination, long_asc_node, arg_per, true_anomaly

def GR1(m1, m2, a, e, dt, tend, P):
    t = 0. | units.yr
    omega = []
    w=0
    GM = constants.G*(m1 + m1)
    p = a*(1-e**2)*constants.c**2

    while t < tend:
        omega.append(w)
        w += dt/P*206265*3*GM*math.pi / p
        t += dt
    print("Einstein precession = {} [as]".format(omega[-1]))
    return omega

def sun_and_planets():
    particles = Particles(9)
    sun = particles[0] 
    mercury = particles[1]
    venus = particles[2]
    earth = particles[3]
    mars = particles[4]
    jupiter = particles[5]
    saturn = particles[6]
    uranus = particles[7]
    neptune = particles[8]
        
    sun.mass = 1047.517| units.MJupiter  
    sun.radius = 1.0 | units.RSun      
    sun.position = ( 0.005717 , -0.00538 , -2.130e-5 ) | units.AU
    sun.velocity = ( 0.007893 , 0.01189 , 0.0002064 )  | units.kms 
            
    mercury.mass = 0.000174 | units.MJupiter  
    mercury.radius =  0  | units.RSun    
    mercury.position = ( -0.31419 , 0.14376 , 0.035135 ) | units.AU 
    mercury.velocity = ( -30.729 , -41.93 , -2.659 )  | units.kms  
        
    venus.mass = 0.002564 | units.MJupiter     
    venus.radius =   0    | units.RSun 
    venus.position = ( -0.3767 , 0.60159 , 0.0393 ) | units.AU 
    venus.velocity = ( -29.7725 , -18.849 , 0.795 )  | units.kms 
        
    earth.mass = 0.003185 | units.MJupiter     
    earth.radius =   0    | units.RSun 
    earth.position = ( -0.98561 , 0.0762 , -7.847e-5 ) | units.AU 
    earth.velocity = ( -2.927 , -29.803 , -0.000533 )  | units.kms 
            
    mars.mass = 0.000338 | units.MJupiter     
    mars.radius =     0  | units.RSun 
    mars.position = ( -1.2895 , -0.9199 , -0.048494 ) | units.AU 
    mars.velocity = ( 14.9 , -17.721 , 0.2979 )  | units.kms 
            
    jupiter.mass = 1 | units.MJupiter     
    jupiter.radius =  0    | units.RSun  
    jupiter.position = ( -4.9829 , 2.062 , -0.10990 ) | units.AU 
    jupiter.velocity = ( -5.158 , -11.454 , -0.13558 )  | units.kms 
        
    saturn.mass = 0.29947 | units.MJupiter     
    saturn.radius =   0    | units.RSun 
    saturn.position = ( -2.075 , 8.7812 , 0.3273 ) | units.AU 
    saturn.velocity = ( -9.9109 , -2.236 , -0.2398 )  | units.kms 
        
    uranus.mass = 0.045737 | units.MJupiter     
    uranus.radius =    0   | units.RSun 
    uranus.position = ( -12.0872 , -14.1917 , 0.184214 ) | units.AU 
    uranus.velocity = ( 5.1377 , -4.7387 , -0.06108 )  | units.kms 
            
    neptune.mass = 0.053962 | units.MJupiter 
    neptune.radius =   0    | units.RSun 
    neptune.position = ( 3.1652 , 29.54882 , 0.476391 ) | units.AU 
    neptune.velocity = ( -5.443317 , 0.61054 , -0.144172 )  | units.kms 
        
    particles.move_to_center()
    return particles

#-----------------------------------------------------------------------

if __name__ == "__main__":
    m1 = 1 | units.MSun
    m2 = 3.285e23 | units.kg
    a = 0.38709893 | units.AU
    e = 0.20563069
    i  = 0.0
    O  = 0.0
    w  = 0.0
    TA = 0.0

    binary = new_binary_from_orbital_elements(m1, m2,
                                                   a, e, TA,
                                                   i, O, w,
                                                   G=constants.G)
    solar_syst = sun_and_planets()
    R1 = constants.G*m1/constants.c**2
    R2 = constants.G*m2/constants.c**2
    binary[0].radius = R1
    binary[1].radius = R2

    M_scale = m1+m2
    R_scale = a
    converter = nbody_system.nbody_to_si(M_scale, R_scale)

    grav = Brutus(converter, number_of_workers=4)
    grav.particles.add_particles(binary)
    epsilon = 1e-16
    Lw = 4 * abs(math.log10(epsilon)) + 32
    grav.set_PN_terms(1,0,0,0)
    grav.set_bs_tolerance(1e-16)
    grav.calculate_word_length()
    print(grav.get_bs_tolerance())
    
    mp.prec = Lw
    
    P = 2*math.pi*(a**3/(constants.G*(m1+m2))).sqrt()

    t    = 0. | units.yr
    tend = 100 | units.yr
    dt = 1 | units.yr

    channel = grav.particles.new_channel_to(binary)
        
    ts = []
    omega = []

    try:    
        while t < tend:
            try:
                grav.evolve_model(t)
            except:
                break
            print(t)
            
            channel.copy()
            x_check = binary[1].position[0]
            #print(x_check)
    
            mass1, mass2, semimajor_axis, eccentricity, inclination, long_asc_node, arg_per, true_anomaly = get_elements(binary)
 
            ts.append(t.value_in(units.yr))
            omega.append(3600*arg_per)
            
            t += dt

            sys.stdout.flush()
            
    except KeyboardInterrupt:
        print("terminated.")
        pass

    grav.cleanup_code()
    grav.stop()
    print("Brutus precession = {} [as]".format(omega[-1]))
    omegaE = GR1(m1, m2, a, e, dt, tend, P)
    
    fig, ax = plt.subplots(1)

    ax.plot(ts, omega, 'bo')
    ax.plot(ts, omegaE, 'r-')
    
    ax.set_xlabel("t / [yr]")
    ax.set_ylabel(r"$\omega$ / [as]")

    plt.tight_layout()
        
    plt.savefig('Erwan_Edit.pdf', dpi = 300)
