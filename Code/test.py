import numpy
from amuse.lab import *
from amuse.units import quantities
from amuse.ext.rotating_bridge import Rotating_Bridge
from amuse.community.galaxia.interface import BarAndSpirals3D
from amuse.ext.composition_methods import *

class IntegrateOrbit(object):
    """
    This class makes the integration of the Sun in the Milky Way
    by using BarAndSpirals3D. 
    galaxy(): Function that sets the desired Galactic model. Any question on the parameters, contact me
    creation_particles_noinertial(): creates a parti le set in a rotating frame
    noinertial_to_inertial(): converts data from rotating to inertial frame
    get_pos_vel_and_orbit(): Makes the evolution of the particle set
    """
    
    def __init__(self, t_end= 10 |units.Myr, dt_bridge=0.5 |units.Myr, method= SPLIT_6TH_SS_M13, phase_bar= 0, phase_spiral= 0, omega_spiral= -20 |(units.kms/units.kpc), amplitude= 650|(units.kms**2/units.kpc), m=4, omega_bar= -50 |(units.kms/units.kpc), mass_bar= 1.1e10 |units.MSun ):
        # Simulation parameters
        self.t_end= t_end
        self.dt_bridge= dt_bridge
        self.method= method
        self.time= 0 |units.Myr
        #galaxy parameters
        self.omega= 0 | (units.kms/units.kpc)
        self.initial_phase= 0
        self.bar_phase= phase_bar
        self.spiral_phase= phase_spiral
        self.omega_spiral= omega_spiral
        self.amplitude= amplitude
        self.rsp= 3.12 |units.kpc
        self.m= m
        self.tan_pitch_angle= 0.227194425
        self.omega_bar= omega_bar
        self.mass_bar= mass_bar
        self.aaxis_bar= 3.12 |units.kpc
        self.axis_ratio_bar= 0.37
        return
    
    def galaxy(self):
        global I
        galaxy= BarAndSpirals3D(redirection='file', redirect_stdout_file="GAL{0}.log".format(I))
        I = I + 1
        galaxy.kinetic_energy=quantities.zero
        galaxy.potential_energy=quantities.zero
        galaxy.parameters.bar_contribution= True
        galaxy.parameters.bar_phase= self.bar_phase
        galaxy.parameters.omega_bar= self.omega_bar
        galaxy.parameters.mass_bar= self.mass_bar
        galaxy.parameters.aaxis_bar= self.aaxis_bar
        galaxy.parameters.axis_ratio_bar= self.axis_ratio_bar 
        galaxy.parameters.spiral_contribution= False
        galaxy.parameters.spiral_phase= self.spiral_phase
        galaxy.parameters.omega_spiral= self.omega_spiral
        galaxy.parameters.amplitude= self.amplitude
        galaxy.parameters.rsp= self.rsp
        galaxy.parameters.m= self.m
        galaxy.parameters.tan_pitch_angle= self.tan_pitch_angle
        galaxy.commit_parameters()
        self.omega= galaxy.parameters.omega_system
        self.initial_phase= galaxy.parameters.initial_phase
        print("INITIAL_PHASE:", self.initial_phase)

        galaxy.kinetic_energy=quantities.zero
        galaxy.potential_energy=quantities.zero
        return galaxy 
        
    def creation_particles_noinertial(self, particles):
        """
        makes trans in a counterclockwise frame.
        If the Galaxy only has bar or only spiral arms, the frame corotates with
        the bar or with the spiral arms. If the Galaxy has bar and spiral arms, the frame corotates with the bar
        """
        no_inertial_system= particles.copy()
        angle= self.initial_phase + self.omega*self.time
        C1= particles.vx + self.omega*particles.y
        C2= particles.vy - self.omega*particles.x
        no_inertial_system.x = particles.x*numpy.cos(angle) + particles.y*numpy.sin(angle)
        no_inertial_system.y = -particles.x*numpy.sin(angle) + particles.y*numpy.cos(angle) 
        no_inertial_system.z = particles.z
        no_inertial_system.vx = C1*numpy.cos(angle) + C2*numpy.sin(angle) 
        no_inertial_system.vy = C2*numpy.cos(angle) - C1*numpy.sin(angle)
        no_inertial_system.vz = particles.vz
        return no_inertial_system    

    def noinertial_to_inertial(self, part_noin, part_in):
        #makes trans in a counterclockwise frame
        angle= self.initial_phase + self.omega*self.time
        C1= part_noin.vx - part_noin.y*self.omega
        C2= part_noin.vy + part_noin.x*self.omega
        part_in.x= part_noin.x*numpy.cos(angle)-part_noin.y*numpy.sin(angle)
        part_in.y= part_noin.x*numpy.sin(angle)+part_noin.y*numpy.cos(angle)
        part_in.z= part_noin.z
        part_in.vx= C1*numpy.cos(angle) - C2*numpy.sin(angle)
        part_in.vy= C1*numpy.sin(angle) + C2*numpy.cos(angle)
        part_in.vz= part_noin.vz
        return

    
    def testing_potential_and_force(self, galaxy, x, y, z):
        dx, dy, dz = 0.001 |units.kpc, 0.001 |units.kpc, 0.001 |units.kpc
        phi1x= galaxy.get_potential_at_point(0 |units.kpc, (x+dx), y, z)
        phi2x= galaxy.get_potential_at_point(0 |units.kpc, (x-dx), y, z)
        f1x= -(phi1x-phi2x)/(2*dx)
        phi1y= galaxy.get_potential_at_point(0 |units.kpc, x, (y+dy), z)
        phi2y= galaxy.get_potential_at_point(0 |units.kpc, x, (y-dy), z)
        f1y= -(phi1y-phi2y)/(2*dy)
        phi1z= galaxy.get_potential_at_point(0 |units.kpc, x, y, (z+dz))
        phi2z= galaxy.get_potential_at_point(0 |units.kpc, x, y, (z-dz))
        f1z= -(phi1z-phi2z)/(2*dz)
        fx,fy,fz= galaxy.get_gravity_at_point(0 |units.kpc, x, y, z)
        print("analytic", "numerical") 
        print(fx.value_in(100*units.kms**2/units.kpc) , f1x.value_in(100*units.kms**2/units.kpc))
        print(fy.value_in(100*units.kms**2/units.kpc) , f1y.value_in(100*units.kms**2/units.kpc))
        print(fz.value_in(100*units.kms**2/units.kpc) , f1z.value_in(100*units.kms**2/units.kpc))
        return