#include "stdbool.h"

int get_mass(int index_of_the_particle, double * mass);

int commit_particles();

int get_time(double * time);

int set_mass(int index_of_the_particle, double mass);

int get_index_of_first_particle(int * index_of_the_particle);

int new_particle_string(int * identity_of_the_particle, char * mass, char * x, char * y, char * z, char * vx, char * vy, char * vz, char * radius);

int get_total_radius(double * radius);

int get_total_mass_string(char * * M);

int get_bs_tolerance_string(char * * epsilon);

int set_bs_tolerance_string(char * epsilon);

int get_bs_tolerance(double * epsilon);

int get_kinetic_energy_string(char * * kinetic_energy_string);

int get_lightspeed_flag(bool * lightspeed_flag);

int get_total_mass(double * mass);

int get_position_string(int id, char * * x, char * * y, char * * z);

int get_lightspeed(double * zzeta);

int evolve_model(double time);

int set_eta_string(char * dt_param);

int set_word_length(int numBits);

int get_energy(double * energy);

int get_radius_string(int id, char * * radius);

int set_position_string(int id, char * x, char * y, char * z);

int set_eps2(double epsilon_squared);

int set_bs_tolerance(double epsilon);

int set_lightspeed_string(char * zzeta);

int get_begin_time(double * time);

int get_eps2(double * epsilon_squared);

int get_t(double * time);

int get_index_of_next_particle(int index_of_the_particle, int * index_of_the_next_particle);

int calculate_word_length();

int get_word_length(int * numBits);

int delete_particle(int index_of_the_particle);

int get_eta(double * dt_param);

int get_potential(int index_of_the_particle, double * potential);

int synchronize_model();

int set_state(int index_of_the_particle, double mass, double x, double y, double z, double vx, double vy, double vz, double radius);

int get_state(int index_of_the_particle, double * mass, double * x, double * y, double * z, double * vx, double * vy, double * vz, double * radius);

int get_time_step(double * time_step);

int recommit_particles();

int get_kinetic_energy(double * kinetic_energy);

int get_number_of_particles(int * number_of_particles);

int new_particle_float64(int * identity_of_the_particle, double mass, double x, double y, double z, double vx, double vy, double vz, double radius);

int set_eta(double dt_param);

int set_mass_string(int id, char * mass);

int set_acceleration(int index_of_the_particle, double ax, double ay, double az);

int get_center_of_mass_position(double * x, double * y, double * z);

int get_velocity_string(int id, char * * vx, char * * vy, char * * vz);

int set_t(double time);

int set_brutus_output_directory(char * brutus_output_directory);

int get_center_of_mass_velocity(double * vx, double * vy, double * vz);

int set_velocity_string(int id, char * vx, char * vy, char * vz);

int get_radius(int index_of_the_particle, double * radius);

int set_PN_terms(int PN1p, int PN1c, int PN2, int PN2_5);

int set_begin_time(double time);

int set_state_string(int id, char * mass, char * x, char * y, char * z, char * vx, char * vy, char * vz, char * radius);

int get_potential_energy_string(char * * potential_energy_string);

int get_total_energy_string(char * * total_energy_string);

int set_radius(int index_of_the_particle, double radius);

int cleanup_code();

int get_eta_string(char * * dt_param);

int set_lightspeed(double zzeta);

int recommit_parameters();

int get_state_string(int id, char * * m, char * * x, char * * y, char * * z, char * * vx, char * * vy, char * * vz, char * * radius);

int initialize_code();

int get_potential_energy(double * potential_energy);

int get_velocity(int index_of_the_particle, double * vx, double * vy, double * vz);

int get_brutus_output_directory(char * * brutus_output_directory);

int set_t_string(char * time);

int get_position(int index_of_the_particle, double * x, double * y, double * z);

int set_position(int index_of_the_particle, double x, double y, double z);

int set_radius_string(int id, char * radius);

int get_mass_string(int id, char * * m);

int get_t_string(char * * time);

int get_acceleration(int index_of_the_particle, double * ax, double * ay, double * az);

int commit_parameters();

int set_velocity(int index_of_the_particle, double vx, double vy, double vz);

int get_lightspeed_string(char * * zzeta);

