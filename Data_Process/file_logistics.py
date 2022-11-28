from amuse.lab import *
import pickle as pkl
import numpy as np
import glob
import fnmatch
import natsort
import os
import matplotlib.ticker as mtick

class plotter_setup(object):
    def moving_average(self, array, smoothing):
        """
        Function to remove the large fluctuations in various properties by taking the average
        
        Inputs:
        array:     Array consisting of variable for which to produce running average of
        smoothing: Number of elements to average over
        output:    Smoothened array
        """

        value = np.cumsum(array, dtype=float)
        value[smoothing:] = value[smoothing:] - value[:-smoothing]

        return value[smoothing-1:]/smoothing
         
    def tickers(self, ax, plot_type):
        """
        Function to setup axis
        """

        ax.yaxis.set_ticks_position('both')
        ax.xaxis.set_ticks_position('both')
        ax.xaxis.set_minor_locator(mtick.AutoMinorLocator())
        ax.yaxis.set_minor_locator(mtick.AutoMinorLocator())
        if plot_type == 'plot':
            ax.tick_params(axis="y", which = 'both', direction="in")
            ax.tick_params(axis="x", which = 'both', direction="in")

        return ax

    def tickers_pop(self, ax, pop):
        """
        Function to setup axis for population plots
        """

        xints = [i for i in range(1+int(max(pop))) if i % 10 == 0]

        ax.set_xlabel(r'IMBH Population [$N$]')
        ax.yaxis.set_ticks_position('both')
        ax.xaxis.set_ticks_position('both')
        ax.yaxis.set_minor_locator(mtick.AutoMinorLocator())
        ax.tick_params(axis="y", which = 'both', direction="in")
        ax.tick_params(axis="x", which = 'both', direction="in")     
        ax.set_xticks(xints)
        ax.set_xlim(min(xints)-5, max(xints)+5)
        ax.set_xlim(5, 105)

        return ax
        
    def val_filter(self, arr):
        """
        Function which removes the excessive terms
        """
        
        arr[:][abs(arr[:]) > 10**7] = np.NaN
        return arr
    
def bulk_stat_extractor(file_string, rewrite):
    """
    Function which extracts all files in a given dir.
    
    Inputs:
    file_string: The directory wished to extract files from
    rewrite:     (Y|N) string to dictate whether to compress the file based on needed data
    """

    filename = glob.glob(file_string)
    filename = natsort.natsorted(filename)
    data = [ ]

    if rewrite == 'N':
        for file_ in range(len(filename)):
            with open(filename[file_], 'rb') as input_file:
                print('Reading file: ', input_file)
                data.append(pkl.load(input_file))
                
    else:
        for file_ in range(len(filename)):
            with open(filename[file_], 'rb') as input_file:
                print('Reading file: ', input_file)
                rewrite_file = pkl.load(input_file)
            if np.shape(rewrite_file)[1] > 50 and np.shape(rewrite_file)[1] < 50:
                data_pts = round((np.shape(rewrite_file)[1])/25)
            elif np.shape(rewrite_file)[1] > 500 and np.shape(rewrite_file)[1] < 5000:
                data_pts = round((np.shape(rewrite_file)[1])/250)
            elif np.shape(rewrite_file)[1] > 5000:
                data_pts = round((np.shape(rewrite_file)[1])/2500)
            else:
                data_pts = np.shape(rewrite_file)[1]
            rewrite_file = rewrite_file.drop(rewrite_file.iloc[:, data_pts:-1*data_pts], axis = 1) 
            data.append(rewrite_file)

    return data

def ejected_extract_final(set, ejected, ejec_merge):
    """
    Extracts the final info on the ejected particle into an array
    
    Inputs:
    set:        The complete particle set plotting
    ejected:    The ejected particle
    ejec_merge: String whether we care for ejection simulations (E) or want to account for all simulations (else)
    """

    Nclose = 0
    if ejec_merge == 'E':
        for i in range(len(set)):
            if set.iloc[i][0][0] == ejected.iloc[0][4]: 
                print('Ejection detected for pop.: ', len(set)-1)
                ejec_data = set.iloc[i]   #Make data set of only ejected particle
                ejec_data = ejec_data.replace(np.NaN, "[Np.NaN, [np.NaN, np.NaN, np.NaN], [np.NaN, np.NaN, np.NaN]")
                ejec_vel = []
                tot_steps = round(len(ejec_data)/2)
                
                for steps_ in range(tot_steps):
                    vel_ = ejec_data.iloc[(-steps_)][3].value_in(units.kms)
                    vx = vel_[0] - set.iloc[0][(-steps_)][3][0].value_in(units.kms)
                    vy = vel_[1] - set.iloc[0][(-steps_)][3][1].value_in(units.kms)
                    vz = vel_[2] - set.iloc[0][(-steps_)][3][2].value_in(units.kms)
                    ejec_vel.append(np.sqrt(vx**2+vy**2+vz**2))
                idx = np.where(ejec_vel == max(ejec_vel))[0]
                idx -= tot_steps
                ejec_vel = np.asarray(ejec_vel)
                esc_vel = ejec_vel[idx]

                xpos = (ejec_data.iloc[idx][0][2][0]-set.iloc[0][idx][0][2][0]).value_in(units.pc)
                ypos = (ejec_data.iloc[idx][0][2][1]-set.iloc[0][idx][0][2][1]).value_in(units.pc)
                zpos = (ejec_data.iloc[idx][0][2][2]-set.iloc[0][idx][0][2][2]).value_in(units.pc)

                KE = ejec_data.iloc[-1][4].value_in(units.J)
                PE = ejec_data.iloc[-1][5].value_in(units.J)

                for j in range(len(ejec_data)):
                    if abs(ejec_data.iloc[j][-1]) < 1e-2:
                        Nclose += 1
                Nmerge = ejected.iloc[0][10]
                            
                return xpos, ypos, zpos, esc_vel, KE, PE, Nclose, Nmerge
        else:
            print('Simulation ended with merger')

    else:
        for i in range(len(set)):
            if set.iloc[i][0][0] == ejected.iloc[0][4]: 
                ejec_data = set.iloc[i]   #Make data set of only ejected particle
                ejec_data = ejec_data.replace(np.NaN, "[Np.NaN, [np.NaN, np.NaN, np.NaN], [np.NaN, np.NaN, np.NaN]")
                ejec_vel = []
                tot_steps = round(len(ejec_data)/2)
                
                for steps_ in range(tot_steps):
                    vel_ = ejec_data.iloc[(-steps_)][3].value_in(units.kms)
                    vx = vel_[0] - set.iloc[0][(-steps_)][3][0].value_in(units.kms)
                    vy = vel_[1] - set.iloc[0][(-steps_)][3][1].value_in(units.kms)
                    vz = vel_[2] - set.iloc[0][(-steps_)][3][2].value_in(units.kms)
                    ejec_vel.append(np.sqrt(vx**2+vy**2+vz**2))
                idx = np.argwhere(ejec_vel == max(ejec_vel))[0]
                idx -= tot_steps
                ejec_vel = np.asarray(ejec_vel)
                esc_vel = ejec_vel[idx]

                xpos = (ejec_data.iloc[idx][0][2][0]-set.iloc[0][idx][0][2][0]).value_in(units.pc)
                ypos = (ejec_data.iloc[idx][0][2][1]-set.iloc[0][idx][0][2][1]).value_in(units.pc)
                zpos = (ejec_data.iloc[idx][0][2][2]-set.iloc[0][idx][0][2][2]).value_in(units.pc)

                KE = ejec_data.iloc[idx][0][4].value_in(units.J)
                PE = ejec_data.iloc[idx][0][5].value_in(units.J)

                for j in range(len(ejec_data)):
                    if abs(ejec_data.iloc[j][-1]) < 1e-2:
                        Nclose += 1
                Nmerge = ejected.iloc[0][10]

                return xpos, ypos, zpos, esc_vel, KE, PE, Nclose, Nmerge
    return np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN

def ejected_stat_extractor(chaos_dir, int):
    """
    Function to extract a cropped data set of any given simulation to analyse ejected particles

    Inputs:
    chaos_dir:  The directory to extract data from
    int:        String dictating whether it is Hermite or GRX data
    """

    chaos_data = glob.glob(chaos_dir)
    chaos_data = natsort.natsorted(chaos_data)

    filt_IMBH = []
    filt_Chaotic = []

    if int == 'Hermite':
        lbound = 95
        ubound = 105
    else:
        lbound = 91
        ubound = 101

    for file_ in range(len(chaos_data)):
            with open(chaos_data[file_], 'rb') as input_file:
                data = pkl.load(input_file)
                if data.iloc[0][-4] > 0:
                    pass
                else:
                    filt_Chaotic.append(data)
                    input_file = str(input_file)
                    IMBH_data = glob.glob('data/'+str(int)+'/particle_trajectory/*'+str(input_file[lbound:ubound])+'*')
                    with open(IMBH_data[0], 'rb') as input_file:
                        data = pkl.load(input_file)
                        data_pts = round((np.shape(data)[1])/15)
                        data = data.drop(data.iloc[:, data_pts:-1*data_pts], axis = 1)
                        filt_IMBH.append(data)

    return filt_IMBH, filt_Chaotic

def ejected_index(pset, ejected, int):
    """
    Extracts index of the ejected particle
    
    Inputs:
    pset:     The complete particle pset plotting
    ejected:  The ejected particle
    int:      The integrator used
    """

    merger = False
    eject = False
    end = False

    if int == 'Hermite':
        for i in range(len(pset)):
            part_mass = pset.iloc[i][0][1]
            SMBH_mass = pset.iloc[0][0][1]
            tot_mass = part_mass + SMBH_mass
            if not (isinstance(pset.iloc[i][-1][0], np.uint64)) or pset.iloc[i][-1][1] == tot_mass:
                ejec_idx = i
                merger = True
                string = 'merger'

                if ejec_idx == 0:
                    print('SMBH gets absorbed')

                    com_velx = 0 | units.MSun * units.kms
                    com_vely = 0 | units.MSun * units.kms
                    com_velz = 0 | units.MSun * units.kms
                    mass = 0 | units.MSun
                    for j in range(len(pset)):
                        com_velx += pset.iloc[j][0][1] * pset.iloc[j][-2][3][0]
                        com_vely += pset.iloc[j][0][1] * pset.iloc[j][-2][3][1]
                        com_velz += pset.iloc[j][0][1] * pset.iloc[j][-2][3][2]
                        mass += pset.iloc[j][0][1]
                    com_velx /= mass
                    com_vely /= mass
                    com_velz /= mass

                    SMBH_mass = pset.iloc[0][0][1]
                    velx_SMBH = pset.iloc[0][-2][3][0]
                    vely_SMBH = pset.iloc[0][-2][3][1]
                    velz_SMBH = pset.iloc[0][-2][3][2]
                    vel_SMBH = ((velx_SMBH-com_velx)**2 + (vely_SMBH-com_vely)**2 + (velz_SMBH-com_velz)**2).sqrt()
                    
                    for i in range(len(pset)):
                        if pset.iloc[i][-1][1] > 4*10**6 | units.MSun:
                            ejec_idx = i

                            IMBH_mass = pset.iloc[i][0][1]
                            velx_IMBH = pset.iloc[i][-2][3][0]
                            vely_IMBH = pset.iloc[i][-2][3][1]
                            velz_IMBH = pset.iloc[i][-2][3][2]
                            vel_IMBH = ((velx_IMBH-com_velx)**2 + (vely_IMBH-com_vely)**2 + (velz_IMBH-com_velz)**2**2).sqrt()
                            vel_merger = (IMBH_mass*vel_IMBH + SMBH_mass*vel_SMBH)/(IMBH_mass+SMBH_mass)

                            print('SMBH com velocity:    ', vel_SMBH.in_(units.kms))
                            print('remnant com velocity: ', vel_merger.in_(units.kms))

        if not (merger) and np.shape(pset)[1] < 100000: 
            for i in range(len(pset)):   
                if pset.iloc[i][-1][0] == ejected.iloc[0][4]:
                    ejec_idx = i
                    eject = True
                    string = 'ejected'
                    if i == 0:
                        com_velx = 0 | units.MSun * units.kms
                        com_vely = 0 | units.MSun * units.kms
                        com_velz = 0 | units.MSun * units.kms
                        mass = 0 | units.MSun
                        for j in range(len(pset)):
                            com_velx += pset.iloc[j][0][1] * pset.iloc[j][-2][3][0]
                            com_vely += pset.iloc[j][0][1] * pset.iloc[j][-2][3][1]
                            com_velz += pset.iloc[j][0][1] * pset.iloc[j][-2][3][2]
                            mass += pset.iloc[j][0][1]
                        com_velx /= mass
                        com_vely /= mass
                        com_velz /= mass

                        SMBH_mass = pset.iloc[0][0][1]
                        velx_SMBH = pset.iloc[0][-2][3][0]
                        vely_SMBH = pset.iloc[0][-2][3][1]
                        velz_SMBH = pset.iloc[0][-2][3][2]
                        vel_SMBH = ((velx_SMBH-com_velx)**2 + (vely_SMBH-com_vely)**2 + (velz_SMBH-com_velz)**2).sqrt()
                    
                        print('SMBH ejected with velocity: ', vel_SMBH.in_(units.kms))
        
        if not (eject) and not (merger):
            ejec_idx = 5       #Replace this with the most sustaining binary
            string = 'over time'

    else:
        if np.shape(pset)[1] > 100000: 
            ejec_idx = 5 
            string = 'over time'
            end = True

        if not (end):
            for i in range(len(pset)):
                part_mass = pset.iloc[i][0][1]
                SMBH_mass = pset.iloc[0][0][1]
                tot_mass = part_mass + SMBH_mass
                if not (isinstance(pset.iloc[i][-1][0], np.uint64)) or pset.iloc[i][-1][1] == tot_mass:
                    ejec_idx = i
                    merger = True
                    string = 'merger'

                    if ejec_idx == 0:
                        print('SMBH gets absorbed')

                        com_velx = 0 | units.MSun * units.kms
                        com_vely = 0 | units.MSun * units.kms
                        com_velz = 0 | units.MSun * units.kms
                        mass = 0 | units.MSun
                        for j in range(len(pset)):
                            com_velx += pset.iloc[j][0][1] * pset.iloc[j][-2][3][0]
                            com_vely += pset.iloc[j][0][1] * pset.iloc[j][-2][3][1]
                            com_velz += pset.iloc[j][0][1] * pset.iloc[j][-2][3][2]
                            mass += pset.iloc[j][0][1]
                        com_velx /= mass
                        com_vely /= mass
                        com_velz /= mass

                        SMBH_mass = pset.iloc[0][0][1]
                        velx_SMBH = pset.iloc[0][-2][3][0]
                        vely_SMBH = pset.iloc[0][-2][3][1]
                        velz_SMBH = pset.iloc[0][-2][3][2]
                        vel_SMBH = ((velx_SMBH-com_velx)**2 + (vely_SMBH-com_vely)**2 + (velz_SMBH-com_velz)**2).sqrt()
                    
                        for i in range(len(pset)):
                            if pset.iloc[i][-1][1] > 4*10**6 | units.MSun:
                                ejec_idx = i

                                IMBH_mass = pset.iloc[i][0][1]
                                velx_IMBH = pset.iloc[i][-2][3][0]
                                vely_IMBH = pset.iloc[i][-2][3][1]
                                velz_IMBH = pset.iloc[i][-2][3][2]
                                vel_IMBH = ((velx_IMBH-com_velx)**2 + (vely_IMBH-com_vely)**2 + (velz_IMBH-com_velz)**2).sqrt()
                                vel_merger = (IMBH_mass*vel_IMBH + SMBH_mass*vel_SMBH)/(IMBH_mass+SMBH_mass)

                                print('SMBH com velocity:    ', vel_SMBH.in_(units.kms))
                                print('remnant com velocity: ', vel_merger.in_(units.kms))

        if not (merger) and not (end): 
            for i in range(len(pset)):
                if pset.iloc[i][-1][0] == ejected.iloc[0][4]:
                    ejec_idx = i
                    eject = True
                    string = 'ejected'
                    if i == 0:
                        vel = [0, 0, 0] | units.MSun * units.kms
                        mass = 0 | units.MSun
                        for j in range(len(pset)):
                            vel += pset.iloc[j][0][1] * pset.iloc[j][-2]
                            mass += pset.iloc[j][0][1]
                        vel /= mass
                        
                        velx = pset.iloc[i][-1][3][0]
                        vely = pset.iloc[i][-1][3][1]
                        velz = pset.iloc[i][-1][3][2]
                        vel = ((velx-vel[0])**2 + (vely-vel[1])**2 + (velz-vel[2])**2).sqrt().in_(units.kms)
                        print('SMBH ejected with velocity: ', vel)

    print(string)
    return ejec_idx, string


def file_counter(int_string):
    """
    Function which counts the number of files in a directory.
    """
    
    dir_path = 'data/'+str(int_string)+'/simulation_stats/'
    return len(fnmatch.filter(os.listdir(dir_path), '*.*'))
    
def file_opener(dir):
    """
    Function which opens and reads the most recent pickle file
    
    Input:
    dir: The directory for which to access the file.
    """

    filename = glob.glob(dir)
    with open(os.path.join(max(filename, key=os.path.getctime)), 'rb') as input_file:
        temp_data = pkl.load(input_file)

    return temp_data

def file_reset(dir):
    """
    Function to remove all files from a directory
    """

    filelist = glob.glob(os.path.join(dir, "*"))
    for f in filelist:
        os.remove(f)

def stats_chaos_extractor(dir):
    steadytime_data = bulk_stat_extractor(dir, 'N')
    no_Data = len(steadytime_data)

    inj_mass_data  = np.empty(no_Data)
    cum_merge_mass = np.empty(no_Data)
    eje_mass_data  = np.empty(no_Data)
    ejected_parti  = np.empty(no_Data)
    fin_parti_data = np.empty(no_Data)
    init_dist_data = np.empty(no_Data)
    init_mass_data = np.empty((no_Data, 2))
    ini_parti_data = np.empty(no_Data)
    number_mergers = np.empty(no_Data)
    simulated_end  = np.empty(no_Data)
    stab_time_data = np.empty(no_Data)

    faults = 0
    for i in range(no_Data):
        sim_data = steadytime_data[i]
        
        if isinstance(sim_data.iloc[0][9], float):
            inj_mass_data[i]  = sim_data.iloc[0][0].value_in(units.MSun)
            cum_merge_mass[i] = sim_data.iloc[0][2].value_in(units.MSun)
            eje_mass_data[i]  = sim_data.iloc[0][3].value_in(units.MSun)
            ejected_parti[i]  = sim_data.iloc[0][4]
            fin_parti_data[i] = sim_data.iloc[0][6] + sim_data.iloc[0][10]
            init_dist_data[i] = sim_data.iloc[0][7].value_in(units.parsec)
            init_mass_data[i] = [int(min(sim_data.iloc[0][8].value_in(units.MSun))), int(max(sim_data.iloc[0][8].value_in(units.MSun)))]
            ini_parti_data[i] = sim_data.iloc[0][9]
            number_mergers[i] = sim_data.iloc[0][10]
            simulated_end[i]  = sim_data.iloc[0][12].value_in(units.Myr)
            stab_time_data[i] = sim_data.iloc[0][13].value_in(units.Myr)
        else:
            faults += 1
        
    return ini_parti_data, fin_parti_data, number_mergers, cum_merge_mass, simulated_end, ejected_parti, \
           stab_time_data, init_dist_data, init_mass_data, inj_mass_data, eje_mass_data

def stats_stable_extractor(dir):
    steadytime_data = bulk_stat_extractor(dir, 'N')
    no_Data = len(steadytime_data)
    
    ini_parti_data = np.empty(no_Data)
    inj_event_data = np.empty(no_Data)
    merge_no_data  = np.empty(no_Data)
    mergerm_data   = np.empty(no_Data)
    simulated_end  = np.empty(no_Data)
    initial_dist   = np.empty(no_Data)
    cluster_rad    = np.empty(no_Data)
    init_parti_m   = np.empty((no_Data, 2))
    trelax_data    = np.empty(no_Data)

    for i in range(no_Data):
        sim_data = steadytime_data[i]
        ini_parti_data[i] = sim_data.iloc[0][0]
        inj_event_data[i] = sim_data.iloc[0][1]
        merge_no_data[i]  = sim_data.iloc[0][2]
        mergerm_data[i]   = sim_data.iloc[0][3].value_in(units.MSun)
        simulated_end[i]  = sim_data.iloc[0][4].value_in(units.Myr)
        initial_dist[i]   = sim_data.iloc[0][5].value_in(units.pc)
        cluster_rad[i]    = sim_data.iloc[0][6].value_in(units.pc)
        init_parti_m[i]   = [int(min(sim_data.iloc[0][7].value_in(units.MSun))), int(max(sim_data.iloc[0][7].value_in(units.MSun)))]
        trelax_data[i]    = sim_data.iloc[0][8].value_in(units.Myr)

    return ini_parti_data, inj_event_data, merge_no_data, mergerm_data, simulated_end, initial_dist, \
           cluster_rad, init_parti_m, trelax_data

def simulation_stats_checker(int_string):
    """
    Function to check the final outcomes of all the simulations
    """

    filename = glob.glob('data/'+str(int_string)+'/simulation_stats/*')
    filename = natsort.natsorted(filename)
    SMBH_merger = 0
    IMBH_merger = 0
    ejection = 0
    tot_sims = 0
    complete = 0

    for file_ in range(len(filename)):
        tot_sims += 1
        with open(filename[file_]) as f:
            line = f.readlines()
            line1 = line[-9][:-7]
            line2 = line[3][17:]
            data = line1.split()
            data2 = line2.split()
            data = [float(i) for i in data]
            if 4001000 in data:
                SMBH_merger += 1
            if 2000 in data:
                IMBH_merger += 1
            if 4001000 not in data and 2000 not in data and '100000000.0' not in data2:
                ejection += 1
            if '100000000.0' in data2:
                complete += 1

    with open('figures/'+str(int_string)+'summary.txt', 'w') as file:
        file.write('\nSimulation outcomes for '+str(int_string))
        file.write('\nTotal simulations:   '+str(tot_sims))
        file.write('\nSMBH merging events: '+str(SMBH_merger))
        file.write('\nIMBH merging events: '+str(IMBH_merger))
        file.write('\nEjection events:     '+str(ejection))
        file.write('\nCompleted sims:      '+str(complete))
        file.write('\n========================================')

simulation_stats_checker('GRX')
simulation_stats_checker('Hermite')