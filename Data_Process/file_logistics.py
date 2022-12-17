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

def ejected_extract_final(set, ejected):
    """
    Extracts the final info on the ejected particle into an array
    
    Inputs:
    set:        The complete particle set plotting
    ejected:    The ejected particle
    """

    Nclose = 0
    for parti_ in range(len(set)):
        if set.iloc[parti_][0][0] == ejected.iloc[0][4]: 
            ejec_data = set.iloc[parti_]   #Make data set of only ejected particle
            ejec_data = ejec_data.replace(np.NaN, "[Np.NaN, [np.NaN, np.NaN, np.NaN], [np.NaN, np.NaN, np.NaN]")
            ejec_vel = []
            tot_steps = min(round(len(ejec_data)**0.5), 10)
            
            for steps_ in range(tot_steps):
                vel_ = ejec_data.iloc[(-steps_)][3].value_in(units.kms)
                vx = vel_[0] - set.iloc[0][(-steps_)][3][0].value_in(units.kms)
                vy = vel_[1] - set.iloc[0][(-steps_)][3][1].value_in(units.kms)
                vz = vel_[2] - set.iloc[0][(-steps_)][3][2].value_in(units.kms)
                ejec_vel.append(np.sqrt(vx**2+vy**2+vz**2))
                
            idx = np.where(ejec_vel == np.nanmax(ejec_vel))[0]
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
        lbound = 97
        ubound = 107
        
    else:
        lbound = 91
        ubound = 101

    for file_ in range(len(chaos_data)):
        with open(chaos_data[file_], 'rb') as input_file:
            print('Reading file : ', input_file)
            data = pkl.load(input_file)
            if data.iloc[0][-4] < 1:
                filt_Chaotic.append(data)
                input_file = str(input_file)
                IMBH_data = glob.glob('/media/erwanh/Elements/'+str(int)+'/particle_trajectory/*'+str(input_file[lbound:ubound])+'*')
                with open(IMBH_data[0], 'rb') as input_file:
                    data = pkl.load(input_file)
                    data_pts = round((np.shape(data)[1])/15)
                    data = data.drop(data.iloc[:, data_pts:-1*data_pts], axis = 1)
                    filt_IMBH.append(data)

    return filt_IMBH, filt_Chaotic

def ejected_index(pset, ejected):
    """
    Extracts index of the ejected particle
    
    Inputs:
    pset:     The complete particle pset plotting
    ejected:  The ejected particle
    """

    merger = False
    eject = False

    for parti_ in range(len(pset)):
        part_mass = pset.iloc[parti_][0][1]
        SMBH_mass = pset.iloc[0][0][1]
        tot_mass = part_mass + SMBH_mass

        if not (isinstance(pset.iloc[parti_][-1][0], np.uint64)) or pset.iloc[parti_][-1][1] == tot_mass:
            ejec_idx = parti_
            merger = True
            string = 'merger'

            if ejec_idx == 0:
                print('SMBH gets absorbed')
                com_px = 0 | units.MSun * units.kms
                com_py = 0 | units.MSun * units.kms
                com_pz = 0 | units.MSun * units.kms
                mass = 0 | units.MSun

                for j in range(len(pset)):
                    com_px += pset.iloc[j][0][1] * pset.iloc[j][-2][3][0]
                    com_py += pset.iloc[j][0][1] * pset.iloc[j][-2][3][1]
                    com_pz += pset.iloc[j][0][1] * pset.iloc[j][-2][3][2]
                    mass += pset.iloc[j][0][1]
                com_px /= mass
                com_py /= mass
                com_pz /= mass

                velx_SMBH = pset.iloc[0][-2][3][0]
                vely_SMBH = pset.iloc[0][-2][3][1]
                velz_SMBH = pset.iloc[0][-2][3][2]
                vel_SMBH = ((velx_SMBH-com_px)**2 + (vely_SMBH-com_py)**2 + (velz_SMBH-com_pz)**2).sqrt()
                
                for j in range(len(pset)):
                    if pset.iloc[j][-1][1] > 4*10**6 | units.MSun:
                        ejec_idx = parti_

                        IMBH_mass = pset.iloc[j][0][1]
                        velx_IMBH = pset.iloc[j][-2][3][0]
                        vely_IMBH = pset.iloc[j][-2][3][1]
                        velz_IMBH = pset.iloc[j][-2][3][2]
                        vel_IMBH = ((velx_IMBH-com_px)**2 + (vely_IMBH-com_py)**2 + (velz_IMBH-com_pz)**2).sqrt()
                        vel_merger = (IMBH_mass*vel_IMBH + SMBH_mass*vel_SMBH)/(IMBH_mass+SMBH_mass)

                        print('SMBH com velocity:    ', vel_SMBH.in_(units.kms))
                        print('remnant com velocity: ', vel_merger.in_(units.kms))

    if not (merger) and np.shape(pset)[1] <= 10**6: 
        for parti_ in range(len(pset)):   
            if pset.iloc[parti_][-1][0] == ejected.iloc[0][4]:
                ejec_idx = parti_
                eject = True
                string = 'ejected'
                if parti_ == 0:
                    com_px = 0 | units.MSun * units.kms
                    com_py = 0 | units.MSun * units.kms
                    com_pz = 0 | units.MSun * units.kms
                    mass = 0 | units.MSun
                    for j in range(len(pset)):
                        com_px += pset.iloc[j][0][1] * pset.iloc[j][-2][3][0]
                        com_py += pset.iloc[j][0][1] * pset.iloc[j][-2][3][1]
                        com_pz += pset.iloc[j][0][1] * pset.iloc[j][-2][3][2]
                        mass += pset.iloc[j][0][1]
                    com_px /= mass
                    com_py /= mass
                    com_pz /= mass

                    SMBH_mass = pset.iloc[0][0][1]
                    velx_SMBH = pset.iloc[0][-2][3][0]
                    vely_SMBH = pset.iloc[0][-2][3][1]
                    velz_SMBH = pset.iloc[0][-2][3][2]
                    vel_SMBH = ((velx_SMBH-com_px)**2 + (vely_SMBH-com_py)**2 + (velz_SMBH-com_pz)**2).sqrt()
                
                    print('SMBH ejected with velocity: ', vel_SMBH.in_(units.kms))
    
    if not (eject) and not (merger):
        ejec_idx = 5 
        string = 'over time'

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

def stats_chaos_extractor(dir):
    steadytime_data = bulk_stat_extractor(dir, 'N')
    no_Data = len(steadytime_data)

    fin_parti_data = np.empty(no_Data)
    init_dist_data = np.empty(no_Data)
    init_mass_data = np.empty((no_Data, 2))
    stab_time_data = np.empty(no_Data)

    for i in range(no_Data):
        sim_data = steadytime_data[i]
        
        if isinstance(sim_data.iloc[0][9], float):
            fin_parti_data[i] = sim_data.iloc[0][6] + sim_data.iloc[0][10]
            init_dist_data[i] = sim_data.iloc[0][7].value_in(units.parsec)
            init_mass_data[i] = [int(min(sim_data.iloc[0][8].value_in(units.MSun))), 
                                 int(max(sim_data.iloc[0][8].value_in(units.MSun)))]
            stab_time_data[i] = sim_data.iloc[0][13].value_in(units.Myr)
        
    return fin_parti_data, stab_time_data, init_dist_data, init_mass_data

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