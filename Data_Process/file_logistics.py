from amuse.lab import *
import pickle as pkl
import numpy as np
import glob
import fnmatch
import natsort
import os

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
                data.append(pkl.load(input_file))
                
    else:
        for file_ in range(len(filename)):
            with open(filename[file_], 'rb') as input_file:
                rewrite_file = pkl.load(input_file)
            data_pts = round((np.shape(rewrite_file)[1])/15)
            rewrite_file = rewrite_file.drop(rewrite_file.iloc[:, data_pts:-1*data_pts], axis = 1) 
            data.append(rewrite_file)

    return data

def ejected_stat_extractor(chaos_dir):

    chaos_data = glob.glob(chaos_dir)
    chaos_data = natsort.natsorted(chaos_data)

    filt_IMBH = []
    filt_Chaotic = []

    for file_ in range(len(chaos_data)):
            with open(chaos_data[file_], 'rb') as input_file:
                data = pkl.load(input_file)
                if data.iloc[0][-4] > 0:
                    pass
                else:
                    filt_Chaotic.append(data)
                    input_file = str(input_file)
                    IMBH_data = glob.glob('data/Hermite/particle_trajectory/*'+str(input_file[95:105])+'*')
                    with open(IMBH_data[0], 'rb') as input_file:
                        data = pkl.load(input_file)
                        filt_IMBH.append(data)
    return filt_IMBH, filt_Chaotic

def ejected_extract_final(set, ejected, ejec_merge):
    """
    Extracts the final positional info on the ejected particle into an array
    
    Inputs:
    set:        The complete particle set plotting
    ejected:    The ejected particle
    ejec_merge: String whether we care for ejection simulations (E) or want to account for all simulations (else)
    """

    Nclose = 0
    if ejec_merge == 'E':
        for i in range(len(set)):
            esc_vel = [ ]
            
            if set.iloc[i][0][0] == ejected.iloc[0][4]: 
                ejec_data = set.iloc[i]   #Make data set of only ejected particle
                ejec_data = ejec_data.replace(np.NaN, "[Np.NaN, [np.NaN, np.NaN, np.NaN], [np.NaN, np.NaN, np.NaN]")
                if len(ejec_data) > 3:
                    for vel_ in [ejec_data.iloc[-3][3].value_in(units.kms), 
                                ejec_data.iloc[-2][3].value_in(units.kms), 
                                ejec_data.iloc[-1][3].value_in(units.kms)]: #Last three time steps ~the typical crossing time based on cluster param
                        esc_vel.append(np.sqrt(vel_[0]**2+vel_[1]**2+vel_[2]**2))
                    idx = np.argwhere(esc_vel == max(esc_vel))
                    idx = np.asarray([i-3 for i in idx])[0]
                    esc_vel = esc_vel[idx[0]]

                else:
                    for vel_ in [ejec_data.iloc[-1][3].value_in(units.kms)]: #Last three time steps ~the typical crossing time based on cluster param
                        esc_vel.append(np.sqrt(vel_[0]**2+vel_[1]**2+vel_[2]**2))
                    idx = [-1]

                xpos = (ejec_data.iloc[idx][0][2][0]-set.iloc[0][idx][0][2][0]).value_in(units.pc)
                ypos = (ejec_data.iloc[idx][0][2][1]-set.iloc[0][idx][0][2][1]).value_in(units.pc)
                zpos = (ejec_data.iloc[idx][0][2][2]-set.iloc[0][idx][0][2][2]).value_in(units.pc)

                KE = ejec_data.iloc[-1][4].value_in(units.J)
                PE = ejec_data.iloc[-1][5].value_in(units.J)
                for j in range(len(ejec_data)):
                    if j == 0:
                        pass
                    else:
                        deltaKE = ejec_data.iloc[j][4]/ejec_data.iloc[j-1][4]
                        if abs(deltaKE) > 10:
                            Nclose += 1
                Nmerge = ejected.iloc[0][10]
                            
                return xpos, ypos, zpos, esc_vel, KE, PE, Nclose, Nmerge
        else:
            print('Simulation ended with merger')

    else:
        for i in range(len(set)):
            esc_vel = [ ]
            if set.iloc[i][0][0] == ejected.iloc[0][4]: 
                ejec_data = set.iloc[i]   #Make data set of only ejected particle
                ejec_data = ejec_data.replace(np.NaN, "[Np.NaN, [np.NaN, np.NaN, np.NaN], [np.NaN, np.NaN, np.NaN]")

                if len(ejec_data) > 3:
                    for vel_ in [ejec_data.iloc[-3][3].value_in(units.kms), 
                                ejec_data.iloc[-2][3].value_in(units.kms), 
                                ejec_data.iloc[-1][3].value_in(units.kms)]: #Last three time steps ~the typical crossing time based on cluster param
                        esc_vel.append(np.sqrt(vel_[0]**2+vel_[1]**2+vel_[2]**2))
                    idx = np.argwhere(esc_vel == max(esc_vel))
                    idx = np.asarray([i-3 for i in idx])[0]
                    esc_vel = esc_vel[idx[0]]

                else:
                    for vel_ in [ejec_data.iloc[-1][3].value_in(units.kms)]: #Last three time steps ~the typical crossing time based on cluster param
                        esc_vel.append(np.sqrt(vel_[0]**2+vel_[1]**2+vel_[2]**2))
                    idx = [-1]

                xpos = (ejec_data.iloc[idx][0][2][0]-set.iloc[0][idx][0][2][0]).value_in(units.pc)
                ypos = (ejec_data.iloc[idx][0][2][1]-set.iloc[0][idx][0][2][1]).value_in(units.pc)
                zpos = (ejec_data.iloc[idx][0][2][2]-set.iloc[0][idx][0][2][2]).value_in(units.pc)

                KE = ejec_data.iloc[idx][0][4].value_in(units.J)
                PE = ejec_data.iloc[idx][0][5].value_in(units.J)

                for j in range(len(ejec_data)):
                    if j == 0:
                        pass
                    else:
                        deltaKE = ejec_data.iloc[j][4]/ejec_data.iloc[j-1][4]
                        if abs(deltaKE) > 10:
                            Nclose += 1
                Nmerge = ejected.iloc[0][10]

                return xpos, ypos, zpos, esc_vel, KE, PE, Nclose, Nmerge
    return np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN

def ejected_index(set, ejected):
    """
    Extracts index of the ejected particle
    
    Inputs:
    set:     The complete particle set plotting
    ejected: The ejected particle
    """

    merger = False
    eject = False
    
    for i in range(len(set)): #Loop for particle that is merged
        if isinstance(set.iloc[i][-1][0], float):
            print('Simulation ended in merger')
            merger = True
            ejec_idx = i

    for i in range(len(set)):      
        if not (merger): #Loop for particle ejected but NOT bound
            if set.iloc[i][-1][0] == ejected.iloc[0][4]: 
                print('Simulation ended with ejection')
                eject = True
                ejec_idx = i
    
    if not (eject) and not (merger):
        print('Simulation ended over the time limit')
        for i in range(len(set)):
            ejec_idx = 5       #Replace this with the most sustaining binary

    return ejec_idx

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
            pass
        
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