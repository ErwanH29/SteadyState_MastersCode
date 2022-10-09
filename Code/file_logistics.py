from amuse.lab import *
import pickle as pkl
import numpy as np
import glob
import fnmatch
import os

def bulk_stat_extractor(file_string):
    """
    Function which extracts all files in a given dir.
    
    Inputs:
    file_string: The directory wished to extract files from
    """

    filename = glob.glob(file_string)
    filename = sorted(filename, key=lambda t: os.stat(t).st_mtime)
    data = [ ]

    for file_ in range(len(filename)):
        with open(filename[file_], 'rb') as input_file:
            data.append(pkl.load(input_file))
    
    return data

def ejected_extract(set, ejected, col_len):
    """
    Extracts positional info on the ejected particle into an array
    
    Inputs:
    set:     The complete particle set plotting
    ejected: The ejected particle
    col_len: The number of time-steps simulated
    """

    line_x = np.empty((1, col_len, 1))
    line_y = np.empty((1, col_len, 1))
    line_z = np.empty((1, col_len, 1))

    line_vx = np.empty((1, col_len))
    line_vy = np.empty((1, col_len))
    line_vz = np.empty((1, col_len))

    for i in range(len(set)):
        if set.iloc[i,0] == ejected.iloc[0][5]:
            temp_comp = set.iloc[i]
            temp_comp = temp_comp.replace(np.NaN, "[Np.NaN, [np.NaN, np.NaN, np.NaN], [np.NaN, np.NaN, np.NaN]")
            for j in range(col_len):
                coords = temp_comp.iloc[j+1][1]
                vel = temp_comp.iloc[-1][2]  #Only the final velocity is important
                
                if len(coords) == 1:
                    pass
                else:
                    line_x[0][j][0] = coords[0].value_in(units.pc)
                    line_y[0][j][0] = coords[1].value_in(units.pc)
                    line_z[0][j][0] = coords[2].value_in(units.pc)

                    line_vx[0][j] = vel[0].value_in(units.kms)
                    line_vy[0][j] = vel[1].value_in(units.kms)
                    line_vz[0][j] = vel[2].value_in(units.kms)

    return line_x, line_y, line_z, line_vx, line_vy, line_vz

def file_counter():
    """
    Function which counts the number of files in a directory.
    """

    dir_path = r'data/simulation_stats/'
    count = len(fnmatch.filter(os.listdir(dir_path), '*.*'))
    return count

def file_manipulator(col_len, data):
    """
    Function to read POSITIONAL arrays of SINGLE components.
    Manipulates them so they can be read by plotters.
    
    Input:
    col_len: Number of time-steps simulated
    data:    File with the data
    """

    temp_x = np.empty((col_len, 1))
    temp_y = np.empty((col_len, 1))
    temp_z = np.empty((col_len, 1))

    for i in range(col_len):
        temp_vals = data.iloc[i]
        temp_x[i][0] = temp_vals[0].value_in(units.pc)
        temp_y[i][0] = temp_vals[1].value_in(units.pc)
        temp_z[i][0] = temp_vals[2].value_in(units.pc)

    return temp_x, temp_y, temp_z
    
def file_opener(file_string):
    """
    Function which opens and reads the most recent pickle file
    
    Input:
    file_string: The directory for which to access the file.
    """

    filename = glob.glob(file_string)
    with open(os.path.join(max(filename, key=os.path.getctime)), 'rb') as input_file:
        temp_data = pkl.load(input_file)
    temp_length = len(temp_data)

    return temp_data, temp_length

def file_reset(dir):
    """
    Function to remove all files from a directory
    """

    filelist = glob.glob(os.path.join(dir, "*"))
    for f in filelist:
        os.remove(f)

def stats_chaos_extractor(dir):
    steadytime_data = bulk_stat_extractor(dir)
    no_Data = len(steadytime_data)
    
    ini_parti_data = np.empty(no_Data)
    fin_parti_data = np.empty(no_Data)
    number_mergers = np.empty(no_Data)
    cum_merge_mass = np.empty(no_Data)
    simulated_end  = np.empty(no_Data)
    ejected_parti  = np.empty(no_Data)
    stab_time_data = np.empty(no_Data)
    init_dist_data = np.empty(no_Data)
    cluster_radius = np.empty(no_Data)
    init_mass_data = np.empty((no_Data, 2))
    inj_mass_data  = np.empty(no_Data)
    eje_mass_data  = np.empty(no_Data)
    reltime_data   = np.empty(no_Data)

    for i in range(no_Data):
        sim_data = steadytime_data[i]
        ini_parti_data[i] = sim_data.iloc[0][0]
        fin_parti_data[i] = sim_data.iloc[0][1]
        number_mergers[i] = sim_data.iloc[0][2]
        cum_merge_mass[i] = sim_data.iloc[0][3].value_in(units.MSun)
        simulated_end[i]  = sim_data.iloc[0][4].value_in(units.Myr)
        ejected_parti[i]  = sim_data.iloc[0][5]
        stab_time_data[i] = sim_data.iloc[0][6].value_in(units.Myr)
        init_dist_data[i] = sim_data.iloc[0][7].value_in(units.parsec)
        cluster_radius[i] = sim_data.iloc[0][8].value_in(units.parsec)
        init_mass_data[i] = [int(min(sim_data.iloc[0][9].value_in(units.MSun))), int(max(sim_data.iloc[0][9].value_in(units.MSun)))]
        inj_mass_data[i]  = sim_data.iloc[0][10].value_in(units.MSun)
        eje_mass_data[i]  = sim_data.iloc[0][11].value_in(units.MSun)
        #reltime_data[i]   = sim_data.iloc[0][11].value_in(units.yr)

    return ini_parti_data, fin_parti_data, number_mergers, cum_merge_mass, simulated_end, ejected_parti, stab_time_data, \
           init_dist_data, cluster_radius, init_mass_data, inj_mass_data, eje_mass_data, reltime_data

def stats_stable_extractor(dir):
    steadytime_data = bulk_stat_extractor(dir)
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