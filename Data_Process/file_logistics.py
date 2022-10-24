from amuse.lab import *
import pickle as pkl
import numpy as np
import glob
import fnmatch
import os

def bulk_stat_extractor(file_string, rewrite):
    """
    Function which extracts all files in a given dir.
    
    Inputs:
    file_string: The directory wished to extract files from
    rewrite:     (Y|N) string to dictate whether to compress the file based on needed data
    """

    filename = glob.glob(file_string)
    filename = sorted(filename, key=lambda t: os.stat(t).st_mtime)
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
            rewrite_file = rewrite_file.drop(rewrite_file.iloc[:, :-1*data_pts], axis = 1) 
            data.append(rewrite_file)

    return data

def ejected_extract_traj(set, ejected, col_len):
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
        if set.iloc[i][0][0] == ejected.iloc[0][5]: 
            ejec_data = set.iloc[i]
            ejec_data = ejec_data.replace(np.NaN, "[Np.NaN, [np.NaN, np.NaN, np.NaN], [np.NaN, np.NaN, np.NaN]")
            for j in range(col_len):
                coords = ejec_data.iloc[j+1][2]
                vel = ejec_data.iloc[-1][3]  #Only the final velocity is important

                line_x[0][j][0] = coords[0].value_in(units.pc)
                line_y[0][j][0] = coords[1].value_in(units.pc)
                line_z[0][j][0] = coords[2].value_in(units.pc)

                line_vx[0][j] = vel[0].value_in(units.kms)
                line_vy[0][j] = vel[1].value_in(units.kms)
                line_vz[0][j] = vel[2].value_in(units.kms)

    return line_x, line_y, line_z, line_vx, line_vy, line_vz

def ejected_extract_final(set, ejected):
    """
    Extracts positional info on the ejected particle into an array
    
    Inputs:
    set:     The complete particle set plotting
    ejected: The ejected particle
    col_len: The number of time-steps simulated
    """

    Nclose = 0
    for i in range(len(set)):
        esc_vel = [ ]
        if set.iloc[i][0][0] == ejected.iloc[0][5]: 
            ejec_data = set.iloc[i]   #Make data set of only ejected particle
            ejec_data = ejec_data.replace(np.NaN, "[Np.NaN, [np.NaN, np.NaN, np.NaN], [np.NaN, np.NaN, np.NaN]")
            if len(ejec_data) > 3:
                for vel_ in [ejec_data.iloc[-3][3].value_in(units.kms), 
                             ejec_data.iloc[-2][3].value_in(units.kms), 
                             ejec_data.iloc[-1][3].value_in(units.kms)]: #Last three time steps ~the typical crossing time based on cluster param
                    esc_vel.append(np.sqrt(vel_[0]**2+vel_[1]**2+vel_[2]**2))
                idx = np.argwhere(esc_vel == max(esc_vel))
                idx = np.asarray([i-3 for i in idx])[0]

            else:
                for vel_ in [ejec_data.iloc[-1][3].value_in(units.kms)]: #Last three time steps ~the typical crossing time based on cluster param
                    esc_vel.append(np.sqrt(vel_[0]**2+vel_[1]**2+vel_[2]**2))
                idx = [-1]

            xpos = (ejec_data.iloc[idx][0][2][0]-set.iloc[0][idx][0][2][0]).value_in(units.pc)
            ypos = (ejec_data.iloc[idx][0][2][1]-set.iloc[0][idx][0][2][1]).value_in(units.pc)
            zpos = (ejec_data.iloc[idx][0][2][2]-set.iloc[0][idx][0][2][2]).value_in(units.pc)
            distance = np.sqrt(xpos**2+ypos**2+zpos**2)
            
            vx = ejec_data.iloc[idx][0][3][0].value_in(units.kms)
            vy = ejec_data.iloc[idx][0][3][1].value_in(units.kms)
            vz = ejec_data.iloc[idx][0][3][2].value_in(units.kms)

            incl = np.degrees(np.arcsin(zpos/distance))
            KE = ejec_data.iloc[idx][0][4].value_in(units.J)
            PE = ejec_data.iloc[idx][0][5].value_in(units.J)
            for j in range(len(ejec_data)):
                if j == 0:
                    pass
                else:
                    deltaKE = ejec_data.iloc[j][4]/ejec_data.iloc[j-1][4]
                    if abs(deltaKE) > 10:
                        Nclose += 1

            return xpos, ypos, zpos, vx, vy, vz, KE, PE, incl, Nclose
    
    return print('Nope')

def file_counter(int_string):
    """
    Function which counts the number of files in a directory.
    """
    
    dir_path = r'data/Hermite/GC/simulation_stats/' #Hard-coded change for HErmtieGRX -> GRX
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
    steadytime_data = bulk_stat_extractor(dir, 'N')
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

    faults = 0
    for i in range(no_Data):
        sim_data = steadytime_data[i]
        if isinstance(sim_data.iloc[0][9], float):
            ini_parti_data[i] = sim_data.iloc[0][9]
            fin_parti_data[i] = sim_data.iloc[0][6]
            number_mergers[i] = sim_data.iloc[0][10]
            cum_merge_mass[i] = sim_data.iloc[0][3].value_in(units.MSun)
            simulated_end[i]  = sim_data.iloc[0][-2].value_in(units.Myr)
            ejected_parti[i]  = sim_data.iloc[0][5]
            stab_time_data[i] = sim_data.iloc[0][-1].value_in(units.Myr)
            init_dist_data[i] = sim_data.iloc[0][7].value_in(units.parsec)
            cluster_radius[i] = sim_data.iloc[0][1].value_in(units.parsec)
            init_mass_data[i] = [int(min(sim_data.iloc[0][8].value_in(units.MSun))), int(max(sim_data.iloc[0][8].value_in(units.MSun)))]
            inj_mass_data[i]  = sim_data.iloc[0][0].value_in(units.MSun)
            eje_mass_data[i]  = sim_data.iloc[0][4].value_in(units.MSun)
            #reltime_data[i]   = sim_data.iloc[0][11].value_in(units.yr)
        else:
            faults += 1
            pass
            
    return ini_parti_data, fin_parti_data, number_mergers, cum_merge_mass, simulated_end, ejected_parti, stab_time_data, \
           init_dist_data, cluster_radius, init_mass_data, inj_mass_data, eje_mass_data, reltime_data

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