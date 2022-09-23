from amuse.lab import *
import pickle as pkl
import numpy as np
import glob
import os

def ejected_extract(complete, ejected, col_len):
    """
    Extracts positional info on the ejected particle into an array
    
    Inputs:
    complete: The complete particle set plotting
    ejected:  The ejected particle
    col_len:  The number of time-steps simulated
    """

    line_x = np.empty((1, col_len, 1))
    line_y = np.empty((1, col_len, 1))
    line_z = np.empty((1, col_len, 1))

    for i in range(len(complete)):
        if complete.iloc[i,0] == ejected.iloc[0][-2]:
            temp_comp = complete.iloc[i]
            temp_comp = temp_comp.replace(np.NaN, "[Np.NaN, [np.NaN, np.NaN, np.NaN]")
            for j in range(col_len):
                coords = temp_comp.iloc[j+1][1]
                if len(coords) == 1:
                    pass
                else:
                    line_x[0][j][0] = coords[0].value_in(units.pc)
                    line_y[0][j][0] = coords[1].value_in(units.pc)
                    line_z[0][j][0] = coords[2].value_in(units.pc)
    
    return line_x, line_y, line_z

def file_manipulator(col_len, data):
    """
    Function to read POSITIONAL arrays of SINGLE components.
    Manipulates them so they can be read by plotters.
    
    Input:
    col_len: Number of time-steps simulated
    data:    File with the data
    """

    temp_x = np.empty((1, col_len, 1))
    temp_y = np.empty((1, col_len, 1))
    temp_z = np.empty((1, col_len, 1))

    for i in range(col_len):
        temp_vals = data.iloc[i]
        temp_x[0][i][0] = temp_vals[0].value_in(units.pc)
        temp_y[0][i][0] = temp_vals[1].value_in(units.pc)
        temp_z[0][i][0] = temp_vals[2].value_in(units.pc)

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

def file_reset(directory):
    filelist = glob.glob(os.path.join(directory, "*"))
    for f in filelist:
        os.remove(f)