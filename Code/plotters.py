from amuse.lab import *
from initialiser import *
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import matplotlib.animation as animation
from distinct_colours import get_distinct
import pickle as pkl
import numpy as np
import glob
import os
from datetime import datetime

def file_opener(file_string):
    """
    Function which opens and reads the pickle files
    outputs: Data frames
    """

    filename = glob.glob(file_string)
    with open(os.path.join(max(filename, key=os.path.getctime)), 'rb') as input_file:
        temp_data = pkl.load(input_file)
    temp_length = len(temp_data)

    return temp_data, temp_length

def file_manipulator(col_len, data):
    """
    Function to read POSITIONAL arrays of SINGLE components.
    Manipulates them so they can be read by plotters.
    outputs: Manipulatable data
    """

    temp_x = np.empty((1, col_len, 1))
    temp_y = np.empty((1, col_len, 1))
    temp_z = np.empty((1, col_len, 1))

    for i in range(col_len):
        temp_vals = data.iloc[i]
        temp_x[0][i][0] = temp_vals[0].value_in(units.AU)
        temp_y[0][i][0] = temp_vals[1].value_in(units.AU)
        temp_z[0][i][0] = temp_vals[2].value_in(units.AU)

    return temp_x, temp_y, temp_z


def spatial_plotter():
    """
    Function to plot the evolution of the system
    outputs: Plots of the evolution of the system in Cartesian coordinates
    """

    SMBH_code = MW_SMBH()

    com_tracker, col_len  = file_opener('data/center_of_mass/*')
    IMBH_tracker, col_len = file_opener('data/positions_IMBH/*')
    GC_tracker, col_len   = file_opener('data/positions_GC/*')
    
    line_x = np.empty((len(IMBH_tracker), col_len, 1))
    line_y = np.empty((len(IMBH_tracker), col_len, 1))
    line_z = np.empty((len(IMBH_tracker), col_len, 1))

    com_x, com_y, com_z = file_manipulator(col_len, com_tracker)
    GC_x, GC_y, GC_z = file_manipulator(col_len, GC_tracker)

    for i in range(len(IMBH_tracker)):
        tIMBH_tracker = IMBH_tracker.iloc[i]
        
        for j in range(col_len):
            coords = tIMBH_tracker.iloc[j+1][0]
            line_x[i][j][0] = coords[0].value_in(units.AU)
            line_y[i][j][0] = coords[1].value_in(units.AU)
            line_z[i][j][0] = coords[2].value_in(units.AU)

    fig = plt.figure(figsize=(12.5, 8))
    ax1 = fig.add_subplot(221)
    ax2 = fig.add_subplot(222)
    ax3 = fig.add_subplot(212)

    colors = get_distinct(len(IMBH_tracker))
    ax1.set_title('Overall System')
    ax2.set_title('Focus On One IMBH')
    
    ax1.set_xlabel(r'$x$-Coordinate [AU]')
    ax1.set_ylabel(r'$z$-Coordinate [AU]')
    for ax_ in [ax1, ax2]:
        ax_.set_xlabel(r'$x$-Coordinate [AU]')
        ax_.set_ylabel(r'$y$-Coordinate [AU]')

    xy_lim = 0.001 | units.parsec # Hard-coded value but corresponds to the shift
    xy_lim = 1.5* xy_lim.value_in(units.AU)
    ax3.set_xlim(-xy_lim, xy_lim)
    ax3.set_ylim(-xy_lim, xy_lim)

    lim_x = [ ]
    lim_y = [ ]

    min_x = min((line_x[:]-com_x[0])[0][0])
    max_x = max((line_x[:]-com_x[0])[0][0])
    lim_x.append(max(abs(min_x), abs(max_x)))

    min_y = min((line_y[:]-com_y[0])[0][0])
    max_y = max((line_y[:]-com_y[0])[0][0])
    lim_y.append(max(abs(min_y), abs(max_y)))
    
    ax3.plot(GC_x[0][:], GC_y[0][:], color = 'black', label = 'Globular Cluster', linestyle = '--')
    ax3.scatter(SMBH_code.position.x.value_in(units.AU), SMBH_code.position.y.value_in(units.AU), 
                color = 'black', s = 1000*SMBH_code.bh_rad.value_in(units.AU), label = r'SMBH [$M=4$e$6 M_{\odot}$]' )
    for i in range(len(IMBH_tracker)):
        ax3.plot(line_x[i][:], line_y[i][:], color = colors[i])
        ax1.plot(line_x[i][:]-com_x[0][:], line_y[i][:]-com_y[0][:], c = colors[i], lw = 1.4)
        ax1.scatter(line_x[i][0]-com_x[0][0], line_y[i][0]-com_y[0][0], 
                    alpha = 0.7, c = colors[i], edgecolors = 'black', s = 50)
        ax1.scatter(line_x[i][-1]-com_x[0][-1], line_y[i][-1]-com_y[0][-1], 
                    c = colors[i], edgecolors = 'black', s = 50)
        ax2.plot(line_x[i][:]-line_x[0][:], line_y[i][:]-line_y[0][:], c = colors[i], lw = 1.4)
        ax2.scatter(line_x[i][0]-line_x[0][0], line_y[i][0]-line_y[0][0], 
                    alpha = 0.7, c = colors[i], edgecolors = 'black', s = 50)
        ax2.scatter(line_x[i][-1]-line_x[0][-1], line_y[i][-1]-line_y[0][-1], 
                    c = colors[i], edgecolors = 'black', s = 50)
    ax2.scatter(0, 0, c = colors[0], edgecolors = 'black', s = 50)
    plt.legend()
    plt.savefig('figures/spatial_tracker'+str(datetime.now())+'.pdf', dpi=300, bbox_inches='tight')
    plt.clf()
    plt.close()

    return

def energy_plotter():
    """
    Function to plot the evolution of the system
    outputs: Plots of the evolution of the system in Cartesian coordinates
    """

    energy_tracker, col_len = file_opener('data/energy/*')

    iteration = np.empty((1, col_len, 1))
    Et_array  = np.empty((1, col_len, 1))
    dE_array  = np.empty((1, col_len, 1))
    dEs_array = np.empty((1, col_len, 1))

    for i in range(col_len):
        vals = energy_tracker.iloc[i]
        iteration[0][i][0] = vals[0]
        Et_array[0][i][0]  = vals[1].value_in(units.J)
        dE_array[0][i][0]  = vals[2]
        dEs_array[0][i][0] = vals[3]

    fig = plt.figure(figsize=(12, 4))
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)
    ax1.set_title('Energy Conservation')
    ax2.set_title('Stabilised Energy Conservation')

    ax1.set_xlabel(r'Time Step [$\eta$]')
    ax1.set_ylabel(r'$\frac{|E(t)-E_0|}{|E_0|}$')
    ax2.set_xlabel(r'Time Step [$\eta$]')
    ax1.set_yscale('log')
    ax2.set_yscale('log')
        
    ax1.plot(iteration[0][:],   dE_array[0][:])
    ax2.plot(iteration[0][17:], dEs_array[0][17:])
    plt.savefig('figures/energy_tracker'+str(datetime.now())+'.pdf', dpi=300, bbox_inches='tight')
    plt.clf()
    plt.close()

    return


def animator(tend, eta):
    """
    Function to produce animations. WARNING: SLOW
    
    Inputs:
    tend:    The total time the simulation evolves till
    eta:     The time-step used for integration
    outputs: Movie of the animation
    """

    print('You have chosen to animate. This will take a while...')

    plt.clf()
    SMBH_code = MW_SMBH()

    com_tracker, col_len  = file_opener('data/center_of_mass/*')
    IMBH_tracker, col_len = file_opener('data/positions_IMBH/*')
    GC_tracker,   col_len = file_opener('data/positions_GC/*')
    energy_tracker, col_len = file_opener('data/energy/*')

    iteration = np.empty((1, col_len, 1))
    dE_array  = np.empty((1, col_len, 1))

    for i in range(col_len):
        vals = energy_tracker.iloc[i]
        iteration[0][i][0] = vals[0]
        dE_array[0][i][0]  = vals[2]

    col_len = len(IMBH_tracker.iloc[0])-2
    
    line_x = np.empty((len(IMBH_tracker), col_len, 1))
    line_y = np.empty((len(IMBH_tracker), col_len, 1))
    line_z = np.empty((len(IMBH_tracker), col_len, 1))
    
    com_x, com_y, com_z = file_manipulator(col_len, com_tracker)
    GC_x, GC_y, GC_z = file_manipulator(col_len, GC_tracker)

    for i in range(len(IMBH_tracker)):
        tIMBH_tracker = IMBH_tracker.iloc[i]
        
        for j in range(col_len):
            coords = tIMBH_tracker.iloc[j+1][0]
            line_x[i][j][0] = coords[0].value_in(units.AU)
            line_y[i][j][0] = coords[1].value_in(units.AU)
            line_z[i][j][0] = coords[2].value_in(units.AU)

    colors = get_distinct(len(IMBH_tracker))

    fig = plt.figure(figsize=(12.5, 8))
    ax1 = fig.add_subplot(221)
    ax2 = fig.add_subplot(222)
    ax3 = fig.add_subplot(223)
    ax4 = fig.add_subplot(224)

    def animate(col_len):
        if col_len %100 == 0:
            print('Current frame: ', col_len)

        xy_lim = 0.001 | units.parsec # Hard-coded value but corresponds to the shift
        xy_lim = 1.1* xy_lim.value_in(units.AU)

        ax1.clear()
        ax2.clear()
        ax3.clear()
        ax4.clear()

        for i in range(len(IMBH_tracker)):
            ax1.plot(line_x[i][0:col_len+1], line_y[i][0:col_len+1], c = colors[i])

            ax3.plot(line_x[i][max(0,col_len+1-70):col_len+1]-com_x[0][max(0,col_len+1-70):col_len+1],
                     line_y[i][max(0,col_len+1-70):col_len+1]-com_y[0][max(0,col_len+1-70):col_len+1],
                     c = colors[i], lw = 2)

            ax3.scatter(line_x[i][col_len]-com_x[0][col_len], line_y[i][col_len]-com_y[0][col_len],
                        c = colors[i], edgecolors = 'black', s = 40)

            ax4.plot(line_x[i][max(0,col_len+1-70):col_len+1]-com_x[0][max(0,col_len+1-70):col_len+1],
                     line_z[i][max(0,col_len+1-70):col_len+1]-com_z[0][max(0,col_len+1-70):col_len+1],
                     c = colors[i], lw = 2)
            ax4.scatter(line_x[i][col_len]-com_x[0][col_len], line_z[i][col_len]-com_z[0][col_len],
                        c = colors[i], edgecolors = 'black', s = 40)

        ax1.plot(GC_x[0][0:col_len+1], GC_y[0][0:col_len+1], c = 'black', linestyle = '--')

        ax1.set_title(str("{:.3f}".format(iteration[0][col_len][0]*tend.number*eta*365)+" Days"), loc = 'left')
        ax1.scatter(SMBH_code.position.x.value_in(units.AU), SMBH_code.position.y.value_in(units.AU), 
                    color = 'black', s = SMBH_code.bh_mass.number**0.3)
        ax1.set_xlim(-xy_lim, xy_lim)
        ax1.set_ylim(-xy_lim, xy_lim)     
        ax1.set_xlabel(r'$x$-Coordinate [AU]')
        ax1.set_ylabel(r'$y$-Coordinate [AU]')

        ax2.set_xlabel(r'Time Step [$\eta$]')
        ax2.set_ylabel(r'$\frac{|E(t)-E_0|}{|E_0|}$')
        ax2.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.0e'))
        ax2.plot(iteration[0][:col_len+1], abs(dE_array[0][:col_len+1]))
        ax2.set_yscale('log')

        ax3.set_xlabel(r'$x$-Coordinate [AU]')
        ax3.set_ylabel(r'$y$-Coordinate [AU]')

        ax4.set_xlabel(r'$x$-Coordinate [AU]')
        ax4.set_ylabel(r'$z$-Coordinate [AU]')

    anim = animation.FuncAnimation(fig, animate, interval = 100, frames=col_len)
    anim.save('figures/animation'+str(datetime.now())+'.gif', writer='pillow')
    return 
