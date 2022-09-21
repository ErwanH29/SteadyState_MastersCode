from amuse.lab import *
from parti_initialiser import *
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import matplotlib.animation as animation
from amuse.ext.distinct_colours import get_distinct
import pickle as pkl
import numpy as np
import glob
import os
from itertools import cycle
from evol_func import file_counter

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

def animator():
    """
    Function to produce animations. WARNING: SLOW
    
    Inputs:
    tend:    The total time the simulation evolves till
    eta:     The time-step used for integration
    outputs: Movie of the animation
    """

    print('You have chosen to animate. This will take a while...')

    plt.clf()
    count = file_counter()

    Lag_tracker, col_len  = file_opener('data/lagrangians/*')
    com_tracker, col_len  = file_opener('data/center_of_mass/*')
    IMBH_tracker, col_len = file_opener('data/positions_IMBH/*')
    energy_tracker, col_len = file_opener('data/energy/*')

    time = np.empty((1, col_len, 1))
    dE_array  = np.empty((1, col_len, 1))

    for i in range(col_len):
        vals = energy_tracker.iloc[i]
        time[0][i][0] = vals[0].value_in(units.Myr)
        dE_array[0][i][0]  = vals[2]

    col_len = len(IMBH_tracker.iloc[0])-2
    
    com_x, com_y, com_z = file_manipulator(col_len, com_tracker)
    line_x = np.empty((len(IMBH_tracker), col_len, 1))
    line_y = np.empty((len(IMBH_tracker), col_len, 1))
    line_z = np.empty((len(IMBH_tracker), col_len, 1))

    for i in range(len(IMBH_tracker)-1):
        tIMBH_tracker = IMBH_tracker.iloc[i+1]
        tIMBH_tracker = tIMBH_tracker.replace(np.NaN, "[Np.NaN, [np.NaN, np.NaN, np.NaN]")
        for j in range(col_len):
            coords = tIMBH_tracker.iloc[j+1][1]
            if len(coords) == 1:
                pass
            else:
                line_x[i][j][0] = coords[0].value_in(units.AU)
                line_y[i][j][0] = coords[1].value_in(units.AU)
                line_z[i][j][0] = coords[2].value_in(units.AU)

    line_x[:][abs(line_x[:]) < 10**-5] = np.NaN
    line_y[:][abs(line_y[:]) < 10**-5] = np.NaN
    line_z[:][abs(line_z[:]) < 10**-5] = np.NaN

    line_x[:][abs(line_x[:]) > 1e8] = np.NaN
    line_y[:][abs(line_y[:]) > 1e8] = np.NaN
    line_z[:][abs(line_z[:]) > 1e8] = np.NaN 

    LG25_array  = np.empty((1, col_len, 1))
    LG50_array  = np.empty((1, col_len, 1))
    LG75_array  = np.empty((1, col_len, 1))

    for i in range(col_len):
        vals = Lag_tracker.iloc[i]
        LG25_array[0][i][0] = vals[1].value_in(units.AU)
        LG50_array[0][i][0] = vals[2].value_in(units.AU)
        LG75_array[0][i][0] = vals[3].value_in(units.AU)

    fig = plt.figure(figsize=(12.5, 8))
    ax1 = fig.add_subplot(221)
    ax2 = fig.add_subplot(222)
    ax3 = fig.add_subplot(223)
    ax4 = fig.add_subplot(224)

    for ax_ in [ax1, ax2, ax3, ax4]:
        ax_.xaxis.set_ticks_position('both')
        ax_.yaxis.set_ticks_position('both')
        ax_.tick_params(axis="y",direction="in")
        ax_.tick_params(axis="x",direction="in")

    def animate(col_len):
        if col_len %100 == 0:
            print('Current frame: ', col_len)

        xy_lim = 0.001 | units.parsec # Hard-coded value but corresponds to the shift
        xy_lim = 1.1* xy_lim.value_in(units.AU)

        ax1.clear()
        ax2.clear()
        ax3.clear()
        ax4.clear()

        ax1.set_xlim(-10000, 10000)
        ax1.set_ylim(-10000, 10000)

        for i in range(len(IMBH_tracker)):
            ax1.plot(line_x[i][max(0,col_len+1-1):col_len+1]-line_x[1][max(0,col_len+1-1):col_len+1],
                     line_y[i][max(0,col_len+1-1):col_len+1]-line_x[1][max(0,col_len+1-1):col_len+1],
                     c = 'black', lw = 2)
            ax1.scatter(line_x[i][col_len]-line_x[1][max(0,col_len+1-1):col_len+1], line_y[i][col_len]-line_x[1][max(0,col_len+1-1):col_len+1],
                        c = 'black', edgecolors = 'black', s = 40)

            ax3.plot(line_x[i][max(0,col_len+1-1):col_len+1], line_y[i][max(0,col_len+1-1):col_len+1], c = 'black', lw = 2)
            ax3.scatter(line_x[i][col_len], line_y[i][col_len], c = 'black', edgecolors = 'black', s = 40)

            ax4.plot(time[0][:col_len+1], LG25_array[0][:col_len+1], color = 'blue', label = r'$r_{25,L}$')
            ax4.plot(time[0][:col_len+1], LG75_array[0][:col_len+1], color = 'red',  label = r'$r_{75,L}$')

        ax2.plot(time[0][:col_len+1], abs(dE_array[0][:col_len+1]))
        ax1.set_title(str("{:.3f}".format(time[0][col_len][0])+" Myr"), loc = 'left')   
        ax1.set_xlabel(r'$y$-Coordinate [AU]')
        ax1.set_ylabel(r'$z$-Coordinate [AU]')

        ax2.set_xlabel(r'Time [Myr]')
        ax2.set_ylabel(r'$\frac{|E(t)-E_0|}{|E_0|}$')
        ax2.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.0e'))
        ax2.set_yscale('log')

        ax3.set_xlabel(r'$x$-Coordinate [AU]')
        ax3.set_ylabel(r'$y$-Coordinate [AU]')

        ax4.set_xlabel(r'Time [Myr]')
        ax4.set_ylabel(r'Lagrangian Radius [AU]')

    anim = animation.FuncAnimation(fig, animate, interval = 100, frames=col_len)
    anim.save('figures/animation'+str(count)+'.gif', writer='pillow')
    return 
    
def energy_plotter():
    """
    Function to plot the evolution of the system
    outputs: Plots of the evolution of the system in Cartesian coordinates
    """

    count = file_counter()

    energy_tracker, col_len = file_opener('data/energy/*')
    IMBH_energy_tracker, col_len = file_opener('data/particle_energies/*')

    time = np.empty((1, col_len, 1))
    Et_array    = np.empty((1, col_len, 1))
    dE_array    = np.empty((1, col_len, 1))
    dEs_array   = np.empty((1, col_len, 1))
    IMBH_birth  = np.empty((1, col_len, 1))
    collisions  = np.empty((1, col_len, 1)) 
    merger_mass = np.empty((1, col_len, 1)) 

    BE_array   = np.empty((1, col_len, 1))
    KE_array   = np.empty((1, col_len, 1))
    TotE_array = np.empty((1, col_len, 1))

    for i in range(col_len):
        vals = energy_tracker.iloc[i]
        time[0][i][0] = vals[0].value_in(units.Myr)
        Et_array[0][i][0]    = vals[1].value_in(units.J)
        dE_array[0][i][0]    = vals[2]
        dEs_array[0][i][0]   = vals[3]
        IMBH_birth[0][i][0]  = vals[4].value_in(units.Myr)
        collisions[0][i][0]  = vals[5].value_in(units.Myr)
        merger_mass[0][i][0] = vals[6].value_in(units.MSun)

        vals = IMBH_energy_tracker.iloc[i]
        BE_array[0][i][0]    = vals[1].value_in(units.J)
        KE_array[0][i][0]    = vals[2].value_in(units.J)
        TotE_array[0][i][0]  = vals[3].value_in(units.J)    

    coll_id = np.argwhere(collisions > time[0][1])
    coll_id[:,0] = [i-1 for i in coll_id[:,0]]
    app_id  = np.argwhere(IMBH_birth[0]  > time[0][1])
    app_id[:,0]  = [i-1 for i in app_id[:,0]]

    fig = plt.figure(figsize=(12.5, 8))
    ax1 = fig.add_subplot(221)
    ax2 = fig.add_subplot(222)
    ax1.set_title('Energy Conservation')
    ax2.set_title('Total Energy Evolution')
    for ax_ in [ax1, ax2]:
        ax_.yaxis.set_ticks_position('both')
        ax_.xaxis.set_ticks_position('both')
        ax_.tick_params(axis="y",direction="in")
        ax_.tick_params(axis="x",direction="in")
        ax_.set_xlabel(r'Time [Myr]')
        ax_.set_yscale('log')
    ax1.set_ylabel(r'$\frac{|E(t)-E_0|}{|E_0|}$')
    ax2.set_ylabel(r'Energy [J]')
    ax2.set_ylim(0.1*min(abs(BE_array[0][:])), 5*max(KE_array[0][:]))

    ax1.plot(time[0][:], dE_array[0][:], color = 'black', zorder = 1)
    if len(coll_id) > 0:
        dE_coll = dE_array[0][coll_id[:,1]-1]
        collisions = collisions[0][coll_id[:,1]]
        merger_mass = merger_mass[0][coll_id[:,1]]
        color_axes = ax1.scatter(collisions, dE_coll, c = merger_mass, zorder=3)
        plt.colorbar(color_axes, ax=ax1, label = r'Merger Mass [$M_{\odot}$]')
    ax2.plot(time[0][:], abs(BE_array[0][:]), color = 'blue', label = 'Potential Energy (abs)', zorder = 1)
    ax2.plot(time[0][:], KE_array[0][:], color = 'red', label = 'Kinetic Energy', zorder = 2)
    ax2.plot(time[0][:], abs(TotE_array[0][:]), color = 'black', label = 'Total Energy (abs)', linestyle = '--', zorder = 3)
    ax2.legend()

    if len(app_id) > 0:
        dE_app = dE_array[0][app_id[:,0]+1]
        IMBH_app = IMBH_birth[0][app_id[:,0]+1]
        totE_app = abs(TotE_array[0][app_id[:,0]+1])
        ax1.scatter(IMBH_app, dE_app, marker = 'x', color = 'red', zorder = 2)
        ax2.scatter(IMBH_app, totE_app, color = 'black', zorder = 4)
    plt.savefig('figures/energy_tracker'+str(count)+'.pdf', dpi=300, bbox_inches='tight')
    plt.clf()
    plt.close()

    return

def spatial_plotter(init_dist):
    """
    Function to plot the evolution of the system
    outputs: Plots of the evolution of the system in Cartesian coordinates
    """

    count = file_counter()
    com_tracker, col_len  = file_opener('data/center_of_mass/*')
    IMBH_tracker, col_len = file_opener('data/positions_IMBH/*')
    tdyn_tracker, col_len = file_opener('data/dynamical_time/*')
    Lag_tracker, col_len  = file_opener('data/lagrangians/*')
    com_x, com_y, com_z   = file_manipulator(col_len, com_tracker)

    time = np.empty((1, col_len, 1))
    LG25_array  = np.empty((1, col_len, 1))
    LG50_array  = np.empty((1, col_len, 1))
    LG75_array  = np.empty((1, col_len, 1))

    for i in range(col_len):
        vals = Lag_tracker.iloc[i]
        time[0][i][0] = vals[0].value_in(units.Myr)
        LG25_array[0][i][0] = vals[1].value_in(units.pc)
        LG50_array[0][i][0] = vals[2].value_in(units.pc)
        LG75_array[0][i][0] = vals[3].value_in(units.pc)

    line_x = np.empty((len(IMBH_tracker), col_len, 1))
    line_y = np.empty((len(IMBH_tracker), col_len, 1))
    line_z = np.empty((len(IMBH_tracker), col_len, 1))
    tdyn   = np.empty((len(IMBH_tracker), col_len, 1))

    colors = get_distinct(12)

    for i in range(len(IMBH_tracker)-1):
        tIMBH_tracker = IMBH_tracker.iloc[i+1]
        tIMBH_tracker = tIMBH_tracker.replace(np.NaN, "[Np.NaN, [np.NaN, np.NaN, np.NaN]")
        tDyntime_trk  = tdyn_tracker.iloc[i+1]
        for j in range(col_len):
            coords = tIMBH_tracker.iloc[j+1][1]
            if len(coords) == 1:
                pass
            else:
                line_x[i][j][0] = coords[0].value_in(units.AU)
                line_y[i][j][0] = coords[1].value_in(units.AU)
                line_z[i][j][0] = coords[2].value_in(units.AU)

                tdynval = tDyntime_trk.iloc[j+1][0]
                tdyn[i][j][0] = tdynval[0].value_in(units.yr)

    line_x[:][abs(line_x[:]) < 10**-5] = np.NaN
    line_y[:][abs(line_y[:]) < 10**-5] = np.NaN
    line_z[:][abs(line_z[:]) < 10**-5] = np.NaN

    line_x[:][abs(line_x[:]) > 1e8] = np.NaN
    line_y[:][abs(line_y[:]) > 1e8] = np.NaN
    line_z[:][abs(line_z[:]) > 1e8] = np.NaN

    tdyn[:][tdyn[:] < 10**-30] = np.NaN
    tdyn[:][tdyn[:] > 10**10] = np.NaN

    fig = plt.figure(figsize=(12.5, 8))
    ax1 = fig.add_subplot(221)
    ax2 = fig.add_subplot(222)
    ax3 = fig.add_subplot(223)
    ax4 = fig.add_subplot(224)

    colorscycler = cycle(colors)
    
    ax1.set_title('Center of Mass TO FIX')
    ax2.set_title('Ejected IMBH Focus')

    for ax_ in [ax1, ax2, ax3]:
        ax_.set_xlabel(r'$x$-Coordinate [AU]')
        ax_.set_ylabel(r'$y$-Coordinate [AU]')
    for ax_ in [ax1, ax2, ax3, ax4]:
        ax_.yaxis.set_ticks_position('both')
        ax_.xaxis.set_ticks_position('both')
        ax_.tick_params(axis="y",direction="in")
        ax_.tick_params(axis="x",direction="in")
 
    ax2.set_xlim(-100,100)
    ax2.set_ylim(-100,100)        
    ax3.set_xlim(-init_dist.value_in(units.AU)**1.01, init_dist.value_in(units.AU)**1.01)
    ax3.set_ylim(-init_dist.value_in(units.AU)**1.01, init_dist.value_in(units.AU)**1.01)        

    ax4.set_xlabel(r'Time [Myr]')
    ax4.set_ylabel(r'Lagrangian Radius [pc]')

    for i in range(len(IMBH_tracker)):
        colors = next(colorscycler)
        ax1.scatter(line_x[i][:]-com_x[0][:], line_y[i][:]-com_y[0][:], c = colors, zorder = 1)
        ax1.scatter(line_x[i][0]-com_x[0][0], line_y[i][0]-com_y[0][0], 
                    alpha = 0.7, c = colors, edgecolors = 'black', s = 50, zorder = 2)
        ax1.scatter(line_x[i][-1]-com_x[0][-1], line_y[i][-1]-com_y[0][-1], 
                    c = colors, edgecolors = 'black', s = 50, zorder = 3)
        ax3.scatter(line_x[i][:-1], line_y[i][:-1], c = colors, zorder = 1, s = 5)
        ax3.scatter(line_x[i][-1], line_y[i][-1], c = colors, edgecolors = 'black', s = 50, zorder = 3)
        ax3.scatter(line_x[i][0], line_y[i][0], alpha = 0.7, c = colors, edgecolors = 'black', s = 50, zorder = 2)
        ax2.plot(line_x[i][:]-line_x[1][:], line_y[i][:]-line_y[1][:], c = colors, lw = 3)
        ax2.scatter(line_x[i][0]-line_x[1][0], line_y[i][0]-line_y[1][0], 
                    alpha = 0.7, c = colors, edgecolors = 'black', s = 50)
        ax2.scatter(line_x[i][-1]-line_x[1][-1], line_y[i][-1]-line_y[1][-1], 
                    c = colors, edgecolors = 'black', s = 50)
    ax3.scatter(0,0, color = 'black', s = 50, zorder = 2)
    ax4.plot(time[0][:], LG25_array[0][:], color = 'red',   label = r'$r_{25,L}$')
    ax4.plot(time[0][:], LG50_array[0][:], color = 'black', label = r'$r_{25,L}$')
    ax4.plot(time[0][:], LG75_array[0][:], color = 'blue',  label = r'$r_{75,L}$')
    ax4.legend()
    plt.savefig('figures/spatial_tracker'+str(count)+'.pdf', dpi=300, bbox_inches='tight')
    plt.clf()
    plt.close()

    return

spatial_plotter(0.1 | units.parsec)