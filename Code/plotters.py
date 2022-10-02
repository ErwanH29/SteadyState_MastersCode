from sqlite3 import adapt
from amuse.lab import *
from parti_initialiser import *
from file_logistics import *
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import matplotlib.animation as animation
import numpy as np
from itertools import cycle

def colour_picker():
    """
    Colour chooser for the spatial plots and animation
    """

    colors = ['red', 'blue', 'orange', 'purple', 'salmon', 'slateblue', 
              'gold', 'darkviolet', 'cornflowerblue',  'cyan',
              'lightskyblue', 'magenta',  'dodgerblue']

    return colors

def animator(init_dist):
    """
    Function to produce animations. WARNING: SLOW
    
    Inputs:
    init_dist: The initial distance of the cluster to the central SMBH
    """

    print('!!!!! You have chosen to animate. This will take a while !!!!!')
    plt.clf()
    count = file_counter()

    Lag_tracker, col_len = file_opener('data/lagrangians/*')
    com_tracker, col_len = file_opener('data/center_of_mass/*')
    IMBH_tracker, col_len = file_opener('data/positions_IMBH/*')
    energy_tracker, col_len = file_opener('data/energy/*')
    com_x, com_y, com_z = file_manipulator(col_len, com_tracker)

    time = np.empty((col_len, 1))
    dE_array = np.empty((col_len, 1))

    for i in range(col_len):
        vals = energy_tracker.iloc[i]
        time[i][0] = vals[0].value_in(units.Myr)
        dE_array[i][0] = vals[2]

    col_len = len(IMBH_tracker.iloc[0])-2
    line_x = np.empty((len(IMBH_tracker), col_len, 1))
    line_y = np.empty((len(IMBH_tracker), col_len, 1))
    line_z = np.empty((len(IMBH_tracker), col_len, 1))

    for i in range(len(IMBH_tracker)):
        tIMBH_tracker = IMBH_tracker.iloc[i]
        tIMBH_tracker = tIMBH_tracker.replace(np.NaN, "[Np.NaN, [np.NaN, np.NaN, np.NaN]")
        for j in range(col_len):
            coords = tIMBH_tracker.iloc[j+1][1]
            if len(coords) == 1:
                pass
            else:
                line_x[i][j][0] = coords[0].value_in(units.pc)
                line_y[i][j][0] = coords[1].value_in(units.pc)
                line_z[i][j][0] = coords[2].value_in(units.pc)

    temp_coord_array = [line_x, line_y, line_z]
    for arr_ in temp_coord_array:
        arr_[:][abs(arr_[:]) < 10**-5] = np.NaN
        arr_[:][abs(arr_[:]) > 1e8] = np.NaN

    LG25_array = np.empty((col_len, 1))
    LG50_array = np.empty((col_len, 1))
    LG75_array = np.empty((col_len, 1))
    rtid_array = np.empty((col_len, 1))

    for i in range(col_len):
        vals = Lag_tracker.iloc[i]
        LG25_array[i][0] = vals[1].value_in(units.pc)
        LG50_array[i][0] = vals[2].value_in(units.pc)
        LG75_array[i][0] = vals[3].value_in(units.pc)
        rtid_array[i][0] = vals[4].value_in(units.pc)

    totsys_lim = init_dist.value_in(units.pc)

    plt.ioff()
    fig = plt.figure()
    ax3D1 = fig.add_subplot(121, projection="3d")
    ax3D2 = fig.add_subplot(122, projection="3d")
    colours = colour_picker()

    def animate_3D(col_len):
        skip_zeroth = col_len + 1
        iter = -1

        if skip_zeroth %100 == 0:
            print('Current 3D frame: ', skip_zeroth)

        ax3D1.clear()
        ax3D2.clear()
        ax3D1.set_title(str("{:.3f}".format(time[skip_zeroth][0])+" Myr"), loc = 'left') 
        ax3D1.set_xlim3d(-1.05*totsys_lim, 1.05*totsys_lim)
        ax3D1.set_ylim3d(-1.05*totsys_lim, 1.05*totsys_lim)
        ax3D1.set_zlim3d(-1.05*totsys_lim, 1.05*totsys_lim)

        for i in range(len(IMBH_tracker)):
            iter += 1
            if iter > len(colours):
                iter = 0
                
            if i == 0:
                pass
            if i == 1:
                ax3D1.scatter((line_x[i][max(0,skip_zeroth-20):skip_zeroth+1]-line_x[0][max(0,skip_zeroth-20):skip_zeroth+1]), 
                              (line_y[i][max(0,skip_zeroth-20):skip_zeroth+1]-line_y[0][max(0,skip_zeroth-20):skip_zeroth+1]), 
                              (line_z[i][max(0,skip_zeroth-20):skip_zeroth+1]-line_z[0][max(0,skip_zeroth-20):skip_zeroth+1]), 
                               ls=':', c = 'black', lw = 2)

            if i > 1:
                ax3D1.scatter((line_x[i][max(0,skip_zeroth-20):skip_zeroth+1]-line_x[0][max(0,skip_zeroth-20):skip_zeroth+1]), 
                              (line_y[i][max(0,skip_zeroth-20):skip_zeroth+1]-line_y[0][max(0,skip_zeroth-20):skip_zeroth+1]), 
                              (line_z[i][max(0,skip_zeroth-20):skip_zeroth+1]-line_z[0][max(0,skip_zeroth-20):skip_zeroth+1]), 
                              s =1, c = colours[iter-2], lw = 2)
                ax3D1.scatter((line_x[i][skip_zeroth]-line_x[0][skip_zeroth]), 
                              (line_y[i][skip_zeroth]-line_y[0][skip_zeroth]), 
                              (line_z[i][skip_zeroth]-line_z[0][skip_zeroth]), 
                               c = colours[iter-2], edgecolors = 'black', s = 50)

                ax3D2.scatter((line_x[i][max(0,skip_zeroth-20):skip_zeroth+1]-com_x[max(0,skip_zeroth-20):skip_zeroth+1]),
                              (line_y[i][max(0,skip_zeroth-20):skip_zeroth+1]-com_y[max(0,skip_zeroth-20):skip_zeroth+1]),
                              (line_z[i][max(0,skip_zeroth-20):skip_zeroth+1]-com_z[max(0,skip_zeroth-20):skip_zeroth+1]), 
                               s =1, c = colours[iter-2], lw = 2)
                ax3D2.scatter((line_x[i][skip_zeroth]-com_x[skip_zeroth]), 
                              (line_y[i][skip_zeroth]-com_y[skip_zeroth]), 
                              (line_z[i][skip_zeroth]-com_z[skip_zeroth]),
                               c = colours[iter-2], edgecolors = 'black', s = 40)

        for ax_ in [ax3D1, ax3D2]:
            ax_.xaxis.set_major_formatter(mtick.FormatStrFormatter('%0.2f'))
            ax_.yaxis.set_major_formatter(mtick.FormatStrFormatter('%0.2f'))
            ax_.zaxis.set_major_formatter(mtick.FormatStrFormatter('%0.2f'))
            
        ax3D1.scatter(0,0, color = 'black', s = 100 )
        ax3D1.set_xlabel(r'$x$-Coordinate [pc]')
        ax3D2.set_xlabel(r'$x$-Coordinate [pc]')
        ax3D2.set_ylabel(r'$y$-Coordinate [pc]')

    anim3D = animation.FuncAnimation(fig, animate_3D, interval = 100, frames=(col_len-1))
    anim3D.save('figures/animation3D'+str(count)+'.gif', writer='pillow')

    fig = plt.figure(figsize=(12.5, 8))
    ax1 = fig.add_subplot(221)
    ax2 = fig.add_subplot(222)
    ax3 = fig.add_subplot(223)
    ax4 = fig.add_subplot(224)

    for ax_ in [ax1, ax2, ax3, ax4]:
        ax_.xaxis.set_ticks_position('both')
        ax_.yaxis.set_ticks_position('both')
        ax_.xaxis.set_minor_locator(mtick.AutoMinorLocator())
        ax_.yaxis.set_minor_locator(mtick.AutoMinorLocator())

    def animate(col_len):
        skip_zeroth = col_len + 1
        if skip_zeroth %100 == 0:
            print('Current frame: ', skip_zeroth)
            
        for ax_ in [ax1, ax2, ax3, ax4]:
            ax_.clear()
            ax_.yaxis.set_ticks_position('both')
            ax_.xaxis.set_ticks_position('both')
            ax_.xaxis.set_minor_locator(mtick.AutoMinorLocator())
            ax_.yaxis.set_minor_locator(mtick.AutoMinorLocator())
            ax_.tick_params(axis="y", which = 'both', direction="in")
            ax_.tick_params(axis="x", which = 'both', direction="in")    
        
        for i in range(len(IMBH_tracker)):
            ax1.scatter(line_x[i][max(0,skip_zeroth-20):skip_zeroth+1]-com_x[max(0,skip_zeroth-20):skip_zeroth+1],
                        line_y[i][max(0,skip_zeroth-20):skip_zeroth+1]-com_y[max(0,skip_zeroth-20):skip_zeroth+1],
                        c = colours[i], lw = 2, alpha = 0.7, s = 1)
            ax1.scatter(line_x[i][skip_zeroth]-com_x[skip_zeroth], line_y[i][skip_zeroth]-com_y[skip_zeroth],
                        c = colours[i], edgecolors = 'black', s = 40)

            ax3.scatter((line_x[i][max(0,skip_zeroth-20):skip_zeroth+1]-line_x[0][max(0,skip_zeroth-20):skip_zeroth+1]), 
                        (line_y[i][max(0,skip_zeroth-20):skip_zeroth+1]-line_y[0][max(0,skip_zeroth-20):skip_zeroth+1]), 
                        s =1, c = colours[i], lw = 2)
            ax3.scatter((line_x[i][skip_zeroth]-line_x[0][skip_zeroth]), 
                        (line_y[i][skip_zeroth]-line_y[0][skip_zeroth]), 
                        c = colours[i], edgecolors = 'black', s = 50)

            ax4.plot(time[3:skip_zeroth+1], (rtid_array[3:skip_zeroth+1]), color = 'black', linestyle = ':', label = r'$r_{tidal}$')
            ax4.plot(time[3:skip_zeroth+1], (LG25_array[3:skip_zeroth+1]), color = 'red',  label = r'$r_{25,L}$')
            ax4.plot(time[3:skip_zeroth+1], (LG50_array[3:skip_zeroth+1]), color = 'black', label = r'$r_{50,L}$')
            ax4.plot(time[3:skip_zeroth+1], (LG75_array[3:skip_zeroth+1]), color = 'blue',   label = r'$r_{75,L}$')

        ax3.scatter(0,0, color = 'black', s = 100 )
        ax2.plot(time[:skip_zeroth+1], abs(dE_array[:skip_zeroth+1]), color = 'black')

        for ax_ in [ax1, ax3]:
            ax_.set_xlabel(r'$x$-Coordinate [pc]')
            ax_.set_ylabel(r'$y$-Coordinate [pc]')

        for ax_ in [ax2, ax4]:
            ax_.set_xlabel(r'Time [Myr]')

        ax1.set_title(str("{:.0f}".format(time[skip_zeroth][0])+" Myr"), loc = 'left')   
        ax1.yaxis.set_major_formatter(mtick.FormatStrFormatter('%0.3f'))
        ax1.xaxis.set_major_formatter(mtick.FormatStrFormatter('%0.3f'))
        ax1.set_xlim(min(line_x[i][max(0,skip_zeroth-80):skip_zeroth+1]-com_x[max(0,skip_zeroth-80):skip_zeroth+1]),
                     max(line_x[i][max(0,skip_zeroth-80):skip_zeroth+1]-com_x[max(0,skip_zeroth-80):skip_zeroth+1]))
        ax1.set_ylim(min(line_y[i][max(0,skip_zeroth-80):skip_zeroth+1]-com_y[max(0,skip_zeroth-80):skip_zeroth+1]),
                     max(line_y[i][max(0,skip_zeroth-80):skip_zeroth+1]-com_y[max(0,skip_zeroth-80):skip_zeroth+1]))

        ax2.set_ylabel(r'$\frac{|E(t)-E_0|}{|E_0|}$')
        ax2.xaxis.set_major_formatter(mtick.FormatStrFormatter('%0.2f'))
        ax2.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.e'))
        ax2.set_yscale('log')

        ax3.yaxis.set_major_formatter(mtick.FormatStrFormatter('%0.1f'))
        ax3.xaxis.set_major_formatter(mtick.FormatStrFormatter('%0.1f'))
        ax3.set_xlim(-1.05*totsys_lim, 1.05*totsys_lim)
        ax3.set_ylim(-totsys_lim, totsys_lim)  

        ax4.set_ylabel(r'Distance [pc]')

    anim = animation.FuncAnimation(fig, animate, interval = 100, frames=(col_len-1))
    anim.save('figures/animation2D'+str(count)+'.gif', writer='pillow')

    return 
        
def energy_plotter():
    """
    Function to plot the energy evolution of the system
    """

    count = file_counter()
    energy_tracker, col_len = file_opener('data/energy/*')
    IMBH_energy_tracker, col_len = file_opener('data/particle_energies/*')

    time = np.empty((1, col_len, 1))
    Et_array = np.empty((1, col_len, 1))
    dE_array = np.empty((1, col_len, 1))
    dEs_array = np.empty((1, col_len, 1))
    IMBH_birth = np.empty((1, col_len, 1))
    collisions = np.empty((1, col_len, 1)) 
    merger_mass = np.empty((1, col_len, 1)) 

    BE_array = np.empty((1, col_len, 1))
    KE_array = np.empty((1, col_len, 1))
    TotE_array = np.empty((1, col_len, 1))

    for i in range(col_len):
        vals = energy_tracker.iloc[i]
        time[0][i][0] = vals[0].value_in(units.Myr)
        Et_array[0][i][0] = vals[1].value_in(units.J)
        dE_array[0][i][0] = vals[2]
        dEs_array[0][i][0] = vals[3]
        IMBH_birth[0][i][0] = vals[4].value_in(units.Myr)
        collisions[0][i][0] = vals[5].value_in(units.Myr)
        merger_mass[0][i][0] = vals[6].value_in(units.MSun)

        vals = IMBH_energy_tracker.iloc[i]
        BE_array[0][i][0] = vals[1].value_in(units.J)
        KE_array[0][i][0] = vals[2].value_in(units.J)
        TotE_array[0][i][0] = vals[3].value_in(units.J)    

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
        ax_.xaxis.set_minor_locator(mtick.AutoMinorLocator())
        ax_.yaxis.set_minor_locator(mtick.AutoMinorLocator())
        ax_.tick_params(axis="y", which = 'both', direction="in")
        ax_.tick_params(axis="x", which = 'both', direction="in")    
        ax_.set_xlabel(r'Time [Myr]')
        ax_.set_yscale('log')
    ax1.set_xlim(time[0][5], time[0][-1])
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
    ax2.plot(time[0][5:], abs(BE_array[0][5:]), color = 'blue', label = 'Potential Energy (abs)', zorder = 1)
    ax2.plot(time[0][5:], KE_array[0][5:], color = 'red', label = 'Kinetic Energy', zorder = 2)
    ax2.plot(time[0][5:], abs(TotE_array[0][5:]), color = 'black', label = 'Total Energy (abs)', linestyle = '--', zorder = 3)
    ax2.legend()

    if len(app_id) > 0:
        dE_app = dE_array[0][app_id[:,0]+1]
        IMBH_app = IMBH_birth[0][app_id[:,0]+1]
        totE_app = abs(TotE_array[0][app_id[:,0]+1])

        ax1.scatter(IMBH_app, dE_app, marker = 'x', color = 'red', zorder = 2)
        ax2.scatter(IMBH_app, totE_app, marker = 'x', color = 'red', zorder = 4)
    plt.savefig('figures/energy_tracker'+str(count)+'.pdf', dpi=300, bbox_inches='tight')
    plt.clf()
    plt.close()

    return

def spatial_plotter(init_dist):
    """
    Function to plot the evolution of the system
    outputs: Plots of the evolution of the system in Cartesian coordinates
    """

    gc_code = globular_cluster()
    count = file_counter()
    ejec_parti, col_len = file_opener('data/chaotic_simulation/*')
    com_tracker, col_len = file_opener('data/center_of_mass/*')
    IMBH_tracker, col_len = file_opener('data/positions_IMBH/*')
    Lag_tracker, col_len = file_opener('data/lagrangians/*')
    com_x, com_y, com_z = file_manipulator(col_len, com_tracker)

    time = np.empty((1, col_len, 1))
    LG25_array  = np.empty((1, col_len, 1))
    LG50_array  = np.empty((1, col_len, 1))
    LG75_array  = np.empty((1, col_len, 1))
    rtide_array = np.empty((1, col_len, 1))
    relax_time  = np.empty((1, col_len, 1))

    for i in range(col_len):
        vals = Lag_tracker.iloc[i]
        time[0][i][0] = vals[0].value_in(units.Myr)
        LG25_array[0][i][0]  = vals[1].value_in(units.pc)
        LG50_array[0][i][0]  = vals[2].value_in(units.pc)
        LG75_array[0][i][0]  = vals[3].value_in(units.pc)
        rtide_array[0][i][0] = vals[4].value_in(units.pc)
        relax_time[0][i][0]  = vals[5].value_in(units.Myr)

    line_x = np.empty((len(IMBH_tracker), col_len, 1))
    line_y = np.empty((len(IMBH_tracker), col_len, 1))
    line_z = np.empty((len(IMBH_tracker), col_len, 1))
    tdyn   = np.empty((len(IMBH_tracker), col_len, 1))

    colours = colour_picker()

    for i in range(len(IMBH_tracker)):
        tIMBH_tracker = IMBH_tracker.iloc[i]
        tIMBH_tracker = tIMBH_tracker.replace(np.NaN, "[Np.NaN, [np.NaN, np.NaN, np.NaN]")
        for j in range(col_len):
            coords = tIMBH_tracker.iloc[j+1][1]
            tdynval = tIMBH_tracker.iloc[j+1][4]
            if len(coords) == 1:
                pass
            else:
                line_x[i][j][0] = coords[0].value_in(units.pc)
                line_y[i][j][0] = coords[1].value_in(units.pc)
                line_z[i][j][0] = coords[2].value_in(units.pc)
                tdyn[i][j][0] = tdynval.value_in(units.Myr)

    ejected_x, ejected_y, ejected_z = ejected_extract(IMBH_tracker, ejec_parti, col_len)
    ejected_x[:][abs(ejected_x[:]) < 10**-5] = np.NaN
    ejected_y[:][abs(ejected_y[:]) < 10**-5] = np.NaN
    ejected_z[:][abs(ejected_z[:]) < 10**-5] = np.NaN
    ejected_x[:][abs(ejected_x[:]) > 10**5] = np.NaN
    ejected_y[:][abs(ejected_y[:]) > 10**5] = np.NaN
    ejected_z[:][abs(ejected_z[:]) > 10**5] = np.NaN

    ejected_dist = np.sqrt((ejected_x-com_x)**2+(ejected_y-com_y)**2+(ejected_z-com_z)**2)

    line_x[:][abs(line_x[:]) < 10**-5] = np.NaN
    line_y[:][abs(line_y[:]) < 10**-5] = np.NaN
    line_z[:][abs(line_z[:]) < 10**-5] = np.NaN
    line_x[:][abs(line_x[:]) > 1e8] = np.NaN
    line_y[:][abs(line_y[:]) > 1e8] = np.NaN
    line_z[:][abs(line_z[:]) > 1e8] = np.NaN 

    tdyn[:][tdyn[:] < 10**-30] = np.NaN
    tdyn[:][tdyn[:] > 10**10] = np.NaN

    fig = plt.figure(figsize=(12.5, 15))
    ax1 = fig.add_subplot(321)
    ax2 = fig.add_subplot(322)
    ax3 = fig.add_subplot(323)
    ax4 = fig.add_subplot(324)
    ax5 = fig.add_subplot(325)
    ax6 = fig.add_subplot(326)
    
    ax1.set_title('Cluster Center of Mass')
    ax3.set_title('Ejected IMBH Focus')
    ax4.set_title('Complete System')
    ax5.set_title('Distance Evolution of IMBH Particles')
    ax6.set_title('Cluster Relaxation Timescale TO FIX')

    for ax_ in [ax1, ax2, ax3, ax4]:
        ax_.set_xlabel(r'$x$-Coordinate [pc]')
        ax_.set_ylabel(r'$y$-Coordinate [pc]')

    ax5.set_xlabel(r'Time [Myr]')
    ax5.set_ylabel(r'Distance [pc]')
    ax6.set_xlabel(r'Time [Myr]')
    ax6.set_ylabel(r'Relaxation Time [Myr]')

    for ax_ in [ax1, ax2, ax3, ax4]:
        ax_.yaxis.set_ticks_position('both')
        ax_.xaxis.set_ticks_position('both')
        ax_.xaxis.set_minor_locator(mtick.AutoMinorLocator())
        ax_.yaxis.set_minor_locator(mtick.AutoMinorLocator())
        ax_.tick_params(axis="y", which = 'both', direction="in")
        ax_.tick_params(axis="x", which = 'both', direction="in")       

    iter = -1
    for i in range(len(IMBH_tracker)):
        iter += 1
        if i == 0:
            adapt_c = 'black'
            ax4.scatter((line_x[i][:-1]-line_x[0][:-1]), (line_y[i][:-1]-line_y[0][:-1]), c = adapt_c, zorder = 1, s = 21)
            ax4.scatter((line_x[i][-1]-line_x[0][-1]), (line_y[i][-1]-line_y[0][-1]), c = adapt_c, s = 250, zorder = 3)
            ax4.scatter((line_x[i][0]-line_x[0][0]), (line_y[i][0]-line_y[0][0]), alpha = 0.7, c = adapt_c, s = 250, zorder = 2)
        if i == 1:
            ax4.plot((line_x[i][:-1]-line_x[0][:-1]), (line_y[i][:-1]-line_y[0][:-1]), c = adapt_c, ls = '--', zorder = 1)
        if iter > len(colours):
            iter = 0
        if i > 1:
            xval = (line_x[i][1:]-line_x[1][1:])
            yval = (line_y[i][1:]-line_y[1][1:])

            for ax_ in [ax1, ax2]:
                ax_.scatter(xval, yval, c = colours[iter-2], zorder = 1, s = 5)
                ax_.scatter(line_x[i][-1]-line_x[1][-1], line_y[i][-1]-line_y[1][-1], 
                            c = colours[iter-2], edgecolors = 'black', s = 50, zorder = 3)
            ax1.scatter(line_x[i][1]-line_x[1][1], line_y[i][1]-line_y[1][1], 
                        alpha = 0.7, c = colours[iter-2], edgecolors = 'black', s = 30, zorder = 2)

            ax3.scatter(line_x[i][:-1]-ejected_x[0][:-1], line_y[i][:-1]-ejected_y[0][:-1], 
                        c = colours[iter-2], s = 5, zorder = 1)
            ax3.scatter(line_x[i][1]-ejected_x[0][1], line_y[i][1]-ejected_y[0][1], 
                        alpha = 0.7,c = colours[iter-2], edgecolors = 'black', s = 30, zorder=2)
            ax3.scatter(line_x[i][-2]-ejected_x[0][-2], line_y[i][-2]-ejected_y[0][-2], 
                        c = colours[iter-2], edgecolors = 'black', s = 50, zorder=3)

            ax4.scatter((line_x[i][:-1]-line_x[0][:-1]), (line_y[i][:-1]-line_y[0][:-1]), c = colours[iter-2], zorder = 1, s = 1)
            ax4.scatter((line_x[i][-1]-line_x[0][-1]), (line_y[i][-1]-line_y[0][-1]), c = colours[iter-2], edgecolors = 'black', s = 50, zorder = 3)
            ax4.scatter((line_x[i][1]-line_x[0][1]), (line_y[i][1]-line_y[0][1]), alpha = 0.7, c = colours[iter-2], edgecolors = 'black', s = 50, zorder = 2)

    cluster_rad1 = plt.Circle((0,0), gc_code.gc_rad.value_in(units.pc), edgecolor = 'black', fill=False, ls = ':', label = 'Cluster Radius')
    cluster_rad2 = plt.Circle((0,0), gc_code.gc_rad.value_in(units.pc), edgecolor = 'black', fill=False, ls = ':', label = 'Cluster Radius')
    tidal_rad   = plt.Circle((0,0), rtide_array[0][1], edgecolor = 'black', fill=False, ls = '-.', label = 'Tidal Radius')
    ax1.add_patch(cluster_rad1)
    ax2.add_patch(tidal_rad)
    ax2.add_patch(cluster_rad2)
    ax1.legend(loc = 1)
    ax2.legend(loc = 1)

    ax5.plot(time[0][3:], rtide_array[0][3:], color = 'black',  label = r'$r_{tidal}$')
    ax5.plot(time[0][3:], LG25_array[0][3:],  color = 'red',   label = r'$r_{25,L}$')
    ax5.plot(time[0][3:], LG75_array[0][3:],  color = 'blue',  label = r'$r_{75,L}$')
    ax5.plot(time[0][3:], ejected_dist[0][3:], color = 'orange', label = 'Ejected Particle')
    ax5.legend()

    ax6.plot(time[0][3:], relax_time[0][3:], color = 'black',  label = r'$r_{tidal}$', linestyle = ":")
    ax6.legend()

    ax1.set_xlim(-1.25*gc_code.gc_rad.value_in(units.pc), 1.25*gc_code.gc_rad.value_in(units.pc))
    ax1.set_ylim(-1.25*gc_code.gc_rad.value_in(units.pc), 1.25*gc_code.gc_rad.value_in(units.pc)) 
    for ax_ in [ax2, ax3]:
        ax_.set_xlim(-1.3*rtide_array[0][0], 1.3*rtide_array[0][0])
        ax_.set_ylim(-1.3*rtide_array[0][0], 1.3*rtide_array[0][0])
    ax4.set_xlim(-1.05*init_dist.value_in(units.pc), 1.05*init_dist.value_in(units.pc))
    ax4.set_ylim(-1.05*init_dist.value_in(units.pc), 1.05*init_dist.value_in(units.pc)) 

    ax5.xaxis.set_major_formatter(mtick.FormatStrFormatter('%0.3f'))
    ax5.yaxis.set_major_formatter(mtick.FormatStrFormatter('%0.2f'))
    ax5.set_xlim(time[0][3], time[0][-1])
    ax5.set_ylim(0, 2*max(rtide_array[0]))
    plt.savefig('figures/spatial_tracker'+str(count)+'.pdf', dpi=300, bbox_inches='tight')
    plt.clf()
    plt.close()

    return

class stability_time(object):
    def __init__(self, dirC = 'data/chaotic_simulation/*'):
        self.chaos_ini_parti_data, self.chaos_fin_parti_data, self.chaos_number_mergers, self.chaos_cumulative_mm, self.chaos_simulated_end, self.chaos_ejected_parti, \
        self.chaos_stab_time_data, self.chaos_init_dist_data, self.chaos_cluster_radius, self.chaos_init_mass_data, self.chaos_inj_mass_data, \
        self.chaos_eje_mass_data, self.chaos_reltime_data = steadytime_extractor(dirC)  

    def stab_plotter_logistics(self, ax, pop, xints, ydata):
        """
        Function to set up the various stability time plotters
        
        Inputs:
        ax:     The axis
        pop:    The population, N, of the data
        xints:  The integer spacing for which to place the x-labels
        ydata:  The y-values for the plot
        """

        ax.yaxis.set_ticks_position('both')
        ax.xaxis.set_ticks_position('both')
        ax.xaxis.set_minor_locator(mtick.AutoMinorLocator())
        ax.yaxis.set_minor_locator(mtick.AutoMinorLocator())
        ax.tick_params(axis="y", which = 'both', direction="in")
        ax.tick_params(axis="x", which = 'both', direction="in")    
        ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%0.3f'))
        ax.xaxis.labelpad = 20
        plt.xticks(xints)
        plt.xlim(2.5,max(pop)+1)
        plt.ylim(0, 1.2*max(ydata))
        plt.ylabel(r'Ejection Time [Myr]')
        plt.xlabel(r'Number of IMBH [$N$]')
        plt.title(r'Chaotic Black Hole Population vs. Stability Time')

    def steadytime_distdep_plotter(self):
        """
        Function to plot the steady time with various lines corresponding to different
        SMBH Distances:
        """
        in_mass = np.unique(self.chaos_init_mass_data[:,0], axis = 0) #Get the unique initial mass and distances array
        in_dist = np.unique(self.chaos_init_dist_data)

        for mass_ in in_mass:           #For every initial mass, we will plot a separate graph showing different coloured data points
            fig, ax = plt.subplots()    #depending on their distance to the central SMBH
            colourcycler = cycle(colour_picker())
            tot_pop = [ ]
            y_max = []  
            iter = -1

            for dist_ in in_dist:
                iter += 1
                init_mass_idx = np.where((self.chaos_init_mass_data == mass_).all(1))[0] #Find indices where data files correspond to the correct initial masses
                fin_parti = self.chaos_fin_parti_data[init_mass_idx]                     #Filter the needed data based on the data files satisfying condition
                stab_time = self.chaos_stab_time_data[init_mass_idx]
                init_dist = self.chaos_init_dist_data[init_mass_idx]
                dist_arrays = np.where(init_dist == dist_)                     #Find indices corresponding to the correct distances. Here we split

                colours = next(colourcycler)
                N_parti_avg = [ ]
                
                fin_parti = fin_parti[dist_arrays]                             #The data relative to their distances.
                stab_time = stab_time[dist_arrays]

                pop_size, pop_samples = np.unique(fin_parti, return_counts=True)  #Count the number of unique final populations
                if len(pop_size) == 0:
                    pass
                else:
                    pop_id = np.argwhere(pop_size > 2)                                #Only look at samples with N>2
                    pop_size = pop_size[pop_id]

                    pop_samples = pop_samples[pop_id]
                    tot_pop.append(max(pop_size))

                    for pop_, samp_ in zip(pop_size, pop_samples):
                        N_parti = np.argwhere(fin_parti == pop_)                      #For every N, we will gather their indices where they are in the pickle file
                        N_parti_avg.append(np.mean(stab_time[N_parti]))               #This enables us to compute their average between the whole data set                     
                    
                    y_max.append(max(N_parti_avg))
                    ax.scatter(pop_size, N_parti_avg, color = colours, edgecolor = 'black', label = r'$r_{SMBH} =$'+str(dist_)+' pc')
                    for j, xpos in enumerate(pop_size):
                        if iter == 0:
                            ax.text(xpos, -2500, 'Simulations\n'+'Set '+str(iter)+': '+str('{:.0f}'.format(pop_samples[j][0])), fontsize = 'xx-small', ha = 'center' )
                        else:
                            ax.text(xpos, -2500*(1+0.4*iter),  'Set '+str(iter)+': '+str(pop_samples[j][0]), fontsize = 'xx-small', ha = 'center' )
                xints = [i for i in range(1+int(max(tot_pop)))]

            order = int('{:.0f}'.format(np.log10(mass_)))
            ax.text(0.34*max(tot_pop), 1.125*max(y_max), r'$m_{i} \in$ '+str('{:.0f}'.format(mass_/10**order))+r'$\times 10^{}$'.format(order)+r'$M_\odot$')
            self.stab_plotter_logistics(ax, tot_pop, xints, y_max)
            plt.legend()
            plt.savefig('figures/chaotic_stab_time_equal_mass_'+str(mass_)+'_mean.pdf', dpi = 300, bbox_inches='tight')

        for mass_ in in_mass:      #For every initial mass and distance we will plot a separate graph. This time it includes std spread
            for dist_ in in_dist:  #For three different distances, every mass has 3 unique plots with errors.
                y_max = []
                
                plt.clf()
                fig, ax = plt.subplots()
                
                init_mass_idx = np.where((self.chaos_init_mass_data == mass_).all(1))[0]     #Make sure to use only data who has corresponding mass
                fin_parti = self.chaos_fin_parti_data[init_mass_idx]                         #Filter data through indices satisfying the mass requirement
                stab_time = self.chaos_stab_time_data[init_mass_idx]
                init_dist = self.chaos_init_dist_data[init_mass_idx]
                dist_arrays = np.where(init_dist == dist_)                        #Find indices corresponding to the correct distances. Here we split

                colours = next(colourcycler)
                N_parti_avg = [ ]
                std = [ ]
                
                fin_parti = fin_parti[dist_arrays]                                #Filtering the filtered data with the correct distances
                stab_time = stab_time[dist_arrays]

                pop_size, pop_samples = np.unique(fin_parti, return_counts=True)  #Count the number of unique final populations
                
                if len(pop_size) == 0:
                    pass
                else:
                    pop_id = np.argwhere(pop_size > 2)
                    pop_size = pop_size[pop_id]
                    pop_samples = pop_samples[pop_id]

                    for pop_, samp_ in zip(pop_size, pop_samples):                    #For the unique populations, compute their average stab time
                        N_parti = np.argwhere(fin_parti == pop_)                      #and average spread in data.
                        N_parti_avg.append(np.mean(stab_time[N_parti]))
                        std.append(np.std(stab_time[N_parti]))

                    y_max.append(max(np.add(N_parti_avg, std)))
                    ax.errorbar(pop_size, N_parti_avg, color = 'black', yerr=std, fmt = 'o')
                    ax.scatter(pop_size, np.add(N_parti_avg, std), marker = '_', color = 'black')
                    ax.scatter(pop_size, np.subtract(N_parti_avg, std), marker = '_', color = 'black')
                    ax.scatter(pop_size, N_parti_avg, color = 'black')

                    for j, xpos in enumerate(pop_size):
                        ax.text(xpos, -0.115*max(np.add(N_parti_avg, std)), 'Simulations\n'+str(pop_samples[j][0]), fontsize = 'xx-small', ha = 'center' )
                    
                    order = int('{:.0f}'.format(np.log10(mass_)))
                    xints = [i for i in range(1+int(max(tot_pop)))]
                    ax.text(0.34*max(tot_pop), 1.05*(max(np.add(N_parti_avg, std))), r'$m_{IMBH}=$'+str('{:.0f}'.format(mass_/10**order))+r'$\times 10^{}$'.format(order)+r' $M_\odot$'+'\n'+r'$r_{SMBH}=$'+str(dist_)+' pc')
                    self.stab_plotter_logistics(ax, tot_pop, xints, np.add(N_parti_avg, std))
                    plt.savefig('figures/chaotic_stab_time_equal_mass_'+str(mass_)+'_err_dist_'+str(dist_)+'.pdf', dpi = 300, bbox_inches='tight')

    def steadytime_massdep_plotter(self):

        in_dist = np.unique(self.chaos_init_dist_data)
        in_mass = np.unique(self.chaos_init_mass_data, axis=0)

        for dist_ in in_dist:
            fig, ax = plt.subplots()
            tot_pop = []
            y_max = []     
            iter = -1   
            init_dist_idx = np.where((self.chaos_init_dist_data == dist_))
            coloursycler = cycle(colour_picker())

            for mass_ in in_mass:
                iter += 1

                init_mass = self.chaos_init_mass_data[init_dist_idx]
                init_mass = self.chaos_init_mass_data[init_dist_idx[0]]
                fin_parti = self.chaos_fin_parti_data[init_dist_idx[0]]
                stab_time = self.chaos_stab_time_data[init_dist_idx[0]]

                colours = next(coloursycler)
                N_parti_avg = [ ]
                std = [ ]
                mass_arrays = np.where((init_mass == mass_).all(1))[0]  #Indices with the correct mass column
                fin_parti = fin_parti[mass_arrays]
                stab_time = stab_time[mass_arrays]
                pop_size, pop_samples = np.unique(fin_parti, return_counts=True)

                if len(pop_size) == 0:
                    pass
                else:
                    pop_id = np.argwhere(pop_size > 2)
                    pop_size = pop_size[pop_id]
                    pop_samples = pop_samples[pop_id]
                    tot_pop.append(max(pop_size))

                    for pop_, samp_ in zip(pop_size, pop_samples):
                        N_parti = np.argwhere(fin_parti == pop_)
                        N_parti_avg.append(np.mean(stab_time[N_parti]))
                        std.append(np.std(stab_time[N_parti]))

                    y_max.append(max(N_parti_avg))
                    ax.scatter(pop_size, N_parti_avg, color = colours, edgecolor = 'black',
                            label = r'$m_{i} \in$ ['+str(mass_[0])+', '+str(mass_[1])+r'] $M_\odot$')
                    
                    for j, xpos in enumerate(pop_size):
                        if iter == 0:
                            ax.text(xpos, -0.115*max(N_parti_avg), 'Simulations\n'+'Set '+str(iter)+': '+str('{:.0f}'.format(pop_samples[j][0])), fontsize = 'xx-small', ha = 'center' )
                        else:
                            ax.text(xpos, -0.115*max(N_parti_avg)*(1+0.6*iter), 'Set '+str(iter)+': '+str(pop_samples[j][0]), fontsize = 'xx-small', ha = 'center' )
                xints = [i for i in range(1+int(max(tot_pop)))]

            ax.text(0.34*max(tot_pop), 1.125*max(y_max), r'$r_{SMBH}=$'+str(dist_)+' pc')
            self.stab_plotter_logistics(ax, tot_pop, xints, y_max)
            plt.legend()
            plt.savefig('figures/chaotic_stab_time_equal_dist_'+str(dist_)+'_mean.pdf', dpi = 300, bbox_inches='tight')

        for dist_ in in_dist:
            for mass_ in in_mass: 
                y_max = []

                plt.clf()
                fig, ax = plt.subplots()

                init_dist_idx = np.where((self.chaos_init_dist_data == dist_))
                init_mass = self.chaos_init_mass_data[init_dist_idx]
                fin_parti = self.chaos_fin_parti_data[init_dist_idx]
                stab_time = self.chaos_stab_time_data[init_dist_idx]

                tot_pop = [ ]
                colours = next(coloursycler)
                N_parti_avg = [ ]
                std = [ ]
                mass_arrays = np.where((init_mass == mass_).all(1))[0]  #Indices with the correct mass column

                fin_parti = fin_parti[mass_arrays]
                stab_time = stab_time[mass_arrays]
                pop_size, pop_samples = np.unique(fin_parti, return_counts=True)
                pop_id = np.argwhere(pop_size > 2)
                pop_size = pop_size[pop_id]

                if len(pop_size) == 0:
                    pass
                else:
                    pop_samples = pop_samples[pop_id]
                    tot_pop.append(max(pop_size))

                    for pop_, samp_ in zip(pop_size, pop_samples):
                        N_parti = np.argwhere(fin_parti == pop_)
                        N_parti_avg.append(np.mean(stab_time[N_parti]))
                        std.append(np.std(stab_time[N_parti]))

                    y_max.append(max(np.add(N_parti_avg, std)))
                    ax.errorbar(pop_size, N_parti_avg, color = 'black', yerr=std, fmt = 'o')
                    ax.scatter(pop_size, np.add(N_parti_avg, std), marker = '_', color = 'black')
                    ax.scatter(pop_size, np.subtract(N_parti_avg, std), marker = '_', color = 'black')
                    ax.scatter(pop_size, N_parti_avg, color = 'black')

                    for j, xpos in enumerate(pop_size):
                        ax.text(xpos, -0.13*max(np.add(N_parti_avg, std)), 'Simulations\n'+str('{:.0f}'.format(pop_samples[j][0])), fontsize = 'xx-small', ha = 'center' )
                    xints = [i for i in range(1+int(max(tot_pop)))]

                    ax.text(0.34*max(tot_pop), 1.05*(max(np.add(N_parti_avg, std))), r'$r_{SMBH}=$'+str(dist_)+' pc\n'+r'$m_{i} \in$ ['+str(mass_[0])+', '+str(mass_[1])+r'] $M_\odot$')
                    self.stab_plotter_logistics(ax, tot_pop, xints, np.add(N_parti_avg, std))
                    plt.savefig('figures/chaotic_stab_time_equal_dist_'+str(dist_)+'_std_mass_'+str(mass_)+'.pdf', dpi = 300, bbox_inches='tight')

#gc_code = globular_cluster()
#spatial_plotter(1.15*gc_code.gc_dist)
#energy_plotter()
#animator(1.5 | units.parsec)

st = stability_time()
st.steadytime_massdep_plotter()
st.steadytime_distdep_plotter()
