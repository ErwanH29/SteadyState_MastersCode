from amuse.lab import *
from parti_initialiser import *
from file_logistics import *
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import matplotlib.animation as animation
import numpy as np

class plotter_setup():
    def val_filter(self, arr):
        arr[:][abs(arr[:]) < 10**-5] = np.NaN
        arr[:][abs(arr[:]) > 10**5] = np.NaN
        return arr
         
    def tickers(self, ax):
        """
        Function to give outlay for the axis
        """

        ax.yaxis.set_ticks_position('both')
        ax.xaxis.set_ticks_position('both')
        ax.xaxis.set_minor_locator(mtick.AutoMinorLocator())
        ax.yaxis.set_minor_locator(mtick.AutoMinorLocator())
        ax.tick_params(axis="y", which = 'both', direction="in")
        ax.tick_params(axis="x", which = 'both', direction="in")    
        return ax

    def xy_pc(self, ax):
        """
        Function which gives the x,y labels if they are in [pc]
        """

        ax.set_xlabel(r'$x$-Coordinate [pc]')
        ax.set_ylabel(r'$x$-Coordinate [pc]')
        return ax

def colour_picker():
    """
    Colour chooser for the spatial plots and animation
    """

    colors = ['red', 'blue', 'orange', 'purple', 'salmon', 'slateblue', 
              'gold', 'darkviolet', 'cornflowerblue',  'cyan',
              'lightskyblue', 'magenta',  'dodgerblue']

    return colors

def animator(init_dist, int_string):
    """
    Function to produce animations. WARNING: SLOW
    
    Inputs:
    init_dist: The initial distance of the cluster to the central SMBH
    """

    print('!!!!! You have chosen to animate. This will take a while !!!!!')
    
    plot_ini = plotter_setup()
    count = file_counter(int_string)

    Lag_tracker, col_len = file_opener('data/'+str(int_string)+'/lagrangians/*')
    IMBH_tracker, col_len = file_opener('data/'+str(int_string)+'/particle_trajectory/*')
    energy_tracker, col_len = file_opener('data/'+str(int_string)+'/energy/*')

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
            coords = tIMBH_tracker.iloc[j+1][2]
            if len(coords) == 1:
                pass
            else:
                line_x[i][j][0] = coords[0].value_in(units.pc)
                line_y[i][j][0] = coords[1].value_in(units.pc)
                line_z[i][j][0] = coords[2].value_in(units.pc)

    temp_coord_array = [line_x, line_y, line_z]
    for arr_ in temp_coord_array:
        plot_ini.val_filter(arr_)

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

    plt.clf()
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

                ax3D2.scatter((line_x[i][max(0,skip_zeroth-20):skip_zeroth+1]-line_x[1][max(0,skip_zeroth-20):skip_zeroth+1]),
                              (line_y[i][max(0,skip_zeroth-20):skip_zeroth+1]-line_y[1][max(0,skip_zeroth-20):skip_zeroth+1]),
                              (line_z[i][max(0,skip_zeroth-20):skip_zeroth+1]-line_z[1][max(0,skip_zeroth-20):skip_zeroth+1]), 
                               s =1, c = colours[iter-2], lw = 2)
                ax3D2.scatter((line_x[i][skip_zeroth]-line_x[1][skip_zeroth]), 
                              (line_y[i][skip_zeroth]-line_y[1][skip_zeroth]), 
                              (line_z[i][skip_zeroth]-line_z[1][skip_zeroth]),
                               c = colours[iter-2], edgecolors = 'black', s = 40)

        for ax_ in [ax3D1, ax3D2]:
            ax_.xaxis.set_major_formatter(mtick.FormatStrFormatter('%0.2f'))
            ax_.yaxis.set_major_formatter(mtick.FormatStrFormatter('%0.2f'))
            ax_.zaxis.set_major_formatter(mtick.FormatStrFormatter('%0.2f'))
            
        ax3D1.scatter(0,0, color = 'black', s = 100 )
        ax3D1.set_xlabel(r'$x$-Coordinate [pc]')
        plot_ini.xy_pc(ax3D2)

    anim3D = animation.FuncAnimation(fig, animate_3D, interval = 100, frames=(col_len-1))
    anim3D.save('figures/animation3D'+str(count)+'.gif', writer='pillow')

    fig = plt.figure(figsize=(12.5, 8))
    ax1 = fig.add_subplot(221)
    ax2 = fig.add_subplot(222)
    ax3 = fig.add_subplot(223)
    ax4 = fig.add_subplot(224)

    def animate(col_len):
        skip_zeroth = col_len + 1
        if skip_zeroth %100 == 0:
            print('Current frame: ', skip_zeroth)
            
        for ax_ in [ax1, ax2, ax3, ax4]:
            plot_ini.tickers(ax_) 
        
        for i in range(len(IMBH_tracker)):
            ax1.scatter(line_x[i][max(0,skip_zeroth-20):skip_zeroth+1]-line_x[1][max(0,skip_zeroth-20):skip_zeroth+1],
                        line_y[i][max(0,skip_zeroth-20):skip_zeroth+1]-line_y[1][max(0,skip_zeroth-20):skip_zeroth+1],
                        c = colours[i], lw = 2, alpha = 0.7, s = 1)
            ax1.scatter(line_x[i][skip_zeroth]-line_x[1][skip_zeroth], line_y[i][skip_zeroth]-line_t[1][skip_zeroth],
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
            plot_ini.xy_pc(ax_)

        for ax_ in [ax2, ax4]:
            ax_.set_xlabel(r'Time [Myr]')

        ax1.set_title(str("{:.0f}".format(time[skip_zeroth][0])+" Myr"), loc = 'left')   
        ax1.yaxis.set_major_formatter(mtick.FormatStrFormatter('%0.3f'))
        ax1.xaxis.set_major_formatter(mtick.FormatStrFormatter('%0.3f'))
        ax1.set_xlim(min(line_x[i][max(0,skip_zeroth-80):skip_zeroth+1]-line_x[1][max(0,skip_zeroth-80):skip_zeroth+1]),
                     max(line_x[i][max(0,skip_zeroth-80):skip_zeroth+1]-line_x[1][max(0,skip_zeroth-80):skip_zeroth+1]))
        ax1.set_ylim(min(line_y[i][max(0,skip_zeroth-80):skip_zeroth+1]-line_y[1][max(0,skip_zeroth-80):skip_zeroth+1]),
                     max(line_y[i][max(0,skip_zeroth-80):skip_zeroth+1]-line_y[1][max(0,skip_zeroth-80):skip_zeroth+1]))

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
        
def energy_plotter(int_string):
    """
    Function to plot the energy evolution of the system
    """

    plot_ini = plotter_setup()
    count = file_counter(int_string)
    energy_tracker, col_len = file_opener('data/'+str(int_string)+'/energy/*')
    #IMBH_energy_tracker, col_len = file_opener('data/'+str(int_string)+'/particle_energies/*')

    time = np.empty((1, col_len, 1))
    Et_array = np.empty((1, col_len, 1))
    dE_array = np.empty((1, col_len, 1))
    dEs_array = np.empty((1, col_len, 1))
    IMBH_birth = np.empty((1, col_len, 1))
    collisions = np.empty((1, col_len, 1)) 
    merger_mass = np.empty((1, col_len, 1)) 

    """BE_array = np.empty((1, col_len, 1))
    KE_array = np.empty((1, col_len, 1))
    TotE_array = np.empty((1, col_len, 1))"""

    for i in range(col_len):
        vals = energy_tracker.iloc[i]
        time[0][i][0] = vals[0].value_in(units.Myr)
        Et_array[0][i][0] = vals[1].value_in(units.J)
        dE_array[0][i][0] = vals[2]
        dEs_array[0][i][0] = vals[3]
        IMBH_birth[0][i][0] = vals[4].value_in(units.Myr)
        collisions[0][i][0] = vals[5].value_in(units.Myr)
        merger_mass[0][i][0] = vals[6].value_in(units.MSun)

        """vals = IMBH_energy_tracker.iloc[i]
        BE_array[0][i][0] = vals[1].value_in(units.J)
        KE_array[0][i][0] = vals[2].value_in(units.J)
        TotE_array[0][i][0] = vals[3].value_in(units.J)"""    

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
        plot_ini.tickers(ax_)  
        ax_.set_xlabel(r'Time [Myr]')
        ax_.set_yscale('log')
    ax1.set_xlim(time[0][5], time[0][-1])
    ax1.set_ylabel(r'$\frac{|E(t)-E_0|}{|E_0|}$')
    ax2.set_ylabel(r'Energy [J]')
    #ax2.set_ylim(0.1*min(abs(BE_array[0][:])), 5*max(KE_array[0][:]))

    ax1.plot(time[0][:], dE_array[0][:], color = 'black', zorder = 1)
    if len(coll_id) > 0:
        dE_coll = dE_array[0][coll_id[:,1]-1]
        collisions = collisions[0][coll_id[:,1]]
        merger_mass = merger_mass[0][coll_id[:,1]]
        color_axes = ax1.scatter(collisions, dE_coll, c = merger_mass, zorder=3)
        plt.colorbar(color_axes, ax=ax1, label = r'Merger Mass [$M_{\odot}$]')
    """ax2.plot(time[0][5:], abs(BE_array[0][5:]), color = 'blue', label = 'Potential Energy (abs)', zorder = 1)
    ax2.plot(time[0][5:], KE_array[0][5:], color = 'red', label = 'Kinetic Energy', zorder = 2)
    ax2.plot(time[0][5:], abs(TotE_array[0][5:]), color = 'black', label = 'Total Energy (abs)', linestyle = '--', zorder = 3)"""
    ax2.legend()

    if len(app_id) > 0:
        dE_app = dE_array[0][app_id[:,0]+1]
        IMBH_app = IMBH_birth[0][app_id[:,0]+1]
        #totE_app = abs(TotE_array[0][app_id[:,0]+1])

        ax1.scatter(IMBH_app, dE_app, marker = 'x', color = 'red', zorder = 2)
        ax2.scatter(IMBH_app, totE_app, marker = 'x', color = 'red', zorder = 4)
    plt.savefig('figures/energy_tracker'+str(count)+'.pdf', dpi=300, bbox_inches='tight')
    plt.clf()
    plt.close()

    return

def spatial_plotter(init_dist, int_string):
    """
    Function to plot the evolution of the system
    """

    plot_ini = plotter_setup()
    gc_code = globular_cluster()
    count = file_counter(int_string)
    IMBH_tracker, col_len = file_opener('data/'+str(int_string)+'/particle_trajectory/*')
    ejec_parti, col_len = file_opener('data/'+str(int_string)+'/no_addition/chaotic_simulation/*')
    #Lag_tracker, col_len = file_opener('data/'+str(int_string)+'/lagrangians/*')

    time = np.empty((1, col_len, 1))
    """LG25_array  = np.empty((1, col_len, 1))
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
        relax_time[0][i][0]  = vals[5].value_in(units.Myr)"""

    col_len = np.shape(IMBH_tracker)[1]-1

    line_x = np.empty((len(IMBH_tracker), col_len, 1))
    line_y = np.empty((len(IMBH_tracker), col_len, 1))
    line_z = np.empty((len(IMBH_tracker), col_len, 1))
    tdyn   = np.empty((len(IMBH_tracker), col_len, 1))

    colours = colour_picker()

    for i in range(len(IMBH_tracker)):
        tIMBH_tracker = IMBH_tracker.iloc[i]
        tIMBH_tracker = tIMBH_tracker.replace(np.NaN, "[Np.NaN, [np.NaN, np.NaN, np.NaN]")
        for j in range(col_len):
            coords = tIMBH_tracker.iloc[j+1][2]
            tdynval = tIMBH_tracker.iloc[j+1][6]
            if len(coords) == 1:
                pass
            else:
                line_x[i][j][0] = coords[0].value_in(units.pc)
                line_y[i][j][0] = coords[1].value_in(units.pc)
                line_z[i][j][0] = coords[2].value_in(units.pc)
                #tdyn[i][j][0] = tdynval.value_in(units.Myr)

    ejected_x, ejected_y, ejected_z, evx, evy, evz = ejected_extract(IMBH_tracker, ejec_parti, col_len)
    for arr_ in [ejected_x, ejected_y, ejected_z]:
        plot_ini.val_filter(arr_)
    ejected_dist = np.sqrt((ejected_x-line_x[1])**2+(ejected_y-line_y[1])**2+(ejected_z-line_z[1])**2)

    for arr_ in [line_x, line_y, line_z]:
        plot_ini.val_filter(arr_)

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
        plot_ini.xy_pc(ax_)

    ax5.set_xlabel(r'Time [Myr]')
    ax5.set_ylabel(r'Distance [pc]')
    ax6.set_xlabel(r'Time [Myr]')
    ax6.set_ylabel(r'Relaxation Time [Myr]')

    for ax_ in [ax1, ax2, ax3, ax4]:
        plot_ini.tickers(ax_)      
        
    iter = -1
    for i in range(len(IMBH_tracker)):
        iter += 1
        if iter > len(colours):
            iter = 0
        if i == 0:
            adapt_c = 'black'
            ax4.scatter((line_x[i][:-1]-line_x[0][:-1]), (line_y[i][:-1]-line_y[0][:-1]), c = adapt_c, zorder = 1, s = 21)
            ax4.scatter((line_x[i][-1]-line_x[0][-1]), (line_y[i][-1]-line_y[0][-1]), c = adapt_c, s = 250, zorder = 3)
            ax4.scatter((line_x[i][0]-line_x[0][0]), (line_y[i][0]-line_y[0][0]), alpha = 0.7, c = adapt_c, s = 250, zorder = 2)
        #if i == 1:
         #   ax4.plot((line_x[i][:-1]-line_x[0][:-1]), (line_y[i][:-1]-line_y[0][:-1]), c = adapt_c, ls = '--', zorder = 1)
        
        else:
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
    #tidal_rad   = plt.Circle((0,0), rtide_array[0][1], edgecolor = 'black', fill=False, ls = '-.', label = 'Tidal Radius')
    ax1.add_patch(cluster_rad1)
    #ax2.add_patch(tidal_rad)
    ax2.add_patch(cluster_rad2)
    ax1.legend(loc = 1)
    ax2.legend(loc = 1)

    """ax5.plot(time[0][3:], rtide_array[0][3:], color = 'black',  label = r'$r_{tidal}$')
    ax5.plot(time[0][3:], LG25_array[0][3:],  color = 'red',   label = r'$r_{25,L}$')
    ax5.plot(time[0][3:], LG75_array[0][3:],  color = 'blue',  label = r'$r_{75,L}$')
    ax5.plot(time[0][4:], ejected_dist[0][3:], color = 'orange', label = 'Ejected Particle')
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
    ax5.set_ylim(0, 2*max(rtide_array[0]))"""
    plt.savefig('figures/spatial_tracker'+str(count)+'.pdf', dpi=300, bbox_inches='tight')
    plt.clf()
    plt.close()

    return

class vejec_mass(object):
    """
    Class to plot the indep_variable vs. ejection velocity plots
    """
    def __init__(self, int_string = 'Hermite'):
        """
        Extracts the required data
        """
        self.ejec_data = bulk_stat_extractor('data/'+str(int_string)+'/temp_plot/no_addition/chaotic_simulation/*')
        self.IMBH_tracker = bulk_stat_extractor('data/'+str(int_string)+'/temp_plot/particle_trajectory/*')
        self.data_entries = 1
        
        self.ejec_vx = np.empty((len(self.ejec_data)))
        self.ejec_vy = np.empty((len(self.ejec_data)))
        self.ejec_vz = np.empty((len(self.ejec_data)))
        self.tot_mass = np.empty((len(self.ejec_data)))
        self.tot_pop = np.empty((len(self.ejec_data)))
        self.surv_time = np.empty((len(self.ejec_data)))

        for i in range(len(self.ejec_data)):
            self.ex, self.ey, self.ez, self.ejec_vx[i], self.ejec_vy[i], self.ejec_vz[i] = ejected_extract(self.IMBH_tracker[i], 
                                                                                                           self.ejec_data[i], 
                                                                                                           self.data_entries)
            
            self.vals_df = self.ejec_data[i].iloc[0]
            self.mass_vals = np.sum(self.vals_df[9].value_in(units.MSun))
            self.tot_pop[i] = self.vals_df[0]
            self.tot_mass[i] = self.mass_vals

            self.surv_df = self.ejec_data[i].iloc[0]
            self.surv_time[i] = self.surv_df[4].value_in(units.Myr)

        self.ejec_vel = [np.sqrt(i**2+j**2+k**2) for i,j,k in zip(self.ejec_vx, self.ejec_vy, self.ejec_vz)]

    def plotter_func(self, xdata, ydata, cdata, clabel):
        """
        Function to plot the wanted data
        
        Inputs:
        xdata, ydata: The x and y variables wanting to plot
        cdata:        The data to be used as the colour code
        clabel:       Label describing the colour plot values
        """

        plot_ini = plotter_setup()
        fig, ax = plt.subplots()
        color_axes = ax.scatter(xdata, ydata, c = cdata)
        plt.colorbar(color_axes, ax=ax, label = clabel)
        plot_ini.tickers(ax)
        return ax

    def vejec_sysmass(self):
        """
        Function to plot how the total initial mass of the system influences the ejection velocity
        """

        in_mass = np.unique(self.tot_mass)
        avg_vel = np.empty((len(in_mass)))
        avg_surv = np.empty((len(in_mass)))

        iter = -1
        for mass_ in in_mass:
            iter += 1
            indices = np.where((self.tot_mass == mass_))[0]
            temp_evel = [self.ejec_vel[i] for i in indices]
            temp_surv = [self.surv_time[i] for i in indices]
            avg_vel[iter] = np.mean(temp_evel)
            avg_surv[iter] = np.mean(temp_surv)

        ax = self.plotter_func(in_mass, avg_vel, avg_surv, r'Ejection Time [Myr]')

        ax.set_xlabel(r'Total IMBH Mass [$\frac{M}{M_{\odot}}$]')
        ax.set_ylabel(r'Ejection Velocity [km/s]')
        #ax.set_ylim(0, 200)
        plt.show()
        plt.clf()

    def vejec_syspop(self):
        """
        Plot to show how the total population of the system influences the ejection velocity
        """
        
        ax = self.plotter_func(self.tot_pop, self.ejec_vel, (self.tot_mass), r'Total IMBH Mass [$\frac{M}{M_{\odot}}$]')
        ax.set_xlabel(r'IMBH Population [$N$]')
        ax.set_ylabel(r'Ejection Velocity [km/s]')
        #ax.set_ylim(0, 200)
        plt.show()
        plt.clf()

spatial_plotter(60 | units.parsec, 'Hermite')