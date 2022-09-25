from amuse.lab import *
from parti_initialiser import *
from file_logistics import *
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import matplotlib.animation as animation
import pickle as pkl
import numpy as np
import glob
from evol_func import file_counter

def colour_picker():
    """
    Colour chooser for the spatial plots and animation
    """

    colors = ['red', 'blue', 'orange', 'purple', 'salmon', 'slateblue', 
              'darkviolet', 'cornflowerblue', 'cyan', 'sienna', 
              'lightskyblue', 'magenta', 'gold', 'dodgerblue', 'red',
              'blue', 'red', 'orange', 'purple', 'salmon', 'slateblue', 
              'darkviolet', 'cornflowerblue', 'cyan', 'sienna', 
              'lightskyblue', 'magenta', 'gold', 'dodgerblue', 'red', 
              'blue', 'red', 'orange', 'purple', 'salmon', 'slateblue', 
              'darkviolet', 'cornflowerblue', 'cyan', 'sienna', 
              'lightskyblue', 'magenta', 'gold', 'dodgerblue', 'red', 
              'blue', 'red', 'orange', 'purple', 'salmon', 'slateblue', 
              'darkviolet', 'cornflowerblue', 'cyan', 'sienna', 
              'lightskyblue', 'magenta', 'gold', 'dodgerblue', 'red']
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

    Lag_tracker, col_len  = file_opener('data/lagrangians/*')
    com_tracker, col_len  = file_opener('data/center_of_mass/*')
    IMBH_tracker, col_len = file_opener('data/positions_IMBH/*')
    energy_tracker, col_len = file_opener('data/energy/*')

    time = np.empty((1, col_len, 1))
    dE_array  = np.empty((1, col_len, 1))

    for i in range(col_len):
        vals = energy_tracker.iloc[i]
        time[0][i][0] = vals[0].value_in(units.kyr)
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
                line_x[i][j][0] = coords[0].value_in(units.pc)
                line_y[i][j][0] = coords[1].value_in(units.pc)
                line_z[i][j][0] = coords[2].value_in(units.pc)

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

    totsys_lim = init_dist.value_in(units.pc)

    plt.ioff()
    fig = plt.figure()
    ax3D1 = fig.add_subplot(121, projection="3d")
    ax3D2 = fig.add_subplot(122, projection="3d")
    colours = colour_picker()

    def animate_3D(col_len):
        skip_zeroth = col_len + 1
        if skip_zeroth %100 == 0:
            print('Current 3D frame: ', skip_zeroth)

        ax3D1.clear()
        ax3D2.clear()
        ax3D1.set_title(str("{:.3f}".format(time[0][skip_zeroth][0])+" kyr"), loc = 'left') 
        ax3D1.set_xlim3d(-1.05*totsys_lim, 1.05*totsys_lim)
        ax3D1.set_ylim3d(-1.05*totsys_lim, 1.05*totsys_lim)
        ax3D1.set_zlim3d(-1.05*totsys_lim, 1.05*totsys_lim)

        for i in range(len(IMBH_tracker)):
            ax3D1.scatter(line_x[i][max(0,skip_zeroth-20):skip_zeroth+1], 
                          line_y[i][max(0,skip_zeroth-20):skip_zeroth+1], 
                          line_z[i][max(0,skip_zeroth-20):skip_zeroth+1], 
                          s =1, c = colours[i], lw = 2)
            ax3D1.scatter(line_x[i][skip_zeroth], line_y[i][skip_zeroth], line_z[i][skip_zeroth], 
                          c = colours[i], edgecolors = 'black', s = 50)

            ax3D2.scatter(line_x[i][max(0,skip_zeroth-20):skip_zeroth+1]-com_x[0][max(0,skip_zeroth-20):skip_zeroth+1],
                          line_y[i][max(0,skip_zeroth-20):skip_zeroth+1]-com_y[0][max(0,skip_zeroth-20):skip_zeroth+1],
                          line_z[i][max(0,skip_zeroth-20):skip_zeroth+1]-com_z[0][max(0,skip_zeroth-20):skip_zeroth+1], 
                          s =1, c = colours[i], lw = 2)
            ax3D2.scatter(line_x[i][skip_zeroth]-com_x[0][skip_zeroth], 
                          line_y[i][skip_zeroth]-com_y[0][skip_zeroth], 
                          line_z[i][skip_zeroth]-com_z[0][skip_zeroth],
                          c = colours[i], edgecolors = 'black', s = 40)
        ax3D1.scatter(0,0, color = 'black', s = 100 )
        ax3D1.set_xlabel(r'$x$-Coordinate [pc]')
        ax3D1.set_ylabel(r'$y$-Coordinate [pc]')

    anim3D = animation.FuncAnimation(fig, animate_3D, interval = 100, frames=(col_len-1))
    anim3D.save('figures/animation'+str(count)+'.gif', writer='pillow')

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
        skip_zeroth = col_len + 1
        if skip_zeroth %100 == 0:
            print('Current frame: ', skip_zeroth)

        ax1.clear()
        ax2.clear()
        ax3.clear()
        ax4.clear()
        
        colours = colour_picker()
        for i in range(len(IMBH_tracker)):
            ax1.plot(line_x[i][max(0,skip_zeroth-20):skip_zeroth+1]-com_x[0][max(0,skip_zeroth-20):skip_zeroth+1],
                     line_y[i][max(0,skip_zeroth-20):skip_zeroth+1]-com_y[0][max(0,skip_zeroth-20):skip_zeroth+1],
                     c = colours[i], lw = 2)
            ax1.scatter(line_x[i][skip_zeroth]-com_x[0][skip_zeroth], line_y[i][skip_zeroth]-com_y[0][skip_zeroth],
                        c = colours[i], edgecolors = 'black', s = 40)

            ax3.scatter(line_x[i][max(0,skip_zeroth-20):skip_zeroth+1], line_y[i][max(0,skip_zeroth-20):skip_zeroth+1], 
                        s =1, c = colours[i], lw = 2)
            ax3.scatter(line_x[i][skip_zeroth], line_y[i][skip_zeroth], c = colours[i], edgecolors = 'black', s = 50)

            ax4.plot(time[0][:skip_zeroth+1], LG25_array[0][:skip_zeroth+1], color = 'blue', label = r'$r_{25,L}$')
            ax4.plot(time[0][:skip_zeroth+1], LG75_array[0][:skip_zeroth+1], color = 'red',  label = r'$r_{75,L}$')

        ax3.scatter(0,0, color = 'black', s = 100 )
        ax2.plot(time[0][:skip_zeroth+1], abs(dE_array[0][:skip_zeroth+1]), color = 'black')

        ax1.set_title(str("{:.3f}".format(time[0][skip_zeroth][0])+" kyr"), loc = 'left')   
        ax1.set_xlabel(r'$x$-Coordinate [pc]')
        ax1.set_ylabel(r'$y$-Coordinate [pc]')
        ax1.yaxis.set_major_formatter(mtick.FormatStrFormatter('%0.2f'))
        ax1.xaxis.set_major_formatter(mtick.FormatStrFormatter('%0.2f'))

        ax2.set_xlabel(r'Time [kyr]')
        ax2.set_ylabel(r'$\frac{|E(t)-E_0|}{|E_0|}$')
        ax2.xaxis.set_major_formatter(mtick.FormatStrFormatter('%0.2f'))
        ax2.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.e'))
        ax2.set_yscale('log')

        ax3.set_xlabel(r'$x$-Coordinate [pc]')
        ax3.set_ylabel(r'$y$-Coordinate [pc]')
        ax3.yaxis.set_major_formatter(mtick.FormatStrFormatter('%0.1f'))
        ax3.xaxis.set_major_formatter(mtick.FormatStrFormatter('%0.1f'))
        ax3.set_xlim(-1.05*totsys_lim, 1.05*totsys_lim)
        ax3.set_ylim(-totsys_lim, totsys_lim)  

        ax4.set_xlabel(r'Time [kyr]')
        ax4.set_ylabel(r'Lagrangian Radius [AU]')

    anim = animation.FuncAnimation(fig, animate, interval = 100, frames=(col_len-1))
    anim.save('figures/animation'+str(count)+'.gif', writer='pillow')

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
        time[0][i][0] = vals[0].value_in(units.kyr)
        Et_array[0][i][0] = vals[1].value_in(units.J)
        dE_array[0][i][0] = vals[2]
        dEs_array[0][i][0] = vals[3]
        IMBH_birth[0][i][0] = vals[4].value_in(units.kyr)
        collisions[0][i][0] = vals[5].value_in(units.kyr)
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
        ax_.tick_params(axis="y",direction="in")
        ax_.tick_params(axis="x",direction="in")
        ax_.set_xlabel(r'Time [kyr]')
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
    ejec_parti, col_len   = file_opener('data/stability_time/*')
    com_tracker, col_len  = file_opener('data/center_of_mass/*')
    IMBH_tracker, col_len = file_opener('data/positions_IMBH/*')
    tdyn_tracker, col_len = file_opener('data/dynamical_time/*')
    Lag_tracker, col_len  = file_opener('data/lagrangians/*')
    com_x, com_y, com_z   = file_manipulator(col_len, com_tracker)

    time = np.empty((1, col_len, 1))
    LG25_array  = np.empty((1, col_len, 1))
    LG50_array  = np.empty((1, col_len, 1))
    LG75_array  = np.empty((1, col_len, 1))

    for i in range(col_len-1):
        vals = Lag_tracker.iloc[i+1]
        time[0][i+1][0] = vals[0].value_in(units.kyr)
        LG25_array[0][i+1][0] = vals[1].value_in(units.AU)
        LG50_array[0][i+1][0] = vals[2].value_in(units.AU)
        LG75_array[0][i+1][0] = vals[3].value_in(units.AU)

    line_x = np.empty((len(IMBH_tracker), col_len, 1))
    line_y = np.empty((len(IMBH_tracker), col_len, 1))
    line_z = np.empty((len(IMBH_tracker), col_len, 1))
    tdyn   = np.empty((len(IMBH_tracker), col_len, 1))

    colours = colour_picker()
    for i in range(len(IMBH_tracker)-1):
        tIMBH_tracker = IMBH_tracker.iloc[i+1]
        tIMBH_tracker = tIMBH_tracker.replace(np.NaN, "[Np.NaN, [np.NaN, np.NaN, np.NaN]")
        tDyntime_trk  = tdyn_tracker.iloc[i+1]
        for j in range(col_len):
            coords = tIMBH_tracker.iloc[j+1][1]
            if len(coords) == 1:
                pass
            else:
                line_x[i][j][0] = coords[0].value_in(units.pc)
                line_y[i][j][0] = coords[1].value_in(units.pc)
                line_z[i][j][0] = coords[2].value_in(units.pc)

                tdynval = tDyntime_trk.iloc[j+1][0]
                tdyn[i][j][0] = tdynval[0].value_in(units.yr)


    ejected_x, ejected_y, ejected_z = ejected_extract(IMBH_tracker, ejec_parti, col_len)

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
    
    ax1.set_title('Center of Mass TO FIX')
    ax2.set_title('Ejected IMBH Focus')

    for ax_ in [ax1, ax2]:
        ax_.set_xlabel(r'$x$-Coordinate [AU]')
        ax_.set_ylabel(r'$y$-Coordinate [AU]')
    ax3.set_xlabel(r'$x$-Coordinate [pc]')
    ax3.set_ylabel(r'$y$-Coordinate [pc]')

    for ax_ in [ax1, ax2, ax3, ax4]:
        ax_.yaxis.set_ticks_position('both')
        ax_.xaxis.set_ticks_position('both')
        ax_.tick_params(axis="y",direction="in")
        ax_.tick_params(axis="x",direction="in")
    
    ax2.set_xlim(-3e2, 3e2)
    ax2.set_ylim(-3e2, 3e2)

    ax3.set_xlim(-1.15*init_dist.value_in(units.pc), 1.15*init_dist.value_in(units.pc))
    ax3.set_ylim(-1.05*init_dist.value_in(units.pc), 1.05*init_dist.value_in(units.pc))        

    ax4.set_xlabel(r'Time [kyr]')
    ax4.set_ylabel(r'Lagrangian Radius [AU]')

    for i in range(len(IMBH_tracker)):
        ax1.scatter(206265*(line_x[i][1:]-com_x[0][1:]), 206265*(line_y[i][1:]-com_y[0][1:]), 
                    s = 5, c = colours[i], zorder = 1)
        ax1.scatter(206265*(line_x[i][1]-com_x[0][1]), 206265*(line_y[i][1]-com_y[0][1]), 
                    alpha = 0.7, c = colours[i], edgecolors = 'black', s = 50, zorder = 2)
        ax1.scatter(206265*(line_x[i][-1]-com_x[0][-1]), 206265*(line_y[i][-1]-com_y[0][-1]), 
                    c = colours[i], edgecolors = 'black', s = 50, zorder = 3)
        ax3.scatter(line_x[i][:-1], line_y[i][:-1], c = colours[i], zorder = 1, s = 1)
        ax3.scatter(line_x[i][-1], line_y[i][-1], c = colours[i], edgecolors = 'black', s = 50, zorder = 3)
        ax3.scatter(line_x[i][1], line_y[i][1], alpha = 0.7, c = colours[i], edgecolors = 'black', s = 50, zorder = 2)
        ax2.scatter(206265*(line_x[i][:-1]-ejected_x[0][:-1]), 206265*(line_y[i][:-1]-ejected_y[0][:-1]), 
                    c = colours[i], s = 5, zorder = 1)
        ax2.scatter(206265*(line_x[i][0]-ejected_x[0][0]), 206265*(line_y[i][0]-ejected_y[0][0]), 
                    alpha = 0.7,c = colours[i], edgecolors = 'black', s = 50, zorder=2)
        ax2.scatter(206265*(line_x[i][-2]-ejected_x[0][-2]), 206265*(line_y[i][-2]-ejected_y[0][-2]), 
                    c = colours[i], edgecolors = 'black', s = 50, zorder=3)
    ax3.scatter(0,0, color = 'black', s = 250, zorder = 2)
    ax4.plot(time[0][:-5], LG25_array[0][:-5], color = 'red',   label = r'$r_{25,L}$')
    ax4.plot(time[0][:-5], LG50_array[0][:-5], color = 'black', label = r'$r_{50,L}$')
    ax4.plot(time[0][:-5], LG75_array[0][:-5], color = 'blue',  label = r'$r_{75,L}$')
    ax4.set_ylim(0, max(LG75_array[LG25_array < 5000]))
    ax4.legend()
    ax1.yaxis.set_major_formatter(mtick.FormatStrFormatter('%0.2f'))
    plt.savefig('figures/spatial_tracker'+str(count)+'.pdf', dpi=300, bbox_inches='tight')
    plt.clf()
    plt.close()

    return

def bulk_stat_extractor(file_string):
    
    filename = glob.glob(file_string)
    data = [ ]
    
    for file_ in range(len(filename)):
        with open(filename[file_], 'rb') as input_file:
            data.append(pkl.load(input_file))

    return data

def steadytime_dist_plotter():
    steadytime_data = bulk_stat_extractor('data/stability_time/*')
    no_Data = len(steadytime_data)
    
    fin_parti  = np.empty(no_Data)
    stab_time  = np.empty(no_Data)
    init_dist  = np.empty(no_Data)
    init_mass  = np.empty(no_Data)

    for i in range(no_Data):
        sim_data = steadytime_data[i]
        fin_parti[i]  = sim_data.iloc[0][1]
        stab_time[i]  = sim_data.iloc[0][5].value_in(units.kyr)
        init_dist[i]  = sim_data.iloc[0][6].value_in(units.parsec)
        init_mass[i]  = sim_data.iloc[0][8][0].value_in(units.MSun)

    in_mass = np.unique(init_mass)
    in_dist = np.unique(init_dist)
    
    colours = colour_picker()

    for exp_ in in_mass:
        fig, ax = plt.subplots()
        tot_pop = [ ]

        for i in range(len(in_dist)): 
            N_parti_avg = [ ]
            std = [ ]
            dist_arrays = np.argwhere(init_dist == in_dist[i])
            fin_parti = fin_parti[dist_arrays]
            stab_time = stab_time[dist_arrays]
            pop_size, pop_samples = np.unique(fin_parti, return_counts=True)
            tot_pop.append(max(pop_size))

            for pop_, samp_ in zip(pop_size, pop_samples):
                N_parti = np.argwhere(fin_parti == pop_)
                N_parti_avg.append(np.mean(stab_time[N_parti]))
                std.append(np.std(stab_time[N_parti]))

            ax.errorbar(pop_size, N_parti_avg, color = colours[i], yerr=std, 
                        markeredgecolor = 'black', fmt = 'o', label = r'$r_{SMBH}=$'+str(in_dist[i])+' pc')
            ax.scatter(pop_size, np.add(N_parti_avg, std), marker = '_', color = colours[i])
            ax.scatter(pop_size, np.subtract(N_parti_avg, std), marker = '_', color = colours[i])
            ax.scatter(pop_size, N_parti_avg, color = colours[i])

            for j, xpos in enumerate(pop_size):
                if i == 0:
                    ax.text(xpos, -0.09, 'Simulations\n'+'Set '+str(i)+': '+str(pop_samples[j]), fontsize = 'xx-small', ha = 'center' )
                else:
                    ax.text(xpos, -0.09*(1+0.3*i), 'Set '+str(i)+': '+str(pop_samples[j]), fontsize = 'xx-small', ha = 'center' )
        
        order = int('{:.0f}'.format(np.log10(exp_)))
        xints = [i for i in range(1+int(max(tot_pop)))]
        ax.text(2.8, (0.05), r'$m_{IMBH}$: '+str('{:.0f}'.format(exp_/10**order))+r'$\times 10^{}$'.format(order)+r' $M_\odot$')
        ax.yaxis.set_ticks_position('both')
        ax.xaxis.set_ticks_position('both')
        ax.tick_params(axis="y",direction="in")
        ax.tick_params(axis="x",direction="in")
        ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%0.3f'))
        ax.xaxis.labelpad = 20
        plt.xticks(xints)
        plt.xlim(2.5,max(pop_size)+1)
        plt.ylim(0, 1)
        plt.ylabel(r'Ejection Time [kyr]')
        plt.xlabel(r'Number of IMBH [$N$]')
        plt.title(r'Black Hole Population vs. Stability Time')
        plt.legend()
        plt.savefig('figures/stab_time_equal_mass_'+str(exp_)+'.pdf', dpi = 300, bbox_inches='tight')

def steadytime_initmass_plotter():
    steadytime_data = bulk_stat_extractor('data/stability_time/*')
    no_Data = len(steadytime_data)
    
    fin_parti  = np.empty(no_Data)
    stab_time  = np.empty(no_Data)
    init_dist  = np.empty(no_Data)
    init_mass  = np.empty((no_Data, 2))

    for i in range(no_Data):
        sim_data = steadytime_data[i]
        fin_parti[i]  = sim_data.iloc[0][1]
        stab_time[i]  = sim_data.iloc[0][5].value_in(units.kyr)
        init_dist[i]  = sim_data.iloc[0][6].value_in(units.parsec)
        init_mass[i]  = [min(sim_data.iloc[0][8].value_in(units.MSun)), max(sim_data.iloc[0][8].value_in(units.MSun))]

    in_dist = np.unique(init_dist)
    in_mass = np.unique(init_mass, axis=0)

    colours = colour_picker()

    for exp_ in in_dist:
        fig, ax = plt.subplots()
        tot_pop = [ ]

        for i in range(len(in_mass)): 
            N_parti_avg = [ ]
            std = [ ]
            mass_arrays = np.argwhere(init_mass == in_mass[1][i])[:,0]  #Indices with the correct mass column
            fin_parti = fin_parti[mass_arrays]
            stab_time = stab_time[mass_arrays]
            pop_size, pop_samples = np.unique(fin_parti, return_counts=True)
            tot_pop.append(max(pop_size))

            for pop_, samp_ in zip(pop_size, pop_samples):
                N_parti = np.argwhere(fin_parti == pop_)
                N_parti_avg.append(np.mean(stab_time[N_parti]))
                std.append(np.std(stab_time[N_parti]))

            ax.errorbar(pop_size, N_parti_avg, color = colours[i], yerr=std, 
                        markeredgecolor = 'black', fmt = 'o', label = r'$m_{i} \in$ ['+str(in_mass[i][0])+', '+str(in_mass[i][1])+r'] $M_\odot$')
            ax.scatter(pop_size, np.add(N_parti_avg, std), marker = '_', color = colours[i])
            ax.scatter(pop_size, np.subtract(N_parti_avg, std), marker = '_', color = colours[i])
            ax.scatter(pop_size, N_parti_avg, color = colours[i])

            for j, xpos in enumerate(pop_size):
                if i == 0:
                    ax.text(xpos, -0.09*(1+0.05*i), 'Simulations\n'+'Set '+str(i)+': '+str(pop_samples[j]), fontsize = 'xx-small', ha = 'center' )
                else:
                    ax.text(xpos, -0.09*(1+0.3*i), 'Set '+str(i)+': '+str(pop_samples[j]), fontsize = 'xx-small', ha = 'center' )
        xints = [i for i in range(1+int(max(tot_pop)))]

        ax.text(2.8, 0.05, 'SMBH Distance: '+str(exp_)+' pc')
        ax.yaxis.set_ticks_position('both')
        ax.xaxis.set_ticks_position('both')
        ax.tick_params(axis="y",direction="in")
        ax.tick_params(axis="x",direction="in")
        ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%0.2f'))
        ax.xaxis.labelpad = 40
        plt.xticks(xints)
        plt.xlim(2.5,max(tot_pop)+1)
        plt.ylim(0, 1)
        plt.ylabel(r'Ejection Time [kyr]')
        plt.xlabel(r'Number of IMBH [$N$]')
        plt.title(r'Black Hole Population vs. Stability Time')
        plt.legend()
        plt.savefig('figures/stab_time_equal_dist_'+str(exp_)+'.pdf', dpi = 300, bbox_inches='tight')

steadytime_dist_plotter()
#steadytime_initmass_plotter()