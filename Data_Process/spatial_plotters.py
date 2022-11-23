from amuse.lab import *
from file_logistics import *
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import matplotlib.animation as animation
import numpy as np

class plotter_setup(object):
    def ejected_idx(self, int_string):
        """
        Function to extract the index of the 'ejected/merged' particle
        
        Inputs:
        int_string: Choice of GRX/Hermite based on simulation
        output:     Row index of the ejected particle
        """

        IMBH_tracker = file_opener('data/'+str(int_string)+'/spatial_plotters/particle_trajectory/*')
        ejec_parti = file_opener('data/'+str(int_string)+'/spatial_plotters/chaotic_simulation/*')

        return ejected_index(IMBH_tracker, ejec_parti)

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
        ax.set_xlim(min(xints)-5, max(xints)+5)
        ax.set_xlim(5, 105)

        return ax
        
    def val_filter(self, arr):
        """
        Function which removes the excessive terms
        """
        
        arr[:][abs(arr[:]) > 10**7] = np.NaN
        return arr

def colour_picker():
    """
    Colour chooser for the various plots
    """

    colors = ['red', 'blue', 'orange', 'purple', 'salmon', 'slateblue', 
              'gold', 'darkviolet', 'cornflowerblue',  'cyan',
              'lightskyblue', 'magenta',  'dodgerblue']

    return colors

def animator(init_dist, int_string):
    """
    Function to produce animations. WARNING: SLOW
    
    Inputs:
    init_dist:  The initial distance of the cluster to the central SMBH
    int_string: String to delineate between GRX and Hermite data files
    """

    print('!!!!! You have chosen to animate. This will take a while !!!!!')
    
    plot_ini = plotter_setup()
    count = file_counter(int_string)

    Lag_tracker = file_opener('data/'+str(int_string)+'/lagrangians/*')
    IMBH_tracker = file_opener('data/'+str(int_string)+'/particle_trajectory/*')
    energy_tracker = file_opener('data/'+str(int_string)+'/energy/*')
    col_len = np.shape(IMBH_tracker)[1]

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
        ax3D1.set_xlabel(r'$x$ [pc]')
        plot_ini.xy_pc(ax3D2)

    anim3D = animation.FuncAnimation(fig, animate_3D, interval = 100, frames=(col_len-1))
    anim3D.save('figures/system_evolution/animation3D'+str(count)+'.gif', writer='pillow')

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
            plot_ini.tickers(ax_, 'plot') 
        
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
    anim.save('figures/system_evolution/animation2D'+str(count)+'.gif', writer='pillow')

    return 
        
def chaos_deviate():
    """
    Function to plot a comparison of two identical simulations where in one there is a 1e-15 
    [pc] perturbation in the particles x coordinate.
    """

    plot_ini = plotter_setup()
    IMBH_tracker = file_opener('data/GRX/chaos_example/system/particle_trajectory/*')
    energy_tracker = file_opener('data/GRX/chaos_example/energy/*')
    IMBH_tracker_pert = file_opener('data/GRX/chaos_example/pert_system/particle_trajectory/*')
    
    col_len = np.shape(IMBH_tracker)[1] - 1
    col_len_pert = np.shape(IMBH_tracker_pert)[1] - 1
    col_len = min(col_len, col_len_pert)

    focus_idx = 1
    particle = IMBH_tracker.iloc[focus_idx]
    SMBH_data = IMBH_tracker.iloc[0]
    particle_pert = IMBH_tracker_pert.iloc[focus_idx]
    SMBH_data_pert = IMBH_tracker_pert.iloc[0]

    time = np.empty(col_len)
    neigh_key = np.empty((3, col_len))
    nearest_neigh = np.empty(col_len)
    SMBH_dist = np.empty(col_len)
    neigh_key_pert = np.empty((3, col_len))
    nearest_neigh_pert = np.empty(col_len)
    SMBH_dist_pert = np.empty(col_len)

    semi_SMBH = np.empty(col_len)
    semi_SMBH_pert = np.empty(col_len)
    
    line_x = np.empty((col_len))
    line_y = np.empty((col_len))
    line_x_pert = np.empty((col_len))
    line_y_pert = np.empty((col_len))
    rel_phase = np.empty((col_len))
    for j in range(col_len):
        SMBH_coords = SMBH_data.iloc[j][2]
        SMBH_coords_pert = SMBH_data_pert.iloc[j][2]
        SMBH_vel = SMBH_data.iloc[j][3]
        SMBH_vel_pert = SMBH_data_pert.iloc[j][3]
        time_arr = energy_tracker.iloc[j]
        time[j] = time_arr[6].value_in(units.Myr)

        coords = particle.iloc[j][2]
        vel = particle.iloc[j][3]
        coords_pert = particle_pert.iloc[j][2]
        vel_pert = particle_pert.iloc[j][3]
        if len(coords) != 1:
            line_x[j] = (coords[0]-SMBH_coords[0]).value_in(units.pc)
            line_y[j] = (coords[1]-SMBH_coords[1]).value_in(units.pc)
            line_z = (coords[2]-SMBH_coords[2]).value_in(units.pc)
            SMBH_dist[j] = np.sqrt(line_x[j]**2+line_y[j]**2+line_z**2)
            velx = (vel[0]-SMBH_vel[0]).value_in(units.kms)
            vely = (vel[1]-SMBH_vel[1]).value_in(units.kms)
            velz = (vel[2]-SMBH_vel[2]).value_in(units.kms)

            line_x_pert[j] = (coords_pert[0]-SMBH_coords_pert[0]).value_in(units.pc)
            line_y_pert[j] = (coords_pert[1]-SMBH_coords_pert[1]).value_in(units.pc)
            line_z_pert = (coords_pert[2]-SMBH_coords_pert[2]).value_in(units.pc)
            SMBH_dist_pert[j] = np.sqrt(line_x_pert[j]**2+line_y_pert[j]**2+line_z_pert**2)
            velx_pert = (vel_pert[0]-SMBH_vel_pert[0]).value_in(units.kms)
            vely_pert = (vel_pert[1]-SMBH_vel_pert[1]).value_in(units.kms)
            velz_pert = (vel_pert[2]-SMBH_vel_pert[2]).value_in(units.kms)

            rel_x = (line_x_pert[j]-line_x[j])**2
            rel_y = (line_y_pert[j]-line_y[j])**2
            rel_z = (line_z_pert-line_z)**2
            rel_dist = np.sqrt(rel_x+rel_y+rel_z)

            rel_vx = (velx_pert-velx)**2
            rel_vy = (vely_pert-vely)**2
            rel_vz = (velz_pert-velz)**2
            rel_vel = np.sqrt(rel_vx+rel_vy+rel_vz)

            rel_phase[j] = 0.5 * np.log(rel_dist + rel_vel)

        time_step = particle.iloc[j]
        time_step_pert = particle_pert.iloc[j]
        semi_SMBH[j] = time_step[7][0].value_in(units.parsec)
        semi_SMBH_pert[j] = time_step_pert[7][0].value_in(units.parsec)
        for i in range(3):
            neigh_key[i][j] = time_step[6][i]
            neigh_key_pert[i][j] = time_step_pert[6][i]
        nearest_neigh[j] = abs(time_step[-1])
        nearest_neigh_pert[j] = abs(time_step_pert[-1])

    smoothing = round(0.05 * col_len)
    normalise = semi_SMBH[0]
    semi_SMBH /= (normalise)
    semi_SMBH_pert /= (normalise)
    semi_SMBH_smooth = plot_ini.moving_average(semi_SMBH, smoothing)
    semi_SMBH_smooth_pert = plot_ini.moving_average(semi_SMBH_pert, smoothing)

    time_smooth = plot_ini.moving_average(time, smoothing)
    nearest_smooth = plot_ini.moving_average(nearest_neigh, smoothing)
    nearest_smooth_pert = plot_ini.moving_average(nearest_neigh_pert, smoothing)
    SMBH_dist_smooth = plot_ini.moving_average(SMBH_dist, smoothing)
    SMBH_dist_smooth_pert = plot_ini.moving_average(SMBH_dist_pert, smoothing)
    phase_smooth = plot_ini.moving_average(rel_phase, round(0.1*smoothing))

    fig = plt.figure(figsize=(12.5, 10))
    ax1 = fig.add_subplot(221)
    ax2 = fig.add_subplot(222)
    ax3 = fig.add_subplot(223)
    ax4 = fig.add_subplot(224)

    plot_ini.tickers(ax1, 'plot')
    ax2.set_title(r'Time vs. $\ln [\delta w(\vec r,\vec v)]$')
    ax3.set_title(r'Time vs. $a_{SMBH}$')
    ax4.set_title('Time vs. Distance to Nearest Neighbour')
    ax1.set_xlabel(r'$x$ [pc]')
    ax1.set_ylabel(r'$y$ [pc]')
    ax2.set_ylabel(r'$\ln [\delta w(\vec r,\vec v)]$')
    ax3.set_ylabel(r'$a(t)/a_{\rm{SMBH},0}^{\rm{Hermite}}$')
    ax4.set_ylabel(r'$r_{\rm {nn}}$ [pc]')
    for ax_ in [ax2, ax3, ax4]:
        ax_.set_xlabel('Time [Myr]')
        plot_ini.tickers(ax_, 'plot') 
        ax_.set_xlim(0,max(time_smooth))
    ax1.set_xlim(-0.5,0.55)
    ax1.set_ylim(-0.5,0.5)
    ax2.set_ylim(0,1.05*max(phase_smooth))

    steps = round(0.35*col_len)
    ax1.scatter(line_x[:steps], line_y[:steps], s = 5, color = 'red', label = r'$\delta = 0$')
    ax1.scatter(line_x_pert[:steps], line_y_pert[:steps], s = 5, color = 'blue', label = r'$\delta = 10^{-13}$ ')
    ax1.scatter(0, 0, s = 250, color = 'black')

    ax2.plot(time[:len(phase_smooth)], phase_smooth, color = 'black')
    
    ax3.plot(time_smooth, semi_SMBH_smooth, color = 'red')
    ax3.plot(time_smooth, semi_SMBH_smooth_pert, color = 'blue')

    ax4.plot(time_smooth, SMBH_dist_smooth, color = 'red', label = r'$r_{\rm{SMBH}}$')
    ax4.plot(time_smooth, SMBH_dist_smooth_pert, color = 'blue', )
    ax4.plot(time_smooth, nearest_smooth, color = 'red', linestyle = ':', label = r'$r_{\rm{NN}}$')
    ax4.plot(time_smooth, nearest_smooth_pert, color = 'blue', linestyle = ':')

    ax1.legend()
    ax4.legend()

    plt.savefig('figures/system_evolution/GRX_vs_Hermite_sys_evol.pdf', dpi=300, bbox_inches='tight')
    
def ejected_evolution(int_string):
    """
    Function which plots various Kepler elements of the IMBH particle stopping the simulation
    
    output: Plots of the inclination,semi-major axis, nearest neighbour and eccentricity of the 
            merged/ejected particle relative to the SMBH, nearest neighbour and second-nearest neighbour
    """
    
    plot_ini = plotter_setup()
    count = file_counter(int_string)
    IMBH_tracker = file_opener('data/'+str(int_string)+'/spatial_plotters/particle_trajectory/*')
    energy_tracker = file_opener('data/'+str(int_string)+'/spatial_plotters/energy/*')

    col_len = np.shape(IMBH_tracker)[1] - 1
    smoothing = round(0.1 * col_len)
    
    focus_idx = plot_ini.ejected_idx(int_string)
    ejec_particle = IMBH_tracker.iloc[focus_idx]
    SMBH_data = IMBH_tracker.iloc[0]

    time = np.empty(col_len)
    neigh_key = np.empty((3, col_len))
    nearest_neigh = np.empty(col_len)
    SMBH_dist = np.empty(col_len)

    ejec_semi_SMBH = np.empty(col_len)
    ejec_semi_bin = np.empty(col_len)
    ejec_semi_ter = np.empty(col_len)
    
    ejec_ecc_SMBH = np.empty(col_len)
    ejec_ecc_bin = np.empty(col_len)
    ejec_ecc_ter = np.empty(col_len)

    ejec_incl_SMBH = np.empty(col_len)
    ejec_incl_bin = np.empty(col_len)
    ejec_incl_ter = np.empty(col_len)

    for j in range(col_len):
        time_step = ejec_particle.iloc[j]
        time_arr = energy_tracker.iloc[j]
        SMBH_coords = SMBH_data.iloc[j][2]

        ejec_semi_SMBH[j] = time_step[7][0].value_in(units.parsec)
        ejec_semi_bin[j]  = time_step[7][1].value_in(units.parsec)
        ejec_semi_ter[j]  = time_step[7][2].value_in(units.parsec)

        ejec_ecc_SMBH[j] = (1-time_step[8][0])
        ejec_ecc_bin[j]  = (1-time_step[8][1])
        ejec_ecc_ter[j]  = (1-time_step[8][2])

        ejec_incl_SMBH[j] = time_step[9][0]
        ejec_incl_bin[j]  = time_step[9][1]
        ejec_incl_ter[j]  = time_step[9][2]

        for i in range(3):
            neigh_key[i][j] = time_step[6][i]

        nearest_neigh[j] = abs(time_step[-1])
        time[j] = time_arr[6].value_in(units.Myr)
        
        line_x = (time_step[2][0] - SMBH_coords[0])
        line_y = (time_step[2][1] - SMBH_coords[1])
        line_z = (time_step[2][2] - SMBH_coords[2])
        SMBH_dist[j] = np.sqrt(line_x**2+line_y**2+line_z**2).value_in(units.pc)

    normalise = ejec_semi_SMBH[0]
    ejec_semi_SMBH /= (normalise)
    ejec_semi_bin  /= (normalise)
    ejec_semi_ter  /= (normalise)

    ejec_semi_SMBH_smooth = plot_ini.moving_average(ejec_semi_SMBH, smoothing)

    ejec_ecc_SMBH_smooth = plot_ini.moving_average(ejec_ecc_SMBH, smoothing)
    ejec_ecc_bin_smooth = plot_ini.moving_average(ejec_ecc_bin, smoothing)
    ejec_ecc_ter_smooth = plot_ini.moving_average(ejec_ecc_ter, smoothing)

    ejec_incl_SMBH_smooth = plot_ini.moving_average(ejec_incl_SMBH, smoothing)
    ejec_incl_bin_smooth = plot_ini.moving_average(ejec_incl_bin, smoothing)
    ejec_incl_ter_smooth = plot_ini.moving_average(ejec_incl_ter, smoothing)

    time_smooth = plot_ini.moving_average(time, smoothing)
    nearest_smooth = plot_ini.moving_average(nearest_neigh, smoothing)
    SMBH_dist_smooth = plot_ini.moving_average(SMBH_dist, smoothing)

    fig = plt.figure(figsize=(12.5, 10))
    ax1 = fig.add_subplot(221)
    ax2 = fig.add_subplot(222)
    ax3 = fig.add_subplot(223)
    ax4 = fig.add_subplot(224)

    ax1.set_title('Time vs. Eccentricity')
    ax2.set_title('Time vs. Semi-Major Axis')
    ax3.set_title('Time vs. Inclination')
    ax4.set_title('Time vs. Distance to Nearest Neighbour')
    ax1.set_ylabel(r'$(1-e)$')
    ax2.set_ylabel(r'$a(t)/a_{SMBH,0}$')
    ax3.set_ylabel(r'$i$ [deg]')
    ax4.set_ylabel(r'$r_{nn}$ [pc]')
    for ax_ in [ax1, ax2, ax3, ax4]:
        ax_.set_xlabel('Time [Myr]')
        plot_ini.tickers(ax_, 'plot') 
        ax_.set_xlim(0,max(time_smooth))
    ax1.set_ylim(0,1)
    ax3.set_ylim(-180,180)

    ax1.plot(time, ejec_ecc_SMBH, color = 'black', alpha = 0.3, linestyle = ':')
    ax1.plot(time_smooth, ejec_ecc_SMBH_smooth, color = 'black', label = 'w.r.t SMBH')
    ax1.plot(time_smooth, ejec_ecc_bin_smooth, color = 'red', label = 'w.r.t Binary')
    ax1.plot(time_smooth, ejec_ecc_ter_smooth, color = 'blue', label = 'w.r.t Tertiary')

    ax2.plot(time, ejec_semi_SMBH, color = 'black', alpha = 0.3, linestyle = ':')
    ax2.plot(time_smooth, ejec_semi_SMBH_smooth, color = 'black')

    ax3.plot(time_smooth, ejec_incl_SMBH_smooth, color = 'black')
    ax3.plot(time_smooth, ejec_incl_bin_smooth, color = 'red')
    ax3.plot(time_smooth, ejec_incl_ter_smooth, color = 'blue')

    ax4.plot(time, SMBH_dist, color = 'black', alpha = 0.3, linestyle = ':')
    ax4.plot(time, nearest_neigh, color = 'red', alpha = 0.3, linestyle = ':')
    ax4.plot(time_smooth, nearest_smooth, color = 'red')
    ax4.plot(time_smooth, SMBH_dist_smooth, color = 'black')

    ax1.legend()
    plt.savefig('figures/system_evolution/ejec_bin_trip_evol_'+str(count)+'.pdf', dpi=300, bbox_inches='tight')
    plt.clf()

def energy_plotter(int_string):
    """
    Function to plot the energy evolution of the system

    output: Energy evolution plot of the system
    """

    plot_ini = plotter_setup()
    count = file_counter(int_string)
    energy_tracker = file_opener('data/'+str(int_string)+'/spatial_plotters/energy/*')
    col_len = np.shape(energy_tracker)[0]

    time = np.empty((col_len - 1))
    Et_array = np.empty((col_len - 1))
    dE_array = np.empty((col_len - 1))
    KE_array = np.empty((col_len - 1))
    PE_array = np.empty((col_len - 1))
    IMBHapp_array = np.empty((col_len - 1))
    colltime_array = np.empty((col_len - 1)) 
    merger_mass = np.empty((col_len - 1))

    for i in range(col_len):
        if i != 0:
            vals = energy_tracker.iloc[i]
            IMBHapp_array[i-1] = vals[0].value_in(units.Myr)
            merger_mass[i-1] = vals[1].value_in(units.MSun)
            colltime_array[i-1] = vals[2].value_in(units.Myr)
            Et_array[i-1] = vals[3].value_in(units.J)
            KE_array[i-1] = vals[4].value_in(units.J)
            time[i-1] = vals[6].value_in(units.Myr)
            dE_array[i-1] = vals[7]
            PE_array[i-1] = vals[9].value_in(units.J)

    fig = plt.figure(figsize=(12.5, 8))
    ax1 = fig.add_subplot(221)
    ax2 = fig.add_subplot(222)
    ax1.set_title('Energy Evolution in Time')
    ax2.set_title('Energy Error')

    for ax_ in [ax1, ax2]:
        plot_ini.tickers(ax_, 'plot')  
        ax_.set_xlabel(r'Time [Myr]')
        ax_.set_yscale('log')
    ax1.set_ylabel(r'Energy [J]')
    ax2.set_ylabel(r'$\frac{|E(t)-E_0|}{|E_0|}$')

    ax1.plot(time, abs(Et_array), color = 'black', label = r'$|E_{tot}|$', zorder = 3)
    ax1.plot(time, KE_array, color = 'red', label = r'$K_E$', zorder = 1)
    ax1.plot(time, abs(PE_array), color = 'blue', label = r'$|P_E|$', zorder = 2)
    ax2.plot(time[:-5], dE_array[:-5], color = 'black')

    ax1.set_ylim(0.5*min(abs(KE_array)), 5*max(abs(PE_array)))
    ax1.legend()

    plt.savefig('figures/system_evolution/energy_tracker'+str(count)+'.pdf', dpi=300, bbox_inches='tight')
    plt.clf()
    plt.close()

    return

def nearest_neigh(int_string):
    """
    Function which plots the evolution of all nearest neighbour distances in a specific simulation.
    """
    
    plot_ini = plotter_setup()
    count = file_counter(int_string)
    IMBH_tracker = file_opener('data/'+str(int_string)+'/spatial_plotters/particle_trajectory/*')
    energy_tracker = file_opener('data/'+str(int_string)+'/spatial_plotters/energy/*')

    col_len = np.shape(IMBH_tracker)[1] - 1
    no_parti = np.shape(IMBH_tracker)[0] - 1
    smoothing = round(10**-2 * col_len)
    

    time = np.empty(col_len)
    nn_dist = np.empty((no_parti, col_len))
    nn_dist_smooth = np.empty((no_parti, (col_len-smoothing+1)))

    for parti_ in range(no_parti):
        time_step = IMBH_tracker.iloc[parti_+1]
        for j in range(col_len):
            nn_dist[parti_][j] = time_step.iloc[j][-1]
        nn_dist_smooth[parti_] = plot_ini.moving_average(nn_dist[parti_], smoothing)
        nn_dist_smooth[parti_] = np.asarray(nn_dist_smooth[parti_])

    for i in range(col_len):
        time_arr = energy_tracker.iloc[i]
        time[i] = time_arr[6].value_in(units.Myr)
    time_smooth = plot_ini.moving_average(time, smoothing)
    time = np.asarray(time)

    fig, ax = plt.subplots()
    ax.set_title('Time vs. Distance')
    ax.set_ylabel(r'$r_{NN}$')
    ax.set_xlabel(r'$t$ [Myr]')
    plot_ini.tickers(ax, 'plot') 
    
    for parti_ in range(no_parti):
        ax.plot(time_smooth, np.log10(nn_dist_smooth[parti_]), color = 'black', alpha = 0.5)
        #ax.plot(time, np.log10(nn_dist[parti_]), color = 'black')
    plt.savefig('figures/system_evolution/NN_dist_Evolution'+str(count)+'.pdf', dpi=300, bbox_inches='tight')
    plt.clf()

def spatial_plotter(int_string):
    """
    Function to plot the evolution of the system

    output: The spatial evolution of the system
    """

    plot_ini = plotter_setup()
    count = file_counter(int_string)

    IMBH_tracker = file_opener('data/'+str(int_string)+'/spatial_plotters/particle_trajectory/*')
    energy_tracker = file_opener('data/'+str(int_string)+'/spatial_plotters/energy/*')
    col_len_raw = np.shape(IMBH_tracker)[1]
    col_len = round(col_len_raw**0.6)
    parti_size = 20+len(IMBH_tracker)**-0.5

    line_x = np.empty((len(IMBH_tracker), col_len))
    line_y = np.empty((len(IMBH_tracker), col_len))
    line_z = np.empty((len(IMBH_tracker), col_len))

    for i in range(len(IMBH_tracker)):
        tIMBH_tracker = IMBH_tracker.iloc[i]
        for j in range(col_len):
            coords = tIMBH_tracker.iloc[j][2]
            if len(coords) == 1:
                pass
            else:
                line_x[i][j] = coords[0].value_in(units.pc)
                line_y[i][j] = coords[1].value_in(units.pc)
                line_z[i][j] = coords[2].value_in(units.pc)

    focus_idx = plot_ini.ejected_idx(int_string)
    focus_particle = IMBH_tracker.iloc[focus_idx]

    focus_x = np.empty((col_len))
    focus_y = np.empty((col_len))
    focus_z = np.empty((col_len))

    for j in range(col_len):
        coords = focus_particle.iloc[j][2]
        focus_x[j] = coords[0].value_in(units.pc)
        focus_y[j] = coords[1].value_in(units.pc)
        focus_z[j] = coords[2].value_in(units.pc)

    for arr_ in [focus_x, focus_y, focus_z, line_x, line_y, line_z]:
        plot_ini.val_filter(arr_)

    time = np.empty((col_len_raw - 1))
    dE_array = np.empty((col_len_raw - 1))

    for i in range(col_len_raw):
        if i == 0:
            pass
        else:
            vals = energy_tracker.iloc[i]
            time[i-1] = vals[6].value_in(units.Myr)
            dE_array[i-1] = vals[7]

    colours = colour_picker()
    fig = plt.figure(figsize=(12.5, 15))
    ax1 = fig.add_subplot(321)
    ax2 = fig.add_subplot(322)
    ax3 = fig.add_subplot(323)
    ax4 = fig.add_subplot(324)
    
    ax1.set_title('Overall System')
    ax2.set_title('Energy Error vs. Time')
    ax1.xaxis.set_major_locator(plt.MaxNLocator(3))
    ax1.yaxis.set_major_locator(plt.MaxNLocator(3))
    for ax_ in [ax1, ax2, ax3, ax4]:
        plot_ini.tickers(ax_, 'plot') 

    xaxis_lim = 1.05*np.nanmax(abs(line_x-line_x[0]))
    yaxis_lim = 1.05*np.nanmax(abs(line_y-line_y[0]))
    zaxis_lim = 1.05*np.nanmax(abs(line_z-line_z[0]))

    ax1.set_xlim(-abs(xaxis_lim), abs(xaxis_lim))
    ax1.set_ylim(-abs(yaxis_lim), yaxis_lim)
    ax3.set_xlim(-abs(xaxis_lim), abs(xaxis_lim))
    ax3.set_ylim(-abs(zaxis_lim), zaxis_lim)
    ax4.set_xlim(-abs(yaxis_lim), abs(yaxis_lim))
    ax4.set_ylim(-abs(zaxis_lim), zaxis_lim)
    ax2.set_yscale('log')

    ax1.set_xlabel(r'$x$ [pc]')
    ax1.set_ylabel(r'$y$ [pc]')
    ax2.set_xlabel(r'Time [Myr]')
    ax2.set_ylabel(r'$\frac{|E(t)-E_0|}{|E_0|}$')
    ax3.set_xlabel(r'$x$ [pc]')
    ax3.set_ylabel(r'$z$ [pc]')
    ax4.set_xlabel(r'$y$ [pc]')
    ax4.set_ylabel(r'$z$ [pc]')
    iter = -1
    
    for i in range(len(IMBH_tracker)):
        iter += 1
        if iter > len(colours):
            iter = 0

        if i == 0:
            adapt_c = 'black'
            ax1.scatter((line_x[i]-line_x[0]), (line_y[i]-line_y[0]), 
                         c = adapt_c, zorder = 1, s = 250)
            ax3.scatter((line_x[i]-line_x[0]), (line_z[i]-line_z[0]), 
                         c = adapt_c, zorder = 1, s = 250)
            ax4.scatter((line_z[i]-line_z[0]), (line_y[i]-line_y[0]), 
                         c = adapt_c, zorder = 1, s = 250)
        else:
            ax1.scatter(line_x[i][-1]-line_x[0][-1], line_y[i][-1]-line_y[0][-1], 
                        c = colours[iter-2], edgecolors = 'black', s = parti_size, zorder = 3)
            ax1.scatter(line_x[i]-line_x[0], line_y[i]-line_y[0], 
                        c = colours[iter-2], s = 1, zorder = 1) 

            ax3.scatter(line_x[i][-1]-line_x[0][-1], line_z[i][-1]-line_z[0][-1], 
                        c = colours[iter-2], edgecolors = 'black', s = parti_size, zorder = 3)
            ax3.scatter(line_x[i]-line_x[0], line_z[i]-line_z[0], 
                        c = colours[iter-2], s = 1, zorder = 1) 

            ax4.scatter(line_y[i][-1]-line_y[0][-1], line_z[i][-1]-line_z[0][-1], 
                        c = colours[iter-2], edgecolors = 'black', s = parti_size, zorder = 3)
            ax4.scatter(line_y[i]-line_y[0], line_z[i]-line_z[0], 
                        c = colours[iter-2], s = 1, zorder = 1) 
    ax2.plot(time[:-5], dE_array[:-5], color = 'black')
    plt.savefig('figures/system_evolution/simulation_evolution_'+str(count)+'.pdf', dpi=300, bbox_inches='tight')
    plt.clf()
    plt.close()     

    fig = plt.figure(figsize=(8, 8))
    ax3D = fig.add_subplot(121, projection="3d")
    iter = -1
    for i in range(len(IMBH_tracker)):
        iter += 1
        if iter > len(colours):
            iter = 0
        if i == 0:
            pass
        else:
            ax3D.scatter(line_x[i]-line_x[0], 
                         line_y[i]-line_y[0], 
                         line_z[i]-line_z[0], 
                         c = colours[iter-2], s = 1, zorder = 1)
    ax3D.scatter(0, 0, 0, color = 'black', s = 150, zorder = 2)
    ax3D.xaxis.set_major_formatter(mtick.FormatStrFormatter('%0.2f'))
    ax3D.yaxis.set_major_formatter(mtick.FormatStrFormatter('%0.2f'))
    ax3D.zaxis.set_major_formatter(mtick.FormatStrFormatter('%0.2f'))
    ax3D.set_xlabel(r'$x$ [pc]')
    ax3D.set_ylabel(r'$y$ [pc]')
    ax3D.set_zlabel(r'$z$ [pc]')
    ax3D.view_init(30, 160)
    plt.savefig('figures/system_evolution/simulation_evolution_3D_'+str(count)+'.pdf', dpi=300, bbox_inches='tight')

    return

#direct_comparison('Hermite')
#spatial_plotter('GRX')
#energy_plotter('Hermite')
#ejected_evolution('Hermite')
#chaos_deviate()