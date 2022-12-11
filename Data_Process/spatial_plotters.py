from amuse.lab import *
from file_logistics import *
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import matplotlib.animation as animation
import numpy as np
import warnings

np.seterr(divide='ignore')
warnings.filterwarnings("ignore", category=RuntimeWarning) 

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
    ptracker = file_opener('/media/erwanh/Elements/'+(int_string)+'/particle_trajectory/*')
    etracker = file_opener('data/'+str(int_string)+'/energy/*')

    col_len = np.shape(ptracker)[1]
    time = np.empty((col_len, 1))
    dE_array = np.empty((col_len, 1))

    for i in range(col_len):
        vals = etracker.iloc[i]
        time[i][0] = vals[0].value_in(units.Myr)
        dE_array[i][0] = vals[2]

    col_len = len(ptracker.iloc[0])-2
    line_x = np.empty((len(ptracker), col_len, 1))
    line_y = np.empty((len(ptracker), col_len, 1))
    line_z = np.empty((len(ptracker), col_len, 1))

    for i in range(len(ptracker)):
        tptracker = ptracker.iloc[i]
        tptracker = tptracker.replace(np.NaN, "[Np.NaN, [np.NaN, np.NaN, np.NaN]")
        for j in range(col_len):
            coords = tptracker.iloc[j+1][2]
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
    c = colour_picker()

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

        for i in range(len(ptracker)):
            iter += 1
            if iter > len(c):
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
                              s =1, c = c[iter-2], lw = 2)
                ax3D1.scatter((line_x[i][skip_zeroth]-line_x[0][skip_zeroth]), 
                              (line_y[i][skip_zeroth]-line_y[0][skip_zeroth]), 
                              (line_z[i][skip_zeroth]-line_z[0][skip_zeroth]), 
                               c = c[iter-2], edgecolors = 'black', s = 50)

                ax3D2.scatter((line_x[i][max(0,skip_zeroth-20):skip_zeroth+1]-line_x[1][max(0,skip_zeroth-20):skip_zeroth+1]),
                              (line_y[i][max(0,skip_zeroth-20):skip_zeroth+1]-line_y[1][max(0,skip_zeroth-20):skip_zeroth+1]),
                              (line_z[i][max(0,skip_zeroth-20):skip_zeroth+1]-line_z[1][max(0,skip_zeroth-20):skip_zeroth+1]), 
                               s =1, c = c[iter-2], lw = 2)
                ax3D2.scatter((line_x[i][skip_zeroth]-line_x[1][skip_zeroth]), 
                              (line_y[i][skip_zeroth]-line_y[1][skip_zeroth]), 
                              (line_z[i][skip_zeroth]-line_z[1][skip_zeroth]),
                               c = c[iter-2], edgecolors = 'black', s = 40)

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
        
        for i in range(len(ptracker)):
            ax1.scatter(line_x[i][max(0,skip_zeroth-20):skip_zeroth+1]-line_x[1][max(0,skip_zeroth-20):skip_zeroth+1],
                        line_y[i][max(0,skip_zeroth-20):skip_zeroth+1]-line_y[1][max(0,skip_zeroth-20):skip_zeroth+1],
                        c = c[i], lw = 2, alpha = 0.7, s = 1)
            ax1.scatter(line_x[i][skip_zeroth]-line_x[1][skip_zeroth], line_y[i][skip_zeroth]-line_t[1][skip_zeroth],
                        c = c[i], edgecolors = 'black', s = 40)

            ax3.scatter((line_x[i][max(0,skip_zeroth-20):skip_zeroth+1]-line_x[0][max(0,skip_zeroth-20):skip_zeroth+1]), 
                        (line_y[i][max(0,skip_zeroth-20):skip_zeroth+1]-line_y[0][max(0,skip_zeroth-20):skip_zeroth+1]), 
                        s =1, c = c[i], lw = 2)
            ax3.scatter((line_x[i][skip_zeroth]-line_x[0][skip_zeroth]), 
                        (line_y[i][skip_zeroth]-line_y[0][skip_zeroth]), 
                        c = c[i], edgecolors = 'black', s = 50)

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
    Function to plot a comparison of two identical simulations.
    """

    plot_ini = plotter_setup()
    ptracker = file_opener('data/GRX/chaos_example/system/particle_trajectory/*')
    etracker = file_opener('data/GRX/chaos_example/energy/*')
    ptracker_pert = file_opener('data/GRX/chaos_example/pert_system/particle_trajectory/*')
    
    col_len = np.shape(ptracker)[1] - 1
    col_len_pert = np.shape(ptracker_pert)[1] - 1
    col_len = min(col_len, col_len_pert)

    focus_idx = 1
    particle = ptracker.iloc[focus_idx]
    SMBH_data = ptracker.iloc[0]
    particle_pert = ptracker_pert.iloc[focus_idx]
    SMBH_data_pert = ptracker_pert.iloc[0]

    neigh_key = np.empty((3, col_len))
    NN = np.empty(col_len)
    SMBH_dist = np.empty(col_len)
    semi_SMBH = np.empty(col_len)
    line_x = np.empty((col_len))
    line_y = np.empty((col_len))

    neigh_key_pert = np.empty((3, col_len))
    NN_pert = np.empty(col_len)
    SMBH_dist_pert = np.empty(col_len)
    semi_SMBH_pert = np.empty(col_len)
    line_x_pert = np.empty((col_len))
    line_y_pert = np.empty((col_len))
    rel_phase = np.empty((col_len))

    time = np.empty(col_len)
    for j in range(col_len):
        SMBH_coords = SMBH_data.iloc[j][2]
        SMBH_coords_pert = SMBH_data_pert.iloc[j][2]
        SMBH_vel = SMBH_data.iloc[j][3]
        SMBH_vel_pert = SMBH_data_pert.iloc[j][3]
        time_arr = etracker.iloc[j]
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
        NN[j] = abs(time_step[-1])
        NN_pert[j] = abs(time_step_pert[-1])

    smoothing = round(0.05 * col_len)
    normalise = semi_SMBH[0]
    semi_SMBH /= (normalise)
    semi_SMBH_pert /= (normalise)
    semi_SMBH_smooth = plot_ini.moving_average(semi_SMBH, smoothing)
    semi_SMBH_smooth_pert = plot_ini.moving_average(semi_SMBH_pert, smoothing)

    time_smooth = plot_ini.moving_average(time, smoothing)
    nearest_smooth = plot_ini.moving_average(NN, smoothing)
    nearest_smooth_pert = plot_ini.moving_average(NN_pert, smoothing)
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
    ax1.set_xlim(-0.45,0.9)
    ax1.set_ylim(-0.25,1)
    ax2.set_ylim(0,1.05*max(phase_smooth))

    steps = round(0.6*col_len)

    ax1.scatter(line_x[:steps], line_y[:steps], s = 5, color = 'red')
    ax1.scatter(line_x_pert[:steps], line_y_pert[:steps], s = 5, color = 'blue')
    ax1.scatter(0, 0, s = 250, color = 'black')
    ax2.plot(time[:len(phase_smooth)], phase_smooth, color = 'black')
    ax3.plot(time_smooth, semi_SMBH_smooth, color = 'red')
    ax3.plot(time_smooth, semi_SMBH_smooth_pert, color = 'blue')
    ax4.plot(time_smooth, SMBH_dist_smooth, color = 'red', label = r'$r_{\rm{SMBH}}$')
    ax4.plot(time_smooth, SMBH_dist_smooth_pert, color = 'blue', )
    ax4.plot(time_smooth, nearest_smooth, color = 'red', linestyle = ':', label = r'$r_{\rm{NN}}$')
    ax4.plot(time_smooth, nearest_smooth_pert, color = 'blue', linestyle = ':')

    ax4.legend()
    plt.savefig('figures/system_evolution/GRX_vs_Hermite_sys_evol.pdf', dpi=300, bbox_inches='tight')

def ejected_evolution():
    """
    Function which plots various Kepler elements of the IMBH particle stopping the simulation
    """
    
    plot_ini = plotter_setup()
    
    integrator = ['Hermite', 'GRX']
    ejec_ecc_arr = [[ ], [ ]]
    ejec_vel_arr = [[ ], [ ]]
    ejec_NNdist_arr = [[ ], [ ]]
    ejec_distSMBH_arr = [[ ], [ ]]
    iter = -1
    for int_ in integrator:   
        iter += 1
        data = natsort.natsorted(glob.glob('/media/erwanh/Elements/'+(int_)+'/particle_trajectory/*'))
        if int_ == 'GRX':
            chaotic = ['data/GRX/no_addition/chaotic_simulation/'+str(i[47:]) for i in data]
            energy = ['data/GRX/energy/'+str(i[47:]) for i in data]
        else:
            chaotic = ['data/Hermite/no_addition/chaotic_simulation/'+str(i[51:]) for i in data]
            energy = ['data/Hermite/energy/'+str(i[51:]) for i in data]
        filter = 0
        high_ecc = 0
        
        ejec_ecc = []
        distSMBH = []
        ejec_time = []
        vel_arr = []
        NNdist = []
        ejec_ecc_smooth = []
        distSMBH_smooth = []
        ejec_time_smooth = []
        vel_smooth = []
        NNdist_smooth = []
        cropped_idx = []
        all_idx = []
        SMBHx = []
        SMBHy = []
        SMBHz = []
        SMBH_sample = [[], [], []]
        
        for file_ in range(len(data)):
            with open(chaotic[file_], 'rb') as input_file:
                print('Reading file : ', input_file)
                chaotic_tracker = pkl.load(input_file)
                if chaotic_tracker.iloc[0][6] <= 50:     # TURN OFF TO FIND THE COMPLETE ECCENTRIC STATISTICS
                    if chaotic_tracker.iloc[0][-4] > 0 or chaotic_tracker.iloc[0][5] > 0:
                        with open(data[file_], 'rb') as input_file:
                            ptracker = pkl.load(input_file)

                        with open(energy[file_], 'rb') as input_file:
                            etracker = pkl.load(input_file)
                            
                        col_len = np.shape(ptracker)[1]-1
                        smoothing = round(0.1 * col_len)
                        focus_idx, res = ejected_index(ptracker, chaotic_tracker, int_)
                        if res == 'merger' or res == 'ejected' and focus_idx != 0:
                            filter += 1
                            ejec_particle = ptracker.iloc[focus_idx]
                            SMBH_data = ptracker.iloc[0]

                            time = np.empty(col_len)
                            neigh_key = np.empty((3, col_len))
                            NN = np.empty(col_len)
                            SMBH_dist = np.empty(col_len)
                            pvel = np.empty(col_len)
                            ejec_KE = np.empty(col_len)
                            ejec_ecc_SMBH = np.empty(col_len)

                            for j in range(col_len):
                                time_step = ejec_particle.iloc[j]
                                time_arr = etracker.iloc[j]
                                SMBH_coords = SMBH_data.iloc[j]

                                ejec_KE[j] = time_step[4].value_in(units.J)
                                ejec_ecc_SMBH[j] = np.log10(1-time_step[8][0])

                                for i in range(3):
                                    neigh_key[i][j] = time_step[6][i]

                                NN[j] = np.log10(time_step[-1])
                                time[j] = time_arr[6].value_in(units.Myr)
                                
                                line_x = (time_step[2][0] - SMBH_coords[2][0])
                                line_y = (time_step[2][1] - SMBH_coords[2][1])
                                line_z = (time_step[2][2] - SMBH_coords[2][2])
                                SMBH_dist[j] = np.log10(np.sqrt(line_x**2+line_y**2+line_z**2).value_in(units.pc))

                                vel_x = (time_step[3][0] - SMBH_coords[3][0])
                                vel_y = (time_step[3][1] - SMBH_coords[3][1])
                                vel_z = (time_step[3][2] - SMBH_coords[3][2])
                                vel = np.log10(np.sqrt(vel_x**2+vel_y**2+vel_z**2).value_in(units.kms))
                                pvel[j] = vel

                                SMBHx.append(SMBH_coords[2][0].value_in(units.pc))
                                SMBHy.append(SMBH_coords[2][1].value_in(units.pc))
                                SMBHz.append(SMBH_coords[2][2].value_in(units.pc))
                                if file_ == 0:
                                    SMBH_sample[0].append(SMBH_coords[2][0].value_in(units.pc))
                                    SMBH_sample[1].append(SMBH_coords[2][1].value_in(units.pc))
                                    SMBH_sample[2].append(SMBH_coords[2][2].value_in(units.pc))

                            ejec_KE /= ejec_KE[0]
                            ejec_ecc_SMBH = np.asarray(ejec_ecc_SMBH)

                            ejec_KE_smooth = plot_ini.moving_average(ejec_KE, smoothing)
                            ejec_ecc_SMBH_smooth = plot_ini.moving_average(ejec_ecc_SMBH, smoothing)
                            pvel_smooth = plot_ini.moving_average(pvel, smoothing)
                            time_smooth = plot_ini.moving_average(time, smoothing)
                            nearest_smooth = plot_ini.moving_average(NN, smoothing)
                            SMBH_dist_smooth = plot_ini.moving_average(SMBH_dist, smoothing)

                            fig = plt.figure(figsize=(12.5, 10))
                            ax1 = fig.add_subplot(221)
                            ax2 = fig.add_subplot(222)
                            ax3 = fig.add_subplot(223)
                            ax4 = fig.add_subplot(224)

                            ax1.set_title('Time vs. Eccentricity')
                            ax2.set_title('Time vs. Kinetic Energy')
                            ax3.set_title('Time vs. Inclination')
                            ax4.set_title('Time vs. Distance to Nearest Neighbour')
                            ax1.set_ylabel(r'$\log_{10}(1-e)$')
                            ax2.set_ylabel(r'$K_E / K_{E,0}$ [J]')
                            ax3.set_ylabel(r'$\log_{10} |v|$ [km s$^{-1}$]')
                            ax4.set_ylabel(r'$r_{nn}$ [pc]')
                            ax1.set_ylim(-6.5, 0)
                            for ax_ in [ax1, ax2, ax3, ax4]:
                                ax_.set_xlabel('Time [Myr]')
                                plot_ini.tickers(ax_, 'plot')

                            ax1.plot(time, ejec_ecc_SMBH, color = 'black', alpha = 0.35, linestyle = ':')
                            ax1.plot(time_smooth, ejec_ecc_SMBH_smooth, color = 'black')
                            ax2.plot(time, ejec_KE, color = 'black', alpha = 0.35, linestyle = ':')
                            ax2.plot(time_smooth, ejec_KE_smooth, color = 'black')
                            ax3.plot(time, pvel, color = 'black', alpha = 0.35, linestyle = ':')
                            ax3.plot(time_smooth, pvel_smooth, color = 'black')
                            ax4.plot(time, SMBH_dist, color = 'black', alpha = 0.35, linestyle = ':')
                            ax4.plot(time_smooth, SMBH_dist_smooth, color = 'black', label = 'w.r.t SMBH')
                            ax4.plot(time, NN, color = 'purple', alpha = 0.35, linestyle = ':')
                            ax4.plot(time_smooth, nearest_smooth, color = 'purple', label = 'w.r.t NN')

                            ax4.legend()
                            if focus_idx == 0:
                                plt.savefig('figures/system_evolution/'+str(int_)+'_Merger/ejec_bin_trip_evol_'+str(file_)
                                            +'_pop_'+str(np.shape(ptracker)[0])+'_'+str(res)+'_SMBH_merged.pdf', dpi=300, bbox_inches='tight')
                            else:
                                plt.savefig('figures/system_evolution/'+str(int_)+'_Merger/ejec_bin_trip_evol_'+str(file_)
                                            +'_pop_'+str(np.shape(ptracker)[0])+'_'+str(res)+'.pdf', dpi=300, bbox_inches='tight')
                            plt.clf()

                            ejec_ecc.append(ejec_ecc_SMBH)
                            distSMBH.append(SMBH_dist)
                            ejec_time.append(time)
                            vel_arr.append(pvel)
                            NNdist.append(NN)
                            ejec_ecc_smooth.append(ejec_ecc_SMBH_smooth)
                            distSMBH_smooth.append(SMBH_dist_smooth)
                            ejec_time_smooth.append(time_smooth)
                            vel_smooth.append(pvel_smooth)
                            NNdist_smooth.append(nearest_smooth)
                            if res == 'merger':
                                if len(ejec_ecc_SMBH[ejec_ecc_SMBH < 10**-6]) != 0:
                                    high_ecc += 1
                                    if max(time) < 5:
                                        cropped_idx.append(filter-1)
                                    all_idx.append(filter-1)
                                ejec_ecc_arr[iter].append(ejec_ecc_SMBH)
                                ejec_distSMBH_arr[iter].append(SMBH_dist)
                                ejec_vel_arr[iter].append(pvel)
                                ejec_NNdist_arr[iter].append(NN)
                        
        save_file = ['All', 'Cropped']
        idx = [all_idx, cropped_idx]
        c = colour_picker()
        for j in range(2):        
            colour_iter = 0
            fig = plt.figure(figsize=(12.5, 10))
            ax1 = fig.add_subplot(221)
            ax2 = fig.add_subplot(222)
            ax3 = fig.add_subplot(223)
            ax4 = fig.add_subplot(224) 

            ax1.set_title('Time vs. Eccentricity')
            ax2.set_title('Time vs. Distance to SMBH')
            ax3.set_title('Time vs. Velocity')
            ax4.set_title('Time vs. Distance to Nearest Neighbour')
            for ax_ in [ax1, ax2, ax3, ax4]:
                ax_.set_xlabel('Time [Myr]')
                plot_ini.tickers(ax_, 'plot') 
            for i in idx[j]:
                if colour_iter > 10:
                    colour_iter = 0
                ax1.plot(ejec_time[i], ejec_ecc[i], color = c[colour_iter], alpha = 0.35, linestyle = ':', zorder = 1)
                ax1.plot(ejec_time_smooth[i], ejec_ecc_smooth[i], color = c[colour_iter], zorder = 2)
                ax2.plot(ejec_time[i], distSMBH[i], color = c[colour_iter], alpha = 0.35, linestyle = ':', zorder = 1)
                ax2.plot(ejec_time_smooth[i], distSMBH_smooth[i], color = c[colour_iter], zorder = 2)
                ax3.plot(ejec_time[i], vel_arr[i], color = c[colour_iter], alpha = 0.35, linestyle = ':', zorder = 1)
                ax3.plot(ejec_time_smooth[i], vel_smooth[i], color = c[colour_iter], zorder = 2)
                ax4.plot(ejec_time[i], NNdist[i], color = c[colour_iter], alpha = 0.35, linestyle = ':', zorder = 1)
                ax4.plot(ejec_time_smooth[i], NNdist_smooth[i], color = c[colour_iter], zorder = 2)
                colour_iter += 1

            colour_iter = 0
            for i in idx[j]:
                if colour_iter > 10:
                    colour_iter = 0
                ax1.scatter(ejec_time[i][-1], ejec_ecc[i][-1], color = c[colour_iter], edgecolors='black', zorder = 3)
                ax2.scatter(ejec_time[i][-1], distSMBH[i][-1], color = c[colour_iter], edgecolors='black', zorder = 3)
                ax2.scatter(ejec_time[i][distSMBH[i] == min(distSMBH[i])], min(distSMBH[i]), color = c[colour_iter], s = 5, zorder = 3)
                ax3.scatter(ejec_time[i][-1], vel_arr[i][-1], color = c[colour_iter], edgecolors='black', zorder = 3)
                ax3.scatter(ejec_time[i][vel_arr[i] == max(vel_arr[i])], max(vel_arr[i]), color = c[colour_iter], s = 5, zorder = 3)
                ax4.scatter(ejec_time[i][-1], NNdist[i][-1], color = c[colour_iter], edgecolors='black', zorder = 3)
                ax4.scatter(ejec_time[i][NNdist[i] == min(NNdist[i])], min(NNdist[i]), color = c[colour_iter], s = 5, zorder = 3)
                colour_iter += 1
            ax1.set_ylabel(r'$\log_{10}(1-e)$')
            ax2.set_ylabel(r'$\log_{10}r_{SMBH}$ [pc]')
            ax3.set_ylabel(r'$\log_{10}|v|$ [km s$^{-1}$]')
            ax4.set_ylabel(r'$\log_{10}r_{NN}$ [pc]')
            plt.savefig('figures/system_evolution/'+str(int_)+'_Merger/ejec_bin_trip_evol'+str(save_file[j])+'.pdf', dpi=300, bbox_inches='tight')
            plt.clf()

        with open('figures/system_evolution/output/'+str(int_)+'_frac_ecc.txt', 'w') as file:
            file.write(str(int_)+' has high eccentricity in '+str(high_ecc)+'/'+str(filter)+' of its merging simulations.')
        
    c_hist = ['red', 'blue']

    ejec_ecc_flat = [[ ], [ ]]
    ejec_distSMBH_flat = [[ ], [ ]]
    ejec_vel_flat = [[ ], [ ]]
    ejec_NNdist_flat = [[ ], [ ]]
    for j in range(2):
        for sublist in ejec_ecc_arr[j]:
            for item in sublist:
                ejec_ecc_flat[j].append(item)
        for sublist in ejec_distSMBH_arr[j]:
            for item in sublist:
                ejec_distSMBH_flat[j].append(item)
        for sublist in ejec_vel_arr[j]:
            for item in sublist:
                ejec_vel_flat[j].append(item)
        for sublist in ejec_NNdist_arr[j]:
            for item in sublist:
                ejec_NNdist_flat[j].append(item)
    
    fig = plt.figure(figsize=(12.5, 10))
    ax1 = fig.add_subplot(221)
    ax2 = fig.add_subplot(222)
    ax3 = fig.add_subplot(223)
    ax4 = fig.add_subplot(224) 
    ax1.set_xlabel(r'$\log_{10}(1-e)$')
    ax2.set_xlabel(r'$\log_{10}r_{SMBH}$ [pc]')
    ax3.set_xlabel(r'$\log_{10}|v|$ [km s$^{-1}$]')
    ax4.set_xlabel(r'$\log_{10}r_{NN}$ [pc]')
    for ax_ in [ax1, ax2, ax3, ax4]:
        plot_ini.tickers(ax_, 'plot')
    for j in range(2):
        ejec_ecc_sort = np.sort(ejec_ecc_flat[j])
        ejec_ecc_index = np.asarray([i for i in enumerate(ejec_ecc_sort)])
        ax1.plot(ejec_ecc_sort, ejec_ecc_index[:,0]/ejec_ecc_index[-1,0], color = c_hist[j], label = integrator[j])

        ejec_distSMBH_sort = np.sort(ejec_distSMBH_flat[j])
        ejec_distSMBH_index = np.asarray([i for i in enumerate(ejec_distSMBH_sort)])
        ax2.plot(ejec_distSMBH_sort, ejec_distSMBH_index[:,0]/ejec_distSMBH_index[-1,0], color = c_hist[j])

        ejec_vel_sort = np.sort(ejec_vel_flat[j])
        ejec_vel_index = np.asarray([i for i in enumerate(ejec_vel_sort)])
        ax3.plot(ejec_vel_sort, ejec_vel_index[:,0]/ejec_vel_index[-1,0], color = c_hist[j])

        ejec_NNdist_sort = np.sort(ejec_NNdist_flat[j])
        ejec_NNdist_index = np.asarray([i for i in enumerate(ejec_NNdist_sort)])
        ax4.plot(ejec_NNdist_sort, ejec_NNdist_index[:,0]/ejec_NNdist_index[-1,0], color = c_hist[j])
    ax1.plot(ejec_ecc_sort, [10**(2*i) for i in ejec_ecc_sort], color = 'black', linestyle = ':', label = 'Thermal Distribution')
    ax1.legend(loc='upper left')
    plt.savefig('figures/system_evolution/merged_properties_cdf_histogram.pdf', dpi=300, bbox_inches='tight')

def energy_plotter(int_string):
    """
    Function to plot the energy evolution of the system

    output: Energy evolution plot of the system
    """

    plot_ini = plotter_setup()
    count = file_counter(int_string)
    etracker = file_opener('data/'+str(int_string)+'/spatial_plotters/energy/*')
    col_len = np.shape(etracker)[0]

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
            vals = etracker.iloc[i]
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

    plt.savefig('figures/system_evolution/etracker'+str(count)+'.pdf', dpi=300, bbox_inches='tight')
    plt.clf()
    plt.close()

    return

def global_properties():
    """
    Function which plots various Kepler elements of ALL particles simulated
    
    output: Plots of the inclination,semi-major axis, nearest neighbour and eccentricity of the 
            merged/ejected particle relative to the SMBH, nearest neighbour and second-nearest neighbour
    """
    
    plot_ini = plotter_setup()
    
    integrator = ['Hermite', 'GRX']
    ecc_arr = [[ ], [ ]]
    distSMBH_arr = [[ ], [ ]]
    vel_arr = [[ ], [ ]]
    NNdist_arr = [[ ], [ ]]
    KE_arr = [[ ], [ ]]

    iter = -1
    for int_ in integrator:   
        iter += 1
        data = natsort.natsorted(glob.glob('/media/erwanh/Elements/'+(int_)+'/particle_trajectory/*'))
        if int_ == 'GRX':
            energy = ['data/GRX/energy/'+str(i[47:]) for i in data]
            chaotic = ['data/GRX/no_addition/chaotic_simulation/'+str(i[47:]) for i in data]
        else:
            energy = ['data/Hermite/energy/'+str(i[51:]) for i in data]
            chaotic = ['data/Hermite/no_addition/chaotic_simulation/'+str(i[51:]) for i in data]

        SMBHx = []
        SMBHy = []
        SMBHz = []
        SMBH_sample = [[], [], []]
        
        for file_ in range(len(data)):
            with open(chaotic[file_], 'rb') as input_file:
                chaotic_tracker = pkl.load(input_file)
                if chaotic_tracker.iloc[0][6] <= 50:
                    with open(data[file_], 'rb') as input_file:
                        print('Reading File :', input_file)
                        ptracker = pkl.load(input_file)

                    with open(energy[file_], 'rb') as input_file:
                        etracker = pkl.load(input_file)

                    col_len = np.shape(ptracker)[1]-1
                    for parti_ in range(np.shape(ptracker)[0]):
                        if parti_ == 0:
                            pass
                        else:
                            particle = ptracker.iloc[parti_]
                            SMBH_data = ptracker.iloc[0]

                            time = np.empty(col_len)
                            neigh_key = np.empty((3, col_len))
                            NN = np.empty(col_len)
                            SMBH_dist = np.empty(col_len)
                            pvel = np.empty(col_len)
                            KE = np.empty(col_len)
                            ecc_SMBH = np.empty(col_len)

                            for j in range(col_len):
                                time_step = particle.iloc[j]
                                time_arr = etracker.iloc[j]
                                SMBH_coords = SMBH_data.iloc[j]

                                KE[j] = np.log10(time_step[4].value_in(units.J))
                                ecc_SMBH[j] = np.log10(1-time_step[8][0])

                                for i in range(3):
                                    neigh_key[i][j] = time_step[6][i]

                                NN[j] = np.log10(time_step[-1])
                                time[j] = time_arr[6].value_in(units.Myr)
                                
                                line_x = (time_step[2][0] - SMBH_coords[2][0])
                                line_y = (time_step[2][1] - SMBH_coords[2][1])
                                line_z = (time_step[2][2] - SMBH_coords[2][2])
                                SMBH_dist[j] = np.log10(np.sqrt(line_x**2+line_y**2+line_z**2).value_in(units.pc))

                                vel_x = (time_step[3][0] - SMBH_coords[3][0])
                                vel_y = (time_step[3][1] - SMBH_coords[3][1])
                                vel_z = (time_step[3][2] - SMBH_coords[3][2])
                                pvel[j] = np.log10(np.sqrt(vel_x**2+vel_y**2+vel_z**2).value_in(units.kms))

                                SMBHx.append(SMBH_coords[2][0].value_in(units.pc))
                                SMBHy.append(SMBH_coords[2][1].value_in(units.pc))
                                SMBHz.append(SMBH_coords[2][2].value_in(units.pc))
                                if file_ == 0:
                                    SMBH_sample[0].append(SMBH_coords[2][0].value_in(units.pc))
                                    SMBH_sample[1].append(SMBH_coords[2][1].value_in(units.pc))
                                    SMBH_sample[2].append(SMBH_coords[2][2].value_in(units.pc))

                            ecc_SMBH = np.asarray(ecc_SMBH)

                            ecc_arr[iter].append(ecc_SMBH)
                            distSMBH_arr[iter].append(SMBH_dist)
                            NNdist_arr[iter].append(NN)
                            vel_arr[iter].append(pvel)
                            KE_arr[iter].append(KE)
            
    c_hist = ['red', 'blue']

    ecc_flat = [[ ], [ ]]
    distSMBH_flat = [[ ], [ ]]
    vel_flat = [[ ], [ ]]
    NNdist_flat = [[ ], [ ]]
    KE_flat = [[ ], [ ]]
    for j in range(2):
        for sublist in ecc_arr[j]:
            for item in sublist:
                ecc_flat[j].append(item)
        for sublist in distSMBH_arr[j]:
            for item in sublist:
                distSMBH_flat[j].append(item)
        for sublist in vel_arr[j]:
            for item in sublist:
                vel_flat[j].append(item)
        for sublist in NNdist_arr[j]:
            for item in sublist:
                NNdist_flat[j].append(item)
        for sublist in KE_arr[j]:
            for item in sublist:
                KE_flat[j].append(item)
    
    fig = plt.figure(figsize=(12.5, 10))
    ax1 = fig.add_subplot(221)
    ax2 = fig.add_subplot(222)
    ax3 = fig.add_subplot(223)
    ax4 = fig.add_subplot(224) 
    ax1.set_xlabel(r'$\log_{10}(1-e)$')
    ax2.set_xlabel(r'$\log_{10}r$ [pc]')
    ax3.set_xlabel(r'$\log_{10}|v|$ [km s$^{-1}$]')
    ax4.set_xlabel(r'$\log_{10} K_E$ [J]')
    for ax_ in [ax1, ax2, ax3, ax4]:
        plot_ini.tickers(ax_, 'plot')
    for j in range(2):
        ecc_sort = np.sort(ecc_flat[j])
        ecc_index = np.asarray([i for i in enumerate(ecc_sort)])
        ax1.plot(ecc_sort, ecc_index[:,0]/ecc_index[-1,0], color = c_hist[j], label = integrator[j])

        distSMBH_sort = np.sort(distSMBH_flat[j])
        distSMBH_index = np.asarray([i for i in enumerate(distSMBH_sort)])
        NNdist_sort = np.sort(NNdist_flat[j])
        NNdist_index = np.asarray([i for i in enumerate(NNdist_sort)])
        if j == 0:
            ax2.plot(distSMBH_sort, distSMBH_index[:,0]/distSMBH_index[-1,0], color = c_hist[j], linestyle = ':', label = r'$r_{\rm{SMBH}}$')
            ax2.plot(NNdist_sort, NNdist_index[:,0]/NNdist_index[-1,0], color = c_hist[j], linestyle = '-.', label =r'$r_{\rm{IMBH}}$')
        else:
            ax2.plot(distSMBH_sort, distSMBH_index[:,0]/distSMBH_index[-1,0], color = c_hist[j], linestyle = ':')
            ax2.plot(NNdist_sort, NNdist_index[:,0]/NNdist_index[-1,0], color = c_hist[j], linestyle = '-.')
        
        vel_sort = np.sort(vel_flat[j])
        vel_index = np.asarray([i for i in enumerate(vel_sort)])
        ax3.plot(vel_sort, vel_index[:,0]/vel_index[-1,0], color = c_hist[j])

        KE_sort = np.sort(KE_flat[j])
        KE_index = np.asarray([i for i in enumerate(KE_sort)])
        ax4.plot(KE_sort, KE_index[:,0]/KE_index[-1,0], color = c_hist[j])

    ax1.plot(ecc_sort, [10**(2*i) for i in ecc_sort], color = 'black', linestyle = ':', label = 'Thermal Distribution')
    ax1.legend(loc='upper left')
    ax2.legend(loc='upper left')
    plt.savefig('figures/system_evolution/global_properties_N=20_cdf_histogram.pdf', dpi=300, bbox_inches='tight')

def spatial_plotter(int_string):
    """
    Function to plot the evolution of the system

    output: The spatial evolution of the system
    """

    plot_ini = plotter_setup()

    ptracker_files = natsort.natsorted(glob.glob('/media/erwanh/Elements/'+(int_string)+'/particle_trajectory/*'))
    etracker_files = natsort.natsorted(glob.glob('data/'+str(int_string)+'/energy/*'))
    ctracker_files = natsort.natsorted(glob.glob('data/'+str(int_string)+'/no_addition/chaotic_simulation/*'))
    iter_file = -1
    for file_ in range(len(ptracker_files)):
        iter_file += 1
        with open(ctracker_files[file_], 'rb') as input_file:
            print('Reading File ', file_, ' : ', input_file)
            ctracker = pkl.load(input_file)
            if ctracker.iloc[0][5] > 0:
                with open(ptracker_files[file_], 'rb') as input_file:
                    ptracker = pkl.load(input_file)
                with open(etracker_files[file_], 'rb') as input_file:
                    etracker = pkl.load(input_file)

                col_len_raw = np.shape(ptracker)[1]
                col_len = round(col_len_raw**0.6)
                parti_size = 20+len(ptracker)**-0.5

                line_x = np.empty((len(ptracker), col_len))
                line_y = np.empty((len(ptracker), col_len))
                line_z = np.empty((len(ptracker), col_len))

                for i in range(len(ptracker)):
                    tptracker = ptracker.iloc[i]
                    for j in range(col_len):
                        coords = tptracker.iloc[j][2]
                        if len(coords) == 1:
                            pass
                        else:
                            line_x[i][j] = coords[0].value_in(units.pc)
                            line_y[i][j] = coords[1].value_in(units.pc)
                            line_z[i][j] = coords[2].value_in(units.pc)

                time = np.empty((col_len_raw - 1))
                dE_array = np.empty((col_len_raw - 1))

                for i in range(col_len_raw):
                    if i == 0:
                        pass
                    else:
                        vals = etracker.iloc[i]
                        time[i-1] = vals[6].value_in(units.Myr)
                        dE_array[i-1] = vals[7]

                c = colour_picker()
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
                
                for i in range(len(ptracker)):
                    iter += 1
                    if iter > len(c):
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
                                    c = c[iter-2], edgecolors = 'black', s = parti_size, zorder = 3)
                        ax1.scatter(line_x[i]-line_x[0], line_y[i]-line_y[0], 
                                    c = c[iter-2], s = 1, zorder = 1) 

                        ax3.scatter(line_x[i][-1]-line_x[0][-1], line_z[i][-1]-line_z[0][-1], 
                                    c = c[iter-2], edgecolors = 'black', s = parti_size, zorder = 3)
                        ax3.scatter(line_x[i]-line_x[0], line_z[i]-line_z[0], 
                                    c = c[iter-2], s = 1, zorder = 1) 

                        ax4.scatter(line_y[i][-1]-line_y[0][-1], line_z[i][-1]-line_z[0][-1], 
                                    c = c[iter-2], edgecolors = 'black', s = parti_size, zorder = 3)
                        ax4.scatter(line_y[i]-line_y[0], line_z[i]-line_z[0], 
                                    c = c[iter-2], s = 1, zorder = 1) 
                        
                ax2.plot(time[:-5], dE_array[:-5], color = 'black')
                plt.savefig('figures/system_evolution/Overall_System/simulation_evolution_'+str(iter_file)+'.pdf', dpi=300, bbox_inches='tight')
                plt.clf()
                plt.close()     

                fig = plt.figure(figsize=(8, 8))
                ax3D = fig.add_subplot(121, projection="3d")
                iter = -1
                for i in range(len(ptracker)):
                    iter += 1
                    if iter > len(c):
                        iter = 0
                    if i == 0:
                        pass
                    else:
                        ax3D.scatter(line_x[i]-line_x[0], 
                                    line_y[i]-line_y[0], 
                                    line_z[i]-line_z[0], 
                                    c = c[iter-2], s = 1, zorder = 1)
                ax3D.scatter(0, 0, 0, color = 'black', s = 150, zorder = 2)
                ax3D.xaxis.set_major_formatter(mtick.FormatStrFormatter('%0.2f'))
                ax3D.yaxis.set_major_formatter(mtick.FormatStrFormatter('%0.2f'))
                ax3D.zaxis.set_major_formatter(mtick.FormatStrFormatter('%0.2f'))
                ax3D.set_xlabel(r'$x$ [pc]')
                ax3D.set_ylabel(r'$y$ [pc]')
                ax3D.set_zlabel(r'$z$ [pc]')
                ax3D.view_init(30, 160)
                plt.savefig('figures/system_evolution/Overall_System/simulation_evolution_3D_'+str(iter_file)+'.pdf', dpi=300, bbox_inches='tight')

    return