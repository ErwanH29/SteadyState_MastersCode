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

def spat_limits(xcoords, ycoords):
    """
    Function which defines the x, y limits for spatial plots
    
    Inputs:
    spatial_data:  The positional data of the various IMBH tracked
    outputs:       The x, y limits of the spatial plots
    """
    
    x_lim1 = abs(spatial_data[:][:,0].min())
    x_lim2 = abs(spatial_data[:][:,0].max())
    y_lim1 = abs(spatial_data[:][:,1].min())
    y_lim2 = abs(spatial_data[:][:,1].max())
    xy_choice = 1.2*max(x_lim1, x_lim2, y_lim1, y_lim2)

    return xy_choice

def spatial_plotter():
    """
    Function to plot the evolution of the system
    outputs: Plots of the evolution of the system in Cartesian coordinates
    """

    SMBH_code = MW_SMBH()
    filename = glob.glob('data/center_of_mass/*')
    with open(os.path.join(max(filename, key=os.path.getctime)), 'rb') as input_file:
        com_tracker = pkl.load(input_file)

    filename = glob.glob('data/positions/*')
    with open(os.path.join(max(filename, key=os.path.getctime)), 'rb') as input_file:
        spatial_tracker = pkl.load(input_file)

    for i in range(len(spatial_tracker)):
        IMBH_tracker = spatial_tracker.iloc[i]
        col_len = len(IMBH_tracker)-1

        line_x = np.empty((len(spatial_tracker)+1, col_len, 1))
        line_y = np.empty((len(spatial_tracker)+1, col_len, 1))
        line_z = np.empty((len(spatial_tracker)+1, col_len, 1))
        
        for j in range(col_len-1):
            coords = IMBH_tracker.iloc[j+1][0] # * (1 | units.m)**-1
            line_x[i][j][0] = coords[0].value_in(units.AU)
            line_y[i][j][0] = coords[1].value_in(units.AU)
            line_z[i][j][0] = coords[2].value_in(units.AU)
        print(line_x[0])
        print(i)

    fig = plt.figure(figsize=(12, 4))
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)
    colors = get_distinct(len(spatial_tracker))
    ax1.set_title('Close-Encounters')
    ax2.set_title('Overall System')

    xy_lim= spat_limits(spatial_tracker)
    ax1.set_xlabel(r'$x$-Coordinate [AU]')
    ax1.set_ylabel(r'$y$-Coordinate [AU]')
    ax2.set_xlabel(r'$x$-Coordinate [AU]')
    ax2.set_ylabel(r'$y$-Coordinate [AU]')
    #ax2.yaxis.set_major_formatter(mtick.FormatStrFormatter('%1.0e'))
    #ax2.xaxis.set_major_formatter(mtick.FormatStrFormatter('%1.0e'))

    ax2.set_xlim(-xy_lim.value_in(units.AU), xy_lim.value_in(units.AU))
    ax2.set_ylim(-xy_lim.value_in(units.AU), xy_lim.value_in(units.AU))
    ax2.plot(spatial_tracker[-1][:,0].value_in(units.AU), 
             spatial_tracker[-1][:,1].value_in(units.AU), 
                c = 'black', linestyle = '--')

    lim_x = [ ]
    lim_y = [ ]

    for i in range(len(spatial_tracker)-1):
        ax1.plot(spatial_tracker[i][:,0].value_in(units.AU)-com_tracker[0][:,0].value_in(units.AU), 
                 spatial_tracker[i][:,1].value_in(units.AU)-com_tracker[0][:,1].value_in(units.AU), 
                 c = colors[i], lw = 1)
        ax2.plot(spatial_tracker[i][:,0].value_in(units.AU), 
                 spatial_tracker[i][:,1].value_in(units.AU), 
                 c = colors[i])
    min_x = min(spatial_tracker[i][:,0].value_in(units.AU)-com_tracker[0][:,1].value_in(units.AU))
    max_x = max(spatial_tracker[i][:,0].value_in(units.AU)-com_tracker[0][:,1].value_in(units.AU))
    lim_x.append(max(abs(min_x), abs(max_x)))

    min_y = min(spatial_tracker[i][:,1].value_in(units.AU)-com_tracker[0][:,1].value_in(units.AU))
    max_y = max(spatial_tracker[i][:,1].value_in(units.AU)-com_tracker[0][:,1].value_in(units.AU))
    lim_y.append(max(abs(min_y), abs(max_y)))

    ax1.set_xlim(max(-2, -(1.1*np.max(lim_x))), 
                 min(2,    1.1*np.max(lim_x)))   
    ax1.set_ylim(max(-2, -(1.1*np.max(lim_y))),
                 min(2,    1.1*np.max(lim_y)))


    for i in range(len(spatial_tracker)-1):
        ax1.scatter(spatial_tracker[i][0,0].value_in(units.AU)-com_tracker[0][0,0].value_in(units.AU),
                        spatial_tracker[i][0,1].value_in(units.AU)-com_tracker[0][0,1].value_in(units.AU),
                        c = colors[i], edgecolors = 'black')


    ax1.set_ylim()
    ax2.scatter(SMBH_code.position.x.value_in(units.AU), 
                SMBH_code.position.y.value_in(units.AU), 
                color = 'black', s = 500*SMBH_code.bh_rad.value_in(units.AU) )
    plt.savefig('figures/spatial_tracker'+str(datetime.now())+'.pdf', dpi=300, bbox_inches='tight')
    plt.clf()
    plt.close()
    return

def energy_plotter():
    """
    Function to plot the evolution of the system
    outputs: Plots of the evolution of the system in Cartesian coordinates
    """

    filename = glob.glob('data/energy/*')
    with open(os.path.join(max(filename, key=os.path.getctime)), 'rb') as input_file:
        energy_tracker = pkl.load(input_file)

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
        
    ax1.plot(energy_tracker[0][0], energy_tracker[1][0])
    ax2.plot(energy_tracker[0][0][15:], energy_tracker[2][0])
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
    outputs: Movie of the animation"""

    print('You have chosen to animate. This will take a while...')

    plt.clf()
    SMBH_code = MW_SMBH()

    filename = glob.glob('data/center_of_mass/*')
    with open(os.path.join(max(filename, key=os.path.getctime)), 'rb') as input_file:
        com_tracker = pkl.load(input_file)

    filename_energy = glob.glob('data/energy/*')
    with open(os.path.join(max(filename_energy, key=os.path.getctime)), 'rb') as input_file:
        energy_tracker = pkl.load(input_file)

    filename_pos = glob.glob('data/positions/*')
    with open(os.path.join(max(filename_pos, key=os.path.getctime)), 'rb') as input_file:
        spatial_tracker = pkl.load(input_file)

    spatial_tracker[spatial_tracker.number==0] = np.nan | units.parsec

    colors = get_distinct(len(spatial_tracker))
    numData = len(energy_tracker[0][0])

    fig = plt.figure(figsize=(12.5, 8))
    
    ax1 = fig.add_subplot(221)
    ax2 = fig.add_subplot(222)
    ax3 = fig.add_subplot(223)
    ax4 = fig.add_subplot(224)

    def animate(num):
        if num %100 == 0:
            print('Current frame: ', num)

        xy_lim= spat_limits(spatial_tracker)
        ax1.clear()
        ax2.clear()
        ax3.clear()
        ax4.clear()

        for i in range(len(spatial_tracker)-1):
            ax1.plot(spatial_tracker[i][0:num+1,0].value_in(units.AU),
                     spatial_tracker[i][0:num+1,1].value_in(units.AU),
                     c = colors[i])

            ax3.plot(spatial_tracker[i][max(0,num+1-70):num+1,0].value_in(units.AU)-com_tracker[0][max(0,num+1-70):num+1,0].value_in(units.AU),
                     spatial_tracker[i][max(0,num+1-70):num+1,1].value_in(units.AU)-com_tracker[0][max(0,num+1-70):num+1,1].value_in(units.AU),
                     c = colors[i], lw = 2)
            ax3.scatter(spatial_tracker[i][num,0].value_in(units.AU)-com_tracker[0][num,0].value_in(units.AU),
                        spatial_tracker[i][num,1].value_in(units.AU)-com_tracker[0][num,1].value_in(units.AU),
                        c = colors[i], edgecolors = 'black', s = 40)

            ax4.plot(spatial_tracker[i][max(0,num+1-70):num+1,0].value_in(units.AU)-com_tracker[0][max(0,num+1-70):num+1,0].value_in(units.AU),
                     spatial_tracker[i][max(0,num+1-70):num+1,2].value_in(units.AU)-com_tracker[0][max(0,num+1-70):num+1,2].value_in(units.AU),
                     c = colors[i], lw = 2)
            ax4.scatter(spatial_tracker[i][num,0].value_in(units.AU)-com_tracker[0][num,0].value_in(units.AU),
                        spatial_tracker[i][num,2].value_in(units.AU)-com_tracker[0][num,2].value_in(units.AU),
                        c = colors[i], edgecolors = 'black', s = 40)

        ax1.plot(spatial_tracker[-1][0:num+1,0].value_in(units.AU),
                 spatial_tracker[-1][0:num+1,1].value_in(units.AU),
                 c = 'black', linestyle = '--')
        ax1.set_title(str("{:.3f}".format(energy_tracker[0][0][num]*tend.number*eta*365)+" Days"), loc = 'left')
        ax1.scatter(SMBH_code.position.x.value_in(units.AU), 
                    SMBH_code.position.y.value_in(units.AU), 
                    color = 'black', s = SMBH_code.bh_mass.number**0.3)
        ax1.set_xlim(-xy_lim.value_in(units.AU), xy_lim.value_in(units.AU))
        ax1.set_ylim(-xy_lim.value_in(units.AU), xy_lim.value_in(units.AU))        
        ax1.set_xlabel(r'$x$-Coordinate [AU]')
        ax1.set_ylabel(r'$y$-Coordinate [AU]')

        ax2.set_xlabel(r'Time Step [$\eta$]')
        ax2.set_ylabel(r'$\frac{|E(t)-E_0|}{|E_0|}$')
        ax2.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.0e'))
        ax2.plot(energy_tracker[0][0][:num+1], energy_tracker[1][0][:num+1])
        ax2.set_yscale('log')

        ax3.set_xlabel(r'$x$-Coordinate [AU]')
        ax3.set_ylabel(r'$y$-Coordinate [AU]')

        ax4.set_xlabel(r'$x$-Coordinate [AU]')
        ax4.set_ylabel(r'$z$-Coordinate [AU]')

    anim = animation.FuncAnimation(fig, animate, interval = 100, frames=numData)
    #plt.show()
    anim.save('figures/animation'+str(datetime.now())+'.gif', writer='pillow')
    return 


spatial_plotter()