
class vejec_mass(object):
    """
    Class to plot the indep_variable vs. ejection velocity plots
    """
    def __init__(self, int_string = 'Hermite'):
        """
        Extracts the required data
        """

        self.IMBH_tracker = bulk_stat_extractor('data/'+str(int_string)+'/vej_data/particle_trajectory/*', 'Y')
        self.ejec_data = bulk_stat_extractor('data/'+str(int_string)+'/vej_data/no_addition/chaotic_simulation/*', 'N')

        self.ex = np.empty((len(self.ejec_data)))
        self.ey = np.empty((len(self.ejec_data)))
        self.ez = np.empty((len(self.ejec_data)))

        self.ejec_vx = np.empty((len(self.ejec_data)))
        self.ejec_vy = np.empty((len(self.ejec_data)))
        self.ejec_vz = np.empty((len(self.ejec_data)))

        self.ejec_KE = np.empty((len(self.ejec_data)))
        self.ejec_PE = np.empty((len(self.ejec_data)))
        self.incl = np.empty((len(self.ejec_data)))
        self.Nclose = np.empty((len(self.ejec_data)))

        self.tot_mass = np.empty((len(self.ejec_data)))
        self.tot_pop = np.empty((len(self.ejec_data)))
        self.surv_time = np.empty((len(self.ejec_data)))

        for i in range(len(self.ejec_data)):
            self.ex[i], self.ey[i], self.ez[i], self.ejec_vx[i], self.ejec_vy[i], \
            self.ejec_vz[i], self.ejec_KE[i], self.ejec_PE[i], self.incl[i], \
            self.Nclose[i] = ejected_extract_final(self.IMBH_tracker[i], self.ejec_data[i])
            
            self.vals_df = self.ejec_data[i].iloc[0]          
            self.tot_mass[i] = np.sum(self.vals_df[8].value_in(units.MSun))
            self.tot_pop[i] = self.vals_df[6]
            self.surv_time[i] = self.vals_df[-2].value_in(units.Myr)

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
        colour_axes = ax.scatter(xdata, ydata, c = cdata, zorder = 3)
        plt.colorbar(colour_axes, ax=ax, label = clabel)
        plot_ini.tickers(ax)
        ax.set_xlabel(r'Total IMBH Mass [$\frac{M}{10^3 M_{\odot}}$]')
        ax.set_ylabel(r'Ejection Velocity [km/s]')
        return ax

    def vejec_sysmass(self):
        """
        Function to plot how the total initial mass of the system influences the ejection velocity
        """

        in_mass = np.unique(self.tot_mass[self.tot_mass < 11*10**3])
        avg_vel = np.empty((len(in_mass)))
        avg_surv = np.empty((len(in_mass)))

        iter = -1
        popul = [ ]
        for mass_ in in_mass:
            iter += 1
            indices = np.where((self.tot_mass == mass_))[0]
            temp_evel = [self.ejec_vel[i] for i in indices]
            temp_surv = [self.surv_time[i] for i in indices]
            temp_pop  = [self.tot_pop[i] for i in indices]
            popul.append(np.unique(temp_pop))
            avg_vel[iter] = np.mean(temp_evel)
            avg_surv[iter] = np.mean(temp_surv)

        popul = np.asarray(popul).flatten()
        avg_vel = np.asarray(avg_vel)
        yint, slope = polyfit(popul, avg_vel, 1)
        xlist = np.linspace(2.5,10.5)
        ylist = [yint+slope*i for i in xlist]

        ax = self.plotter_func(popul, avg_vel, avg_surv, r'Average Ejection Time [Myr]')
        ax.plot(xlist, ylist, linestyle = '--', zorder = 1)
        ax.text(3, 50, r'$v_{ejec} = '+str('{:.2f}'.format(slope))+r'$N +'+str('{:.2f}'.format(yint)))
        ax.set_ylim(0, 600)
        ax.set_xlim(2.5, 10.5)
        plt.savefig('figures/mean_vej.pdf', dpi=300, bbox_inches='tight')

        """ax = self.plotter_func(popul, avg_vel, avg_surv, r'Average Ejection Time [Myr]')
        ax.set_ylim(0, 600)
        ax.set_xlim(10.5, 35.5)
        plt.savefig('figures/mean_vej_highm.pdf', dpi=300, bbox_inches='tight')"""

    def vejec_syspop(self):
        """
        Plot to show how the total population of the system influences the ejection velocity
        """

        idx = np.where(self.tot_mass < 11 * 10**3)[0]
        ejec_vel = np.asarray(self.ejec_vel)[idx]
        mass = self.tot_mass[idx]
        
        pop_eject = np.where(ejec_vel > 690)[0]
        print('Simulations where v > vej MW:', len(pop_eject), '/', len(ejec_vel))

        ax = self.plotter_func(self.tot_pop[idx], 
                               ejec_vel, 
                               (mass), r'Total IMBH Mass [$\frac{M}{10^3 M_{\odot}}$]')
        plt.axhline(y=690, linestyle = '-.', color = 'black', zorder = 1)
        ax.text(2.1, 705, r'$v_{esc, MW}$')
        ax.set_ylim(0, 1500)
        ax.set_xlim(2, 11)
        plt.savefig('figures/scatter_vej.pdf', dpi=300, bbox_inches='tight')

        """plt.axhline(y=690, linestyle = '-.', color = 'black', zorder = 1)
        ax.text(10.6, 705, r'$v_{esc, MW}$')
        ax.set_ylim(0, 1500)
        ax.set_xlim(10.5, 35.5)
        plt.savefig('figures/scatter_vej_highmass.pdf', dpi=300, bbox_inches='tight')"""

    def finalpos_energy_data_extract(self, filter):
        """
        Function to extract particle energy data
        
        Inputs:
        filter: String 'B' or None which decides whether you are looking (B)elow a certain threshold or not
        """
        
        if filter == 'B':
            idx = np.where(self.ejec_KE < 1e45)
            ejected_KE1 = np.asarray([i for i in self.ejec_KE[idx]])
            ejected_PE1 = np.asarray([-i for i in self.ejec_PE[idx]])
            final_pos = np.asarray([(i**2+j**2+z**2)**0.5 for i, j, z in zip(self.ex, self.ey, self.ez)])[idx]

            idx = np.where(abs(ejected_PE1) < 1e45)
            ejected_KE = np.asarray([i/max(ejected_KE1[idx]) for i in ejected_KE1[idx]])
            ejected_PE = np.asarray([i/max(ejected_KE1[idx]) for i in ejected_PE1[idx]])
            final_pos = np.asarray(final_pos[idx])

        else:
            idx = np.where(self.ejec_KE > -0.01)
            ejected_KE = [i/max(self.ejec_KE[idx]) for i in self.ejec_KE[idx]]
            ejected_PE = [-i/max(self.ejec_KE[idx]) for i in self.ejec_PE[idx]]
            final_pos = np.asarray([(i**2+j**2+z**2)**0.5 for i, j, z in zip(self.ex, self.ey, self.ez)])[idx]

        linex = np.linspace(-0.01*np.min(ejected_KE), 10)
        liney = -1*linex

        return ejected_KE, ejected_PE, final_pos, linex, liney

    def finalpos_energy_plotter(self, KE, PE, linex, liney, cdata, clabel, no_axis, save_file, hist):
        """
        Function to plot the required data
        
        Inputs:
        KE/PE:       The x, y data of the plot
        linex/liney: The line of best fit to distinguish bound from unbound
        cdata:       The data for which the colour code is based off
        no_axis:     The number of axis plotting (zoom or not)
        save_file:   String to ensure no overwritten plots
        hist:        String to state whether plotting histogram or other coloured variable
        """

        plot_ini = plotter_setup()
        plt.figure(figsize=(10, 4))
        if no_axis == 2:
            gs = gridspec.GridSpec(1, 2)
            ax1 = plt.subplot(gs[0,0])
            ax2 = plt.subplot(gs[0,1])
            ax1.set_ylabel(r'$E_P/E_{{K, {}}}$'.format(str('max')))
           
            for ax_ in [ax1, ax2]:
                ax_.plot(linex, liney, color = 'black', linestyle = '-.')
                ax_.set_xlabel(r'$E_K/E_{{K, {}}}$'.format(str('max')))
                plot_ini.tickers(ax_)
                
            if hist == 'Y':
                bin2d_sim, xedges_s, yedges_s, image = ax2.hist2d(KE, PE, bins=(200, 200), range=([0,1.1],[1.1*min(PE),0]))
                bin2d_sim /= np.max(bin2d_sim)
                extent = [0, 1.1, 1.1*min(PE), 0]

                contours = ax2.imshow(np.log10(bin2d_sim), extent = extent, aspect='auto', origin = 'upper')
                #ax2.scatter(KE, PE, s = 15, color = 'red', edgecolors = 'black')
                #ax2.set_ylim(m(PE),0)
                ax2.set_xlim(0,max(KE))
                bin2d_sim, xedges_s, yedges_s, image = ax1.hist2d(KE, PE, bins=(80, 80), range=([0,0.35],[-0.35,0]))
                bin2d_sim /= np.max(bin2d_sim)
                extent = [0, 0.3, -0.3,0]
                contours = ax1.imshow(np.log10(bin2d_sim), extent = extent, aspect='auto', origin = 'upper')
                #ax1.scatter(KE, PE, s = 15, color = 'red', edgecolors = 'black')
                cbar = plt.colorbar(contours, ax=ax2)
                cbar.set_label(label = r'$\log_{10}(N/N_{max})$',rotation=270,labelpad=15)

                ax1.set_xlim(0,0.3)
                ax1.set_ylim(-0.3,0)

            if hist == 'N':
                colour_axes = ax1.scatter(KE, PE, c=cdata, s=2)
                colour_axes = ax2.scatter(KE, PE, c=cdata, s=2)
                plt.colorbar(colour_axes, ax=ax2, label = clabel)
                ax2.set_ylim(-1.05,0)
                ax2.set_xlim(0,1.05)
                ax1.set_xlim(0,0.3)
                ax1.set_ylim(-0.3,0)
                
            plt.savefig('figures/energy_diagram_'+str(save_file)+'_.pdf', dpi=500, bbox_inches='tight')
            return
        
        else:
            if hist == 'N':
                fig, ax = plt.subplots()
                ax.set_title('Ejected Particles \nKinetic Energy vs. Potential Energy')
                ax.plot(linex, liney, color = 'black', linestyle = '-.')
                colour_axes = ax.scatter(KE, PE, c = cdata, s= 2)
                plt.colorbar(colour_axes, ax = ax, label = r'Final Distance to Core [pc]')
                plt.plot(linex, liney, color = 'black', linestyle = '-.')
                plt.xlim(1e-5,1.05)
                plt.ylim(-1, -10**-15)
                #plt.yscale('symlog')
                #plt.xscale('log')
                plt.xlabel(r'$E_K/E_{{K, {}}}$'.format(str('max')))
                plt.ylabel(r'$E_P/E_{{K, {}}}$'.format(str('max')))
                plot_ini.tickers(ax)
                plt.savefig('figures/energy_diagram_'+str(save_file)+'.pdf', dpi=300, bbox_inches='tight')
                plt.clf()
            else: 
                fig, ax = plt.subplots()

                bin2d_sim, xedges_s, yedges_s, image = ax.hist2d(KE, PE, bins=(100,300), range=([np.min(KE), 2],[1.1*min(PE),0.1]))
                bin2d_sim /= np.max(bin2d_sim)
                extent = [0, 1.1, 1.1*min(PE),0]

                ax.set_title('Ejected Particles \nKinetic Energy vs. Potential Energy')
                ax.set_xlim(1e-6, 1.1)
                ax.plot(linex, liney, color = 'black', linestyle = '-.')
                ax.set_ylim(-1, -10**-6)
                ax.set_xlabel(r'$E_K/E_{{K, {}}}$'.format(str('max')))
                ax.set_ylabel(r'$E_P/E_{{K, {}}}$'.format(str('max')))
                #ax.set_yscale('symlog')
                #ax.set_xscale('log')
                contours = ax.imshow(np.log10(bin2d_sim), extent = extent, aspect='auto', origin = 'upper')
                ax.plot(linex, liney, color = 'black', linestyle = '-.')
                #ax.scatter(KE, PE, s=15, color = 'red', edgecolors ='black')
                cbar = plt.colorbar(contours, ax=ax)
                cbar.set_label(r"$\log_{10}(N/N_{max})$", rotation=270,labelpad=15)
                
                ax.set_yscale('symlog')
                plot_ini.tickers(ax)
                plt.savefig('figures/energy_diagram_'+str(save_file)+'.pdf', dpi=300, bbox_inches='tight')
                plt.clf()
            return

    def finalpos_energy(self):
            """
            Function to plot the KE vs. PE plot coloured with final positions
            """

            save_file = ['fpos', 'fpos_crop']
            data_filt = [None, 'B']
            for i in range(2):
                axis = i+1
                ejected_KE, ejected_PE, final_pos, unbounded_x, unbounded_y = self.finalpos_energy_data_extract(data_filt[i])
                self.finalpos_energy_plotter(ejected_KE, ejected_PE, unbounded_x, unbounded_y, final_pos, 'Final Distance to Core [pc]', axis, save_file[i], 'N')     

    def finalpos_energy_hist(self):
            """
            Function to plot the KE vs. PE plot coloured with final positions
            """
            ejected_KE, ejected_PE, final_pos, unbounded_x, unbounded_y = self.finalpos_energy_data_extract('B')
            self.finalpos_energy_plotter(ejected_KE, ejected_PE, unbounded_x, unbounded_y, final_pos, r'$\log_{10}n$', 2, 'histogram_crop', 'Y')    

            ejected_KE, ejected_PE, final_pos, unbounded_x, unbounded_y = self.finalpos_energy_data_extract(None)
            self.finalpos_energy_plotter(ejected_KE, ejected_PE, unbounded_x, unbounded_y, final_pos, r'$\log_{10}n$', 1, 'histogram_All', 'Y')   

    def finalpos_incl_plotter(self):
            """
            Function to plot the KE vs. PE plot coloured with final positions
            """

            dist = np.sqrt(self.ex**2+self.ey**2+self.ez**2)
            plot_ini = plotter_setup()

            for i in range(2):
                fig, ax = plt.subplots()
                plt.xlabel(r'Distance to Core [pc]')
                plt.ylabel(r'Inclination [deg]'.format(str('max')))
                ax.set_ylim(-95,95)

                if i == 0:
                    colour_axes = ax.scatter(dist, self.incl, c = self.Nclose, s= 2)
                    plot_ini.tickers(ax)
                    ax.set_xscale('log')
                    plt.title('Ejected Particles \nFinal Distance vs. Orbital Inclination')
                    plt.colorbar(colour_axes, ax = ax, label = r'Number of Close Encounters')
                    plt.savefig('figures/inclination.pdf', dpi=300, bbox_inches='tight')
                    plt.clf()

                if i == 1:
                    bin2d_sim, xedges_s, yedges_s, image = ax.hist2d(dist, self.incl, bins = (60,60), range = ([0.9*min(dist),4],[-90,90]))
                    plot_ini.tickers(ax)
                    #ax.set_xscale('log')
                    bin2d_sim /= np.max(bin2d_sim)
                    extent = [0.9*min(dist), 5, -90, 90]
                    contours = ax.imshow(np.log10(bin2d_sim), extent = extent, aspect='auto', origin = 'upper')
                    #ax.scatter(dist, self.incl, s = 15, color = 'red', edgecolors = 'black')
                    plt.title('Ejected Particles \nFinal Distance vs. Orbital Inclination')
                    plt.colorbar(contours, ax=ax, label = r'$\log_{10}(N/N_{max})$')
                    plt.savefig('figures/inclination_histogram.pdf', dpi=300, bbox_inches='tight')
                    plt.clf()

class event_tracker(object):
    """
    Class to take stats of ejection vs. mergers
    """
    def __init__(self, int_string = 'Hermite'):
        self.merge_events = bulk_stat_extractor('data/'+str(int_string)+'/collision_events/*')

        merge_events = np.asarray(self.merge_events)
        no_merge = 0 
        for i in range(len(merge_events)):
            if merge_events[i][0][2].number != 0:
                no_merge += 1
        print('Total Simulations: ', len(merge_events),
              '\nNumber of mergers: ', no_merge)