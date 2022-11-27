from amuse.lab import *
from file_logistics import *
from spatial_plotters import *
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import quad
from matplotlib.ticker import FormatStrFormatter
from scipy.interpolate import interp1d

def calc_loop(iter, redshift, range_a, range_b, phi0, phi0_diff, mc, mc_diff, alpha, alpha_diff, phen):
    """
    Function to calculate the events per cubic Gpc per year for both a constant GC formation
    and a redshift-dependent one
    
    Inputs:
    iter:        Current iteration of the for loop
    redshift:    The redshift calculating at
    range_a:     The number of data points in the current computing redshift interval
    range_b:     The number of data points in the previous redshift interval
    phi0:        The phi0 constant
    phi0_diff:   The change in phi0 over the redshift interval
    mc:          The mc constant
    mc_diff:     The change in mc over the redshift interval
    alpha:       The power-law relating how mass evolves over redshift
    alpha_diff:  The change in the power-law over the redshift interval
    phen:        The phenomena (ULX | SMBH-IMBH merger) wishing to compute rate of
    """
    
    if phen == 'ULX':
        eff = 0.6
        const_rate = (11/3)
        NIMBH = 10
    else:
        eff = 1
        const_rate = (4/3)
        NIMBH = 1

    phi0_val = phi0 + phi0_diff * (iter-range_b)/(range_a-range_b)
    mc_val = mc + mc_diff * (iter-range_b)/(range_a-range_b)
    alpha_val = alpha + alpha_diff * (iter-range_b)/(range_a-range_b)
    PS_val = quad(PS_function, 10**8, 10**14, args=(phi0_val, mc_val, alpha_val))[0]

    Nevents = (4/3)*np.pi*eff*NIMBH*const_rate*10**-12 * PS_val * GC_formation(redshift)
    Nevents_const = (4/3)*np.pi*eff*NIMBH*const_rate*10**-12 * PS_val

    return Nevents, Nevents_const

def cmove_dist(z):
    """
    The distance based on redshift.
    The constants are taken from Planck 2018.
    Value is in Mpc.

    Inputs:
    z:      The redshift to integrate over
    """

    H0 = 67
    c = constants.c.value_in(units.kms)
    omegaL = 0.673
    omegaM = 0.315
    return c/H0 * (np.sqrt(omegaL+omegaM*(1+z)**3))**-1

def GC_formation(zGC):
    """
    The GC formation rate over cosmic time.
    Taken from Gratton et al. 1997, van den Berg et al. 2013 and El-Badry et al. 2019
    Constant parameters are taken from Antonini et al. 2015.
    Value is dimensionless.

    Inputs:
    zGC:    The redshift integrating over
    """

    z = 3.2
    sigmaGC = 1.5
    return np.exp(-(z-zGC)/sigmaGC)

def phenomena_event(phen):
    """
    Function which manipulates the data and plots it depending on the
    phenomena (phen) asked.
    """

    plot_ini = plotter_setup()

    # Data from Furlong et al. (2015)
    redshift_range = [0, 0.5, 1, 2, 3]
    phi0_range = np.asarray([8.4, 8.4, 7.4, 4.5, 2.2]) 
    phi0_range *= 10**-4
    mc_range = [10**11.14, 10**11.11, 10**11.06, 10**10.91, 10**10.78]
    alpha_range = [-1.43, -1.45, -1.48, -1.57, -1.66]

    ps_int = []
    GC_form_int = []
    cmove_dist_int = []
    for i in range(len(redshift_range)-1):
        ps_int.append(quad(PS_function, 10**8, 10**14, args=(phi0_range[i], mc_range[i], alpha_range[i]))[0])
        GC_form_int.append(quad(GC_formation, redshift_range[i], redshift_range[i+1])[0])
        cmove_dist_int.append(quad(cmove_dist, redshift_range[i], redshift_range[i+1])[0])

    phi0_diffs = []
    mc_diffs = []
    alpha_diffs = []
    for i in range(len(phi0_range)-1):
        phi0_diffs.append(phi0_range[i+1] - phi0_range[i])
        mc_diffs.append(mc_range[i+1] - mc_range[i])
        alpha_diffs.append(alpha_range[i+1] - alpha_range[i])
    phi0_diffs = np.asarray(phi0_diffs)
    mc_diffs = np.asarray(mc_diffs)
    alpha_diffs = np.asarray(alpha_diffs)

    redshift_vals = np.linspace(0, 3, 1000)
    range_1 = len(redshift_vals[(redshift_vals<=0.5)])
    range_2 = len(redshift_vals[(redshift_vals<=1.0)])
    range_3 = len(redshift_vals[(redshift_vals<=2.0)])
    range_4 = len(redshift_vals[(redshift_vals<=3.0)])

    rate = []
    rate_const = []

    iter = 0
    for z_ in redshift_vals:
        iter += 1
        if z_ <= 0.5:
            Nev, Nevc = calc_loop(iter, z_, range_1, 0, phi0_range[0], phi0_diffs[0], 
                                mc_range[0], mc_diffs[0], alpha_range[0], 
                                alpha_diffs[0], str(phen))
            rate.append(Nev)
            rate_const.append(Nevc)
            
        elif z_ > 0.5 and z_ <= 1.0:
            Nev, Nevc = calc_loop(iter, z_, range_2, range_1, phi0_range[1], phi0_diffs[1], 
                                  mc_range[1], mc_diffs[1], alpha_range[1], 
                                  alpha_diffs[1], str(phen))
            rate.append(Nev)
            rate_const.append(Nevc)
            
        elif z_ > 1.0 and z_ <= 2.0:
            Nev, Nevc = calc_loop(iter, z_, range_3, range_2, phi0_range[2], phi0_diffs[2], 
                                  mc_range[2], mc_diffs[2], alpha_range[2], 
                                  alpha_diffs[2], str(phen))
            rate.append(Nev)
            rate_const.append(Nevc)

        elif z_ > 2.0 and z_ <= 3.0:
            Nev, Nevc = calc_loop(iter, z_, range_4, range_3, phi0_range[3], phi0_diffs[3], 
                                  mc_range[3], mc_diffs[3], alpha_range[3], 
                                  alpha_diffs[3], str(phen))
            rate.append(Nev)
            rate_const.append(Nevc)
            
    rate = np.asarray(rate)
    rate_const = np.asarray(rate_const)

    rate *= 10**3
    rate_const *= 10**3

    cum_GCz = []
    cum_GCc = []
    value_temp = 0
    value_const_temp = 0

    if phen == 'ULX':
        eff = 0.6
        const_rate = (11/3)
        NIMBH = 10
    else:
        eff = 1
        const_rate = (4/3)
        NIMBH = 1

    for i in range(len(ps_int)):
        value_temp += (4/3)*np.pi*eff*NIMBH*const_rate*10**-12*(ps_int[i] * GC_form_int[i] * cmove_dist_int[i]**3) * 10**-6
        value_const_temp += (4/3)*np.pi*eff*NIMBH*const_rate*10**-12*(ps_int[i] * cmove_dist_int[i]**3) * 10**-6
        cum_GCz.append(value_temp)
        cum_GCc.append(value_const_temp)
    zvals = [0, 0.5, 1, 2, 3]
    cum_GCz.insert(0,0)
    cum_GCc.insert(0,0)

    fz = interp1d(zvals, cum_GCz, kind='quadratic')
    fc = interp1d(zvals, cum_GCc, kind='quadratic')
    xnew = np.linspace(0,3,1000)

    fz = fz(xnew) + rate[0]
    fc = fc(xnew) + rate_const[0]
    for i in range(len(fz)):
        fz[i] -= rate[0] * (i/len(fz))
        fc[i] -= rate_const[0] * (i/len(fc))

    if phen != 'ULX':
        zvals = [0, 0.5, 1, 1.38, 2, 2.591, 3]
        cum_GCz.insert(3, 50.8)
        cum_GCc.insert(3, 414)
        cum_GCz.insert(5, 103)
        cum_GCc.insert(5, 551)

    fig = plt.figure(figsize=(12.5, 5))
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)
    ax1.set_ylabel(r'Events [Gpc$^{-3}$ yr$^{-1}$]')
    ax2.set_ylabel(r'Events [yr$^{-1}$]')
    for ax_ in [ax1, ax2]:
        ax_.set_xlim(0,3)
        ax_.set_xlabel(r'Redshift')
        ax_.xaxis.set_ticks_position('bottom')
        ax_.yaxis.set_ticks_position('left')
        plot_ini.tickers(ax_, 'plot')
    ax1.set_yscale('log')
    ax2.set_yscale('log')
    ax1.plot(redshift_vals, rate_const, color = 'red')
    ax1.plot(redshift_vals, rate, color = 'black')
    ax2.plot(xnew, fc, color = 'red', label = r'$N_{GC,const}$')
    ax2.plot(xnew, fz, color = 'black', label = r'$N_{GC}(z)$')
    ax1.plot(redshift_vals, 0.1*rate, color = 'black', linestyle = 'dashdot')
    ax1.plot(redshift_vals, 0.01*rate, color = 'black', linestyle = 'dashed')
    ax2.plot(xnew, 0.1*fz, color = 'black', linestyle = 'dashdot', label = r'$f_{IMBH} = 0.1$')
    ax2.plot(xnew, 0.01*fz, color = 'black', linestyle = 'dashed', label = r'$M_{IMBH} \leq 0.1M_{GC}$')
    if phen == 'ULX':
        ax1.set_ylim(5e-1, 800) 
        ax2.set_ylim(5e-1, 12000)
    else:   
        ax1.set_ylim(3e-2, 60) 
        ax2.set_ylim(3e-2, 1000)
    ax2.legend()
    plt.savefig('figures/forecast/'+str(phen)+'_rate.pdf', dpi=300, bbox_inches='tight')

def PS_function(zmass, phi0, mc, alpha):
    """
    The Press-Schecter stellar mass function (Press, Schechter 1974).
    Value is in (Mpc)^-3
    
    Inputs:
    zmass:  The redshift-dependent mass integrating over
    phi0:   The redshift-dependent normalisation factor
    mc:     The redshift-dependent characteristic mass. Describes the knee
    alpha:  The redshift-dependent power-law mass slope
    """

    ps_func = phi0 * (zmass/mc)**alpha*np.exp(-zmass/mc)
    return ps_func

phen = ['ULX', 'SMBH-IMBH_merger']
for p_ in phen:
    phenomena_event(p_)