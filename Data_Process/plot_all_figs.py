import time as cpu_time
from loss_cone_plotters import *
from tGW_plotters import *
from sustainable_bintert_plotters import *
from steady_plotters import *
from ejection_stat_plotters import *
from spatial_plotters import *

start_time = cpu_time.time()


print('...steady_plotter...')
cst = stability_plotters()
cst.overall_steady_plotter()
end_stab = cpu_time.time()
print('Plotting time [mins]:', (end_stab - start_time)/60)
"""
print('...spatial_plotters...')
#ejected_evolution()
global_properties()
energy_plotter('Hermite')
spatial_plotter('Hermite')
#chaos_deviate()
energy_plotter('GRX')
spatial_plotter('GRX')
end_spat = cpu_time.time()
print('Plotting time [mins]:', (end_spat - end_stab)/60)"""


print('...sustainable_bintert_plotters...')
cst = sustainable_sys()
cst.system_formation_plotter()
cst.GW_emissions()
end_sys = cpu_time.time()
#print('Plotting time [mins]:', (end_sys - end_spat)/60)

print('...loss_cone_plotters...')
cst = loss_cone()
cst.lcone_plotter()
cst.lcone_timescale()
cst.lcone_fininit_plotter()
end_LC = cpu_time.time()
#print('Plotting time [mins]:', (end_LC - end_sys)/60)

print('... ejection_Stat_plotters ...')
cst = event_tracker()
cst = ejection_stats()
cst.new_data_extractor()
cst.combine_data()
cst.energy_plotters()
cst.vejec_plotters()
end_ejec = cpu_time.time()
print('Plotting time [mins]:', (end_ejec - end_LC)/60)

print('...tGW_plotters...')
cst = gw_calcs()
#cst.new_data_extractor()
cst.orbital_hist_plotter()
cst.Nenc_tgw_plotter()
cst.strain_freq_plotter()
cst.transient_events()
end_tGW = cpu_time.time()

print('Plotting time [mins]:', (end_tGW - end_ejec)/60)

print('Plotting time [mins]:', (end_tGW - start_time)/60)