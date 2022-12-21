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

print('... ejection_Stat_plotters ...')
cst = event_tracker()
cst = ejection_stats()
cst.new_data_extractor()
cst.combine_data()
cst.energy_plotters()
cst.vejec_plotters()

print('...spatial_plotters...')
ejected_evolution()
global_properties()
energy_plotter('Hermite')
spatial_plotter('Hermite')
#chaos_deviate()

print('...loss_cone_plotters...')
cst = loss_cone()
cst.lcone_plotter()
cst.lcone_timescale()
cst.lcone_fininit_plotter()

print('...sustainable_bintert_plotters...')
cst = sustainable_sys()
cst.system_formation_plotter()

print('...tGW_plotters...')
cst = gw_calcs()
#cst.new_data_extractor()
cst.orbital_hist_plotter()
cst.Nenc_tgw_plotter()
cst.strain_freq_plotter()
cst.transient_events()

print('Plotting time [mins]:', (end_time - start_time)/60)