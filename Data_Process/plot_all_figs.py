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
end_time = cpu_time.time()
print('Plotting time [mins]:', (end_time - start_time)/60)
start_time2 = end_time
print('...loss_cone_plotters...')
cst = loss_cone()
cst.lcone_plotter()
cst.lcone_timescale()
cst.lcone_fininit_plotter()
end_time = cpu_time.time()
print('Plotting time [mins]:', (end_time - start_time)/60)
start_time2 = end_time

print('...sustainable_bintert_plotters...')
cst = sustainable_sys()
cst.system_formation_plotter()
end_time = cpu_time.time()
print('Plotting time [mins]:', (end_time - start_time2)/60)
start_time2 = end_time

print('...tGW_plotters...')
cst = coupled_systems()
#cst.new_data_extractor()
cst.combine_data()
cst.IMBH_tgw_plotter()
cst.SMBH_tgw_plotter()
cst.strain_freq_plotter()
cst.IMBH_tgw_plotter()
cst.transient_events()
end_time = cpu_time.time()
print('Plotting time [mins]:', (end_time - start_time2)/60)


##########################
#### NO UPDATE NEEDED ####
##########################


print('...spatial_plotters...')
#ejected_evolution()
global_properties()
#spatial_plotter('Hermite')
energy_plotter('Hermite')
#chaos_deviate()
end_time = cpu_time.time()
print('Plotting time [mins]:', (end_time - start_time2)/60)
start_time2 = end_time

print('Plotting time [mins]:', (end_time - start_time)/60)