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

print('...spatial_plotters...')
spatial_plotter('Hermite')
ejected_evolution()
energy_plotter('GRX')
nearest_neigh('Hermite')
#chaos_deviate()

print('...sustainable_bintert_plotters...')
cst = sustainable_sys()
cst.system_formation_plotter()


print('...loss_cone_plotters...')
cst = loss_cone()
cst.lcone_plotter()
cst.lcone_timescale()
cst.lcone_fininit_plotter()

print('...tGW_plotters...')
cst = coupled_systems()
cst.transient_events()
cst.SMBH_tgw_plotter()
cst.strain_freq_plotter()
cst.IMBH_tgw_plotter()

print('...ejection_stat_plotters...')
cst = event_tracker()
cst = vejection()
cst.vejec_plotter()
cst = KE_PE_plotters()
cst.KEPE_plotter()


end_time = cpu_time.time()
print('Plotting time [mins]:', (end_time - start_time)/60)