import time as cpu_time
from loss_cone_plotters import *
from tGW_plotters import *
from sustainable_bintert_plotters import *
from steady_plotters import *
from ejection_stat_plotters import *
from spatial_plotters import *

start_time = cpu_time.time()

print('...loss_cone_plotters...')
cst = loss_cone()
cst.lcone_evolution_plotter()
cst.lcone_fininit_plotter()

print('...tGW_plotters...')
cst = bin_tert_systems()
cst.semi_ecc_SMBH_gw_plotter()
cst.semi_ecc_bin_gw_plotter()

print('...sustainable_bintert_plotters...')
cst = sustainable_sys()
cst.system_formation()

print('...steady_plotter...')
cst = stability_plotters()
cst.massdep_plotter()

print('...ejection_stat_plotters...')
cst = KE_PE_plotters()
cst.KEPE_plotter()
cst = vejection()
cst.vejec_histogram()
cst.vejec_mean_plotter()
cst = event_tracker()

print('...spatial_plotters...')
ejected_evolution('Hermite')
spatial_plotter('Hermite')
energy_plotter('Hermite')

end_time = cpu_time.time()
print('Plotting time:', start_time - end_time)