from tGW_plotters import *
from sustainable_bintert_plotters import *
from steady_plotters import *
from ejection_stat_plotters import *
from spatial_plotters import *

print('...GW plotters...')
cst = bin_tert_systems()
cst.semi_ecc_SMBH_gw_plotter()
cst.semi_ecc_bin_gw_plotter()

print('...System plotters...')
cst = sustainable_sys()
cst.system_formation()

print('...Energy plotters...')
cst = KE_PE_plotters()
cst.KEPE_plotter()

print('...Ejection plotters...')
cst = vejection()
cst.vejec_histogram()
cst.vejec_mean_plotter()

print('...Merging ratio plotter...')
cst = event_tracker()

print('...Stability time plotters...')
cst = stability_plotters()
cst.massdep_plotter()

print('...Indiv. simulation plotters...')
spatial_plotter('Hermite')
energy_plotter('Hermite')
ejected_evolution('Hermite')