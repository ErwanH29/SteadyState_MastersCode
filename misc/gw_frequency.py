import numpy as np
from amuse.lab import *

def freqISCO(mass1,mass2):
    """
    Function based on Maggiore (2008) to calculate the GW frequency at ISCO
    """
    return (6**1.5*2*np.pi)**-1 * (constants.c)**3/(constants.G*((mass1+mass2)*(1 | units.MSun)))

fISCO_SMBH = freqISCO(4*10**6, 10**3)
print('IMBH-SMBH fISCO:', fISCO_SMBH.in_(units.Hz))
LISAz = fISCO_SMBH/(10**-4 | units.Hz) - 1
print('LISA will observe till z < :', LISAz)

fISCO_IMBH = freqISCO(10**3, 10**3)
print('IMBH-SMBH fISCO:', fISCO_IMBH.in_(units.Hz))
LISAz = fISCO_IMBH/(10**-4 | units.Hz) - 1
print('LISA will observe till z < :', LISAz)