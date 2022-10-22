from astropy.coordinates import SkyCoord
import astropy.coordinates as coord
import numpy as np
from amuse.lab import units
from astropy import units as u

# Data is from https://ui.adsabs.harvard.edu/abs/1996AJ....112.1487H/abstract

def conv_coord(right_ascension, declination, dist):

    coordinates = SkyCoord(ra=right_ascension*u.degree, 
                           dec=declination*u.degree, 
                           distance=dist*u.kpc)
    coordinates_co = coordinates.transform_to(coord.Galactocentric)
    x_co = np.array(coordinates_co.cartesian.x)
    y_co = np.array(coordinates_co.cartesian.y)
    z_co = np.array(coordinates_co.cartesian.z)
    dist = np.sqrt(x_co**2+y_co**2+z_co**2)

    return dist

ra_dec = []
helio_r = []

with open('GC_datafile.txt', 'r') as inp_:
    for i, line in enumerate(inp_):
        if i == 0:
            continue
        line = line.split("|")
        c = SkyCoord(line[1], line[2], unit=(u.hourangle, u.deg))
        ra_dec.append(c)
        helio_r.append(float(line[3])) #in kpc
distances = []
for i in range(len(ra_dec)):
    distances.append(conv_coord(ra_dec[i].ra.radian, ra_dec[i].dec.radian, helio_r[i]))

print(distances)
print('Minimum distances to core of MW', min(distances))