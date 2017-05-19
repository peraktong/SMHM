import numpy as np
from astropy.io import fits
from astropy.table import Table


path = "/Users/caojunzhi/Desktop/apVisit-r6-4812-55725-003.fits"

ap_star = ""

star = fits.open(path)

VHELIO = star[0].header["VHELIO"]

print(VHELIO)

