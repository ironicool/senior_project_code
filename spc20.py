
import matplotlib.pyplot as plt
import numpy as np
import galpy

from astropy import units as u

from galpy.potential import HernquistPotential
from galpy.potential import MiyamotoNagaiPotential
from galpy.potential import LogarithmicHaloPotential
from galpy.potential import MWPotential2014

# Units
km_s = u.km / u.s

km_kpc_conversion = 3.086e16 # m/kpc, multiply to convert kpc to m
kg_msun_conversion = 1.989e30 # kg/Msun, multipy to convert Msun to kg

# Universal parameters

G = (6.67e-11 * (1 / 1000)**3 * kg_msun_conversion) * (u.km)**3 / (1*u.solMass * (1*u.s)**2)

r_sun = 8.31*u.kpc # radial distance of the sun for the galactic core
dr_sun = 0.08*u.kpc # r_sun error
v_sun_phi = 251.3*km_s # speed of the sun in the phi direction

a = 6.5 * u.kpc # shape parameter for disk
b = 0.26 * u.kpc # shape parameter for disk
c = 0.7 * u.kpc # shape parameter for bulge
r_h = 12.0 * u.kpc # shape parameter for halo, == core in galpy

vi0_bulge = 121.9*km_s
vi0_disk = 154.9*km_s

mass_interior = r_sun * vi0_bulge**2 / G # amount of bulge mass within the Sun's radius
mass_bulge = mass_interior * (r_sun + c)**2 / r_sun**2 # total mass of the bulge potential
mass_disk = vi0_disk**2 * (r_sun**2 + (a + b)**2)**(3/2) / (r_sun**2 * G) # total mass of the disk potential

# Standard Model parameters

# Bulge, Hernquist Potential
f_bulge_Stan = 0.307
# Disc, MiyamotoNagai Potential
f_disk_Stan = 0.496
# Halo, Logarithmic Halo Potential
f_halo_Stan = 0.197
vi0_halo_Stan = 97.65*km_s

vc_halo_Stan = vi0_halo_Stan * np.sqrt(r_sun**2 + r_h**2) / r_sun # amp of Log potential, circular speed far away from r=0

'''
print(mass_Standard.to(u.solMass), mass_test.to(u.solMass))
'''

# Massive Model parameters

# Bulge, Hernquist Potential
f_bulge_Mass = 0.270
# Disc, MiyamotoNagai Potential
f_disk_Mass = 0.4345
# Halo, Logarithmic Halo Potential
f_halo_Mass = 0.2955
vi0_halo_Mass = 127.75*km_s

vc_halo_Mass = vi0_halo_Mass * np.sqrt(r_sun**2 + r_h**2) / r_sun # amp of Log potential, circular speed far away from r=0


# Small Model parameters

# Bulge, Hernquist Potential
f_bulge_Smal = 0.337
# Disc, MiyamotoNagai Potential
f_disk_Smal = 0.555
# Halo, Logarithmic Halo Potential
f_halo_Smal = 0.119
vi0_halo_Smal = 72.44*km_s

vc_halo_Smal = vi0_halo_Smal * np.sqrt(r_sun**2 + r_h**2) / r_sun # amp of Log potential, circular speed far away from r=0

###########################

from galpy.potential import plotRotcurve


bulge_pot = HernquistPotential(amp=2*mass_bulge, a=c)
disk_pot = MiyamotoNagaiPotential(amp=G*mass_disk, a=a, b=b)

Stan_halo_pot = LogarithmicHaloPotential(amp=vc_halo_Stan**2, core=r_h)
Standard = bulge_pot + disk_pot + Stan_halo_pot

Mass_halo_pot = LogarithmicHaloPotential(amp=vc_halo_Mass**2, core=r_h)
Massive = bulge_pot + disk_pot + Mass_halo_pot

Smal_halo_pot = LogarithmicHaloPotential(amp=vc_halo_Smal**2, core=r_h)
Small = bulge_pot + disk_pot + Smal_halo_pot

plotRotcurve(Standard, overplot=True)
plotRotcurve(Massive, overplot=True)
plotRotcurve(Small, overplot=True)

plt.xlim(0, 25)
plt.ylim(100, 280)
plt.grid()
plt.show()

'''
# Testing
hp = HernquistPotential(amp=2*mass_Standard, a=c)
mp = MiyamotoNagaiPotential(amp=G*mass_test, a=a, b=b)
lp = LogarithmicHaloPotential(amp=vc_halo_Stan**2, core=r_h)

#plotRotcurve(hp)
#plotRotcurve(mp)
#plotRotcurve(lp)

plotRotcurve(hp + mp + lp)

plt.xlim(0, 25)
plt.ylim(100, 280)
plt.grid()
plt.show()
'''