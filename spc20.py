
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
#v_sun_phi = 251.3*km_s # speed of the sun in the phi direction

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

# Initializing Potentials

bulge_pot = HernquistPotential(amp=2*mass_bulge, a=c)
disk_pot = MiyamotoNagaiPotential(amp=G*mass_disk, a=a, b=b)

Stan_halo_pot = LogarithmicHaloPotential(amp=vc_halo_Stan**2, core=r_h)
Standard = bulge_pot + disk_pot + Stan_halo_pot

Mass_halo_pot = LogarithmicHaloPotential(amp=vc_halo_Mass**2, core=r_h)
Massive = bulge_pot + disk_pot + Mass_halo_pot

Smal_halo_pot = LogarithmicHaloPotential(amp=vc_halo_Smal**2, core=r_h)
Small = bulge_pot + disk_pot + Smal_halo_pot

no_Halo = bulge_pot + disk_pot

###########################

# Finding circular velocities due to potentials


from galpy.potential import vcirc

ov_Stan = vcirc(Standard, r_sun)
ov_Mass = vcirc(Massive, r_sun)
ov_Smal = vcirc(Small, r_sun)

#print("Stan vel", ov_Stan, "upper error", vcirc(Standard, r_sun+dr_sun)-ov_Stan, "lower error", ov_Stan-vcirc(Standard, r_sun-dr_sun))
#print("Mass vel", ov_Mass, "upper error", vcirc(Massive, r_sun+dr_sun)-ov_Mass, "lower error", ov_Mass-vcirc(Massive, r_sun-dr_sun))
#print("Smal vel", ov_Smal, "upper error", vcirc(Small, r_sun+dr_sun)-ov_Smal, "lower error", ov_Smal-vcirc(Small, r_sun-dr_sun))
#print("bulge", mass_bulge.to(u.solMass), "disk", mass_disk.to(u.solMass))
#print("Stan vo", vc_halo_Stan.to(km_s), "Mass vo", vc_halo_Mass.to(km_s), "Smal vo", vc_halo_Smal.to(km_s))

ov_nohalo = vcirc(no_Halo, r_sun)
print("nohalo vel", ov_nohalo, "upper error", vcirc(no_Halo, r_sun+dr_sun)-ov_nohalo, "lower error", ov_nohalo-vcirc(no_Halo, r_sun-dr_sun))

###########################

# Plotting Rotation Curves to compare with original paper

from galpy.potential import plotRotcurve


'''
plotRotcurve(no_Halo)
plt.title = "Rotation Curve of Bulge and Disk Potentials"
plt.xlabel = r'$r\,(\mathrm{kpc})$'
plt.ylabel = r"$v_c(R)\,(\mathrm{km\,s}^{-1})$"
plt.xlim(0, 25)
plt.ylim(100, 280)
plt.grid()
plt.show()
'''
'''
plotRotcurve(bulge_pot)
plt.xlabel = r'$r\,(\mathrm{kpc})$'
plt.ylabel = r"$v_c(R)\,(\mathrm{km\,s}^{-1})$"
plt.xlim(0, 25)
plt.ylim(100, 280)
plt.grid()
plt.show()
'''
'''
plotRotcurve(disk_pot)
plt.xlabel = r'$r\,(\mathrm{kpc})$'
plt.ylabel = r"$v_c(R)\,(\mathrm{km\,s}^{-1})$"
plt.xlim(0, 25)
plt.ylim(100, 280)
plt.grid()
plt.show()
'''

'''
Stan_plot = plotRotcurve(Standard, overplot=True)
Mass_plot = plotRotcurve(Massive, overplot=True)
Smal_plot = plotRotcurve(Small, overplot=True)

plt.xlim(0, 25)
plt.ylim(100, 280)
plt.xlabel = r'$r\,(\mathrm{kpc})$'
plt.ylabel = r"$v_c(R)\,(\mathrm{km\,s}^{-1})$"
plt.grid()
plt.show()
'''

###########################

# Initializing orbit of Palomar
'''
from galpy.orbit import Orbit

mas_yr = 1e-3 * u.arcsec / u.yr

RA = 229.018*u.deg
Dec = -0.124*u.deg
Dist = (20.9 + (23.2-20.9)/2)*u.kpc
mu_ra = -2.296*mas_yr
mu_dec = -2.257*mas_yr
v_Pal = -58.7*km_s

init_cond = [RA, Dec, Dist, mu_ra, mu_dec, v_Pal]

V_sun_r = -11.1#*km_s
V_sun_phi = 251.3#*km_s
V_sun_z = 7.3#*km_s

#V_sun = [V_sun_r, V_sun_phi, V_sun_z]

Palomar = Orbit(init_cond, radec=True)

ts= np.linspace(0.,300.,1001)

Palomar.integrate(ts, Standard)
Palomar.plot3d(alpha=0.4)
plt.xlim(-100.,100.)
plt.ylim(-100.,100)
plt.show()
'''

###########################

# Initializing Stream

'''
from galpy.df import streamdf
from galpy.actionAngle import actionAngleIsochroneApprox

rad_disp = 0.4*km_s

aAI = actionAngleIsochroneApprox(b=20.9+(23.2-20.9)/2, pot=Standard, ro=r_sun, vo=V_sun_phi)
sdf = streamdf(sigv=rad_disp, projenitor=Palomar, pot=Standard, aA=aAI, tdisrupt=5e9*u.yr)

sdf.plotTrack(d1='ll',d2='dist',interp=True,color='k',spread=2,overplot=False,lw=2.)
sdf.plotTrack(d1='ll',d2='dist',interp=False,color='k',spread=0,overplot=True,ls='none',marker='o')
sdf.plotProgenitor(d1='ll',d2='dist',color='r',overplot=True,ls='--')
xlim(155.,255.); ylim(7.5,14.8)
'''