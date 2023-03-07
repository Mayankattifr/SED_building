from astropy.table import Table
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import scipy as ss

radius=5



RA=63.4779


Dec=28.1922

Name='NAME IRAS 04108+2803B'

target=str(RA)+'%20'+str(Dec)

#target='IRAS%04154+2823'
sed=Table.read(f"https://vizier.cds.unistra.fr/viz-bin/sed?-c={target}&-c.rs={radius}&-c.u=arcsec")
# print(sed)
dd = sed.to_pandas()

m=3e14/(sed["sed_freq"]*1e9) >1
dd=dd[m]

m1=dd['sed_flux']/dd['sed_eflux']>3
dd=dd[m1]
dd=dd.sort_values(by=['sed_freq'])
dd=dd.drop_duplicates(subset=['sed_filter'])
dis=140
# 
#%%
from scipy import interpolate
fig = plt.figure()
ax = fig.add_subplot()
ax.set_xlabel('Wave (microns)')
ax.set_ylabel('Flux (Jy)')
ax.scatter(3e14/(dd["sed_freq"]*1e9), dd["sed_flux"])
ax.set_yscale('log')
ax.set_xscale('log')
plt.show()
x=np.linspace(min(dd["sed_freq"]*1e9), max(dd["sed_freq"]*1e9),20000)
f= interpolate.interp1d(dd["sed_freq"]*1e9, dd["sed_flux"])
plt.plot(3e14/x,f(x))
dd.to_csv(str(Name)+".csv",index=False)

Lbol=4*3.14*((3.086e+18*dis)**2)*1e-23*np.trapz(dd["sed_flux"],x=dd["sed_freq"]*1e9)/(3.8e33)

Tbol=1.25e-11*np.trapz(dd["sed_flux"]*dd["sed_freq"]*1e9,x=dd["sed_freq"]*1e9)/np.trapz(dd["sed_flux"],x=dd["sed_freq"]*1e9)

print(Lbol,Tbol)

