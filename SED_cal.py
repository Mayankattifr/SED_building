from astropy.table import Table
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import scipy as ss
import pandas as pd

import glob





def BLT(file):

    dd=pd.read_csv(file)
    m1=dd['sed_flux']/dd['sed_eflux']>3
    dd=dd[m1]
    dd=dd.sort_values(by=['sed_freq'])
    dd=dd.drop_duplicates(subset=['sed_filter'])
    dis=140
    
    fig = plt.figure()
    ax = fig.add_subplot()
    ax.set_xlabel('Wave (microns)')
    ax.set_ylabel('Flux (Jy)')
    ax.plot(3e14/(dd["sed_freq"]*1e9), dd["sed_flux"],'--o')
    ax.set_yscale('log')
    ax.set_xscale('log')
    plt.show()
    
    
    
    Lbol=4*3.14*((3.086e+18*dis)**2)*1e-23*np.trapz(dd["sed_flux"],x=dd["sed_freq"]*1e9)/(3.8e33)
    
    Tbol=1.25e-11*np.trapz(dd["sed_flux"]*dd["sed_freq"]*1e9,x=dd["sed_freq"]*1e9)/np.trapz(dd["sed_flux"],x=dd["sed_freq"]*1e9)
    
    print(Lbol,Tbol)
    
    return Lbol,Tbol



i=0

Lbol=[]
Tbol=[]
Star=[]

for name in glob.glob('*csv'):
    xx=BLT(name)
    Lbol.append(xx[0])
    Tbol.append(xx[1])
    Star.append(name)
    i=i+1
    
df=pd.DataFrame()

df['Star'] =Star
df['Lbol'] = Lbol
df['Tbol'] = Tbol


df.to_csv('Taurus_protostar.txt')
    