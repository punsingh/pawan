#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 28 10:50:07 2019

@author: SumeetKumar
"""
import matplotlib.pyplot as plt
#plt.ioff()
from scipy.interpolate import interp1d
import numpy as np
import math
import os
import tikzplotlib

def plot_simple(Fig,Ax,Title,Data,Limits):
    # f = interp1d(Xdata, Ydata, kind='cubic')
    # Xdata_smooth = np.linspace(min(Xdata), max(Xdata), 3*len(Xdata))
    # Ydata_smooth = f(Xdata_smooth)
    # Ax.plot(Xdata_smooth, Ydata_smooth, '-', label=Legend)
    Xdata = Data['x']
    Ydata = Data['y']
    
    Ax.plot(Xdata, Ydata, '-')
    Ax.plot(Xdata,Ydata,'blue',label='Dymore')
    if Data['expt']:
        Ydata_expt = Data['expt']
        Ax.plot(Xdata,Ydata_expt,'red',label='Expt')
    Ax.set_title(Title, fontsize=9,fontweight='bold',y=1.0) 
    Ax.legend()
    # Ax.set_xlabel(Xlabel,fontweight='bold')
    # Ax.set_ylabel(Ylabel,fontweight='bold')
    # Ax.set_xlim(0.0, 360.0)    #This needs to be changed if output of a non-azimuthal quantity is being plotted
    if Limits: Ax.set_ylim(Limits[0], Limits[1])   
    Ax.ticklabel_format(useOffset=False, style='plain')
    #plt.show()

def plotvalidation(Ax,Titel,data,limits):
    # f = interp1d(Xdata, Ydata, kind='cubic')
    # Xdata_smooth = np.linspace(min(Xdata), max(Xdata), 3*len(Xdata))
    # Ydata_smooth = f(Xdata_smooth)
    # Ax.plot(Xdata_smooth, Ydata_smooth, '-', label=Legend)
    linestyles=['-','--','-.',':','','--','-.',':']
    for ((label,arr),ls) in zip(data['data'].items(),linestyles):
        Ax.plot(arr[:,0],arr[:,1],label=label,linestyle=ls)
        
    Ax.legend()
    Ax.set_xlabel(data['xlabel'])
    Ax.set_ylabel(data['ylabel'])#,fontweight='bold')
    Ax.set_xlim(0.0, 5.0)    #This needs to be changed if output of a non-azimuthal quantity is being plotted
    Ax.set_ylim(limits[0], limits[1]) 
    Ax.grid(True)
    # if limits: Ax.set_ylim(limits[0], limits[1])   
    # Ax.ticklabel_format(useOffset=False, style='plain')
    #plt.show()


def plot_and_savesubplot(Figname,Data,Limits,Plot_dir,conf='AIAA_scitech',savetikz=False):

    Axislabel_size=17
    Ticklabel_size=14
    Legend_size=15
    Fig = plt.figure()
    Ax = Fig.add_subplot(1,1,1)
    # for Ydata,label in zip(Data['y'],Data['labels']):
    #     if len(Ydata)!=0: #works for both np arrays and lists
    #         Ax.plot(Xdata,Ydata,'blue',label='Dymore')
    #     else:
    #         print(f'no data for {label}')
    # if Data['expt']:
    #     Ax.plot(Xdata,Data['expt'],'red',label='Expt')
    # Ax.legend(frameon=False, prop={'size': Legend_size})
    # Ax.set_ylabel(Data['ylabel'],fontsize=Axislabel_size)
    # Ax.set_xlabel(Data['xlabel'], fontsize=Axislabel_size)
    # Ax.set_xlim(0.0, 5.0)    #This needs to be changed if output of a non-azimuthal quantity is being plotted
    # Ax.set_ylim(Limits[0], Limits[1]) 
    plotvalidation(Ax,'',Data,Limits)
    
    # Ax.set_xticks(np.linspace(min(Xdata),max(Xdata),7,endpoint=True))
    # Ax.tick_params(axis='both', labelsize=Ticklabel_size)
    # Ax.ticklabel_format(useOffset=False, style='plain')
    try:
        os.makedirs(f"{Plot_dir}/Savedfigs")
    except FileExistsError:
        pass    
    Fig.savefig(f"{Plot_dir}/Savedfigs/{Figname}.png", bbox_inches='tight')
    print('saved: ', Figname)
    if savetikz:
        tikzplotlib.save(f"{Plot_dir}/Savedfigs/{Figname}.pgf")    
    plt.close(Fig)
 
































































































