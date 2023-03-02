#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 13 15:28:30 2022

@author: Sumeet Kumar
"""

import pickle
import numpy as np
import tecplot as tp
import readWake as rw
from Tecplot import utils as tecutils
from Tecplot import plot
from tecplot.constant import *
import subprocess

print('running...')
directry = "../data"
filename = "temp.wake"
#filename = "wing_rectangularCOARSE_pitchingPawan.wake"
wake = rw.readWake(f"{directry}/{filename}")

# tecutils.gen_plt(directry,filename,wake) #not working for now
# generate *.dat and then save as plt
datfilepath = tecutils.gen_dat(directry,filename,wake) 
subprocess.run(["/usr/local/tecplot/360ex_2019r1/bin/preplot", datfilepath])

pltfilepath = datfilepath.replace(".dat",".plt")
# plot.scatter(pltfilepath)
tp.session.connect()
tp.new_layout()
dataset = tp.data.load_tecplot(pltfilepath)

frame = tp.active_frame()
frame.plot_type = PlotType.Cartesian3D
plot = frame.plot()
plot.contour(0).variable = dataset.variable('vor_x')
plot.show_scatter = True

plot.scatter.variable = dataset.variable('strength')

for z in dataset.zones():
    scatter = plot.fieldmap(z).scatter
    scatter.symbol_type = SymbolType.Geometry
    scatter.symbol().shape = GeomShape.Sphere
    scatter.fill_mode = FillMode.UseSpecificColor
    scatter.fill_color = plot.contour(0)
    scatter.color = plot.contour(0)
    scatter.size_by_variable = False

# frame.add_text('Size of dots indicate relative pressure', (20, 80))

# # ensure consistent output between interactive (connected) and batch
# plot.contour(0).levels.reset_to_nice()
    
