#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 27 16:07:41 2020

@author: saus
"""

import pylab
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy

data = numpy.loadtxt("data.dat")
#zgrid функции 
def makeData():
    x = numpy.arange(-10, 10, 0.05)
    y = numpy.arange(-10, 10, 0.05)
    xgrid, ygrid = numpy.meshgrid(x, y)
    zgrid = (xgrid*numpy.exp((7+xgrid)/2.)+pow(ygrid,2)*numpy.exp((7+xgrid)/2.)-8.*ygrid*numpy.exp((7+xgrid)/2.)+23.*numpy.exp((7+xgrid)/2.))
    return xgrid, ygrid, zgrid


if __name__ == '__main__':
    x, y, z = makeData()
    fig, ax = plt.subplots()
    fig2 = pylab.figure()
    axes = Axes3D(fig)
    axes.plot_surface(x, y, z)
    cs = pylab.contour(x, y, z, levels = 20,linestyles='dashed', linewidths=1)
    ax.clabel(cs,inline=1, fontsize=10)
    plt.plot(data[:,0],data[:,1],marker = '*')
    cs = pylab.contourf(x, y, z)
    
    fig2.set_figwidth(16)    
    fig2.set_figheight(8)
    pylab.show()