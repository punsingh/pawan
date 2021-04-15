"""@plotWake docstring
Documentation for plotWake package
@author Puneet Singh
@date 04/15/2021
"""
import struct
import numpy as np
import readWake as rw
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

class plotWake:
    """ Documentation for plotWake
    """
    
    def __init__(self,fileName):
        ''" Constructor """
        data = rw.readWake(fileName)
        self.nParticles = data.nParticles
        self.position = data.position
        self.x = self.position[:,0]
        self.y = self.position[:,1]
        self.z = self.position[:,2]
        self.vorticity = data.vorticity
        self.p = self.vorticity[:,0]
        self.q = self.vorticity[:,1]
        self.r = self.vorticity[:,2]
        self.radius = data.radius
        self.volume = data.volume
        self.birthstrength = data.birthstrength
        self.strength = np.sqrt(self.vorticity[:,0]**2+self.vorticity[:,1]**2+self.vorticity[:,2]**2)
        self.maxstrength = max(self.strength)
        self.strengthcolors = self.strength/np.sum(self.strength)
    
    def plot3DQuiver(self):
        """ plotDQuiver Plots the vortex particles as colored arrows 
        """
        norm = matplotlib.colors.Normalize()
        norm.autoscale(self.strengthcolors)
        cm = matplotlib.cm.viridis
        sm = matplotlib.cm.ScalarMappable(cmap=cm, norm=norm)
        sm.set_array([])
        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')
        ax.scatter(self.x,self.y,self.z,color=cm(norm(self.strengthcolors)))
        ax.quiver(self.x,self.y,self.z,self.p,self.q,self.r,length=0.1,color=cm(norm(self.strengthcolors)),normalize=True)
        plt.show()


    def printData(self):
        """ printData prints wake properties
        """
        print(self.nParticles)
        print(self.position)
        print(self.vorticity)
        print(self.radius)
        print(self.volume)
        print(self.birthstrength)
        print(self.strength)
        print(self.maxstrength)
        print(self.strengthcolors)


if __name__ == "__main__":
    pw = plotWake("data/temp.wake")
    pw.plot3DQuiver()
    # pw.printData()

