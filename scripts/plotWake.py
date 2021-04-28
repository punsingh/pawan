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
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D

class plotWake:
    """ Documentation for plotWake
    """
    
    def __init__(self,fileName):
        ''" Constructor """
        data = rw.readWake(fileName)
        self.nTimesteps = data.nTimesteps
        self.nParticles = data.nParticles
        self.position = data.position
        self.vorticity = data.vorticity
        self.radius = data.radius
        self.volume = data.volume
        self.birthstrength = data.birthstrength
        
        self.position_arr = np.asarray(self.position)        
        multiplier = 1.25
        X = multiplier
        xmax = X*max(np.amax(self.position_arr[:,:,0]),np.amax(self.position_arr[:,:,1]))
        xmin = X*min(np.amin(self.position_arr[:,:,0]),np.amin(self.position_arr[:,:,1]))
        ymax = X*xmax
        ymin = X*xmin
        zmax = X*np.amax(self.position_arr[:,:,2])
        zmin = X*np.amin(self.position_arr[:,:,2])
        self.plot_limits = {'x':[xmax,xmin],'y':[ymax,ymin],'z':[zmax,zmin]}
        
    def getTimeStep(self,n=0):
        """ 
        getTimeStep Loads particle position and vorticity at nth time step 
        
        Parameters
        ----------
        n: int
            Time step index
        """
        self.x = self.position[n][:,0]
        self.y = self.position[n][:,1]
        self.z = self.position[n][:,2]
        self.p = self.vorticity[n][:,0]
        self.q = self.vorticity[n][:,1]
        self.r = self.vorticity[n][:,2]
        self.strength = np.sqrt(self.p**2+self.q**2+self.r**2)
        self.maxstrength = max(self.strength)
        self.strengthcolors = self.strength/np.sum(self.strength)
    
    def plot3DQuiver(self,n=0):
        """ 
        plot3DQuiver Plots the vortex particles as colored arrows 
        
        Parameters
        ----------
        n: int
            Time step index
        """
        self.getTimeStep(n)
        y = self.setColors()
        fig = plt.figure()
        print('plotting static figure')
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(self.x,self.y,self.z,color=y)
        ax.quiver(self.x,self.y,self.z,self.p,self.q,self.r,length=0.1,color=y,normalize=True)
        ax.set_xlim([-2,2])
        ax.set_ylim([-2,2])
        ax.set_zlim([-1,3])
        plt.show()

    def plotAnimate(self):
        """ 
        plotAnimate Plots the vortex particles as colored spheres 
        """
        fig = plt.figure(figsize=(12,10))
        print('plotting animation')
        ax = fig.add_subplot(111, projection='3d')
        self.getTimeStep(0)
        self.scplot = ax.scatter(self.x,self.y,self.z)
        self.qvplot = ax.quiver(self.x,self.y,self.z,self.p,self.q,self.r,length=0.2)
        ax.set_xlim(self.plot_limits['x'])
        ax.set_ylim(self.plot_limits['y'])
        ax.set_zlim(self.plot_limits['z'])
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        
        def update(num):
            self.getTimeStep(n=num)
            y = self.setColors()
            self.qvplot.remove()
            self.scplot.remove()
            self.qvplot = ax.quiver(self.x,self.y,self.z,self.p,self.q,self.r,length=0.1,color=y,normalize=True)
            self.scplot = ax.scatter(self.x,self.y,self.z,color=y)
        ani = animation.FuncAnimation( fig, update, self.nTimesteps, interval = 200)
        plt.show()
        ani.save('wake.mp4',writer='ffmpeg',fps=2, bitrate = 1000)

    def setColors(self):
        norm = matplotlib.colors.Normalize()
        norm.autoscale(self.strengthcolors)
        cm = matplotlib.cm.viridis
        sm = matplotlib.cm.ScalarMappable(cmap=cm, norm=norm)
        sm.set_array([])
        y = cm(norm(self.strengthcolors))
        return y


    def printData(self):
        """ printData prints wake properties
        """
        print(self.nParticles)
        print(self.position)
        print(self.vorticity)
        print(self.radius)
        print(self.volume)
        print(self.birthstrength)

if __name__ == "__main__":
    pw = plotWake("../data/temp.wake")
    # pw.plotAnimate()
    pw.plotAnimate()
    # pw.plot3DQuiver(0)
    # pw.plot3DQuiver(1)
    # pw.plot3DQuiver(2)
    # pw.plot3DQuiver(3)
    # pw.plot3DQuiver(4)
    # pw.plot3DQuiver(5)
    # pw.printData()

