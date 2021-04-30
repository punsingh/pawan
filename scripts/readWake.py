"""@readWake docstring
Documentation for readWakeFile package
@author Puneet Singh
@date 04/14/2021
"""
import struct
import numpy as np

class readWake:
    """ Documentation for readWake
    Wake file includes:
        double  time
        size_t  _numParticles
        matrix  _position
        matrix  _vorticity
        vector  _radius
        vector  _volume
        vector  _birthstrength
    """
    def __init__(self,fileName):
        ''" Constructor """
        with open(fileName, mode='rb') as file:
            fileContent = file.read()
        self.fileSize = len(fileContent)
        self.time = []
        self.nParticles = []
        self.position = []
        self.vorticity = []
        self.radius = []
        # self.volume = []
        # self.birthstrength = []
        while len(fileContent)>0:
            # [t,n,p,q,r,v,b,fileContent] = self.getAllWakeState(fileContent)
            [t,n,p,q,r,fileContent] = self.getAllWakeState(fileContent)
            self.time.append(t)
            self.nParticles.append(n)
            self.position.append(p)
            self.vorticity.append(q)
            self.radius.append(r)
            # self.volume.append(v)
            # self.birthstrength.append(b)
        self.nTimesteps = len(self.time)
    
    def getInteger(self,fileContent):
        """ getInteger returns an integer 
        """
        ans = struct.unpack("i", fileContent[:4])[0]
        fileContent = fileContent[4:]
        return [ans, fileContent]
    
    def getDouble(self,fileContent):
        """ getDouble returns a double 
        """
        ans = struct.unpack("d", fileContent[:8])[0]
        fileContent = fileContent[8:]
        return [ans, fileContent]
    
    def getVector(self,fileContent,n):
        """ getVector returns a double vector of length n
        """
        ans = np.array(struct.unpack("d"*n,fileContent[:8*n]));
        fileContent = fileContent[8*n:]
        return [ans, fileContent]

    def getMatrix(self,fileContent,m,n=3):
        """ getMatrix returns a double matrix with m rows and n columns
        """
        ans = np.array(struct.unpack("d"*m*n,fileContent[:8*m*n])).reshape((m,n))
        fileContent = fileContent[8*m*n:]
        return [ans, fileContent]
    
    def getAllWakeState(self, fileContent):
        """ getAllWakeState returns wake data from file
        """
        [time, fileContent] = self.getDouble(fileContent)
        [nWake, fileContent] = self.getInteger(fileContent)
        fileContent = fileContent[4:]   # Skipping 4 bit blank space
        nParticles = 0 
        position = []
        vorticity = []
        radius = []
        # volume = []
        # birthstrength = []
        for n in range(nWake):
            # [npar,pos,vort,rad,vol,bs,fileContent] = self.getWakeState(fileContent)
            [npar,pos,vort,rad,fileContent] = self.getWakeState(fileContent)
            nParticles = nParticles + npar
            position.append(pos)
            vorticity.append(vort)
            radius.append(rad)
            # volume.append(vol)
            # birthstrength.append(bs)
        position = np.concatenate(position)
        vorticity = np.concatenate(vorticity)
        radius = np.concatenate(radius)
        # volume = np.concatenate(volume)
        # birthstrength = np.concatenate(birthstrength)
        # return [time, nParticles, position, vorticity, radius, volume, birthstrength, fileContent]
        return [time, nParticles, position, vorticity, radius, fileContent]

    def getWakeState(self, fileContent):
        """ getWakeState returns wake data from file
        """
        [nParticles, fileContent] = self.getInteger(fileContent)
        fileContent = fileContent[4:]   # Skipping 4 bit blank space
        [position, fileContent] = self.getMatrix(fileContent,nParticles)
        [vorticity, fileContent] = self.getMatrix(fileContent,nParticles)
        [radius, fileContent] = self.getVector(fileContent,nParticles)
        # [volume, fileContent] = self.getVector(fileContent,nParticles)
        # [birthstrength, fileContent] = self.getVector(fileContent,nParticles)
        # return [nParticles, position, vorticity, radius, volume, birthstrength, fileContent]
        return [nParticles, position, vorticity, radius, fileContent]

    def printData(self):
        """ printData prints wake properties
        """
        print(self.nTimesteps)
        print(self.time)
        print(self.nParticles)
        print(self.position)
        print(self.vorticity)
        print(self.radius)
        print(self.volume)
        print(self.birthstrength)

if __name__ == "__main__":
    rw = readWake("data/temp.wake")
    rw.printData()

