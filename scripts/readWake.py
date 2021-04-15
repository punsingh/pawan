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
        [self.nParticles, fileContent] = self.getInteger(fileContent)
        fileContent = fileContent[4:]   # Skipping 4 bit blank space
        [self.position, fileContent] = self.getMatrix(fileContent,self.nParticles)
        [self.vorticity, fileContent] = self.getMatrix(fileContent,self.nParticles)
        [self.radius, fileContent] = self.getVector(fileContent,self.nParticles)
        [self.volume, fileContent] = self.getVector(fileContent,self.nParticles)
        [self.birthstrength, fileContent] = self.getVector(fileContent,self.nParticles)
    
    def getInteger(self,fileContent):
        """ getInteger returns an integer 
        """
        ans = struct.unpack("i", fileContent[:4])[0]
        fileContent = fileContent[4:]
        return [ans, fileContent]
        # return struct.unpack("i", fileContent[:4])[0]
    
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
    rw = readWake("data/temp.wake")
    rw.printData()

