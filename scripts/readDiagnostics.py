#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 27 16:11:13 2023

@author: ge56beh
"""
import struct
import numpy as np
import argparse as ap


class readDiagnostics:

    def __init__(self,fileName):
        ''" Constructor """
        with open(fileName, mode='rb') as file:
            fileContent = file.read()
        self.fileSize = len(fileContent)
        self.time = []
        self.totalvorticity = []
        self.linearimpulse = []
        self.angularimpulse = []
        self.helicity = []
        self.enstrophy = []
        self.enstrophyF = []
        self.kineticenergy = []
        self.kineticenergyF = []
        self.centroidpos = []
        
        self.nu = struct.unpack("d", fileContent[:8])[0]
        fileContent = fileContent[8:]        
        while len(fileContent)>0:
            [t,tv,li,ai,h,e,ef,ke,kef,cp,fileContent] = self.getSimulationDiagnostics(fileContent)
            self.time.append(t)
            self.totalvorticity.append(tv)
            self.linearimpulse.append(li)
            self.angularimpulse.append(ai)
            self.helicity.append(h)
            self.enstrophy.append(e)
            self.enstrophyF.append(ef)
            self.kineticenergy.append(ke)
            self.kineticenergyF.append(kef)
            self.centroidpos.append(cp)
        self.nTimesteps = len(self.time)
        
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
    
    def getSimulationDiagnostics(self, fileContent):
        [time, fileContent] = self.getDouble(fileContent)
        [totalvorticity, fileContent] = self.getVector(fileContent,3)
        [linearimpulse, fileContent] = self.getVector(fileContent,3)
        [angularimpulse, fileContent] = self.getVector(fileContent,3)
        [helicity, fileContent] = self.getDouble(fileContent)
        [enstrophy, fileContent] = self.getDouble(fileContent)
        [enstrophyF, fileContent] = self.getDouble(fileContent)
        [kineticenergy, fileContent] = self.getDouble(fileContent)
        [kineticenergyF, fileContent] = self.getDouble(fileContent)
        [centroidpos, fileContent] = self.getDouble(fileContent)

        return [time, totalvorticity, linearimpulse, angularimpulse, helicity, 
                enstrophy, enstrophyF, kineticenergy, kineticenergyF, centroidpos, 
                fileContent]

if __name__ == "__main__":

    directry = "../data"
    filename = "temp.diagnosis"
    wake = readDiagnostics(f"{directry}/{filename}")
