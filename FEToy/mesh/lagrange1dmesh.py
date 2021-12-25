__author__="Yang Bai"
__copyright__= "Copyright (C) 2021-present by M3 Group"
__version__ = "1.0"
__maintainer__ = "Yang Bai"
__email__ = "yangbai90@outlook.com"
__status__ = "development"
__date__ = "Dec 19, 2021"

import numpy as np
import matplotlib.pyplot as plt
import sys

class mesh1d:
    def __init__(self,xmin=0.0,xmax=1.0,nx=5,meshtype='edge2'):
        """
        Initialize the mesh class 
        :param xmin: the left point along x-axis
        :type xmin: double
        :param xmax: the right point along x-axis
        :type xmax: double
        :param nx: the number of element
        :type nx: int
        :param meshtype: the type of mesh you want to generated, it should be edge2, edge3, edge4
        :type meshtype: string
        """
        self.meshtype=meshtype
        self.dim=1
        self.xmin=xmin
        self.xmax=xmax
        self.nx=nx
        self.order=1
        self.nodes=0
        self.nodesperelement=2
        self.elements=0
        self.meshtype=meshtype
        self.setmeshtype(meshtype)
    def setnx(self,nx):
        """
        set up the number of element we want to generate
        :param nx: the number of element
        :type nx: int
        """
        self.nx=nx
    def setmeshtype(self,meshtype):
        """
        set up the type of mesh we want to use
        :param meshtype: the type name should be edge2,edge3,edge4
        :type meshtype: string
        """
        if 'edge2' in meshtype:
            self.order=1
            self.nodesperelement=2
        elif 'edge3' in meshtype:
            self.order=2
            self.nodesperelement=3
        elif 'edge4' in meshtype:
            self.order=3
            self.nodesperelement=4
        else:
            sys.exit('sorry, unsuported 1d mesh type in FEToy!')
    def setdomainsize(self,xmin,xmax):
        """
        set up the domain size
        :param xmin: the left point of x-axis
        :type xmin: double
        :param xmax: the right point of x-axis
        :type xmax: double
        """
        self.xmin=xmin
        self.xmax=xmax
    def update(self):
        self.createmesh()
    def createmesh(self):
        """
        generate the lagrange mesh in 1d case
        """
        self.elements=self.nx
        self.nodes=self.elements*self.order+1
        dx=(self.xmax-self.xmin)/(self.nodes-1)
        self.nodecoords=np.zeros(self.nodes)
        for i in range(self.nodes):
            self.nodecoords[i]=self.xmin+i*dx

        self.elementconn=np.zeros((self.elements,self.nodesperelement),dtype=np.int16)
        for e in range(self.elements):
            for j in range(self.nodesperelement):
                self.elementconn[e,j]=e*self.order+j+1
    #####################################################
    def printnodes(self):
        print('*** node coordinates of the mesh (total nodes=%d, nodes per element=%d)'%(self.nodes,self.nodesperelement))
        for i in range(self.nodes):
            str='%6d-th node: x=%14.6e'%(i+1,self.nodecoords[i])
            print(str)
    def printelements(self):
        print('*** element connectivity information(bulk elements=%d)'%(self.elements))
        for e in range(self.elements):
            str='%6d-th element'%(e+1)
            for i in range(self.nodesperelement):
                str+='%5d '%(self.elementconn[e,i])
            print(str)
    ######################################################
    def plotmesh(self,withnode=False,withnodeid=False):
        y=np.zeros(self.nodes)
        plt.figure()
        plt.plot(self.nodecoords,y,'k')
        if withnode:
            plt.plot(self.nodecoords,y,'r*')
        if withnodeid:
            for i in range(self.nodes):
                x=self.nodecoords[i]
                plt.text(x,0.0,'%d'%(i+1))
        plt.xlabel('X',fontsize=14)
        plt.ylabel('Y',fontsize=14)
