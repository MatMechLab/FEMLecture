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

class mesh2d:
    def __init__(self,xmin=0.0,xmax=1.0,ymin=0.0,ymax=1.0,nx=5,ny=5,meshtype='quad4'):
        """
        Initialize the mesh class 
        :param xmin: the left point along x-axis
        :type xmin: double
        :param xmax: the right point along x-axis
        :type xmax: double
        :param ymin: the bottom point along y-axis
        :type ymin: double
        :param ymax: the top point along y-axis
        :type ymax: double
        :param nx: the number of element in x direction
        :type nx: int
        :param ny: the number of element in y direction
        :type ny: int
        :param meshtype: the type of mesh you want to generated, it should be quad4, edge8, edge9
        :type meshtype: string
        """
        self.meshtype=meshtype
        self.dim=2
        self.xmin=xmin
        self.xmax=xmax
        self.ymin=ymin
        self.ymax=ymax 
        self.nx=nx
        self.ny=ny
        self.order=1
        self.nodes=0
        self.nodesperelement=4
        self.elements=0
        self.meshtype=meshtype
        self.setmeshtype(meshtype)
    def setnx(self,nx):
        """
        set up the number of element we want to generate along x direction
        :param nx: the number of element
        :type nx: int
        """
        self.nx=nx
    def setny(self,ny):
        """
        set up the number of element we want to generate along y direction
        :param ny: the number of element
        :type ny: int
        """
        self.ny=ny
    def setmeshtype(self,meshtype):
        """
        set up the type of mesh we want to use
        :param meshtype: the type name should be qud4,quad8,quad9
        :type meshtype: string
        """
        if 'quad4' in meshtype:
            self.order=1
            self.nodesperelement=4
            self.meshtype='quad4'
        elif 'quad9' in meshtype:
            self.order=2
            self.nodesperelement=9
            self.meshtype='quad9'
        else:
            sys.exit('sorry, unsuported 2d mesh type in FEToy!')
    def setdomainsize(self,xmin,xmax,ymin,ymax):
        """
        set up the domain size
        :param xmin: the left point of x-axis
        :type xmin: double
        :param xmax: the right point of x-axis
        :type xmax: double
        :param ymin: the bottom point of y-axis
        :type ymin: double
        :param ymax: the top point of y-axis
        :type ymax: double
        """
        self.xmin=xmin
        self.xmax=xmax
        self.ymin=ymin
        self.ymax=ymax
    def update(self):
        self.createmesh()
    def createmesh(self):
        """
        generate the lagrange mesh in 2d case
        """
        if 'quad4' in self.meshtype:
            self.elements=self.nx*self.ny
            self.nodes=(self.nx+1)*(self.ny+1)
            dx=(self.xmax-self.xmin)/(self.nx-1)
            dy=(self.ymax-self.ymin)/(self.ny-1)
            self.nodecoords=np.zeros((self.nodes,2))

            # for bc nodes
            leftnodes=np.zeros(self.ny+1,dtype=np.int16)
            rightnodes=np.zeros(self.ny+1,dtype=np.int16)
            bottomnodes=np.zeros(self.nx+1,dtype=np.int16)
            topnodes=np.zeros(self.nx+1,dtype=np.int16)
            for j in range(self.ny+1):
                for i in range(self.nx+1):
                    k=j*(self.nx+1)+i
                    self.nodecoords[k,0]=self.xmin+i*dx
                    self.nodecoords[k,1]=self.ymin+j*dx
                    if i==0:
                        # for left nodes
                        leftnodes[j]=k
                    if i==self.nx:
                        # for right nodes
                        rightnodes[j]=k
                    if j==0:
                        # for bottom nodes
                        bottomnodes[i]=k
                    if j==self.ny:
                        # for top nodes
                        topnodes[i]=k
            self.bcnodeids={'left':leftnodes,'right':rightnodes,'bottom':bottomnodes,'top':topnodes}
            ###########################################
            self.elementconn=np.zeros((self.elements,self.nodesperelement),dtype=np.int16)
            # for bc elements
            leftconn=np.zeros((self.ny,2),dtype=np.int16)
            rightconn=np.zeros((self.ny,2),dtype=np.int16)
            bottomconn=np.zeros((self.nx,2),dtype=np.int16)
            topconn=np.zeros((self.nx,2),dtype=np.int16)
            for j in range(1,self.ny+1):
                for i in range(1,self.nx+1):
                    e=(j-1)*self.nx+i-1
                    i1=(j-1)*(self.nx+1)+i
                    i2=i1+1
                    i3=i2+self.nx+1
                    i4=i3-1
                    self.elementconn[e,0]=i1-1
                    self.elementconn[e,1]=i2-1
                    self.elementconn[e,2]=i3-1
                    self.elementconn[e,3]=i4-1
                    # 4 +-----+ 3
                    #   |     |
                    #   |     |
                    # 1 +-----+ 2
                    if i==1:
                        # for left side
                        leftconn[j-1,0]=i4
                        leftconn[j-1,1]=i1
                    if i==self.nx:  
                        # for right side
                        rightconn[j-1,0]=i2
                        rightconn[j-1,1]=i3
                    if j==1:  
                        # for bottom side
                        bottomconn[i-1,0]=i1
                        bottomconn[i-1,1]=i2
                    if j==self.ny:
                        # for top side
                        topconn[i-1,0]=i3
                        topconn[i-1,1]=i4
            self.bcconn={'left':leftconn,'right':rightconn,'bottom':bottomconn,'top':topconn}
        elif 'quad9' in self.meshtype:
            self.elements=self.nx*self.ny
            self.nodes=(2*self.nx+1)*(2*self.ny+1)
            dx=(self.xmax-self.xmin)/(2*self.nx)
            dy=(self.ymax-self.ymin)/(2*self.ny)
            self.nodecoords=np.zeros((self.nodes,2))
            # for bc nodes
            leftnodes=np.zeros(2*self.ny+1,dtype=np.int16)
            rightnodes=np.zeros(2*self.ny+1,dtype=np.int16)
            bottomnodes=np.zeros(2*self.nx+1,dtype=np.int16)
            topnodes=np.zeros(2*self.nx+1,dtype=np.int16)
            for j in range(2*self.ny+1):
                for i in range(2*self.nx+1):
                    k=j*(2*self.nx+1)+i
                    self.nodecoords[k,0]=self.xmin+i*dx
                    self.nodecoords[k,1]=self.ymin+j*dy
                    if i==0:
                        # for left nodes
                        leftnodes[j]=k
                    if i==2*self.nx:
                        # for right nodes
                        rightnodes[j]=k
                    if j==0:
                        # for bottom nodes
                        bottomnodes[i]=k
                    if j==2*self.ny:
                        # for top nodes
                        topnodes[i]=k
            self.bcnodeids={'left':leftnodes,'right':rightnodes,'bottom':bottomnodes,'top':topnodes}
            #######################################        
            self.elementconn=np.zeros((self.elements,self.nodesperelement),dtype=np.int16)
            # for bc elements
            leftconn=np.zeros((self.ny,3),dtype=np.int16)
            rightconn=np.zeros((self.ny,3),dtype=np.int16)
            bottomconn=np.zeros((self.nx,3),dtype=np.int16)
            topconn=np.zeros((self.nx,3),dtype=np.int16)
            for j in range(1,self.ny+1):
                for i in range(1,self.nx+1):
                    e=(j-1)*self.nx+i-1
                    i1=(j-1)*2*(2*self.nx+1)+2*i-1
                    i2=i1+2
                    i3=i2+2*(2*self.nx+1)
                    i4=i3-2
                    i5=i1+1
                    i6=i2+(2*self.nx+1)
                    i7=i3-1
                    i8=i1+(2*self.nx+1);
                    i9=i8+1
                    self.elementconn[e,1-1]=i1-1
                    self.elementconn[e,2-1]=i2-1
                    self.elementconn[e,3-1]=i3-1
                    self.elementconn[e,4-1]=i4-1
                    self.elementconn[e,5-1]=i5-1
                    self.elementconn[e,6-1]=i6-1
                    self.elementconn[e,7-1]=i7-1
                    self.elementconn[e,8-1]=i8-1
                    self.elementconn[e,9-1]=i9-1
                    ##
                    # 4 +---7---+ 3
                    #   |       |
                    # 8 +   9   + 6
                    #   |       |
                    # 1 +---5---+ 2
                    if i==1:
                        # for left side
                        leftconn[j-1,0]=i4
                        leftconn[j-1,1]=i8
                        leftconn[j-1,2]=i1
                    if i==self.nx:
                        # for right side
                        rightconn[j-1,0]=i2
                        rightconn[j-1,1]=i6
                        rightconn[j-1,2]=i3
                    if j==1:
                        # for bottom side
                        bottomconn[i-1,0]=i1
                        bottomconn[i-1,1]=i5
                        bottomconn[i-1,2]=i2
                    if j==self.ny:
                        # for top side
                        topconn[i-1,0]=i3
                        topconn[i-1,1]=i7
                        topconn[i-1,2]=i4
            self.bcconn={'left':leftconn,'right':rightconn,'bottom':bottomconn,'top':topconn}
    #####################################################
    def printnodes(self):
        print('*** node coordinates of the mesh (total nodes=%d, nodes per element=%d)'%(self.nodes,self.nodesperelement))
        for i in range(self.nodes):
            str='%6d-th node: x=%14.6e,y=%14.6e'%(i+1,self.nodecoords[i,0],self.nodecoords[i,1])
            print(str)
    def printelements(self):
        print('*** element connectivity information(bulk elements=%d)'%(self.elements))
        for e in range(self.elements):
            str='%6d-th element :'%(e+1)
            for i in range(self.nodesperelement):
                str+='%5d '%(self.elementconn[e,i])
            print(str)
    ######################################################
    def plotmesh(self,withnode=False,withnodeid=False):
        plt.plot()
        for e in range(self.elements):
            conn=self.elementconn[e,:]
            if 'quad4' in self.meshtype:
                conn=np.append(conn,conn[0])
            elif 'quad9' in self.meshtype:
                conn=np.array([conn[1-1],conn[5-1],conn[2-1],conn[6-1],conn[3-1],conn[7-1],conn[4-1],conn[8-1],conn[1-1]])
            x=self.nodecoords[conn,0]
            y=self.nodecoords[conn,1]
            plt.plot(x,y,'k')
        if withnode:
            plt.plot(self.nodecoords[:,0],self.nodecoords[:,1],'r*')
        if withnodeid:
            for i in range(self.nodes):
                x=self.nodecoords[i,0];y=self.nodecoords[i,1]
                str='%d'%(i+1)
                plt.text(x,y,str)
        plt.xlabel('X',fontsize=14)
        plt.ylabel('Y',fontsize=14)
