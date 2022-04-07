__author__="Yang Bai"
__copyright__= "Copyright (C) 2021-present by M3 Group"
__version__ = "1.0"
__maintainer__ = "Yang Bai"
__email__ = "yangbai90@outlook.com"
__status__ = "development"
__date__ = "Dec 25, 2021"

import numpy as np
import matplotlib.pyplot as plt
import sys


class shape1d:
    def __init__(self,meshtype='edge2'):
        """
        Parameters
        ----------
        meshtype :  string
            the type of mesh, the default one is 'edge2'
        """
        self.nNodes=2
        self.shape_val=np.zeros(self.nNodes)
        self.shape_grad=np.zeros(self.nNodes)
        self.jacdet=1.0
        self.meshtype=meshtype
        self.setmeshtype(meshtype)
        self.update()
    def setmeshtype(self,meshtype):
        """
        set up the mesh type for shape1d class

        Parameters
        ----------
        meshtype : string
            the type of mesh you plan to use, it could be edge2, edge3, edge4
        """
        if 'edge2' in meshtype:
            self.nNodes=2
            self.meshtype='edge2'
        elif 'edge3' in meshtype:
            self.nNodes=3
            self.meshtype='edge3'
        elif 'edge4' in meshtype:
            self.nNodes=4
            self.meshtype='edge4'
        else:
            sys.exit('unsupported mesh type in shape1d->setmeshtype')
    def update(self):
        """
        update the shape function
        """
        self.shape_val=np.zeros(self.nNodes)
        self.shape_grad=np.zeros(self.nNodes)
        self.jacdet=1.0
    def getshpnumber(self):
        return self.nNodes
    ###########################################################
    def calc(self,xi,x,flag=True):
        """
        calculate the shape function value and its derivative w.r.t xi(local) or x(global)
        
        Parameters
        ----------
        xi : scalar
            the local cooridinate
        x : vector
            the global coordinate
        flag : boolean 
            True for the calculation based on global coordinate, otherwise, it use the local one
        """
        if self.nNodes==2:
            self.shape_val[1-1] = 0.5*(1.0-xi)
            self.shape_grad[1-1]=-0.5

            self.shape_val[2-1] = 0.5*(1.0+xi)
            self.shape_grad[2-1]= 0.5
        elif self.nNodes==3:
            self.shape_val[1-1] =0.5*xi*(xi-1.0)
            self.shape_grad[1-1]=0.5*(2*xi-1)

            self.shape_val[2-1] =-(xi+1.0)*(xi-1.0)
            self.shape_grad[2-1]=-2.0*xi

            self.shape_val[3-1] =0.5*xi*(xi+1.0)
            self.shape_grad[3-1]=0.5*(2.0*xi+1.0)
        elif self.nNodes==4:
            self.shape_val[1-1]=-(3.0*xi+1.0)*(3.0*xi-1.0)*(    xi-1.0)/16.0
            self.shape_grad[1-1]=-27.0*xi*xi/16.0+9.0*xi/8.0+ 1.0/16.0

            self.shape_val[2-1]=(3.0*xi+3.0)*(3.0*xi-1.0)*(3.0*xi-3.0)/16.0
            self.shape_grad[2-1]=-27.0*xi*xi/16.0+9.0*xi/8.0+ 1.0/16.0

            self.shape_val[3-1]=-(3.0*xi+3.0)*(3.0*xi+1.0)*(3.0*xi-3.0)/16.0
            self.shape_grad[3-1]=-81.0*xi*xi/16.0-9.0*xi/8.0+27.0/16.0

            self.shape_val[4-1]=(    xi+1.0)*(3.0*xi+1.0)*(3.0*xi-1.0)/16.0
            self.shape_grad[4-1]=27.0*xi*xi/16.0+9.0*xi/8.0- 1.0/16.0
        else:
            sys.exit('unsupported shape function calculation in shape1d')

        dxdxi=0.0
        for i in range(self.nNodes):
            dxdxi+=self.shape_grad[i]*x[i]
            
        self.jacdet=np.abs(dxdxi)
        if self.jacdet<1.0e-16:
            sys.exit('error: you have one singular 1d mesh !!!')
        if flag==True:
            for i in range(self.nNodes):
                self.shape_grad[i]=self.shape_grad[i]/self.jacdet

        return self.shape_val,self.shape_grad,self.jacdet
    def plot(self):
        """
        plot the 1d shape function
        """
        n=50
        xivec=np.linspace(-1.0,1.0,n)
        x=np.linspace(0.0,1.0,self.nNodes)
        shp=np.zeros((self.nNodes,n))
        for i in range(n):
            xi=xivec[i]
            shp[:,i],shpgrad,jac=self.calc(xi,x)
        plt.figure()
        for i in range(self.nNodes):
            plt.plot(xivec,shp[i,:],label=r'$N_{%d}$'%(i+1))
        plt.legend(fontsize=13)
        plt.xlabel(r'$\xi$',fontsize=15)
        plt.ylabel('shape function value',fontsize=15)
#################################################################################
class shape2d:
    def __init__(self,meshtype='quad4'):
        """
        Initialize the 2d shape function class 
        """
        self.nNodes=4
        self.shape_val=np.zeros(self.nNodes)
        self.shape_grad=np.zeros((self.nNodes,2))
        self.jacdet=1.0
        self.meshtype=meshtype
        self.setmeshtype(meshtype)
        self.update()
    def setmeshtype(self,meshtype):
        """
        set up the mesh type for shape2d class

        Parameters
        ----------
        meshtype: string
            the type of mesh, it should be quad4,quad9
        """
        if 'quad4' in meshtype:
            self.nNodes=4
            self.meshtype='quad4'
        elif 'quad9' in meshtype:
            self.nNodes=9
            self.meshtype='quad9'
        else:
            sys.exit('unsupported mesh type in shape2d->setmeshtype')
    def update(self):
        """
        update the shape function vector
        """
        self.shape_val=np.zeros(self.nNodes)
        self.shape_grad=np.zeros((self.nNodes,2))
        self.jacdet=1.0
    def getshpnumber(self):
        """
        return the number of shape function
        """
        return self.nNodes
    ###########################################################
    def calc(self,xi,eta,x,y,flag=True):
        """
        calculate the shape function value and its derivative w.r.t xi(local) or x(global)

        Parameters
        ----------
        xi : double
            the local cooridinate
        eta : double
            the local cooridinate
        x : double 
            the global coordinate
        y : double
            the global coordinate
        flag : boolean
            True for the calculation based on global coordinate, otherwise, it use the local one
        """
        if self.nNodes==4:
            self.shape_val[0]=(1.0-xi)*(1.0-eta)/4.0
            self.shape_val[1]=(1.0+xi)*(1.0-eta)/4.0
            self.shape_val[2]=(1.0+xi)*(1.0+eta)/4.0
            self.shape_val[3]=(1.0-xi)*(1.0+eta)/4.0
            
            self.shape_grad[0,1-1]= (eta-1.0)/4.0  
            self.shape_grad[0,2-1]= (xi -1.0)/4.0  
  
            self.shape_grad[1,1-1]= (1.0-eta)/4.0  
            self.shape_grad[1,2-1]=-(1.0+xi )/4.0  
  
            self.shape_grad[2,1-1]= (1.0+eta)/4.0  
            self.shape_grad[2,2-1]= (1.0+xi )/4.0  
  
            self.shape_grad[3,1-1]=-(1.0+eta)/4.0
            self.shape_grad[3,2-1]= (1.0-xi )/4.0
        elif self.nNodes==9:
            self.shape_val[0]=(xi*xi-xi )*(eta*eta-eta)/4.0
            self.shape_val[1]=(xi*xi+xi )*(eta*eta-eta)/4.0
            self.shape_val[2]=(xi*xi+xi )*(eta*eta+eta)/4.0
            self.shape_val[3]=(xi*xi-xi )*(eta*eta+eta)/4.0
            self.shape_val[4]=(1.0-xi*xi)*(eta*eta-eta)/2.0
            self.shape_val[5]=(xi*xi+xi )*(1.0-eta*eta)/2.0
            self.shape_val[6]=(1.0-xi*xi)*(eta*eta+eta)/2.0
            self.shape_val[7]=(xi*xi-xi )*(1.0-eta*eta)/2.0
            self.shape_val[8]=(1.0-xi*xi)*(1.0-eta*eta)   

            self.shape_grad[0,1-1]=(2.0*xi-1.0)*(eta*eta-eta)/4.0  
            self.shape_grad[0,2-1]=(xi*xi-xi  )*(2.0*eta-1.0)/4.0  
  
            self.shape_grad[1,1-1]=(2.0*xi+1.0)*(eta*eta-eta)/4.0  
            self.shape_grad[1,2-1]=(xi*xi+xi  )*(2.0*eta-1.0)/4.0  
     
            self.shape_grad[2,1-1]=(2.0*xi+1.0)*(eta*eta+eta)/4.0  
            self.shape_grad[2,2-1]=(xi*xi+xi  )*(2.0*eta+1.0)/4.0  
   
            self.shape_grad[3,1-1]=(2.0*xi-1.0)*(eta*eta+eta)/4.0  
            self.shape_grad[3,2-1]=(xi*xi-xi  )*(2.0*eta+1.0)/4.0  

            self.shape_grad[4,1-1]=-xi*(eta*eta-eta)
            self.shape_grad[4,2-1]=(1.0-xi*xi )*(2.0*eta-1.0)/2.0
        
            self.shape_grad[5,1-1]=(2.0*xi+1.0)*(1.0-eta*eta)/2.0
            self.shape_grad[5,2-1]=-(xi*xi+xi )*eta
      
            self.shape_grad[6,1-1]=-xi*(eta*eta+eta)
            self.shape_grad[6,2-1]=(1.0-xi*xi )*(2.0*eta+1.0)/2.0
  
            self.shape_grad[7,1-1]=(2.0*xi-1.0)*(1.0-eta*eta)/2.0
            self.shape_grad[7,2-1]=-(xi*xi-xi )*eta
  
            self.shape_grad[8,1-1]=-2.0*xi*(1.0-eta*eta)
            self.shape_grad[8,2-1]=-2.0*eta*(1.0-xi*xi)
        else:
            sys.exit('unsupported shape function calculation in shape1d')

        dxdxi=0.0;dxdeta=0.0
        dydxi=0.0;dydeta=0.0
        for i in range(self.nNodes):
            dxdxi +=self.shape_grad[i,0]*x[i]
            dxdeta+=self.shape_grad[i,1]*x[i]

            dydxi +=self.shape_grad[i,0]*y[i]
            dydeta+=self.shape_grad[i,1]*y[i]
        
        jac=np.array([[dxdxi,dydxi],[dxdeta,dydeta]])
        self.jacdet=np.linalg.det(jac)
        xjac=np.linalg.inv(jac)
        if self.jacdet<1.0e-16:
            sys.exit('error: you have one singular 1d mesh !!!')
        if flag==True:
            for i in range(self.nNodes):
                temp1=self.shape_grad[i,1-1]*xjac[0,0]+self.shape_grad[i,2-1]*xjac[0,1]
                temp2=self.shape_grad[i,1-1]*xjac[1,0]+self.shape_grad[i,2-1]*xjac[1,1]
                self.shape_grad[i,1-1]=temp1
                self.shape_grad[i,2-1]=temp2

        return self.shape_val,self.shape_grad,self.jacdet
    def plot(self):
        """
        plot the 2d shape function
        """
        n=20
        xivec=np.linspace(-1.0,1.0,n)
        etavec=np.linspace(-1.0,1.0,n)
        x=np.array([0.0,1.0,1.0,0.0, 0.5,1.0,0.5,0.0,0.5])
        y=np.array([0.0,0.0,1.0,1.0, 0.0,0.5,1.0,0.5,0.5])
        shp=np.zeros((self.nNodes,n,n))
        X,Y=np.meshgrid(xivec,etavec)
        for i in range(n):
            for j in range(n):
                xi=xivec[i];eta=etavec[j]
                shp[:,i,j],shpgrad,jac=self.calc(xi,eta,x,y)
        fig=plt.figure()
        ax = fig.add_subplot(111, projection="3d")
        for i in range(self.nNodes):
            Z=shp[i,:,:]
            ax.plot_surface(X,Y,Z,label=r'$N_{%d}$'%(i+1))
            #ax.contour(X, Y, Z, 10, lw=3, cmap="autumn_r", linestyles="solid", offset=-1)
            #ax.contour(X, Y, Z, 10, lw=3, colors="k", linestyles="solid")
        ax.set_xlabel(r'$\xi$',fontsize=15)
        ax.set_ylabel(r'$\eta$',fontsize=15)

