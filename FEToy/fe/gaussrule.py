__author__="Yang Bai"
__copyright__= "Copyright (C) 2021-present by M3 Group"
__version__ = "1.0"
__maintainer__ = "Yang Bai"
__email__ = "yangbai90@outlook.com"
__status__ = "development"
__date__ = "Dec 21, 2021"


import numpy as np
import matplotlib.pyplot as plt
import sys

class gausspoint1d:
    def __init__(self,ngp=1):
        """
        Initialize the gauss point class
        :param ngp: the number of gauss points
        :type ngp: int
        """
        self.ngp=ngp
        self.gpcoords=np.zeros((ngp,2)) # 0-> weight, 1-> xi

    def setnum(self,num):
        """
        Set up the gauss point number
        :param num: number of gauss point one want to generate
        :type num: int
        """
        self.ngp=num
    def getgausspointsnumber(self):
        """
        Get the number of gauss points
        """
        return self.ngp
    def creategausspoint(self):
        """
        Generate the gauss points 
        """
        self.gpcoords=np.zeros((self.ngp,2))
        if self.ngp==1:
            self.gpcoords[0,0]=2.0
            self.gpcoords[0,1]=0.0
        elif self.ngp==2:
            self.gpcoords[0,0]=1.0
            self.gpcoords[0,1]=-np.sqrt(1.0/3.0)

            self.gpcoords[1,0]=1.0
            self.gpcoords[1,1]= np.sqrt(1.0/3.0)
        elif self.ngp==3:
            self.gpcoords[0,0]=5.0/9.0
            self.gpcoords[0,1]=-np.sqrt(3.0/5.0);

            self.gpcoords[1,0]=8.0/9.0
            self.gpcoords[1,1]= 0.0

            self.gpcoords[2,0]=5.0/9.0
            self.gpcoords[2,1]= np.sqrt(3.0/5.0)
        elif self.ngp==4:
            t=np.sqrt(4.8)
            w=1.0/3.0/t

            self.gpcoords[0,1]=-np.sqrt((3.0+t)/7.0)
            self.gpcoords[1,1]=-np.sqrt((3.0-t)/7.0)
            self.gpcoords[2,1]= np.sqrt((3.0-t)/7.0)
            self.gpcoords[3,1]= np.sqrt((3.0+t)/7.0)

            self.gpcoords[0,0]=0.5-w
            self.gpcoords[1,0]=0.5+w
            self.gpcoords[2,0]=0.5+w
            self.gpcoords[3,0]=0.5-w
        elif self.ngp==5:
            t=np.sqrt(1120.0)
            self.gpcoords[1-1, 1] = -np.sqrt((70.0 + t) / 126.0)
            self.gpcoords[2-1, 1] = -np.sqrt((70.0 - t) / 126.0)
            self.gpcoords[3-1, 1] = 0.0;
            self.gpcoords[4-1, 1] = np.sqrt((70.0 - t) / 126.0)
            self.gpcoords[5-1, 1] = np.sqrt((70.0 + t) / 126.0)
     
            self.gpcoords[1-1, 0] = (21.0 * t + 117.60) / (t * (70.0 + t))
            self.gpcoords[2-1, 0] = (21.0 * t - 117.60) / (t * (70.0 - t))
            self.gpcoords[3-1, 0] = 2.0 * (1.0 - self.gpcoords[1-1, 0] -self.gpcoords[2-1, 0])
            self.gpcoords[4-1, 0] = self.gpcoords[2-1, 0]
            self.gpcoords[5-1, 0] = self.gpcoords[1-1, 0]
        else:
            sys.exit("unsuported gauss point number in 1d case")
    def update(self):
        """
        Update/re-generate the gauss point number
        """
        self.creategausspoint()
    def print(self):
        for i in range(self.ngp):
            str='%d-th gauss point: xi=%14.5e, weight=%14.5e'%(i+1,self.gpcoords[i,1],self.gpcoords[i,0])
            print(str)
    def getIthweight(self,i):
        """
        Get the i-th gauss point's weight
        """
        if i>self.ngp:
            sys.exit('%d is out of 1d gauss points range !'%(i))
        return self.gpcoords[i,0]
    def getIthcoord(self,i):
        """
        Get the i-th gauss point's coordinate
        """
        if i>self.ngp:
            sys.exit('%d is out of 1d gauss points range !'%(i))
        return self.gpcoords[i,1]
    def plot(self):
        plt.plot(self.gpcoords[:,1],self.gpcoords[:,0],'k*')
        plt.xlabel(r'$\xi$',fontsize=16)
        plt.ylabel(r'$w$',fontsize=16)
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
###########################################################################
class gausspoint2d:
    def __init__(self,ngp):
        """
        Initialize the 2d gauss points
        """
        self.ngp=ngp
        self.ngp2=ngp*ngp
        self.gpcoords=np.zeros((self.ngp2,3)) # 0->weight,1->xi,2->eta
    def setnum(self,ngp):
        """
        Set up the gauss point number in single direction(here we use the tensor product)
        :param ngp: number of gauss point in single direction
        :type ngp: int
        """
        self.ngp=ngp
        self.ngp2=ngp*ngp
        if ngp>5:
            sys.exit("unsupported gauss points number in 2d case!")
    def creategausspoint(self):
        """
        Generate the 2d gauss points
        """
        if self.ngp>5:
            sys.exit('can not generate 2d gauss point, unsupported gauss point number (must <=5) !!!')
        self.gpcoords=np.zeros((self.ngp2,3))
        gp1d=gausspoint1d(self.ngp)
        gp1d.creategausspoint()
        k=0
        for i in range(self.ngp):
            for j in range(self.ngp):
                self.gpcoords[k,0]=gp1d.gpcoords[i,0]*gp1d.gpcoords[j,0]
                self.gpcoords[k,1]=gp1d.gpcoords[i,1]
                self.gpcoords[k,2]=gp1d.gpcoords[j,1]
                k+=1
    def print(self):
        for i in range(self.ngp2):
            str='%d-th gauss point: xi=%14.5e, eta=%14.5e, weight=%14.5e'%(i+1,self.gpcoords[i,1],self.gpcoords[i,2],self.gpcoords[i,0])
            print(str)
    def update(self):
        """
        Update the gauss points 
        """
        self.creategausspoint()
    def getgausspointsnumber(self):
        return self.ngp2
    def getIthcoords(self,i):
        """
        get the i-th gauss point's coordinate
        :param i: i the node index
        :type i: int
        """
        if i>self.ngp2:
            sys.exit('i=%d is out of 2d gauss point range'%(i))
        return self.gpcoords[i,1],self.gpcoords[i,2]
    def getIthweight(self,i):
        """
        Get the i-th gauss point's weight
        :param i: the node index
        :type i: int
        """
        if i>self.ngp:
            sys.exit('i=%d is out of 2d gauss point range'%(i))
        return self.gpcoords[i,0]
    def plot(self):
        fig=plt.figure()
        plt.clf()
        ax=fig.add_subplot(projection='3d')
        ax.scatter(self.gpcoords[:,1],self.gpcoords[:,2],self.gpcoords[:,0])
        ax.set_xlabel(r'$\xi$',fontsize=14)
        ax.set_ylabel(r'$\eta$',fontsize=14)
        ax.set_zlabel(r'$w$',fontsize=14)
