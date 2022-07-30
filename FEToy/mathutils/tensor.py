__author__="Yang Bai"
__copyright__= "Copyright (C) 2022-present by M3 Group"
__version__ = "1.0"
__maintainer__ = "Yang Bai"
__email__ = "yangbai90@outlook.com"
__status__ = "development"
__date__ = "July 30, 2022"

import imp
import numpy as np
import matplotlib.pyplot as plt
import sys
# from vector import vector
from FEToy.mathutils.vector import vector

class rank2tensor:
    def __init__(self,dim=2,value=0.0):
        """
        Initialize the tensor

        Parameters
        ----------
        dim : int
            the dimension or size of current tensor
        value : double
            default value for current tensor
        """
        self.dim=dim
        if dim>3:
            sys.exit('error: dim=%d is invalid for a rank-2 tensor'%(dim))
        self.vals=np.zeros((self.dim,self.dim))
        self.vals[:,:]=1.0*value
    def setToRandom(self):
        """
        Set current rank-2 tensor to random value

        Parameters
        ----------
        """
        for i in range(self.dim):
            for j in range(self.dim):
                self.vals[i,j]=np.random.rand()
    def setToIdentity(self):
        """
        Set current rank-2 tensor to identity tensor

        Parameters
        ----------
        """
        for i in range(self.dim):
            for j in range(self.dim):
                self.vals[i,j]=1.0*(i==j)
    def det(self):
        """
        get the determinte of current tensor

        Parameters
        ----------
        """
        return np.linalg.det(self.vals)
    def trace(self):
        """
        get the trace of current tensor

        Parameters
        ----------
        """
        if self.dim==2:
            return self.vals[0,0]+self.vals[1,1]
        else:
            return self.vals[0,0]+self.vals[1,1]+self.vals[2,2]
    def dev(self):
        """
        get the deviatoric part of current tensor

        Parameters
        ----------
        """
        temp=rank2tensor(self.dim)
        tr=self.trace()
        for i in range(self.dim):
            for j in range(self.dim):
                temp.vals[i,j]=self.vals[i,j]-(1.0/3.0)*tr*(i==j)
        return temp
    def doublecontraction(self,a):
        """
        get the double contraction result between two tensor

        Parameters
        ----------
        a : rank2 tensor
            right hand side tensor
        """
        if isinstance(a,rank2tensor):
            if self.dim==a.dim:
                sum=0.0
                for i in range(self.dim):
                    for j in range(self.dim):
                        sum+=self.vals[i,j]*a.vals[i,j]
                return sum
            else:
                sys.exit('error: the right hand side value should be either a tensor or a scalar for \'-\' operator')
        else:
            sys.exit('error: double contraction should be used for two tensors !!!')
    def ithRow(self,i):
        """
        get i-th row of current tensor

        Parameters
        ----------
        i : integer
            row index
        """
        temp=vector(self.dim)
        for j in range(self.dim):
            temp.vals[j]=self.vals[i,j]
        return temp
    def ithCol(self,i):
        """
        get i-th column of current tensor

        Parameters
        ----------
        i : integer
            col index
        """
        temp=vector(self.dim)
        for j in range(self.dim):
            temp.vals[j]=self.vals[j,i]
        return temp
    def print(self,str=None):
        """
        Print the rank-2 tensor

        Parameters
        ----------
        str : string
            the content to be printed
        """
        if self.dim==2:
            if str is None:
                print('tensor components are: |%14.5e, %14.5e|'%(self.vals[0,0],self.vals[0,1]))
                print('                       |%14.5e, %14.5e|'%(self.vals[1,0],self.vals[1,1]))
            else:
                print(str+', tensor components are: ')
                print('  |%14.5e, %14.5e|'%(self.vals[0,0],self.vals[0,1]))
                print('  |%14.5e, %14.5e|'%(self.vals[1,0],self.vals[1,1]))
        else:
            if str is None:
                print('tensor components are: |%14.5e, %14.5e %14.5e|'%(self.vals[0,0],self.vals[0,1],self.vals[0,2]))
                print('                       |%14.5e, %14.5e %14.5e|'%(self.vals[1,0],self.vals[1,1],self.vals[1,2]))
                print('                       |%14.5e, %14.5e %14.5e|'%(self.vals[2,0],self.vals[2,1],self.vals[2,2]))
            else:
                print(str+', tensor components are: ')
                print('  |%14.5e, %14.5e %14.5e|'%(self.vals[0,0],self.vals[0,1],self.vals[0,2]))
                print('  |%14.5e, %14.5e %14.5e|'%(self.vals[1,0],self.vals[1,1],self.vals[1,2]))
                print('  |%14.5e, %14.5e %14.5e|'%(self.vals[2,0],self.vals[2,1],self.vals[2,2]))
    #################################################
    ### for operators
    #################################################
    def __getitem__(self,pos):
        """
        Get the i,j-th element of current vector

        Parameters
        ----------
        pos : int
            the integer index pair
        """
        i,j=pos
        if isinstance(i,int) and isinstance(j,int):
            if i<0 or i>self.dim-1:
                sys.exit('i=%d is out of range for a tensor element access !!!'%(i))
            if j<0 or j>self.dim-1:
                sys.exit('j=%d is out of range for a tensor element access !!!'%(j))
            return self.vals[i,j]
    def __setitem__(self,pos,val):
        """
        Set the i,j-th element of current tensor

        Parameters
        ----------
        pos : int pair
            the integer index pair
        val : float
            the float scalar 
        """
        i,j=pos
        if isinstance(i,int) and isinstance(j,int):
            if i<0 or i>self.dim-1:
                sys.exit('i=%d is out of range for a tensor element access !!!'%(i))
            if j<0 or j>self.dim-1:
                sys.exit('j=%d is out of range for a tensor element access !!!'%(i))
            self.vals[i,j]=val
    def __str__(self):
        """
        Print() method overload
        ----------
        """
        if self.dim==2:
            return "(%14.5e, %14.5e)\n(%14.5e, %14.5e)"%(self.vals[0,0],self.vals[0,1],
                                                         self.vals[1,0],self.vals[1,1])
        else:
            return "(%14.5e, %14.5e, %14.5e)\n(%14.5e, %14.5e, %14.5e)\n(%14.5e, %14.5e, %14.5e)"%(self.vals[0,0],self.vals[0,1],self.vals[0,2],
                                                                                                   self.vals[1,0],self.vals[1,1],self.vals[1,2],
                                                                                                   self.vals[2,0],self.vals[2,1],self.vals[2,2])
    # for '+'
    def __add__(self,a):
        """
        + operator overload

        Parameters
        ----------
        a : tensor or scalar
            right hand side value
        """
        if isinstance(a,rank2tensor):
            if self.dim==a.dim:
                temp=rank2tensor(self.dim,0.0)
                for i in range(self.dim):
                    for j in range(self.dim):
                        temp.vals[i,j]=self.vals[i,j]+a.vals[i,j]
                return temp
            else:
                sys.exit('error: \'+\' operator must be applied to two tensors with the same size !!!')
        else:
            if isinstance(a,int) or isinstance(a,float):
                # for scalar
                temp=rank2tensor(self.dim,0.0)
                for i in range(self.dim):
                    for j in range(self.dim):
                        temp.vals[i,j]=self.vals[i,j]+a
                return temp
            else:
                sys.exit('error: the right hand side value should be either a tensor or a scalar for \'+\' operator')
    # for '-'
    def __sub__(self,a):
        """
        - operator overload

        Parameters
        ----------
        a : temspr or scalar
            right hand side value
        """
        if isinstance(a,rank2tensor):
            if self.dim==a.dim:
                temp=rank2tensor(self.dim,0.0)
                for i in range(self.dim):
                    for j in range(self.dim):
                        temp.vals[i,j]=self.vals[i,j]-a.vals[i,j]
                return temp
            else:
                sys.exit('error: \'-\' operator must be applied to two tensors with the same size !!!')
        else:
            if isinstance(a,int) or isinstance(a,float):
                # for scalar
                temp=rank2tensor(self.dim,0.0)
                for i in range(self.dim):
                    for j in range(self.dim):
                        temp.vals[i,j]=self.vals[i,j]-a
                return temp
            else:
                sys.exit('error: the right hand side value should be either a tensor or a scalar for \'-\' operator')
    # for '*'
    def __mul__(self,a):
        """
        * operator overload

        Parameters
        ----------
        a : tensor or scalar
            right hand side value
        """
        if isinstance(a,rank2tensor):
            if self.dim==a.dim:
                temp=rank2tensor(self.dim,0.0)
                for i in range(self.dim):
                    for j in range(self.dim):
                        for k in range(self.dim):
                            temp.vals[i,j]+=self.vals[i,k]*a.vals[k,j]
                return temp
            else:
                sys.exit('error: \'*\' operator must be applied to two tensors with the same size !!!')
        elif isinstance(a,vector):
            if self.dim==a.dim:
                temp=vector(self.dim,0.0)
                for i in range(self.dim):
                    for j in range(self.dim):
                        temp.vals[i]+=self.vals[i,j]*a.vals[j]
                return temp
            else:
                sys.exit('error: \'*\' operator must be applied to tensor/vector with the same size !!!')
        else:
            if isinstance(a,int) or isinstance(a,float):
                # for scalar
                temp=rank2tensor(self.dim,0.0)
                for i in range(self.dim):
                    for j in range(self.dim):
                        temp.vals[i,j]=self.vals[i,j]*a
                return temp
            else:
                sys.exit('error: the right hand side value should be either a tensor or a scalar for \'*\' operator')
    # for '/'
    def __truediv__(self,a):
        """
        / operator overload

        Parameters
        ----------
        a : scalar
            right hand side value
        """
        if isinstance(a,rank2tensor):
            sys.exit('error: \'/\' operator can only be applied to a scalar !!!')
        else:
            if isinstance(a,int) or isinstance(a,float):
                # for scalar
                if np.abs(a)<1.0e-15:
                    sys.exit('error: \'/\' operator can only be applied to a non-singular scalar !!!')
                temp=rank2tensor(self.dim,0.0)
                for i in range(self.dim):
                    for j in range(self.dim):
                        temp.vals[i,j]=self.vals[i,j]/a
                return temp
            else:
                sys.exit('error: the right hand side value should be either a tensor or a scalar for \'/\' operator')

