__author__="Yang Bai"
__copyright__= "Copyright (C) 2022-present by M3 Group"
__version__ = "1.0"
__maintainer__ = "Yang Bai"
__email__ = "yangbai90@outlook.com"
__status__ = "development"
__date__ = "June 30, 2022"

import numpy as np
import matplotlib.pyplot as plt
import sys


class vector:
    def __init__(self,dim=2,value=0.0):
        """
        Initialize the vector

        Parameters
        ----------
        dim : int
            the dimension or size of current vector
        value : double
            default value for current vector
        """
        self.dim=dim
        if dim>3:
            sys.exit('error: dim=%d is invalid for a vector'%(dim))
        self.vals=np.zeros(self.dim)
        self.vals[:]=1.0*value
    def setToRandom(self):
        """
        Set current vector to random value

        Parameters
        ----------
        """
        for i in range(self.dim):
            self.vals[i]=np.random.rand()
    def print(self,str=None):
        """
        Print the vector

        Parameters
        ----------
        str : string
            the content to be printed
        """
        if self.dim==2:
            if str is None:
                print('vector components are: %14.5e, %14.5e'%(self.vals[0],self.vals[1]))
            else:
                print(str+', vector components are: %14.5e, %14.5e'%(self.vals[0],self.vals[1]))
        else:
            if str is None:
                print('vector components are: %14.5e, %14.5e %14.5e'%(self.vals[0],self.vals[1],self.vals[2]))
            else:
                print(str+',vector components are: %14.5e, %14.5e %14.5e'%(self.vals[0],self.vals[1],self.vals[2]))
    #################################################
    ### for operators
    #################################################
    def __getitem__(self,i):
        """
        Get the i-th element of current vector

        Parameters
        ----------
        i : int
            the integer index
        """
        if isinstance(i,int):
            if i<0 or i>self.dim-1:
                sys.exit('i=%d is out of range for a vector element access !!!'%(i))
            return self.vals[i]
        elif isinstance(i,slice):
            if i.start is not None and i.stop is not None:
                m = max(i.start, i.stop)
                return [self[ii] for ii in xrange(*i.indices(m+1))]
            else:
                return self.vals[:]
            # if np.min(i)<1 or np.max(i)>self.dim:
            #     sys.exit(i+' is out of range for a vector element access !!!')
            # return self.vals[i-1]
    def __setitem__(self,i,val):
        """
        Set the i-th element of current vector

        Parameters
        ----------
        i : int
            the integer index
        val : float
            the float scalar 
        """
        if isinstance(i,int):
            if i<0 or i>self.dim-1:
                sys.exit('i=%d is out of range for a vector element access !!!'%(i))
            self.vals[i]=val
        elif isinstance(i,slice):
            if i.start is not None and i.stop is not None:
                m = max(i.start, i.stop)
                for ii in xrange(*i.indices(m+1)):
                    self.vals[ii]=val
            else:
                self.vals[:]=val
    def __str__(self):
        """
        Print() method overload
        ----------
        """
        if self.dim==2:
            return "(%14.5e, %14.5e)"%(self.vals[0],self.vals[1])
        else:
            return "(%14.5e, %14.5e, %14.5e)"%(self.vals[0],self.vals[1],self.vals[2])
    # for '+'
    def __add__(self,a):
        """
        + operator overload

        Parameters
        ----------
        a : vector or scalar
            right hand side value
        """
        if isinstance(a,vector):
            if self.dim==a.dim:
                temp=vector(self.dim,0.0)
                temp.vals[:]=self.vals[:]+a.vals[:]
                return temp
            else:
                sys.exit('error: \'+\' operator must be applied to two vector with the same size !!!')
        else:
            if isinstance(a,int) or isinstance(a,float):
                # for scalar
                temp=vector(self.dim,0.0)
                temp.vals[:]=self.vals[:]+a
                return temp
            else:
                sys.exit('error: the right hand side value should be either a vector or a scalar for \'+\' operator')
    # for '-'
    def __sub__(self,a):
        """
        - operator overload

        Parameters
        ----------
        a : vector or scalar
            right hand side value
        """
        if isinstance(a,vector):
            if self.dim==a.dim:
                temp=vector(self.dim,0.0)
                temp.vals[:]=self.vals[:]-a.vals[:]
                return temp
            else:
                sys.exit('error: \'-\' operator must be applied to two vector with the same size !!!')
        else:
            if isinstance(a,int) or isinstance(a,float):
                # for scalar
                temp=vector(self.dim,0.0)
                temp.vals[:]=self.vals[:]-a
                return temp
            else:
                sys.exit('error: the right hand side value should be either a vector or a scalar for \'-\' operator')
    # for '*'
    def __mul__(self,a):
        """
        * operator overload

        Parameters
        ----------
        a : vector or scalar
            right hand side value
        """
        if isinstance(a,vector):
            if self.dim==a.dim:
                if self.dim==2:
                    sum=self.vals[0]*a.vals[0]\
                       +self.vals[1]*a.vals[1]
                else:
                    sum=self.vals[0]*a.vals[0]\
                       +self.vals[1]*a.vals[1]\
                       +self.vals[2]*a.vals[2]
                return sum
            else:
                sys.exit('error: \'*\' operator must be applied to two vector with the same size !!!')
        else:
            if isinstance(a,int) or isinstance(a,float):
                # for scalar
                temp=vector(self.dim,0.0)
                temp.vals[:]=self.vals[:]*a
                return temp
            else:
                sys.exit('error: the right hand side value should be either a vector or a scalar for \'*\' operator')
    # for '/'
    def __truediv__(self,a):
        """
        / operator overload

        Parameters
        ----------
        a : scalar
            right hand side value
        """
        if isinstance(a,vector):
            sys.exit('error: \'/\' operator can only be applied to a scalar !!!')
        else:
            if isinstance(a,int) or isinstance(a,float):
                # for scalar
                if np.abs(a)<1.0e-15:
                    sys.exit('error: \'/\' operator can only be applied to a non-singular scalar !!!')
                temp=vector(self.dim,0.0)
                temp.vals[:]=self.vals[:]/a
                return temp
            else:
                sys.exit('error: the right hand side value should be either a vector or a scalar for \'/\' operator')

