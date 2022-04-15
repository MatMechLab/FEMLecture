__author__="Yang Bai"
__copyright__= "Copyright (C) 2021-present by M3 Group"
__version__ = "1.0"
__maintainer__ = "Yang Bai"
__email__ = "yangbai90@outlook.com"
__status__ = "development"
__date__ = "April 15, 2022"


import numpy as np

class ResultIO:
    def __init__(self,prefixname=''):
        """
        Initialize the ResultIO class 

        Parameters
        ----------
        prefixname : string
            the name of the output file(only prefix)
        """
        self.prefixname=prefixname
        self.filename=''
        if len(prefixname)<1:
            self.prefixname='myresult'
    def save2csv(self,mesh,solution,varnamelist,step):
        """
        Save results to csv file

        Parameters
        ----------
        mesh : mesh class
            the mesh class
        solution : numpy array
            the solution of each node
        varnamelist : list
            the name list of your dofs
        step : int
            the current time step
        """
        filename='%06d.csv'%(step)
        self.filename=self.prefixname+'-'+filename
        inp=open(self.filename,'w+')
        str=''
        if mesh.dim==1:
            str='x'
        elif mesh.dim==2:
            str='x,y'
        for i in varnamelist:
                str+=','+i
        str+='\n'
        inp.write(str)
        nodedofs=len(varnamelist)
        for i in range(mesh.nodes):
            if mesh.dim==1:
                x=mesh.nodecoords[i]
                str='%14.5e'%(x)
            elif mesh.dim==2:
                x=mesh.nodecoords[i,0]
                y=mesh.nodecoords[i,1]
                str='%14.5e,%14.5e'%(x,y)
            for j in range(nodedofs):
                iInd=i*nodedofs+j
                str+=',%14.5e'%(solution[iInd])
            str+='\n'
            inp.write(str)

        inp.close()
        print('write result to %s'%(self.filename))

