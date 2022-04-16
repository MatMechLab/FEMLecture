__author__="Yang Bai"
__copyright__= "Copyright (C) 2021-present by M3 Group"
__version__ = "1.0"
__maintainer__ = "Yang Bai"
__email__ = "yangbai90@outlook.com"
__status__ = "development"
__date__ = "April 15, 2022"


import numpy as np
import sys

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
        if not mesh.nodes*len(varnamelist)==len(solution):
            sys.exit('your varnamelist length*nodes does not match with your solution!')
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

    def save2vtu(self,mesh,solution,varnamelist,step):
        """
        Save results to vtu file

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
        if not mesh.nodes*len(varnamelist)==len(solution):
            sys.exit('your varnamelist length*nodes dose not match with your solution!')
        filename='%06d.vtu'%(step)
        self.filename=self.prefixname+'-'+filename
        inp=open(self.filename,'w+')
        vtkcelltype=mesh.vtkcelltype

        inp.write("<?xml version=\"1.0\"?>\n")
        inp.write("<VTKFile type=\"UnstructuredGrid\" version=\"1.0\">\n")
        inp.write("<UnstructuredGrid>\n")

        str="<Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n"%(mesh.nodes,mesh.elements)
        inp.write(str)
        inp.write("<Points>\n")
        inp.write("<DataArray type=\"Float64\" Name=\"nodes\"  NumberOfComponents=\"3\"  format=\"ascii\">\n")

        # write out nodal coordinates
        for i in range(mesh.nodes):
            x=0.0;y=0.0
            if mesh.dim==1:
                x=mesh.nodecoords[i]
            elif mesh.dim==2:
                x=mesh.nodecoords[i,0]
                y=mesh.nodecoords[i,1]
            str='%14.6e %14.6e 0.0\n'%(x,y)
            inp.write(str)
        inp.write("</DataArray>\n")
        inp.write("</Points>\n")

        # write out cell info
        inp.write("<Cells>\n")
        inp.write("<DataArray type=\"Int32\" Name=\"connectivity\" NumberOfComponents=\"1\" format=\"ascii\">\n")
        for e in range(mesh.elements):
            str=''
            for i in range(mesh.nodesperelement):
                str+='%6d '%(mesh.elementconn[e,i])
            str+='\n'
            inp.write(str)
        inp.write("</DataArray>\n")

        inp.write("<DataArray type=\"Int32\" Name=\"offsets\" NumberOfComponents=\"1\" format=\"ascii\">\n")
        offset=0
        for e in range(mesh.elements):
            offset+=mesh.nodesperelement
            str='%6d\n'%(offset)
            inp.write(str)
        inp.write("</DataArray>\n")

        # for cell type
        inp.write("<DataArray type=\"Int32\" Name=\"types\"  NumberOfComponents=\"1\"  format=\"ascii\">\n")
        for e in range(mesh.elements):
            str='%6d\n'%(vtkcelltype)
            inp.write(str)
        inp.write("</DataArray>\n")
        inp.write("</Cells>\n")

        str="<PointData Scalar=\""
        for j in range(len(varnamelist)):
            dofname=varnamelist[j]
            str+=dofname+' '
        str+="\" >\n"
        inp.write(str)

        for j in range(len(varnamelist)):
            dofname=varnamelist[j]
            str="<DataArray type=\"Float64\" Name=\"" +dofname+ "\"  NumberOfComponents=\"1\" format=\"ascii\">\n"
            inp.write(str)
            for i in range(mesh.nodes):
                iInd=i*len(varnamelist)+j
                str='%14.6e\n'%(solution[iInd])
                inp.write(str)
            inp.write("</DataArray>\n\n")
        inp.write("</PointData>\n")

        inp.write("</Piece>\n")
        inp.write("</UnstructuredGrid>\n")
        inp.write("</VTKFile>")

        inp.close()

        print('write result to %s'%(self.filename))