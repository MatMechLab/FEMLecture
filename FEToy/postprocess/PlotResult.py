__author__="Yang Bai"
__copyright__= "Copyright (C) 2021-present by M3 Group"
__version__ = "1.0"
__maintainer__ = "Yang Bai"
__email__ = "yangbai90@outlook.com"
__status__ = "development"
__date__ = "Jan 23, 2022"


import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
#from matplotlib.mlab import griddata

class Plot2D:
    def Contour2D(mesh,sol,savefig=False,figname=''):
        fig, ax = plt.subplots(constrained_layout=True)
        x,y=mesh.nodecoords[:,0],mesh.nodecoords[:,1]
        X,Y=np.meshgrid(x,y)
        Z=griddata((x,y),sol,(X,Y),method='linear')
        #cs=plt.contourf(X,Y,Z,cmap=plt.cm.hsv,levels=200)
        #cs=plt.contourf(X,Y,Z,cmap=plt.cm.viridis,levels=200,antialiased=True,extend='both')
        cs=plt.tricontourf(x,y,sol,levels=20)
        fig.colorbar(cs)
        ax.axis('equal')

        if savefig:
            if len(figname)>4:
                fig.savefig(figname,dpi=300,bbox_inches='tight')
                print('save results to ',figname)
            else:
                fig.savefig('result.jpg',dpi=300,bbox_inches='tight')
                print('save result to result.jpg')
    def Contour2DDeform(mesh,disp,sol,savefig=False,figname=''):
        fig, ax = plt.subplots(constrained_layout=True)
        x=np.zeros(mesh.nodes)
        y=np.zeros(mesh.nodes)
        for i in range(mesh.nodes):
            x[i]=mesh.nodecoords[i,0]+disp[2*i]
            y[i]=mesh.nodecoords[i,1]+disp[2*i+1]
        #cs=plt.contourf(X,Y,Z,cmap=plt.cm.hsv,levels=200)
        #cs=plt.contourf(X,Y,Z,cmap=plt.cm.viridis,levels=200,antialiased=True,extend='both')
        X,Y=np.meshgrid(x,y)
        Z=griddata((x,y),sol,(X,Y),method='linear')
        #cs=plt.contourf(X,Y,Z,cmap=plt.cm.viridis,levels=200,antialiased=True,extend='both')
        cs=plt.tricontourf(x,y,sol,levels=20)
        #cs=plt.contourf(x,y,sol,200)
        fig.colorbar(cs)
        ax.axis('equal')

        if savefig:
            if len(figname)>4:
                fig.savefig(figname,dpi=300,bbox_inches='tight')
                print('save results to ',figname)
            else:
                fig.savefig('result.jpg',dpi=300,bbox_inches='tight')
                print('save result to result.jpg')
