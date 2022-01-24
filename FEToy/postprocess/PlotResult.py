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

class Plot2D:
    def Contour2D(x,y,sol):
        fig, ax = plt.subplots(constrained_layout=True)
        X,Y= np.meshgrid(x,y)
        Z = griddata((x,y), sol, (X, Y),method='nearest')
        #cs=plt.contourf(X,Y,Z,cmap=plt.cm.hsv,levels=200)
        cs=plt.contourf(X,Y,Z,cmap=plt.cm.viridis,levels=400)
        fig.colorbar(cs)
