import numpy as np
import pylab as pl


f=open("reduced_data.csv")
data=[[float(el) for el in line.strip().split(",")] for line in f]
data=np.array(data)
f=open("means.csv")
means=[[float(el) for el in line.strip().split()] for line in f]
means=np.array(means)
pl.figure()
if(len(data[0])>1):
    pl.scatter(data[:,0],data[:,1], color="y")
    pl.scatter(means[:,0],means[:,1],color="g")
else:
    pl.scatter(data[:,0],np.zeros(len(data)), color="y")
    pl.scatter(means[:,0],np.zeros(len(means)),color="g")
pl.show()