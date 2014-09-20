import matplotlib.pyplot as plt
from matplotlib import cm
import numpy
import pickle

fig, ax = plt.subplots(1)
f = open('Desktop/working/hamiltonian.pkl','r')
h1 = numpy.array(pickle.load(f))
h1 = h1[:h1.shape[0]/2 , :h1.shape[1]/2]
f.close()
f = open('Desktop/working/Ham_flow_M1.pkl','r')
h2 = numpy.array(pickle.load(f))
f.close()
f = open('Desktop/working/Ham_reduce_M1.pkl','r')
h3 = numpy.array(pickle.load(f))
f.close()


normalizeObj1 = plt.Normalize(vmin = -0.1, vmax = 0.1, clip = True)
normalizeObj2 = plt.matplotlib.colors.SymLogNorm(linthresh=0.1, linscale=0.1, vmin=-8.0, vmax=8.0, clip=True)


plt.imshow(h2, interpolation='nearest' , norm=normalizeObj2 , cmap=cm.bwr)
plt.colorbar()
plt.show()
print h3.shape