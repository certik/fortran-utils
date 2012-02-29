from numpy import array
from h5py import File

f = File("dsetf2.h5")
b = array(f["/g/f/test5"])
c = array(f["/test6"])
print b
print c
