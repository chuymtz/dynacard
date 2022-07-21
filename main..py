from scipy.io import loadmat
from matplotlib import pyplot as plt
mat = loadmat('Everitt/PUnit.mat')

type(mat)
mat.keys()
type(mat['Xi'])
mat['Xi'].shape

Xi = mat['Xi'][:,0]
plt.plot(Xi)

case1 = loadmat('Everitt/case1s.mat')
case1.keys()
 
u = case1['X'][:,0]
f = case1['F'][:,0]
 
plt.plot(u, f)