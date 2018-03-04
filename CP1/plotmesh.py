import numpy as np
import matplotlib.pyplot as plt

#Read Data
x,y,u,v,p = np.loadtxt('velocity.plt',delimiter=',',dtype=np.float,skiprows=1,unpack=True)

#Reshape data from vector to a matrix (mesh)
N = int(np.sqrt(len(x)))
x = np.reshape(x,(N,N))
y = np.reshape(y,(N,N))
u = np.reshape(u,(N,N))
v = np.reshape(v,(N,N))
p = np.reshape(p,(N,N))

#Plot data
#Plot u
plt.figure(facecolor='white')
plt.contourf(x,y,u,levels=np.arange(np.min(u),np.max(u),0.01),cmap=plt.cm.jet)
plt.colorbar()
plt.title('u')
plt.xlabel('x')
plt.ylabel('y')
plt.xlim(0,1)
plt.ylim(0,1)
plt.savefig('u.jpg')

#Plot v
plt.figure(facecolor='white')
plt.contourf(x,y,v,levels=np.arange(np.min(v),np.max(v),0.01),cmap=plt.cm.jet)
plt.colorbar()
plt.title('v')
plt.xlabel('x')
plt.ylabel('y')
plt.xlim(0,1)
plt.ylim(0,1)
plt.savefig('v.jpg')

#Plot p
plt.figure(facecolor='white')
plt.contourf(x,y,p,levels=np.arange(np.min(p),np.max(p),0.5),cmap=plt.cm.jet)
plt.colorbar()
plt.title('p')
plt.xlabel('x')
plt.ylabel('y')
plt.xlim(0,1)
plt.ylim(0,1)
plt.savefig('p.jpg')

#Plot streamlines
plt.figure(facecolor='white')
plt.streamplot(x,y,u,v)
plt.title('Streamlines')
plt.xlabel('x')
plt.ylabel('y')
plt.xlim(0,1)
plt.ylim(0,1)
plt.axis('equal')
plt.savefig('streamlines.jpg')

plt.show()
