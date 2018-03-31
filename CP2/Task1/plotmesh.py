import numpy as np
import matplotlib.pyplot as plt

#Read Velocity Data
x,y,u,v = np.loadtxt('velocity.plt',delimiter=',',dtype=np.float,skiprows=1,unpack=True)

#Reshape into a mesh
N = int(np.sqrt(len(x)))
x = np.reshape(x,(N,N))
y = np.reshape(y,(N,N))
u = np.reshape(u,(N,N))
v = np.reshape(v,(N,N))

#Plot data

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

#Read and plot phi
for iter in range(0,641,80):
    x,y,phi = np.loadtxt('phi'+str(iter)+'.plt',delimiter=',',dtype=np.float,skiprows=1,unpack=True)
    N = int(np.sqrt(len(x)))
    x = np.reshape(x,(N,N))
    y = np.reshape(y,(N,N))
    phi = np.reshape(phi,(N,N))
    
    plt.figure(facecolor='white')
    CS=plt.contour(x,y,phi,levels=[-0.25,0,0.25],colors='k')
    plt.clabel(CS,fontsize=10,inline=1)
    plt.title('Phi at Iteration '+str(iter))
    plt.xlabel('x')
    plt.ylabel('y')
    plt.xlim(0,8)
    plt.ylim(0,8)
    plt.axis('equal')
	plt.savefig('phi'+str(iter)+'.jpg')
