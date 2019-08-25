import numpy as np
import matplotlib.pyplot as plt
def scatterplot_city(c,N,numblocks,positions,colorlist=None,sizes=30,filename=None):
    colors=c/np.max(c,0,keepdims=True)
    block_colors=np.zeros([numblocks,N,4])
    if colorlist is None:
        colorlist=np.array([[1,0,0],[0,1,0],[0,0,1],\
                            [1,0,1],[1,1,0],[0,1,1],[0,0,0]])
    print(numblocks)
    fig=plt.figure(figsize=(15,15), dpi= 80, facecolor='w', edgecolor='k')
    for i in range(0,numblocks):
        block_colors[i,:,0:3]=colorlist[i%7,:]
        block_colors[i,:,3]=colors[:,i]
    for i in range(0,numblocks):
        plt.scatter(positions[:,1]
                    , positions[:,0]\
                    ,c=block_colors[i,:,:],s=sizes,edgecolors='k')
    if not(filename is None):
        plt.savefig(filename+"scatter.png")

def block_omegas(r,numblocks,omega_max=None,filename=None):
    if omega_max is None:
        omega_max=np.max(r)
    f,axes=plt.subplots(nrows=numblocks,ncols=numblocks,figsize=(20,10))
    
    colorlist=np.array([[1,0,0],[0,1,0],[0,0,1],\
                            [1,0,1],[1,1,0],[0,1,1],[0,0,0]])
    for block1 in range(0,numblocks):
        for block2 in range(0,numblocks):
            axes[block1,block2].plot(range(0,24),r[block1,block2,:])
            axes[block1,block2].set_ylim(0,omega_max)
            axes[block1,block2].set_xlim(0,23)
    #f.tight_layout()
    for ax, col in zip(axes[0], range(1,numblocks+1)):
        ax.set_title("Block-"+str(col),color=colorlist[(col-1)%7])
    
    for ax, row in zip(axes[:,0], range(1,numblocks+1)):
        ax.set_ylabel("Block-"+str(row), rotation=0, size='large',color=colorlist\
                      [(row-1)%7])
    
    if not(filename is None):
        plt.savefig(filename+”block_omegas.png")
def save_model(c,r,N,numiter,ll,numblocks,positions,cityname):
    numblocks=int(numblocks)
    f1=open(cityname+'Blocks', 'w+')
    f1.write(",lat,lon,")
    for k in range(0,numblocks):
        f1.write("block"+str(k))
        if k<numblocks-1:
                f1.write(",")
    f1.write("\n")
    for i in range(0,N):
        f1.write(str(i)+","+str(positions[i][0])+","+str(positions[i][1])+",")
        for k in range(0,numblocks):
            f1.write(str(c[i][k]))
            if k<numblocks-1:
                f1.write(",")
        f1.write("\n")
            
    f1.close()
    f2=open(cityname+'Blocks', 'w+')

    for i in range(0,numblocks):
        for j in range(0,numblocks):
            f2.write(str(i)+" "+str(j)+"-time profile:")
            if (i!=numblocks-1)|(j!=numblocks-1):
                    f2.write(",")
    f2.write("\n")
    for t in range(0,24):
        for i in range(0,numblocks):
            for j in range(0,numblocks):        
                f2.write(str(r[i][j][t]))
                if (i!=numblocks-1)|(j!=numblocks-1):
                    f2.write(",")
        f2.write("\n")
    f2.close()
    f3=open(cityname+'loglikelihood','w+')
    f3.write('number of iterations,log-likelihood \n')
    f3.write(str(numiter)+","+str(ll))
    f3.close()
