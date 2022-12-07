import numpy as np
import sys
import matplotlib.pyplot as plt

def plotImage(iMap):
    plt.imshow(iMap, cmap='jet')
    plt.show()
def plotxVelbySect(u):
    fig1, ax1 = plt.subplots()
    nx= len(u[0,:])
    ny=len(u[:,0])
    uo=u[ny//2,0]
    #uo=0.2
    Tm1 = np.zeros(ny)
    Tm2 = np.zeros(ny)
    Tm3 = np.zeros(ny)
    Tm4 = np.zeros(ny)
    Tm5 = np.zeros(ny)
    sec1 = nx//6
    sec2 = nx//6*2
    sec3 = nx//6*3
    sec4 = nx//6*4
    sec5 = nx//6*5
    Tm1=u[1:-1,sec1]/uo
    Tm2=u[1:-1,sec2]/uo
    Tm3=u[1:-1,sec3]/uo
    Tm4=u[1:-1,sec4]/uo
    Tm5=u[1:-1,sec5]/uo
    y = np.arange(1,ny-1,1)
    ax1.plot(y,Tm1,linewidth=1.5,label='Section {}'.format(sec1))
    ax1.plot(y,Tm2,linewidth=1.5,label='Section {}'.format(sec2))
    ax1.plot(y,Tm3,linewidth=1.5,label='Section {}'.format(sec3))
    ax1.plot(y,Tm4,linewidth=1.5,label='Section {}'.format(sec4))
    ax1.plot(y,Tm5,linewidth=1.5,label='Section {}'.format(sec5))
    ax1.set_xlabel('Y')
    ax1.set_ylabel('U')
    ax1.set_title('X component of velocity along channel')
    ax1.legend()
    plt.show()

if __name__ == '__main__':
    fname = sys.argv[1]
    u = np.genfromtxt(fname, delimiter='\t')
    nrow = u.shape[0]
    ncol= u.shape[1]
    if (nrow > ncol):
        u = np.transpose(u)
        u = -u
    plotImage(u)
    plotxVelbySect(u)
