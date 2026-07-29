import os
import numpy as np
from matplotlib import pyplot as plt
from pathlib import Path 
# functions
def getTurbFiles():
    # get the file dir first!
    global data_dir
    path_dir=Path(__file__).resolve().parent
    run_dir=path_dir.parent/ "run"
    data_dir=run_dir/"data"

    # now lets look for the files
    f_whole=os.listdir(data_dir)
    if not f_whole:
        raise Exception("there is no file")
    f_turb=[]
    for fname in f_whole:
        if "turb" in fname:
            f_turb.append(fname)
    if not f_turb:
        raise Exception("there is no turbulence statistics file")
    f_turb.sort()
    return f_turb

def up_zp(fturb,hmean,re_tau):
    data=data_load(fturb)
    fig,ax = plt.subplots()
    nu=1/re_tau
    condition=(data[:,0]>=hmean) & (data[:,0]<=1.0+hmean)
    data=data[condition]
    zp=(data[:,0]-hmean)/nu
    up=data[:,1]
    ax.semilogx(zp,up,label="u+ over z+")
    #ax.plot(zp,up,label="u+ over z+")
    ax.set_xlabel("z+")
    ax.set_ylabel("u+")
    ax.set_title("u+ over z+")
    ax.legend()
    plt.show()
    return

def meanVelocity(fturb,hmean,re_tau):
    data=data_load(fturb)
    fig,ax = plt.subplots()
    zp=data[:,0]
    um=data[:,1]
    ax.plot(um,zp,label="Mean X Velocity")
    #ax.plot(data[:,2],data[:,0],label="Mean Y Velocity")
    #ax.plot(data[:,3],data[:,0],label="Mean Z Velocity")
    ax.set_xlabel("Um")
    ax.set_ylabel("z")
    ax.set_title("Mean Velocity")
    ax.legend()
    plt.show()
    return  

def statiscialStationary(fturb):
    # CaNS already provides us with Bulk Velocity
    # information and time information under "forcing.out" file
    # for checking the statiscal stationary we need to sample each 
    # ETT (ETT=h/u_tau) w/ our non-dimensionalization its ~1
    data=data_load(fturb)
    fig,ax=plt.subplots()
    # bulk velocity
    ub=data[:,4]
    zg=data[:,0]
    ax.plot(zg,ub,label="Ub")
    #ax.plot(data[:,0],data[:,5],label="Vb")
    #ax.plot(data[:,0],data[:,6],label="Wb")
    # pressure grad
    #ax.plot(data[:,0],data[:,1],label="Pressure Gradient")
    ax.set_xlabel("ETT")
    ax.set_ylabel("Bulk Velocity- Pressure Gradient")
    ax.legend()
    ax.set_title("Statistical Stationary Map")
    fig.tight_layout()
    plt.show()
    return 

def getWallShearSt(fturb,hmean,re_tau):
    data=data_load(fturb)

    condition=(data[:,0]>=hmean) & (data[:,0]<=1+hmean) 
    data=data[condition]
    fig,ax=plt.subplots()
    nu=1.0/re_tau
    z_g=data[:,0]-hmean
    uw=-data[:,7]
    um=data[:,1]
    dudz=np.gradient(um,z_g,edge_order=2)
    sh_visc=nu*dudz
    sh_tot=sh_visc+uw
    sh_tot_an=1-z_g[:]
    ax.plot(z_g,uw,label="Reynolds Shear Stress")
    ax.plot(z_g,sh_visc,label="Viscous Shear Stress")
    ax.plot(z_g,sh_tot,label="Total Shear Stress")
    ax.plot(z_g,sh_tot_an,label="Analytical Total Shear Stress")

    ax.legend()
    ax.set_xlabel("z/h")
    ax.set_ylabel("Shear St")
    ax.set_title("Shear Stress Distribution")
    plt.show()
    return

def getReynoldsStComparisson(fturb,hmean,re_tau):
    data=data_load(fturb)
    condition=(data[:,0]>=hmean) & (data[:,0]<=1+hmean) 
    data=data[condition]
    fig,ax=plt.subplots()
    z_g=data[:,0]-hmean
    u2=data[:,4]
    v2=data[:,5]
    w2=data[:,6]
    uw=data[:,7]
    uv=data[:,8]
    vw=data[:,9]
    k=np.zeros_like(u2)
    k=0.5*(u2+v2+w2)
    # z+ = z*u_tau/nu
    # nu = 1/Re; u_tau=sqrt(tau_w/rho)=1  
    nu=1/re_tau
    zp=z_g[:]/nu
    ax.plot(zp,u2,label="u2")
    ax.plot(zp,v2,label="v2")
    ax.plot(zp,w2,label="w2")
    ax.plot(zp,uw,label="uw")
    ax.plot(zp,uv,label="uv")
    ax.plot(zp,vw,label="vw")
    ax.plot(zp,k,label="k")
    ax.legend()
    ax.set_xlabel("z+")
    ax.set_ylabel("Reynolds Stresses and T.K.E \n uiuj/u_tau")
    ax.set_title("Reynolds Stresses over z+")
    fig.tight_layout()
    plt.show()
    return

def plotHelper(p1,p2,title=None,xlabel=None,ylabel=None,\
               p_label=None,*args):
    fig,ax=plt.subplots()
    if not args:
        ax.plot(p1,p2,label=p_label)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.legend()
        fig.tight_layout
        plt.show()
    else:
        for arg in args:
            ax.plot(p1,arg,label=str(arg))
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)
            ax.set_title(str(arg)+" - "+str(p1))
            ax.legend()
            fig.tight_layout
            plt.show()
    return
def data_load(fturb):
    #either data location or data it self
    if isinstance(fturb,str):
        file_dir=f"{data_dir}/{fturb}"
        data=np.loadtxt(file_dir)
    else:
        data=fturb
    return data
def getTimeAve(sequence:list):
    #here we define the sequence
    if not sequence:
        raise Exception("No file sequence are given")
    n=np.shape(sequence)[0] #amount of files we have
    #lets define the vars
    sizing_data=data_load(sequence[0])
    time_aveData=np.zeros_like(sizing_data)
    um=0.0
    vm=0.0
    wm=0.0
    u2=0.0
    v2=0.0
    w2=0.0
    uw=0.0
    uv=0.0
    vw=0.0
    for file in sequence:
        data=data_load(file)
        zg=data[:,0]
        um=um+data[:,1]
        vm=vm+data[:,2] 
        wm=wm+data[:,3]
        u2=u2+data[:,4]
        v2=v2+data[:,5]
        w2=w2+data[:,6]
        uw=uw+data[:,7]
        uv=uv+data[:,8]
        vw=vw+data[:,9]
    time_aveData[:,0]=zg
    time_aveData[:,1]=um/n
    time_aveData[:,2]=vm/n
    time_aveData[:,3]=wm/n
    time_aveData[:,4]=u2/n
    time_aveData[:,5]=v2/n
    time_aveData[:,6]=w2/n
    time_aveData[:,7]=uw/n
    time_aveData[:,8]=uv/n
    time_aveData[:,9]=vw/n
    meanVelocity(time_aveData,0.162,180)
    up_zp(time_aveData,0.162,180)
    getWallShearSt(time_aveData,0.162,180)
    getReynoldsStComparisson(time_aveData,0.162,180)
    return time_aveData
#
#
fturb=getTurbFiles()
print(fturb)
meanVelocity(fturb[-1],0.162,180)
up_zp(fturb[-1],0.162,180)
statiscialStationary("forcing.out")
getWallShearSt(fturb[-1],0.162,180)
getReynoldsStComparisson(fturb[-1],0.162,180)
time_aveData = getTimeAve(fturb[0:3])