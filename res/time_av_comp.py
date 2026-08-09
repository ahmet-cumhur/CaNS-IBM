import os
import numpy as np
from matplotlib import pyplot as plt
from pathlib import Path 
# functions
def getTurbFiles(floc:str):
    # get the file dir first!
    path_dir=Path(__file__).resolve().parent
    data_dir=path_dir/f"{floc}"

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
    return f_turb,data_dir


def data_load(fturb,data_dir=None):
    #either data location or data it self
    if isinstance(fturb,str):
        file_dir=f"{data_dir}/{fturb}"
        data=np.loadtxt(file_dir)
    else:
        data=fturb
    return data
def getTimeAve(sequence:list,data_dir,hmean):
    #here we define the sequence
    if not sequence:
        raise Exception("No file sequence are given")
    n=np.shape(sequence)[0] #amount of files we have
    #lets define the vars
    sizing_data=data_load(sequence[0],data_dir)
    #that +1 is for Uc
    row,col=np.shape(sizing_data)
    time_aveData=np.zeros((row,col+1))
    dz=abs(sizing_data[1,0]-sizing_data[2,0])
    um=0.0
    vm=0.0
    wm=0.0
    u2=0.0
    v2=0.0
    w2=0.0
    uw=0.0
    uv=0.0
    vw=0.0
    uc=0.0
    for file in sequence:
        data=data_load(file,data_dir)
        condCenter=(data[:,0]>=+1+hmean-dz/2)&(data[:,0]<=+1+hmean+dz/2)
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
        data=data[condCenter]
        uc=uc+np.mean(data[:,1])
        
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
    time_aveData[:,10]=uc/n
    return time_aveData
def getBulkMean(st,fin,datadir):
    data = np.loadtxt(datadir)
    Ub_mean = np.mean(data[st:fin, 4])
    return Ub_mean
def up_zp(dataRough,dataSmooth,hmean,re_tau):
    dataR=data_load(dataRough)
    dataS=data_load(dataSmooth)
    fig,ax = plt.subplots()
    nu=1/re_tau
    conditionR=(dataR[:,0]>=hmean) & (dataR[:,0]<=1.0+hmean)
    conditionS=(dataS[:,0]>=0.0) & (dataS[:,0]<=1.0)
    dataR=dataR[conditionR]
    dataS=dataS[conditionS]
    zpR=(dataR[:,0]-hmean)/nu
    upR=dataR[:,1]
    zpS=(dataS[:,0])/nu
    upS=dataS[:,1]
    ax.semilogx(zpR,upR,label="Rough Data")
    ax.semilogx(zpS,upS,label="Smooth Data")
    ax.set_xlabel("z+")
    ax.set_ylabel("u+")
    ax.set_title("u+ over z+")
    ax.legend()
    fig.tight_layout()
    plt.show()
    fig.savefig("up.png")
    return

def uc_zd(dataRough,dataSmooth,hmean,re_tau):
    dataR=data_load(dataRough)
    dataS=data_load(dataSmooth)
    fig,ax = plt.subplots()
    nu=1/re_tau
    conditionR=(dataR[:,0]>=hmean) & (dataR[:,0]<=1.0+hmean)
    conditionS=(dataS[:,0]>=0.0) & (dataS[:,0]<=1.0)
    dataR=dataR[conditionR]
    dataS=dataS[conditionS]
    zpR=(dataR[:,0]-hmean)
    upR=dataR[0,10]-dataR[:,1]
    zpS=(dataS[:,0])
    upS=dataS[:,10]-dataS[:,1]
    ax.plot(zpR,upR,label="Rough Data")
    ax.plot(zpS,upS,label="Smooth Data")
    ax.set_xlabel("z/delta")
    ax.set_ylabel("Uc-u(z)/u_tau")
    ax.set_title("Uc-u(z)")
    ax.legend()
    fig.tight_layout()
    plt.show()
    fig.savefig("Ucuz.png")
    return
def wrTurbData(ubS,ubR,UcS,UcR,dUc,dub):
    header="ubS,ubR,UcS,UcR,dUc,dUb"
    data=np.array([[ubS,ubR,UcS,UcR,dUc,dub]])
    np.savetxt("caseResult.out",data,delimiter=" ",\
               fmt="%.18e",header=header)
    return 
dir_smooth="smooth"
dir_rough="rough"
fsmooth,ddirS=getTurbFiles(dir_smooth)
frough,ddirR=getTurbFiles(dir_rough)

fSmoothTA=getTimeAve(fsmooth[40:170],ddirS,0.0)
fRoughTA=getTimeAve(frough[40:170],ddirR,0.162)

up_zp(fRoughTA,fSmoothTA,0.162,180)
uc_zd(fRoughTA,fSmoothTA,0.162,180)

ubS=getBulkMean(40,170,f"./{dir_smooth}/forcing.out")
ubR=getBulkMean(40,170,f"./{dir_rough}/forcing.out")
UcS=fSmoothTA[0,-1]
UcR=fRoughTA[0,-1]
dUc=UcS-UcR
dub=ubS-ubR
#smooth Ub
print("Smooth Ub time averaged:",ubS)
#rough Ub
print("Rough Ub time averaged:",ubR)
#smooth Uc
print("Smooth Uc time averaged:",UcS)
# Rough Uc
print("Rough Ub time averaged:",UcR)
#DeltaUc+
print("delta Uc+:",dUc)
#DeltaUb+
print("delta Ub+:",dub)

wrTurbData(ubS,ubR,UcS,UcR,dUc,dub)
