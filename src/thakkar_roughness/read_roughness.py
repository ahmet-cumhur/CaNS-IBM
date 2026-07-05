# data from: https://www.tandfonline.com/doi/full/10.1080/14685248.2016.1258119
#before importing the data enter the data set number from thakkar paper!
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
#read data from thakkar roughness
#choose the data number 1-17
nmbr = 8
delta=1.0#normalization

path_data_t3="./data_2016JoT/Table3.csv"
path_data_t2="./data_2016JoT/Table2.csv"
path = f'./data_2016JoT/RoughnessHeightMaps/s{nmbr}.csv'
sample_name=f"s{nmbr}"
#load and edit the data of table 3 and table 2  
data_t3 = pd.read_csv(path_data_t3, skiprows=2)
data_t2 = pd.read_csv(path_data_t2, skiprows=2)

data_t3_im=data_t3.iloc[0],data_t3.iloc[nmbr]
data_t2_im=data_t2.iloc[:,1],data_t2.iloc[:,nmbr]

# load data
data = np.loadtxt(path,delimiter=',')
data=np.array(data)
ny_d,nx_d=np.shape(data)
nx=data_t3.iloc[nmbr-1,4]
ny=data_t3.iloc[nmbr-1,5]
nz=data_t3.iloc[nmbr-1,6]
Lx=np.astype(data_t3.iloc[nmbr-1,1],np.float64)
Ly=np.astype(data_t3.iloc[nmbr-1,2],np.float64)
Lz=np.astype(data_t3.iloc[nmbr-1,3],np.float64)
y = Ly*np.arange(data.shape[0])/data.shape[0]
x = Lx*np.arange(data.shape[1])/data.shape[1]


# roughness.bin (BINARY) -> the actual data:    x, y, height_map(y,x)    WARNING C/Python-Order of indices used here
# roughness.nfo (ASCII)  -> length(x), length(y), decription
with open("roughness.bin",'wb') as f:
    data.T.tofile(f)
    f.close()
with open("roughness.nfo","w") as f:
    f.write(f"&roughnessinfo\n")
    f.write(f"lx={Lx*delta}\n")
    f.write(f"ly={Ly*delta}\n")
    f.write(f"lz={Lz*delta}\n")
    f.write(f"nx={nx}\n")
    f.write(f"ny={ny}\n")
    f.write(f"nz={nz}\n")
    f.write(f"nx_data={nx_d}\n")
    f.write(f"ny_data={ny_d}\n")
    f.write("/\n")
    f.close()

