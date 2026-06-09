import numpy as np 
import os
import subprocess
import shutil  

def reset_case():
    subprocess.run(["bash","-c","rm -rf /home/cumhur/codes/CaNS-IBM/run/data/*"])
    shutil.copy2("write_xdmf.py","/home/cumhur/codes/CaNS-IBM/run/data/")
    subprocess.run(["bash","-c","rm -f /home/cumhur/codes/CaNS-IBM/run/input.nml"])
    shutil.copy2("input.nml","/home/cumhur/codes/CaNS-IBM/run/") 

reset_case()