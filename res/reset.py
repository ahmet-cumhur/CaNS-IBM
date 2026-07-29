import numpy as np 
import os
import subprocess
from pathlib import Path
def reset_case():
    file_dir=Path(__file__).resolve().parent
    main_dir=file_dir.parent
    run_dir=main_dir / "run"
    data_dir=run_dir / "data" 

    subprocess.run(["bash","-c",f"rm -rf {data_dir}/*"])
    subprocess.run(["bash","-c",f"cp {file_dir}/write_xdmf.py {data_dir}"])
    subprocess.run(["bash","-c",f"rm -f {run_dir}/input.nml"])
    subprocess.run(["bash","-c",f"cp {file_dir}/input.nml {run_dir}"])

reset_case()