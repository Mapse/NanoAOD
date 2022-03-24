import os, subprocess, sys
import datetime

today = datetime.date.today()

njobs = 1 # Change for the wanted number of jobs
template = "crab_config_NANO_MC"

for path in subprocess.check_output("ls paths_monte_carlo/", shell=True).decode("utf-8").splitlines():
    if not (path.endswith(".txt")):
        continue
    dataset = path[0: path.find(".")]
    with open(template + ".py", 'r') as f:
        new_file = f.read().replace("DATASET", dataset)
        new_file = new_file.replace("DATE", today.strftime("%d-%m-%Y"))
        new_file = new_file.replace("FILE", path)
        new_file = new_file.replace("NJOBS", str(njobs))

        with open(template + "_" + dataset + ".py", 'w') as nf:
            nf.write(new_file)
    os.system("crab submit -c " + template + "_" + dataset + ".py")
