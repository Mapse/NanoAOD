import os, subprocess, sys, time

ndays = 5
nhours = ndays * 24
nminutes = nhours * 60
nsecs = nminutes * 60

def resub():
    for project in subprocess.check_output("ls CharmoniumRun2017UL/", shell=True).decode("utf-8").splitlines():
        output = subprocess.check_output("crab status -d CharmoniumRun2017UL/" + project, shell=True).decode("utf-8")
        print(output)
        if output.find("failed") >-1:
            os.system("crab resubmit -d CharmoniumRun2017UL/" + project)

for i in range(0, nhours):
    resub()
    print("Processing...")
    time.sleep(1800)
    print("Processing...")
