#!/usr/bin/env python3
import sys 
from subprocess import Popen, PIPE
from snakemake.utils import read_job_properties
from time import sleep

jobscript = sys.argv[-1]  
dependencies = set(sys.argv[1:-1])

jobprops = read_job_properties(jobscript)

cmdline = ["sbatch"]

cmdline.append("--parsable")

if dependencies:
    cmdline.append("--dependency")
    cmdline.append( "afterok:" + ",".join(dependencies))

if jobprops["cluster"]:
    if jobprops["cluster"]["time"]:
        cmdline.append("--time=" + jobprops["cluster"]["time"])
    if jobprops["cluster"]["partition"]:
        cmdline.append("--partition=" + jobprops["cluster"]["partition"])
    if jobprops["cluster"]["cpus"]:
        cmdline.append("--cpus-per-task=" + jobprops["cluster"]["cpus"])
    if jobprops["cluster"]["mem"]:
        cmdline.append("--mem=" + jobprops["cluster"]["mem"])
    if jobprops["cluster"]["jobname"]:
        cmdline.append("--job-name=" + jobprops["cluster"]["jobname"])
    if jobprops["cluster"]["output"]:
        cmdline.append("--output=" + jobprops["cluster"]["output"])
    if jobprops["cluster"]["error"]:
        cmdline.append("--error=" + jobprops["cluster"]["error"])
   

cmdline.append(jobscript)

# Constructs and submits
#cmdline = " ".join(cmdline)
#out = open("1.txt", "a")
#out.write(" ".join(cmdline) + "\n")
#out.close()

sleepTime = 0.25
timeMultiplier = 4
for i in range(7):
    submitOut = Popen(cmdline, stdout = PIPE, stderr = PIPE)
    jobid, joberr = submitOut.communicate()
    
    jobid = jobid.decode("utf-8")
    
    if submitOut.returncode == 0:
        print(jobid, end = "")
        break
    
    sleep(sleepTime)
    sleepTime *= timeMultiplier
