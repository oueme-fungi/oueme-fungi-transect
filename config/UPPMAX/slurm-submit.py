#!/usr/bin/env python3
import sys
import os
import re
import argparse
import subprocess

from snakemake.utils import read_job_properties

parser = argparse.ArgumentParser(add_help=False)
parser.add_argument(
    "--help", help="Display help message.", action="store_true")
parser.add_argument(
    "positional", action="append",
    nargs="?", metavar="POS",
    help="additional arguments not in slurm parser group to pass to sbatch")

# A subset of SLURM-specific arguments
slurm_parser = parser.add_argument_group("slurm-specific arguments")
slurm_parser.add_argument(
    "-a", "--array", help="job array index values")
slurm_parser.add_argument(
    "-A", "--account", help="charge job to specified account")
slurm_parser.add_argument(
    "--begin", help="defer job until HH:MM MM/DD/YY")
slurm_parser.add_argument(
    "-c", "--cpus-per-task", help="number of cpus required per task")
slurm_parser.add_argument(
    "-M", "--clusters", help="which cluster to run the job on")
slurm_parser.add_argument(
    "-d", "--dependency",
    help="defer job until condition on jobid is satisfied")
slurm_parser.add_argument(
    "-D", "--workdir", help="set working directory for batch script")
slurm_parser.add_argument(
    "-e", "--error", help="file for batch script's standard error",
    default="logs/snakemake-%j.log" if "logs/snakemake-%j.log" else None)
slurm_parser.add_argument(
    "-J", "--job-name", help="name of job")
slurm_parser.add_argument(
    "--mail-type", help="notify on state change: BEGIN, END, FAIL or ALL")
slurm_parser.add_argument(
    "--mail-user", help="who to send email notification for job state changes")
slurm_parser.add_argument(
    "-n", "--ntasks", help="number of tasks to run")
slurm_parser.add_argument(
    "-N", "--nodes", help="number of nodes on which to run (N = min[-max])")
slurm_parser.add_argument(
    "-o", "--output", help="file for batch script's standard output",
    default="logs/snakemake-%j.log" if "logs/snakemake-%j.log" else None)
slurm_parser.add_argument(
    "-p", "--partition", help="partition requested")
slurm_parser.add_argument(
    "-Q", "--quiet", help="quiet mode (suppress informational messages)")
slurm_parser.add_argument(
    "-t", "--time", help="time limit")
slurm_parser.add_argument(
    "--wrap", help="wrap command string in a sh script and submit")
slurm_parser.add_argument(
    "-C", "--constraint", help="specify a list of constraints")
slurm_parser.add_argument(
    "--mem", help="minimum amount of real memory")

args = parser.parse_args()

if args.help:
    parser.print_help()
    sys.exit(0)

jobscript = sys.argv[-1]
job_properties = read_job_properties(jobscript)

extras = ""
if args.positional:
    for m in args.positional:
        if m is not None:
            extras = extras + " " + m

arg_dict = dict(args.__dict__)


# Process resources
if "resources" in job_properties:
    resources = job_properties["resources"]
    if arg_dict["time"] is None:
        if "runtime" in resources:
            arg_dict["time"] = resources["runtime"]
        elif "walltime" in resources:
            arg_dict["time"] = resources["walltime"]
    if arg_dict["mem"] is None:
        if "mem" in resources:
            arg_dict["mem"] = resources["mem"]
        elif "mem_mb" in resources:
            arg_dict["mem"] = resources["mem_mb"]


# Threads
if "threads" in job_properties:
    arg_dict["cpus_per_task"] = job_properties["threads"]

# Dependencies
if arg_dict["dependency"] == "afterok:":
    arg_dict.pop("dependency")
else:
    arg_dict["dependency"] = re.sub("\s+", ":", arg_dict["dependency"])    

opt_keys = ["array", "account", "begin", "cpus_per_task",
            "dependency", "workdir", "error", "job_name", "mail_type",
            "mail_user", "ntasks", "nodes", "output", "partition",
            "quiet", "time", "wrap", "constraint", "mem", "clusters"]

# Set default partition
if arg_dict["partition"] is None:
    if not "core":
        # partitions and SLURM - If not specified, the default behavior is to
        # allow the slurm controller to select the default partition as
        # designated by the system administrator.
        opt_keys.remove("partition")
    else:
        arg_dict["partition"] = "core"

# Set default account
if arg_dict["account"] is None:
    if "snic2018-8-131" != "":
        arg_dict["account"] = "snic2018-8-131"

# Ensure output folder for Slurm log files exist.
# This is a bit hacky; will run for every Slurm submission...
if arg_dict["output"] is not None:
    os.makedirs(os.path.dirname(arg_dict["output"]), exist_ok=True)
if arg_dict["error"] is not None:
    os.makedirs(os.path.dirname(arg_dict["error"]), exist_ok=True)

opts = ""
for k, v in arg_dict.items():
    if k not in opt_keys:
        continue
    if v is not None:
        opts += " --{} \"{}\" ".format(k.replace("_", "-"), v)

if arg_dict["wrap"] is not None:
    cmd = "sbatch {opts}".format(opts=opts)
else:
    cmd = "sbatch {opts} {extras}".format(opts=opts, extras=extras)

try:
    res = subprocess.run(cmd, check=True, shell=True, stdout=subprocess.PIPE)
except subprocess.CalledProcessError as e:
    raise e

# Get jobid
res = res.stdout.decode()
try:
    m = re.search("Submitted batch job (\d+)", res)
    jobid = m.group(1)
    print(jobid)
except Exception as e:
    print(e)
    raise
