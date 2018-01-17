import argparse 
import json
import numpy as np 
import os
from os import path 
import subprocess as proc
import re

# Settings 

# - Current path for relative directories
cdir = path.dirname(path.abspath(__file__))

# - path to gmx executable (Gromacs 5.x)
gmx = path.join(cdir, 'gmx')

# - Number of walkers to use 
nwalkers = 4 

# - Task file containing prefixes for JSON files
taskfile = 'args.list'

# - Prefix for template files
prefix = 'template'

parser = argparse.ArgumentParser()

parser.add_argument(
    "-t",
    help="Template directory",
    type=str,
    default="template"
)

parser.add_argument(
    "-s",
    help="DAT files source directory",
    type=str,
    default="npt"
)

parser.add_argument(
    "-o",
    help="Destination directory", 
    type=str,
    default="preelastic"
)

parser.add_argument(
    "--nproc",
    help="Number of processors per walker",
    type=int,
    default=8
)

parser.add_argument(
    "--percent",
    help="Percentage of data to average",
    type=float,
    default=0.8
)

parser.add_argument(
    "--runtime",
    help="MD runtime",
    type=int,
    default=100e6
)

args = parser.parse_args()

# Change the setting of an option in an mdp file
def change_opt(buff, opt, val):
    for l in range(len(buff)):
        s = re.sub(\
            r"({0}\s+[=])(\s+)(\s?)([a-zA-Z0-9-_.]+)".format(opt),\
            r"\1\3 {0}".format(val),\
            buff[l])
        buff[l] = s


# Get the current setting of an option in an mdp file 
def get_opt(buff, opt):
    for l in range(len(buff)):
        s = re.search(\
            r"({0}\s+[=])(\s+?)([a-zA-Z0-9-_.]+)".format(opt),\
            buff[l])
        if s:
            return s.group(3)


# Load JSON template 
root = {} 
with open(path.join(args.t, "{0}.json".format(prefix))) as f:
    root = json.load(f)

# Generate appropriate number of drivers. 
root["driver"] = []
for i in range(nwalkers):
    driver = {}
    driver["number processors"] = int(args.nproc)
    driver["type"] = "Gromacs"
    driver["MDSteps"] = int(args.runtime)
    driver["logfile"] = "node-{0}".format(i)
    root["driver"].append(driver)

# Delete first CV because we're not interested for pre-elastic. 
del root["CVs"][0]
del root["method"]["ksprings"][0]
del root["method"]["centers"][0]

# Get list of dat file prefixes from source folder 
fnames = [path.splitext(f)[0] for f in os.listdir(path.join(cdir,args.s)) if f.endswith('.dat')]

# Task list 
tlist = []

# Topology file 
topf = path.join(cdir, "forcefield", "topol.top")

# Loop through file prefixes 
for fname in fnames:
    # Load mdp file for given prefix.
    mdp = open(path.join(args.s, "{0}.mdp".format(fname)), "r").readlines()
    
    # Make sure it's NVT.
    change_opt(mdp, "pcoupl", "No")

    # Change output frequency.
    change_opt(mdp, "nstxout", int(1e5))
    change_opt(mdp, "nstvout", int(1e5))
    change_opt(mdp, "nstfout", int(1e5))
    change_opt(mdp, "nstxtcout", int(1e5))
    change_opt(mdp, "nstenergy", int(1e5))
    change_opt(mdp, "gen-vel", "yes")

    # Write output mdp
    mdpf = path.join(args.o, "{0}.mdp".format(fname))
    with open(mdpf, "w") as f:
        f.writelines(mdp)
    
    # Load dat file and get average box length. 
    dat = np.loadtxt(path.join(args.s, "{0}.dat".format(fname)), skiprows=1)

    # Get beginning row number.
    n = int(dat.shape[0]*(1-args.percent))

    # Mean box length. Convert from angstrom to nm.
    L = 0.1*np.mean(dat[n:,1])

    # Generate new gro file.
    tpri = path.join(args.s, "{0}.tpr".format(fname))
    groi = path.join(args.s, "{0}.gro".format(fname))
    grof = path.join(args.o, "{0}-init.gro".format(fname))
    p = proc.Popen(\
        [gmx, "trjconv", "-f", groi, "-s", tpri, "-box", "{0:.5f}".format(L), "{0:.5f}".format(L), "{0:.5f}".format(L), "-o", grof],\
        stdin=proc.PIPE)
    p.communicate(input=b'0\n')

    # Create tpr file.
    tprf = path.join(args.o, "{0}.tpr".format(fname))
    proc.call(\
        [gmx, "grompp", "-f", mdpf, "-c", grof, "-p", topf, "-o", tprf],\
        shell=False)
    
    # Create individual tpr files for each walker 
    for i in range(nwalkers):
        tprf = path.join(args.o, "{0}{1:d}.tpr".format(fname, i))
        proc.call(\
            [gmx, "grompp", "-f", mdpf, "-c", grof, "-p", topf, "-o", tprf],\
            shell=False)

    # Add to task list 
    tlist.append("{0}\n".format(fname))

    # Set JSON parameters
    root["observers"][0]["file name"] = "{0}.chkpt".format(fname)
    root["inputfile"] = "{0}.tpr".format(fname)

    # Set restriction region based on box dimension.
    # We want 10% on either side.
    root["CVs"][0]["restriction"]["min"] = 0.9*L
    root["CVs"][0]["restriction"]["max"] = 0.1*L

    # Set output file name.
    root["method"]["file name"] = "{0}-umbrella.dat".format(fname)

    # Write json file.
    jsonf = path.join(args.o, "{0}.json".format(fname))
    with open(jsonf, "w") as f:
        json.dump(root, f, indent=4, separators=(',', ': '))

# Write task list 
with open(path.join(args.o, taskfile), "w") as f:
    f.writelines(tlist)

