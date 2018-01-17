import re
import argparse
import subprocess as proc
from os import path
import numpy as np

# Settings 

# - Current path for relative directories
cdir = path.dirname(path.abspath(__file__))

# - path to gmx executable (Gromacs 5.x)
gmx = path.join(cdir, 'gmx')

# - Prefix for source files
prefix = 'template'

# - Task file containing tpr prefixes for deffnm call
taskfile = 'args.list'

# Parser
parser = argparse.ArgumentParser()

parser.add_argument(
    "Tmin",
    help="Minimum temperature",
    type=float
)

parser.add_argument(
    "Tmax",
    help="Maximum temperature",
    type=float
)

parser.add_argument(
    "dt",
    help="Temperature increment",
    type=float
)

parser.add_argument(
    "-o",
    help="Output directory",
    type=str,
    default="npt"
)

parser.add_argument(
    "-s",
    help="Source directory",
    type=str,
    default="template"
)

parser.add_argument(
    "-p", 
    help="Prefix",
    type=str,
    default="npt"
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

# Define temperature interval 
Ts = np.arange(args.Tmin, args.Tmax, args.dt)

# Load template mdp and gro files
with open(path.join(args.s, "{0}.mdp".format(prefix)), "r") as f:
    mdp = f.readlines()

# Turn on barostat
change_opt(mdp, "pcoupl", "Parrinello-Rahman")

# Task list
tlist = []

# Topology file 
topf = path.join(cdir, "../forcefield", "topol.top")

# gro file 
grof = path.join(args.s, "{0}.gro".format(prefix))

# Go through temperatures and adjust settings
for T in Ts:
    # Create file name prefix
    pref = '{0}-{1:.2f}'.format(args.p, T)
    tlist.append(pref+"\n")

    # Set target temperature
    change_opt(mdp, "ref-t", T)

    # Write files
    mdpf = path.join(args.o, "{0}.mdp".format(pref))
    with open(mdpf, "w") as f:
        f.writelines(mdp)

    # tpr file 
    tprf = path.join(args.o, "{0}.tpr".format(pref))

    # Call Gromacs preprocessor
    proc.call(\
        [gmx, "grompp", "-f", mdpf, "-c", grof, "-p", topf, "-o", tprf],\
        shell=False)

# Write task list 
with open(path.join(args.o, taskfile), "w") as f:
    f.writelines(tlist)

    



