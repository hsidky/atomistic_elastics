import argparse
import os 
from os import path
import json
import re 
import MDAnalysis as md
import subprocess as proc

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
    help="Gro files source directory",
    type=str,
    default="preelastic"
)

parser.add_argument(
    "--nproc",
    help="Number of processors per walker",
    type=int,
    default=6
)

parser.add_argument(
    "--runtime",
    help="MD runtime",
    type=int,
    default=700e6
)

parser.add_argument(
    "--frac",
    help="Fraction of box dimension for restriction",
    type=float,
    default=0.2
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

# Add constraints

root["constraints"] = []
constraint = {}
constraint["type"] = "Umbrella"
constraint["ksprings"] = [0, 1e6]
constraint["centers"] = [0, 1.0]
constraint["log_frequency"] = int(1000)
root["constraints"].append(constraint)

root["method"] = {}
root["method"]["type"] = "Basis"
root["method"]["cycle_frequency"] = int(10000000)
root["method"]["frequency"] = int(1)
root["method"]["weight"] = 1.0 
root["method"]["CV_coefficients"] = [int(2)]
root["method"]["CV_restraint_spring_constants"] = [0]
root["method"]["CV_restraint_minimums"] = [-0.52]
root["method"]["CV_restraint_maximums"] = [0.52]
root["method"]["tolerance"] = 1e-6
root["method"]["convergence_exit"] = True

# Add grid 
root["grid"] = {}
root["grid"]["lower"] = [-0.3]
root["grid"]["upper"] = [0.3]
root["grid"]["number_points"] = [100]
root["grid"]["periodic"] = [False]

# Get list of file prefixes from taskfile.
fnames = sorted([line.rstrip('\n') for line in open(path.join(args.s, taskfile), "r").readlines()])

# New task list 
tlist = []

# Topology file 
topf = path.join(cdir, "forcefield", "topol.top")

# Elastic modes 
modes = ["splay", "twist", "bend"]

# Director for each mode.
directors = [
    [1, 0, 0],
    [0, 1, 0],
    [1, 0, 0]
]

# Restriction dimension for each mode 
dims = ["x", "x", "z"]

# Create directory for each mode 
for mode in modes:
    if not path.exists(mode):
        os.makedirs(mode)

# Disable GMX backups.
os.environ["GMX_MAXBACKUP"] = "-1"

# Loop through file prefixes 
for fname in fnames:
    # Load mdp file for given prefix.
    mdp = open(path.join(args.s, "{0}.mdp".format(fname)), "r").readlines()
    
    # Make sure it's NVT.
    change_opt(mdp, "pcoupl", "No")

    # Change output frequency.
    change_opt(mdp, "nstxout", int(1e7))
    change_opt(mdp, "nstvout", int(1e7))
    change_opt(mdp, "nstfout", int(1e7))
    change_opt(mdp, "nstxtcout", int(1e7))
    change_opt(mdp, "nstenergy", int(1e6))
    change_opt(mdp, "gen-vel", "No")

    # Read initial gro file file to get coordinates.
    u = md.Universe(path.join(args.s, "{0}.tpr".format(fname)), path.join(args.s, "{0}-init.gro".format(fname)))
    L = 0.1*u.trajectory[0].triclinic_dimensions[0][0] # nm

    # Add to task list 
    tlist.append("{0}\n".format(fname))

    # Configure non mode dependent JSON.	
    root["observers"][0]["file name"] = "{0}.chkpt".format(fname)
    root["inputfile"] = "{0}.tpr".format(fname)

    # Go through each mode 
    for i in range(len(modes)):
        mode = modes[i]
        director = directors[i]
        dim = dims[i]
        
        # Write output mdp
        mdpf = path.join(mode, "{0}.mdp".format(fname))
        with open(mdpf, "w") as f:
            f.writelines(mdp)
        
        # Create main tpr file.
        tprf = path.join(mode, "{0}.tpr".format(fname))
        groi = path.join(args.s, "{0}0.gro".format(fname))
        proc.call(\
            [gmx, "grompp", "-f", mdpf, "-c", groi, "-p", topf, "-o", tprf],\
            shell=False)

        # Create tpr file for each walker from final gro.
        for j in range(nwalkers):
            tprf = path.join(mode, "{0}{1:1d}.tpr".format(fname, j))
            groi = path.join(args.s, "{0}{1:1d}.gro".format(fname, j))
            proc.call(\
                [gmx, "grompp", "-f", mdpf, "-c", groi, "-p", topf, "-o", tprf],\
                shell=False)
        
        # Mode specific JSON.
        root["constraints"][0]["file_name"] = "{0}-umbrella.dat".format(fname)
        root["method"]["basis_filename"] = fname
        root["method"]["coeff_filename"] = fname

        # - First CV is the central region. 
        root["CVs"][0]["director"] = director
        root["CVs"][0]["restriction"]["dimension"] = dim
        root["CVs"][0]["restriction"]["min"] = 0.5*L*(1 - args.frac)
        root["CVs"][0]["restriction"]["max"] = 0.5*L*(1 + args.frac)
        # - Second CV is the edge restriction.
        root["CVs"][1]["restriction"]["dimension"] = dim
        root["CVs"][1]["restriction"]["min"] = (1-0.5*args.frac)*L
        root["CVs"][1]["restriction"]["max"] = 0.5*args.frac*L

        # Write JSON FileExistsError
        jsonf = path.join(mode, "{0}.json".format(fname))
        with open(jsonf, "w") as f:
            json.dump(root, f, indent=4, separators=(',', ': '))
        
# Write task list
for mode in modes: 
    with open(path.join(mode, taskfile), "w") as f:
        f.writelines(tlist)



        


