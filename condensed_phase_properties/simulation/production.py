
from openmm.app import (
    StateDataReporter,
    DCDReporter,
    CharmmPsfFile,
    CharmmParameterSet,
    Simulation,
    CharmmCrdFile,
    PME
)
from openmm import (
    unit,
    LangevinIntegrator,
    NoseHooverIntegrator,
    VerletIntegrator,
    Platform,
    XmlSerializer,
    MonteCarloBarostat,
)
from openmm.openmm import CMMotionRemover

from openmmml import MLPotential
import nutmegpotentials
import sys
import argparse
import os
import torch
import argparse

torch._C._jit_set_nvfuser_enabled(False)

def gen_box(psf, crd):
    coords = crd.positions

    min_crds = [coords[0][0], coords[0][1], coords[0][2]]
    max_crds = [coords[0][0], coords[0][1], coords[0][2]]

    for coord in coords:
        min_crds[0] = min(min_crds[0], coord[0])
        min_crds[1] = min(min_crds[1], coord[1])
        min_crds[2] = min(min_crds[2], coord[2])
        max_crds[0] = max(max_crds[0], coord[0])
        max_crds[1] = max(max_crds[1], coord[1])
        max_crds[2] = max(max_crds[2], coord[2])

    boxlx = max_crds[0] - min_crds[0]
    boxly = max_crds[1] - min_crds[1]
    boxlz = max_crds[2] - min_crds[2]

    psf.setBox(boxlx, boxly, boxlz)
    return psf

parser = argparse.ArgumentParser()
parser.add_argument('-cnt', dest='counter', help='Counter variable fpr consecutive runs', required=True)
parser.add_argument('-i', dest='inpfile', help='Input parameter file', required=False)
parser.add_argument('-p', dest='psffile', help='Input CHARMM PSF file', required=True)
parser.add_argument('-c', dest='crdfile', help='Input CHARMM CRD file (EXT verison)', required=True)
parser.add_argument('-t', dest='toppar', help='Input CHARMM-GUI toppar stream file', required=False)
parser.add_argument('-n', dest='name', help='Input base name of dcd and restart files', required=True)
parser.add_argument('-npt', dest='n_npt', help='Total number of npt runs', required=False)
parser.add_argument('-nvt', dest='n_nvt', help='Total number of nvt runs', required=False)
parser.add_argument('-m', dest='theory', help='Integration method', required=False)
parser.add_argument('-e', dest='ensemble', help='Thermodynamic ensemble', required=False)
parser.add_argument('-int', dest='int_method', help='Integrator for simulation', required=False)
parser.add_argument('-s', dest='step', help='Number of steps per run', required=False)
parser.add_argument('-sys', dest='system', help='Simulation system', required=False)

parser.add_argument('-b', dest='sysinfo', help='Input CHARMM-GUI sysinfo stream file (optional)', default=None)
parser.add_argument('-icrst', metavar='RSTFILE', dest='icrst', help='Input CHARMM RST file (optional)', default=None)
parser.add_argument('-irst', metavar='RSTFILE', dest='irst', help='Input restart file (optional)', default=None)
parser.add_argument('-ichk', metavar='CHKFILE', dest='ichk', help='Input checkpoint file (optional)', default=None)
parser.add_argument('-opdb', metavar='PDBFILE', dest='opdb', help='Output PDB file (optional)', default=None)
parser.add_argument('-minpdb', metavar='MINIMIZEDPDBFILE', dest='minpdb', help='Output PDB file after minimization (optional)', default=None)
parser.add_argument('-orst', metavar='RSTFILE', dest='orst', help='Output restart file (optional)', default=None)
parser.add_argument('-ochk', metavar='CHKFILE', dest='ochk', help='Output checkpoint file (optional)', default=None)
parser.add_argument('-odcd', metavar='DCDFILE', dest='odcd', help='Output trajectory file (optional)', default=None)
parser.add_argument('-ocsv', metavar='CSVFILE', dest='ocsv', help='Output csv file (optional)', default=None)
args = parser.parse_args()


f_name      = args.name
theory      = args.theory #sys.argv[1]          # mm/ani2x/mace-s/mace-m/mace-l
system_name = args.system #sys.argv[2]          # mono/tip125/tip572
ensemble    = args.ensemble #sys.argv[3]
int_method  = args.int_method
step        = int(args.step) #int(sys.argv[4])

type = "production" # production/NNP_equi
path = "/site/raid7/anna/cp_prop/methanol/simulation/trajectories/"


if not any([args.orst, args.ochk]):
    args.orst = path + theory + "/traj_" + ensemble + "_" + int_method + "/" + f_name + '_' + str(args.counter) + '_' + ensemble + "_" + int_method + '.rst'
if not args.odcd:
    args.odcd = path + theory + "/traj_" + ensemble + "_" + int_method + "/" + f_name + '_' + str(args.counter) + '_' + ensemble + "_" + int_method + '.dcd'
if not args.ocsv:
    args.ocsv = path + theory + "/traj_" + ensemble + "_" + int_method + "/" + f_name + '_' + str(args.counter) + '_' + ensemble + "_" + int_method + '.csv'

pcnt = int(args.counter) - 1
# hier hin schreiben: wenn der counter = 1 ist -> .rst file von der mm equilibrierung verwenden
if pcnt == 0:
    args.irst = f'/site/raid7/anna/cp_prop/methanol/simulation/rst/equi_{ensemble}_{system_name}_mm.rst'
if int(args.counter) > 1 and not any([args.icrst, args.irst, args.ichk]):
    args.irst = path + theory + "/traj_" + ensemble + "_" + int_method + "/" + f_name + '_' + str(pcnt) + '_' + ensemble + "_" + int_method + '.rst'


temp = 300 
dt   = 0.0005  #0.5 fs

# for mm runs----------------------
if system_name=="moh73":
    r_off = 0.7
    r_on = r_off - 0.1
else:
    r_off = 1.2
    r_on = r_off - 0.2
# ---------------------------------



# combination mono and mm makes no sense, pot energy = 0!

# input files for the waterbox

psf = CharmmPsfFile(args.psffile)
crd = CharmmCrdFile(args.crdfile)


if system_name!="mono":
    if ensemble=="NVT":
        box_size = 27.97898408 
        psf.setBox(box_size * unit.angstrom, box_size * unit.angstrom, box_size * unit.angstrom) #     water at 300 Kelvin = 26.85 grad celsius is 0.996571 g/mL
    else:
        psf=gen_box(psf, crd)



if theory=="mm":
    parms = ()
    parms += (f"/site/raid7/anna/cp_prop/acetone/charmm-gui-box_aco178/toppar/par_all36_cgenff.prm",)
    params = CharmmParameterSet(*parms)
    if system_name!="mono":
        system = psf.createSystem(
            params,
            nonbondedMethod=PME,
            nonbondedCutoff=r_off * unit.nanometers,
            switchDistance=r_on * unit.nanometers,
            temperature=temp,
            rigidWater=True,
        )
    else:
        system = psf.createSystem(
        params,
        # nonbondedMethod=PME,
        nonbondedCutoff=r_off * unit.nanometers,
        switchDistance=r_on * unit.nanometers,
        temperature=temp,
        rigidWater=True,
    )
else:
    if theory=="mace-s":
        potential = MLPotential('mace-off23-small')
    if theory=="mace-m":
        potential = MLPotential('mace-off23-medium')
    if theory=="mace-l":
        potential = MLPotential('mace-off23-large') 
    if theory=="ani2x":
        potential = MLPotential("ani2x")
    if theory=="nutmeg-s":
        potential = MLPotential("nutmeg-small")


    if theory=="ani2x":
        system = potential.createSystem(
            psf.topology,
            implementation='torchani',
        )
    else:
        if theory=="nutmeg-s":
            system = potential.createSystem(psf.topology, total_charge=0)
        else:
            system = potential.createSystem(
                psf.topology,
            )
    system.addForce(CMMotionRemover())
    

# only for NPT simulations
if ensemble == "NPT" and system_name!="mono":
    barostat = MonteCarloBarostat(1.0 * unit.bar, temp * unit.kelvin)
    system.addForce(barostat)

if int_method == "L":
    integrator = LangevinIntegrator(
        temp * unit.kelvin, 1 / unit.picosecond, dt * unit.picoseconds
    )
else:
    if int_method=="NH":
        integrator = NoseHooverIntegrator(
            temp * unit.kelvin, 1 / unit.picosecond, dt * unit.picoseconds
        )
    else:
        print("No valid integrator chosen! Choose L for Langevin or NH for Nose-Hoover!")


platform = Platform.getPlatformByName("CUDA")
prop = dict(CudaPrecision="mixed")
simulation = Simulation(psf.topology, system, integrator, platform, prop)

simulation.context.setPositions(crd.positions)
simulation.context.setVelocitiesToTemperature(temp * unit.kelvin)

if args.irst:
    with open(args.irst, 'r') as f:
        simulation.context.setState(XmlSerializer.deserialize(f.read()))


# # Calculate initial system energy
print(f"\nInitial system energy for {theory} box:")
print(simulation.context.getState(getEnergy=True).getPotentialEnergy())


simulation.reporters.append(DCDReporter(args.odcd, 100, enforcePeriodicBox=True))
simulation.reporters.append(
    StateDataReporter(
        args.ocsv,
        reportInterval=100,
        step=True,
        time=True,
        potentialEnergy=True,
        totalEnergy=True,
        temperature=True,
        volume=True,
        density=True,
        speed=True,
        separator="\t",
    )
)

if step > 0:
    print("\nMD run: %s steps" % step)
    
    simulation.step(step)


# Write restart file
if not (args.orst or args.ochk): args.orst = 'output.rst'
if args.orst:
    state = simulation.context.getState( getPositions=True, getVelocities=True )
    with open(args.orst, 'w') as f:
        f.write(XmlSerializer.serialize(state))



