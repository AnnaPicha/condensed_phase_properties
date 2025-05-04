
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

theory      = sys.argv[1]          # mm/ani2x/mace-s/mace-m/mace-l
system_name = sys.argv[2]          # mono/moh323
ensemble    = sys.argv[3]
step        = int(sys.argv[4])

type = "equi"
path = "/site/raid7/anna/cp_prop/methanol/simulation"


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

psf_file = f"/site/raid7/anna/cp_prop/methanol/coor_moh/{system_name}.psf"
psf      = CharmmPsfFile(psf_file)
crd      = CharmmCrdFile(f"/site/raid7/anna/cp_prop/methanol/coor_moh/{system_name}.crd")

if system_name!="mono":
    # if ensemble=="NVT":
    box_size = 27.97898408 
    psf.setBox(box_size * unit.angstrom, box_size * unit.angstrom, box_size * unit.angstrom) #     water at 300 Kelvin = 26.85 grad celsius is 0.996571 g/mL
    # else:
    #     psf=gen_box(psf, crd)


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

integrator = LangevinIntegrator(
    temp * unit.kelvin, 1 / unit.picosecond, dt * unit.picoseconds
)


platform = Platform.getPlatformByName("CUDA")
prop = dict(CudaPrecision="mixed")
simulation = Simulation(psf.topology, system, integrator, platform, prop)

simulation.context.setPositions(crd.positions)
simulation.context.setVelocitiesToTemperature(temp * unit.kelvin)


# # Calculate initial system energy
print("\nInitial system energy")
print(simulation.context.getState(getEnergy=True).getPotentialEnergy())


simulation.reporters.append(DCDReporter(f"{path}/trajectories/{theory}/{system_name}_{theory}_{ensemble}_{type}.dcd", 100, enforcePeriodicBox=True))
simulation.reporters.append(
    StateDataReporter(
        f"{path}/trajectories/{theory}/{system_name}_{theory}_{ensemble}_{type}.csv",
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

if system_name != "mono":
    # Write restart file
    restart_file = f"{path}/rst/equi_{ensemble}_{system_name}_{theory}.rst"
    system_file = f"{path}/rst/equi_{ensemble}_state_{system_name}_{theory}.txt"
    if not (restart_file): restart_file = f"{path}/rst/equi_{ensemble}_{system_name}_{theory}.rst"
    if restart_file:
        state = simulation.context.getState( getPositions=True, getVelocities=True )
        with open(restart_file, 'w') as f:
            f.write(XmlSerializer.serialize(state))
        with open(system_file, 'w') as f:
            f.write(XmlSerializer.serialize(system))