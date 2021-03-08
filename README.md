# LIGGGHTS-INL
Granular Flow 

**LIGGGHTS-INL is a capability-extended version of LIGGGHTS**. LIGGGHTS is an Open Source Discrete Element Method Particle Simulation Software ([https://www.cfdem.com/media/DEM/docu/Manual.html](https://www.cfdem.com/media/DEM/docu/Manual.html)). LIGGGHTS stands for LAMMPS improved for general granular and granular heat transfer simulations. LAMMPS is a classical molecular dynamics simulator. It is widely used in the field of Molecular Dynamics ([https://lammps.sandia.gov/doc/Manual.html](https://www.cfdem.com/media/DEM/docu/Manual.html)). Thanks to physical and algorithmic analogies, LAMMPS is a very good platform for DEM simulations. LAMMPS offers a GRANULAR package to perform these kinds of simulations. LIGGGHTS aims to improve those capability with the goal to apply it to industrial applications.

## Citing LIGGGHTS-INL

If you are using results generated with LIGGGHTS-INL, please cite the following articles:

Yidong Xia, Zhengshou Lai, Tyler Westover, Jordan Klinger, Hai Huang, Qiushi Chen. "Discrete element modeling of deformable pinewood chips in cyclic loading test." *Powder Technology* 345 (2019): 1-14. [https://doi.org/10.1016/j.powtec.2018.12.072](https://doi.org/10.1016/j.powtec.2018.12.072)

Yuan Guo, Qiushi Chen, Yidong Xia, Tyler Westover, Sandra Eksioglu, Mohammad Roni. "Discrete element modeling of switchgrass particles under compression and rotational shear," *Biomass & Bioenergy* 141 (2020): 105649. [https://doi.org/10.1016/j.biombioe.2020.105649](https://doi.org/10.1016/j.biombioe.2020.105649)

## Documentation

Notice that the website for LIGGGHTS-PUBLIC documentation ([https://www.cfdem.com/media/DEM/docu/Manual.html](https://www.cfdem.com/media/DEM/docu/Manual.html)) does not contain the extended capabilities in LIGGGHTS-INL, e.g., the bonded-sphere model. 

Documentation of the bonded-sphere model:

`<LIGGGHTS-INL repository>/doc/gran_cohesion_bond.html`

User examples of the bonded-sphere model:

`<LIGGGHTS-INL repository>/examples/LIGGGHTS/INL/cohesive_bond `

## Instructions for Code Compilation

### Linux

#### Ubuntu 20.04

Install OpenMPI: `sudo apt-get install libopenmpi-dev openmpi-bin`

Install VTK 7.1: `sudo apt-get install libvtk7-dev`

Navigate to `<LIGGGHTS-INL repository>/src/`

Compile the code with VTK: `make -j4 mpi_ubuntu-20.04`

If compilation is successful, an executable will be generated: `lmp_mpi_ubuntu-20.04` 


### macOS

Use Homebrew to install OpenMPI.

**Special notice**: Do not use Homebrew to install VTK. The VTK libraries from Homebrew were pre-built without MPI enabled. So a number of VTK libraries needed by LIGGGHTS-INL are missing, e.g., vtkIOParallel, vtkIOParallel, vtkParallelCore, vtkParallelMPI, and vtkIOParallelXML. We suggest users build VTK 7.1 from source (with MPI enabled and in Release mode) and install the libraries in a common directory. To build VTK from source, additional prerequisite libraries are needed and can be installed from Homebrew.

An example of Makefile for compiling LIGGGHTS-INL on macOS can be found at:

`<LIGGGHTS-INL repository>/src/MAKE/Makefile.mpi_vtk-7.1`

(Detailed instruction to be updated ...)

## Examples of LIGGGHTS-INL Capabilities

### Strain-hardening nonlinear normal contact

An example of controlled collision between two spherical particles can be found at:

`<LIGGGHTS-INL repository>/examples/LIGGGHTS/INL/normal_contact_hysteretic_nonlinear1`

The simulation result is shown in the following plot:

<img src="figs/fig_strain_hardening_nonlinear_normal_contact.png">

### Elastoplastic bond normal stiffness

An example of microfiber made of five bonded spheres can be found at:

`<LIGGGHTS-INL repository>/examples/LIGGGHTS/INL/elastoplastic_bond/chain_bending_um_2`

The simulation result is shown in the following plot:

<img src="figs/fig_microfiber.png">

## Other Software
Idaho National Laboratory is a cutting edge research facility which is a constantly producing high quality research and software. Feel free to take a look at our other software and scientific offerings at:

[Primary Technology Offerings Page](https://www.inl.gov/inl-initiatives/technology-deployment)

[Supported Open Source Software](https://github.com/idaholab)

[Raw Experiment Open Source Software](https://github.com/IdahoLabResearch)

[Unsupported Open Source Software](https://github.com/IdahoLabCuttingBoard)

## License

Copyright 2020 Battelle Energy Alliance, LLC

Licensed under the GPL v2 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

 https://opensource.org/licenses/GPL-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

Licensing
-----
This software is licensed under the terms you may find in the file named "LICENSE" in this directory.


Developers
-----
By contributing to this software project, you are agreeing to the following terms and conditions for your contributions:

You agree your contributions are submitted under the GPL v2 license. You represent you are authorized to make the contributions and grant the license. If your employer has rights to intellectual property that includes your contributions, you represent that you have received permission to make contributions and grant the required license on behalf of that employer.

# Content below is from LIGGGHTS
----------------------------------------------------------------------
This is the

██╗     ██╗ ██████╗  ██████╗  ██████╗ ██╗  ██╗████████╗███████╗
██║     ██║██╔════╝ ██╔════╝ ██╔════╝ ██║  ██║╚══██╔══╝██╔════╝
██║     ██║██║  ███╗██║  ███╗██║  ███╗███████║   ██║   ███████╗
██║     ██║██║   ██║██║   ██║██║   ██║██╔══██║   ██║   ╚════██║
███████╗██║╚██████╔╝╚██████╔╝╚██████╔╝██║  ██║   ██║   ███████║
╚══════╝╚═╝ ╚═════╝  ╚═════╝  ╚═════╝ ╚═╝  ╚═╝   ╚═╝   ╚══════╝®

DEM simulation engine, released by 
DCS Computing GmbH, Linz, Austria
www.dcs-computing.com, office@dcs-computing.com

LIGGGHTS® is open-source, distributed under the terms of the GNU Public 
License, version 2 or later.

LIGGGHTS® is part of CFDEM®project: 
www.liggghts.com | www.cfdem.com

Core developer and main author:
Christoph Kloss, christoph.kloss@dcs-computing.com

LIGGGHTS® and CFDEM® are registered trade marks of DCS Computing GmbH, 
the producer of the LIGGGHTS® software and the CFDEM®coupling software
See http://www.cfdem.com/terms-trademark-policy for details.

----------------------------------------------------------------------
Copyright 2012-     DCS Computing GmbH, Linz
Copyright 2009-2015 JKU Linz
Some parts of LIGGGHTS® are based on LAMMPS and Copyright on these
parts is held by Sandia Corporation and other parties. Info on LAMMPS below
Some parts of LIGGGHTS® are contributied by other parties, which are
holding the Copyright. This is listed in each file of the distribution.
----------------------------------------------------------------------

The LIGGGHTS® distribution includes the following files and directories:

README          this file
LICENSE         the GNU General Public License (GPL)
doc             documentation
examples        simple example simulation setups
lib             libraries LIGGGHTS® can be linked with
python          Python wrapper on LIGGGHTS® as a library
src             source files

Point your browser at any of these files to get started:

doc/Manual.html	           the manual
doc/Section_intro.html	   hi-level introductio
doc/Section_start.html	   how to build and use

----------------------------------------------------------------------

Some parts of LIGGGHTS® are based on LAMMPS
LAMMPS stands for Large-scale Atomic/Molecular Massively Parallel
Simulator. 

LAMMPS is Copyright (2003) Sandia Corporation.  Under the terms of Contract
DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
certain rights in this software.  This software is distributed under
the GNU General Public License.

LAMMPS is a classical molecular dynamics simulation code designed to
run efficiently on parallel computers.  It was developed at Sandia
National Laboratories, a US Department of Energy facility, with
funding from the DOE.  It is an open-source code, distributed freely
under the terms of the GNU Public License (GPL).

The primary author of LAMMPS is Steve Plimpton, who can be emailed
at sjplimp@sandia.gov.  The LAMMPS WWW Site at lammps.sandia.gov has
more information about the code and its uses.

----------------------------------------------------------------------

Core developer and main author of LIGGGHTS®:
Christoph Kloss, christoph.kloss@dcs-computing.com
