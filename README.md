# LIGGGHTS-INL

<img src="figs/fig_liggghts_inl_banner.png">

**LIGGGHTS-INL** is a capability-extended adaptation of the [LIGGGHTS](https://www.cfdem.com/media/DEM/docu/Manual.html) Open Source Discrete Element Method (DEM) Particle Simulation Software based on LIGGGHTS release version 4.0.0.

**Why LIGGGHTS-INL?** The public version of LIGGGHTS DEM software, [**LIGGGHTS-PUBLIC**](https://github.com/CFDEMproject/LIGGGHTS-PUBLIC), is maintained by LIGGGHTS community users. LIGGGHTS-INL had been but is no longer a fork of LIGGGHTS-PUBLIC. LIGGGHTS-INL adopts some updates in LIGGGHTS-PUBLIC to fix bugs in the code and offers extended capabilities such as bonded-sphere model (or bonded-particle model in some literature), elastoplastic bond normal stiffness, and strain-hardening nonlinear normal contact model. Those capabilities have been implemented in LIGGGHTS-INL originally for biomass granular flow and structural biomechanics simulations, and can be used for broader applications upon adaptation.  For many common scenarios of granular flow simulations that do not require those specialized features, LIGGGHTS-INL can be used interchangeably with LIGGGHTS-PUBLIC, though some input command syntaxes have become different due to our independent maintenance. LIGGGHTS-INL will routinely release user tutorials and examples to improve user experience.

## Installing LIGGGHTS-INL

Supported operating systems: **Linux** and **macOS**. We provide [detailed instructions for installing LIGGGHTS-INL on the Linux Ubuntu LTS and macOS releases](/compile/README.md).

## Citing LIGGGHTS-INL

Reference articles with results generated using LIGGGHTS-INL are listed below:

**Static angle of repose**

* A. Hamed et al. Flowability of Crumbler rotary shear size-reduced granular biomass: An experiment-informed modeling study on the angle of repose. [*Frontiers in Energy Research, section Bioenergy and Biofuel* (2022).](https://doi.org/10.3389/fenrg.2022.859248) (<span style="color: green;">**Open Access**</span>)

**Uniaxial loading test**

* F. Chen et al. A set of hysteretic nonlinear contact models for DEM: Theory, formulation, and application for lignocellulosic biomass. [*Powder Technology* 397 (2022): 117100.](https://doi.org/10.1016/j.powtec.2021.117100) (Download [Authors' manuscript PDF](https://www.researchgate.net/publication/357393650_A_set_of_hysteretic_nonlinear_contact_models_for_DEM_Theory_formulation_and_application_for_lignocellulosic_biomass))
* Y. Xia et al. Discrete element modeling of deformable pinewood chips in cyclic loading test. [*Powder Technology* 345 (2019): 1-14.](https://doi.org/10.1016/j.powtec.2018.12.072) (Download [Authors' manuscript PDF](https://www.researchgate.net/publication/329869479_Discrete_element_modeling_of_deformable_pinewood_chips_in_cyclic_loading_test))

**Ring shear test**

* Y. Guo et al. Discrete element modeling of switchgrass particles under compression and rotational shear. [*Biomass & Bioenergy* 141 (2020): 105649.](https://doi.org/10.1016/j.biombioe.2020.105649) (Download [Authors' manuscript PDF](https://www.researchgate.net/publication/342170187_Discrete_element_modeling_of_switchgrass_particles_under_compression_and_rotational_shear))

**Hopper discharge flow**

* Z. Lai et al. Discrete element modeling of granular hopper flow of irregular-shaped deformable particles. [*Advanced Powder Technology* 34 (2023): 104106.](https://doi.org/10.1016/j.apt.2023.104106) (Download [Authors' manuscript PDF](https://www.researchgate.net/publication/371447845_Discrete_element_modeling_of_granular_hopper_flow_of_irregular-shaped_deformable_particles))
* F. Chen et al. Hopper discharge flow dynamics of milled pine and prediction of process upsets using the discrete element method. [*Powder Technology* 415 (2023): 118165.](https://doi.org/10.1016/j.powtec.2022.118165) (Download [Authors' manuscript PDF](https://www.researchgate.net/publication/366219764_Hopper_discharge_flow_dynamics_of_milled_pine_and_prediction_of_process_upsets_using_the_discrete_element_method))

**Screw feeder**

* A. Hamed et al. Particle size and shape effect of Crumbler® rotary shear-milled granular woody biomass on the performance of Acrison® screw feeder: A computational and experimental investigation. [*Powder Technology* 427 (2023): 118707.](https://doi.org/10.1016/j.powtec.2023.118707) (<span style="color: green;">**Open Access**</span>)

**FEM vs. DEM: granular flow benchmark case studies**

* W. Jin et al. On the fidelity of computational models for the flow of milled loblolly pine: A benchmark study on continuum-mechanics models and discrete-particle models. [*Frontiers in Energy Research, section Bioenergy and Biofuel* (2022).](https://doi.org/10.3389/fenrg.2022.855848) (<span style="color: green;">**Open Access**</span>)

**Structural biomechanics**

* Q. Sun et al. Reverse scaling of a bonded-sphere DEM model: Formulation and application to lignocellulosic biomass microstructures. [*Powder Technology* 409 (2022): 117797.](https://doi.org/10.1016/j.powtec.2022.117797) (Download [Authors' manuscript PDF](https://www.researchgate.net/publication/362419823_Reverse_scaling_of_a_bonded-sphere_DEM_model_Formulation_and_application_to_lignocellulosic_biomass_microstructures))
* Y. Guo et al. A nonlinear elasto-plastic bond model for the discrete element modeling of woody biomass particles. [*Powder Technology* 385 (2021): 557-571.](https://doi.org/10.1016/j.powtec.2021.03.008) (Download [Authors' manuscript PDF](https://www.researchgate.net/publication/350048423_A_nonlinear_elasto-plastic_bond_model_for_the_discrete_element_modeling_of_woody_biomass_particles))

## Tutorials for Granular Material Flow Characterization

We are in the process of releasing tutorials for beginners to use DEM simulation as a numerical tool for granular material flow characterization. Examples can be found here: [LIGGGHTS-INL tutorials](/examples/LIGGGHTS/INL_tutorials).

## LIGGGHTS-INL Documentation and Extended Capabilities

Refer to the [LIGGGHTS-PUBLIC documentation](https://www.cfdem.com/media/DEM/docu/Manual.html) for the common **LIGGGHTS** features. Documentation of some of the **LIGGGHTS-INL** extended capabilities is in this repository, e.g., the bonded-sphere model ([HTML documentation](/doc/gran_cohesion_bond.html) and [user examples](/examples/LIGGGHTS/INL/cohesive_bond)).


**LIGGGHTS-INL** provides original **nonlinear contact and bond stiffness models**. No user manual has yet been created for these models. Refer to the listed articles for more information. User examples are introduced below:

* **Strain-hardening nonlinear normal contact**: [An example of collision between two spherical particles](/examples/LIGGGHTS/INL/normal_contact_hysteretic_nonlinear1)
<img src="figs/fig_nonlinear_contact.png">

* **Elastoplastic bond normal stiffness** (a): [An example of macro-fiber made of five bonded spheres](/examples/LIGGGHTS/INL/cohesive_bond_nonlinear_compression/chain_bending_mm_2)
<img src="figs/fig_string_mm.png">

* **Elastoplastic bond normal stiffness** (b): [An example of microfiber made of five bonded spheres](/examples/LIGGGHTS/INL/cohesive_bond_nonlinear_compression/chain_bending_um_2)
<img src="figs/fig_string_um.png">

* **Micro-biomechanics**: [An example of compression on a pine wood particle microstructure](/examples/LIGGGHTS/INL/microstructure_compression)
<img src="figs/fig_microstructure_compression.png">

## Auxiliary 3D Image Processing and Analysis Tools

- A set of MATLAB codes for 3D image-based porosity analysis: [tools/PorosityAnalysis3D](tools/PorosityAnalysis3D)
- A concise FIJI user tutorial for 3D image binarization: [tools/FIJI](tools/FIJI)

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
