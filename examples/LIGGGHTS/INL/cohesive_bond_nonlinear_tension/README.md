This is work in progress. Use caution.

In order to use the nonlinear tension bond model, which is still under development, calibration, and validation, follow the following steps.

* Navigate to the repository.
* `cp examples/LIGGGHTS/INL/cohesive_bond_nonlinear_tension/src/cohesion_model_bond_nonlinear_tension.h src/`
* `cd src/`
* `mv cohesion_model_bond_nonlinear_tension.h cohesion_model_bond_nonlinear.h`
* Compile or recompile the code.

The above action will overwrite the existing file `cohesion_model_bond_nonlinear.h`. To restore the original file, you can do `git checkout cohesion_model_bond_nonlinear.h` 
