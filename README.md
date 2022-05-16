##### SlopeInit
a tool for editing SALEc checkpoint files.
###### COMPONENT
- [SlopeInit] calculate the pressure and stress distribution in a slope topography with given gravity
- [PostProcess] interpolate the .vtu output into SALEc checkpoint
###### REQUIRMENTS
- vtk package(>=9.1.0)
- dealii package(>9.3.0)
- cmake(>=3.16)
- gcc and mpicc
###### INPUTS
- [SlopeInit]'s input is similar with SALEc, but with additional parameters
in *cycle* field as follows.
  ```shell
   [numerical]
    chkload    = off
    chklpre    = export_check.0000  /* prefix of checkpoint files to be loaded */
    chkstore   = on
    chkspre    = export_check /* prefix of checkpoint produced in runtime */
    chkstep    = 5000 /* step interval of checkpoint */
  ```
- [PostProcess] named postprocess.inp in default
```shell
# this is input for SlopeInit/PostProcess
[process]
    pvtu = solution_00.pvtu
    chkprefix  = export_check.0000
    outprefix  = reset_check.0000
    nrank = 8
    nchk  = 8
```
##### COMPILE
```shell
mkdir -p <directory-of-source-files>/build
cd <directory-of-source-files>/build
cmake .. && make 
# the executable files can be found in <directory-of-source-files>/build/bin
```