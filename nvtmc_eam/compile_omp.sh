gfortran -fopenmp variables.f90 fvector.f90 filehandle.f90 gen_pos.f90 func_mc.f90 domdec.f90 dcdmod.f90 mc_eam.f90 -o mc_eam  -fcheck=bounds -ffree-line-length-256
gfortran -fopenmp variables.f90 fvector.f90 filehandle.f90 gen_pos.f90 func_anal.f90 domdec.f90 dcdmod.f90 anal_htraj.f90 -o anal_htraj  -fcheck=bounds -ffree-line-length-256
gfortran -fopenmp anal_hmap.f90 -o anal_hmap -ffree-line-length-256 -fcheck=all -Wall -g -fbacktrace
ulimit -s 1048576
