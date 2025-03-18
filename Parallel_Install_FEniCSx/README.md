# Installation

This procedure was used to install FEniCSx on the HPC Cassiopee (RED HAT Enterprise Linux 8.10).

* Get version 0.23.0 of spack:
```bash
wget https://github.com/spack/spack/releases/download/v0.23.0/spack-0.23.0.tar.gz
```

* untar the archive:
```bash
tar -xf spack-0.23.0.tar.gz
```

* Optional: disable the local configuration of spack (current user home directory). Useful if you have multiple versions of spack
```bash
export SPACK_DISABLE_LOCAL_CONFIG=true
```

* activate spack:
```bash
cd spack-0.23.0
source ./share/spack/setup-env.sh
```

* add system compiler to spack (need to be installed in your os):
```bash
spack compiler find 
```

* optional: add a decent compiler (useful if you have an old distribution) and add it as compiler to spack:
```bash
spack install gcc@12
spack compiler find $(spack location -i gcc@12)
```

* create a spack environment (useful if you want to load several packages at one time):
```bash
spack env create MyFenicsxEnv
```

* activate the environment, add some packages, and install them:
```bash
spack env activate MyFenicsxEnv
spack add petsc +fortran +mumps +trilinos +superlu-dist
spack add python
spack add py-fenics-dolfinx@0.9.0
spack add py-numpy
spack add py-matplotlib
spack add py-pandas
spack add py-gmsh
spack install
```
Congratulation, you have installed fenicsx 0.9.0, working with mpi!

* Create the porous_fenicsx package
```bash
cd $PATH/porous_fenicsx
python -m pip install .
```
* deactivate your environment:
```bash
spack env deactivate
```

# Utilisation

* activate spack in your session (can be added to your .bashrc):
```bash
source FullPathToSpack/share/spack/setup-env.sh
```
* activate your environment:
```bash
spack env activate MyFenicsxEnv
```
* launch your script in parallel (replace NCPUs with the number of CPUs you want):
```bash
mpirun -n NCPUs python myscript.py
```
* deactivate your environment:
```bash
spack env deactivate
```




