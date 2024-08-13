# Elastic Multimaterial Beam

## Description of the problem

This example aims in providing an example of a 3D beam of dimensions $`40\times5\times5`$ composed of two materials.

The beam is subdivided into two subdomains of dimensions $`20\times5\times5`$ with different material properties. Both sides respect a same constitutive law and a mapping is applied on the material properties.

The beam is clamped on its left face (Dirichlet boundary condition) and a vertical traction force is applied on its right face (Neumann Boundary condition).

The objective is to find the resulting displacement. Even though the proposed problem is linear, a non-linear solver is used for the example. Linear solvers have been proposed for the Stokes and Thermodynamic problems.

## Implementation

Computing the Finite element problem within FEniCSx in python requires to load the libraries:

```python
import dolfinx
from dolfinx.fem.petsc import NonlinearProblem
from dolfinx.nls.petsc import NewtonSolver
import ufl
import basix
import petsc4py
import mpi4py
import numpy
import pyvista
``` 
One can assess the version of FEniCSx with the following:
```python
print("Dolfinx version is:",dolfinx.__version__)
```
