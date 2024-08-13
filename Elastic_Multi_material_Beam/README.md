# Elastic Multimaterial Beam

## Description of the problem

This example aims in providing an example of a 3D beam of dimensions $`40\times5\times5`$ composed of two materials.

The beam is subdivided into two subdomains of dimensions $`20\times5\times5`$ with different material properties. Both sides respect a same constitutive law and a mapping is applied on the material properties.

The beam is clamped on its left face (Dirichlet boundary condition) and a vertical traction force is applied on its right face (Neumann Boundary condition).

Even though the proposed problem is linear, a non-linear solver is used for the example. Linear solvers have been proposed for the Stokes and Thermodynamic problems.

## Implementation

