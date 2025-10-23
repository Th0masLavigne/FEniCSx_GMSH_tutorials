import dolfinx
import basix.ufl
import ufl
import random as random
import numpy as np
import os
import gmsh
import sys
import petsc4py
import matplotlib.pyplot as plt
#
from mpi4py import MPI
from petsc4py import PETSc
from dolfinx import fem, mesh, io, plot
from dolfinx.nls.petsc import NewtonSolver
from dolfinx.fem.petsc import NonlinearProblem
def Heaviside(H_, sigma, subspace_number, sigma_min, sigma_max):
    """
    Computes the smoothed Heaviside function for a given subspace of sigma.

    Parameters:
        H_ (Function): Function to store the Heaviside values.
        sigma (Function): Input function to apply the Heaviside transformation.
        subspace_number (int): Index of the subspace.
        sigma_min (float): Minimum threshold for sigma.
        sigma_max (float): Maximum threshold for sigma.

    Returns:
        H_ (Function): Updated Heaviside function.
    """
    fsub = sigma.sub(subspace_number).collapse()
    tmp = np.clip(fsub.x.array, sigma_min, sigma_max)
    H_.x.petsc_vec[:] = 0.5 - 0.5 * np.cos(np.pi * (tmp - sigma_min) / (sigma_max - sigma_min))
    return H_

def axisymmetric_grad(u,r):
    """
    Gradient in cylindrical coordinates: (r, θ, z).
    u: 3D vector field with components [u_r, u_θ, u_z]
    [[∂u_r/∂r   , (1/r)*(∂u_r/∂θ) - u_θ/r           , ∂u_r/∂z)     ],
     [∂u_θ/∂r   , (1/r)*(∂u_θ/∂θ) - u_r/r           , ∂u_θ/∂z      ],
     [∂u_z/∂r   , (1/r)*(∂u_z/∂θ)                   , ∂u_z/∂z      ]]
    #
    In our case : Custom gradient in axisymmetric cylindrical coordinates: (r, z).
    u: 2D vector field with components [u_r, u_z]
    [[∂u_r/∂r   , 0             , ∂u_r/∂z)      ],
     [0         , - u_r/r       , 0             ],
     [∂u_z/∂r   , 0             , ∂u_z/∂z       ]]
    """
    return ufl.as_tensor([[ufl.Dx(u[0], 0)  , 0               , ufl.Dx(u[0], 1)],
                          [0                , -u[0]/r         , 0              ],
                          [ufl.Dx(u[1], 0)  , 0               , ufl.Dx(u[1], 1)]])                          
#
def axisymmetric_div(u, r):
    """
    Divergence in cylindrical coordinates: (r, θ, z).
    u: 3D vector field with components [u_r, u_θ, u_z]
    div(u) = (1/r)*[∂(u_r*r)/∂r] + (1/r)*[∂(u_θ)/∂θ] +  ∂u_z/∂z
    #
    In our case : Custom divergence in axisymmetric cylindrical coordinates: (r, z).
    u: 2D vector field with components [u_r, u_z]
    div(u) = (1/r)*[∂(u_r*r)/∂r] + ∂u_z/∂z
    """
    r_safe = ufl.max_value(r, 1e-8)  # To avoid dividing by 0 the radius minumum value is 100 times smaller than the size of 1 element
    return ufl.Dx(u[0], 0) + u[0]/r_safe + ufl.Dx(u[1], 1)
#
def axisymmetric_inner_grad(u, v,r):
    """
    Inner product of gradients in axisymmetric setting.
    Equivalent to: inner(grad(u), grad(v)) + u_r * v_r / r^2
    """
    return ufl.inner(axisymmetric_grad(u,r), axisymmetric_grad(v,r))
#
def axi(integrand,r, dx):
    """
    Axisymmetric volume integration: ∫ 2πr · integrand dx
    This accounts for rotation around the axis (volume of revolution).
    """
    return r * integrand * dx

def readme_O():
    return None