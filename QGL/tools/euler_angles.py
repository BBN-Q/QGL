"""
Tools for creating Euler-angle based gates. 

Original Author: Guilhem Ribeill, Luke Govia

Copyright 2020 Raytheon BBN Technologies

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

   http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
"""

import numpy as np 
from scipy.linalg import expm
from .Cliffords import C1
from .PulsePrimitives import *
from .tools.matrix_tools import *

def xyx_unitary(α, β, γ):
    """Unitary decomposed as Rx, Ry, Rx rotations. 
        Angles are in matrix order, not in circuit order!
    """
    return expm(-0.5j*α*pX)@expm(-0.5j*β*pY)@expm(-0.5j*γ*pX)

def zyz_unitary(ϕ, θ, λ):
    """Unitary decomposed as Rz, Ry, Rz rotations. 
        Angles are in matrix order, not in circuit order!
    """
    return expm(-0.5j*ϕ*pZ)@expm(-0.5j*θ*pY)@expm(-0.5j*λ*pZ)

def diatomic_unitary(a, b, c):
    """Unitary decomposed as a diatomic gate of the form  
        Ztheta + X90 + Ztheta + X90 + Ztheta
    """
    X90 = expm(-0.25j*np.pi*pX)
    return expm(-0.5j*a*pZ)@X90@expm(-0.5j*b*pZ)@X90@expm(-0.5j*c*pZ)

def zyz_angles(U):
    """Euler angles for a unitary matrix U in the sequence Z-Y-Z.
        Note that angles are returned in matrix multiplication, not circuit order.
    """
    assert U.shape == (2,2), "Must use a 2x2 matrix!"
    k = 1.0/np.sqrt(np.linalg.det(U))
    SU = k*U 
    θ = 2 * np.arctan2(np.abs(SU[1,0]), np.abs(SU[0,0]))
    a = 2 * np.angle(SU[1,1])
    b = 2 * np.angle(SU[1,0])
    ϕ = (a + b) * 0.5
    λ = (a - b) * 0.5
    return (ϕ, θ, λ)

def _mod_2pi(angle):
    if angle > np.pi:
        angle -= 2*np.pi
    if angle < -np.pi:
        angle += 2*np.pi
    return angle

def xyx_angles(U):
    """Euler angles for a unitary matrix U in the sequence X-Y-X.
        Note that angles are returned in matrix multiplication, not circuit order.
        We make use of the identity:
            Rx(a)Ry(b)Rx(c) = H Rz(a) Ry(-b) Rz(c) H
    """
    H = np.array([[1., 1.], [1., -1.]], dtype=np.complex128)/np.sqrt(2)
    ϕ, θ, λ = zyz_angles(H@U@H)
    return (_mod_2pi(ϕ), _mod_2pi(-1.0*θ), _mod_2pi(λ))

def diatomic_angles(U):
    ϕ, θ, λ = zyz_angles(U)
    a = ϕ
    b = np.pi - θ
    c = λ - np.pi
    return (a, b, c)