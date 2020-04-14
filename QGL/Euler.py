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

#Pauli matrices
pX = np.array([[0, 1], [1, 0]], dtype=np.complex128)
pZ = np.array([[1, 0], [0, -1]], dtype=np.complex128)
pY = 1j * pX @ pZ
pI = np.eye(2, dtype=np.complex128)

#Machine precision
_eps = np.finfo(np.complex128).eps 

#### FUNCTIONS COPIED FROM PYGSTI
#### See: https://github.com/pyGSTio/pyGSTi

#### PYGSTI NOTICE
# Python GST Implementation (PyGSTi) v. 0.9
# Copyright 2015, 2019 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
# Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains certain rights in this software.
#### END PYGSTI NOTICE

#### PYGSTI COPYRRIGHT
# Copyright 2015, 2019 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
# Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains certain rights
# in this software.
# Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except
# in compliance with the License.  You may obtain a copy of the License at
# http://www.apache.org/licenses/LICENSE-2.0 or in the LICENSE file in the root pyGSTi directory.
#### END PYGSTI COPYRIGHT

def tracenorm(A, tol=np.sqrt(_eps)):
    """Compute the trace norm of a matrix A given by:
        Tr(sqrt{A^dagger * A})
        From: https://github.com/pyGSTio/pyGSTi/blob/master/pygsti/tools/optools.py
    """
    if np.linalg.norm(A - np.conjugate(A.T)) < tol:
        #Hermitian, so just sum eigenvalue magnitudes
        return np.sum(np.abs(np.linalg.eigvals(A)))
    else:
        #Sum of singular values (positive by construction)
        return np.sum(np.linalg.svd(A, compute_uv=False))

def tracedist(A, B, tol=np.sqrt(_eps)):
    """Compute the trace distance between matrices A and B given by:
        0.5 * Tr(sqrt{(A-B)^dagger * (A-B)})
        From: https://github.com/pyGSTio/pyGSTi/blob/master/pygsti/tools/optools.py
    """
    return 0.5 * tracenorm(A - B)

#### END FUNCTIONS COPIED FROM PYGSTI

def is_close(A, B, tol=np.sqrt(_eps)):
    """Check if two matrices are close in the sense of trace distance. 
    """
    if tracedist(A, B) < tol:
        return True 
    else:
        A /= np.exp(1j*np.angle(A[0,0]))
        B /= np.exp(1j*np.angle(B[0,0]))
        return tracedist(A, B) < tol

def haar_unitary(d):
    """Generate a Haar-random unitary matrix of dimension d.
       Algorithm from:
         F. Medrazzi. "How to generate random matrices from the classical compact groups"
           arXiv: math-ph/0609050
    """
    assert d > 1, 'Dimension must be > 1!'
    re_Z = np.random.randn(d*d).reshape((d,d))
    im_Z = np.random.randn(d*d).reshape((d,d))
    Z = (re_Z + 1j*im_Z)/np.sqrt(2.0)
    Q, R = np.linalg.qr(Z)
    L =  np.diag(np.diag(R) / np.abs(np.diag(R)))
    return Q @ L @ Q

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

def xyx_angles(U):
    """Euler angles for a unitary matrix U in the sequence X-Y-X.
        Note that angles are returned in matrix multiplication, not circuit order.
        We make use of the identity:
            Rx(a)Ry(b)Rx(c) = H Rz(a) Ry(-b) Rz(c) H
    """
    H = np.array([[1., 1.], [1., -1.]], dtype=np.complex128)/np.sqrt(2)
    ϕ, θ, λ = zyz_angles(H@U@H)
    return (ϕ, -1.0*θ, λ)

def diatomic_angles(U):
    ϕ, θ, λ = zyz_angles(U)
    a = ϕ
    b = np.pi - θ
    c = λ - np.pi
    return (a, b, c)

# C1 = {}
# C1[0] = pI
# C1[1] = expm(-1j * (pi / 4) * pX)
# C1[2] = expm(-2j * (pi / 4) * pX)
# C1[3] = expm(-3j * (pi / 4) * pX)
# C1[4] = expm(-1j * (pi / 4) * pY)
# C1[5] = expm(-2j * (pi / 4) * pY)
# C1[6] = expm(-3j * (pi / 4) * pY)
# C1[7] = expm(-1j * (pi / 4) * pZ)
# C1[8] = expm(-2j * (pi / 4) * pZ)
# C1[9] = expm(-3j * (pi / 4) * pZ)
# C1[10] = expm(-1j * (pi / 2) * (1 / np.sqrt(2)) * (pX + pY))
# C1[11] = expm(-1j * (pi / 2) * (1 / np.sqrt(2)) * (pX - pY))
# C1[12] = expm(-1j * (pi / 2) * (1 / np.sqrt(2)) * (pX + pZ))
# C1[13] = expm(-1j * (pi / 2) * (1 / np.sqrt(2)) * (pX - pZ))
# C1[14] = expm(-1j * (pi / 2) * (1 / np.sqrt(2)) * (pY + pZ))
# C1[15] = expm(-1j * (pi / 2) * (1 / np.sqrt(2)) * (pY - pZ))
# C1[16] = expm(-1j * (pi / 3) * (1 / np.sqrt(3)) * (pX + pY + pZ))
# C1[17] = expm(-2j * (pi / 3) * (1 / np.sqrt(3)) * (pX + pY + pZ))
# C1[18] = expm(-1j * (pi / 3) * (1 / np.sqrt(3)) * (pX - pY + pZ))
# C1[19] = expm(-2j * (pi / 3) * (1 / np.sqrt(3)) * (pX - pY + pZ))
# C1[20] = expm(-1j * (pi / 3) * (1 / np.sqrt(3)) * (pX + pY - pZ))
# C1[21] = expm(-2j * (pi / 3) * (1 / np.sqrt(3)) * (pX + pY - pZ))
# C1[22] = expm(-1j * (pi / 3) * (1 / np.sqrt(3)) * (-pX + pY + pZ))
# C1[23] = expm(-2j * (pi / 3) * (1 / np.sqrt(3)) * (-pX + pY + pZ))

def XYXClifford(qubit, cliff_num):
    """
    The set of 24 Diatomic Clifford single qubit pulses. Each pulse is decomposed
    as Rx(α)Ry(β)Rx(γ).

    Parameters
    ----------
    qubit : logical channel to implement sequence (LogicalChannel)
    cliffNum : the zero-indexed Clifford number

    Returns
    -------
    pulse object
    """
    α, β, γ = xyx_angles(C1[cliff_num])

    p1 =  Id(qubit) if np.isclose(γ, 0.0) else Xtheta(qubit, angle=γ)
    p2 =  Id(qubit) if np.isclose(β, 0.0) else Ytheta(qubit, angle=β)
    p3 =  Id(qubit) if np.isclose(α, 0.0) else Xtheta(qubit, angle=α)

    return p1 + p2 + p3

