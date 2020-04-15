"""
Tools for manipulating matrices

Original Author: Guilhem Ribeill

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
        A[np.abs(A) < tol] = 0.0
        B[np.abs(B) < tol] = 0.0
        A /= np.exp(1j*np.angle(A[0,0]))
        B /= np.exp(1j*np.angle(B[0,0]))
        return ((tracedist(A, B) < tol) or (tracedist(A, -1.0*B) < tol))

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