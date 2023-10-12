# ----------------------------------------------------------------------
# Newton 2nd Law Utilities
# Created: Wed Oct 11 2023 Harrison B. Prosper
# ----------------------------------------------------------------------
import os, sys
import numpy as np
# ----------------------------------------------------------------------
# The following functions can operate on an array of vectors
# ----------------------------------------------------------------------
# compute magnitudes of one or more vectors
def magnitude(v):
    return np.sqrt(np.sum(v*v, axis=-1))

# compute unit vectors from one or more vectors
def norm(v):
    magv = magnitude(v)
    magv = np.where(magv < 1.e-10, 1, magv) # handle zero-length vectors
    return v / magv[:, np.newaxis]
# ----------------------------------------------------------------------
# Given the masses `m` and positions `r` of the objects (planets, etc.), 
# compute the net gravitational force on each object by performing the 
# double loop over positions ourselves. 
# ----------------------------------------------------------------------
def compute_grav_forces1(m, r):
    
    # this list will contain the net forces acting on
    # each particle
    f = []
    
    # loop over the positions of the particles.
    # use the "enumerate" construct so that we 
    # have access to the index "i" of each position (i = 0, 1,..., n-1)
    
    for i, ri in enumerate(r):
        
        # "gi" will be the net gravitational field at the position
        # of particle i
        gi = np.zeros(3)
        
        # loop over the positions of all other particles
        for j, rj in enumerate(r):
            
            if i == j: continue  # skip self-interactions
            
            # by convention, we compute the vector that points
            # from particle "j" to particle "i"
            rij = ri - rj
            
            # find its magnitude, that is, the scalar distance
            # between particles "i" and "j"
            magrij = magnitude(rij)
            
            # sum the gravitational fields at the location of particle "i".
            gi += m[j] * rij / magrij**3
            
        gi = -G * gi
        
        # compute net force on particle "i"
        fnet = m[i] * gi
        
        # cache net force
        f.append(fnet)
    
    # convert list of forces to a numpy array of forces
    return np.array(f)
# ----------------------------------------------------------------------
# Given the masses m and positions r of the objects (planets, etc.), 
# compute the net gravitational force on each object, this time using 
# the numpy broadcasting mechanism. 
# ----------------------------------------------------------------------
def compute_grav_forces2(m, r):
    
    # initial shapes of arrays m (masses) and r (position vectors)
    # m.shape: (n, )
    # r.shape: (n, 3)
    
    ri = r[np.newaxis, :]    # change shape from (n, 3) to (1, n, 3)
    rj = r[:, np.newaxis]    # change shape from (n, 3) to (n, 1, 3)
    
    # compute all possible vector differences using broadcasting:
    # 1. the row (1, n, 3) with n columns, each element of which 
    #    is a 3-vector, is replicated vertically n times to form 
    #    an array of shape (n, n, 3)
    #
    # 2. the column (n, 1, 3) with n rows, each element of which 
    #    is a 3-vector, is replicated horizontally n times to form
    #    an array of shape (n, n, 3).
    #
    # 3. now that we have arrays of the same shape, we can perform 
    #    element-by-element subtractions.
    #
    # note: the vectors along the diagonal are all zero-length vectors (see above)
    rij= ri - rj            # shape: (n, n, 3)
    
    # magrij has shape (n, n) 
    magrij = magnitude(rij)
    
    # for the magnitudes, replace the zeros along the diagonal
    # with ones. why must we do this in this function but not
    # in compute_forces1?
    np.fill_diagonal(magrij, 1)
        
    # why must we change the shape of the magnitudes array?
    magrij = magrij[:, :, np.newaxis]  # change shape to (n, n, 1)
    
    # why must we change the shape of the mass array?
    mj = m[:, np.newaxis, np.newaxis]  # change shape from (n, ) to (n, 1, 1)

    # for each particle, compute all gravitational fields
    gi = mj * rij / magrij**3         # shape of f: (n, n, 3)
    
    # compute net gravitational fields by summing the fields along axis 0, that is,
    # "vertically".
    gi = -G * gi.sum(axis=0)
    
    # compute the net force acting on each particle 
    mi = m[:, np.newaxis]
    fnet = mi * gi # shape: (n, 3)
    
    return fnet
