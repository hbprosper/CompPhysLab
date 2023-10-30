# ----------------------------------------------------------------------
# Vector utilities based on numpy
# Created: Wed Oct 25 2023 Harrison B. Prosper
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
    magv = np.sqrt(np.sum(v*v, axis=-1))
    magv = np.where(magv < 1.e-15, 1, magv) # handle zero-length vectors
    try: 
        u = v / magv[:, np.newaxis]
    except:
        u = v / magv
    return u

def dot(a, b):
    '''
    
    Return the dot product of vectors a and b or an array of vectors a and b, where
    the vectors are modeled as numpy arrays.
    
    ab = dot(a, b)
    
    a: a vector or array of vectors
    b: a  
    '''
    try:
        c = (a*b).sum(axis=1)
    except:
        c = np.dot(a, b)
    return c

def tangent(u, n):
    '''
    Given the incident unit vector u and normal unit vector n, return the unit 
    vector n x u x n / |u x n| that is at right angles to n and it lies 
    in the plane defined by u and n.
    
    nt = tangent(u, n)
    
    u: a vector or array of vectors
    n: a vector of array of vectors
    '''
    return norm(np.cross(n, np.cross(u, n)))

def reflection(u, n):
    '''
    Given the incident unit vector u and normal unit vector n that defines the orientation
    of the boundary between two media, return the unit vector corresponding to a reflection 
    at the boundary. 
    
    ur = reflection(u, n)
    
    u: a vector or array of vectors
    n: a vector of array of vectors
    '''
    udotn = dot(u, n)
    try:
        udotn = udotn[:, np.newaxis] # why do we try this?
    except:
        pass
    return u - 2*udotn*n

def transmission(u, n, n1, n2):
    '''
    Given the incident unit vector u, normal unit vector n that defines the orientation
    of the boundary between two media of refractive indices n1 and n2, 
    return the unit vector corresponding to a transmission at the boundary.  
    
    ut = transmission(u, n, n1, n2)
    
    u: a vector or array of vectors
    n: a vector of array of vectors
    '''
    
    # n x u x n
    nun   = np.cross(n, np.cross(u, n))
    
    udotn = dot(u, n)
    
    n12   = n1/n2
    scale = np.sign(udotn)*np.sqrt(1-n12**2*(1-udotn**2))
    
    try:
        scale = scale[:, np.newaxis] # again, why try this?
    except:
        pass
 
    try:
        n12 = n12[:, np.newaxis]
    except:
        pass
    
    return scale * n + n12 * nun
