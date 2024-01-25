import numpy as np
from typing import Tuple
import rotation_functions as rFunc

def MultQuat(q1:np.ndarray[4,1],q2:np.ndarray[4,1]) -> np.ndarray[4,1]:

    qs = q1[0,0]
    ps = q2[0,0]
    qv = q1[1:,0:]
    pv = q2[1:,0:]
    
    scalar = qs*ps - qv.T@pv

    vector = qs*pv + ps*qv + rFunc.sk(qv)@pv

    result = np.array([[scalar[0,0]],
                       [vector[0,0]],
                       [vector[1,0]],
                       [vector[2,0]]],dtype='float64')

    return result

def Normalized(q:np.ndarray[4,1]) -> np.ndarray[4,1]:
    ret = q.copy()
    ret = ret/np.linalg.norm(ret)
    return ret

def ConjQuat(q:np.ndarray[4,1]) -> np.ndarray[4,1]:
    qConj = q.copy()
    qConj[1:,0:] *= -1
    return qConj

def UnitQuat(q:np.ndarray[4,1]) -> np.ndarray[4,1]:
    return q.copy()/np.linalg.norm(q)

def Eaa2Quat(axis:np.ndarray[3,1],angle:float) -> np.ndarray[4,1]:
    sin = np.sin(angle)
    cos = np.cos(angle)
    u = axis.copy()
    return np.array([[cos],[sin*u[0,0]],[sin*u[1,0]],[sin*u[2,0]]])

def Rotate3D(v:np.ndarray[3,1], q:np.ndarray[4,1]) -> np.ndarray[3,1]:
    '''
    Rotates a vector in 3D space. quaternion must be normalized
    '''
    qv = np.array([[0],[v[0,0]],[v[1,0]],[v[2,0]]]) #convert the vector to quaternion for calculation
    v2 = MultQuat(MultQuat(q,qv),ConjQuat(q)) #v2 = q*v*~q    inner functions first, 
    return np.array(v2[1:,0:])


def Quat2RotM(q:np.ndarray[4,1]) -> np.ndarray[3,3]:

    q_copy = q.copy()
    qs = q_copy[0,0]
    qv = q_copy[1:,0:]

    #Referente al documento de teoria sobre conversion entre quaternion y matriz de rotacion:

    #Calcula por raices cuadradas. Usa el mayor valor para determinar el orden de calculos para evitar valores erroneos

    R1 = (qs**2-qv.T@qv)*np.eye(3)
    R2 = 2*(qv@qv.T)
    R3 = 2*qs*rFunc.sk(q[1:,:])

    R = R1+R2+R3

    return R

def DVec2Quat(v1:np.ndarray[3,1],v2:np.ndarray[3,1])->np.ndarray[4,1]:
    '''
    Returns a quaternion that represents the rotation between two vectors
    '''
    #First, normalize both vectors
    vec1 = v1.flatten()/np.linalg.norm(v1)
    vec2 = v2.flatten()/np.linalg.norm(v2)

    #Then we obtain the rotation axis between those two vectors by using cross product
    crossProd = np.cross(vec1,vec2)

    #Quaternion is formed by the dot product of both vectors +1 as the real part, and the cross product as the imaginary part
    q = np.ones((4,1))
    q[0,:] = 1 + vec1.dot(vec2)
    q[1:,0] = crossProd[:]

    #given two parallel vectors, an identity quaternion will be returned
    return Normalized(q)