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
    return q.copy()/np.linalg.norm(q)

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

def Rotate3D(v:np.ndarray[3,1], q:np.ndarray[4,1]) -> Tuple[np.ndarray[3,1],float]:
    '''
    Rotates a vector in 3D space. quaternion must be normalized
    '''
    qv = np.array([[0],[v[0,0]],[v[1,0]],[v[2,0]]]) #convert the vector to quaternion for calculation
    v2 = MultQuat(q,MultQuat(qv,ConjQuat(q))) #v2 = q*v*~q
    return np.array(v2[1:,0:]), v2[0,0]


def Quat2RotM(q:np.ndarray[4,1]) -> np.ndarray[3,3]:
    R = np.zeros((3,3))
    #Referente al documento de teoria sobre conversion entre quaternion y matriz de rotacion:
    #Calcula por raices cuadradas. Usa el mayor valor para determinar el orden de calculos para evitar valores erroneos
    

    return R