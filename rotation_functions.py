from typing import Tuple
import numpy as np

def rotation_angles(matrix, order):
    """
    input
        matrix = 3x3 rotation matrix (numpy array)
        oreder(str) = rotation order of x, y, z : e.g, rotation XZY -- 'xzy'
    output
        theta1, theta2, theta3 = rotation angles in rotation order
    """
    r11, r12, r13 = matrix[0]
    r21, r22, r23 = matrix[1]
    r31, r32, r33 = matrix[2]

    if order == 'xzx':
        theta1 = np.arctan(r31 / r21)
        theta2 = np.arctan(r21 / (r11 * np.cos(theta1)))
        theta3 = np.arctan(-r13 / r12)

    elif order == 'xyx':
        theta1 = np.arctan(-r21 / r31)
        theta2 = np.arctan(-r31 / (r11 *np.cos(theta1)))
        theta3 = np.arctan(r12 / r13)

    elif order == 'yxy':
        theta1 = np.arctan(r12 / r32)
        theta2 = np.arctan(r32 / (r22 *np.cos(theta1)))
        theta3 = np.arctan(-r21 / r23)

    elif order == 'yzy':
        theta1 = np.arctan(-r32 / r12)
        theta2 = np.arctan(-r12 / (r22 *np.cos(theta1)))
        theta3 = np.arctan(r23 / r21)

    elif order == 'zyz':
        theta1 = np.arctan(r23 / r13)
        theta2 = np.arctan(r13 / (r33 *np.cos(theta1)))
        theta3 = np.arctan(-r32 / r31)

    elif order == 'zxz':
        theta1 = np.arctan(-r13 / r23)
        theta2 = np.arctan(-r23 / (r33 *np.cos(theta1)))
        theta3 = np.arctan(r31 / r32)

    elif order == 'xzy':
        theta1 = np.arctan(r32 / r22)
        theta2 = np.arctan(-r12 * np.cos(theta1) / r22)
        theta3 = np.arctan(r13 / r11)

    elif order == 'xyz':
        theta1 = np.arctan(-r23 / r33)
        theta2 = np.arctan(r13 * np.cos(theta1) / r33)
        theta3 = np.arctan(-r12 / r11)

    elif order == 'yxz':
        theta1 = np.arctan(r13 / r33)
        theta2 = np.arctan(-r23 * np.cos(theta1) / r33)
        theta3 = np.arctan(r21 / r22)

    elif order == 'yzx':
        theta1 = np.arctan(-r31 / r11)
        theta2 = np.arctan(r21 * np.cos(theta1) / r11)
        theta3 = np.arctan(-r23 / r22)

    elif order == 'zyx':
        theta1 = np.arctan(r21 / r11)
        theta2 = np.arctan(-r31 * np.cos(theta1) / r11)
        theta3 = np.arctan(r32 / r33)

    elif order == 'zxy':
        theta1 = np.arctan(-r12 / r22)
        theta2 = np.arctan(r32 * np.cos(theta1) / r22)
        theta3 = np.arctan(-r31 / r33)

    theta1 = theta1 * 180 / np.pi
    theta2 = theta2 * 180 / np.pi
    theta3 = theta3 * 180 / np.pi

    return (theta1, theta2, theta3)



def rotation_angles_with_zeros_zyx(matrix:np.ndarray) -> Tuple[float, float, float]:
    """
      input
          matrix = 3x3 rotation matrix (numpy array)
      output
          psi, theta, phi = rotation angles in rotation order
      """
    r11, r12, r13 = matrix[0]
    r21, r22, r23 = matrix[1]
    r31, r32, r33 = matrix[2]

    psi, theta, phi = 0, 0, 0

    if (abs(r31) == 1):
      if (r31 == -1):
        theta = np.pi/2
        psi = phi + np.arctan2(r12,r13)
      else:
        theta = -np.pi/2
        psi = -phi + np.arctan2(-r12,-r13)

    psi = psi * 180 / np.pi
    theta = theta * 180 / np.pi
    phi = phi * 180 / np.pi


    return (psi,theta,phi)



def sk(v:np.ndarray[3,1]) -> np.ndarray[3,3]:
    '''
    returns skew-symmetrical matrix of vector cross product with itself
    devuelve matriz antisimetrica del producto vectorial del vector con si mismo
    '''
    v = v.flatten()
    vsk = np.array([[0,  -v[2], v[1]],
                    [v[2],  0, -v[0]],
                    [-v[1], v[0], 0]])

    return vsk


def Eaa2rotM(angle: float, axis: np.ndarray[3,1]) -> np.ndarray[3,3]:
    '''
    Returns the rotation matrix R able to rotate vectors an angle 'angle' (in rads) about the axis 'axis'
    '''
    axis = axis / np.linalg.norm(axis)

    sin = np.sin(angle)
    cos = np.cos(angle)
    R1 = np.eye(3)*cos
    R2 = (1-cos)*(axis @ axis.T)
    R3 = sk(axis)*sin

    R = R1+R2+R3

    return R


def RotM2Eaa(matrix:np.ndarray) ->Tuple[np.ndarray,float]:
    '''
    WIP: Converts a rotation matrix to an axis and angle tuple
    '''
    axis = 0
    angle = 0


    return axis, angle












