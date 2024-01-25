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


def Eaa2rotM(angle, axis):
    '''
    Returns the rotation matrix R able to rotate vectors an angle 'angle' (in rads) about the axis 'axis'
    Axis = X Y Z
    '''    
    axis_norm = np.linalg.norm(axis)
    

    if axis_norm > 1:
        axis = axis / axis_norm

    if axis.ndim == 1:
        axis = axis.reshape((-1, 1))
    

    R = np.eye(3) * np.cos(np.radians(angle)) + (1 - np.cos(np.radians(angle))) * np.outer(axis, axis) + np.sin(np.radians(angle)) * np.array([[0, -axis[2, 0], axis[1, 0]], [axis[2, 0], 0, -axis[0, 0]], [-axis[1, 0], axis[0, 0], 0]])

    return R

def RotM2Eaa(R:np.ndarray) ->Tuple[np.ndarray,float]:
    '''
    WIP: Converts a rotation matrix to an axis and angle tuple
    '''
    angle = np.arccos((np.trace(R)-1)/2)
    if angle == 0:
        axis = np.atleast_2d([1,1,1]).T
        axis = axis / np.linalg.norm(axis)
        return axis, angle

    u_sk = (R-R.T)/(2*np.sin(angle))
    axis = np.atleast_2d([u_sk[1,2],u_sk[2,0],u_sk[0,1]]).T

    return axis, angle


def RotVec2RotM(vector:np.ndarray) -> np.ndarray:
    '''
    Returns the rotation matrix R able to rotate vectors an angle 'angle' (in rads) about the axis 'axis'
    '''
    angle = np.linalg.norm(vector)
    axis = vector / angle

    sin = np.sin(angle)
    cos = np.cos(angle)
    R1 = np.eye(3)*cos
    R2 = (1-cos)*(axis @ axis.T)
    R3 = sk(axis)*sin

    R = R1+R2+R3

    return R

def Rotate3D_RV(v:np.ndarray[3,1], vector:np.ndarray) -> np.ndarray[3,1]:
    '''
    Rotates a vector in 3D space. quaternion must be normalized
    '''
  
    q = RotVec2RotM(vector)
    return q @ v


def Rotate3D_AA(v:np.ndarray[3,1], axis:np.ndarray[3,1], angle:float) -> np.ndarray[3,1]:
    '''
    Rotates a vector in 3D space. quaternion must be normalized
    '''
    q = Eaa2rotM(angle,axis)
    return q@v

#Rotation vector to axis angle
def RV2AA(vector:np.ndarray) -> Tuple[np.ndarray,float]:
    '''
    Converts a rotation vector to an axis and angle tuple
    '''
    angle = np.linalg.norm(vector)
    axis = vector / angle

    return axis, angle




def euler2rotM(angles):
    '''
    Returns the 3x3 rotation matrix R for a given set of Euler angles (in radians).
    Euler angles are specified as [phi, theta, psi], representing rotations around Z, Y, and X axes, respectively.
    '''

    phi, theta, psi = angles

    cos_phi, sin_phi = np.cos(phi), np.sin(phi)
    cos_theta, sin_theta = np.cos(theta), np.sin(theta)
    cos_psi, sin_psi = np.cos(psi), np.sin(psi)

    R = np.array([[cos_theta*cos_psi, -cos_phi*sin_psi + sin_phi*sin_theta*cos_psi, sin_phi*sin_psi + cos_phi*sin_theta*cos_psi],
                  [cos_theta*sin_psi, cos_phi*cos_psi + sin_phi*sin_theta*sin_psi, -sin_phi*cos_psi + cos_phi*sin_theta*sin_psi],
                  [-sin_theta, sin_phi*cos_theta, cos_phi*cos_theta]])

    return R

