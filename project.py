import customtkinter
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
# Implement the default Matplotlib key bindings.

from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib.pyplot as plt

import numpy as np

import quaternion_functions as quatFunc
import rotation_functions as rotFunc

customtkinter.set_appearance_mode("Dark")  # Modes: "System" (standard), "Dark", "Light"
customtkinter.set_default_color_theme("dark-blue")  # Themes: "blue" (standard), "green", "dark-blue"


class Arcball(customtkinter.CTk):

    def __init__(self):
        super().__init__()

        # Mouse position from previous drag frame
        self.mouseX = 0
        self.mouseY = 0

        # Orientation vars. Initialized to represent 0 rotation
        self.quat = np.array([[1.0],[0.0],[0.0],[0.0]])
        self.rotM = np.eye(3)
        self.AA = {"axis": np.array([[0.0],[0.0],[0.0]]), "angle":0.0}
        self.rotv = np.array([[0.0],[0.0],[0.0]])
        self.euler = np.array([[0.0],[0.0],[0.0]])

        # configure window
        self.title("Holroyd's arcball")
        self.geometry(f"{1100}x{580}")
        self.resizable(False, False)

        self.grid_columnconfigure((0,1), weight=0   )
        self.grid_rowconfigure((0,1), weight=1)
        self.grid_rowconfigure(2, weight=0)

        # Cube plot
        self.init_cube()

        self.canvas = FigureCanvasTkAgg(self.fig, self)  # A tk.DrawingArea.
        self.bm = BlitManager(self.canvas,[self.facesObj])
        self.canvas.draw()
        self.canvas.get_tk_widget().grid(row=0, column=0, rowspan=2, padx=(20, 20), pady=(20, 20), sticky="nsew")

        self.pressed = False #Bool to bypass the information that mouse is clicked
        self.canvas.mpl_connect('button_press_event', self.onclick)
        self.canvas.mpl_connect('motion_notify_event', self.onmove)
        self.canvas.mpl_connect('button_release_event', self.onrelease)
        
        # Reset button
        self.resetbutton = customtkinter.CTkButton(self, text="Reset", command=self.resetbutton_pressed)
        self.resetbutton.grid(row=3, column=0, padx=(0, 0), pady=(5, 20), sticky="ns")
        
        # Selectable atti
        self.tabview = customtkinter.CTkTabview(self, width=150, height=150)
        self.tabview.grid(row=0, column=1, padx=(0, 20), pady=(20, 0), sticky="nsew")
        self.tabview.add("Axis angle")
        self.tabview.add("Rotation vector")
        self.tabview.add("Euler angles")
        self.tabview.add("Quaternion")

        # Selectable atti: AA
        self.tabview.tab("Axis angle").grid_columnconfigure(0, weight=0)  # configure grid of individual tabs
        self.tabview.tab("Axis angle").grid_columnconfigure(1, weight=0)  # configure grid of individual tabs

        self.label_AA_axis= customtkinter.CTkLabel(self.tabview.tab("Axis angle"), text="Axis:")
        self.label_AA_axis.grid(row=0, column=0, rowspan=3, padx=(80,0), pady=(45,0), sticky="e")

        self.entry_AA_ax1 = customtkinter.CTkEntry(self.tabview.tab("Axis angle"))
        self.entry_AA_ax1.insert(0,"1.0")
        self.entry_AA_ax1.grid(row=0, column=1, padx=(5, 0), pady=(50, 0), sticky="ew")

        self.entry_AA_ax2 = customtkinter.CTkEntry(self.tabview.tab("Axis angle"))
        self.entry_AA_ax2.insert(0,"0.0")
        self.entry_AA_ax2.grid(row=1, column=1, padx=(5, 0), pady=(5, 0), sticky="ew")

        self.entry_AA_ax3 = customtkinter.CTkEntry(self.tabview.tab("Axis angle"))
        self.entry_AA_ax3.insert(0,"0.0")
        self.entry_AA_ax3.grid(row=2, column=1, padx=(5, 0), pady=(5, 10), sticky="ew")

        self.label_AA_angle = customtkinter.CTkLabel(self.tabview.tab("Axis angle"), text="Angle:")
        self.label_AA_angle.grid(row=3, column=0, padx=(120,0), pady=(10, 20),sticky="w")
        self.entry_AA_angle = customtkinter.CTkEntry(self.tabview.tab("Axis angle"))
        self.entry_AA_angle.insert(0,"0.0")
        self.entry_AA_angle.grid(row=3, column=1, padx=(5, 0), pady=(0, 10), sticky="ew")

        self.button_AA = customtkinter.CTkButton(self.tabview.tab("Axis angle"), text="Apply", command=self.apply_AA, width=180)
        self.button_AA.grid(row=5, column=0, columnspan=2, padx=(0, 0), pady=(5, 0), sticky="e")

        # Selectable atti: rotV
        self.tabview.tab("Rotation vector").grid_columnconfigure(0, weight=1)
        self.tabview.tab("Rotation vector").grid_columnconfigure(1, weight=0)
        
        self.label_rotV= customtkinter.CTkLabel(self.tabview.tab("Rotation vector"), text="rot. Vector:")
        self.label_rotV.grid(row=0, column=0, rowspan=3, padx=(2,0), pady=(45,0), sticky="e")

        self.entry_rotV_1 = customtkinter.CTkEntry(self.tabview.tab("Rotation vector"))
        self.entry_rotV_1.insert(0,"0.0")
        self.entry_rotV_1.grid(row=0, column=1, padx=(5, 60), pady=(50, 0), sticky="ew")

        self.entry_rotV_2 = customtkinter.CTkEntry(self.tabview.tab("Rotation vector"))
        self.entry_rotV_2.insert(0,"0.0")
        self.entry_rotV_2.grid(row=1, column=1, padx=(5, 60), pady=(5, 0), sticky="ew")

        self.entry_rotV_3 = customtkinter.CTkEntry(self.tabview.tab("Rotation vector"))
        self.entry_rotV_3.insert(0,"0.0")
        self.entry_rotV_3.grid(row=2, column=1, padx=(5, 60), pady=(5, 10), sticky="ew")

        self.button_rotV = customtkinter.CTkButton(self.tabview.tab("Rotation vector"), text="Apply", command=self.apply_rotV, width=180)
        self.button_rotV.grid(row=5, column=0, columnspan=2, padx=(0, 60), pady=(5, 0), sticky="e")

        # Selectable atti: Euler angles
        self.tabview.tab("Euler angles").grid_columnconfigure(0, weight=1)
        self.tabview.tab("Euler angles").grid_columnconfigure(1, weight=0)
        
        self.label_EA_roll= customtkinter.CTkLabel(self.tabview.tab("Euler angles"), text="roll:")
        self.label_EA_roll.grid(row=0, column=0, padx=(2,0), pady=(50,0), sticky="e")

        self.label_EA_pitch= customtkinter.CTkLabel(self.tabview.tab("Euler angles"), text="pitch:")
        self.label_EA_pitch.grid(row=1, column=0, padx=(2,0), pady=(5,0), sticky="e")

        self.label_EA_yaw= customtkinter.CTkLabel(self.tabview.tab("Euler angles"), text="yaw:")
        self.label_EA_yaw.grid(row=2, column=0, rowspan=3, padx=(2,0), pady=(5,10), sticky="e")

        self.entry_EA_roll = customtkinter.CTkEntry(self.tabview.tab("Euler angles"))
        self.entry_EA_roll.insert(0,"0.0")
        self.entry_EA_roll.grid(row=0, column=1, padx=(5, 60), pady=(50, 0), sticky="ew")

        self.entry_EA_pitch = customtkinter.CTkEntry(self.tabview.tab("Euler angles"))
        self.entry_EA_pitch.insert(0,"0.0")
        self.entry_EA_pitch.grid(row=1, column=1, padx=(5, 60), pady=(5, 0), sticky="ew")

        self.entry_EA_yaw = customtkinter.CTkEntry(self.tabview.tab("Euler angles"))
        self.entry_EA_yaw.insert(0,"0.0")
        self.entry_EA_yaw.grid(row=2, column=1, padx=(5, 60), pady=(5, 10), sticky="ew")

        self.button_EA = customtkinter.CTkButton(self.tabview.tab("Euler angles"), text="Apply", command=self.apply_EA, width=180)
        self.button_EA.grid(row=5, column=0, columnspan=2, padx=(0, 60), pady=(5, 0), sticky="e")

        # Selectable atti: Quaternion
        self.tabview.tab("Quaternion").grid_columnconfigure(0, weight=1)
        self.tabview.tab("Quaternion").grid_columnconfigure(1, weight=0)
        
        self.label_quat_0= customtkinter.CTkLabel(self.tabview.tab("Quaternion"), text="q0:")
        self.label_quat_0.grid(row=0, column=0, padx=(2,0), pady=(50,0), sticky="e")

        self.label_quat_1= customtkinter.CTkLabel(self.tabview.tab("Quaternion"), text="q1:")
        self.label_quat_1.grid(row=1, column=0, padx=(2,0), pady=(5,0), sticky="e")

        self.label_quat_2= customtkinter.CTkLabel(self.tabview.tab("Quaternion"), text="q2:")
        self.label_quat_2.grid(row=2, column=0, padx=(2,0), pady=(5,0), sticky="e")

        self.label_quat_3= customtkinter.CTkLabel(self.tabview.tab("Quaternion"), text="q3:")
        self.label_quat_3.grid(row=3, column=0, padx=(2,0), pady=(5,10), sticky="e")

        self.entry_quat_0 = customtkinter.CTkEntry(self.tabview.tab("Quaternion"))
        self.entry_quat_0.insert(0,"1.0")
        self.entry_quat_0.grid(row=0, column=1, padx=(5, 60), pady=(50, 0), sticky="ew")

        self.entry_quat_1 = customtkinter.CTkEntry(self.tabview.tab("Quaternion"))
        self.entry_quat_1.insert(0,"0.0")
        self.entry_quat_1.grid(row=1, column=1, padx=(5, 60), pady=(5, 0), sticky="ew")

        self.entry_quat_2 = customtkinter.CTkEntry(self.tabview.tab("Quaternion"))
        self.entry_quat_2.insert(0,"0.0")
        self.entry_quat_2.grid(row=2, column=1, padx=(5, 60), pady=(5, 0), sticky="ew")

        self.entry_quat_3 = customtkinter.CTkEntry(self.tabview.tab("Quaternion"))
        self.entry_quat_3.insert(0,"0.0")
        self.entry_quat_3.grid(row=3, column=1, padx=(5, 60), pady=(5, 10), sticky="ew")

        self.button_quat = customtkinter.CTkButton(self.tabview.tab("Quaternion"), text="Apply", command=self.apply_quat, width=180)
        self.button_quat.grid(row=4, column=0, columnspan=2, padx=(0, 60), pady=(5, 0), sticky="e")

        # Rotation matrix info
        self.RotMFrame = customtkinter.CTkFrame(self, width=150)
        self.RotMFrame.grid(row=1, column=1, rowspan=3, padx=(0, 20), pady=(20, 20), sticky="nsew")

        self.RotMFrame.grid_columnconfigure((0,1,2,3,4), weight=1)

        self.label_RotM= customtkinter.CTkLabel(self.RotMFrame, text="RotM = ")
        self.label_RotM.grid(row=0, column=0, rowspan=3, padx=(2,0), pady=(20,0), sticky="e")

        self.entry_RotM_11= customtkinter.CTkEntry(self.RotMFrame, width=50, border_width=0)
        self.entry_RotM_11.insert(0,"1.0")
        self.entry_RotM_11.configure(state="disabled")
        self.entry_RotM_11.grid(row=0, column=1, padx=(2,0), pady=(20,0), sticky="ew")

        self.entry_RotM_12= customtkinter.CTkEntry(self.RotMFrame, width=50, border_width=0)
        self.entry_RotM_12.insert(0,"0.0")
        self.entry_RotM_12.configure(state="disabled")
        self.entry_RotM_12.grid(row=0, column=2, padx=(2,0), pady=(20,0), sticky="ew")

        self.entry_RotM_13= customtkinter.CTkEntry(self.RotMFrame, width=50, border_width=0)
        self.entry_RotM_13.insert(0,"0.0")
        self.entry_RotM_13.configure(state="disabled")
        self.entry_RotM_13.grid(row=0, column=3, padx=(2,0), pady=(20,0), sticky="ew")

        self.entry_RotM_21= customtkinter.CTkEntry(self.RotMFrame, width=50, border_width=0)
        self.entry_RotM_21.insert(0,"0.0")
        self.entry_RotM_21.configure(state="disabled")
        self.entry_RotM_21.grid(row=1, column=1, padx=(2,0), pady=(2,0), sticky="ew")

        self.entry_RotM_22= customtkinter.CTkEntry(self.RotMFrame, width=50, border_width=0)
        self.entry_RotM_22.insert(0,"1.0")
        self.entry_RotM_22.configure(state="disabled")
        self.entry_RotM_22.grid(row=1, column=2, padx=(2,0), pady=(2,0), sticky="ew")

        self.entry_RotM_23= customtkinter.CTkEntry(self.RotMFrame, width=50, border_width=0)
        self.entry_RotM_23.insert(0,"0.0")
        self.entry_RotM_23.configure(state="disabled")
        self.entry_RotM_23.grid(row=1, column=3, padx=(2,0), pady=(2,0), sticky="ew")

        self.entry_RotM_31= customtkinter.CTkEntry(self.RotMFrame, width=50, border_width=0)
        self.entry_RotM_31.insert(0,"0.0")
        self.entry_RotM_31.configure(state="disabled")
        self.entry_RotM_31.grid(row=2, column=1, padx=(2,0), pady=(2,0), sticky="ew")

        self.entry_RotM_32= customtkinter.CTkEntry(self.RotMFrame, width=50, border_width=0)
        self.entry_RotM_32.insert(0,"0.0")
        self.entry_RotM_32.configure(state="disabled")
        self.entry_RotM_32.grid(row=2, column=2, padx=(2,0), pady=(2,0), sticky="ew")

        self.entry_RotM_33= customtkinter.CTkEntry(self.RotMFrame, width=50, border_width=0)
        self.entry_RotM_33.insert(0,"1.0")
        self.entry_RotM_33.configure(state="disabled")
        self.entry_RotM_33.grid(row=2, column=3, padx=(2,0), pady=(2,0), sticky="ew")
    
    def setRotMatrix(self, rotM):
        """
        Funcion que hace cambiar la rotM de abajo
        """

        rotM = rotM.copy()
       
        self.entry_RotM_11.configure(state="normal")
        self.entry_RotM_12.configure(state="normal")
        self.entry_RotM_13.configure(state="normal")
        self.entry_RotM_21.configure(state="normal")
        self.entry_RotM_22.configure(state="normal")
        self.entry_RotM_23.configure(state="normal")
        self.entry_RotM_31.configure(state="normal")
        self.entry_RotM_32.configure(state="normal")
        self.entry_RotM_33.configure(state="normal")

        self.entry_RotM_11.delete(0, 99)
        self.entry_RotM_12.delete(0, 99)
        self.entry_RotM_13.delete(0, 99)
        self.entry_RotM_21.delete(0, 99)
        self.entry_RotM_22.delete(0, 99)
        self.entry_RotM_23.delete(0, 99)
        self.entry_RotM_31.delete(0, 99)
        self.entry_RotM_32.delete(0, 99)
        self.entry_RotM_33.delete(0, 99)

        
        if(abs(rotM[0][0]) < 1.0e-16):
            rotM[0][0] = 0
        if(abs(rotM[0][1]) < 1.0e-16):
            rotM[0][1] = 0
        if(abs(rotM[0][2]) < 1.0e-16):
            rotM[0][2] = 0
        if(abs(rotM[1][0]) < 1.0e-16):
            rotM[1][0] = 0
        if(abs(rotM[1][1]) < 1.0e-16):
            rotM[1][1] = 0
        if(abs(rotM[1][2]) < 1.0e-16):
            rotM[1][2] = 0
        if(abs(rotM[2][0]) < 1.0e-16):
            rotM[2][0] = 0
        if(abs(rotM[2][1]) < 1.0e-16):
            rotM[2][1] = 0
        if(abs(rotM[2][2]) < 1.0e-16):
            rotM[2][2] = 0
        

        

        self.entry_RotM_11.insert(0,rotM[0][0])
        self.entry_RotM_12.insert(0,rotM[0][1])
        self.entry_RotM_13.insert(0,rotM[0][2])
        self.entry_RotM_21.insert(0,rotM[1][0])
        self.entry_RotM_22.insert(0,rotM[1][1])
        self.entry_RotM_23.insert(0,rotM[1][2])
        self.entry_RotM_31.insert(0,rotM[2][0])
        self.entry_RotM_32.insert(0,rotM[2][1])
        self.entry_RotM_33.insert(0,rotM[2][2])

        self.entry_RotM_11.configure(state="disabled")
        self.entry_RotM_12.configure(state="disabled")
        self.entry_RotM_13.configure(state="disabled")
        self.entry_RotM_21.configure(state="disabled")
        self.entry_RotM_22.configure(state="disabled")
        self.entry_RotM_23.configure(state="disabled")
        self.entry_RotM_31.configure(state="disabled")
        self.entry_RotM_32.configure(state="disabled")
        self.entry_RotM_33.configure(state="disabled")

    def rotMToAngleAxis(self, rotM):

        
        R = rotM

        angle = np.arccos((np.trace(R) - 1) / 2)
        if(angle != 0):
            axis = (R - R.T) / (2 * np.sin(angle))
            axis_org = np.array([axis[2, 1], axis[0, 2], axis[1, 0]])
        else:
            axis_org = np.array([1,0,0])

        angle = np.rad2deg(angle)

        self.entry_AA_angle.delete(0,99)
        self.entry_AA_angle.insert(0, angle)

        self.entry_AA_ax1.delete(0,99)
        self.entry_AA_ax2.delete(0,99)
        self.entry_AA_ax3.delete(0,99)
        self.entry_AA_ax1.insert(0,axis_org[0])
        self.entry_AA_ax2.insert(0,axis_org[1])
        self.entry_AA_ax3.insert(0,axis_org[2])

    
    def rotMToRotationVector(self, rotM):
        R = rotM

        angle = np.arccos((np.trace(R) - 1) / 2)

        


        if(angle != 0):
            axis = (R - R.T) / (2 * np.sin(angle))
            axis_org = np.array([axis[2, 1], axis[0, 2], axis[1, 0]])
        else:
            axis_org = np.array([1,0,0])

        rotV = axis_org * np.rad2deg(angle)


       
        self.entry_rotV_1.delete(0,99)
        self.entry_rotV_2.delete(0,99)
        self.entry_rotV_3.delete(0,99)
        self.entry_rotV_1.insert(0,rotV[0])
        self.entry_rotV_2.insert(0,rotV[1])
        self.entry_rotV_3.insert(0,rotV[2])
        
     
    def resetbutton_pressed(self):
      """
      Event triggered function on the event of a push on the button Reset
      """
      # Restaurar la posición inicial del cubo
      self.M = np.array([
          [-1, -1, 1],
          [-1, 1, 1],
          [1, 1, 1],
          [1, -1, 1],
          [-1, -1, -1],
          [-1, 1, -1],
          [1, 1, -1],
          [1, -1, -1]
      ], dtype=float).transpose()

      # Actualizar el cubo en la interfaz gráfica
      self.update_cube()

    
    def apply_AA(self):
        """
        Event triggered function on the event of a push on the button button_AA
        """
        self.M = np.array(
            [[ -1,  -1, 1],   #Node 0
            [ -1,   1, 1],    #Node 1
            [1,   1, 1],      #Node 2
            [1,  -1, 1],      #Node 3
            [-1,  -1, -1],    #Node 4
            [-1,  1, -1],     #Node 5
            [1,   1, -1],     #Node 6
            [1,  -1, -1]], dtype=float).transpose()

        angle = float(self.entry_AA_angle.get())
        axisX = float(self.entry_AA_ax1.get())
        axisY = float(self.entry_AA_ax2.get())
        axisZ = float(self.entry_AA_ax3.get())
        rotM = rotFunc.Eaa2rotM(angle, np.array([axisX, axisY, axisZ]))
    
        #self.M = rotM@self.M 
        
        
        self.M = np.dot(rotM, self.M)
        


        self.setRotMatrix(rotM)
        self.rotMToRotationVector(rotM)
        self.rotM = rotM
        self.update_cube()

        pass
        
    
    def convert_to_axis_angle(rotv):
        """
        Convert a rotation vector to axis and angle representation.
        """
        angle = np.linalg.norm(rotv)
        axis = rotv / angle
        return axis, angle

    def apply_rotV(self):
        """
        Event triggered function on the event of a push on the button button_rotV 
        """
        rV = np.array([[float(self.entry_rotV_1.get())],[float(self.entry_rotV_2.get())],[float(self.entry_rotV_3.get())]])
        axisSacado = rV / np.linalg.norm(rV)
        angleSacado = np.linalg.norm(rV)

        rotM = rotFunc.Eaa2rotM(angleSacado, axisSacado)
        self.M = np.dot(rotM, self.M)
        self.setRotMatrix(rotM)
        self.rotMToAngleAxis(rotM)
        self.rotM = rotM
        self.update_cube()

        
        entries = [
            (self.entry_RotM_11, self.rotM[0, 0]),
            (self.entry_RotM_12, self.rotM[0, 1]),
            (self.entry_RotM_13, self.rotM[0, 2]),
            (self.entry_RotM_21, self.rotM[1, 0]),
            (self.entry_RotM_22, self.rotM[1, 1]),
            (self.entry_RotM_23, self.rotM[1, 2]),
            (self.entry_RotM_31, self.rotM[2, 0]),
            (self.entry_RotM_32, self.rotM[2, 1]),
            (self.entry_RotM_33, self.rotM[2, 2])
        ]

        for entry, value in entries:
            entry.configure(state="normal")
            entry.delete(0, "end")
            entry.insert(0, value)
            entry.configure(state="disabled")
                
                
    def apply_EA(self):
     """
     Event triggered function on the event of a push on the button button_euler_angles
     """
     angles = np.zeros(3)
     angles[0] = self.entry_EA_roll.get()
     angles[1] = self.entry_EA_pitch.get()
     angles[2] = self.entry_EA_yaw.get()

     for i in range(self.M.shape[1]): 
        v = np.array(self.M[:,i], ndmin=2).T
        
        yaw, pitch, roll = angles

        # Matrices de rotación en los ejes respectivos
        R_yaw = np.array([[np.cos(yaw), -np.sin(yaw), 0],
                          [np.sin(yaw), np.cos(yaw), 0],
                          [0, 0, 1]])

        R_pitch = np.array([[np.cos(pitch), 0, np.sin(pitch)],
                            [0, 1, 0],
                            [-np.sin(pitch), 0, np.cos(pitch)]])

        R_roll = np.array([[1, 0, 0],
                           [0, np.cos(roll), -np.sin(roll)],
                           [0, np.sin(roll), np.cos(roll)]])
    
        # Aplicar rotaciones en el orden yaw-pitch-roll
        rotated_vector = R_yaw @ R_pitch @ R_roll @ v

        self.M[:,i] = rotated_vector[:,0].T
 
     self.rotM = rotFunc.euler2rotM(angles)
    
     self.update_rotM_visual()

     self.update_cube()

    
    def apply_quat(self):
        """
        Event triggered function on the event of a push on the button button_quat
        """
        self.quat[0,0] = self.entry_quat_0.get()
        self.quat[1,0] = self.entry_quat_1.get()
        self.quat[2,0] = self.entry_quat_2.get()
        self.quat[3,0] = self.entry_quat_3.get()

        self.quat = self.quat/np.linalg.norm(self.quat) #normalize the quaternion
        
        for i in range(self.M.shape[1]):
            v = np.array(self.M[:,i],ndmin=2).T
            vr = quatFunc.Rotate3D(v, self.quat)
            self.M[:,i] = vr[:,0].T
        
        #update quaternion shown on screen
        self.entry_quat_0.delete(0,'end')
        self.entry_quat_1.delete(0,'end')
        self.entry_quat_2.delete(0,'end')
        self.entry_quat_3.delete(0,'end')
        self.entry_quat_0.insert(0,self.quat[0,0])
        self.entry_quat_1.insert(0,self.quat[1,0])
        self.entry_quat_2.insert(0,self.quat[2,0])
        self.entry_quat_3.insert(0,self.quat[3,0])

        #convert quaternion into rotation matrix to update the shown rotation matrix
        self.rotM = quatFunc.Quat2RotM(q)
        self.update_rotM_visual()
        
        self.update_cube()

    def update_rotM_visual(self):
        """
        Updates the rotation matrix shown on the bottom-right corner of the application
        """
        entries = [
            (self.entry_RotM_11,self.rotM[0,0]),
            (self.entry_RotM_12,self.rotM[0,1]),
            (self.entry_RotM_13,self.rotM[0,2]),
            (self.entry_RotM_21,self.rotM[1,0]),
            (self.entry_RotM_22,self.rotM[1,1]),
            (self.entry_RotM_23,self.rotM[1,2]),
            (self.entry_RotM_31,self.rotM[2,0]),
            (self.entry_RotM_32,self.rotM[2,1]),
            (self.entry_RotM_33,self.rotM[2,2]),
        ]

        for entry,value in entries:
            entry.configure(state="normal")
            entry.delete(0,"end")
            entry.insert(0,value)
            entry.configure(state="disabled")
    
    def onclick(self, event):
        """
        Event triggered function on the event of a mouse click inside the figure canvas
        """
        print("Pressed button", event.button)

        if event.button:
            self.pressed = True # Bool to control(activate) a drag (click+move)
            self.mouseX, self.mouseY = self.canvas_coordinates_to_figure_coordinates(event.x,event.y) #Set initial mouse drag coordinates

    def GetVectorFromSurface(self,x,y,r):
        """
        Returns a vector with Holroyd's arcball method, given the figure coordinates and a virtual sphere radius
        """
        ret = np.ones((3,1))

        if (((x**2) + (y**2))< 0.5*r**2): #calculate from sphere's surface
            ret[0,0] = x
            ret[1,0] = y
            ret[2,0] = np.sqrt((r**2)-(x**2)-(y**2))
        else:
            ret[0,0] = x
            ret[1,0] = y
            ret[2,0] = ((r**2)/(2*np.sqrt(x**2+y**2)))
        
        norm = np.linalg.norm(ret)
        ret = ret/norm #normalize the vector

        return ret

    def onmove(self,event):
        """
        Event triggered function on the event of a mouse motion
        """
        
        #Example
        if self.pressed: #Only triggered if previous click
            x_fig,y_fig= self.canvas_coordinates_to_figure_coordinates(event.x,event.y) #Extract viewport coordinates
            
            r2 = x_fig*x_fig+y_fig*y_fig

            #print("x: ", x_fig)
            #print("y", y_fig)
            #print("r2", r2)

            #Roger: no tengo ni idea de por que el cubo termina rotando tan rapido, pero esto es lo mas cerca que he podido llegar de hacer que gire correctamente

            movX,movY = x_fig,y_fig#x_fig-self.mouseX, y_fig-self.mouseY
            self.mouseX, self.mouseY = x_fig, y_fig

            prevAxis, prevAngle = rotFunc.RotM2Eaa(self.rotM)

            #use vertex distance from origin as virtual sphere radius
            sVec = self.GetVectorFromSurface(movX,movY,np.linalg.norm(self.M[:,0]))

            newQ = quatFunc.DVec2Quat(prevAxis,sVec)

            print("quat:",newQ.flatten())

            self.rotM = R = quatFunc.Quat2RotM(newQ)
            #R = np.array([[0, 1, 0], [-1, 0, 0], [0, 0, 1]])
                    
            self.M = R.dot(self.M) #Modify the vertices matrix with a rotation matrix M

            self.update_rotM_visual()
            self.update_cube() #Update the cube


    def onrelease(self,event):
        """
        Event triggered function on the event of a mouse release
        """
        self.pressed = False # Bool to control(deactivate) a drag (click+move)


    def init_cube(self):
        """
        Initialization function that sets up cube's geometry and plot information
        """

        self.M = np.array(
            [[ -1,  -1, 1],   #Node 0
            [ -1,   1, 1],    #Node 1
            [1,   1, 1],      #Node 2
            [1,  -1, 1],      #Node 3
            [-1,  -1, -1],    #Node 4
            [-1,  1, -1],     #Node 5
            [1,   1, -1],     #Node 6
            [1,  -1, -1]], dtype=float).transpose() #Node 7

        self.con = [
            [0, 1, 2, 3], #Face 1
            [4, 5, 6, 7], #Face 2
            [3, 2, 6, 7], #Face 3
            [0, 1, 5, 4], #Face 4
            [0, 3, 7, 4], #Face 5
            [1, 2, 6, 5]] #Face 6

        faces = []

        for row in self.con:
            faces.append([self.M[:,row[0]],self.M[:,row[1]],self.M[:,row[2]],self.M[:,row[3]]])

        self.fig = plt.figure()
        ax = self.fig.add_subplot(111, projection='3d')

        for item in [self.fig, ax]:
            item.patch.set_visible(False)

        self.facesObj = Poly3DCollection(faces, linewidths=.2, edgecolors='k',animated = True)
        self.facesObj.set_facecolor([(0,0,1,0.9), #Blue
        (0,1,0,0.9), #Green
        (.9,.5,0.13,0.9), #Orange
        (1,0,0,0.9), #Red
        (1,1,0,0.9), #Yellow
        (0,0,0,0.9)]) #Black

        #Transfering information to the plot
        ax.add_collection3d(self.facesObj)

        #Configuring the plot aspect
        ax.azim=-90
        ax.roll = -90
        ax.elev=0   
        ax.set_xlim3d(-2, 2)
        ax.set_ylim3d(-2, 2)
        ax.set_zlim3d(-2, 2)
        ax.set_aspect('equal')
        ax.disable_mouse_rotation()
        ax.set_axis_off()

        self.pix2unit = 1.0/60 #ratio for drawing the cube 


    def update_cube(self):
        """
        Updates the cube vertices and updates the figure.
        Call this function after modifying the vertex matrix in self.M to redraw the cube
        """

        faces = []

        for row in self.con:
            faces.append([self.M[:,row[0]],self.M[:,row[1]],self.M[:,row[2]], self.M[:,row[3]]])

        self.facesObj.set_verts(faces)
        self.bm.update()


    def canvas_coordinates_to_figure_coordinates(self,x_can,y_can):
        """
        Remap canvas coordinates to cube centered coordinates
        """

        (canvas_width,canvas_height)=self.canvas.get_width_height()
        figure_center_x = canvas_width/2+14
        figure_center_y = canvas_height/2+2
        x_fig = (x_can-figure_center_x)*self.pix2unit
        y_fig = (y_can-figure_center_y)*self.pix2unit

        return(x_fig,y_fig)


    def destroy(self):
        """
        Close function to properly destroy the window and tk with figure
        """
        try:
            self.destroy()
        finally:
            exit()


class BlitManager:
    def __init__(self, canvas, animated_artists=()):
        """
        Parameters
        ----------
        canvas : FigureCanvasAgg
            The canvas to work with, this only works for sub-classes of the Agg
            canvas which have the `~FigureCanvasAgg.copy_from_bbox` and
            `~FigureCanvasAgg.restore_region` methods.

        animated_artists : Iterable[Artist]
            List of the artists to manage
        """
        self.canvas = canvas
        self._bg = None
        self._artists = []

        for a in animated_artists:
            self.add_artist(a)
        # grab the background on every draw
        self.cid = canvas.mpl_connect("draw_event", self.on_draw)

    def on_draw(self, event):
        """Callback to register with 'draw_event'."""
        cv = self.canvas
        if event is not None:
            if event.canvas != cv:
                raise RuntimeError
        self._bg = cv.copy_from_bbox(cv.figure.bbox)
        self._draw_animated()

    def add_artist(self, art):
        """
        Add an artist to be managed.

        Parameters
        ----------
        art : Artist

            The artist to be added.  Will be set to 'animated' (just
            to be safe).  *art* must be in the figure associated with
            the canvas this class is managing.

        """
        if art.figure != self.canvas.figure:
            raise RuntimeError
        art.set_animated(True)
        self._artists.append(art)

    def _draw_animated(self):
        """Draw all of the animated artists."""
        fig = self.canvas.figure
        for a in self._artists:
            fig.draw_artist(a)

    def update(self):
        """Update the screen with animated artists."""
        cv = self.canvas
        fig = cv.figure
        # paranoia in case we missed the draw event,
        if self._bg is None:
            self.on_draw(None)
        else:
            # restore the background
            cv.restore_region(self._bg)
            # draw all of the animated artists
            self._draw_animated()
            # update the GUI state
            cv.blit(fig.bbox)
        # let the GUI event loop process anything it has to do
            cv.draw_idle()


if __name__ == "__main__":
    app = Arcball()
    app.mainloop()
    exit()
