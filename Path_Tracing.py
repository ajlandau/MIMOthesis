#Class for path-loss
#Need to include
#Line of Sight
#Reflections off walls
#Perhaps this just needs to be part of room design
#Must have:     
#Room shape  and size, base station position and antenna objects.
import numpy as np
import scipy as sci
import scipy.constants as scic
import random

class Path_Tracing:
    def __init__(self, room_size, node_px, node_py, base_station, frequency, antenna_pat = [0,0], wall_attenuation = 1 ):
        self.room_size = room_size
        self.base_station = base_station
        self.ant_pat = antenna_pat
        self.wall_attenuation = wall_attenuation
        self.node_px = node_px
        self.node_py = node_py
        self.wave_len = 3e8/frequency
    
    def LOS_path(self):
        LOS_traces = np.zeros((np.size(self.base_station,0), np.size(self.node_px), np.size(self.node_py), 2), dtype = complex)
    #need to fix points return since px and py can be different length vectors
        for i in range(len(self.base_station)):
            for j in range(len(self.node_px)):
                for k in range(len(self.node_py)):
                    dist = np.linalg.norm(self.base_station[i] - (self.node_px[j],self.node_py[k]))
                    LOS_traces[i,j,k,1] = dist/(3e8) #tof of signal
                    LOS_traces[i,j,k,0] = (1/(((4*sci.pi*dist)/self.wave_len)**2))*np.exp(1.j*(2*scic.pi*dist/self.wave_len))

        return LOS_traces

    def Reflect_LW(self):
        ref_LW = np.zeros((np.size(self.base_station,0), np.size(self.node_px), np.size(self.node_py), 2), dtype = complex)
        for i in range(len(self.base_station)):
            BS_point =[ -1* self.base_station[i, 0], self.base_station[i,1]]
            for j in range(len(self.node_px)):
                for k in range(len(self.node_py)):
                    dist = np.linalg.norm(BS_point - (self.node_px[j],self.node_py[k]))
                    slope = (BS_point[1] - self.node_py[k])/(BS_point[0] - self.node_px[j])
                    angle = np.tan(slope)
                    ref_LW[i,j,k,1] = dist/(3e8)
                    ref_LW[i,j,k,0] = (1/(((4*scic.pi*dist)/self.wave_len)**2))*np.exp(1.j*(2*scic.pi*dist)/self.wave_len)
        return ref_LW



    def Reflect_BW(self):
        ref_BW = np.zeros((np.size(self.base_station,0), np.size(self.node_px), np.size(self.node_py), 2), dtype = complex)
        for i in range(len(self.base_station)):
            BS_point =[self.base_station[i, 0], -1*self.base_station[i,1]]
            for j in range(len(self.node_px)):
                for k in range(len(self.node_py)):
                    dist = np.linalg.norm(BS_point - (self.node_px[j],self. node_py[k]))

        return ref_BW

    def Reflect_RW(self):
        ref_RW = np.zeros((np.size(self.base_station,0), np.size(self.node_px), np.size(self.node_py), 2), dtype = complex)
        for i in range(len(self.base_station)):
            BS_point =[(2*(self.room_size[0] - self.base_station[i, 0]) + self.base_station[i,0] ), self.base_station[i,1]]
            for j in range(len(self.node_px)):
                for k in range(len(self.node_py)):
                    dist = np.linalg.norm(BS_point - (self.node_px[j],self.node_py[k]))

        return ref_RW
"""
    def Reflect_TW(self):
        ref_TW = np.zeros((np.size(self.base_station,0), np.size(self.node_px), np.size(self.node_py), 2), dtype = complex)
        for i in range(len(self.base_station)):
            BS_point =[self.base_station[i, 0], (self.base_station[i,1] * self.room_size[1]]
            for j in range(len(self.node_px)):
                for k in range(len(self.node_py)):
                    dist = np.linalg.norm(BS_point - (self.node_px[j],self. node_py[k]))

        return ref_TW

"""

    #def Reflect_Obstacle(self, rand, size):
"""
    def scattering_point(self, location = [0,0], attenuation = 1):
        if location == [0,0]:
            ptx = random.randint(edge)
            """