#Class for building base stations
#Need to include:

#Number of antennas
#base station location
#antenna patterns
import numpy as np
import math as m
import matplotlib.pyplot as plt
import scipy.constants as scic
import random

class Base_Station:

    def __init__(self, room_size, num_antennas, frequency, spacing, edge = 0.3, angle = 45, antenna_pattern =0):
        #All the shared variables
        self.ant = num_antennas
        self.ant_pat = antenna_pattern
        self.room_size = room_size
        self.spacing = spacing
        self.angle = angle
        self.edge = edge
        self.frequency = frequency
        self.wave_len = (3e8)/frequency
        self.BS_sidewall = []
        self.BS_corner = []
        self.BS_corner_ref = []
        self.BS_sidewall_ref = []

        self.k = (2*scic.pi)/self.wave_len
        self.W =self.wave_len/2
        self.L = self.W
        self.blind_corner = []
        self.blind_sidewall =[]
        self.reflection_points = []
        self.reflectedCornerRef = []
        self.reflectedSideRef = []
    # def magnitude(self,vector):

    #     mag = np.sqrt(np.sum(vector*vector, axis = 1))
    #     mag = np.reshape(mag, (len(mag), -1))
    #     return mag

    # def norm(self,vector):
    #     return np.array(vector)/self.magnitude(np.array(vector))

    # def lineRayIntersectionPoint(self, rayOrigin, rayDirection, point1, point2):
        # rayOrigin = np.array(rayOrigin, dtype=np.float)
        # rayDirection = np.array(self.norm(rayDirection), dtype=np.float)
        # point1 = np.array(point1, dtype=np.float)
        # point2 = np.array(point2, dtype=np.float)
        
        # # Ray-Line Segment Intersection Test in 2D
        # v1 = rayOrigin - point1 

        # v2 = point2 - point1
        # v3 = np.array([-rayDirection[:,1], rayDirection[:,0]])
        # t1 = np.cross(v2, v1, axisc = 1) / np.dot(v2, v3)
        # v3 = np.rot90(v3)
        # t2 = np.sum(v1*v3, axis = 1)/np.sum(v2*v3, axis = 1)
        # #reflected = np.zeros(len(t1), dtype=bool)
        # reflected =((t1 >= 0) & (t2 >= 0.0) & (t2 <= 1.0))

        # return reflected
    #Reflects across a vertical line x= whatever
    def _horizon_reflect(self, py, reflection_pt): #from third to second quadrant
        new_py = []
        for j in range(len(py)):
            new_py.append((reflection_pt-py[j])*2 + py[j])
        return new_py
    #Reflects across a horizontal line y = whatever
    def _vertical_reflect(self, px, reflection_pt):# from third to fourth quadrant
        new_px = []
        for j in range(len(px)):
            new_px.append((reflection_pt-px[j])*2 + px[j])
        return new_px

    def patchAntenna(self, angle_deg,wave_length, phi = 0):
        #Normalized radiation model for a patch antenna 
        W = 0.5*wave_length
        L = W
        k = (2*m.pi)/wave_length
        phi = np.deg2rad(phi)
        angle = []
        for i in range(len(angle_deg)):
            angle.append(m.radians(angle_deg[i]))
        angle = np.array(angle)
        E_t = np.zeros(len(angle))
        for j in range(len(angle)):
            if phi == 0:
                if abs(angle[j]) > np.deg2rad(90):
                        E_t[j] = 0
                else:
                    E_t[j] = m.cos(((k*L)/2)*m.sin(angle[j]))
            else:
                # if abs(angle[j]) > np.deg2rad(90):
                #     E_t [j] = 0
                # else:
                    E_t[j] = (np.sin((k*W*np.sin(angle[j])*np.sin(phi))/2)/((k*W*np.sin(angle[j])*np.sin(phi))/2))*np.cos((k*L/4)*(np.sin(angle[j]+phi)+np.sin(angle[j]-phi)))*np.cos(angle[j])

        return E_t
        
    def calcAngle(self, BS_point, BS_ref, xxyy_points):
        lineA = np.array([[BS_point[0],BS_point[1]],[BS_ref[0],BS_ref[1]]] )
        lineB = np.array([[BS_point[0],BS_point[1]],[xxyy_points[0],xxyy_points[1]]])
        vA = [(lineA[0][0]-lineA[1][0]), (lineA[0][1]-lineA[1][1])]
        vB = [(lineB[0][0]-lineB[1][0]), (lineB[0][1]-lineB[1][1])]

        angle1 = m.atan2(vA[1], vA[0])
        angle2 = m.atan2(vB[1], vB[0])
        ang_deg = (angle1-angle2) * 360 /(2*m.pi)
        if ang_deg - 180 >= 0:
            return 360 - ang_deg
        return ang_deg
    def attenuation(self, loss_array, bs_array, ref_points,wave_len, xxyy, px, py):
        # Get angles of points 
        angles = np.zeros((len(xxyy)))
        lobe_effect = []
        for i in range(len(bs_array)):
            for j in range(len(xxyy)):
                angles[j] = self.calcAngle(bs_array[i,:], ref_points[i,:], xxyy[j,:])
            lobe_effect.append(self.patchAntenna(angles[:], wave_len, phi = 0))

        lobe_effect = np.reshape(lobe_effect, (len(bs_array), np.size(px, axis = 0)))
        antenna_applied = loss_array*lobe_effect
        #antenna_summed = sum_arrays(antenna_applied,antennas*4)
        return antenna_applied

    def corner_points(self,bs_1):
        end_array = []
        bs_1 = np.reshape(bs_1,(self.ant,2))
        bs_2y = np.reshape(self._horizon_reflect(bs_1[:,1],self.room_size[1]/2),(len(bs_1[:,0]),1))
        bs_2x = np.reshape(bs_1[:,0],(len(bs_1[:,0]),1))
        bs_2 = np.concatenate((bs_2x, bs_2y), axis = 1)

        bs_4x = np.reshape((self._vertical_reflect(bs_1[:,0], self.room_size[0]/2)),(len(bs_1),1))
        bs_4y = np.reshape(bs_1[:,1],(len(bs_1[:,1]),1))
        bs_4 = np.concatenate((bs_4x, bs_4y), axis = 1)

        bs_3 = np.concatenate((bs_4x,bs_2y), axis = 1)

        end_array = np.concatenate((bs_1, bs_2, bs_3, bs_4), axis = 0)
        return end_array


    def sidewall_points(self, offset_dist):
        end_array = []
        len_cen = self.room_size[0]/2
        wid_cen = self.room_size[1]/2
        ant1 = []
        ant2 = []
        ant3 = []
        ant4 = []
        if self.ant%2 == 1:  #Odd number of antennas
            offset_x = len_cen - (((self.ant -1)/2) * self.spacing)
            offset_y = wid_cen - (((self.ant -1) /2) * self.spacing)
        else: 
        #Even number of antennas
            offset_x = len_cen - ((((self.ant/2)-1) * self.spacing) + self.spacing*0.5)
            offset_y = wid_cen - ((((self.ant/2)-1) *self.spacing) + self.spacing*0.5)
        for i in range(self.ant):
            ant1.append([(offset_x+(i*self.spacing)),offset_dist])
        end_array.append(ant1)
        #y-axis 
        for f in range(self.ant):
            ant2.append([offset_dist,(offset_y+(f*self.spacing))])
        end_array.append(ant2) 
        for k in range(self.ant):
            ant3.append([self.room_size[0] -offset_dist,(offset_y+(k*self.spacing))])
        end_array.append(ant3)
        for h in range(self.ant):
            ant4.append([(offset_x+(h*self.spacing)),self.room_size[1]-offset_dist])
        end_array.append(ant4)
        end_array = np.array(np.reshape(end_array, (self.ant*4, 2)))

        return end_array
 
    def get_location(self):
        #Corner BS point creation
        x_incep = [0,0]
        y_incep = [0,0]
        y_side = 0
        x_side = 0
        slope = 0
        ant_array = self.edge*2 + (self.ant - 1)*self.spacing #0.3 added for tolerance on side
        if self.angle == 45:
            incep = ant_array/np.sqrt(2)
            y_incep = np.array([0, incep])
            x_incep = np.array([incep,0])
            y_side = incep
            x_side = incep
            slope = -1
            print(y_incep)
            print(x_incep)
        else:
            y_side = (ant_array*m.sin(m.radians(self.angle)))
            x_side = (ant_array*m.cos(m.radians(self.angle)))
            y_incep =np.array([0, y_side])
            x_incep = np.array([x_side, 0 ])
            slope = -(y_side/x_side)
        #Array with edges of line segments behind the arrays
        # self.blind_corner = [[[0,self.room_size[1]], [x_incep[0],self.room_size[1]], [(self.room_size[0]-x_incep[0]),self.room_size[1]], [self.room_size[0],self.room_size[1]]],
        #                      [[0,0], [x_incep[0],0], [(self.room_size[0]-x_incep[0]),0], [self.room_size[0],0]],   
        #                     [[self.room_size[0],0], [self.room_size[0], y_incep[1]], [self.room_size[0],self.room_size[1]-y_incep[1]], [self.room_size[0],self.room_size[0]]],
        #                     [[0,0],[0,y_incep[1]],[0,self.room_size[1]-y_incep[1]],[0,self.room_size[1]]]]
                    
 
        # self.blind_corner = np.array(self.blind_corner)

        f_antx = self.edge/np.sqrt((1+slope**2))
        f_ant =  [f_antx, (slope*f_antx+y_side)]
        bs_1 = []
        bs_1.append(f_ant)
        for i in range(self.ant -1):
            loc_x = bs_1[i][0] + self.spacing/np.sqrt(1+slope**2)
            bs_1.append([loc_x,(slope*loc_x+y_side)])
        self.BS_corner = self.corner_points(bs_1)
        self.BS_corner = np.array(self.BS_corner)
        

        ref_1 = np.ndarray(np.shape(bs_1))
        for i in range(np.size(bs_1, axis = 0 )):
            slope_bs = slope*-1
            inter_bs= bs_1[i][1] - slope_bs*bs_1[i][0]
            ref_1[i,0] = bs_1[i][0]+1
            ref_1[i,1] = slope_bs*ref_1[i,0] +inter_bs
        self.BS_corner_ref =self.corner_points(ref_1)
        self.BS_corner_ref = np.array(self.BS_corner_ref)
        
        
        
        #Sidewall BS creation
        wall_offset = 0.3
        ref_offset = 1
        self.BS_sidewall = self.sidewall_points(wall_offset)
        self.BS_sidewall = np.array(self.BS_sidewall)

        self.BS_sidewall_ref = self.sidewall_points(ref_offset)
        #Line segments for blind spots for side arrays
        self.blind_sidewall = np.array([[[self.BS_sidewall[0,0],self.room_size[1]],[self.BS_sidewall[self.ant-1,0],self.room_size[1]]],
                               [[self.BS_sidewall[0,0],0],[self.BS_sidewall[(self.ant-1),0],0]],
                               [[self.room_size[0],self.BS_sidewall[self.ant-1,1]],[self.room_size[0],self.BS_sidewall[self.ant*2-1,1]]],
                                [[0,self.BS_sidewall[self.ant-1,1]],[0,self.BS_sidewall[self.ant*2-1,1]]]])
             
    
    def los_path(self, node_px, node_py, xxyy, location = 'sidewall'):
        bs = ''
        if location == 'sidewall':
            bs = np.array(self.BS_sidewall)
        else:
            bs = np.array(self.BS_corner)

        dist = []
        for x in range(len(bs)):
            dist.append(np.linalg.norm(bs[x,:]- xxyy, axis = 1))
        dist = np.array(dist)
        tof = dist/(3e8)
        LOS_traces = (1/((2*scic.pi*dist)/self.wave_len))*np.exp(1.j*(2*scic.pi*dist/self.wave_len))
        #LOS_traces = 1*np.exp(1.j*(2*scic.pi*dist/self.wave_len))
        LOS_traces = np.reshape(LOS_traces, (len(bs), len(node_px), len(node_py)))
        return LOS_traces

    def mirroring(self, bs):
        walls = [[self._vertical_reflect(bs[:,0], 0),self._vertical_reflect(bs[:,0], self.room_size[0]),bs[:,0],bs[:,0]],
                  [ bs[:,1],bs[:,1],self._horizon_reflect(bs[:,1],0), self._horizon_reflect(bs[:,1], self.room_size[1])] ]
        walls = np.array(walls)
        walls = np.reshape(walls, (2,4*4*self.ant))
        walls = np.rot90(walls)
        self.reflection_points = walls
        return walls
    def wall_reflect(self, node_px, node_py, xxyy, location = 'sidewall'):
        node_px = np.array(node_px)
        node_py = np.array(node_py)
        
        bs = []
        if location == 'sidewall':
            bs = np.array(self.BS_sidewall)
        else:
            bs = np.array(self.BS_corner)
        walls = self.mirroring(bs)
        self.reflection_points = walls
        dist = []
        for x in range(np.size(walls, axis = 0)):
            dist.append(np.linalg.norm((walls[x,:]-xxyy), axis = 1))
        dist = np.array(dist)

        reflected_rays =  (1/((2*scic.pi*dist)/self.wave_len))*np.exp(1.j*scic.pi*dist/self.wave_len)
        reflected_rays = np.array(reflected_rays)
        
        reflected_rays  = np.reshape(reflected_rays, (len(walls), np.size(node_px), np.size(node_py)))
        return reflected_rays


    def scattering_node(self, node_px, node_py, xxyy, point = [0,0], location = 'sidewall', case = 'LOS', attenuation = 1):
        bs = ''
        xpt = 0
        ypt = 0
        if location == 'sidewall':
            bs = self.BS_sidewall
        else:
            bs = self.BS_corner
        # if point == al[0,0]):
        #     xpt = random.uniform(self.edge,np.int((self.room_size[0] - self.edge)))
        #     ypt = random.uniform(self.edge, (self.room_size[1]- self.edge))
        #     point = [xpt, ypt]
        point = np.reshape(point, (-1,2))

        min_distance = self.wave_len*10
        dist_scatter = np.linalg.norm(point - xxyy, axis = 1)
        dist_bs = []
        if case == 'LOS':
            dist_bs = np.linalg.norm(bs - point, axis = 1)
        else:
            dist_bs = self.wall_reflect(xpt, ypt, point, location)
        
        #path one with attenuation
        loss1= attenuation*(1/(((2*scic.pi*dist_scatter)/self.wave_len)))*np.exp(1.j*(2*scic.pi*dist_scatter/self.wave_len))
        print(np.shape(loss1))
        loss2 =(1/(((2*scic.pi*dist_bs)/self.wave_len)))*np.exp(1.j*(2*scic.pi*dist_bs/self.wave_len))
        
        totalLoss = []
        for i in range(len(loss2)):
            totalLoss.append(loss2[i]*loss1)
        totalLoss = np.array(totalLoss, dtype = 'complex64')
        losses = []
        for i in range(len(loss2)):
           losses.append(np.where(dist_scatter > min_distance, totalLoss, 0))
        #dists = np.array(dists)
        scattering = totalLoss
        
        #scattering = attenuation*np.exp(1.j*(2*scic.pi*dists/self.wave_len))
        if case == 'LOS':
            losses = np.reshape(scattering, (self.ant*4, len(node_px), len(node_py)))
        else:
            losses = np.reshape(scattering, (self.ant*4*4, len(node_px), len(node_py)))
        return losses
    
    def get_cornerBS(self):
        return self.BS_corner

    def get_sidewallBS(self):    
        return self.BS_sidewall

    def get_left_BS(self):
        return self.BS_sidewall[(self.ant -1):(self.ant*2 -1),:]

    def get_right_BS(self):
        return self.BS_sidewall[(self.ant*2 -1):(self.ant*3 -1),:]

    def get_bottom_BS(self):
        return self.BS_sidewall[0:self.ant,:]
        
    def get_top_BS(self):
        return self.BS_sidewall[(self.ant*3 -1): (self.ant*4 -1),:]

    def get_bottom_left(self):
        return self.BS_corner[0:(self.ant -1),:]
    def get_bottom_right(self):
        return self.BS_corner[(self.ant*3 - 1):(self.ant*4 -1),:]
    def get_top_left(self):
        return self.BS_corner[(self.ant-1):(self.ant*2 -1), :]
    def get_top_right(self):
        return self.BS_corner[(self.ant*2 -1):(self.ant*3 -1),:]


