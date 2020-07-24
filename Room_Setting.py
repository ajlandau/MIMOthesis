#Class for creating self
#Recieve:
#Size of self
#Frequency
#antenna spacing
import numpy as np
class Room_Setting:
    def __init__(self, length, width, freq, lambda_step, edge_gap):
        self.length = length
        self.width = width
        self.wave_len = (3e8)/freq
        self.step = lambda_step
        self.edge = edge_gap
        self.px= []
        self.py = []
        return self.points()
    def points(self):
        #fix so that x and y can be different length vectors    
        src_x = (self.length - self.edge*2) / (self.wave_len/self.step)
        src_y = (self.width - self.edge*2) / (self.wave_len/self.step)
        xpts = np.array([(self.edge+((self.wave_len/self.step)*i)) for i in range(int(src_x+1))])
        ypts = np.array([(self.edge+((self.wave_len/self.step)*i)) for i in range(int(src_y+1))])
        
        self.px = xpts #[x,y]
        self.py = ypts
        
