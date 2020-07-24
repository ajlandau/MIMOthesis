#class for defining antenna pattern
import numpy as np
class Antenna_Pattern:
    def __init__(self, fileName):
        #define an array that will hold the degrees values and the gain values
        self.fileName = fileName
        self.pattern = np.zeros((360,2))

    def get_pattern(self):
        data = np.loadtxt(self.fileName, skiprows =2, usecols =(0,1,2))  
        x = 0 
        for i in range(len(data)):
            if data[i,1] >179:
                self.pattern[i,0] = 180+x
                self.pattern[i,1] = data[i,2]
                x = x+1            
            else:
                self.pattern[i,0] = data[i,0]
                self.pattern[i,1] = data[i,2]

        return self.pattern
    def ant_attenuation(self, BS_coord, UE_coord):
        #Slope of BS_coord normal vector to array
        #slope of linebetween BS and UE
        #determine if the UE slope of larger or smaller than BS_normal
        
        return dB_effect

    
