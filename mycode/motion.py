import os
from datetime import datetime, timezone
import matplotlib.pyplot as plt
from numpy import pi, cos, sin, log, exp
import numpy as np

# 7.0: Creating Classes
class Motion:
    """A Class for handling Motion Data"""

    # 7.0.0 Class Initialization and Attributes
    def __init__(self):
        # The data attributes
        self.times = list()
        self.yaw = list()
        self.roll = list()
        self.pitch = list()
        self.heave = list()
        self.metadata = dict()
        self.metadata["angle__units"] = "rad"
        self.metadata["distance_units"] = "m"
        self.metadata["time_basis"] = "UTC"
        
    # 7.0.1 The String Representation Method
    def __str__(self):
        text = ""
        if self.times:
            text += "Start time: %s\n"%(str(self.times[0]))
            text += "End time:   %s\n"%(str(self.times[-1]))
        for key, value in self.metadata.items():
            text += "%s: %s\n"%(key, str(value)) 
                                       
        return text
    
    # 7.8.0 Creating a Read Method for the Motion Class
    def read_jhc_file(self, fullpath):

        # Check the File's existence
        if os.path.exists(fullpath):
            self.metadata["Source File"] = fullpath
            print('Opening motion data file:' + fullpath)
        else:  # Raise a meaningful error
            raise RuntimeError('Unable to locate the input file' + fullpath)

        # Open, read and close the file
        motion_file = open(fullpath)
        motion_content = motion_file.read()
        motion_file.close

        # Tokenize the contents
        motion_lines = motion_content.splitlines()
        count = 0  # initialize the counter for the number of rows read
        for motion_line in motion_lines:
            observations = motion_line.split()  # Tokenize the string
            time = datetime.fromtimestamp(
                float(observations[5]), timezone.utc)
            self.times.append(time)
            self.yaw.append(float(observations[6])*pi/180)
            self.roll.append(float(observations[7])*pi/180)
            self.pitch.append(float(observations[8])*pi/180)
            self.heave.append(float(observations[9]))
            count += 1
            
    # 7.9 Creating a Draw Method for the Motion Class
    def draw(self):
        print('Drawing Motion Data')
        
        # 7.9.0 Defining the plot area
        plt.figure(figsize=(20, 10))
        
        # 7.9.1 Creating subplots
        ax1 = plt.subplot(4, 1, 1)
        ax2 = plt.subplot(4, 1, 2, sharex=ax1)
        ax3 = plt.subplot(4, 1, 3, sharex=ax1)
        ax4 = plt.subplot(4, 1, 4, sharex=ax1, sharey=ax3)

        # 7.9.2 Plot the Data
        ax1.plot(self.times, np.degrees(self.yaw))
        ax2.plot(self.times, self.heave)
        ax3.plot(self.times, np.degrees(self.roll))
        ax4.plot(self.times, np.degrees(self.pitch))
        
        # 7.9.3 Labeling the Axes
        ax4.set_xlabel('Time ('+self.metadata['time_basis']+') →')
        plt.gcf().autofmt_xdate()

        ax1.set_ylabel('Heading [deg] →')
        ax2.set_ylabel('Heave [m] →')
        ax3.set_ylabel('Roll [deg] →')
        ax4.set_ylabel('Pitch [deg] →')

        #Last thing to do
        plt.show()
        
    # 7.10 Getting Motion Data as a Vector for a Specific Epoch
    def get_motion(self, time = datetime.fromtimestamp(0, timezone.utc)):

        # 7.10.1 Allocating Memory
        attitude = np.zeros(4)
        
        # 7.10.2 Map the Input Times to POSIX Times
        times = np.array([e.timestamp() for e in self.times])
        
        attitude[0] = np.interp(time.timestamp(), times, self.roll)
        attitude[1] = np.interp(time.timestamp(), times, self.pitch)
        attitude[2] = np.interp(time.timestamp(), times, self.yaw)
        attitude[3] = np.interp(time.timestamp(), times, self.heave)
        
        return attitude    
    
    # 7.11 Getting Motion Data as a Rotation Matrix for a Specific Epoch
    def get_rotation_matrix(self, time = datetime.fromtimestamp(0, timezone.utc)):
        # 7.11.1 Determining the Attitude
        att = self.get_motion(time)

        Rx = np.array([[1, 0,            0          ],
                       [0, cos(att[0]), -sin(att[0])],
                       [0, sin(att[0]),  cos(att[0])]])

        Ry = np.array([[ cos(att[1]),  0, sin(att[1])],
                       [ 0          ,  1, 0          ],
                       [-sin(att[1]),  0, cos(att[1])]])

        Rz = np.array([[cos(att[2]), -sin(att[2]), 0],
                       [sin(att[2]),  cos(att[2]), 0],
                       [0          ,  0          , 1]])
        

        return Rz@Ry@Rx
    
    def geo_reference_la(self, time = datetime.fromtimestamp(0, timezone.utc), la = np.array([[0],[0],[0]])):
        R = self.get_rotation_matrix(time)
        la_n = R@la
        return la_n
        
#     def pos_to_rp(self, time = datetime.fromtimestamp(0, timezone.utc), pos = np.array([]), la = np.array([])):
#         la_n = self.geo_reference_la(time, la) #Already done in student's case (2021, not after)
#         # Be careful to copy the position as it is passed in by reference
#         rp = pos.copy()
        
#         # Be careful about the alignment of the axes
#         rp[0] -= la_n[1]
#         rp[1] -= la_n[0]
#         rp[2] += la_n[2]
        
#         return rp
    
