# B.1 Reading the Sound Speed Profiles
import os
from datetime import datetime, timezone
import numpy as np
from numpy import pi, cos, sin, log, exp, arccos, tan, arctan, tanh, arctanh
from mycode.position import *

# B.1.0 Creating the SSP Class
class SSP:
    """A Class for handling Sound Speed Profile data"""

    def __init__(self):
        # B3.0.0 Add Attribute Variables to the SSP class
    
    # B3.1
    def read_mvp_file(self, fullpath):
         
            #NOTE: The following is indented at the right level, using the indentation from B3.1 - 
            # adjusting it will create errors (assuming you indented the previous code correctly).
            
            # B3.4 Find the needed SSP metadata
            # Find and extract the time string
            
            # Find and extract the position string
               
            # Find and extract the depth string
                        
            # Find and extract the vessel speed string
            
        # B3.5 Check to see if There is a Date
        
        # B3.6.1 Parse the Time
        
        # B3.6.2 Parse the Position
        
        # B3.6.3 Parse the Depth
        
        # B3.6.4 Parse the Speed
        
        # B.3.7 Parsing the Record Identifier Line
        
        # B.3.7.0 Find the Index of the Depth (m) Records
        
        # B.3.7.1 Find the Index of the `'SV(m/s)'` Column
        
        # B3.8 Parsing the Records
        
        # B3.9 Sorting the Records
        
        # B3.10 Removing the Duplicates
        
        # B3.11 Creating numpy Arrays
        
        # B3.12 Extend the Profiles to the Surface
        
        # B3.13 Calculate the Sound Speed Gradients
        
        # B3.14 Extending the Profiles to Full Ocean Depth
          
        # B3.15 Replacing Zero Gradients
        
        # B3.16 Updating the Profile
        return
    
    # B4.0 Create the `ray_trace_twtt` Method
    def ray_trace_twtt(self, d_start, th_start, ss_start, twtt):
        
        # B4.0.0 Parameter Initialization, Part I
        
        # B4.0.1 Parameter Initialization, Part II
        
        # B4.0.2 Swapping the Depression Angle
        
        # B4.0.3 Determine the Start Layer 
        
        # B4.0.4 Determine the Ray Constant
        
        # B4.0.5 Calculate Ray Path Properties for Each Layer  
        
        # B4.0.6 Calculate Inversion Sound Speed
        
        # B4.0.7 Determine Properties for the First Layer
        
        # B4.0.8 Accumulate From the Start Layer
        
        # B4.0.9 Offset Cumulative Sums by Values From the Start Layer
        
        # B4.0.10  Determine the Number of Boundaries Crossed and the End Layer Index       
        
        # B4.0.11  Determine Travel Time in Final Layer
        
        # B4.0.12 Determining Depth and Radial Distance
        
        # B4.0.13 Determining Depth and Radial Distance
        
        # B4.0.14 Return Results as Tuple
        return (depth, rad_dist, layer_s, layer_e)
    
    # B4.1 Calculating the TWTT for a Given Depth
    def determine_twtt(self, d_start, th_start, ss_start, depth):
        
        # B4.0.0 Parameter Initialization, Part I
       
        # B4.0.1 Parameter Initialization, Part II
        
        # B4.0.2 Swapping the Depression Angle
       
        # B4.0.3 Determine the Start Layer 
        
        # B4.0.4 Determine the Ray Constant
        
        # B4.0.5 Calculate Ray Path Properties for Each Layer  
        
        # B4.0.6 Calculate Reversal Sound Speed
        
        # B4.0.7 Determine Properties for the First Layer
        
        # B4.0.8 Accumulate From the Start Layer
        
        # B4.0.9 Offset Cumulative Sums by Values From the Start Layer
        
        # B4.1.0 initialize twtt
        
        # B4.1.1 Determine the Number of Boundaries Crossed and the End Layer Index
        
        # B4.1.2 Determine the Vertical Distance Traversed in the Last Layer
        
        # B4.1.3 Determine the Sound Speed at the Final Depth
        
        # B4.1.4 Determine the Depression Angle at the Final Depth
        
        # B4.1.5 Determine the final TWTT and dx      
        
        # B4.1.6 Determine the Total Two Way Travel Time
        
        # B4.1.7 Determining Depth and Radial Distance
        
        # B4.1.8 Return Results as Tuple
        return ( twtt, rad_dist, layer_s, layer_e)               



    # B5.0  Plotting the Profiles
    def draw(self, full_profile=False, ax1=False, depth_range=False, ss_range=False, label=True):

        # Get a view to the processed depths and sound speeds
        d = self.proc_depth
        c = self.proc_ss
        
        if ax1 == False:
            fig, ax1 = plt.subplots()   
            
        if full_profile:
            if depth_range == False:
                depth_range = (min(d), max(d))
            if ss_range == False:
                ss_range = (min(c), max(c))
            plt.plot(c[0:], d[0:])
        else:
            if depth_range == False:
                depth_range = (min(d[0:-1]), max(d[0:-1]))
            if ss_range == False:
                ss_range = (min(c[0:-1]), max(c[0:-1]))
            plt.plot(c[1:-1], d[1:-1])
            
        plt.ylim(depth_range)
        plt.xlim(ss_range)
        
        if label:
            plt.ylabel('← Depth [m]')
        else:
            labels = [item.get_text() for item in ax1.get_yticklabels()]
            empty_string_labels = ['']*len(labels)
            ax1.set_yticks(ax1.get_yticks().tolist()) # Bug in Matplotlib requires this
            ax1.set_yticklabels(empty_string_labels)
            
        plt.xlabel('Sound Speed [m/s] →')
        ax1.invert_yaxis()
        ax1.xaxis.tick_top()
        ax1.xaxis.set_label_position('bottom')
        
        # Set the title from the file name that contained the data
        ax1.title.set_text(os.path.splitext(self.metadata['name'])[0])


    # B5.1.1 Add Method to SSP to Determine c at Given Depth
    def determine_c( self, d_interest):
        
        return ss
    