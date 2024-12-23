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
        self.obs_time = None
        self.log_time = None
        self.vessel_speed = None
        self.bot_depth = None
        self.pos_obs = Position()
        self.pos_vessel = Position()
        self.obs_depths = list()
        self.obs_sample = list()
        self.obs_ss = list()
        self.proc_depth = np.array([])        
        self.proc_ss = np.array([])        
        self.twtt_g = np.array([])    
        self.metadata = dict()
        self.metadata["angle_units"] = "rad"
        self.metadata["distance_units"] = "m"
        self.metadata["speed_units"] = "m/s"
        self.metadata["time_units"] = "s"
        self.metadata["count"] = None
        self.metadata["geoid"] = None
        self.metadata["ellipsoid"] = None
        self.metadata["chart_datum"] = None
        self.metadata["time_basis"] = "UTC"
        self.metadata["name"] = None    

    
    # B3.1
    def read_mvp_file(self, fullpath):
        # Check to see whether data already exists in the object
        if self.obs_depths:
            raise RuntimeError('SSP object already contains a profile')

        # Check the File's existence
        print(fullpath)

        if os.path.exists(fullpath):
            self.metadata["Source File"] = fullpath
            print('Opening sound speed profile data file:' + fullpath)
        else:  # Raise a meaningful error
            raise RuntimeError('Unable to locate the input file' + fullpath)

        # Open, read and close the file
        svp_file = open(fullpath)
        svp_content = svp_file.read()
        svp_file.close

        # Save the file name as meta data

        self.metadata["name"] = os.path.basename(fullpath)

        # Tokenize the contents
        svp_lines = svp_content.splitlines()

        n_header_lines = 0
        for line in svp_lines:
            # print(line) # Remove after step B3.1 has been successfully completed

            n_header_lines += 1

            if line == '':
                break

         
            #NOTE: The following is indented at the right level, using the indentation from B3.1 - 
            # adjusting it will create errors (assuming you indented the previous code correctly).
            
            # B3.4 Find the needed SSP metadata
            # Find and extract the time string
            if "gps time" in line.lower():
                # Extract the ZDA record
                obs_time = line.split()[2]
            
            # Find and extract the position string
            if "gps position" in line.lower():
                # Extract the GGA record
                obs_pos = line.split()[2]
               
            # Find and extract the depth string
            if "bottom depth" in line.lower():
                obs_depth = line.split()[2]
                        
            # Find and extract the vessel speed string
            if "ship speed:" in line.lower():
               obs_vessel_speed = line.split()[2]
            
        # B3.5 Check to see if There is a Date
        if 'obs_time' not in locals():
            raise RuntimeError('SSP.read_mvp_file(): Missing date and time info!')
        # B3.6.1 Parse the Time
        self.obs_time = ParseNMEA0183_ZDA(obs_time)
        
        # B3.6.2 Parse the Position
        self.pos_obs.ParseNMEA0183_GGA(
            obs_pos,
            "EGM08",
            "WGS84", 
            "geoid", 
            self.obs_time)
        
        # B3.6.3 Parse the Depth
        self.bot_depth=float(obs_depth)
        
        # B3.6.4 Parse the Speed
        self.vessel_speed=float(obs_vessel_speed)*1852/3600
        
        # B.3.7 Parsing the Record Identifier Line
        rec_type = svp_lines[n_header_lines].split(',')  
        
        # B.3.7.0 Find the Index of the Depth (m) Records
        index_depth = rec_type.index('Dpth(m)') 
        
        # B.3.7.1 Find the Index of the `'SV(m/s)'` Column
        index_ss = rec_type.index('SV(m/s)') 
        
        # B3.8 Parsing the Records
        for line in svp_lines[ n_header_lines + 1:]:
            obs = line.split(',')
            self.obs_depths.append( float(obs[index_depth]))
            self.obs_ss.append( float(obs[index_ss]))
        
        # B3.9 Sorting the Records
        temp = sorted(zip(self.obs_depths, self.obs_ss), key=lambda x: x[0])
        processed_depths, processed_ss = map(list, zip(*temp))
        # B3.10 Removing the Duplicates
        d_p=processed_depths[0]
        index = 0
        unwanted = []
        for d in processed_depths[1:]:
            index += 1
            if d  == d_p:
                unwanted.append(index)
            d_p = d

        for e in sorted( unwanted, reverse = True):
            del processed_depths[e]
            del processed_ss[e]
        
        # B3.11 Creating numpy Arrays
        self.proc_depth = np.array(processed_depths)
        self.proc_ss = np.array(processed_ss)

        # B3.12 Extend the Profiles to the Surface
        if self.proc_depth[0] > 0:
            self.proc_depth = np.insert(self.proc_depth,0,0)
            self.proc_ss = np.insert(self.proc_ss,0,self.proc_ss[0])

        # B3.13 Calculate the Sound Speed Gradients
        self.g = np.divide(np.diff(self.proc_ss),np.diff(self.proc_depth))
        
        # B3.14 Extending the Profiles to Full Ocean Depth    

        if self.proc_depth[-1] < 12000:
            self.proc_depth = np.append(self.proc_depth,12000)            
            self.proc_ss = np.append(self.proc_ss, self.proc_ss[-1] \
                        + 0.017 * (12000 - self.proc_depth[-2]))          
            self.g = np.append(self.g,0.017)   
          
        # B3.15 Replacing Zero Gradients
        self.g[self.g == 0] = 10**-9

        # B3.16 Updating the Profile
        for i in range(1,len(self.g)):
            self.proc_ss[i]=self.proc_ss[i-1]+(self.proc_depth[i] - \
            self.proc_depth[i-1])*self.g[i-1]

        return
    
    # B4.0 Create the `ray_trace_twtt` Method
    def ray_trace_twtt(self, d_start, th_start, ss_start, twtt):
 
        # B4.0.0 Parameter Initialization, Part I
        c = self.proc_ss      
        d = self.proc_depth
        g = self.g
                
        # B4.0.1 Parameter Initialization, Part II
        depth = 0
        rad_dist = 0
        layer_s = 0
        layer_e = 0
        
        # B4.0.2 Swapping the Depression Angle
        swap = False
        if th_start < 0 or th_start > pi:
            raise RuntimeError('SSP.ray_trace_twtt: depression angle th_start out of range')
        elif th_start > pi/2:
                swap = True
                th_start = pi - th_start
        
        # B4.0.3 Determine the Start Layer 
        layer_s = sum( d_start >= d) - 1
        
        # B4.0.4 Determine the Ray Constant
        ray_c = cos(th_start)/ss_start
        
        # B4.0.5 Calculate Ray Path Properties for Each Layer  
        dz = np.diff(d)
        r_curve = -1 / (g[0:] * ray_c)
        th = arccos(c[0:] * ray_c)
        dx = r_curve * (sin(th[1:]) - sin(th[:-1]))
        h_m = 2 / g[0:] * log(c[1:] / c[:-1])
        dt = h_m + (2 / g) * log((1 + sin(th[:-1])) / (1 + sin(th[1:])))      
       
        
        # B4.0.6 Calculate Inversion Sound Speed
        c_invert = 1/ray_c
#         print( "Inversion Sound Speed: %.2f"%(c_invert))  

               
        
        # B4.0.7 Determine Properties for the First Layer
        dx_init = r_curve[layer_s] * (sin(th_start) - sin(th[layer_s]))
        dz_init = d_start - d[layer_s]   
        dt_init = 2 / g[layer_s] * log(ss_start / c[layer_s])        
        dt_init += 2 / g[layer_s] * log((1 + sin(th[layer_s])) / (1 + sin(th_start)))      
        
       
       
        # B4.0.8 Accumulate From the Start Layer
        sum_dx = np.cumsum(dx[layer_s:])
        sum_dt = np.cumsum(dt[layer_s:])
        sum_dz = np.cumsum(dz[layer_s:])
        
        
        # B4.0.9 Offset Cumulative Sums by Values From the Start Layer
        sum_dx -= dx_init
        sum_dz -= dz_init
        sum_dt -= dt_init
        
        
        # B4.0.10  Determine the Number of Boundaries Crossed and the End Layer Index       
        n_bounds =  sum( twtt >= sum_dt[:-1])
        layer_e = n_bounds + layer_s 
            
 
        # B4.0.11  Determine Travel Time in Final Layer
        if n_bounds >= 0:
            t = twtt - sum_dt[n_bounds-1]
        else:
            raise RuntimeError('SSP start depth in same layer as reflector - not yet implemented')
      
        th_end = 2 * arctan(tanh(-t * g[layer_e] / 4 + arctanh(tan(th[layer_e] / 2))))            
        dx_end = r_curve[layer_e] * (sin(th_end) - sin(th[layer_e]))
        dz_end = -r_curve[layer_e] * (cos(th_end) - cos(th[layer_e]))

        # B4.0.12 Determining Depth and Radial Distance
        depth = sum_dz[n_bounds-1] + dz_end + d_start
        rad_dist = sum_dx[n_bounds-1] + dx_end 

        # B4.0.13 Determining Depth and Radial Distance
        if swap:
            rad_dist = -rad_dist
        # # B4.0.14 Return Results as Tuple        
        return (depth, rad_dist, layer_s, layer_e)
    
    # B4.1 Calculating the TWTT for a Given Depth
    def determine_twtt(self, d_start, th_start, ss_start, depth):
        
        # B4.0.0 Parameter Initialization, Part I
        c = self.proc_ss
        d = self.proc_depth
        g = self.g
       
        # B4.0.1 Parameter Initialization, Part II
        dist_z = depth - d_start
        rad_dist = 0
        layer_s = 0
        layer_e = 0
        
        # B4.0.2 Swapping the Depression Angle
        swap = False
        if th_start < 0 or th_start > pi:
            raise RuntimeError('SSP.ray_trace_twtt: depression angle th_start out of range')
        elif th_start > pi/2:
                swap = True
                th_start = pi - th_start
       
        # B4.0.3 Determine the Start Layer 
        layer_s = sum( d_start >= d) - 1
        # B4.0.4 Determine the Ray Constant
        ray_c = cos(th_start)/ss_start

        # B4.0.5 Calculate Ray Path Properties for Each Layer  
        dz = np.diff(d)
        r_curve = -1 / (g[0:] * ray_c)
        th = arccos(c[0:] * ray_c)
        dx = r_curve * (sin(th[1:]) - sin(th[:-1]))
        h_m = (2/g[0:]) * log(c[1:] / c[:-1])
        dt = h_m + (2 / g) * log((1 + sin(th[:-1])) / (1 + sin(th[1:])))

        # B4.0.6 Calculate Reversal Sound Speed
        c_invert = 1/ray_c
        # print( "Inversion Sound Speed: %.2f"%(c_invert))
        

        # B4.0.7 Determine Properties for the First Layer
        dx_init = r_curve[layer_s] * (sin(th_start) - sin(th[layer_s]))
        dz_init = d_start - d[layer_s]   
        dt_init = 2 / g[layer_s] * log(ss_start / c[layer_s])
        dt_init += 2 / g[layer_s] * log((1 + sin(th[layer_s])) / (1 + sin(th_start)))

        # B4.0.8 Accumulate From the Start Layer
        sum_dx = np.cumsum(dx[layer_s:])
        sum_dt = np.cumsum(dt[layer_s:])
        sum_dz = np.cumsum(dz[layer_s:])

        # B4.0.9 Offset Cumulative Sums by Values From the Start Layer
        sum_dx -= dx_init
        sum_dz -= dz_init
        sum_dt -= dt_init

        # B4.1.0 initialize twtt
        twtt = np.nan
        # B4.1.1 Determine the Number of Boundaries Crossed and the End Layer Index
        layer_e = sum(depth >= np.cumsum(dz))
        n_bounds = layer_e - layer_s
        
        # B4.1.2 Determine the Vertical Distance Traversed in the Last Layer
        dz_e = depth - d[layer_e]
       
        # B4.1.3 Determine the Sound Speed at the Final Depths
        c_e = c[layer_e] + g[layer_e] * dz_e
        # B4.1.4 Determine the Depression Angle at the Final Depth
        th_end = arccos(ray_c * c_e)
        
        # B4.1.5 Determine the final TWTT and dx      
        dt_e = (dz_e / sin(th[layer_e])) / c[layer_e] * 2
        dx_end =  r_curve[layer_e] * (sin(th_end) - sin(th[layer_e]))

        # B4.1.6 Determine the Total Two Way Travel Time
        if n_bounds > 0:
            twtt = sum_dt[n_bounds-1] + dt_e
            rad_dist = sum_dx[n_bounds-1] + dx_end 
        else:
            raise RuntimeError('SSP start depth in same layer as reflector - not yet implemented')
        
        # B4.1.7 Determining Depth and Radial Distance
        if swap:
            rad_dist = -rad_dist

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
        c = self.proc_ss
        d = self.proc_depth
        g = self.g
        layer = sum( d_interest >= d) - 1
        ss = c[layer]+(d_interest-d[layer])*g[layer]
        return ss
    # D6 Get the Sound Speed Data
    def read_jhc_file(self, fullpath):
        # Check to see whether data already exists in the object

        if self.obs_depths:
            raise RuntimeError('SSP object already contains a profile')

        # Check the File's existence
        if os.path.exists(fullpath):
            self.metadata["Source File"] = fullpath
            print('Opening sound speed profile data file:' + fullpath)
        else:  # Raise a meaningful error
            raise RuntimeError('Unable to locate the input file' + fullpath)

        # Open, read and close the file
        svp_file = open(fullpath)
        svp_content = svp_file.read()
        svp_file.close

        # Tokenize the contents
        svp_lines = svp_content.splitlines()
        self.obs_time = datetime.fromtimestamp(
            float(svp_lines[1].split()[0]), timezone.utc)
        self.log_time = datetime.fromtimestamp(
            float(svp_lines[2].split()[0]), timezone.utc)
        self.obs_latitude = float(svp_lines[3].split()[0])
        self.obs_longitude = float(svp_lines[3].split()[1])
        self.vessel_latitude = float(svp_lines[4].split()[0])
        self.vessel_longitude = float(svp_lines[4].split()[1])
        self.metadata["count"] = int(svp_lines[5].split()[0])

        count = 0  # initialize the counter for the number of rows read

        for svp_line in svp_lines[16:]:
            observations = svp_line.split()  # Tokenize the stringS
            self.obs_sample.append(float(observations[0]))
            self.obs_depths.append(float(observations[1]))
            self.obs_ss.append(float(observations[2]))
            count += 1

        if self.metadata["count"] != count:
            raise RuntimeError('Nr of Samples read ('+str(count) +
                            ') does not match metadata count (' +
                            str(self.metadata["count"])+')')

        # Process the data - in the jhc data files this is already a one-way profile,
        # this just for illustration
        array_ss = np.zeros((count, 3))

        # Sort the data samples by depth
        sorted_ss = sorted(zip(self.obs_depths, self.obs_ss))

        layer = 0
        for d, ss in sorted_ss:
            array_ss[[layer], [0]] = d
            array_ss[[layer], [1]] = ss
            layer += 1

        # Identify all the depths for which there are multiple observations
        mask = np.full((count, 1), True)
        mask[1:, [0]] = np.diff(array_ss[:, [0]], axis=0) != 0

        # Remove the duplicates - You really should get statistical representations here
        # but to keep this short just remove the duplicates
        array_ss = array_ss[mask[:, 0], ...]

        # Determine the gradients - Note the indexing: the gradient of the first layer 
        # is contained at the same index as the data for the TOP of the layer.
        array_ss[0:-1, [2]] = np.diff(array_ss[:, [1]],
                                        axis=0)/np.diff(array_ss[:, [0]], axis=0)

        # Estimate gradient for last layer assuming that the temperature and salinity remain the same
        # gradient solely a function of pressure (depth)
        array_ss[-1, [2]] = 0.017

        # Extend to 12000 m if necesarry - this is to get around some manufcturers requirements
        if self.obs_depths[-1] < 12000:
            ss = array_ss[-1:, [1]] + array_ss[-1:, [2]] \
            * (12000-array_ss[-1:, [0]])
            array_ss = np.vstack((array_ss, [12000, ss, 0.017]))

        # Make sure that the last gradient is 0.017
        array_ss[-1,2] = 0.017

        # Extend to 0 m if necesarry - assume well mixed
        if self.obs_depths[0] > 0:
            array_ss = np.vstack(
                ([0, array_ss[0, [1]], 0.], array_ss))

        # Step 5 Create a look-up array of twtts for each full layer
        # Allows for great gain in efficiency (do not have to calculate for each ping)
        self.twtt_layer = np.zeros((count, 1))

        for layer in range(0,self.metadata["count"]-1):
            if array_ss[layer, [2]] == 0:
                self.twtt_layer[layer] = 2 * \
                    (array_ss[layer+1, [0]] - array_ss[layer, [0]])/ \
                    array_ss[layer, [1]]
            else:
                self.twtt_layer[layer] = 2 / array_ss[layer, [2]] * \
                log(array_ss[layer+1, [1]]/array_ss[layer, [1]])

        self.proc_ss = array_ss[:,1]
        self.proc_depth = array_ss[:,0]
        self.g = np.diff(self.proc_ss)/np.diff(self.proc_depth)
        self.g[self.g == 0] = .24

        # Updating the Profile
        for i in range(1,len(self.g)):
            self.proc_ss[i]=self.proc_ss[i-1]+(self.proc_depth[i] - self.proc_depth[i-1])*self.g[i-1]

        