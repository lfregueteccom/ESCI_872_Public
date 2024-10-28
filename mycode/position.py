import numpy as np
from numpy import pi
import os
from datetime import datetime, timedelta, timezone

# A4.0.1 Import the geo-referencing modules
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy
import pyproj as proj 

# A2.0 The `Position` Class
class Position:
    """A Class for handling Position Data"""
    
    # A2.1 Class Initialization and Attributes
    def __init__(self):        
        self.times = list()
        self.latitudes = list()
        self.longitudes = list()
        self.heights = list()
        self.qualities = list()
        self.num_sats = list()
        self.hdops = list()        
        self.undulations = list()
        self.corr_ages = list()
        self.corr_stations = list()
        self.data_path = str()
        self.proj_pos = np.array([])             
        self.metadata = {
                "geodetic_units": "rad",
                "height_units": "m",
                "proj_units": "m",
                "geoid_name": None, 
                "ellipsoid_name": None, 
                "height_relative_to": None,
                "time_basis": "UTC",
                "proj_str": None}
    
    # A2.2 The String Representation Method
    def __str__(self): 
        # Loop through the keys and values
        txt = ''
        for key, value in self.metadata.items():
            txt += key + ': ' + str( value) + '\n'
        
        # If there are no latitudes, logically there are no longitudes either
        
        if len(self.latitudes):
            txt += "Minimum latitude       : %.6f\n" % (min(self.latitudes))           
            txt += "Maximum latitude       : %.6f\n" % (max(self.latitudes))       
            txt += "Minimum longitude      : %.6f\n" % (min(self.longitudes))           
            txt += "Maximum longitude      : %.6f\n" % (max(self.longitudes))       
            
        # Test heights separately, as we may just have horizontal positioning
            
        if len(self.heights):
            txt += "Minimum height         : %.2f%s\n" % \
                (min(self.heights), self.metadata["height_units"])         
            txt += "Maximum height         : %.2f%s\n" % \
                (max(self.heights), self.metadata["height_units"])
        else:
            txt += 'No height data present\n'
              
        
        if len(self.times):
            txt += "Start Time             : %s\n" % (min(self.times))
            txt += "End Time               : %s\n" % (max(self.times))
        else:
            txt += 'No time data present\n'
        return txt


    def read_hypack_raw_file(self, fullpath):
        
        # This function will currently only function provided that there are GGA sentences in the records.
        # You may update the function to include other positioning messages as well, but this is 
        # outside the scope of the class
                        
        # Check the File's existence
        if os.path.exists(fullpath):
            self.data_path = fullpath
            print('Opening GNSS data file:' + fullpath)
        else:  # Raise a meaningful error
            raise RuntimeError('Unable to locate the input file' + fullpath)
            
        # Open, read and close the file
        hypack_file = open(fullpath)
        hypack_content = hypack_file.read()
        hypack_file.close    

        # Split the file in lines
        hypack_records = hypack_content.splitlines()
        
        # Go through the header lines to find the date of the survey (not contained in the GGA records)
        
        lines_parsed=0
        for hypack_record in hypack_records:
            
            # Check for the time and date
            
            if hypack_record[:3].lower() == "tnd":
                hypack_datetime=datetime.strptime(hypack_record[4:23], "%H:%M:%S %m/%d/%Y")
                
                print("HYPACK RAW Header start time and date: " + hypack_datetime.ctime())
                
            # Keep track of the lines parsed
            lines_parsed+=1

            # Stop going through the records if the record starts with the string eoh (End Of Header)
            if hypack_record[:3].lower() == "eoh":
                break         
        
        # We are at the end of the header - start looking for the first GGA record and compare its time 
        # to the TND record
        # This is so that we can set the correct date
        
        # Keep track of the number of GGA records found
        
        num_gga_recs=0
        
        for hypack_record in hypack_records[lines_parsed:]:

            if hypack_record[19:22] == "GGA":
                gga_data=hypack_record.split()[3]
                self.ParseNMEA0183_GGA(gga_data,"EGM08","WGS84", "geoid", hypack_datetime)
            
                


    
    def write_hotlink(self, hotlink_path):
        
        fullpath, _ = os.path.splitext(self.data_path)
        
        fullpath = fullpath + "_pos.txt"
       
        # Check the File's existence

        if os.path.exists(fullpath):
            # Let the user now we are overwriting an existing file
            self.data_path = fullpath
            print('Overwriting file: ' + fullpath)
        else:  # Let the user know we are writing to a file
            print('Writing to file: ' + fullpath)
            
        output_file = open(fullpath, mode="w")  # mode="w" to open the file in writing mode
        
        # Write the header

        output_file.write('date time latitude longitude path\n')
        
        # Determine the duration associated to the positions in this object
        
        start_time = min(self.times)
        duration = max(self.times) - start_time
    
        # For each position write the date time longitude latitude path and fraction
        
        for i in range(0,len(self.times)):
            fraction = (self.times[i] - start_time) / duration
            line_content=str(self.times[i]) \
            + " %012.8f %013.8f %s?%f\n" % \
            (self.latitudes[i], self.longitudes[i], hotlink_path, fraction)
            output_file.write(line_content)


    def ParseNMEA0183_GGA(self, dt_str, geoid_name, ellipsoid_name, height_relative_to, date = None):
    
        # Get the GGA string and tokenize it
        gga_data = dt_str.split(',')

        # verify that we have a GGA string
        if not gga_data[0][-3:] == "GGA":
            raise RuntimeError(
                    'ParseNMEA0183_GGA: argument `dt_str` must be a GGA message')

        self.metadata["geoid_name"] = geoid_name
        self.metadata["ellipsoid_name"] = ellipsoid_name
        self.metadata["height_relative_to"] = height_relative_to
        
        if len(gga_data) < 2:
            return

        # Determine the time of day from both the header and the GGA string

        gga_timedelta=timedelta(hours=int(gga_data[1][0:2]), \
                                 minutes = int(gga_data[1][2:4]), \
                                 seconds = int(gga_data[1][4:6]))
        
        # Set the time of the date to midnight
        
        date = datetime(date.year, date.month, date.day, 0, 0, 0)
        self.times.append(date + gga_timedelta)

        # Parse the latitude
        if gga_data[3].lower() == "n":
            lat = float(gga_data[2][0:2])+float(gga_data[2][2:])/60.
        else:
            lat = -float(gga_data[2][0:2])-float(gga_data[2][2:])/60.   
        self.latitudes.append(lat)

        # Parse the longitude
        if gga_data[5].lower == "w":
            lon = float(gga_data[4][0:3])+float(gga_data[4][3:])/60.
        else:
            lon = -float(gga_data[4][0:3])-float(gga_data[4][3:])/60.
        self.longitudes.append(lon)

        # Parse the GNSS Quality indicator
        q = int(gga_data[6])
        self.qualities.append(q)

        # Parse the number of GNSS satellites used for the solution
        n_sats = int(gga_data[7])
        self.num_sats.append(n_sats)

        # Parse the HDOP Quality indicator
        hdop = float(gga_data[8])
        self.hdops.append(hdop)

        # Parse the orthometric height 
        height = float(gga_data[9])
        

        # Generate an error if the units of the orthometric height is not meters

        if not gga_data[10].lower() == "m":
            raise RuntimeError('Orthomeric height units are not meters!')
            
        self.heights.append(height)

        # Parse the geoid ellipsoid separation
        undulation = float(gga_data[11])
        self.undulations.append(undulation)
        
        if gga_data[12].lower() != "m":
            raise RuntimeError('Undulation height units are not meters!') 


        # If there is more data then parse it
        corr_age = None
        corr_station = None
        if not gga_data[13] == "":
            corr_ages.append(float(gga_data[13]))
            corr_stations.append(float(gga_data[14][0:-3]))

        # For now, ignore the checksum (this would become a computer science assignment

    # A2.3 Reading the data
    def read_jhc_file(self, fullpath):

        # 5.0 Set the reference ellipsoid and geoid
        self.metadata["ellipsoid_name"] = "WGS84"
        self.metadata["geoid_name"] = "EGM08"
        self.metadata["height_relative_to"] = "geoid"
        
        # Check the File's existence
        if os.path.exists(fullpath):
            self.data_path = fullpath
            print('Opening GNSS data file:' + fullpath)
        else:  # Raise a meaningful error
            raise RuntimeError('Unable to locate the input file' + fullpath)

        # Open, read and close the file
        gnss_file = open(fullpath)
        gnss_content = gnss_file.read()
        gnss_file.close
        
        times=list();

        # Tokenize the contents
        gnss_lines = gnss_content.splitlines()
        count = 0  # initialize the counter for the number of rows read
        for gnss_line in gnss_lines:
            observations = gnss_line.split()  # Tokenize the string
            time = datetime.fromtimestamp(
                float(observations[5]), timezone.utc)
            times.append(time)
            self.latitudes.append(float(observations[8])*pi/180)
            self.longitudes.append(float(observations[7])*pi/180)
            self.heights.append(float(observations[6]))
            count += 1

        self.times=np.asarray(times)
        

        
    # A4.0.2 Add the draw method to the Position class
    def draw(self, projection='auto'):
        print('Drawing Positioning Data')

        # A6.0 Determine the central latitude and longitude of the data in the Position object
        central_lat = (max(self.latitudes)+min(self.latitudes))/2*180/pi
        central_lon = (max(self.longitudes)+min(self.longitudes))/2*180/pi

        # A6.1 Create a figure and a title
        fig = plt.figure(figsize=(18, 6))
        fig.suptitle("Positioning  Data Plots")

        # A6.2 Create an Orthographic Coordinate Reference System (Orthographic CRS)
        crs_ortho = ccrs.Orthographic(central_lon, central_lat)

        # A6.3 Plot the Orthographic map
        ax1 = fig.add_subplot(1, 2, 1, projection=crs_ortho)
        ax1.set_global()

        # A6.4 Add Oceans, Land and a Graticule
        ax1.add_feature(cartopy.feature.OCEAN)
        ax1.add_feature(cartopy.feature.LAND, edgecolor='black')
        ax1.gridlines()
        ax1.set_title('Orthographic Map of Coverage Area')

        # A6.5 Plot the central point
        plt.plot(central_lon, central_lat, marker='o', markersize=7.0, markeredgewidth=2.5,
                 markerfacecolor='black', markeredgecolor='white',
                 transform=crs_ortho)

        # A6.6 Create a UTM Coordinate Reference System
        zone_number = int((np.floor((central_lon + 180) / 6) % 60) + 1) # calculate the UTM zone         
        if central_lat < 0:
            southern_hemisphere = True
        else:
            southern_hemisphere = False
        crs_utm = ccrs.UTM( zone=zone_number, southern_hemisphere=southern_hemisphere) # creating a CRS in UTUM c=with cartopy
        
        # A6.7 Plot the UTM map
        ax2 = fig.add_subplot( 1, 2, 2, projection=crs_utm)
        e_buffer=(np.max(self.longitudes)-np.min(self.longitudes))/10
        n_buffer=(np.max(self.latitudes)-np.min(self.latitudes))/5
        ax2.set_extent((np.min(self.longitudes)-e_buffer, np.max(self.longitudes)+e_buffer, np.min(
            self.latitudes)-n_buffer, np.max(self.latitudes)+n_buffer), crs=crs_utm)

        # A6.8 Add oceans, land,  graticule and Title
        ax2.add_feature(cartopy.feature.OCEAN)
        ax2.add_feature(cartopy.feature.LAND, edgecolor='black')
        ax2.gridlines()
        ax2.set_title('UTM Map of Positioning Data')    

        # A6.9 Plot the Positions
        for lon, lat in zip(self.longitudes, self.latitudes):
            plt.plot(lon, lat, marker='.', markersize=2,
                markerfacecolor='black', transform=crs_utm)

        # Display the figure
        plt.show()
        
    # A6.9.0 Add Projection Method
    def carto_project(self, projection_name, z_reference):

        # A6.9.1 Test the Projection Name
        # Test whether projection_name is a string
        if not isinstance(projection_name, str):
            raise RuntimeError(
                'Position.project(): argument `projection` must be of type str')

        # A6.9.2 Keep a List of Implemented Projections
        # Keep a list of projection that are implemented in this function
        implemented_projections = list()
        implemented_projections.append('utm')
        
        if projection_name.lower() not in implemented_projections:
            raise RuntimeError(
                'Position.project(): The projection `' + projection_name + '` is not yet implemented')
        
        # A6.9.4 Ensure a Reference Ellipsoid is Defined
        if self.metadata["ellipsoid_name"] == None:
            raise RuntimeError(
                'Position.carto_project(): Requires ellipsoid metadata to be defined!')
        
        # A6.9.5 Create the Projection String
        if projection_name.lower() == 'utm':
            proj_str = '+proj=utm'

            # A6.9.6  Determine Central Latitude and Longitude 
            central_lat = (max(self.latitudes)+min(self.latitudes))/2*180/pi
            central_lon = (max(self.longitudes)+min(self.longitudes))/2*180/pi
            
            # A6.9.7  Determine The UTM Zone Number
            zone_number = int((np.floor((central_lon + 180) / 6) % 60) + 1) # calculate the zone
            proj_str+= " +zone={}".format(zone_number)                
            
            # A6.9.8  Determine the Hemisphere
            if central_lat > 0:
                proj_str += ' +north'
            else:
                proj_str += ' +south'
            
            # A6.9.9 Set the Geodetic Datum for the input coordinates
            proj_str += ' +ellps=' + self.metadata["ellipsoid_name"]

            # A6.9.10 Set the Geodetic Datum for the output coordinates
            proj_str += ' +datum=' + self.metadata["ellipsoid_name"]


            # A6.9.11 Set the Units for the output coordinates
            proj_str += ' +units=' + self.metadata["proj_units"]

            # A6.9.12 Prevent Proj default behavior
            proj_str += ' +no_defs'

            # A6.9.13 Create a Pyproj Object
            proj_obj = proj.Proj(proj_str)

            # A6.9.14 Calculate the Easting and Northings
            E, N = proj_obj([l*180/pi for l in self.longitudes],[l*180/pi for l in self.latitudes])
            
            # Create a matrix of positions as 3D column vectors
            if z_reference.lower()=='ortho' or z_reference.lower()=='geoid':
                self.proj_pos=np.asarray([np.asarray(E),np.asarray(N),np.asarray(self.heights)])
            else:
                raise RuntimeError(\
                  'Position.carto_project(): currently only implemented for orthometric heights')
            
            # Add the string to the metadata
            self.metadata['proj_str'] = proj_str
        
        