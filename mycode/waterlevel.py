import os
from datetime import datetime, timezone
import matplotlib.pyplot as plt

class WaterLevel:
    """A Class for handling Water Level Data"""

    def __init__(self):
        pass
        # The data attributes
        self.times = list()
        self.water_levels = list()
        self.metadata = dict()
        self.data_path=str()
        self.metadata["units"] = "m"
        self.metadata["datum_type"] = None
        self.metadata["datum_name"] = None
        self.metadata["time_basis"] = "UTC"        
        self.metadata["location_name"] = "Unknown"


    def __str__(self): 
        # A1.2 The String Representation Method
        txt = "Vessel Metadata: "
        for key, value in self.metadata.items():
            txt += key + ": " + str( value) + '\n'
        
        txt  = "Location name\t\t: %s\n" % (self.metadata["location_name"])
        txt  += "Reference Surface Type\t: %s\n" % (self.metadata["datum_type"])
        txt  += "Reference Surface Name\t: %s\n" % (self.metadata["datum_name"])
        txt  += "Observation Time Basis\t: %s\n" % (self.metadata["time_basis"])
        txt  += "Observations Units\t: %s\n" % (self.metadata["units"])    

        # A1.3 More on the String Representation Method
        # The code fragment below assumes that there is some string called txt
        if len(self.water_levels):           
            txt += "Start Time\t: %.2f%s\n" % (min(self.times), self.metadata["time_basis"])           
            txt += "End Time\t: %.2f%s\n" % (max(self.times), self.metadata["time_basis"])
            txt += "Minimum Water Level\t: %.2f%s\n" % (min(self.water_levels), self.metadata["units"])           
            txt += "Maximum Water Level\t: %.2f%s\n" % (max(self.water_levels), self.metadata["units"])
        else:
            txt += "No time data present\n"
            txt += "No water level data present"
        return txt
    
    def read_jhc_file(self, fullpath):

        # Check the File's existence
        if os.path.exists(fullpath):
            self.metadata["Source File"] = fullpath
            print('Opening water level data file:' + fullpath)
        else:  # Raise a meaningful error
            raise RuntimeError('Unable to locate the input file' + fullpath)

        # Open, read and close the file
        wl_file = open(fullpath)
        wl_content = wl_file.read()
        wl_file.close

        # Tokenize the contents
        wl_lines = wl_content.splitlines()
        count = 0  # initialize the counter for the number of rows read
        for wl_line in wl_lines:
            observations = wl_line.split()  # Tokenize the string
            epoch=datetime.fromtimestamp(float(observations[5]), timezone.utc)
            self.times.append(epoch)
            self.water_levels.append(float(observations[6]))
            count += 1
            
    def draw(self):
        plt.figure(figsize=(10, 10))
        print('Drawing Water Level Data')
  
        # plotting the points  
        plt.plot(self.times, self.water_levels) 
        plt.title('Water Levels in [m]') 
        plt.ylabel('Water Level in [m] →') 
        plt.xlabel('Time ('+self.metadata['time_basis']+') →') 
        plt.gcf().autofmt_xdate()