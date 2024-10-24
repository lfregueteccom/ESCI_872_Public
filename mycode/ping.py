import os
import numpy as np
from numpy import pi
from datetime import datetime, timezone, timedelta

class Ping:
    """A Class for handling Ping Specific Data"""

    def __init__(self):
        
        # Ping specific
        self.lat        = None
        self.lon        = None
        self.tx_roll    = None
        self.tx_pitch   = None
        self.tx_heading = None
        self.tx_heave   = None
        self.tx_time    = None
        
        # Sector specific
        self.tilt_angle        = np.array([])
        self.focus_range       = np.array([])
        self.signal_length     = np.array([])
        self.tx_time_offset    = np.array([])
        self.center_frequency  = np.array([])
        self.mean_absorp_coeff = np.array([])
        self.waveform_id       = np.array([])
        self.sector_id         = np.array([])
        self.bandwidth         = np.array([])
        
        # Beam specific
        
        self.beam             = np.array([], dtype = np.int64)
        self.across           = np.array([])
        self.along            = np.array([])
        self.depth            = np.array([])
        self.samp_nr          = np.array([])
        self.twtt             = np.array([])
        self.range            = np.array([])
        self.steer_rx         = np.array([])
        self.sector           = np.array([])
        self.steer_tx         = np.array([])
        self.tx_t_offset_beam = np.array([],dtype = timedelta)
        self.ctr_freq_beam    = np.array([])
        
        self.metadata = dict()

    def read(self, fullpath):
        
        # Check to see whether data already exists in the object
        if self.lat:
            raise RuntimeError('Ping object already contains data')

        # Check the File's existence
        if os.path.exists(fullpath):
            self.metadata["Source File"] = fullpath
            print('Opening sound speed profile data file:' + fullpath)
        else:  # Raise a meaningful error
            raise RuntimeError('Unable to locate the input file' + fullpath)

        # Open, read and close the file
        ping_file = open(fullpath)
        ping_content = ping_file.read()
        ping_file.close
        
        # Save the file name as meta data
        
        self.metadata["name"] = os.path.basename(fullpath)

        # Tokenize the contents
        ping_lines = ping_content.splitlines()
        
        # Read the general ping parameters
        
        self.metadata["Ping_ID"]=ping_lines[0].split()[1]
        
        obs = ping_lines[1].split()
        self.lat        = float(obs[1])*pi/180
        self.lon        = float(obs[2])*pi/180
        
        obs = ping_lines[2].split()
        self.tx_roll    = float(obs[2])*pi/180
        self.tx_pitch   = float(obs[4])*pi/180
        self.tx_heading = float(obs[8])*pi/180
        self.tx_heave   = float(obs[6])
        self.tx_time    = datetime.fromtimestamp(float(ping_lines[3].split()[2]), timezone.utc)
        
        # Determine the number of sectors and allocate memory
        
        curr_line = 5
        n_sectors = 0
        while ping_lines[curr_line].split()[0] == 'Transmit':
            n_sectors +=1
            curr_line += 10
            
        # Allocate memory for the arrays
        self.tilt_angle        = np.zeros(n_sectors)
        self.focus_range       = np.zeros(n_sectors)
        self.signal_length     = np.zeros(n_sectors)
        self.tx_time_offset    = np.zeros(n_sectors)
        self.center_frequency  = np.zeros(n_sectors)
        self.mean_absorp_coeff = np.zeros(n_sectors)
        self.waveform_id       = np.zeros(n_sectors, dtype = int)
        self.sector_id         = np.zeros(n_sectors, dtype = int)
        self.bandwidth         = np.zeros(n_sectors)
        
        # Read the sector specific info
       
        curr_line = 5    
        for i in range(n_sectors):
            self.tilt_angle[i]        = float(ping_lines[curr_line+1].split()[-1])*pi/180
            self.focus_range[i]       = float(ping_lines[curr_line+2].split()[-1])
            self.signal_length[i]     = float(ping_lines[curr_line+3].split()[-1])
            self.tx_time_offset[i]    = float(ping_lines[curr_line+4].split()[-1])
            self.center_frequency[i]  = float(ping_lines[curr_line+5].split()[-1])
            self.mean_absorp_coeff[i] = float(ping_lines[curr_line+6].split()[-1])/100
            self.waveform_id[i]       = int(ping_lines[curr_line+7].split()[-1])
            self.sector_id[i]         = int(ping_lines[curr_line+8].split()[-1])
            self.bandwidth[i]         = float(ping_lines[curr_line+9].split()[-1])
            curr_line += 10
                  
        # Determine the number of beams and allocate memory
        curr_line += 1
        num_beams = len(ping_lines)-curr_line
        
        self.beam             = np.zeros(num_beams, dtype = np.int64)
        self.across           = np.zeros(num_beams)
        self.along            = np.zeros(num_beams)
        self.depth            = np.zeros(num_beams)
        self.samp_nr          = np.zeros(num_beams)
        self.twtt             = np.zeros(num_beams)
        self.range            = np.zeros(num_beams)
        self.steer_rx         = np.zeros(num_beams)
        self.sector           = np.zeros(num_beams)
        self.steer_tx         = np.zeros(num_beams)
        self.tx_t_offset_beam = np.zeros(num_beams, dtype = timedelta)
        self.ctr_freq_beam    = np.zeros(num_beams)

        # Read the Beam specific info
       
        for i in range(num_beams):
            obs =  ping_lines[curr_line].split()

            self.beam[i]             = np.int64(obs[ 0])
            self.across[i]           = float(obs[ 1])
            self.along[i]            = float(obs[ 2])
            self.depth[i]            = float(obs[ 3])
            self.samp_nr[i]          = float(obs[ 4])
            self.twtt[i]             = float(obs[ 5])
            self.range[i]            = float(obs[ 6])
            self.steer_rx[i]         = float(obs[ 7])*pi/180
            self.sector[i]           = float(obs[ 8])
            self.steer_tx[i]         = float(obs[ 9])*pi/180
            self.tx_t_offset_beam[i] = timedelta(seconds = float(obs[10]))
            self.ctr_freq_beam[i]    = float(obs[11])
            curr_line += 1
        
        
    def __str__(self):
        
        np.set_printoptions(formatter={'all':lambda x: ' '+str('%.4f'%x)})

        s  = 'Ping     : ' + self.metadata["Ping_ID"] + '\n'
        s += 'Latitude : ' + str(self.lat*180/pi) + u"\N{DEGREE SIGN}" + '\n'
        s += 'Longitude: ' + str(self.lon*180/pi) + u"\N{DEGREE SIGN}" + '\n'
        s += 'Roll     : ' + str(self.tx_roll*180/pi) + u"\N{DEGREE SIGN}" + '\n'
        s += 'Pitch    : ' + str(self.tx_pitch*180/pi) + u"\N{DEGREE SIGN}" + '\n'
        s += 'Headding : ' + str(self.tx_heading*180/pi) + u"\N{DEGREE SIGN}" + '\n'
        s += 'Heave    : ' + str(self.tx_heave) + ' [m]\n'
        s += 'Tx Time  : ' + self.tx_time.strftime("%H:%M:%S  %B %-d, %Y") + '\n'
        
        s += '\n'
        s += str(len(self.tilt_angle)) + ' Sectors --- \n'
       

        for i in range(len(self.tilt_angle)):
            s += '--- Sector ' + str(self.sector_id[i]) + '\n'
            s += 'Tx Tilt Angle   : ' + str(self.tilt_angle[i]*180/pi) + u"\N{DEGREE SIGN}" + '\n'
            s += 'Focus Range     : ' + str(self.focus_range[i]) + ' [m]\n'
            s += 'Signal Length   : ' + str(self.signal_length[i]) + ' [s]\n'
            s += 'Tx time offset  : ' + str(self.tx_time_offset[i]) + ' [s]\n'
            s += 'Tx central freq : ' + str(self.center_frequency[i]) + ' [Hz]\n'
            s += 'Mean att coeff  : ' + str(self.mean_absorp_coeff[i]) + ' [dB/km]\n'
            s += 'Waveform ID     : ' + str(self.waveform_id[i]) + '\n'
            s += 'Tx bandwidth    : ' + str(self.bandwidth[i]) + ' [Hz]\n'

            
            
        s += '\n---\nBeam 0 data:\n'    
            
        # The data for the first beam
        s += 'Beam nr        : ' + str(self.beam[0]) + '\n'
        s += 'Across dist    : ' + str(self.across[0]) + ' [m]\n'
        s += 'Along dist     : ' + str(self.along[0]) + ' [m]\n'
        s += 'Depth          : ' + str(self.depth[0]) + ' [m]\n'
        s += 'Sample nr      : ' + str(self.samp_nr[0]) + '\n'
        s += 'TWTT           : ' + str(self.twtt[0]) + ' [s]\n'
        s += 'Slant range    : ' + str(self.range[0]) + ' [m]\n'
        s += 'Rx steer angle : ' + str(self.steer_rx[0]) + u"\N{DEGREE SIGN}" + '\n'
        s += 'Sector id      : ' + str(self.sector[0]) + '\n'
        s += 'Tx steer angle : ' + str(self.steer_tx[0]) + u"\N{DEGREE SIGN}" + '\n'
        s += 'Time delay     : ' + str(self.tx_t_offset_beam[0]) + ' [m]\n'
        s += 'Central freq   : ' + str(self.ctr_freq_beam[0]) + ' [Hz]\n'
        return s
    
    def get_beam_index(self, beam_nr):
        
        beam_nr = int(beam_nr)
        
        if beam_nr < min(self.beam) or beam_nr > max(self.beam):
            return None
        else:
            return np.where(self.beam == beam_nr)[0][0]

        
        
        
        
