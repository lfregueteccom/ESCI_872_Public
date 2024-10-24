import numpy as np
from scipy.interpolate import interp1d
from numpy import pi, cos, sin, log, exp
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

from mycode.twtt import TWTT
from mycode.motion import Motion
from mycode.waterlevel import WaterLevel
from mycode.vessel import Vessel
from mycode.position import Position

# A8.2 Integration
class Integration:
    """A Class for Integrating Data to Create Soundings"""

    def __init__(self, twtt, pos, motions, sound_speed_profile, water_levels, vessel):
        
        # A8.2.0 Data for integration
        # For now we can only integrate if the Positions have been projected - we will use 
        # either UTM or UPS depending on the latitude
        
        # A8.2.1 Determine the number of pings

        # A8.2.2 Memory Allocation
              
        # A8.2.3 Converting the `datetime` objects to POSIX Times

        # A8.2.4.0 Transmit event interpolation
        
        # Determine the interpolation function for the positions

        # A8.2.4.1  Receive time determination

        # A8.2.4.2 Receive Event interpolation

        # A8.2.5  Allocating Memory for Processed Data
        
        # A8.2.6  Looping through the TWTTs
        ping = 0
        for t in t_twtt:
            # A8.2.7 Calculate the Eucliden Euler Angle Rotation Matrices
            

            # A8.2.8 Compound Rotation Matrix at Transmit
            
            # A8.2.9 Euler Angle Rotation Matrices at Reception
            
            # Calculate the total rotation matrix at receive in the order x, y, z

            # A8.2.10 Calculate the geo-referenced lever arms at Transmit

            # A8.2.11 Calculate the geo-referenced lever arms at Reception

            # A8.2.12 Calculate the Depth Observation
            
            # A8.2.13 Determine the RP and the transducer position at the time of transmit

            # A8.2.14 Determine the RP and the transducer position at the time of receive
            
            # A8.2.15 Virtual Transducer
            
            # A8.2.16 Soundings with Respect to the Geoid

                
            # A8.2.6  Looping through the TWTTs - Leave in place
            ping += 1
            
        # A8.2.17 Soundings with Respect to the MSL
        

        

    def draw(self, **kwargs):

        # Parameters:
        # trange [t_ping_min, t_ping_max] to show on subplot 2. Default: [0, 10]
        # drange [depth_min, depth_max] to show on subplot 3. Default: [min_depth -1, max_depth -1]

        print("Drawing Positions of RP, Positioning Antenna and Transmit Transducer")
        print("Drawing Depths")

        # Depths
        # Heave removed
        depths_corr_heave = self.depth + (self.h_tx + self.h_rx) / 2  # still affected by induced heave
        # Heave removed, induced heave removed, lever arm from trans to RP applied
        depths_corr_heave_indh = depths_corr_heave + \
            (self.lever_arm_rec_rx[2, :] + self.lever_arm_trans_tx[2, :]) / 2
        # Heave removed, induced heave removed, lever arm from trans to RP applied, 
        # lever arm from RP to waterline applied
        soundings2 = depths_corr_heave_indh - self.wl_tx - self.vessel.wl

        if 'drange' in kwargs:
            depth_window = kwargs['drange']
        else:
            depth_window = [-1+min(np.nanmin(self.depth), \
                                   np.nanmin(self.sounding), \
                                   np.nanmin(depths_corr_heave),\
                                   np.nanmin(depths_corr_heave_indh), 
                                   np.nanmin(soundings2)),
                            1+max(np.nanmax(self.depth), \
                                  np.nanmax(self.sounding), \
                                  np.nanmax(depths_corr_heave), \
                                  np.nanmax(depths_corr_heave_indh), \
                                  np.nanmax(soundings2))]

        if 'trange' in kwargs:
            t_ping_min = min(kwargs['trange'])
            t_ping_max = max(kwargs['trange'])
        else:
            t_ping_min = 0
            t_ping_max = 10

        # Finding pings between t_ping_min and t_ping_max (in seconds)
        # Semme's suggestion
        # Less verbose than a cycle with conditions..and probably more efficient
        t = np.array([e.timestamp() for e in self.twtt.times])
        ping_max = len(t[t - t[0] < t_ping_max])
        ping_min = len(t[t - t[0] < t_ping_min])

        ping_window = range(ping_min, ping_max)

        # PLOTTING STARTS HERE
        fig, ax = plt.subplots(nrows=3, ncols=1, figsize=(20, 20))

        # Getting projection information for titles and labels
        proj_str = self.pos.metadata["proj_str"]
        a = proj_str.replace('+', '').split(' ')
        proj = a[0].upper()[5:]
        zone = str(int(float(a[1].upper()[5:])))
        hemisphere = a[2].upper()[0]

        projlabel = proj + zone + hemisphere
        projunits = a[5].lower()[-1]

        # should add code to get sound speed / twtt units to determine depth units

        # Projected positions
        ax[0] = plt.subplot(3, 1, 1)
        plt.plot(self.pos_rp_tx[0, :], self.pos_rp_tx[1, :],'.')
        plt.plot(self.pos_proj_ant_tx[0, :], self.pos_proj_ant_tx[1, :],'.')
        plt.plot(self.pos_trans_tx[0, :], self.pos_trans_tx[1, :])
        plt.legend(['RP_tx', 'PosAntenna_tx', 'TransTX_tx'])
        plt.ylabel("Northing [%s]" % projunits)
        plt.xlabel("Easting [%s]" % projunits)
        plt.title("Projected positions [%s]" % projlabel, fontweight="bold")
        plt.grid(True)
        ax[0].get_xaxis().get_major_formatter().set_useOffset(False)
        ax[0].get_yaxis().get_major_formatter().set_useOffset(False)

        # Projected positions considering a time window.
        ax[1] = plt.subplot(3, 1, 2)
        plt.plot(self.pos_rp_tx[0, ping_window], self.pos_rp_tx[1, ping_window],'.')
        plt.plot(self.pos_proj_ant_tx[0, ping_window], self.pos_proj_ant_tx[1, ping_window],'.')
        plt.plot(self.pos_trans_tx[0, ping_window], self.pos_trans_tx[1, ping_window],'.')
        plt.legend(['RP_tx', 'PosAntenna_tx', 'TransTX_tx'])
        plt.ylabel("Northing [%s]" % projunits)
        plt.xlabel("Easting [%s]" % projunits)
        plt.title("Projected positions [%s]. Time window: %2.2f - %2.2fs" % \
                  (projlabel, t_ping_min, t_ping_max), fontweight="bold")
        plt.grid(True)
        ax[1].get_xaxis().get_major_formatter().set_useOffset(False)
        ax[1].get_yaxis().get_major_formatter().set_useOffset(False)

        ax[2] = plt.subplot(3, 1, 3)
        plt.plot(self.depth, label='Depths')
        plt.plot(self.sounding, label='Soundings wrt EGM08')
        plt.plot(depths_corr_heave, '-r', label='Depths w/Heave Correction', )
        plt.plot(self.sounding_wl, label='Soundings wrt Chart Datum')
        plt.legend()
        plt.ylabel("Depths [m]")
        plt.xlabel("Ping number")
        plt.title("Depths and soundings")
        plt.grid(True)
        ax[2].get_xaxis().get_major_formatter().set_useOffset(False)
        ax[2].get_yaxis().get_major_formatter().set_useOffset(False)
        plt.gca().set_ylim(depth_window)

        if depth_window[0] < depth_window[1]:
            plt.gca().invert_yaxis()

        plt.show()
        fig.tight_layout(pad=5)
        
    def draw_depths(self):
        fig=plt.figure(figsize=(12, 6))
        plt.plot(self.twtt.times, self.depth+(self.h_tx+self.h_rx)/2+self.la_trans_rec_txrx[2,:])

        plt.title('Depths [m]')
        plt.ylabel('Depths [m] →')
        plt.xlabel('Time ('+self.twtt.metadata['time_basis']+') →')
        plt.gca().invert_yaxis()
    
            

#         plt.plot(self.twtt.times, self.pos_ant
#         plt.plot(self.twtt.times, self.depth+(self.h_tx+self.h_rx)/2+self.la_trans_rec_txrx)
        