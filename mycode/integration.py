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
        self.pos = pos
        self.vessel = vessel
        self.twtt = twtt
        self.ssp = sound_speed_profile
        
        # A8.2.1 Determine the number of pings
        n_twtt_times = len(twtt.times)

        # A8.2.2 Memory Allocation
        R_tx = list()
        R_rx = list()
        self.lever_arm_pos_tx=np.zeros([3,n_twtt_times])
        self.lever_arm_pos_rx=np.zeros([3,n_twtt_times])
        self.lever_arm_trans_tx=np.zeros([3,n_twtt_times])
        self.lever_arm_rec_rx=np.zeros([3,n_twtt_times])
        self.pos_rp_tx=np.zeros([3,n_twtt_times])
        self.pos_rp_rx=np.zeros([3,n_twtt_times])
        self.pos_trans_tx=np.zeros([3,n_twtt_times])
        self.pos_rec_rx=np.zeros([3,n_twtt_times])
              
        # A8.2.3 Converting the `datetime` objects to POSIX Times
        t_twtt = np.array([e.timestamp() for e in twtt.times])
        t_pos = np.array([e.timestamp() for e in pos.times])
        t_mru = np.array([e.timestamp() for e in motions.times])
        t_wl = np.array([e.timestamp() for e in water_levels.times])

        # A8.2.4.0 Transmit event interpolation
        
        self.p_tx = np.interp(t_twtt, t_mru, motions.pitch)    
        self.r_tx = np.interp(t_twtt, t_mru, motions.roll)
        self.y_tx = np.interp(t_twtt, t_mru, motions.yaw)
        self.h_tx = np.interp(t_twtt, t_mru, motions.heave)
        self.wl_tx = np.interp(t_twtt, t_wl, water_levels.water_levels)

        # Determine the interpolation function for the positions
        f=interp1d(t_pos,pos.proj_pos,bounds_error=False)
        self.pos_proj_ant_tx=f(t_twtt)
        
        # A8.2.4.1  Receive time determination
        t_twtt += twtt.twtts

        # A8.2.4.2 Receive Event interpolation
        self.p_rx = np.interp(t_twtt, t_mru, motions.pitch)    
        self.r_rx = np.interp(t_twtt, t_mru, motions.roll)
        self.y_rx = np.interp(t_twtt, t_mru, motions.yaw)
        self.h_rx = np.interp(t_twtt, t_mru, motions.heave)
        self.wl_rx = np.interp(t_twtt, t_wl, water_levels.water_levels)
        self.pos_proj_ant_rx=f(t_twtt) #interpolation with the reception time

        # A8.2.5  Allocating Memory for Processed Data
        self.depth = np.zeros((n_twtt_times))
        self.sounding = np.zeros((n_twtt_times))
        self.virtual_txrx = np.zeros((3,n_twtt_times))
        
        # A8.2.6  Looping through the TWTTs
        ping = 0
        for t in t_twtt:
            # A8.2.7 Calculate the Eucliden Euler Angle Rotation Matrices
            Rx_tx = np.array([[1, 0,                     0                   ],
                              [0, cos(self.r_tx[ping]), -sin(self.r_tx[ping])],
                              [0, sin(self.r_tx[ping]),  cos(self.r_tx[ping])]])
            
            Ry_tx = np.array([[cos(self.p_tx[ping]), 0, sin(self.p_tx[ping]) ],
                              [0, 1,                     0                   ],
                              [-sin(self.p_tx[ping]), 0, cos(self.p_tx[ping])]])
            
            Rz_tx = np.array([[cos(self.y_tx[ping]), -sin(self.y_tx[ping]),0  ],
                              [sin(self.y_tx[ping]), cos(self.y_tx[ping]), 0  ],
                              [0, 0,                     1                   ]])
            
            

            # A8.2.8 Compound Rotation Matrix at Transmit
            R_1 = Rz_tx.copy()
            R_2 = Ry_tx.copy()
            R_3 = Rx_tx.copy()
            R_tx.append( R_3@R_2@R_1)
            
            # A8.2.9 Euler Angle Rotation Matrices at Reception
            Rx_rx = np.array([[1, 0,                     0                   ],
                        [0, cos(self.r_rx[ping]), -sin(self.r_rx[ping])],
                        [0, sin(self.r_rx[ping]),  cos(self.r_rx[ping])]])
            
            Ry_rx = np.array([[cos(self.p_rx[ping]), 0, sin(self.p_rx[ping]) ],
                              [0, 1,                     0                   ],
                              [-sin(self.p_rx[ping]), 0, cos(self.p_rx[ping])]])
            
            Rz_rx = np.array([[cos(self.y_rx[ping]), -sin(self.y_rx[ping]),0  ],
                              [sin(self.y_rx[ping]), cos(self.y_rx[ping]), 0  ],
                              [0, 0,                     1                   ]])
            
            # Calculate the total rotation matrix at receive in the order x, y, z
            R_1 = Rx_rx.copy()
            R_2 = Ry_rx.copy()
            R_3 = Rz_rx.copy()
            R_rx.append( R_3@R_2@R_1)

            # A8.2.10 Calculate the geo-referenced lever arms at Transmit
            self.lever_arm_pos_tx[:,[ping]]=R_tx[ping]@vessel.lever_arm_pos
            self.lever_arm_trans_tx[:,[ping]] = R_tx[ping]@vessel.lever_arm_trans

            # A8.2.11 Calculate the geo-referenced lever arms at Reception
            self.lever_arm_pos_rx[:,[ping]]=R_rx[ping]@vessel.lever_arm_pos
            self.lever_arm_rec_rx[:,[ping]] = R_rx[ping]@vessel.lever_arm_rec

            # A8.2.12 Calculate the Depth Observation
            self.depth[ping] =  np.mean(sound_speed_profile)*twtt.twtts[ping]/2
            
            # A8.2.13 Determine the RP and the transducer position at the time of transmit
            self.pos_rp_tx[0,[ping]]=self.pos_proj_ant_tx[0,[ping]]-self.lever_arm_pos_tx[1,[ping]]
            self.pos_rp_tx[1,[ping]]=self.pos_proj_ant_tx[1,[ping]]-self.lever_arm_pos_tx[0,[ping]]
            self.pos_rp_tx[2,[ping]]=self.pos_proj_ant_tx[2,[ping]]+self.lever_arm_pos_tx[2,[ping]]
            self.pos_trans_tx[0,[ping]]=self.pos_rp_tx[0,[ping]]+self.lever_arm_trans_tx[1,[ping]]
            self.pos_trans_tx[1,[ping]]=self.pos_rp_tx[1,[ping]]+self.lever_arm_trans_tx[0,[ping]]
            self.pos_trans_tx[2,[ping]]=self.pos_rp_tx[2,[ping]]-self.lever_arm_trans_tx[2,[ping]]

            # A8.2.14 Determine the RP and the transducer position at the time of receive
            self.pos_rp_rx[0,[ping]]=self.pos_proj_ant_rx[0,[ping]]-self.lever_arm_pos_rx[1,[ping]]
            self.pos_rp_rx[1,[ping]]=self.pos_proj_ant_rx[1,[ping]]-self.lever_arm_pos_rx[0,[ping]]
            self.pos_rp_rx[2,[ping]]=self.pos_proj_ant_rx[2,[ping]]+self.lever_arm_pos_rx[2,[ping]]
            self.pos_rec_rx[0,[ping]]=self.pos_rp_rx[0,[ping]]+self.lever_arm_rec_rx[1,[ping]]
            self.pos_rec_rx[1,[ping]]=self.pos_rp_rx[1,[ping]]+self.lever_arm_rec_rx[0,[ping]]
            self.pos_rec_rx[2,[ping]]=self.pos_rp_rx[2,[ping]]-self.lever_arm_rec_rx[2,[ping]]
            
            # A8.2.15 Virtual Transducer
            self.virtual_txrx[:,[ping]] = (self.pos_trans_tx[:,[ping]]+self.pos_rec_rx[:,[ping]])/2
            
            # A8.2.16 Soundings with Respect to the Geoid
            self.sounding[ping] = self.depth[ping] - self.virtual_txrx[2,[ping]]    
    
            # A8.2.6  Looping through the TWTTs - Leave in place
            ping += 1
            
        # A8.2.17 Soundings with Respect to the MSL
        self.sounding_wl = np.copy(self.depth)
        self.sounding_wl += (self.h_tx + self.h_rx) / 2 + \
            (self.lever_arm_rec_rx[2, :] + self.lever_arm_trans_tx[2, :]) / 2 + \
            - (self.wl_tx + self.wl_rx)/2 - self.vessel.wl

        

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
        