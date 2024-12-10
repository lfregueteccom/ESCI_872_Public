# D0 Create Python Script - Load the Dependecies

import sys
import import_ipynb
import os.path
import matplotlib as plt
import numpy as np
from datetime import datetime, timedelta
from numpy import pi, arctan2, arccos, abs, sin, cos, tan, sqrt, sum, arctan, arcsin
from mycode.position import Position
from mycode.twtt import TWTT
from mycode.integration import Integration
from mycode.integration import Motion
from datetime import datetime, timezone
from mycode.vessel import Vessel
from mycode.ssp import SSP
from mycode.ping import Ping
from mycode.om_math import Rx, Ry, Rz
from numpy.linalg import norm


# D0.1 Add the Current Folder To the Path
sys.path.append(os.getcwd())
abs_path = os.path.abspath(os.path.curdir)

# D1 Defining the Vessel - from A0.1 Class Initialization and Attributes
vessel = Vessel()
vessel.metadata["name"]="USNS Henson"
vessel.metadata["owned_by"]="United States Navy"
vessel.metadata["operated_by"]="United States Navy"
vessel.metadata["pos_source"]="NavCom (C-Nav)"
vessel.metadata["sonar"]="Kongsberg EM710"
vessel.metadata["mru"]="Applanix POS/MV 320"
vessel.metadata["loa"]=100
vessel.wl = -2.59
vessel.lever_arm_trans = np.asarray([16.26, -1.75, 4.15]).reshape((3, 1))
vessel.lever_arm_rec = np.array([14.82, -2.01, 4.17]).reshape((3, 1))
vessel.lever_arm_pos = np.array([-5.73, -0.12, -30.00]).reshape((3, 1))
vessel.lever_arm_mru = np.array([0, 0, 0]).reshape((3, 1))


# D2 Representing the Transmit and Receive Array
tx_ideal = np.array([[1],[0],[0]])
rx_ideal = np.array([[0],[1],[0]])

# D2.0 Representing the Transmit and Receive Array - Reality is Not Ideal
ma_tx=np.array([0.127,1.024,359.957-360])*pi/180  # Tx Mount angles 
ma_rx=np.array([0.101,0.894,0.065])*pi/180    # Rx Mount angles 

# D3 Parameters for Beam Forming
ss_tx=1543.7                       
beam_select=391  
# D4 Get the Positioning Data
positions = Position()
positions.read_jhc_file(os.path.join(abs_path,"Data","Lab_A_GNSS.txt"))

# D4.0 Squaring the Positioning Data Away
positions.carto_project("UTM", "ortho")

# D5 Get the Motion Data
motions = Motion()
motions.read_jhc_file(os.path.join(abs_path,"Data","Lab_A_MRU.txt"))

# D6 Get the Sound Speed Data
ssp = SSP()
ssp.read_jhc_file(os.path.join(abs_path,"Data","Lab_A_SVP.txt"))

# D7 Get the Ping Data
ping = Ping()
ping.read(os.path.join(abs_path,"Data","data_ping_4170.txt"))

# D7.0 Indexing the Ping Data
i = ping.get_beam_index(beam_select)

# D8 Timing it Right
t_tx = ping.tx_time+ping.tx_t_offset_beam[i]
t_rx = t_tx+timedelta(0,ping.twtt[i])

# D9.0 Getting Orientation Vectors
att_tx = motions.get_motion(t_tx)
att_rx = motions.get_motion(t_rx)

# D9.1 Getting Orientation Matrices
R_tx = motions.get_rotation_matrix(t_tx)
R_rx = motions.get_rotation_matrix(t_rx)

# D10 Finding out Where We Are 
geo_tx_ant = positions.get_position(t_tx)
geo_rx_ant = positions.get_position(t_rx)
# D11 Georeferencing the Lever Arms
geo_pos_la_tx = motions.geo_reference_la(t_tx, vessel.lever_arm_pos)
geo_pos_la_rx = motions.geo_reference_la(t_rx, vessel.lever_arm_pos)

# D12 Calculating the Reference Position RP
geo_pos_rp_tx = np.zeros((3,1))
geo_pos_rp_rx = np.zeros((3,1))

geo_pos_rp_tx[0] =  geo_tx_ant[0] - geo_pos_la_tx[1]
geo_pos_rp_tx[1] =  geo_tx_ant[1] - geo_pos_la_tx[0]
geo_pos_rp_tx[2] = -geo_tx_ant[2] - geo_pos_la_tx[2]

geo_pos_rp_rx[0] =  geo_rx_ant[0] - geo_pos_la_rx[1]
geo_pos_rp_rx[1] =  geo_rx_ant[1] - geo_pos_la_rx[0]
geo_pos_rp_rx[2] =  -geo_rx_ant[2] - geo_pos_la_rx[2]

# D13 The Separation of Shot and Receive
geo_pos_diff = geo_pos_rp_rx - geo_pos_rp_tx

# D14 Course Over Ground
cog = arctan2(geo_pos_diff[0],geo_pos_diff[1])

temp = geo_pos_diff.copy()
geo_pos_diff[0] = temp[1]
geo_pos_diff[1] = temp[0]

# D15 Drift Angle
da=(cog - att_tx[2])[0]
while da < -pi:
     da += 2 * pi
while da > pi:
    da -= 2 * pi
    
# D16 Heading Change
d_h=att_rx[2]-att_tx[2]

# D17.0 Aligning the 1-Axis
ref_pos_diff = Rz(-att_tx[2])@geo_pos_diff

# D17.1 New Transmit Transducer Lever Arm
new_la_tx = Ry(att_tx[1])@Rx(att_tx[0])@vessel.lever_arm_trans

# D17.2 New Receive Transducer Lever Arm
new_la_rx = Rz(d_h)@Ry(att_rx[1])@Rx(att_rx[0])@vessel.lever_arm_rec + ref_pos_diff

# D17.3 Virtual Transducer location; Why colocation?
new_virtual_trans = (new_la_rx + new_la_tx)/2

# D17.4 Virtual Transducer location: Vertical Placement
virtual_trans = new_virtual_trans.copy()
virtual_trans[2] += (att_tx[3] + att_rx[3]) / 2 - vessel.wl

# D18.1 Align Transducer Arrays to the Vessel Reference Frame
tx_mount = Rz(ma_tx[2])@Ry(ma_tx[1])@Rx(ma_tx[0])@tx_ideal
rx_mount = Rz(ma_rx[2])@Ry(ma_rx[1])@Rx(ma_rx[0])@rx_ideal

# D.18.2 Correct for Orientation at Transmit
tx=Rz(0)@Ry(att_tx[1])@Rx(att_tx[0])@tx_mount

# D.18.3 Correct for Orientation at Receive
rx=Rz(d_h)@Ry(att_rx[1])@Rx(att_rx[0])@rx_mount

# D18.4 Calculate the Non-Orthogonality angle non_ortho
non_ortho=-arcsin(float(tx.T@rx))

# D18.5 Create a New Orthonormal Basis XYZ'
xp=tx.copy()
zp=np.cross(tx.T,rx.T).T
zp/=norm(zp)
yp=np.cross(zp.T,xp.T).T
yp/=norm(yp)
Tp = np.hstack((xp, yp, zp)) 

# D18.6.0  Formula for Intersection of Non-Orthogonal Arrays
y1=sin(-ping.steer_rx[i]) / cos(non_ortho)
y2=sin(ping.steer_tx[i]) * tan(non_ortho)
radial=sqrt((y1 + y2) ** 2 + sin(ping.steer_tx[i]) ** 2)

# D18.6.1  Formulate the Beam Vector `bv_p` in XYZ'
bv_p=np.array([[sin(ping.steer_tx[i])],[y1 + y2],[sqrt(1-radial**2)]])

# D18.7 Transform the Beam Vector to Geo Referenced Space
bv_g=Tp@bv_p 

# D.18.8 Determine the Depression Angle of the Beam Vector
beam_th=float(arctan(bv_g[2] / sqrt(sum(bv_g[0:2] ** 2))))

# D18.9 Determine the Azimuth of the Beam Vector
beam_az=float(arctan2(bv_g[1],bv_g[0]))

# D19 At last: Ray Tracing
depth,hor_dist,_,_ = ssp.ray_trace_twtt(virtual_trans[2],beam_th,ss_tx,ping.twtt[i])

# D20 Positioning the Bottom Strike Relative to the Virtual Transducer
bot_vx = np.zeros((3,1))
bot_vx[0] = virtual_trans[0] + hor_dist * cos(beam_az)
bot_vx[1] = virtual_trans[1] + hor_dist * sin(beam_az)
bot_vx[2] = depth   

# D21 Positioning the Bottom Strike Relative to the RP at Transmit
bot_rp_tx = bot_vx.copy()
bot_rp_tx[0] -= virtual_trans[0]
bot_rp_tx[1] -= virtual_trans[1]
