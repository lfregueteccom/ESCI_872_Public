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
ping.read(os.path.join(abs_path,"Data","DATA_PING_4170.txt"))

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

# D12 Calculating the Reference Position RP

# D13 The Separation of Shot and Receive

# D14 Course Over Ground

# D15 Drift Angle
    
# D16 Heading Change

# D17.0 Aligning the 1-Axis

# D17.1 New Transmit Transducer Lever Arm

# D17.2 New Receive Transducer Lever Arm

# D17.3 Virtual Transducer location; Why colocation?

# D17.4 Virtual Transducer location: Vertical Placement

# D18.1 Align Transducer Arrays to the Vessel Reference Frame

# D.18.2 Correct for Orientation at Transmit

# D.18.3 Correct for Orientation at Receive

# D18.4 Calculate the Non-Orthogonality angle non_ortho

# D18.5 Create a New Orthonormal Basis XYZ'

# D18.6.0  Formula for Intersection of Non-Orthogonal Arrays

# D18.6.1  Formulate the Beam Vector `bv_p` in XYZ'

# D18.7 Transform the Beam Vector to Geo Referenced Space

# D.18.8 Determine the Depression Angle of the Beam Vector

# D18.9 Determine the Azimuth of the Beam Vector

# D19 At last: Ray Tracing

# D20 Positioning the Bottom Strike Relative to the Virtual Transducer

# D21 Positioning the Bottom Strike Relative to the RP at Transmit
