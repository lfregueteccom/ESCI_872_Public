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


# D1 Defining the Vessel - from A0.1 Class Initialization and Attributes


# D2 Representing the Transmit and Receive Array

# D2.0 Representing the Transmit and Receive Array - Reality is Not Ideal

# D3 Parameters for Beam Forming

# D4 Get the Positioning Data

# D4.0 Squaring the Positioning Data Away

# D5 Get the Motion Data

# D6 Get the Sound Speed Data

# D7 Get the Ping Data

# D7.0 Indexing the Ping Data

# D8 Timing it Right

# D9.0 Getting Orientation Vectors

# D9.1 Getting Orientation Matrices

# D10 Finding out Where We Are 

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
