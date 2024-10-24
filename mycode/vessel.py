import os
import numpy as np


# 0.0 Creating the Class Vessel
class Vessel:
    """A Class for handling Vessel Specific Data"""

    def __init__(self):
        pass # Once you add code you should remove this

        # 0.1.0 Metadata descriptive of the survey
        self.metadata = dict()
        self.metadata["name"]=str()
        self.metadata["owned_by"]=str()
        self.metadata["operated_by"]=str()
        self.metadata["loa"]=float()
        self.metadata["pos_source"]=str()
        self.metadata["sonar"]=str()
        self.metadata["mru"]=list()
        self.metadata["dist_unit"]="m"
        self.metadata["angle_unit"]="rad"      

        # 0.1.2 Adding Other Attributes to the Class
        self.lever_arm_trans =  np.array([])
        self.lever_arm_rec =  np.array([])
        self.lever_arm_pos = np.array([])
        self.lever_arm_mru = np.array([])
        self.wl=np.float()

    # 0.1.1 The String Representation Method
    def __str__(self): 
        
        txt = "Vessel Metadata: "
        for key, value in self.metadata.items():
            txt += key + ": " + str( value) + '\n'

#         txt  = "Vessel Name           : %s\n" % (self.metadata["name"])
#         txt += "Owned by              : %s\n" % (self.metadata["owned_by"])
#         txt += "Operated by           : %s\n" % (self.metadata["operated_by"])
#         txt += "Length Over All       : %.0f%s\n" % (self.metadata["loa"], self.metadata["dist_unit"])
#         txt += "Positioned system     : %s\n" % (self.metadata["pos_source"])
#         txt += "Sonar system          : %s\n" % (self.metadata["sonar"])
#         txt += "Motion Reference Unit : %s\n" % (self.metadata["mru"])
        
        return txt