import os
from datetime import datetime, timezone
from numpy import pi, cos, sin, log, exp
import numpy as np
import matplotlib.pyplot as plt

class AnalyzeSS:
    """A Class for handling Sound Speed Profile data"""
    
    def __init__(self):
        self.SSPs=list()
        
    # B5.0  Plotting the Profiles
    def draw(self):
        # B5.0.0  Determine the Number of Profiles to be Plotted
        n_ssps = len(self.SSPs)
        
        # B5.0.1 Creating the Figure
        fig = plt.figure(figsize=(18, 6))

        # B5.0.2 Set Appropriate Scaling
        # Initialize the ranges
        min_ss = 10000
        max_ss = 0
        max_dep = 0

        # Update the ranges
        for i in range(n_ssps):
            if self.SSPs[i].proc_depth[-2] > max_dep:
                max_dep = self.SSPs[i].proc_depth[-2]
            if min(self.SSPs[i].proc_ss[0:-1]) < min_ss:
                min_ss = min(self.SSPs[i].proc_ss[0:-2])
            elif max(self.SSPs[i].proc_ss) > max_ss:
                max_ss = max(self.SSPs[i].proc_ss[0:-2])

        # Set the ranges as tuples
        dep_range = (0, max_dep+max_dep/20)
        ss_range = (min_ss-3,max_ss+3)
        
        # B5.0.3 Create the Axes for Each SSP and Draw the Data
        ax = plt.subplot(1,n_ssps,1)
        self.SSPs[0].draw(False, ax, dep_range, ss_range, True)
        for i in range(1,n_ssps):
            ax = plt.subplot(1,n_ssps,i+1)
            self.SSPs[i].draw(False, ax, dep_range, ss_range, False)

        plt.show()
    
    # B5.1.3 Select a Profile for Generating Depths from TWTTs
    def select_ssp(self, file_name):
        return 'Not found'
        
    # 5.1.2 Synthesize TWTTs for Given Depth   
    def synthesize_twtt(self, th, tx_depth, depth):
        n_ssps = len(self.SSPs)
        n_beams = len(th)
        twtt = np.zeros((n_ssps,n_beams))
        cross_track_d = np.zeros((n_ssps,n_beams))
        for i in range(n_ssps):
            ss_start = self.SSPs[i].determine_c(tx_depth)
            for j in range(n_beams):
                twtt[i,j],cross_track_d[i,j],_,_= \
                    self.SSPs[i].determine_twtt(tx_depth, th[j], ss_start, depth)

        return twtt, cross_track_d
    
    def select_ssp(self, file_name):
        for i, ssp in enumerate(self.SSPs):
            if os.path.basename(ssp.metadata['name']) == file_name:
                return ssp, i
            i += 1
        return 'Not found'
        
        
        
        
        
        