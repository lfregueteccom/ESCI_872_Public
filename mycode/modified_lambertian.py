from numpy import log10,cos,pi
import numpy as np
import numpy.matlib
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

# C5.0 Create the Function `modified_lambertian`
def modified_lambertian( ..., ..., ... ):
    bs = ...
    n_bot=len(s_b[0])

    # C5.1 Convert the Depression Angles to Incidence Angles
    ...
    
    # C5.2 Determine the Lambertian Component
    rough = s_b[1].reshape(n_bot,1)
    bs_lambertian = rough+10*log10(cos(th)**2)

    # C5.3 Determine the Specular Component
    spec = ... .reshape(...)
    crit = ...
    bs_specular = spec * (1 + cos(pi * th / crit)) / 2
    r_CA=abs(...)>=crit
    bs_specular[...]=0
    
    # C5.4 Total Backscatter Angular Response 
    bs=...+...
    
    if verbose:
        
        # Create a figure
        fig = plt.figure(figsize=(10, 12))
        
        # Add a figure Title
        fig.suptitle("Modified Lambertian Scattering")
        ax=list()
        for i in range(n_bot):

            if i == 0:
                ax.append(fig.add_subplot(n_bot,2,2*(i+1)-1))
            else:
                ax.append(fig.add_subplot(n_bot,2,2*(i+1)-1,sharex=ax[0],sharey=ax[0]))
                
            plt.plot(th[i,:]*180/pi,bs_specular[i,:],'b',linewidth=2)
        
            # Plot a grid
            plt.grid()
            
            # Add a title
            ax[-1].set_title('Specular Response: '+s_b[0][i])
            
            # Add the labels
            plt.xlabel('Angle of Incidence [deg]')
            plt.ylabel('Backscatter Strength BS [dB]')

            if i == 0:
                ax.append(fig.add_subplot(n_bot,2,2*(i+1)))
            else:
                ax.append(fig.add_subplot(n_bot,2,2*(i+1),sharex=ax[1],sharey=ax[1]))
            
            # Plot the La,bertian Scatter
            plt.plot(th[i,:]*180/pi,bs_lambertian[i,:],'g',linewidth=2)
            
            # And the modified Lambertian Scatter
            plt.plot(th[i,:]*180/pi,bs[i,:],'k',linewidth=2)
            
            # Add a legend
            plt.legend(['Lambertian','Combined'])
        
            # Plot a grid
            plt.grid()
            
            # Add a title
            ax[-1].set_title('Modified Lambertian Response: '+s_b[0][i])
            
            # Add the labels
            plt.xlabel('Angle of Incidence [deg]')
            plt.ylabel('Backscatter Strength BS [dB]')

        plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=.5)
        plt.show()
        
    return bs
