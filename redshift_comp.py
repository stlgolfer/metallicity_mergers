import h5py as h5
import numpy as np
import os
import matplotlib.pyplot as plt
import pandas as pd
import string
from compas_python_utils.cosmic_integration import FastCosmicIntegration
from scipy import stats
from compas_python_utils.cosmic_integration.selection_effects import SNRinterpolator
from stellar_types import stellar_types_dictionary
from main import RateProfile, SensitivityProfile

if __name__ == 'main':
    path = './Boesky_sims.h5'

    print('excecuting this code might take a little while (~few min) \n')
    fdata = h5.File(path)

    redshifts = fdata[rate_profile.key]['redshifts'][...].squeeze()
    w_per_z_per_system = fdata[rate_profile.key][profile.sensfile][...].squeeze()
    redshift_fig, redshift_ax = plt.subplots(1,1)
    redshift_ax.plot(redshifts[:-1], np.sum(w_per_z_per_system, axis=0)) # change "2" if you actually ran the weights for more redshifts (here its only for 0 and 1, ie beyond  >2)
    redshift_ax.set_xlabel('redshift z')
    redshift_ax.set_ylabel('Merger rate at z =0 systems/(Gpc^3 * yr)', fontsize=12)
    redshift_ax.set_title('BBH merger rate', fontsize=20)
    redshift_fig.savefig(f'./redshift{profile.name}.png')
    redshift_fig.show()