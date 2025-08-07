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
from main_vis import RateProfile, SensitivityProfile

if __name__ == '__main__':
    path = '/Volumes/Elements/Boesky_sims.h5'

    print('excecuting this code might take a little while (~few min) \n')
    fdata = h5.File(path)
    rate_profile = RateProfile('Rates_mu00.035_muz0.035_alpha-1.778_sigma01.122_sigmaz0.049')
    redshift_fig, redshift_ax = plt.subplots(1,1)

    def plot_profile(profile: SensitivityProfile):
        redshifts = fdata[rate_profile.key]['redshifts'][...].squeeze()
        w_per_z_per_system = fdata[rate_profile.key][profile.sensfile][...].squeeze()
        print(redshifts)

        r_index = np.where(redshifts == profile.up_to_redshift)
        if len(r_index) == 0:
            raise "Redshift not found..consider changing your FastCosmicIntegration settings"
        r_index = r_index[0][0]
        # print(r_index)

        # Make a merger rate vs redshift at z=0 that has
        # 1) no filter
        # 2) no filter + manual cut consistent with detector
        # 3) CE
        # r_pad = np.pad(redshifts[:-1][:r_index], (0,20-r_index))
        truncated = np.pad(np.sum(w_per_z_per_system, axis=0)[:r_index],(0,len(redshifts)-r_index-1))
        print(f'Max for {profile.name} is {np.max(truncated)}')
        # print(truncated)
        # redshift_ax.plot(redshifts[:-1], np.sum(w_per_z_per_system, axis=0))
        redshift_ax.plot(redshifts[:-1], truncated, label=f'{profile.name} up to {profile.up_to_redshift}')

    plot_profile(profile = SensitivityProfile("Nom.",
                                10,
                                sensfile="merger_rate"
                                ))
    plot_profile(profile = SensitivityProfile("Nom.",
                                1,
                                sensfile="merger_rate"
                                ))
    plot_profile(
        profile= SensitivityProfile(
            "design",
            10,
            sensfile='detection_ratedesign'
        )
    )
    plot_profile(
        profile= SensitivityProfile(
            "O1",
            10,
            sensfile='detection_rateO1'
        )
    )

    # plot_profile(
    #     profile= SensitivityProfile(
    #         "CE",
    #         10,
    #         sensfile='detection_rateCE.txt'
    #     )
    # )
    # plot_profile(profile= SensitivityProfile(
    #     'O3',
    #     10,
    #     sensfile='detection_rateO3'
    # ))

    # plot_profile(profile= SensitivityProfile(
    #     'CE',
    #     10,
    #     sensfile='detection_rateCE.txt'
    # ))
    redshift_ax.set_xlabel('redshift z')
    redshift_ax.set_ylabel('Merger rate at z =0 systems/(Gpc^3 * yr)', fontsize=12)
    redshift_ax.set_title('BBH merger rate', fontsize=20)
    redshift_ax.legend()
    redshift_fig.savefig(f'./redshift_comparison.png')
    redshift_fig.show()