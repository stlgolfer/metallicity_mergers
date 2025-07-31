# idea here is to use adam's code and produce the merger rate
# this is basically already done in Floor's code

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

class SensitivityProfile:
    def __init__(self, name: str, up_to_redshift: int, sensfile: str):
        self.name = name
        self.up_to_redshift = up_to_redshift
        self.sensfile = sensfile

class RateProfile:
    def __init__(self, key: str):
        self.key = key

        # now let's see if we can dynamically append all the properties
        properties = key.split('_')[1:] # skip first one
        for p in properties:
            number = ''
            name = ''
            for char in p:
                if char.isdigit() or (char == '.') or (char == '-'):
                    number = number + char
                else:
                    name = name + char
            print(f'setting {p.replace(str(float(number)), '')} to {float(number)}')
            
            # get the name by subtracting the number from the rest of the string
            setattr(self, name, float(number))

# to obtain properties of ALL binaries simulated, do this:
path = './Boesky_sims.h5'

print('excecuting this code might take a little while (~few min) \n')
fdata = h5.File(path)
# shows the different files within the hdf5 folder 

print('the available datasets for this file are:')
print(fdata.keys())

#
# z_index = 0
#
# print('the redshift weights per system at z = ', redshifts[z_index], ' are given by' )
# w_z_index =w_per_z_per_system[:,z_index]
# want to sum up to a certain redshift. For now, sum up to redshift 1

# fdata['BSE_System_Parameters']['Metallicity@ZAMS(1)'] has many more elements than BSE_DCO
# we need to get the seeds from dcomask
# then we get the metallicity from system_params in the spots where the seed is equal

#select BHBH seeds
#  np.where((st1[()] == 13) & (st2[()]==13))[0].shape
# np.where(fdata['BSE_Double_Compact_Objects']['Stellar_Type(2)'][dcomask][()] == 13)
# ^^^^ the above will work since apparently fci doesn't use every dco system, so still have to mask even in dco
# proof: np.isin(fdata['BSE_Double_Compact_Objects']['SEED'], fdata[rate_key]['SEED']) should be all true, but it's not
# ahh it looks like it just pulling BHBH, make sure to set to 'ALL' for dco_type
# actually, instead of writing code to select things after the fact, we can just re-run FCI and select the
# object of interest

# dcomask = fdata[rate_key]['DCOmask'][()]
# fs1 = fdata['BSE_Double_Compact_Objects']['Stellar_Type(1)'][dcomask][()]
# fs2 = fdata['BSE_Double_Compact_Objects']['Stellar_Type(2)'][dcomask][()] # final stellar type
# object_select = np.where((fs1 == 14) & (fs2 == 14)) # look for BHBH. two ways for this to happen
# so must have logical OR
# change this key for the detection rate

# print(np.argwhere(stellar_types_1 == 1))
# print(set(stellar_types_1))
# print(set(stellar_types_2))

# now we would like to overlay the histogram but do it for each type of stellar object
# ie we want to group the rates by progenitors (would be nice to do a fill on this)
# use this as a selector Stellar_Type@ZAMS(2)
# _, bins = np.histogram(m1zams, bins=100, density=True)
# m1zamskde = stats.gaussian_kde(m1zams, weights=w_z_index)
# plt.plot(bins[:-1], m1zamskde(bins[:-1]))
def plot_up_to_redshift(rate_profile: RateProfile, ax, profile: SensitivityProfile):
    dcoseeds = fdata[rate_profile.key]['SEED'][()] # instead of selecting all dcos, get the seeds that are BHBH

    m1zams = fdata['BSE_System_Parameters']['Metallicity@ZAMS(1)'][()]
    search = np.isin(fdata['BSE_System_Parameters']['SEED'], dcoseeds) # find the location of each dco in the larger colleciton
    m1zams = m1zams[search]
    stellar_types_1 = fdata['BSE_System_Parameters']['Stellar_Type@ZAMS(1)'][()][search]
    stellar_types_2 = fdata['BSE_System_Parameters']['Stellar_Type@ZAMS(2)'][()][search]
    # rate_key = 'Rates_mu00.035_muz-0.23_alpha0.0_sigma00.39_sigmaz0.0'
    print(fdata[rate_profile.key].keys())
    redshifts                 = fdata[rate_profile.key]['redshifts'][()] # Redshifts at which the rates were calculated
    # merger_rate_z0_per_system = fdata[rate_key]['detection_rateO3'][()] # detection rate for O3 sensitivity
    # total_rate_z0 = np.sum(merger_rate_z0_per_system) # this is the rate for redshift 0, you can get the rate for all redshifts by summing over merger_rate

    # print some interesting information about the various merger rates per redshift
    # print(total_rate_z0, '[Gpc^-3 yr^_1]')

    w_0 = fdata[rate_profile.key]['merger_rate_z0'][...].squeeze()

    print(w_0)
    redshifts = fdata[rate_profile.key]['redshifts'][...].squeeze()
    print('available redshifts are: ', redshifts, ' this gives %s options'%len(redshifts))

    w_per_z_per_system = fdata[rate_profile.key][f'detection_rate{profile.sensfile}'][...].squeeze()
    r_index = np.where(redshifts == profile.up_to_redshift)
    if len(r_index) == 0:
        raise "Redshift not found..consider changing your FastCosmicIntegration settings"
    r_index = r_index[0][0] # redshifts are monotone increasing so first result is ok
    w_z_summed = np.sum(
        w_per_z_per_system[:,:r_index],
        axis=1
    )

    # create a subroutine that iterates through each star 1 type for now
    def plot_stellar_type_at_zams(type_index, include_histo=False, include_complete_redshift=False):
        bins=1000
        if include_complete_redshift:
            redshift_fig, redshift_ax = plt.subplots(1,1)
            redshift_ax.plot(redshifts[:-1], np.sum(w_per_z_per_system, axis=0)) # change "2" if you actually ran the weights for more redshifts (here its only for 0 and 1, ie beyond  >2)
            redshift_ax.set_xlabel('redshift z')
            redshift_ax.set_ylabel('Merger rate at z =0 systems/(Gpc^3 * yr)', fontsize=12)
            redshift_ax.set_title('BBH merger rate', fontsize=20)
            redshift_fig.savefig(f'./redshift{profile.name}.png')
            redshift_fig.show()
        # stellar_search = np.argwhere(stellar_types_1==type_index)
        # if stellar_search.size <= 2:
        #     print(f'Stellar type {stellar_types_dictionary[type_index]} has <=2 elements, so not painting')
        #     return
        if include_histo:
            ax.hist(
                m1zams, #stellar_search
                bins=bins,
                weights=w_z_summed,
                # label=f'(1) Type {stellar_types_dictionary[type_index]}', histtype='step',
                density=True
            )

        _, bins = np.histogram(m1zams, bins=bins)
        # print(m1zams[stellar_search].shape)
        # print(f'w_z_summed: {w_z_summed[stellar_search].flatten().shape}')
        m1zamskde = stats.gaussian_kde(
            m1zams.flatten(),
            weights=w_z_summed.flatten()
        )
        ax.plot(
            bins[:-1]/0.012, m1zamskde(bins[:-1]),
            label=profile.name
            # label=f'(1) Type {stellar_types_dictionary[type_index]} ({detector})'
        )
        ax.fill_between(
            bins[:-1]/0.012,
            m1zamskde(bins[:-1]),
            interpolate=True,
            alpha=0.3
        )
    # plt.hist(m1zams[np.argwhere(stellar_types_1==16)], bins=100, weights=w_z_index[stellar_search], legend='(1) Type 1')
    # for t in set(stellar_types_1):
    #     plot_stellar_type_at_zams(t)
    plot_stellar_type_at_zams(-1, include_complete_redshift=True)
# plot_stellar_type_at_zams(list(set(stellar_types_1))[0])
fig, ax = plt.subplots(1, 1)

plot_up_to_redshift(
    RateProfile('Rates_mu00.035_muz-0.23_alpha0.0_sigma00.39_sigmaz0.0'),
    ax,
    SensitivityProfile("CE",
                       10, sensfile='CE.txt'
                       )
)

plot_up_to_redshift(
    RateProfile('Rates_mu00.035_muz-0.23_alpha0.0_sigma00.39_sigmaz0.0'),
    ax,
    SensitivityProfile("O3",
                       1, sensfile='O3'
                       )
)
# plot_up_to_redshift(ax, 'CE')
# we eventually want a way to plot this rate for a few different mu0s and sigmas. from van son we want to try a few values
# there are 9 total calculations to complete
# three min values for mu0
# three max values for muz
# sigmas seem to be fixed at 0.036 and 0.006, respectively. should be wary of units

ax.set_title(f'Number of DCO systems/year')
ax.set_xlabel('Metallicity1 at ZAMS Z/Z0') #TODO: check units against COMPAS simulation. Looks like this is just Z
ax.set_ylabel(r'Number of DCO systems/year [simulation weighted]')
ax.set_xscale('log')
ax.legend()
fig.show()

# taking from the example, we also want to be able to get the "primary black hole mass" histogram,
# which from the examples looks like we want to get the more massive one of the two
# will need to get the SEEDs from the DCOs that were redshifted, then select each mass
# dcomask = fdata[rate_key]['DCOmask'][()] # to be able to get seed list
# mass1 = fdata['BSE_Double_Compact_Objects']['Mass(1)'][()][dcomask]
# mass2 = fdata['BSE_Double_Compact_Objects']['Mass(2)'][()][dcomask]
# M_moreMassive = np.maximum(mass1, mass2)

# plt.figure()
# _, bins = np.histogram(M_moreMassive, bins=100, density=True)
# plt.hist(M_moreMassive, bins=100, density=True, weights=w_z_index)
# # now do kde as well
# kde = stats.gaussian_kde(M_moreMassive, weights=w_z_index)
#
# # plt.plot(bins[:-1], counts)
# plt.plot(bins[:-1], kde(bins[:-1]))
# plt.xlabel('BH mass [Msun] ')
# plt.ylabel(r'BBH merger rate at $z =%s [\rm{Gpc}^{-3} \rm{yr}^{-1}]$'%np.round(redshifts[z_index],3), fontsize=12)
# # plt.ylabel(r'Normalized Probability')
# plt.title('BH masses of BH-BH merger population at redshift %s Density'%np.round(redshifts[z_index],3), fontsize=12)
# plt.show()
#

fdata.close()
# then run the cosmic integrator
# !python3 FastCosmicIntegration_large_files.py  --mu0 0.035 --muz -0.23 --sigma0 0.39 --sigmaz 0.0 --alpha 0.0 --weight mixture_weight --zstep 0.5 --sens O3 --m1min 10. --aSF 0.01 --bSF 2.77 --cSF 2.9 --dSF 4.7
# python ~/Documents/Github/COMPAS/compas_python_utils/cosmic_integration/FastCosmicIntegration.py --path ./Boesky_sims.h5 --mu0 0.035 --muz -0.23 --sigma0 0.39 --sigmaz 0.0 --alpha 0.0 --weight mixture_weight --zstep 0.05 --sens O3 --m1min 10. --aSF 0.01 --bSF 2.77 --cSF 2.9 --dSF 4.7 --maxz 20 --maxzdet 20