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

stellar_types_dictionary = ["Main_Sequence_<=_0.7",
                 "Main_Sequence_>_0.7",
                 "Hertzsprung_Gap",
                 "First_Giant_Branch",
                 "Core_Helium_Burning",
                 "Early_Asymptotic_Giant_Branch",
                 "Thermally_Pulsing_Asymptotic_Giant_Branch",
                 "Naked_Helium_Star_MS" ,
             "Naked_Helium_Star_Hertzsprung_Gap" ,
    "Naked_Helium_Star_Giant_Branch" ,
                            "Helium_White_Dwarf" ,
                     "Carbon-Oxygen_White_Dwarf" ,
                       "Oxygen-Neon_White_Dwarf" ,
                                  "Neutron_Star" ,
                                    "Black_Hole" ,
                              "Massless_Remnant" ,
                        "Chemically_Homogeneous" ,
                                          "Star" ,
                                   "Binary_Star" ,
    "Not_a_Star!"]

# to obtain properties of ALL binaries simulated, do this:
path = './Boesky_sims.h5'

print('excecuting this code might take a little while (~few min) \n')
fdata = h5.File(path)
# shows the different files within the hdf5 folder 

print('the available datasets for this file are:')
print(fdata.keys())
# hopefully there are some DCOs in the file to work with
# print(fdata['BSE_Double_Compact_Objects'].keys())

rate_key = 'Rates_mu00.035_muz-0.23_alpha0.0_sigma00.39_sigmaz0.0'
print(fdata[rate_key].keys())
redshifts                 = fdata[rate_key]['redshifts'][()] # Redshifts at which the rates were calculated
merger_rate_z0_per_system = fdata[rate_key]['merger_rate_z0'][()] # detection rate for O3 sensitivity 
total_rate_z0 = np.sum(merger_rate_z0_per_system) # this is the rate for redshift 0, you can get the rate for all redshifts by summing over merger_rate

# print some interesting information about the various merger rates per redshift
print(total_rate_z0, '[Gpc^-3 yr^_1]')

w_0 = fdata[rate_key]['merger_rate_z0'][...].squeeze()

print(w_0)
redshifts = fdata[rate_key]['redshifts'][...].squeeze()
print('available redshifts are: ', redshifts, ' this gives %s options'%len(redshifts))

w_per_z_per_system = fdata[rate_key]['merger_rate'][...].squeeze()

print(np.shape(w_per_z_per_system))

z_index = 0

print('the redshift weights per system at z = ', redshifts[z_index], ' are given by' )
w_z_index =w_per_z_per_system[:,z_index]
# want to sum up to a certain redshift. For now, sum up to redshift 1
w_z_summed = np.sum(w_per_z_per_system[:,:20], axis=1) # here we crudely sum over all redshift weights, but
# eventually we'd like to be able to sum up to a certain redshift
print(w_z_index.shape)

# may not be able to directly mask from dcomask to get metallicities
# fdata['BSE_System_Parameters']['Metallicity@ZAMS(1)'] has many more elements than BSE_DCO
# we need to get the seeds from dcomask
# then we get the metallicity from system_params in the spots where the seed is equal
dcoseeds = fdata[rate_key]['SEED'][()]
m1zams = fdata['BSE_System_Parameters']['Metallicity@ZAMS(1)'][()]
search = np.isin(fdata['BSE_System_Parameters']['SEED'], dcoseeds) # find the location of each dco in the larger colleciton
m1zams = m1zams[search]
stellar_types_1 = fdata['BSE_System_Parameters']['Stellar_Type@ZAMS(1)'][()][search]
stellar_types_2 = fdata['BSE_System_Parameters']['Stellar_Type@ZAMS(2)'][()][search]

# print(np.argwhere(stellar_types_1 == 1))
# print(set(stellar_types_1))
# print(set(stellar_types_2))

plt.figure()
# now we would like to overlay the histogram but do it for each type of stellar object
# ie we want to group the rates by progenitors (would be nice to do a fill on this)
# use this as a selector Stellar_Type@ZAMS(2)
# _, bins = np.histogram(m1zams, bins=100, density=True)
# m1zamskde = stats.gaussian_kde(m1zams, weights=w_z_index)
# plt.plot(bins[:-1], m1zamskde(bins[:-1]))

# create a subroutine that iterates through each star 1 type for now
def plot_stellar_type_at_zams(type_index):
    bins=1000
    stellar_search = np.argwhere(stellar_types_1==type_index)
    # plt.hist(
    #     m1zams[stellar_search],
    #     bins=bins,
    #     weights=w_z_summed[stellar_search],
    #     label=f'(1) Type {stellar_types_dictionary[type_index]}', histtype='step',
    #     density=True
    # )

    _, bins = np.histogram(m1zams[stellar_search], bins=bins)
    # print(m1zams[stellar_search].shape)
    # print(f'w_z_summed: {w_z_summed[stellar_search].flatten().shape}')
    m1zamskde = stats.gaussian_kde(
        m1zams[stellar_search].flatten(),
        weights=w_z_summed[stellar_search].flatten()
    )
    plt.plot(bins[:-1], m1zamskde(bins[:-1]), label=f'(1) Type {stellar_types_dictionary[type_index]}')
# plt.hist(m1zams[np.argwhere(stellar_types_1==16)], bins=100, weights=w_z_index[stellar_search], legend='(1) Type 1')
for t in set(stellar_types_1):
    plot_stellar_type_at_zams(t)
# plot_stellar_type_at_zams(list(set(stellar_types_1))[0])
plt.title('Number of systems/year [Density] summed to z=0')
plt.xlabel('Metallicity1 at ZAMS') #TODO: check units against COMPAS simulation. Is this actually Z/Z0?
plt.ylabel(r'BHNS merger rate [simulation weighted]')
plt.xscale('log')
plt.legend()
plt.show()
assert 0

# taking from the example, we also want to be able to get the "primary black hole mass" histogram,
# which from the examples looks like we want to get the more massive one of the two
# will need to get the SEEDs from the DCOs that were redshifted, then select each mass
dcomask = fdata[rate_key]['DCOmask'][()] # to be able to get seed list
mass1 = fdata['BSE_Double_Compact_Objects']['Mass(1)'][()][dcomask]
mass2 = fdata['BSE_Double_Compact_Objects']['Mass(2)'][()][dcomask]
M_moreMassive = np.maximum(mass1, mass2)

plt.figure()
_, bins = np.histogram(M_moreMassive, bins=100, density=True)
plt.hist(M_moreMassive, bins=100, density=True, weights=w_z_index)
# now do kde as well
kde = stats.gaussian_kde(M_moreMassive, weights=w_z_index)

# plt.plot(bins[:-1], counts)
plt.plot(bins[:-1], kde(bins[:-1]))
plt.xlabel('BH mass [Msun] ')
plt.ylabel(r'BBH merger rate at $z =%s [\rm{Gpc}^{-3} \rm{yr}^{-1}]$'%np.round(redshifts[z_index],3), fontsize=12)
# plt.ylabel(r'Normalized Probability')
plt.title('BH masses of BH-BH merger population at redshift %s Density'%np.round(redshifts[z_index],3), fontsize=12)
plt.show()

plt.plot(redshifts[:-1], np.sum(w_per_z_per_system, axis=0)) # change "2" if you actually ran the weights for more redshifts (here its only for 0 and 1, ie beyond  >2)
plt.xlabel('redshift z')
plt.ylabel(r'BHBH merger rate at $z =%s [\rm{Gpc}^{-3} \rm{yr}^{-1}]$'%np.round(redshifts[z_index],3), fontsize=12)
plt.title('BBH merger rate', fontsize=20)
plt.show()

fdata.close()
# then run the cosmic integrator
# !python3 FastCosmicIntegration_large_files.py  --mu0 0.035 --muz -0.23 --sigma0 0.39 --sigmaz 0.0 --alpha 0.0 --weight mixture_weight --zstep 0.5 --sens O3 --m1min 10. --aSF 0.01 --bSF 2.77 --cSF 2.9 --dSF 4.7