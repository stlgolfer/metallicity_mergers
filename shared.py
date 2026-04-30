import h5py
from syntheticstellarpopconvolve.general_functions import generate_boilerplate_outputfile, extract_unit_dict, temp_dir
import astropy.units as u
import numpy as np
import os

def plot_merger_rate_from_hdf_file(fname, groupname, dco_metallicities):
    with h5py.File(
            fname, "r"
        ) as output_hdf5file:
            # could try and make an actual merger plot. go through each key, sum the values
            yield_locations = output_hdf5file[groupname]
            
            merger_ax_redshifts = []# np.zeros_like(list(yield_locations.keys()))
            
            # merger_ax_rates_total = []
            y_axis_units = u.Unit(extract_unit_dict(output_hdf5file, groupname + list(yield_locations.keys())[0])['yield']).to_string()

            # for each bin, we want a rate vs redshift
            for r, k in enumerate(yield_locations.keys()):
                # get number
                # units = u.Quantity(k).unit
                # print(units)
                # assert units == u.Gyr, "Rest of plotting code assumes Gyr for lookback time. change conv settings"
                t = u.Quantity(k).value
                # if u.Quantity(k) >= u.Quantity('14 Gyr'):
                #     warnings.warn('Skipping a data point that was older than the universe..')
                #     break
                merger_ax_redshifts.append(t)
                # merger_ax_rates_total.append(np.sum(yield_locations[k]['yield'][()]))
                # we want to bin by metallicity
                # print(f'lookback: {t}, sum: {merger_ax_rates[r]}')
                # maybe we have to reweight the counts
                # yield_in_bin = yield_locations[k]['yield'][()]
                # counts, _ = np.histogram(yield_in_bin, weights=dco_mixture_weight)
                # merger_ax_rates[r] = np.sum(counts)
                # merger_ax[r].hist(
                #     yield_in_bin,
                #     weights=dco_mixture_weight,
                #     density=True
                # )
            # at this point we have the lookback times available that have also been filtered for unphysical time
            # now we can sort the times
            merger_ax_redshifts = np.sort(np.array(merger_ax_redshifts).astype(float))
            merger_rate_number_bins = 5 # keep this consistent for both metallicity and lookback
            # _, merger_metallicity_bins = np.histogram(dco_metallicities, bins=merger_rate_number_bins)
            merger_metallicity_bins = np.linspace(-4, -1.5, merger_rate_number_bins)
            merger_rates_binned_by_metallicities = np.zeros((len(merger_metallicity_bins), len(merger_ax_redshifts)))
            # now the times are sorted, go through each lookback time, digitize on metallicity
            m_bin_indices = np.digitize(np.log10(dco_metallicities), merger_metallicity_bins)

            # now we'd also like to bin on delay times well
            # _, merger_delay_bins = np.histogram(delayTimes, bins=merger_rate_number_bins)
            # merger_rates_binned_by_delaytimes = np.zeros((len(merger_delay_bins), len(merger_ax_lookbacks)))
            # dt_bin_indices = np.digitize(delayTimes, merger_delay_bins)

            # we also want to get the rates_per_system so we can re-weight later on
            all_weights = np.zeros((len(dco_metallicities), len(merger_ax_redshifts)))

            for j, lo in enumerate(merger_ax_redshifts):
                yields = yield_locations[u.Quantity(lo).to_string()]['yield'][()]
                all_weights[:, j] = yields
                for bin_index in range(len(merger_metallicity_bins)):
                    # go through each bin, sum yields at the locations where the bin is equal
                    
                    bin_query = np.where(m_bin_indices == bin_index+1) 
                    merger_rates_binned_by_metallicities[bin_index][j] = np.sum(
                        yields[bin_query]
                        # weights=metallicity_distro_at_cur_redshift[
                        #     np.digitize(np.log10(dco_metallicities), options.log_metallicity_bin_centers)-1
                        # ][bin_query]
                    )
                
                    # we want to plug in the metallicities to get the probability
                    # merger_rates_binned_by_delaytimes[bin_index][j] = np.sum(yields[np.where(dt_bin_indices == bin_index+1)])

            # merger_ax_rates_total = np.array(merger_ax_rates_total).astype(float)
            # assert 0
            # merger_ax_redshifts = [lookback_time_to_redshift(l, Planck13) for l in merger_ax_redshifts]

            merger_metallicity_bins = np.round(merger_metallicity_bins,2)
            return {
                'merger_rates_binned_by_metallicities': merger_rates_binned_by_metallicities,
                'merger_metallicity_bins': merger_metallicity_bins,
                'redshifts': np.round(merger_ax_redshifts,2),
                'y_units': y_axis_units,
                'all_weights': all_weights
            }
    
# loading the snr weights is computationally expensive
def load_snr_data(snr_file):
    snrs = {}
    with h5py.File(snr_file, 'r') as f:
        for key in f.keys():
            snrs[key] = np.array(f[key][:]).astype(np.float32)
    return snrs


# as a first pass, let's also do an SNR cut at each redshift
# so first get the redshifts (from the files themselves)
SNRsdir = '/Volumes/Elements/results_SNR'
# source_types = ['BBH', 'BNS', 'NSBH'] note that NSBH is backwards from BHNS, so will need to convert
Cat_base_name = 'SNR_Boesky_fullcat_'
z_bins = np.load('/Volumes/Elements/old_snr/data/redshift_bins_of_interest.npy')

def get_all_snr_weights(dets, snr_cut, dco_metallicities, types):
    # detector = 'KAGRA'
    # snr_cut = 1
    # for this one, we can rescale the count by finding how many objects make the cut and dividing by the total number of systems (not weights!)
    # snr_weights = np.zeros_like(fiducial_plotting_data['redshifts'])
    # this is probably a bad idea for RAM, but let's try a giant array that is per-system
    all_snr_weights = np.zeros((len(z_bins), len(dco_metallicities)))

    # have to iterate
    for i, z in enumerate(z_bins):
        zuse = z
        snrs_at_z = load_snr_data(os.path.join(SNRsdir, Cat_base_name + ('NSBH' if types == 'BHNS' else types) + f'_z{zuse:.2f}' + '_allDetectors.h5'))
        # here is where we need to gather each detector
        for d in dets:
            all_snr_weights[i, :] = np.power(snrs_at_z[d], 2) # sum quadrature
        all_snr_weights[i, :] = np.sqrt(all_snr_weights[i, :])
        # surviving = snrs_at_z[detector][np.where(snrs_at_z[detector] > snr_cut)[0]]
        # total = len(snrs_at_z[detector])

        # snrs_at_z[detector][snrs_at_z[detector] < snr_cut] = 0
        # snr_weights[i] = len(surviving)/total

        # print(len(surviving)/total)

        # store all weights
        # all_snr_weights[i, :] = snrs_at_z[detector]
        # print(snr_weights[i])
    # apply snr cut to all elements
    cut_indices = np.where(all_snr_weights > snr_cut) # do one overall cut
    all_snr_weights[:, :] = 0 # then we just want to count the ones that made the cut
    all_snr_weights[cut_indices] = 1 # zero out the systems that don't have that
    # print(all_snr_weights)
    return all_snr_weights
# since we're going to do only one system at a time, this is a list of the detectors in the network. will eventually make this into a class of some sort or maybe run configuration
networks = {
    # 'ET': {
    #     'array': [
    #         'ETS_0',
    #         'ETS_1',
    #         'ETS_2'
    #     ],
    #     'snr': 1
    # },
    'ETS': {
        'array': ['ETS'],
        'snr': 1
    },
    'ET_CE': {
        'array': ['ETS', 'CE1Id'],
        'snr': 8
    },
    'O4': {
        'array': ['H1_O4a', 'L1_O4a', 'Virgo_O4a'],
        'snr': 8
    },
    'O5': {
        'array': [
            'H1_postO5',
            'KAGRA_O5',
            'L1_O5',
            'Virgo_O5'
        ],
        'snr': 8
    }
}