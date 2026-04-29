import h5py
from syntheticstellarpopconvolve.general_functions import generate_boilerplate_outputfile, extract_unit_dict, temp_dir
import astropy.units as u
import numpy as np

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