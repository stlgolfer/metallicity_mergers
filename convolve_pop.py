import os, copy, h5py
import astropy.units as u
from astropy.cosmology import Planck13
import numpy as np
import pandas as pd
from syntheticstellarpopconvolve import convolve, default_convolution_config, default_convolution_instruction
from syntheticstellarpopconvolve.general_functions import generate_boilerplate_outputfile, extract_unit_dict, temp_dir
from compas_python_utils.cosmic_integration.ClassCOMPAS import COMPASData
from syntheticstellarpopconvolve.starformation_rate_distributions import starformation_rate_distribution_vanSon2023
from syntheticstellarpopconvolve.metallicity_distributions import metallicity_distribution_vanSon2022
from syntheticstellarpopconvolve.general_functions import calculate_bincenters, calculate_bin_edges
from syntheticstellarpopconvolve.SFR_dict_plotting_routines import plot_sfr_dict
import matplotlib.pyplot as plt

# borrowed heavily from the example scipt from sspc

if __name__ == '__main__':

    # Create instance of output
    output_hdf5_filename = '/Volumes/Elements/sspc_output.h5' #os.path.join(TMP_DIR, "output_example.h5")
    generate_boilerplate_outputfile(output_hdf5_filename)

    # SET UP DATA -- from a basic dict, port to an h5 file and then store this into the output file
    # now we want to use some real data
    fdata = h5py.File('/Volumes/Elements/Boesky_sims.h5')
    all_dco_seeds = fdata['BSE_Double_Compact_Objects']['SEED'][()]
    all_seeds = fdata['BSE_System_Parameters']['SEED'][()]
    metallicities = fdata["BSE_System_Parameters"]["Metallicity@ZAMS(1)"][()]
    fdata.close()

    compasdata = COMPASData(
        path='/Volumes/Elements/Boesky_sims.h5'
    )
    compasdata.setCOMPASDCOmask(types='all', withinHubbleTime=True)
    compasdata.setCOMPASData()

    delayTimes = compasdata.delayTimes #/ 1000

    # need to get the metallicities as well
    dco_metallicities = metallicities[np.isin(all_seeds, all_dco_seeds[compasdata.DCOmask])]
    assert len(delayTimes[()]) == len(dco_metallicities), "Something went wrong with masking for dco metallicities"
    dummy_data = {
        'delay_time': delayTimes[()]*u.Myr,
        'probability': (10e-4)*np.ones_like(delayTimes[()]),
        'metallicity': np.log10(dco_metallicities) # I wonder if we need to log these as well?
    }
    dummy_df = pd.DataFrame.from_records(dummy_data) # load as dataframe
    dummy_df.to_hdf(output_hdf5_filename, key="input_data/example") # port pandas to hdf

    # Set up global configuration
    convolution_config = copy.copy(default_convolution_config)
    convolution_config["output_filename"] = output_hdf5_filename

    # Set up SFR, and hopefully also setup metallicity dependence
    num_redshift_bins = 100
    redshift_bin_edges = np.linspace(0, 10, num_redshift_bins)
    redshift_bin_centers = calculate_bincenters(redshift_bin_edges)
    # calculate metallicity distros
    num_metallicities = 500
    log_metallicity_bin_edges = np.linspace(
        np.log10(np.min(dco_metallicities)),
        np.log10(np.max(dco_metallicities)),
        num_metallicities
    )
    
    log_metallicity_bin_centers = calculate_bincenters(log_metallicity_bin_edges)
    # print(np.diff(log_metallicity_bin_centers))
    # assert 0
    
    dpdlogZ = metallicity_distribution_vanSon2022(
        log_metallicity_centers=log_metallicity_bin_centers,
        redshifts=redshift_bin_centers,
    )

    sfr = starformation_rate_distribution_vanSon2023(redshift_bin_centers).to(u.Msun/u.yr/u.Gpc**3)

    high_res_sfr_dict = {
        "redshift_bin_edges": redshift_bin_edges,
        "starformation_rate_array": sfr,
        "metallicity_bin_edges": log_metallicity_bin_edges,
        "metallicity_distribution_array": dpdlogZ,  # We need to transpose!
    }

    axis_dict = plot_sfr_dict(
        high_res_sfr_dict,
        time_type="redshift",
        metallicity_string="logZ",
        metallicity_distribution_multiply_by_metallicity_bin_sizes=False,
        metallicity_distribution_multiply_by_sfr=False,
        metallicity_distribution_scale="log10",
        metallicity_distribution_cmap=copy.copy(plt.cm.viridis),
        return_axis_dict=True,
        figsize=(8,8),
        fontsize=12,
    )
    axis_dict['fig'].savefig('./sfr.png')

    #TODO: assume PLANK13 fit for now
    convolution_config['SFR_info'] = {
        'lookback_time_bin_edges': Planck13.lookback_time(redshift_bin_edges),
        'starformation_rate_array': sfr,
        'metallicity_bin_edges': log_metallicity_bin_edges,
        'metallicity_distribution_array': dpdlogZ
    }

    # set up convolution bin edges
    convolution_config["convolution_lookback_time_bin_edges"] = (
        np.linspace(0,10, 10) * u.Gyr # more or less sets the output bin locations
        # recall that lookback time is related to redshift, so this is basically
        # setting the x-axis for the output
    )

    # Set up the convolution instructions
    convolution_config["convolution_instructions"] = [
        {
            **default_convolution_instruction,
            "input_data_name": "example",
            "output_data_name": "conv_output",
            "data_column_dict": {
                "delay_time": "delay_time",
                "normalized_yield": {"column_name": "probability", "unit": 1/u.Msun},
                'metallicity': 'metallicity'
            }
        }
    ]

    # run convolution
    convolve(convolution_config)

    # read out results
    with h5py.File(
        convolution_config["output_filename"], "r"
    ) as output_hdf5file:
        groupname = "output_data/example/conv_output/convolution_results/"

        # yield_data = output_hdf5file[groupname + "/yield"][()]
        # print(len(yield_data[()]))
        # unit_dict = extract_unit_dict(output_hdf5file, groupname)

        # print(yield_data) # values
        # print(unit_dict) # units
        # could try and make an actual merger plot. go through each key, sum the values
        yield_locations = output_hdf5file['output_data/example/conv_output/convolution_results/']
        merger_fig, merger_ax = plt.subplots(1, 1)
        merger_ax_redshifts = np.zeros_like(list(yield_locations.keys()))
        merger_ax_rates = np.zeros_like(merger_ax_redshifts)
        for r, k in enumerate(yield_locations.keys()):
            # get number
            units = u.Quantity(k).unit
            time = u.Quantity(k).value
            merger_ax_redshifts[r] = time
            merger_ax_rates[r] = np.sum(yield_locations[k]['yield'][()])
        merger_ax.plot(merger_ax_redshifts, merger_ax_rates)
        merger_fig.savefig('merger rates.png')
