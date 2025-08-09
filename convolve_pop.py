import os, copy, h5py
import astropy.units as u
from astropy.cosmology import Planck13
import numpy as np
import pandas as pd
from syntheticstellarpopconvolve import convolve, default_convolution_config, default_convolution_instruction
from syntheticstellarpopconvolve.general_functions import generate_boilerplate_outputfile, extract_unit_dict, temp_dir
from compas_python_utils.cosmic_integration.ClassCOMPAS import COMPASData
from syntheticstellarpopconvolve.starformation_rate_distributions import starformation_rate_distribution_vanSon2023
from syntheticstellarpopconvolve.general_functions import calculate_bincenters, calculate_bin_edges

# borrowed heavily from the example scipt from sspc

if __name__ == '__main__':

    # Create instance of output
    output_hdf5_filename = './output_example.h5' #os.path.join(TMP_DIR, "output_example.h5")
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

    delayTimes = compasdata.delayTimes / 1000

    # need to get the metallicities as well
    # print(np.sum(compasdata.DCOmask))
    # print(np.sum(all_dco_seeds))
    dco_metallicities = metallicities[np.isin(all_seeds, all_dco_seeds[compasdata.DCOmask])]
    assert len(delayTimes[()]) == len(dco_metallicities), "Something went wrong with masking for dco metallicities"
    dummy_data = {
        'delay_time': delayTimes[()],
        'probability': np.ones_like(delayTimes[()])
        # 'metallicity': dco_metallicities
    }
    dummy_df = pd.DataFrame.from_records(dummy_data) # load as dataframe
    dummy_df.to_hdf(output_hdf5_filename, key="input_data/example") # port pandas to hdf

    # Set up global configuration
    convolution_config = copy.copy(default_convolution_config)
    convolution_config["output_filename"] = output_hdf5_filename

    # Set up SFR -- to be determined later
    # convolution_config["SFR_info"] = {
    #     "lookback_time_bin_edges": np.array([0, 1, 2, 3, 4, 5]) * u.yr,
    #     "starformation_rate_array": np.array([1, 2, 3, 4, 5]) * u.Msun / u.yr
    # }
    num_redshift_bins = 100
    redshift_bin_edges = np.linspace(0, 10, num_redshift_bins)
    redshift_bin_centers = calculate_bincenters(redshift_bin_edges)
    # now convert redshift to lookback time
    #TODO: assume PLANK13 fit for now
    convolution_config['SFR_info'] = {
        'lookback_time_bin_edges': Planck13.lookback_time(redshift_bin_edges),
        'starformation_rate_array': starformation_rate_distribution_vanSon2023(redshift_bin_centers).to(u.Msun/u.yr/u.Gpc**3)
    }

    # set up convolution bin edges
    convolution_config["convolution_lookback_time_bin_edges"] = (
        np.array([0, 1, 2]) * u.yr # more or less sets the output bin locations
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
            }
        }
    ]
    # convolution_config['time_type'] = 'redshift'

    # run convolution
    convolve(convolution_config)

    # read out results
    with h5py.File(
        convolution_config["output_filename"], "r"
    ) as output_hdf5file:
        groupname = "output_data/example/conv_output/convolution_results/0.5 yr/"

        yield_data = output_hdf5file[groupname + "/yield"][()]
        print(len(yield_data[()]))
        unit_dict = extract_unit_dict(output_hdf5file, groupname)

        print(yield_data) # values
        print(unit_dict) # units
