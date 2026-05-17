import h5py as h5
import numpy as np
import matplotlib.pyplot as plt
from warnings import warn

if __name__ == '__main__':
    fname = '/Volumes/Elements/Boesky_alpha0.1beta0.5.h5'
    f = h5.File(fname)
    pop_type = 'BBH'

    r1 = str(0.1)
    r2 = str(0.5)
    rates_prefix = f'Rates_mu00.025_muz-0.052_alpha-1.88_sigma01.15_sigmaz0.0477_{pop_type}'
    rates_postfix = '10.0'

    r1_key = f'{rates_prefix}_{r1}_{rates_postfix}'
    r2_key = f'{rates_prefix}_{r2}_{rates_postfix}'
    print(f[r1_key]['merger_rate'][()].shape)
    print(f[r2_key]['merger_rate'][()].shape)
    # now we want to make a new database in the h5 file that has one big table
    # we need to combine the redshifts
    # <KeysViewHDF5 ['DCOmask', 'SEED', 'detection_ratedesign', 'merger_rate', 'merger_rate_z0', 'redshifts']>
    # we also need to copy the other tables over
    redshift_parity = 2

    r1_redshifts = f[r1_key]['redshifts'][()]
    r2_redshifts = f[r2_key]['redshifts'][()]

    r1_rcutoff_index = np.digitize(redshift_parity, r1_redshifts)-1
    r2_rcutoff_index = np.digitize(redshift_parity, r2_redshifts)-1

    # hm, need to deal with the redshift bins might double count certain contributions
    warn("Haven't figured out how to deal with bin overlap")

    new_redshifts = np.concatenate(
        (r1_redshifts[:r1_rcutoff_index], r2_redshifts[r2_rcutoff_index+1:])
    )
    print(new_redshifts)
