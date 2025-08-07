import h5py as h5
import sys

file = h5.File('/Volumes/Elements/Boesky_sims.h5', 'r+')
all_keys = list(file.keys())
for i, k in enumerate(all_keys):
    print(f'[{i}]: {k}')
rate_key = all_keys[int(input('Select rate key: '))]
# rate_key = 'Rates_mu00.035_muz-0.23_alpha0.0_sigma00.39_sigmaz0.0'

if input('Delete this key? [y/n]') == 'y':
    del file[rate_key]
    print('Done')
    sys.exit(0)
# print("rate" in file[rate_key])
print('The following will be deleted (y/n):')
for key in file[rate_key].keys():
    if "detection" in key:
        # del key
        print(key)
if input('Confirm [y/n]: ') == 'y':
    print('delete selected')
    for key in file[rate_key].keys():
        if "detection" in key:
            del file[rate_key][key]
            # print(key)
print(file[rate_key].keys())
file.close()