import h5py as h5
file = h5.File('Boesky_sims.h5', 'r+')
rate_key = 'Rates_mu00.035_muz-0.23_alpha0.0_sigma00.39_sigmaz0.0'

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