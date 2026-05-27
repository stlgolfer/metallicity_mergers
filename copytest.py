import h5py as h5
import numpy as np

def create_dummy_file():
    #create a dummy file with some data
    with h5.File('./dummy.h5', 'w') as fw:
        # print(f.keys())
        fw.create_group('Rates')
        fw['Rates'].create_dataset('merger_rate', data=np.random.rand(20, 100))
        fw['Rates'].create_dataset('redshifts', data=np.linspace(0, 10, num=101))
        
        # print(fw['test'][()])
        print(fw['Rates'].keys())

def combine_rates(rates, red_lefts, r_parity, jump, downsample=True):
    # rates = np.random.randn(20, 100) # simulating 20 systems with 100 redshifts
    # red_lefts = np.linspace(0, 10, num=101) # needs to have 100+1 to simulate left bins
    print(rates.shape)
    print(red_lefts)
    # now, we need to apply a parity and combine
    # r_parity = 2
    z_loc = np.digitize(r_parity, red_lefts) # removing the -1 so we start collecting
    # on the next bin. avoids double counting

    # hm, maybe get the indices manually
    new_bin_indices = []
    j = z_loc #+ jump # this is the step that adjusts the boundary
    # jump = 5 # number of bins
    while j + jump <= len(red_lefts):
        new_bin_indices.append(j)
        j += jump
    # manually add the last bin
    new_bin_indices.append(len(red_lefts) - 1)
    print(red_lefts[new_bin_indices])

    # now we want to actually combine
    new_red_lefts = np.concatenate((red_lefts[:z_loc], red_lefts[new_bin_indices]))
    print(new_red_lefts.shape)
    #hm, maybe we can grab a region, push to 3d by reshape, and then sum?
    # print(rates[np.newaxis, :, :].reshape(5,20,20).shape) eh, just iteratively for now
    new_rates = np.zeros((len(red_lefts[new_bin_indices])-1, rates.shape[0]))
    for r in range(len(new_bin_indices) - 1):
        start = new_bin_indices[r] + np.floor(jump/2).astype(int)

        if downsample:
            new_rates[r,:] = rates[:, start]
        else:
            end = new_bin_indices[r+1]
            block_sum = np.sum(rates[:,start:end],axis=1)
            new_rates[r,:] = block_sum/jump
        # print(block_sum.shape)
        # assert 0
        #block_sum/jump
        # ah, don't do block sum, just downsample the intrinsic rate.
        # this preserves units and still lowers the number of bins
    full_new_rates = np.concatenate((rates[:, :z_loc], new_rates.T), axis=1)
    # full_new_red_lefts = np.concatenate((red_lefts[:z_loc], new_red_lefts))
    print(full_new_rates.shape)
    print(new_red_lefts.shape)
    return new_red_lefts, full_new_rates

# def copy_group_all(f, source, dest, blacklist=[]):

def clean_keys(f, blacklist):
    for b in blacklist:
        if b in f.keys():
            del f[b]
    # del fr[new_key]['merger_rate']
    #     del fr[new_key]['redshifts']
    #     if 'detection_ratedesign' in fr[new_key].keys():
    #         del fr[new_key]

if __name__ == '__main__':
    # create_dummy_file()
    fname = '/Volumes/Elements/Boesky_alpha0.1beta0.5.h5'
    with h5.File(fname, 'r+') as fr:
        species = 'BHNS'
        old_key = f'Rates_mu00.025_muz-0.052_alpha-1.88_sigma01.15_sigmaz0.0477_{species}_0.1_10.0'
        new_key = f'Rates_{species}_mixed'

        assert old_key in fr.keys(), f'Try these keys: {fr.keys()}'

        if new_key in fr.keys():
            del fr[new_key]

        fr.copy(fr[old_key], fr, new_key)
        # then delete these keys
        clean_keys(fr[new_key], [
            'merger_rate',
            'redshifts',
            'detection_ratedesign',
            'merger_rate_z0'
            ])
        
        print(fr[new_key].keys())
        # fr.create_group(new_key) # but now key is already amde
        new_reds, new_rates = combine_rates(
            fr[old_key]['merger_rate'][()],
            fr[old_key]['redshifts'][()],
            2, 5, True
        )
        fr[new_key].create_dataset('merger_rate', data=new_rates)
        fr[new_key].create_dataset('redshifts', data=new_reds)
            
    # print('Now reading')
    # with h5.File('./dummy.h5', 'r') as fr:
    #     print(fr.keys())
    # print('Now creating new group')
    # with h5.File('./dummy.h5', 'r+') as fw: # ok so we need the r+ qualifier
    #     if 'Rates' not in fw.keys():
    #         fw.create_group('Rates')
    #     if 'yuh' not in fw['Rates'].keys():
    #         fw['Rates'].create_dataset('yuh', data=np.random.randn(10))
    #     print(fw['Rates'])

    # now make a dummy dataset
    
    
    # now grab the large
    # blocked = rates[:,z_loc:]
    # course_red_lefts = red_lefts[z_loc:]
    # splitted = np.array_split(course_red_lefts, 5)
    # for i in range(len(splitted)):
    #     splitted[i] = np.array([splitted[i][0],splitted[i][-1]])
    # print(splitted) # this is the list of intervals for the new bins