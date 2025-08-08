# some properties we want to plot for the
# whole distribution -- hoping to explain the bump in the metallicity distribution
import h5py as h5
import numpy as np
import matplotlib.pyplot as plt
from compas_python_utils.cosmic_integration.ClassCOMPAS import COMPASData
from tqdm import tqdm

if __name__ == '__main__':
    # get binary fraction
    fdata = h5.File('/Volumes/Elements/Boesky_sims.h5')
    dco_seeds = fdata['BSE_Double_Compact_Objects']['SEED'][()]
    no_dcos = len(dco_seeds)
    all_seeds = fdata['BSE_System_Parameters']['SEED'][()]
    no_total = len(all_seeds)
    frac = no_dcos / no_total

    m1s = fdata['BSE_System_Parameters']['Mass@ZAMS(1)'][()]
    m2s = fdata['BSE_System_Parameters']['Mass@ZAMS(2)'][()]
    metallicities = fdata["BSE_System_Parameters"]["Metallicity@ZAMS(1)"][()]
    fdata.close()

    compasdata = COMPASData(
        path='/Volumes/Elements/Boesky_sims.h5',
        Mlower=0,
        Mupper=200,
        m2_min=0,
        binaryFraction = frac
        )
    compasdata.setCOMPASDCOmask(types='all', withinHubbleTime=True)
    compasdata.setCOMPASData()

    delayTimes = compasdata.delayTimes / 1000
    # weights = compasdata.weight

    fig, ax = plt.subplots(1, 1)
    ax.hist(delayTimes)
    ax.set_xlabel('Delay time [Gyr]')
    ax.set_ylabel('Weighted rate in COMPAS')
    fig.savefig('./delaytimes.png')

    # now we'd like to plot the formation efficiency
    # there is a function to do this in compas, though it will literally iterate through each
    # unique metallicity, which we might not need
    # compasdata.setGridAndMassEvolved() # this is too slow, even with parallel
    # print(compasdata.totalMassEvolvedPerZ)

    # fraction is 1.77
    # now get the metallicity grad
    zmin = np.min(metallicities)
    zmax = np.max(metallicities)
    metallicity_grid = np.random.choice(metallicities, 1000) # here we have to uniformly choose
    # can't just use linspace since not all of those metallicities are guaranteed to exist
    # for each metallicity, count the total mass
    total = np.zeros(len(metallicity_grid))
    Ndcos_per_metallicity = np.zeros(len(metallicity_grid))
    dco_locs = np.isin(all_seeds, dco_seeds)
    with tqdm(len(metallicity_grid)) as pbar:
        for i, Z in enumerate(metallicity_grid):
                mask = metallicities == Z
                total[i] = np.sum(m1s[mask]) + np.sum(m2s[mask])

                # at the same time, get the number of dcos
                # logical AND the dco mask and the metallicity locations and sum
                # ah no dco mask is for the dco key, so need to match seeds
                
                Ndcos_per_metallicity[i] = np.sum(dco_locs & mask) # total number of dcos with metallicity fixed
                pbar.update(1)
    total / 1.77
    # for each metallicity, get the total number of dcos in that
    f_eff = np.divide(Ndcos_per_metallicity, total)
    eff_fig, eff_ax = plt.subplots(1, 1)
    eff_ax.plot(np.sort(metallicity_grid), f_eff)
    # eff_ax.hist(f_eff, bins = np.sort(metallicity_grid))
    # ax.set_xlabel('Delay time [Gyr]')
    # ax.set_ylabel('Weighted rate in COMPAS')
    eff_fig.savefig('./formation_efficiency.png')
    # instead, maybe use the same code but only for a few metallicities
    