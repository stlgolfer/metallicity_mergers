# some properties we want to plot for the
# whole distribution -- hoping to explain the bump in the metallicity distribution
import h5py as h5
import numpy as np
import matplotlib.pyplot as plt
from compas_python_utils.cosmic_integration.ClassCOMPAS import COMPASData
from tqdm import tqdm
from scipy import stats

# From floor
def analytical_star_forming_mass_per_binary_using_kroupa_imf(
        m1_min, m1_max, m2_min, fbin=1., imf_mass_bounds=[0.01,0.08,0.5,200]
):
    """
    Analytical computation of the mass of stars formed per binary star formed within the
    [m1 min, m1 max] and [m2 min, ..] rage,
    using the Kroupa IMF:

        p(M) \propto M^-0.3 for M between m1 and m2
        p(M) \propto M^-1.3 for M between m2 and m3;
        p(M) = alpha * M^-2.3 for M between m3 and m4;

    @Ilya Mandel's derivation
    """
    m1, m2, m3, m4 = imf_mass_bounds
    if m1_min < m3:
        raise ValueError(f"This analytical derivation requires IMF break m3  < m1_min ({m3} !< {m1_min})")
    alpha = (-(m4**(-1.3)-m3**(-1.3))/1.3 - (m3**(-0.3)-m2**(-0.3))/(m3*0.3) + (m2**0.7-m1**0.7)/(m2*m3*0.7))**(-1)
    # average mass of stars (average mass of all binaries is a factor of 1.5 larger)
    m_avg = alpha * (-(m4**(-0.3)-m3**(-0.3))/0.3 + (m3**0.7-m2**0.7)/(m3*0.7) + (m2**1.7-m1**1.7)/(m2*m3*1.7))
    # fraction of binaries that COMPAS simulates
    fint = -alpha / 1.3 * (m1_max ** (-1.3) - m1_min ** (-1.3)) + alpha * m2_min / 2.3 * (m1_max ** (-2.3) - m1_min ** (-2.3))
    # mass represented by each binary simulated by COMPAS
    m_rep = (1/fint) * m_avg * (1.5 + (1-fbin)/fbin)
    return m_rep

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
    mixture_weights_system_params = fdata['BSE_System_Parameters']['mixture_weight'][()]

    #region get information for kroupa imf -- all from floor
    initial_mass_min = fdata['Run_Details']['initial-mass-min'][()][0]
    initial_mass_max = fdata['Run_Details']['initial-mass-max'][()][0] 
    minimum_secondary_mass = fdata['Run_Details']['minimum-secondary-mass'][()][0] 
    f_binary = 1


    m_rep_per_binary = analytical_star_forming_mass_per_binary_using_kroupa_imf(m1_min=initial_mass_min, m1_max=initial_mass_max,\
                                                                            m2_min=minimum_secondary_mass, fbin=f_binary)

    print('1 binary in COMPAS represents', m_rep_per_binary, ' solar masses formed')
    # now calculate the number of binaries in COMPAS simulation (over the entire simulation)
    n_binaries = np.shape(fdata['BSE_System_Parameters']['SEED'][()])[0]
    print(n_binaries)


    total_mass_evolved_compas = n_binaries * m_rep_per_binary
    print(total_mass_evolved_compas, ' [Msun]')
    #endregion
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

    dco_locs = np.isin(all_seeds, dco_seeds) # TODO: need to mask with dcomask

    # going to try the weights method Floor had suggested
    # first get the metallicities of all the dcos
    # for each metallicity, get the total number of dcos in that
    # f_eff = np.divide(Ndcos_per_metallicity, total)
    eff_fig, eff_ax = plt.subplots(1, 1)
    metallicity_dcos = np.log10(metallicities[dco_locs]/0.012)
    total_bins = 50

    eff_ax.hist(
        metallicity_dcos/total_mass_evolved_compas,
        weights=mixture_weights_system_params[dco_locs],
        bins=total_bins,
        density=True,
        label='Density histogram'
    )

    _, bins = np.histogram(metallicity_dcos, bins=total_bins)
    metallicitykde = stats.gaussian_kde(
        metallicity_dcos.flatten(),
        weights=mixture_weights_system_params[dco_locs].flatten()
    )

    eff_ax.plot(
        bins[:-1], metallicitykde(bins[:-1])/total_mass_evolved_compas,
        label='KDE plot'
        # label=f'(1) Type {stellar_types_dictionary[type_index]} ({detector})'
    )
    eff_ax.fill_between(
        bins[:-1],
        metallicitykde(bins[:-1])/total_mass_evolved_compas,
        interpolate=True,
        alpha=0.3
    )

    eff_ax.set_yscale('log')
    eff_ax.set_xlabel('Log10(Z_dco / Zsun)')
    eff_ax.set_ylabel('R_form 1/M0')
    # eff_ax.set_ylim(min(_), 1)
    eff_ax.set_title('Formation eff. up to constant factor')
    eff_ax.legend()
    eff_fig.savefig('./formation_efficiency.png')
    