# some properties we want to plot for the
# whole distribution -- hoping to explain the bump in the metallicity distribution
import h5py as h5
import numpy as np
import matplotlib.pyplot as plt
from compas_python_utils.cosmic_integration.ClassCOMPAS import COMPASData

compasdata = COMPASData(
    path='./Boesky_sims.h5',
    Mlower=0,
    Mupper=200,
    m2_min=0,
    binaryFraction=1#will need to change
    )
compasdata.setCOMPASDCOmask(types='all', withinHubbleTime=True)
compasdata.setCOMPASData()

delayTimes = compasdata.delayTimes / 1000
# weights = compasdata.weight
compasdata.setGridAndMassEvolved()
print(compasdata.totalMassEvolvedPerZ)

fig, ax = plt.subplots(1, 1)
ax.hist(delayTimes)
ax.set_xlabel('Delay time [Gyr]')
ax.set_ylabel('Weighted rate in COMPAS')
fig.savefig('./delaytimes.png')