import h5py as h5
file = h5.File('Boesky_sims.h5') 
print(file.keys())
print(file['Rates_mu00.035_muz0.035_alpha-1.778_sigma01.122_sigmaz0.049'].keys())
# print(file['Rates_mu00.035_muz-0.23_alpha0.0_sigma00.39_sigmaz0.0']['detection_rateVoyager.txt'])
# print(file['Rates_mu00.035_muz-0.23_alpha0.0_sigma00.39_sigmaz0.0']['detection_rateCE.txt'])