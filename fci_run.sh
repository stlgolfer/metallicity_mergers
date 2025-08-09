read -p "DCO_type: " dcotype
read -p "Z-step: " zstep
read -p "mu0 (default 0.035): " mu0
read -p "muz: " muz

python ~/Documents/Github/COMPAS/compas_python_utils/cosmic_integration/FastCosmicIntegration.py \
--path /Volumes/Elements/Boesky_sims.h5 --mu0 $mu0 --muz $muz --sigma0 1.122 --sigmaz 0.049 --alpha -1.778 \
--zstep $zstep --sens "CE.txt" --m1min 10. --aSF 0.02 --bSF 1.48 --cSF 4.44\
 --dSF 5.9 --maxz 10 --maxzdet 10 --dco_type $dcotype

 python ~/Documents/Github/COMPAS/compas_python_utils/cosmic_integration/FastCosmicIntegration.py \
--path /Volumes/Elements/Boesky_sims.h5 --mu0 $mu0 --muz $muz --sigma0 1.122 --sigmaz 0.049 --alpha -1.778 \
 --zstep $zstep --sens "design" --m1min 10. --aSF 0.02 --bSF 1.48 --cSF 4.44\
 --dSF 5.9 --maxz 10 --maxzdet 10 --dco_type $dcotype
 python ./print_keys.py