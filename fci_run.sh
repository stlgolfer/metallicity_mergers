read -p "DCO_type: " dcotype
read -p "Z-step: " zstep

python ~/Documents/Github/COMPAS/compas_python_utils/cosmic_integration/FastCosmicIntegration.py \
--path ./Boesky_sims.h5 --mu0 0.035 --muz -0.23 --sigma0 0.39 --sigmaz 0.0 --alpha 0.0 \
--weight mixture_weight --zstep $zstep --sens "O3" --m1min 10. --aSF 0.01 --bSF 2.77 --cSF 2.9\
 --dSF 4.7 --maxz 10 --maxzdet 10 --dco_type $dcotype

 python ~/Documents/Github/COMPAS/compas_python_utils/cosmic_integration/FastCosmicIntegration.py \
--path ./Boesky_sims.h5 --mu0 0.035 --muz -0.23 --sigma0 0.39 --sigmaz 0.0 --alpha 0.0 \
--weight mixture_weight --zstep $zstep --sens "CE.txt" --m1min 10. --aSF 0.01 --bSF 2.77 --cSF 2.9\
 --dSF 4.7 --maxz 10 --maxzdet 10 --dco_type $dcotype
 python ./print_keys.py