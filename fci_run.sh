read -p "DCO_type: " dcotype
read -p "Z-step: " zstep

python ~/Documents/Github/COMPAS/compas_python_utils/cosmic_integration/FastCosmicIntegration.py \
--path /Volumes/Elements/Boesky_alpha0.1beta0.5.h5 --mu0 0.025 --muz -0.049 --sigma0 1.125 --sigmaz 0.048 --alpha -1.778 \
--zstep $zstep --sens "design" --aSF 0.02 --bSF 1.48 --cSF 4.44\
 --dSF 5.90 --maxz 10 --maxzdet 10 --dco_type $dcotype --weight mixture_weight

# python ~/Documents/Github/COMPAS/compas_python_utils/cosmic_integration/FastCosmicIntegration.py \
# --path /Volumes/Elements/Boesky_alpha0.1beta0.5.h5 --mu0 0.025 --muz -0.052 --sigma0 1.15 --sigmaz 0.0477 --alpha -1.88 \
# --zstep $zstep --sens "design" --aSF 0.017 --bSF 1.44 --cSF 4.53\
#  --dSF 6.23 --maxz 10 --maxzdet 10 --dco_type $dcotype --weight mixture_weight