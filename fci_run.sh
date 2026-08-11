read -p "DCO_type: " dcotype
# read -p "Z-step: " zstep
zstep=0.1
read -p "SFRD alias: " alias
fidstring="fiducial"
minstring="min"
maxstring="max"

if [ "$alias" == "$fidstring" ]; then
    echo "fiducial config selected"
    mu0=0.025
    muz=-0.049
    sigma0=1.125
    sigmaz=0.048
    alpha=1.778
    aSF=0.02
    bSF=1.48
    cSF=4.44
    dSF=5.90
elif [ "$alias" == "$minstring" ]; then
    echo "min config selected"
    mu0=0.007
    muz=0
    sigma0=0.7
    sigmaz=0
    alpha=-6
    aSF=0.01
    bSF=2.60
    cSF=3.20
    dSF=6.20
elif [ "$alias" == "$maxstring" ]; then
    echo "max config selected"
    mu0=0.035
    muz=-0.5
    sigma0=2.0
    sigmaz=0.1
    alpha=0.0
    aSF=0.03
    bSF=2.6
    cSF=3.3
    dSF=5.9
else
    echo "incorrect alias"
    exit 2
fi

# echo $mu0

python ~/Documents/Github/COMPAS/compas_python_utils/cosmic_integration/FastCosmicIntegration.py \
--path /Volumes/Elements/Boesky_alpha0.1beta0.5.h5 --mu0 $mu0 --muz $muz --sigma0 $sigma0 --sigmaz $sigmaz --alpha $alpha \
--zstep $zstep --sens "design" --aSF $aSF --bSF $bSF --cSF $cSF\
 --dSF $dSF --maxz 10 --maxzdet 10 --dco_type $dcotype --weight mixture_weight --sfrd_alias $alias


# python ~/Documents/Github/COMPAS/compas_python_utils/cosmic_integration/FastCosmicIntegration.py \
# --path /Volumes/Elements/Boesky_alpha0.1beta0.5.h5 --mu0 0.025 --muz -0.049 --sigma0 1.125 --sigmaz 0.048 --alpha -1.778 \
# --zstep $zstep --sens "design" --aSF 0.02 --bSF 1.48 --cSF 4.44\
#  --dSF 5.90 --maxz 10 --maxzdet 10 --dco_type $dcotype --weight mixture_weight --sfrd_alias $alias

# python ~/Documents/Github/COMPAS/compas_python_utils/cosmic_integration/FastCosmicIntegration.py \
# --path /Volumes/Elements/Boesky_alpha0.1beta0.5.h5 --mu0 0.025 --muz -0.052 --sigma0 1.15 --sigmaz 0.0477 --alpha -1.88 \
# --zstep $zstep --sens "design" --aSF 0.017 --bSF 1.44 --cSF 4.53\
#  --dSF 6.23 --maxz 10 --maxzdet 10 --dco_type $dcotype --weight mixture_weight