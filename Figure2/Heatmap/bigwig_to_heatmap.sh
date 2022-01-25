# Author: Sora Yoon
# VahediLab
# Date: 2020. 12. 31
# Goal: (1) Deeptools analysis with B cell data.
#     	(2) A/B compartment analysis
#       (3) Insulation score analysis
#       (4) Stripy/ Non-stripy TAD comparison

# Stripe filter: P<0.1 && Stripiness > 0

hic=/mnt/data0/sora/stripe/hic/Vian_aB_30hrs.mcool::resolutions/5000
stripe_out=/mnt/data0/sora/stripe/stripenn_out/Vian_aB_30hrs_5kb_VCSQRT3/

#stripenn compute --cool $hic --out $stripe_out --norm VC_SQRT -w 16 -m 0.97,0.975,0.98,0.985,0.99 -n 8 -p 0.1

mainDir=/mnt/data0/sora/stripe/
workDir=/mnt/data0/sora/stripe/ChIP_hiC/activated_B/

r=5
p=0.05
s=0

declare -a test_proteins=( CTCF RAD21 NIPBL SMC3 PolII )

stripe=${stripe_out}result_unfiltered.txt
stripe_f=${stripe_out}/result_filtered_p${p}_S${s}.txt

awk -F "\t" '{ if(($11 < "'"$p"'") && ($12 >= "'"$s"'") && NR > 1) {print}}' $stripe > $stripe_f

# Deeptools analysis
mainDir=/mnt/data0/sora/stripe/
workDir=/mnt/data0/sora/stripe/ChIP_hiC/activated_B/
DTdir=${mainDir}'deeptools/2021_01_04_'${r}'kb_p_'${p}'_s'${s}'/'
mkdir $DTdir

for protein in "${test_proteins[@]}"
do
  if [ $protein = CTCF ]
  then
    bw=${workDir}'CTCF/05.bigwig/merged.raw.bw'
    tp=30
  elif [ $protein = RAD21 ]
  then
    bw=${workDir}'Rad21/05.bigwig/merged.raw.bw'
    tp=30
  elif [ $protein = NIPBL ]
  then
    bw=${workDir}'Nipbl/05.bigwig/merged.raw.bw'
    tp=30
  elif [ $protein = SMC3 ]
  then
    bw=${workDir}'Smc3/05.bigwig/merged.raw.bw'
    tp=30
  elif [ $protein = PolII ]
  then
    bw=${workDir}'PolII/05.bigwig/4DNFI7MDXWQU.raw.bw'
    tp=30
  else
    echo "Wrong input"
  fi

  stripe=$stripe_f

  stripe_bed_up=${stripe/.txt/_up.bed}
  stripe_bed_down=${stripe/.txt/_down.bed}

  awk '$2==$5 {print $1"\t"$5"\t"$6}' $stripe > $stripe_bed_down
  awk '$3==$6 {print $1"\t"$5"\t"$6}' $stripe > $stripe_bed_up

  sudo Rscript bigwig_to_heatmap.R $bw $stripe_bed_down ${DTdir}${protein}_down.pdf
  sudo Rscript bigwig_to_heatmap.R $bw $stripe_bed_up ${DTdir}${protein}_up.pdf

done