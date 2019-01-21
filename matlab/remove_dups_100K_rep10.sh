export FILES="EURqc_chr8_chunk125_plink_100K EURqc_chr8_chunk5_plink_100K EURqc_chr15_chunk86_plink_100K EURqc_chr9_chunk68_plink_100K EURqc_chr8_chunk217_plink_100K EURqc_chr2_chunk136_plink_100K EURqc_chr2_chunk15_plink_100K EURqc_chr2_chunk14_plink_100K EURqc_chr1_chunk20_plink_100K EURqc_chr8_chunk126_plink_100K EURqc_chr9_chunk69_plink_100K EURqc_chr8_chunk6_plink_100K EURqc_chr9_chunk70_plink_100K"
export SNPS="rs112741103,rs1677780,rs1671787,rs71501039,rs141927528,rs142701510,rs145926341,rs542232278,rs1610794,rs371891811,rs62011115,rs28763720,rs115010285,rs59201531,rs145970356,rs370528220,rs146021604,rs71218065,rs2558060,rs111827956,rs13374233,rs6658405,rs202107783,rs200269882,rs201680403,rs200141404,rs201551115,rs11261873,rs201093356"
module load plink
for FILE in $FILES; do
    for REPI in `seq 1 10`; do
        #mv /work/users/oleksanf/HAPGEN/"$FILE"_rep"$REPI".bim /work/users/oleksanf/HAPGEN/"$FILE"_rep"$REPI".dups.bim
        #mv /work/users/oleksanf/HAPGEN/"$FILE"_rep"$REPI".bed /work/users/oleksanf/HAPGEN/"$FILE"_rep"$REPI".dups.bed
        #mv /work/users/oleksanf/HAPGEN/"$FILE"_rep"$REPI".fam /work/users/oleksanf/HAPGEN/"$FILE"_rep"$REPI".dups.fam
        #mv /work/users/oleksanf/HAPGEN/"$FILE"_rep"$REPI".nosex /work/users/oleksanf/HAPGEN/"$FILE"_rep"$REPI".dups.nosex
        #mv /work/users/oleksanf/HAPGEN/"$FILE"_rep"$REPI".log /work/users/oleksanf/HAPGEN/"$FILE"_rep"$REPI".dups.log
        plink --bfile /work/users/oleksanf/HAPGEN/"$FILE"_rep"$REPI".dups --make-bed --out /work/users/oleksanf/HAPGEN/"$FILE"_rep"$REPI" --exclude-snps $SNPS
    done
done

