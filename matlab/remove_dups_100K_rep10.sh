export FILES="EURqc_chr8_chunk125_plink_100K EURqc_chr8_chunk5_plink_100K EURqc_chr15_chunk86_plink_100K EURqc_chr9_chunk68_plink_100K EURqc_chr8_chunk217_plink_100K EURqc_chr2_chunk136_plink_100K EURqc_chr2_chunk15_plink_100K EURqc_chr2_chunk14_plink_100K EURqc_chr1_chunk20_plink_100K EURqc_chr8_chunk126_plink_100K EURqc_chr9_chunk69_plink_100K EURqc_chr8_chunk6_plink_100K EURqc_chr9_chunk70_plink_100K"
export SNPS=exclude_snps.txt
module load plink
for FILE in $FILES; do
    for REPI in `seq 1 10`; do
        #mv /work/users/oleksanf/HAPGEN/"$FILE"_rep"$REPI".bim /work/users/oleksanf/HAPGEN/"$FILE"_rep"$REPI".dups.bim
        #mv /work/users/oleksanf/HAPGEN/"$FILE"_rep"$REPI".bed /work/users/oleksanf/HAPGEN/"$FILE"_rep"$REPI".dups.bed
        #mv /work/users/oleksanf/HAPGEN/"$FILE"_rep"$REPI".fam /work/users/oleksanf/HAPGEN/"$FILE"_rep"$REPI".dups.fam
        #mv /work/users/oleksanf/HAPGEN/"$FILE"_rep"$REPI".nosex /work/users/oleksanf/HAPGEN/"$FILE"_rep"$REPI".dups.nosex
        #mv /work/users/oleksanf/HAPGEN/"$FILE"_rep"$REPI".log /work/users/oleksanf/HAPGEN/"$FILE"_rep"$REPI".dups.log
        plink --bfile /work/users/oleksanf/HAPGEN/"$FILE"_rep"$REPI".dups --make-bed --out /work/users/oleksanf/HAPGEN/"$FILE"_rep"$REPI" --exclude $SNPS
    done
done

