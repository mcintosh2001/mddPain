~/ldscore/ldsc/ldsc.py \
 --h2 pain_vonkorff_ldsc.sumstats.gz \
 --ref-ld-chr ~/ldscore/ldsc/refFiles/eur_w_ld_chr/ \
 --w-ld-chr ~/ldscore/ldsc/refFiles/eur_w_ld_chr/   \
 --out pain_vonkorff_h2

~/ldscore/ldsc/ldsc.py \
 --h2 mdd_ldsc.sumstats.gz \
 --ref-ld-chr ~/ldscore/ldsc/refFiles/eur_w_ld_chr/ \
 --w-ld-chr ~/ldscore/ldsc/refFiles/eur_w_ld_chr/   \
 --pop-prev 0.12 \
 --samp-prev 0.5 \
 --out mdd_h2

~/ldscore/ldsc/ldsc.py \
 --rg mdd_ldsc.sumstats.gz,pain_vonkorff_ldsc.sumstats.gz \
 --ref-ld-chr ~/ldscore/ldsc/refFiles/eur_w_ld_chr/ \
 --w-ld-chr ~/ldscore/ldsc/refFiles/eur_w_ld_chr/   \
 --out mdd_pain_vonkorff

~/ldscore/ldsc/ldsc.py \
 --rg mdd_ldsc.sumstats.gz,pain_grade3or4_ldsc.sumstats.gz \
 --ref-ld-chr ~/ldscore/ldsc/refFiles/eur_w_ld_chr/ \
 --w-ld-chr ~/ldscore/ldsc/refFiles/eur_w_ld_chr/   \
 --out mdd_pain_grade3or4

~/ldscore/ldsc/ldsc.py \
 --rg mdd_ldsc.sumstats.gz,pain_grade2or3or4_ldsc.sumstats.gz \
 --ref-ld-chr ~/ldscore/ldsc/refFiles/eur_w_ld_chr/ \
 --w-ld-chr ~/ldscore/ldsc/refFiles/eur_w_ld_chr/   \
 --out mdd_pain_grade2or3or4

