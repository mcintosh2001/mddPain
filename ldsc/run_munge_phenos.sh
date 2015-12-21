~/ldscore/ldsc/munge_sumstats.py \
 --sumstats /sdata/images/projects/GENSCOT/1/PGRS_DATASETS/PAIN_23ME/VonKorff_12thNov.BASIC.txt \
 --N 30000 \
 --out pain_vonkorff_ldsc \
 --merge-alleles ~/ldscore/ldsc/refFiles/w_hm3.snplist

~/ldscore/ldsc/munge_sumstats.py \
 --sumstats /sdata/images/projects/GENSCOT/1/PGRS_DATASETS/PAIN_23ME/Grades2or3or4_12thNov.BASIC.txt \
 --N 30000 \
 --out pain_grade2or3or4_ldsc \
 --merge-alleles ~/ldscore/ldsc/refFiles/w_hm3.snplist

~/ldscore/ldsc/munge_sumstats.py \
 --sumstats /sdata/images/projects/GENSCOT/1/PGRS_DATASETS/PAIN_23ME/Grades3or4_12thNov.BASIC.txt \
 --N 30000 \
 --out pain_grade3or4_ldsc \
 --merge-alleles ~/ldscore/ldsc/refFiles/w_hm3.snplist

~/ldscore/ldsc/munge_sumstats.py \
 --sumstats /sdata/images/projects/GENSCOT/1/PGRS_DATASETS/PGC_MDD/pgc.mdd.full.2012-04.txt \
 --N 18759 \
 --out mdd_ldsc \
 --merge-alleles ~/ldscore/ldsc/refFiles/w_hm3.snplist


