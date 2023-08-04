#rmats-turbo/4.1.1
~
gtf="Human/Homo_sapiens.GRCh38.108.chr.gtf"
outdir="SMNDC1/ENCODE/rMATS"

#SMNDC1 Data inhouse
sbatch -J Pancreas -p general -t 4:00:00 --mem=32g -N 1 -n 16 -e Pancreas.e -o Pancreas.o --wrap="run_rmats --b1 Pancreas_C.txt --b2 Pancreas_T.txt --gtf ${gtf} -t paired --readLength 151 --nthread 16 -t paired --od ${outdir}/Pancreas --tmp ${outdir}/Pancreas_tmp"
sbatch -J Ovarian -p general -t 4:00:00 --mem=32g -N 1 -n 16 -e Ovarian.e -o Ovarian.o --wrap="run_rmats --b1 Ovarian_C.txt --b2 Ovarian_T.txt --gtf ${gtf} -t paired --readLength 151 --nthread 16 -t paired --od ${outdir}/Ovarian --tmp ${outdir}/Ovarian_tmp"

#SMNDC1 Data ENCODE
sbatch -J SMNDC1_HepG2 -p general -t 4:00:00 --mem=120g -N 1 -n 16 -e SMNDC1_HepG2.e -o SMNDC1_HepG2.o --wrap="run_rmats --b1 CtrlSMNDC1_HepG2.txt --b2 shSMNDC1_HepG2.txt --gtf ${gtf} -t paired --readLength 100 --nthread 16 -t paired --od ${outdir}/SMNDC1_HepG2 --tmp ${outdir}/SMNDC1_HepG2_tmp"
sbatch -J SMNDC1_K562 -p general -t 4:00:00 --mem=120g -N 1 -n 16 -e SMNDC1_K562.e -o SMNDC1_K562.o --wrap="run_rmats --b1 CtrlSMNDC1_K562.txt --b2 shSMNDC1_K562.txt --gtf ${gtf} -t paired --readLength 100 --nthread 16 -t paired --od ${outdir}/SMNDC1_K562 --tmp ${outdir}/SMNDC1_K562_tmp"
