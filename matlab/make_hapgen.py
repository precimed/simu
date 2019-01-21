import pandas as pd
import numpy as np
import subprocess
import os
import time

#CHR     CHUNK   TSNP    FROM    TO
#1       1       13116   10616   1264907
#1       2       1265154 1264977 1896931
df = pd.read_table('chunks_EUR_noMHC.txt', delim_whitespace=True)

nK=100    # how many subjects to generate, in 1000
pop='EURqc'

# /work/users/oleksanf/20130502 - copy from /projects/NS9114K/1000Genome/ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502
# /work/users/oleksanf/20130502/recombination_map - recombination map from /projects/NS9114K/MMIL/1000Genome/phase3/build37_released/recombination_map 

dry_run = False
#dry_run=True

def execute_command(command):
    if dry_run:
        return 0 
    print("Execute command: {0}".format(command))
    #print(subprocess.Popen(command.split(), stdout=subprocess.PIPE, stderr=subprocess.STDOUT).communicate()[0].decode("utf-8"))
    #print(subprocess.check_output(command.split()).decode("utf-8"))
    exit_code = subprocess.call(command.split())
    print('Exit code: {}'.format(exit_code))
    return exit_code

class Params(object):
    def __init__(self, chri, chunki, loci_bp, from_bp, to_bp):
        self._chri = chri
        self._chunki = chunki
        self._loci_bp = loci_bp
        self._from_bp = from_bp
        self._to_bp = to_bp
    def out(self):
        return '/work/users/oleksanf/HAPGEN/{}_chr{}_chunk{}'.format(pop, self._chri, self._chunki)
    def describe(self):
        return 'chr{}_chunk{} ({}-{})'.format(self._chri, self._chunki,self._from_bp, self._to_bp)
    def vcftools_out(self):
        return "{out}.impute.hap".format(out=self.out())
    def vcftools_cmd(self):
        return """
module load zlib/1.2.8 && /usit/abel/u1/oleksanf/bin/vcftools \\
        --gzvcf /work/users/oleksanf/20130502/ALL.chr{chri}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz \\
        --chr {chri} \\
        --from-bp {frombp} \\
        --to-bp {tobp} \\
        --keep 1kG_EUR_subjlist.txt \\
        --max-alleles 2 \\
        --min-alleles 2 \\
        --max-missing 1 \\
        --mac 1 --maf 0.005 --max-maf 0.995 \\
        --IMPUTE \\
        --out {out} \\
        --phased \\
""".format(chri=self._chri, tobp=self._to_bp, frombp=self._from_bp, out=self.out())

    def hapgen2_cmd(self, repi):
        return """
/usit/abel/u1/oleksanf/hapgen2/hapgen2 \\
        -h {out}.impute.hap \\
        -l {out}.impute.legend \\
        -m /work/users/oleksanf/20130502/recombination_map/genetic_map_chr{chri}_combined_b37.txt \\
        -dl {loci} 1 2 4 -n {nK}000 1 -int {frombp} {tobp} -no_haps_output \\
        -o {out}_hapgen2_{nK}K_rep{repi}.gz \\
""".format(out=self.out(), repi=repi, nK=nK, chri=self._chri, chunki=self._chunki, loci=self._loci_bp, frombp=self._from_bp, tobp=self._to_bp)
    def hapgen2_out(self, repi):
        return "{out}_hapgen2_{nK}K_rep{repi}.controls.gen.gz".format(out=self.out(), repi=repi, chri=self._chri, chunki=self._chunki, nK=nK)
    def plink_cmd(self, repi):
        return """
module load plink && plink \\
        --gen {out}_hapgen2_{nK}K_rep{repi}.controls.gen.gz \\
        --make-bed --memory 2048 \\
        --out {out}_plink_{nK}K_rep{repi} \\
        --oxford-single-chr {chri} \\
        --sample {out}_hapgen2_{nK}K_rep{repi}.controls.sample \\
""".format(out=self.out(), repi=repi, chri=self._chri, chunki=self._chunki, nK=nK)
    def plink_out(self, repi):
        return "{out}_plink_{nK}K_rep{repi}.fam".format(out=self.out(), repi=repi, chri=self._chri, chunki=self._chunki, nK=nK)

nrep = 10

# generate FAM files
print('generate FAM files...')
fam_df = []
for repi in range(1, 1+nrep):
    with open('/work/users/oleksanf/HAPGEN/rep{repi}.fam'.format(repi=repi), 'w') as f: 
        for id in range(nK*1000*(repi-1), nK*1000*repi):
            f.write("id1_{id} id2_{id} 0 0 0 1\n".format(id=id))
print('done')

tasks_frq_1kg = []
for chri in df.CHR.unique():
    cmd = 'module load zlib/1.2.8 && /usit/abel/u1/oleksanf/bin/vcftools --gzvcf /work/users/oleksanf/20130502/ALL.chr{chri}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --freq --keep 1kG_EUR_subjlist.txt --max-alleles 2 --min-alleles 2 --max-missing 1 --phased --out /work/users/oleksanf/HAPGEN/1kG_EUR_chr{chri}'.format(chri=chri)
    out = '/work/users/oleksanf/HAPGEN/1kG_EUR_chr{chri}.frq'.format(chri=chri)
    tasks_frq_1kg.append(([cmd], [out]))

tasks_merge = []
print('generate merge files...')
for chri in  df.CHR.unique():
    mergelist = '/work/users/oleksanf/HAPGEN/merge_chr{}.txt'.format(chri)
    with open(mergelist, 'w') as f:
        for chunki in df[df.CHR==chri].CHUNK.values:
             for repi in range(1, 1+nrep):
                  fname = '/work/users/oleksanf/HAPGEN/{pop}_chr{chri}_chunk{chunki}_plink_{nK}K_rep{repi}'.format(pop=pop, chri=chri, chunki=chunki, nK=nK, repi=repi)
                  if chunki==1 and repi==1: fname1=fname
                  else: f.write(fname + '\n')
    plinkout = '/work/users/oleksanf/HAPGEN/chr{}'.format(chri)
    tasks_merge.append((['module load plink && plink --bfile {fname} --merge-list {mergelist} --out {out}'.format(fname=fname1, mergelist=mergelist, out=plinkout)], [plinkout + '.fam']))
print('done')

tasks = []; skiplist=[];
for index, row in df.iterrows():
    #if row['CHR'] != 22: continue
    p = Params(row['CHR'], row['CHUNK'], row['TSNP'], row['FROM'], row['TO'])
    cmd_local = []; out_local = []
    
    #cmd_local.append(p.vcftools_cmd()); out_local.append(p.vcftools_out())
    for repi in range(1,1+nrep):
        cmd_local.append(p.hapgen2_cmd(repi)); out_local.append(p.hapgen2_out(repi))
        cmd_local.append(p.plink_cmd(repi)); out_local.append(p.plink_out(repi))
        cmd_local.append("cp /work/users/oleksanf/HAPGEN/rep{repi}.fam {fam}".format(repi=repi, fam=p.plink_out(repi)))

    tasks.append((cmd_local, out_local))
    #break

tasks = tasks_merge
#tasks = tasks_frq_1kg
for (cmd_local, out_local) in tasks:
    if any([os.path.isfile(out) for out in out_local]):
        skiplist.append(cmd_local[0])
        continue

    command = '''#!/bin/bash
# Job name:
#SBATCH --job-name=run3ugTT
#
# Project:
#SBATCH --account=NN9114K
##SBATCH --account=uio
#
# Wall clock limit:
#SBATCH --time=48:00:00
#
## Set up job environment:
source /cluster/bin/jobsetup
module purge   # clear any inherited modules
set -o errexit # exit on errors
{}
'''.format('\n'.join(cmd_local))

    with open('run_script.sh', 'w') as f:
        f.write(command)
    delay = 2
    while True:
       if execute_command('sbatch --cpus-per-task=15  --mem-per-cpu=4096M run_script.sh') == 0:
           if not dry_run:
               subprocess.call('touch {}'.format(out_local[0]).split())
           break
       print('Error, wait {} sec and re-submit'.format(delay))
       time.sleep(delay)
       delay = min(delay * 2, 1200)

print('Submission complete.') 
print('\tTotal number of tasks: {}'.format(len(tasks)))
print('\tNumber of skipped tasks: {}'.format(len(skiplist))) 
