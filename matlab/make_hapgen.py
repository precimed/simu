import pandas as pd
import numpy as np
import subprocess
import os

#CHR     CHUNK   TSNP    FROM    TO
#1       1       13116   10616   1264907
#1       2       1265154 1264977 1896931
df = pd.read_table('chunks_EUR.txt', delim_whitespace=True)


# /work/users/oleksanf/20130502 - copy from /projects/NS9114K/1000Genome/ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502

dry_run = False

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
    def __init__(self, pop, chri, chunki, loci_bp, from_bp, to_bp):
        self._chri = chri
        self._chunki = chunki
        self._loci_bp = loci_bp
        self._from_bp = from_bp
        self._to_bp = to_bp
        self._pop = pop
    def out(self):
        return '/work/users/oleksanf/HAPGEN/{}_chr{}_chunk{}'.format(self._pop, self._chri, self._chunki)
    def describe(self):
        return 'chr{}_chunk{} ({}-{})'.format(self._chri, self._chunki,self._from_bp, self._to_bp)
    def vcftools_out(self):
        return "{out}.impute.hap".format(out=self.out())
    def vcftools_cmd(self):
        return """
/usit/abel/u1/oleksanf/bin/vcftools \\
        --gzvcf /work/users/oleksanf/20130502/ALL.chr{chri}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz \\
        --chr {chri} \\
        --from-bp {frombp} \\
        --to-bp {tobp} \\
        --keep 1kG_EUR_subjlist.txt \\
        --max-alleles 2 \\
        --min-alleles 2 \\
        --max-missing 1 \\
        --IMPUTE \\
        --out {out} \\
        --phased \\
""".format(chri=self._chri, tobp=self._to_bp, frombp=self._from_bp, out=self.out())

tasks = []; skiplist=[];
for index, row in df.iterrows():
    if row['CHR'] != 21: continue
    p = Params('EUR', row['CHR'], row['CHUNK'], row['TSNP'], row['FROM'], row['TO'])
    tasks.append((p.vcftools_cmd(), p.vcftools_out()))
   
for (cmd, out) in tasks:
    if os.path.isfile(out):
        skiplist.append(cmd)
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
#SBATCH --time=4:00:00
#
## Set up job environment:
source /cluster/bin/jobsetup
module purge   # clear any inherited modules
set -o errexit # exit on errors
{}
'''.format(cmd)

    with open('run_script.sh', 'w') as f:
        f.write(command)
    delay = 2
    while True:
       if execute_command('sbatch --cpus-per-task=1  --mem-per-cpu=4096M run_script.sh') == 0:
           if not dry_run:
               subprocess.call('touch {}'.format(out).split())
           break
       print('Error, wait {} sec and re-submit'.format(delay))
       time.sleep(delay)
       delay = min(delay * 2, 1200)

print('Submission complete.') 
print('\tTotal number of tasks: {}'.format(len(tasks)))
print('\tNumber of skipped tasks: {}'.format(len(skiplist))) 
