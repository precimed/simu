#!/usr/bin/env python

# This is a python function to generate and submitt series of HAPGEN2 script
# Currently it only made to fit with mmil cluster environment
# 
# Required database - 1kGP phase 3 functional db 
# Chun C. Fan
#
# Adapter to Abel cluster environment by oleksanf.
# Before running ensure that all data exist on Abel server.
# mkdir /work/users/oleksanf/MMIL/
# mkdir /work/users/oleksanf/MMIL/1000Genome
# mkdir /work/users/oleksanf/MMIL/cfan_SQLdb
# mkdir /work/users/oleksanf/tmp
# rsync -avzP --update oleksandr@login.norstore.uio.no:/projects/NS9114K/MMIL/1000Genome/* /work/users/oleksanf/MMIL/1000Genome
# rsync -avzP --update oleksandr@login.norstore.uio.no:/projects/NS9114K/MMIL/cfan_SQLdb/* /work/users/oleksanf/MMIL/cfan_SQLdb
#
# download and unpack latest HapGen2 from https://mathgen.stats.ox.ac.uk/genetics_software/hapgen/hapgen2.html
# mkdir ~/hapgen2 && cd ~/hapgen2 && wget https://mathgen.stats.ox.ac.uk/genetics_software/hapgen/download/builds/x86_64/v2.2.0/hapgen2_x86_64.tar.gz



import argparse
import datetime
import numpy as np
from string import Template
import os
import re
import sqlite3
import sys

# Input parameters
parser = argparse.ArgumentParser()
parser.add_argument('--hapgen', action="store", dest="bstem", default='/usit/abel/u1/oleksanf/hapgen2/hapgen2', help="executable binary of HAPGEN2")
parser.add_argument('--gen', action="store", dest="gstem",default='/work/users/oleksanf/MMIL/1000Genome/phase3/build37_released/EUR', help="hap and legend file folders")
parser.add_argument('--recomb', action="store", dest="rstem",default='/work/users/oleksanf/MMIL/1000Genome/phase3/build37_released/recombination_map', help="recombination maps folders")
parser.add_argument('--pop', action="store", dest="popf",default='EUR', help="Type of allele frequencies, EUR, EAS, AMR, or AFR")
parser.add_argument('--controls', dest="controls", default=10000, action="store", help="Number of controls we wish to generate")
parser.add_argument('--chunk_int', dest="chunk", default=10000, action="store", help="Chunk size for each output")
parser.add_argument('--db', action="store", dest="db", default="/work/users/oleksanf/MMIL/cfan_SQLdb/Functional_Annot_1kGP_hg19.db", help="SQL database use for allele frequencies and physical positions")
parser.add_argument('--out', action="store", dest="output", default="/work/users/oleksanf/hapmap2_results", help="output absolute path and prefix")
parser.add_argument('--tmp', action="store", dest="tmp", default="/work/users/oleksanf/tmp", help="place for temporary script shells")
args = parser.parse_args()

# immutable parameters
HAPBIN = args.bstem
GSTEM = args.gstem
RSTEM = args.rstem
DBFILE = args.db
CONNUM = args.controls
CHUNKN = int(args.chunk)
OUTSTEM = args.output
TMPDIR = args.tmp
POPTMP = args.popf
NUSIANCEPARA = [1,2,4] # nusiance parameter because no use for case alleles

def get_file(chrnum):
  pattern = '(?=chr' + str(chrnum) + '_)\w+'
  INHAP = [x for x in sorted(os.listdir(GSTEM)) if re.search(pattern,x) and x.endswith('.hap')]
  INLEGEND = [x for x in sorted(os.listdir(GSTEM)) if re.search(pattern,x) and x.endswith('.legend')]
  INRECOMB = [x for x in sorted(os.listdir(RSTEM)) if re.search(pattern,x) and x.endswith('.txt')]
  if not INHAP: INHAP = ['xxx.INHAP.hap']
  if not INLEGEND: INLEGEND=['xxx.INLEGEND.legend']
  if not INRECOMB: INRECOMB=['xxx.INRECOMB.txt']
  return GSTEM + '/' + INHAP[0], GSTEM + '/' + INLEGEND[0], RSTEM + '/' + INRECOMB[0]

def generate_cmd(INHAP, INLEGEND, INRECOMB, pos_SNP, pos_INI, pos_END, OUTFILE):
  params = dict()
  params['HAPBIN'] = HAPBIN
  params['HAP'] = INHAP
  params['LEG'] = INLEGEND
  params['RE'] = INRECOMB
  params['SNP'] = pos_SNP
  params['INI'] = pos_INI
  params['END'] = pos_END
  params['CON'] = CONNUM
  params['OUT'] = OUTFILE
  params['NU'] = str(NUSIANCEPARA[0]) + ' ' + str(NUSIANCEPARA[1]) + ' ' + str(NUSIANCEPARA[2])
  cmd = Template('$HAPBIN -h $HAP -l $LEG -m $RE -dl $SNP $NU -n $CON 1 -int $INI $END -no_haps_output -o $OUT')
  cmd = cmd.substitute(params)
  return cmd

def write_tmp(tmp_fname,cmd):
  outFile = open(tmp_fname, 'w')
  outFile.write('echo\n')
  outFile.write(cmd  + '\n')
  outFile.write('echo\n')
  outFile.close()

class WeightsDB:
  def __init__(self):
    self.conn = sqlite3.connect(DBFILE)

  def query(self, sql, args=None):
    c = self.conn.cursor()
    if args:
      for ret in c.execute(sql, args):
        yield ret
    else:
      for ret in c.execute(sql):
        yield ret

class GetParametersOf:
  def __init__(self):
    self.db = WeightsDB()

  def __call__(self, chrnum):
    print("{now} Generating parameters for chromosome {chr}".format(now=datetime.datetime.now(), chr=chrnum))
    t = (str(chrnum),)
    FNAME = POPTMP + '_AF'
    tmp = [(tup[0],tup[1],tup[2]) for tup in self.db.query("SELECT pos, rsid, %s FROM functionDB WHERE chrnum=?" % (FNAME), t)]
    c_array = range(1, len(tmp), CHUNKN)
    c_array.append(len(tmp))
    for c in range(1,len(c_array)):
      tsnp = [w[0] for w in tmp[int(c_array[c-1]):int(c_array[c])-1] if w[2] >= 0.05 and w[2] < 0.95]
      pos_ini = tmp[int(c_array[c-1])][0]
      pos_end = tmp[int(c_array[c])-1][0]
      yield tsnp[0], pos_ini, pos_end
   
get_para = GetParametersOf()

#
# Execute part
# Currently simulate all chromosomes from 1kGP in the mmil cluster 4

f = open('chunks.txt', 'w')
f.write('CHR\tCHUNK\tTSNP\tFROM\tTO\n')
for i in range(1,23):
  INHAP, INLEGEND, INRECOMB = get_file(i)  
  k = 1
  for chunki in get_para(i):
    simuout = OUTSTEM + '_chr' + str(i) + '_chunk' + str(k) + '.gz' 
    cmd0 = generate_cmd(INHAP, INLEGEND, INRECOMB, chunki[0], chunki[1], chunki[2], simuout)
    tmp_fname = TMPDIR + '/tmp_script_chr' + str(i) + '_chunk' + str(k) + '.sh'
    #write_tmp(tmp_fname, cmd0)
    tmp_log = TMPDIR + '/tmp_script_chr' + str(i) + '_chunk' + str(k)
    cmdsub = 'qsub -hard -l h_vmem=32G -e ' + tmp_log + '_1' + ' -o ' + tmp_log + '_2' + ' -N Simu_chr' + str(i) + '_' + str(k) + ' ' + tmp_fname + '\n'
    #os.system('ssh -XY 169.228.56.67 ' + cmdsub)
    f.write('{}\t{}\t{}\t{}\t{}\n'.format(i, k, chunki[0], chunki[1], chunki[2]))
    print(cmdsub)
    k += 1
f.close()



