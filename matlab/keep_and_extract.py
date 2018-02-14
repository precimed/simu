import os
import glob
import ntpath
import subprocess
#workingpath = 'E:/'
workingpath = '/work/users/oleksanf/MMIL/cfan_SimuPop/EUR'
plink = 'plink'

# 10K individuals, 2M SNPs template
keep(os.path.join(workingpath, 'EUR_10K_2M'), 10000, '/home/oleksandr/2558411_ref.bim')
merge(os.path.join(workingpath, 'EUR_10K_2M'), os.path.join(workingpath, 'EUR_10K_2M_merged'))

# 500 individuals, all SNPs
keep(os.path.join(workingpath, 'EUR_500_80M'), 500)
merge(os.path.join(workingpath, 'EUR_500_80M'), os.path.join(workingpath, 'EUR_500_80M_merged'))

# 10K individuals, all SNPs
keep(os.path.join(workingpath, 'EUR_10K_80M'), 10000)

# 100K individuals, 1M SNPs template
keep(os.path.join(workingpath, 'EUR_100K_1M'), None, '/space/syn03/1/data/oleksandr/w_hm3.snplist/1m.ref')
merge(os.path.join(workingpath, 'EUR_100K_1M'), os.path.join(workingpath, 'EUR_100K_1M_merged'))


# Filter out SNPs and/or individuals
def keep(outputDir, numIndividuals=None, extract=None):
    if not os.path.exists(outputDir):
        os.makedirs(outputDir)
    if numIndividuals:
        with open('tmp.txt', 'w') as text_file:
            for i in range(numIndividuals):
                text_file.write('id1_{0} id2_{0}\n'.format(i))
    #files = glob.glob(os.path.join('/space/syn03/1/data/cfan/SimuPop/EUR', '*.bed'))
    #files = glob.glob(os.path.join(r'E:\EUR_100K_80M', '*.bed'))
    files = glob.glob(os.path.join(r'/work/users/oleksanf/MMIL/cfan_SimuPop/EUR', '*.bed'))
    #~/plink/plink --bfile chr1_chunk123 --keep keep10k.txt --make-bed --out chr1_chunk123_10k --extract ~/2558411_ref.bim
    for i, file in enumerate(files):
        print('Processing {0} out of {1}...\n'.format(i+1, len(files)))
        filename, file_extension = os.path.splitext(file)
        outfilename = os.path.join(outputDir, ntpath.basename(filename));
        command_extract = '--extract {0}'.format(extract) if extract else ''
        command_keep = '--keep tmp.txt' if numIndividuals else ''
        command = '{0} --memory 4096 --bfile {1} {2} {3} --make-bed --out {4}'.format(
            plink, filename, command_keep, command_extract, outfilename)
        print(command)
        subprocess.call(command.split())
    if numIndividuals: os.remove('tmp.txt')

# Merge together all bed files in a folder
def merge(inputDir, outputDir):
    if not os.path.exists(outputDir):
        os.makedirs(outputDir)
    files = glob.glob(os.path.join(inputDir, '*.bed'))
    first, _ = os.path.splitext(files[0])
    with open('mergelist.txt', 'w') as mergelist:
        for file in files[1:]:
            filename, file_extension = os.path.splitext(file)
            mergelist.write('{0}.bed {0}.bim {0}.fam\n'.format(filename))
    command = '{0} --memory 8192 --bfile {1} --merge-list mergelist.txt --allow-no-sex --make-bed --out {2}'.format(
        plink, first, os.path.join(outputDir, 'all'))
    subprocess.call(command.split())
    os.remove('mergelist.txt')

# ToDo: make .mat file with mapping between
# 1. Indices in the original 80M file
# 2. Indices in the resulting file
# 3. Indices in the reference file
# Now, (2) is subset of (1) and (3)

