import os, sys
from posixpath import join

def getlist(targetdir):
    samplesdict = {}
    for i in os.listdir(targetdir):
        if i.endswith('.TSV'):
            target = os.path.join(targetdir, i)
            name = i.split('.')[0]
            samplesdict[target] = name
    keyss = list(samplesdict.keys())
    seprator = ' -I '
    inputs = seprator.join(keyss)
    final = seprator + inputs
    return final, samplesdict

def filterintervals(inputs, PL, annotated, wd):
    output = os.path.join(wd, 'cohort_filtered.interval_list')
    gatk = 'gatk FilterIntervals'
    imr = '-imr OVERLAPPING_ONLY'
    cmd = f'{gatk} -L {PL} --annotated-intervals {annotated} {imr} -O {output}{inputs}'
    os.system(cmd)
    return output

def determinecontig(inputs, filteredL, cpp, wd):
    gatk = 'gatk DetermineGermlineContigPloidy'
    imr = '-imr OVERLAPPING_ONLY'
    optin = f'--output {wd} --output-prefix ploidy --verbosity DEBUG'
    cmd = f'{gatk} -L {filteredL} {imr}{inputs} --contig-ploidy-priors {cpp} {optin}'
    os.system(cmd)
    return os.path,join(wd, 'ploidy-calls')

def GermlineCNVCaller(inputs, filteredL, ploidy, annotated, wd):
    gatk = 'gatk GermlineCNVCaller --run-mode COHORT'
    imr = '-imr OVERLAPPING_ONLY'
    contigcalls = f'--contig-ploidy-calls {ploidy}'
    AnnL = f'--annotated-intervals {annotated}'
    other = f'--output {wd} --output-prefix Cohort --verbosity DEBUG'
    cmd = f'{gatk} -L {filteredL} {imr}{inputs} {contigcalls} {AnnL} {other}'
    os.system(cmd)

def scatterintervals():
    gatk = 'java -jar /home/edi/ngs/gatk/gatk-package-4.2.2.0-local.jar IntervalListTools'
    filteredL = '-I /home/edi/ngs/c/output_run69/cohort_filtered.interval_list'
    other = '--SUBDIVISION_MODE INTERVAL_COUNT --SCATTER_CONTENT 10000 --OUTPUT scatter'
    cmd = f'{gatk} {filteredL} {other}'
    os.system(cmd)


if _name_ == '_main_':
    wd = sys.argv[1]
    PL = '/home/gnks/Desktop/edi/c/bundle/Preproceced.interval_list'
    annotated = '/home/gnks/Desktop/edi/c/bundle/Preproceced_GC_Annotated.interval_list'
    cpp = '/home/gnks/Desktop/edi/c/bundle/contig_ploidy_priors.tsv'
    inputs, sampledicts = getlist(wd)
    print(inputs)
    # filteredL = filterintervals(inputs, PL, annotated, wd)
    # ploidy = determinecontig(inputs, filteredL, cpp, wd)
    # GermlineCNVCaller(inputs, filteredL, ploidy, annotated, wd)