import pandas as pd
import os.path
from glob import glob
import re

# For testing, parse the yaml file (this is automatically done my Snakemake)
import yaml
with open("config/config.yaml", 'r') as ymlfile: config = yaml.load(ymlfile)

configfile: "config/config.yaml"

# single part directory names are given in config.yaml
# this puts together the ones which are composed of previous values
config['moviedir']  = "{seqdir}/rawmovie".format_map(config)
config['ccsdir']    = "{seqdir}/ccs".format_map(config)
config['fastqdir']  = "{seqdir}/rawfastq".format_map(config)
config['demuxdir']  = "{seqdir}/demux".format_map(config)
config['trimdir']   = "{seqdir}/trim".format_map(config)
config['regiondir'] = "{seqdir}/regions".format_map(config)
config['filterdir'] = "{seqdir}/filter".format_map(config)
config['tagdir']    = "{labdir}/tags".format_map(config)

config['dataset']    = '{labdir}/datasets.csv'.format_map(config)
config['regions']    = '{labdir}/regions.csv'.format_map(config)
config['gits7_tags'] = '{labdir}/Hectors_tag_primer_plates.xlsx'.format_map(config)
config['lr5_tags']   = '{labdir}/Brendan_soil2.xlsx'.format_map(config)

# load the dataset and region definitions
datasets = pd.read_csv(config['dataset']).set_index('Seq.Run', drop = False)
regions = pd.read_csv(config['regions']).set_index('Region')

#### PacBio conversions ####
# find the PacBio movie files
seqplates = [seqrun + '_' + '{:03d}'.format(i + 1) for seqrun, runs in datasets['Runs'].iteritems() for i in range(runs)]
moviefiles = {}
for sp in seqplates:
    if sp.startswith('pb'):
        moviefiles[sp] = []
        for m in glob("{rawdir}/**/rawdata/{sp}/**/Analysis_Results/*.bas.h5".format(rawdir = config['rawdir'], sp = sp),
              recursive = True):
            bn = os.path.basename(m)
            moviefiles[sp].append(re.sub(r"\.bas\.h5", "", bn))

datasets2 = (datasets['Runs']
             .map(lambda n: ['{:03d}'.format(i + 1) for i in range(n)])
             .apply(pd.Series, name="Plate")
             .stack())

datasets2 = (datasets.join(pd.DataFrame({"Plate" : datasets2}))
             .assign(SeqPlate=lambda x: x['Seq.Run'].str.cat(x.Plate, sep = "_"))
             .set_index('SeqPlate', drop = False))

# convert a movie from PacBio BAX format to BAM format
# and put output in the raw movie directory

rule convertmovies:
    input:
        expand("{moviedir}/{movie}.{type}.bam",
               moviedir=config['moviedir'],
               movie = [m for sp in moviefiles.values() for m in sp],
               type = ['subreads', 'scraps'])

rule bax2bam:
    input:
        lambda wildcards: glob("{rawdir}/**/rawdata/**/Analysis_Results/{wildcards.movie}.*.h5".format(
                                 wildcards = wildcards,
                                 rawdir = config['rawdir']),
                               recursive = True)
    output:
        expand("{moviedir}/{{movie}}.{movietype}.bam",
               moviedir = config['moviedir'],
               movietype = ["subreads", "scraps"])
    params:
        prefix="{moviedir}/{{movie}}".format_map(config)
    log: "{logdir}/bax2bam_{{movie}}.log".format_map(config)
    shell:
         "bax2bam {input} -o {params.prefix} &> {log}"

rule ccs:
    output:
          "{ccsdir}/{{movie}}.ccs.bam".format_map(config)
    input:
         "{moviedir}/{{movie}}.subreads.bam".format_map(config)
    log: "{logdir}/ccs_{{movie}}.log".format_map(config)
    shell:
         "ccs --polish {input} &>{log}"

rule ccs2fastq:
    output:
        "{fastqdir}/{{seqplate}}.fastq.gz".format_map(config)
    input:
         lambda wildcards: expand("{ccsdir}/{movie}.ccs.bam",
                                  ccsdir = config['ccsdir'],
                                  movie = moviefiles[wildcards.seqplate])
    log: "{logdir}/ccs2fastq_{{seqplate}}.log".format_map(config)
    shell:
         "samtools cat {input} | samtools fastq -0 {output} - &>{log}"

rule pacbio_fastq:
    input:
        expand("{fastqdir}/{seqplate}.fastq.gz",
               fastqdir = config['fastqdir'],
               seqplate = [sp for sp in seqplates if sp.startswith('pb')])

rule tagfiles:
    input:
        gits7_tags = config['gits7_tags'],
        lr5_tags = config['lr5_tags'],
        dataset = config['dataset']
    output:
        expand("{tagdir}/{seqrun}.fasta", tagdir = config['tagdir'], seqrun = datasets['Seq.Run'])
    log: "{logdir}/tagfiles.log".format_map(config)
    script:
        "{config[rdir]}/tags.extract.R"

def find_barcode(wildcards): "{tagdir}/{seqrun}.fasta".format(tagdir = config['tagdir'],
                                                              seqrun = datasets2.loc[wildcards['seqplate'], 'Seq.Run'])

checkpoint pacbio_demux:
    input:
        fastq = "{fastqdir}/{{seqplate}}.fastq.gz".format_map(config),
        barcode = find_barcode
    output:
        directory("{trimdir}/{{seqplate}}/".format_map(config))
    params:
        fpattern = lambda wildcards: "{trimdir}/{seqplate}/{seqplate}-{{name}}f.trim.fastq.gz".format(trimdir = config['trimdir'],
                                                                                           seqplate = wildcards['seqplate']),
        rpattern = lambda wildcards: "{trimdir}/{seqplate}/{seqplate}-{{name}}r.trim.fastq.gz".format(trimdir = config['trimdir'],
                                                                                           seqplate = wildcards['seqplate'])
    log: "{logdir}/pacbio_demux_{{seqplate}}.log".format_map(config)
    shell:
         """
         mkdir -p {output} &&
         cutadapt --quiet\\
           -g file:{input.barcode}\\
           -m 1\\
           --untrimmed-output=-\\
           -o {params.fpattern}\\
           {input.fastq} |
         fastx_reverse_complement |
         cutadapt --quiet\\
           -g file:{input.barcode}\\
           -m 1\\
           --trimmed-only\\
           -o {params.rpattern}\\
           - 2>{log}
         """

def find_ion_bam(wildcards):
    platekey = (pd.read_csv(os.path.join(config['labdir'], datasets.loc[wildcards.seqrun].PlateKey))
                .assign(id = lambda x: x['tag.fwd'].str.rsplit('-').str.get(1)))
    num = platekey.loc[lambda x: x['well'] == wildcards.well, 'id'].get_values()[0]
    dset = datasets.loc[lambda x: x['Seq.Run'] == wildcards.seqrun, 'Dataset'][0]
    return glob("{rawdir}/{dset}/{seqrun}/rawdata/*/IonXpress_{num}*.bam".format(rawdir = config['rawdir'],
                                                                                 dset = dset,
                                                                                 seqrun = wildcards.seqrun,
                                                                                 num = '{:03d}'.format(int(num))))

rule bam2fastq:
    output: "{demuxdir}/{{seqrun}}_{{plate,\d+}}-{{well,[A-H]\d+}}.demux.fastq.gz".format_map(config)
    input: find_ion_bam
    log: "{logdir}/bam2fastq_{{seqrun}}_{{plate}}-{{well}}.log".format_map(config)
    shell:
        """
        mkdir -p {config[demuxdir]}
        samtools fastq {input} | gzip >>{output} 2>{log}
        """

rule ion_trim:
    output: "{trimdir}/{{seqrun}}_{{plate,\d+}}-{{well,[A-H]\d+}}.trim.fastq.gz"
    input:
        fastq = "{demuxdir}/{seqrun}_{plate}-{well}.demux.fastq.gz",
        barcode = find_barcode
    log: "{logdir}/ion_trim_{{seqrun}}_{{plate}}-{{well}}.log".format_map(config)
    shell:
        """
        mkdir -p {config[trimdir]}
        cutadapt --trimmed-only\
                 -m 1\
	             -g file:{input.barcode}\
	             -j {threads}\
	             -o {output}\
	             {input.fastq}\
	             &> {log}
        """



def demux_find(seqplate):
    def subfind(wildcards):
        checkpoints.pacbio_demux.get(seqplate = seqplate)
        return glob("{trimdir}/{seqplate}/{seqplate}-*.trim.fastq.gz".format(trimdir = config['trimdir'],
                                                                   seqplate = seqplate))
    return subfind

rule pacbio_demuxall:
    input:
         demux_find('pb_500_001'),
         demux_find('pb_500_002'),
         demux_find('pb_483_001'),
         demux_find('pb_483_002')