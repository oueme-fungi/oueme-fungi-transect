# Snakemake workflow file
# This workflow identifies raw sequencing reads and performs command-line based
# operations like circular consensus calling and demultiplexing.
# Then it calls R scripts to build a plan in Drake, an R workflow manager.

# Drake developed a lot while I was working on this project, some of the
# more convoluted parts of the combined workflow were intended to use Snakemake
# to compensate for features that weren't yet implemented in Drake.

import pandas as pd
import os.path
from glob import glob
import re
import subprocess
from snakemake.utils import listfiles

# For testing, parse the yaml file (this is automatically done by Snakemake)
#import yaml
#with open("config/config.yaml", 'r') as ymlfile: config = yaml.safe_load(ymlfile)

configfile: "config/config.yaml"

# Find the maximum number of cores available to a single node on SLURM,
# or if we aren't in a SLURM environment, how many we have on the local machine.
try:
    maxthreads = max([int(x) for x in re.findall(r'\d+', subprocess.check_output(["sinfo", "-O", "cpus"]).decode())])
except FileNotFoundError:
    maxthreads = int(subprocess.check_output("nproc").decode())

# load the dataset and region definitions
datasets = pd.read_csv(config['dataset']).set_index('seq_run', drop = False)
regions = pd.read_csv(config['regions']).set_index('region')

#### PacBio conversions ####
# find the PacBio movie files
seqplates = [seqrun + '_' + '{:03d}'.format(i + 1) for seqrun, runs in datasets['runs'].iteritems() for i in range(runs)]
moviefiles = {}
for sp in seqplates:
    if sp.startswith('pb'):
        moviefiles[sp] = []
        for m in glob("{rawdir}/**/rawdata/{sp}/**/Analysis_Results/*.bas.h5".format(rawdir = config['rawdir'], sp = sp),
              recursive = True):
            bn = os.path.basename(m)
            moviefiles[sp].append(re.sub(r"\.bas\.h5", "", bn))

datasets2 = (datasets['runs']
             .map(lambda n: ['{:03d}'.format(i + 1) for i in range(n)])
             .apply(pd.Series, name="plate")
             .stack())

datasets2 = (datasets.join(pd.DataFrame({"plate" : datasets2}))
             .assign(seqplate=lambda x: x['seq_run'].str.cat(x.plate, sep = "_"))
             .set_index('seqplate', drop = False))

# endpoint target
localrules: all
rule all:
    input:
        "{rmddir}/transect_paper.pdf".format_map(config),
        "{rmddir}/transect_supplement.pdf".format_map(config),
        "{outdir}/supp_file_1.tsv".format_map(config),
        "{outdir}/supp_file_2.tsv".format_map(config),
        "{outdir}/supp_file_3.pdf".format_map(config)
        #{outdir}/tech_compare.pdf".format_map(config)

localrules: repair
rule repair:
    input: "{rdir}/repair.R".format_map(config)
    threads: 1
    conda: "config/conda/drake.yaml"
    script: "{rdir}/repair.R".format_map(config)

localrules: clean
rule clean:
    input: "{rdir}/clean.R".format_map(config)
    threads: 1
    conda: "config/conda/drake.yaml"
    script: "{rdir}/clean.R".format_map(config)

# endpoint target: convert all pacbio movies to Sequel format
localrules: convertmovies
rule convertmovies:
    input:
        expand("{moviedir}/{movie}.{type}.bam",
               moviedir=config['moviedir'],
               movie = [m for sp in moviefiles.values() for m in sp],
               type = ['subreads', 'scraps'])

# convert a raw RSII-format (.h5) movie to the Sequel format (.bam)
rule bax2bam:
    output:
        temp(expand("{moviedir}/{{movie}}.{movietype}.bam",
               moviedir = config['moviedir'],
               movietype = ["subreads", "scraps"]))
    input:
        lambda wildcards: glob("{rawdir}/**/rawdata/**/Analysis_Results/{wildcards.movie}.*.bax.h5"
                               .format(wildcards = wildcards,
                                       rawdir = config['rawdir']),
                               recursive = True)
    shadow: "shallow"
    conda: "config/conda/pacbio.yaml"
    group: "pacbio"
    params:
        prefix="{moviedir}/{{movie}}".format_map(config)
    resources:
        walltime=5
    log: "{logdir}/bax2bam_{{movie}}.log".format_map(config)
    shell:
         "bax2bam {input} -o {params.prefix} &> {log}"

# generate a circular consensus sequence from raw PacBio reads
rule ccs:
    output:
          temp("{ccsdir}/{{movie}}.ccs.bam".format_map(config))
    input:
         "{moviedir}/{{movie}}.subreads.bam".format_map(config)
    resources:
        walltime=120
    shadow: "shallow"
    conda: "config/conda/pacbio.yaml"
    group: "pacbio"
    threads: 4
    log: "{logdir}/ccs_{{movie}}.log".format_map(config)
    shell:
         "ccs --numThreads {threads} {input} {output} &>{log}"

# convert all circular consensus files for a sequencing run to fastq.gz format
localrules: ccs2fastq
rule ccs2fastq:
    output:
        "{fastqdir}/{{seqplate}}.fastq.gz".format_map(config)
    input:
         lambda wildcards: expand("{ccsdir}/{movie}.ccs.bam",
                                  ccsdir = config['ccsdir'],
                                  movie = moviefiles[wildcards.seqplate])
    conda: "config/conda/pacbio.yaml"
    params:
        prefix = "{fastqdir}/{{seqplate}}".format_map(config)
    group: "pacbio"
    log: "{logdir}/ccs2fastq_{{seqplate}}.log".format_map(config)
    resources:
        walltime=10
    shell:
         "bam2fastq -o {params.prefix} {input} &>{log}"

# endpoint target: generate all PacBio fastq files
rule pacbio_fastq:
    input:
        expand("{fastqdir}/{seqplate}.fastq.gz",
               fastqdir = config['fastqdir'],
               seqplate = [sp for sp in seqplates if sp.startswith('pb')])

# generate a fasta file of primer sequences for each plate in the format required by cutadapt
localrules: tagfiles
rule tagfiles:
    input:
        gits7_tags = config['gits7_tags'],
        lr5_tags = config['lr5_tags'],
        dataset = config['dataset'],
        script = "{rdir}/tags.extract.R".format_map(config)
    output:
        expand(
            "{tagdir}/{seqrun}.fasta",
            tagdir = config['tagdir'],
            seqrun = datasets["seq_run"][datasets.tech != "Illumina"]
        ),
        expand(
            "{tagdir}/{seqrun}_R{dir}.fasta",
            tagdir = config['tagdir'],
            seqrun = datasets["seq_run"][datasets.tech == "Illumina"],
            dir = [1, 2]
        )
    resources:
        walltime=5
    conda: "config/conda/drake.yaml"
    log: "{logdir}/tagfiles.log".format_map(config)
    script:
        "{rdir}/tags.extract.R".format_map(config)

# look up the primers/barcodes file based on the plate ID.
def find_barcode(wildcards):
    return ("{tagdir}/{seqrun}.fasta"
            .format(tagdir = config['tagdir'],
                    seqrun = datasets2.loc[wildcards['seqplate'], 'seq_run']))

# demultiplex and trim primers from a PacBio circular consensus in fastq.gz format
# The number of output files is not known a priori, so the declared output is a directory and this is a checkpoint
localrules: pacbio_demux
checkpoint pacbio_demux:
    input:
        fastq = "{fastqdir}/{{seqplate}}.fastq.gz".format_map(config),
        barcode = find_barcode
    output:
        directory("{trimdir}/{{seqplate}}/".format_map(config))
    params:
        fpattern = lambda wildcards: ("{trimdir}/{seqplate}/{seqplate}-{{name}}f.trim.fastq.gz"
                                      .format(trimdir = config['trimdir'],
                                              seqplate = wildcards['seqplate'])),
        rpattern = lambda wildcards: ("{trimdir}/{seqplate}/{seqplate}-{{name}}r.trim.fastq.gz"
                                      .format(trimdir = config['trimdir'],
                                              seqplate = wildcards['seqplate']))
    resources:
        walltime=5
    conda: "config/conda/demultiplex.yaml"
    group: "pacbio"
    threads: 3 #2 instances of cutadapt and 1 fastx_reverse_complement
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

# Return a closure which calls checkpoints.pacbio_demux.get() to indicate to Snakemake that this rule is
# dynamically calculated after completion of pacbio_demux, and then find all the demultiplexed pacbio files
# for the given plate.
def pacbio_find(seqplate):
    def subfind(wildcards):
        checkpoints.pacbio_demux.get(seqplate = seqplate)
        return glob("{trimdir}/{seqplate}/{seqplate}-*.trim.fastq.gz"
                    .format(trimdir = config['trimdir'],
                            seqplate = seqplate))
    return subfind

# endpoint rule to demultiplex all the PacBio reads.
localrules: pacbio_demuxall
rule pacbio_demuxall:
    input:
         pacbio_find('pb_500_001'),
         pacbio_find('pb_500_002'),
         pacbio_find('pb_483_001'),
         pacbio_find('pb_483_002')


# put all the demultiplexed long PacBio reads in one file, with labeled samples
localrules: pacbio_singledemux
rule pacbio_singledemux:
    input:
        pacbio_find('pb_500_001'),
        pacbio_find('pb_500_002')
    output:
        fastq = "{seqdir}/long_pacbio.fastq.gz".format_map(config),
        key   = "{seqdir}/long_pacbio.key.tsv".format_map(config)
    shell:
        """
        mkdir -p $(dirname {output}) &&
        rm -f {output.fastq} &&
        touch {output.fastq} &&
        for file in {input}; do
            sample=$(basename $file) &&
            sample=${{sample#pb_500_}} &&
            sample=${{sample%[fr].trim.fastq.gz}} &&
            echo $sample &&
            zcat $file |
            sed '/@/ s/@.*/&;sample:'${{sample}}'/' |I made an LSU alignment in MUSCLE
            gzip - >>{output.fastq}
        done &&
        zcat {output.fastq} | sed -nr '/^@/s/^@([^;]+;sample:(.+))/\1\t\2/ p' >{output.key}
        """


#### IonTorrent files ####
# The IonTorrent files have already been demultiplexed by the IonTorrent pipeline.
# They only need to be renamed, converted to .fastq.gz, and have the primers trimmed away.

# look up the key to convert between barcode numbers and plate wells for an iontorrent run
def ion_platekey(seqrun):
    return (pd.read_csv(os.path.join(config['labdir'], datasets.loc[seqrun].plate_key))
            .assign(id = lambda x: pd.to_numeric(x['tag_fwd'].str.rsplit('-').str.get(1))))

# Find the demultiplexed .bam file (named by barcode index) which corresponds to a particular sample (named by well)
def find_ion_bam(wildcards):
    # load the relevant platekey
    platekey = ion_platekey(wildcards.seqrun)
    # find the id which corresponds to the desired
    num = platekey.loc[lambda x: x['well'] == wildcards.well, 'id'].array[0]
    # find the name of the dataset
    dset = datasets.loc[lambda x: x['seq_run'] == wildcards.seqrun, 'dataset'][0]
    return glob("{rawdir}/{dset}/{seqrun}/rawdata/*/IonXpress_{num}*.bam".format(rawdir = config['rawdir'],
                                                                                 dset = dset,
                                                                                 seqrun = wildcards.seqrun,
                                                                                 num = '{:03d}'.format(int(num))))

localrules: bam2fastq, ion_trim
# convert a demultiplexed IonTorrent .bam file to .fastq.gz
rule bam2fastq:
    output: temp("{demuxdir}/{{seqrun}}_{{plate,\d+}}-{{well,[A-H]\d+}}.demux.fastq.gz".format_map(config))
    input: find_ion_bam
    resources:
        walltime=5
    conda: "config/conda/demultiplex.yaml"
    group: "iontorrent"
    log: "{logdir}/bam2fastq_{{seqrun}}_{{plate}}-{{well}}.log".format_map(config)
    shell:
        """
        mkdir -p {config[demuxdir]}
        samtools fastq {input} | gzip >>{output} 2>{log}
        """

# Trim primers from an IonTorrent .fastq.gz file
rule ion_trim:
    output: "{trimdir}/{{seqplate}}/{{seqplate}}-{{well,[A-H]\d+}}.trim.fastq.gz".format_map(config)
    input:
        fastq = "{demuxdir}/{{seqplate}}-{{well}}.demux.fastq.gz".format_map(config),
        barcode = find_barcode
    resources:
        walltime=5
    threads: 4
    conda: "config/conda/demultiplex.yaml"
    group: "iontorrent"
    log: "{logdir}/ion_trim_{{seqplate}}-{{well}}.log".format_map(config)
    shell:
        """
        mkdir -p {config[trimdir]}
        [ $(wc -c < {input.fastq}) -gt 50 ] &&
        cutadapt --trimmed-only\
                 -m 1\
                     -g file:{input.barcode}\
                     -j {threads}\
                     -o {output}\
                     {input.fastq}\
                     &> {log} ||
        cp {input.fastq} {output}
        """

# Generate the names of all the trimmed IonTorrent files which will be generated from the available .bam files
def ion_find(seqrun, plate):
    platekey = ion_platekey(seqrun).set_index('id')
    # find the name of the dataset
    dset = datasets.loc[lambda x: x['seq_run'] == seqrun, 'dataset'][0]
    ids = [int(i[1].id) for i in listfiles("{rawdir}/{dset}/{seqrun}/rawdata/{{otherdirs}}/IonXpress_{{id}}_rawlib.basecaller.bam"
                                .format(rawdir = config['rawdir'],
                                        dset = dset,
                                        seqrun = seqrun))]
    wells = platekey.loc[ids, 'well'].array
    return [("{trimdir}/{seqrun}_{plate}/{seqrun}_{plate}-{well}.trim.fastq.gz"
             .format(trimdir = config['trimdir'],
                     seqrun = seqrun,
                     plate = plate,
                     well = well))
            for well in wells]

#### Illumina demultiplexing
# look up the primers/barcodes file based on the plate ID.
def find_barcode_illumina(wildcards):
    return {
        'barcode_R1': "{tagdir}/{seq_run}_R1.fasta"
                       .format(tagdir = config['tagdir'],
                               seq_run = wildcards.seq_run),
        'barcode_R2': "{tagdir}/{seq_run}_R2.fasta"
                       .format(tagdir = config['tagdir'],
                               seq_run = wildcards.seq_run)
    }

def find_fastq_illumina(wildcards):
    return {
        'fastq_R1': glob(
            "{rawdir}/**/{seq_run}/**/Sample_{seq_run}-OT{plate}/{seq_run}-OT{plate}_S?_L???_R1_???.fastq.gz"
                .format(
                    rawdir = config['rawdir'],
                    seq_run = wildcards.seq_run,
                    plate = wildcards.plate
                )
        ),
        'fastq_R2': glob(
            "{rawdir}/**/{seq_run}/**/Sample_{seq_run}-OT{plate}/{seq_run}-OT{plate}_S?_L???_R2_???.fastq.gz"
                .format(
                    rawdir = config['rawdir'],
                    seq_run = wildcards.seq_run,
                    plate = wildcards.plate
                )
        )
    }

checkpoint demux_illumina:
    input:
        unpack(find_fastq_illumina),
        unpack(find_barcode_illumina)
    output:
        directory("{trimdir}/{{seq_run}}_00{{plate}}".format_map(config))
    params:
        fpattern_R1 = lambda wildcards: ("{trimdir}/{seq_run}_00{plate}/{seq_run}_00{plate}-{{name}}_R1f.trim.fastq.gz"
                                      .format(trimdir = config['trimdir'],
                                              seq_run = wildcards.seq_run,
                                              plate = wildcards.plate)),
        rpattern_R1 = lambda wildcards: ("{trimdir}/{seq_run}_00{plate}/{seq_run}_00{plate}-{{name}}_R1r.trim.fastq.gz"
                                      .format(trimdir = config['trimdir'],
                                              seq_run = wildcards.seq_run,
                                              plate = wildcards.plate)),
        fpattern_R2 = lambda wildcards: ("{trimdir}/{seq_run}_00{plate}/{seq_run}_00{plate}-{{name}}_R2f.trim.fastq.gz"
                                      .format(trimdir = config['trimdir'],
                                              seq_run = wildcards.seq_run,
                                              plate = wildcards.plate)),
        rpattern_R2 = lambda wildcards: ("{trimdir}/{seq_run}_00{plate}/{seq_run}_00{plate}-{{name}}_R2r.trim.fastq.gz"
                                      .format(trimdir = config['trimdir'],
                                              seq_run = wildcards.seq_run,
                                              plate = wildcards.plate)),
        untrim_R1 = lambda wildcards: ("{seq_run}_{plate}_untrim_R1.fastq.gz"
                                        .format(seq_run = wildcards.seq_run,
                                                plate = wildcards.plate)),
        untrim_R2 = lambda wildcards: ("{seq_run}_{plate}_untrim_R2.fastq.gz"
                                        .format(seq_run = wildcards.seq_run,
                                                plate = wildcards.plate))
                                              
    resources:
        walltime=5
    conda: "config/conda/demultiplex.yaml"
    group: "illumina"
    shadow: "shallow"
    log: "{logdir}/illumina_demux_{{seq_run}}_{{plate}}.log".format_map(config)
    shell:
         """
         mkdir -p {output} &&
         cutadapt \\
           --pair-adapters\\
           -a file:{input.barcode_R1}\\
           -A file:{input.barcode_R2}\\
           -m 1\\
           -o {params.fpattern_R1}\\
           -p {params.fpattern_R2}\\
           --untrimmed-output {params.untrim_R1}\\
           --untrimmed-paired-output {params.untrim_R2}\\
           {input.fastq_R1} \\
           {input.fastq_R2} >>{log} 2>&1 &&
           cutadapt \\
           --pair-adapters\\
           -a file:{input.barcode_R1}\\
           -A file:{input.barcode_R2}\\
           -m 1\\
           --trimmed-only\\
           -o {params.rpattern_R2}\\
           -p {params.rpattern_R1}\\
           {params.untrim_R2} \\
           {params.untrim_R1} >>{log} 2>&1
         """

# Return a closure which calls checkpoints.pacbio_demux.get() to indicate to Snakemake that this rule is
# dynamically calculated after completion of pacbio_demux, and then find all the demultiplexed pacbio files
# for the given plate.
def illumina_find(seq_run, plate):
    def subfind(wildcards):
        checkpoints.demux_illumina.get(seq_run = seq_run, plate = plate)
        return glob("{trimdir}/{seq_run}_00{plate}/{seq_run}_00{plate}-*.trim.fastq.gz"
                    .format(trimdir = config['trimdir'],
                            seq_run = seq_run,
                            plate = plate))
    return subfind

#### Reference databases ####
  
# download taxonomy translation files
localrules: taxon_references
rule taxon_references:
  output: "{ref_root}/{{reference}}.fasta.gz".format_map(config)
  params:
    url = "{reannotate_url}/{{reference}}.fasta.gz".format_map(config)
#  log: "{logdir}/references_{{reference}}.log".format_map(config)
  shell: "cd {config[ref_root]} && wget {params.url}"

#### Drake pipeline ####
# The R-heavy parts of the analysis are organized using the Drake package in R.
# It is very nice for handling dependencies within R.  However, it lacks the capability
# for grouping jobs into SLURM calls, so it ends up having a lot of dead CPU time (or unnecessary usage)
# on the cluster when there are heterogeneous jobs including many short jobs.
# therefore we cut the workflow up into chunks with simple dependency relations and dispatch them to SLURM
# from Snakemake.

# Hash all of the demultiplexed files so Drake will know if they have changed.
localrules: hash_demux
rule hash_demux:
  output:
    config['demux_file']
  input:
    pacbio_find('pb_500_001'),
    pacbio_find('pb_500_002'),
    pacbio_find('pb_483_001'),
    pacbio_find('pb_483_002'),
    ion_find('is_057', '001'),
    illumina_find('SH-2257', '1'),
    illumina_find('SH-2257', '2')
  threads: maxthreads
  resources:
    walltime = 60
  shell:
    """
cat << ENDENDEND | xargs -P {threads} -n 10 md5sum >{output}
{input}
ENDENDEND
    """

# Configure the drake plan
# This does not build any targets, but it does most of the preliminary work,
# including calculating which targets are outdated.
localrules: drake_plan
checkpoint drake_plan:
    output:
        drakedata         = config['drakedata']
    input:
        "{rdir}/dada.R".format_map(config),
        "{rdir}/qstats.R".format_map(config),
        "{rdir}/inferrnal.R".format_map(config),
        "{rdir}/locarrna.R".format_map(config),
        "{rdir}/lsux.R".format_map(config),
        "{rdir}/mantel.R".format_map(config),
        "{rdir}/parallel_helpers.R".format_map(config),
        "{rdir}/plate_check.R".format_map(config),
        "{rdir}/raxml.R".format_map(config),
        "{rdir}/taxonomy.R".format_map(config),
        "{rdir}/utils.R".format_map(config),
        "{rdir}/epa-ng.R".format_map(config),
        "{rdir}/epa_iterate.R".format_map(config),
        "{rdir}/mafft.R".format_map(config),
        demux = rules.hash_demux.output,
        dataset  = config['dataset'],
        regions  = config['regions'],
        methods  = config['methods'],
        cm_5_8S  = config['cm_5_8S'],
        cm_32S   = config['cm_32S'],
        platemap = config['platemap'],
        script = "{rdir}/drake_plan_heavy.R".format_map(config)
    conda: "config/conda/drake.yaml"
    threads: 1
    resources:
        walltime=1440
    log: "{logdir}/drakeplan.log".format_map(config)
    script: "{rdir}/drake_plan_heavy.R".format_map(config)


# Dereplicate the sequences in the fastq.gz files and split them into equal-size groups for region detection using ITSx
rule preITSx:
    output:
        flag = ".preITSx"
    input:
        drakedata = rules.drake_plan.output.drakedata,
        script = "{rdir}/drake-01-preITSx.R".format_map(config)
    conda: "config/conda/drake.yaml"
    threads: maxthreads
    resources:
        walltime=60
    log: "{logdir}/preITSx.log".format_map(config)
    script: "{rdir}/drake-01-preITSx.R".format_map(config)

# Detect rDNA regions using ITSx
# This script will spawn additional jobs.
localrules: ITSx
rule ITSx:
    output:
        flag = ".ITSx"
    input:
        drakedata = rules.drake_plan.output.drakedata,
        flag = rules.preITSx.output.flag,
        script = "{rdir}/drake-02-ITSx.R".format_map(config)
    conda: "config/conda/drake.yaml"
    threads: maxthreads
    resources:
        walltime=360
    log: "{logdir}/ITSx.log".format_map(config)
    script: "{rdir}/drake-02-ITSx.R".format_map(config)

# Recombine the ITSx results, split the fastq files into different regions, and do quality filtering.
rule preDADA:
    output:
        flag = ".preDADA"
    input:
        drakedata = rules.drake_plan.output.drakedata,
        flag      = rules.ITSx.output.flag,
        script    = "{rdir}/drake-03-preDADA.R".format_map(config)
    conda: "config/conda/drake.yaml"
    threads: maxthreads
    resources:
        walltime=120
    log: "{logdir}/preDADA.log".format_map(config)
    script: "{rdir}/drake-03-preDADA.R".format_map(config)

# Fit error model for each sequencing run
rule errmodel:
    output:
        flag = ".errmodel"
    input:
        drakedata = rules.drake_plan.output.drakedata,
        flag   = rules.preDADA.output.flag,
        script    = "{rdir}/drake-04-errmodel.R".format_map(config)
    conda: "config/conda/drake.yaml"
    threads: maxthreads
    resources:
        walltime = 120
    log: "{logdir}/errmodel.log".format_map(config)
    script: "{rdir}/drake-04-errmodel.R".format_map(config)

# Dereplicate, denoise, and remove chimeras for each region/plate combination
rule DADA:
    output:
        flag = ".DADA"
    input:
        drakedata = rules.drake_plan.output.drakedata,
        flag   = rules.errmodel.output.flag,
        script    = "{rdir}/drake-05-DADA.R".format_map(config)
    conda: "config/conda/drake.yaml"
    threads: maxthreads
    resources:
        walltime = 120
    log: "{logdir}/DADA.log".format_map(config)
    script: "{rdir}/drake-05-DADA.R".format_map(config)


regions_table = datasets.assign(regions=datasets.regions.str.split(',')).explode('regions')

# Function to calculate which DADA results represent each region.
def region_inputs(wildcards):
    checkpoints.drake_plan.get()
    region_meta = regions_table.loc[regions_table.regions == wildcards.region]
    return expand('.nochim_{seq_run}_{region}', zip, seq_run = region_meta.seq_run, region = region_meta.regions)

# combine the dada results for each region
localrules: preconsensus
rule preconsensus:
    output:
        flag = ".preconsensus"#,
        #bigfasta = "{datadir}/clusters/{{region}}.fasta.gz".format_map(config)
    input:
        drakedata = rules.drake_plan.output.drakedata,
        flag      = rules.DADA.output.flag,
        script = "{rdir}/drake-06-preconsensus.R".format_map(config)
    conda: "config/conda/drake.yaml"
    threads: maxthreads
    resources:
        walltime = 60
    log: "{logdir}/preconsensus.log".format_map(config)
    script: "{rdir}/drake-06-preconsensus.R".format_map(config)

localrules: cluster
rule cluster:
    output:
        biom     = "{clusterdir}/{{region}}.biom".format_map(config),
        uc       = "{clusterdir}/{{region}}.uc".format_map(config),
        otutable = "{clusterdir}/{{region}}.table".format_map(config)
    input:
        "{clusterdir}/{{region}}.fasta.gz".format_map(config)
    conda: "config/conda/vsearch.yaml"
    threads: 8
    resources:
        walltime = 60
    log: "{logdir}/cluster_{{region}}.log".format_map(config)
    shell:
      """
      vsearch --cluster_smallmem {input}\
              --usersort\
              --biomout {output.biom}\
              --sizein\
              --id 0.97\
              --uc {output.uc}\
              --otutabout {output.otutable}\
              --threads {threads}\
              --log {log}
      """

localrules: tech_compare
rule tech_compare:
  output: "{outdir}/tech_compare.pdf".format_map(config)
  input:
      rmd   = "{rmddir}/tech_compare.Rmd".format_map(config),
      clust = "{clusterdir}/ITS2.table".format_map(config)
  conda: "config/conda/drake.yaml"
  threads: 1
  resources:
      walltime = 5
  log: "{logdir}/tech_compare.log".format_map(config)
  shell:
    """
      mkdir -p {config[outdir]} &&
      R -e 'getwd(); outfile <- file.path(getwd(), "{output}"); print(outfile); rmarkdown::render("{input.rmd}", output_file = outfile)' &>{log}
    """

consensus_table = regions_table.loc[regions_table.seq_run == "pb_500"]

# Calculate a consensus LSU and long amplicon (ITS1--LR5) for each ITS2 ASV,
# and build a guide tree based on LSU
rule consensus:
    output:
        flag    = ".consensus"
    input:
        expand("{ref_root}/{db}.{region}.{method}.fasta.gz",
               ref_root = config['ref_root'],
               db = ['warcup', 'unite'],
               region = ['ITS'],
               method = ['sintax', 'dada2']),
        expand("{ref_root}/{db}.{region}.{method}.fasta.gz",
               ref_root = config['ref_root'],
               db = ['rdp_train'],
               region = ['LSU'],
               method = ['sintax', 'dada2']),
        flag = rules.preconsensus.output.flag,
        drakedata = rules.drake_plan.output.drakedata,
        script = "{rdir}/drake-07-consensus.R".format_map(config)
    conda: "config/conda/drake.yaml"
    threads: maxthreads
    resources:
        walltime = 240
    log: "{logdir}/consensus.log".format_map(config)
    script: "{rdir}/drake-07-consensus.R".format_map(config)

rule raxml:
    output:
        flag = ".raxml"
    input:
        flag = rules.consensus.output.flag,
        drakedata = "{plandir}/drake.Rdata".format_map(config),
        script = "{rdir}/drake-09-raxml.R".format_map(config)
    log: "{logdir}/raxml.log".format_map(config)
    conda: "config/conda/drake.yaml"
    threads: maxthreads
    resources:
        walltime=60*96
    script: "{rdir}/drake-09-raxml.R".format_map(config)


# do all remaining actions in the drake plan:  at the moment, this means making reports.
localrules: heavy
rule heavy:
    output: flag = touch(".heavy")
    input:
        drakedata = rules.drake_plan.output.drakedata,
        flag      = rules.raxml.output.flag,
        script    = "{rdir}/drake-10-finish.R".format_map(config)
    conda: "config/conda/drake.yaml"
    threads: 1
    resources:
        walltime = 60
    log: "{logdir}/drake_finish.log".format_map(config)
    script: "{rdir}/drake-10-finish.R".format_map(config)

localrules: lightplan
rule lightplan:
    output:
        "{rmddir}/transect_paper.pdf".format_map(config),
        "{rmddir}/transect_supplement.pdf".format_map(config),
        "{outdir}/supp_file_1.tsv".format_map(config),
        "{outdir}/supp_file_2.tsv".format_map(config),
        "{outdir}/supp_file_3.pdf".format_map(config)
    input:
        "{rdir}/dada.R".format_map(config),
        "{rdir}/mantel.R".format_map(config),
        "{rdir}/variogram.R".format_map(config),
        "{rdir}/qstats.R".format_map(config),
        "{rdir}/taxonomy.R".format_map(config),
        "{rdir}/plate_check.R".format_map(config),
        "{rdir}/output_functions.R".format_map(config),
        "{clusterdir}/ITS2.table".format_map(config),
        "{clusterdir}/short.table".format_map(config),
        "{rmddir}/preamble.tex".format_map(config),
        "{rmddir}/preamble_supplement.tex".format_map(config),
        "{rmddir}/transect_paper.Rmd".format_map(config),
        "{rmddir}/transect_supplement.Rmd".format_map(config),
        "{rmddir}/all.bib".format_map(config),
        drakedata = rules.drake_plan.output.drakedata,
        flag      = rules.heavy.output.flag,
        script    = "{rdir}/drake_plan_light.R".format_map(config)
    conda: "config/conda/drake.yaml"
    threads: 4
    resources:
        walltime = 60
    log: "{logdir}/drake_light.log".format_map(config)
    script: "{rdir}/drake_plan_light.R".format_map(config)
