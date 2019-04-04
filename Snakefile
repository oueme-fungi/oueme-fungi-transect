import pandas as pd
import os.path
from glob import glob
import re
import subprocess
from snakemake.utils import listfiles

# For testing, parse the yaml file (this is automatically done by Snakemake)
#import yaml
#with open("config/config.yaml", 'r') as ymlfile: config = yaml.load(ymlfile)

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
config['pastadir']  = "{datadir}/pasta".format_map(config)
config['plandir']   = "{datadir}/plan".format_map(config)
config['tagdir']    = "{labdir}/tags".format_map(config)

config['dataset']    = '{labdir}/datasets.csv'.format_map(config)
config['regions']    = '{labdir}/regions.csv'.format_map(config)
config['platemap']   = '{labdir}/Brendan_soil2.xlsx'.format_map(config)
config['gits7_tags'] = '{labdir}/Hectors_tag_primer_plates.xlsx'.format_map(config)
config['lr5_tags']   = '{labdir}/Brendan_soil2.xlsx'.format_map(config)

config['rdp_file']   = '{ref_root}/rdp.fasta.gz'.format_map(config)
config['silva_file']   = '{ref_root}/silva.fasta.gz'.format_map(config)
config['silva_tax_file']   = '{ref_root}/silva_tax.txt'.format_map(config)
config['unite_file']   = '{ref_root}/unite.fasta.gz'.format_map(config)
config['unite_patch_file']   = '{ref_root}/unite_patch.csv'.format_map(config)

# Find the maximum number of cores available to a single node on SLURM,
# or if we aren't in a SLURM environment, how many we have on the local machine.
try:
    maxthreads = max([int(x) for x in re.findall(r'\d+', subprocess.check_output(["sinfo", "-O", "cpus"]).decode())])
except FileNotFoundError:
    maxthreads = int(subprocess.check_output("nproc").decode())

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

# endpoint target
localrules: all
rule all:
    input: ".drake_finish"

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
        expand("{moviedir}/{{movie}}.{movietype}.bam",
               moviedir = config['moviedir'],
               movietype = ["subreads", "scraps"])
    input:
        lambda wildcards: glob("{rawdir}/**/rawdata/**/Analysis_Results/{wildcards.movie}.*.h5"
                               .format(wildcards = wildcards,
                                       rawdir = config['rawdir']),
                               recursive = True)
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
          "{ccsdir}/{{movie}}.ccs.bam".format_map(config)
    input:
         "{moviedir}/{{movie}}.subreads.bam".format_map(config)
    resources:
        walltime=120
    conda: "config/conda/pacbio.yaml"
    group: "pacbio"
    threads: 4
    log: "{logdir}/ccs_{{movie}}.log".format_map(config)
    shell:
         "ccs --polish --numThreads={threads} {input} &>{log}"

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
    group: "pacbio"
    log: "{logdir}/ccs2fastq_{{seqplate}}.log".format_map(config)
    resources:
        walltime=10
    shell:
         "bam2fastq -o {output} {input} &>{log}"

# endpoint target: generate all PacBio fastq files
rule pacbio_fastq:
    input:
        expand("{fastqdir}/{seqplate}.fastq.gz",
               fastqdir = config['fastqdir'],
               seqplate = [sp for sp in seqplates if sp.startswith('pb')])

# generate a fasta file of primer sequences for each platein the format required by cutadapt
localrules: tagfiles
rule tagfiles:
    input:
        gits7_tags = config['gits7_tags'],
        lr5_tags = config['lr5_tags'],
        dataset = config['dataset']
    output:
        expand("{tagdir}/{seqrun}.fasta", tagdir = config['tagdir'], seqrun = datasets['Seq.Run'])
    resources:
        walltime=5
    conda: "config/conda/drake.yaml"
    log: "{logdir}/tagfiles.log".format_map(config)
    script:
        "{config[rdir]}/tags.extract.R"

# look up the primers/barcodes file based on the plate ID.
def find_barcode(wildcards):
    return ("{tagdir}/{seqrun}.fasta"
            .format(tagdir = config['tagdir'],
                    seqrun = datasets2.loc[wildcards['seqplate'], 'Seq.Run']))

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
def demux_find(seqplate):
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
         demux_find('pb_500_001'),
         demux_find('pb_500_002'),
         demux_find('pb_483_001'),
         demux_find('pb_483_002')


#### IonTorrent files ####
# The IonTorrent files have already been demultiplexed by the IonTorrent pipeline.
# They only need to be renamed, converted to .fastq.gz, and have the primers trimmed away.

# look up the key to convert between barcode numbers and plate wells for an iontorrent run
def ion_platekey(seqrun):
    return (pd.read_csv(os.path.join(config['labdir'], datasets.loc[seqrun].PlateKey))
            .assign(id = lambda x: pd.to_numeric(x['tag.fwd'].str.rsplit('-').str.get(1))))

# Find the demultiplexed .bam file (named by barcode index) which corresponds to a particular sample (named by well)
def find_ion_bam(wildcards):
    # load the relevant platekey
    platekey = ion_platekey(wildcards.seqrun)
    # find the id which corresponds to the desired
    num = platekey.loc[lambda x: x['well'] == wildcards.well, 'id'].get_values()[0]
    # find the name of the dataset
    dset = datasets.loc[lambda x: x['Seq.Run'] == wildcards.seqrun, 'Dataset'][0]
    return glob("{rawdir}/{dset}/{seqrun}/rawdata/*/IonXpress_{num}*.bam".format(rawdir = config['rawdir'],
                                                                                 dset = dset,
                                                                                 seqrun = wildcards.seqrun,
                                                                                 num = '{:03d}'.format(int(num))))

localrules: bam2fastq, ion_trim
# convert a demultiplexed IonTorrent .bam file to .fastq.gz
rule bam2fastq:
    output: "{demuxdir}/{{seqrun}}_{{plate,\d+}}-{{well,[A-H]\d+}}.demux.fastq.gz".format_map(config)
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
        cutadapt --trimmed-only\
                 -m 1\
                     -g file:{input.barcode}\
                     -j {threads}\
                     -o {output}\
                     {input.fastq}\
                     &> {log}
        """

# Generate the names of all the trimmed IonTorrent files which will be generated from the available .bam files
def ion_find(seqrun, plate):
    platekey = ion_platekey(seqrun).set_index('id')
    # find the name of the dataset
    dset = datasets.loc[lambda x: x['Seq.Run'] == seqrun, 'Dataset'][0]
    ids = [int(i[1].id) for i in listfiles("{rawdir}/{dset}/{seqrun}/rawdata/{{otherdirs}}/IonXpress_{{id}}_rawlib.basecaller.bam"
                                .format(rawdir = config['rawdir'],
                                        dset = dset,
                                        seqrun = seqrun))]
    wells = platekey.loc[ids, 'well'].get_values()
    return [("{trimdir}/{seqrun}_{plate}/{seqrun}_{plate}-{well}.trim.fastq.gz"
             .format(trimdir = config['trimdir'],
                     seqrun = seqrun,
                     plate = plate,
                     well = well))
            for well in wells]

#### Reference databases ####

# Process an ITS database (e.g., UNITE)
# Split into ITS1, ITS2 using ITSx.

rule itsx_reference:
    output:
        ITS2 = "{ref_root}/{{dbname}}.ITS2.fasta.gz".format_map(config),
        ITS1 = "{ref_root}/{{dbname}}.ITS1.fasta.gz".format_map(config)
    input: "{ref_root}/{{dbname}}.fasta.gz".format_map(config)
    threads: maxthreads
    shadow: "shallow"
    params:
        shards = lambda wildcards, threads: max([threads // 4, 1]),
        cpu_per_shard = lambda wildcards, threads: max([threads // max([threads // 4, 1]), 1])
    conda: "config/conda/itsx.yaml"
    log: "{logdir}/{{dbname}}_ITSx.log".format_map(config)
    shell:
        """
        #if [ {params.shards} -eq 1 ] ; then
        #    zcat {input} |
        #    ITSx --complement F\
        #         --cpu {threads}\
        #         --summary F\
        #         --graphical F\
        #         --preserve T\
        #         --positions F\
        #         --not_found F\
        #         --stdin T\
        #         -o {config[ref_root]}/{wildcards.dbname} &&
        #    gzip {config[ref_root]}/{wildcards.dbname}.full.fasta &&
        #    gzip {config[ref_root]}/{wildcards.dbname}.ITS2.fasta &&
        #    gzip {config[ref_root]}/{wildcards.dbname}.ITS1.fasta &&
        #    mv -f {config[ref_root]}/{wildcards.dbname}.full.fasta.gz {output.ITS} >{log}
        #    exit $?
        #else
        (   zcat {input} >temp.fasta &&
            fasta-splitter --n-parts {params.shards}\
                           temp.fasta &&
            for n in {{1..{params.shards}}}; do
             echo -o temp.part-$n -i temp.part-$n.fasta
            done |
                xargs -n 4 -P {params.shards} \
                    ITSx --complement F\
                 --cpu {threads}\
                 --summary F\
                 --graphical F\
                 --preserve T\
                 --positions F\
                 --not-found F\
                 --fasta F &&
            cat temp.part-{{1..{params.shards}}}.ITS1.fasta | gzip > {output.ITS1} &&
            cat temp.part-{{1..{params.shards}}}.ITS2.fasta | gzip > {output.ITS2} ) &> {log}
        #fi
        """

# Assume the entire database is ITS.
rule its_reference:
    output: "{ref_root}/{{dbname}}.ITS.fasta.gz".format_map(config)
    input: "{ref_root}/{{dbname}}.fasta.gz".format_map(config)
    threads: 1
    log: "{logdir}/{{dbname}}_ITS.log".format_map(config)
    shell: "rm -f {output} && ln -s {input} {output}"

# Process an LSU reference database (e.g., RDP or Silva)
# 1- Take only Eukaryote sequences (Silva)
#    RDP is all Fungi, so everything in it should pass this step.
# 2- Select only sequences with the conserved LR5 primer site, and trim everything in the 3' direction
# 3- Select only sequences with the conserved ITS4 primer site, but don't trim the 5' end.
#    This should select about 900bp of the 5' end of LSU, including the D1(ES5), D2(ES7), and D3(ES9) variable regions.
#
# Ambiguities in the sequences are from the RNA data project and Tedersoo et al 2015
localrules: lsu_reference
rule lsu_reference:
    output: "{ref_root}/{{dbname}}.LSU.fasta.gz".format_map(config)
    input: "{ref_root}/{{dbname}}.fasta.gz".format_map(config)
    threads: 2
    conda: "config/conda/demultiplex.yaml"
    log: "{logdir}/{{dbname}}_LSU.log".format_map(config)
    shell:
        """
        zcat {input} |
        awk '/^>/ {{ p = ($0 ~ /(Eukaryota|Fungi)/)}} p' |
        cutadapt -a CGAAnUUUCnCUCAGGAUAGCnG --trimmed-only --report minimal -e 0.2 -m 50 - |
        cutadapt -g CAUAUnAnUnAGnSSAGGA --trimmed-only --report minimal -e 0.2 - |
        gzip - >{output}
        """

#### Drake pipeline ####
# The R-heavy parts of the analysis are organized using the Drake package in R.
# It is very nice for handling dependencies within R.  However, it lacks the capability
# for grouping jobs into SLURM calls, so it ends up having a lot of dead CPU time (or unnecessary usage)
# on the cluster when there are heterogeneous jobs including many short jobs.
# therefore we cut the workflow up into chunks with simple dependency relations and dispatch them to SLURM
# from Snakemake.

# Configure the drake plan
# This does not build any targets, but it does most of the preliminary work, including calculating which targets
# are outdated.
localrules: drake_plan
checkpoint drake_plan:
    output:
        plan        = "{plandir}/plan.rds".format_map(config),
        itsx_meta   = "{plandir}/itsx_meta.rds".format_map(config),
        predada_meta = "{plandir}/predada_meta.rds".format_map(config),
        dada_meta       = "{plandir}/dada_meta.rds".format_map(config),
        taxonomy_meta       = "{plandir}/taxonomy_meta.rds".format_map(config),
        drakedata   = "{plandir}/drake.Rdata".format_map(config),
        taxonomy_meta_csv    = "{plandir}/taxonomy_meta.csv".format_map(config),
        tids        = "{plandir}/tids.txt".format_map(config)
    input:
        demux_find('pb_500_001'),
        demux_find('pb_500_002'),
        demux_find('pb_483_001'),
        demux_find('pb_483_002'),
        ion_find('is_057', '001'),
        "{rdir}/parallel_helpers.R".format_map(config),
        "{rdir}/combine_derep.R".format_map(config),
        "{rdir}/extract_regions.R".format_map(config),
        "{rdir}/dada.R".format_map(config),
        "{rdir}/plate_check.R".format_map(config),
        "{rdir}/map_LSU.R".format_map(config),
        "{rdir}/taxonomy.R".format_map(config),
        packrat  = "packrat/packrat.lock",
        dataset  = config['dataset'],
        regions  = config['regions'],
        platemap = config['platemap'],
        script   = "{rdir}/drake.R".format_map(config),
        rdp      = config['rdp_file'],
        unite    = config['unite_file'],
        unite_patch = config['unite_patch_file'],
        silva    = config['silva_file'],
        silva_tax= config['silva_tax_file']
    conda: "config/conda/drake.yaml"
    threads: 1
    resources:
        walltime=1440
    log: "{logdir}/drakeplan.log".format_map(config)
    script: "{input.script}"


# Dereplicate the sequences in the fastq.gz files and split them into equal-size groups for region detection using ITSx
rule preITSx:
    output:
        flag = touch(".preITSx")
    input:
        drakedata = rules.drake_plan.output.drakedata,
        script    = "{rdir}/drake-01-preITSx.R".format_map(config)
    conda: "config/conda/drake.yaml"
    threads: 4
    resources:
        walltime=60
    log: "{logdir}/preITSx.log".format_map(config)
    script: "{input.script}"

# Detect rDNA regions using ITSx
# This script will spawn additional jobs.
localrules: ITSx
rule ITSx:
    output:
        flag = touch(".ITSx")
    input:
        drakedata = rules.drake_plan.output.drakedata,
        preITSx = ".preITSx",
        script = "{rdir}/drake-02-ITSx.R".format_map(config)
    conda: "config/conda/drake.yaml"
    threads: 1
    resources:
        walltime=360
    log: "{logdir}/ITSx.log".format_map(config)
    script: "{input.script}"

# Recombine the ITSx results, split the fastq files into different regions, and do quality filtering.
rule preDADA:
    output: touch(".preDADA")
    input:
        drakedata = rules.drake_plan.output.drakedata,
        ITSx      = rules.ITSx.output.flag,
        script = "{rdir}/drake-03-preDADA.R".format_map(config)
    conda: "config/conda/drake.yaml"
    threads: 4
    resources:
        walltime=120
    log: "{logdir}/preDADA.log".format_map(config)
    script: "{input.script}"

# Dereplicate, denoise, and remove chimeras for each region/plate combination
rule DADA:
    output: touch(".nochim_{RID}")
    input:
        drakedata = rules.drake_plan.output.drakedata,
        preDADA   = ".preDADA",
        script    = "{rdir}/drake-04-DADA.R".format_map(config)
    conda: "config/conda/drake.yaml"
    threads: 8
    resources:
        walltime = 120
    log: "{logdir}/DADA_{{RID}}.log".format_map(config)
    script: "{input.script}"

# Function to calculate which DADA results represent each region.
def region_inputs(wildcards):
    checkpoints.drake_plan.get()
    taxonomy_meta = pd.read_csv(rules.drake_plan.output.taxonomy_meta_csv).set_index("Region")
    PIDs = taxonomy_meta.loc[[wildcards.region], 'PlateRegionID'].unique()[0].split(',')
    return expand('.nochim_{PlateRegionID}', PlateRegionID = PIDs)

# combine the dada results for each region
localrules: region_table
rule region_table:
    output: touch(".big_seq_table_{region}")
    input:
        drakedata = rules.drake_plan.output.drakedata,
        dada      = region_inputs,
        script    = "{rdir}/drake-06-pretaxonomy.R".format_map(config)
    conda: "config/conda/drake.yaml"
    threads: 1
    resources:
        walltime = 10
    log: "{logdir}/pretaxonomy_{{region}}.log".format_map(config)
    script: "{input.script}"

# calculate which sequence tables are needed for a taxonomy assignment step
def taxon_inputs(wildcards):
    checkpoints.drake_plan.get()
    taxonomy_meta = pd.read_csv(rules.drake_plan.output.taxonomy_meta_csv).set_index("TaxID")
    regions = taxonomy_meta.loc[wildcards.TaxID, 'Region'].split(sep = ",")
    return expand('.big_seq_table_{region}', region = regions)

# call taxonomy and assign guilds
rule taxonomy:
    output: touch(".guilds_table_{TaxID}")
    input:
        drakedata = rules.drake_plan.output.drakedata,
        dada      = taxon_inputs,
        script    = "{rdir}/drake-07-taxonomy.R".format_map(config)
    conda: "config/conda/drake.yaml"
    threads: 8
    resources:
        walltime = 240
    log: "{logdir}/taxonomy_{{TaxID}}.log".format_map(config)
    script: "{input.script}"

# calculate all the taxonomy outputs which should be calculated
def taxon_outputs(wildcards):
    checkpoints.drake_plan.get()
    file = open(rules.drake_plan.output.tids, mode = 'r')
    TaxIDs = file.read().splitlines()
    file.close()
    return expand(".guilds_table_{TaxID}", TaxID = TaxIDs)

# Calculate a consensus LSU and long amplicon (ITS1--LR5) for each ITS2 ASV,
# and build a guide tree based on LSU
rule consensus:
    output:
        flag    = touch(".consensus"),
        longasv = "{pastadir}/long_ASVs.fasta".format_map(config)
    input:
        drakedata = rules.drake_plan.output.drakedata,
        taxon     = taxon_outputs,
        script    = "{rdir}/drake-08-consensus.R".format_map(config)
    conda: "config/conda/drake.yaml"
    threads: 8
    resources:
        walltime = 240
    log: "{logdir}/consensus.log".format_map(config)
    script: "{input.script}"

# Build a progressive tree using PASTA, with the LSU tree as guide
rule pasta:
    output:
        tree = "{pastadir}/pasta_raxml.tree".format_map(config)
    input:
        consensus = rules.consensus.output.flag,
        longasv   = rules.consensus.output.longasv
    log: "{logdir}/pasta.log".format_map(config)
    threads: 16
    resources:
        walltime=60
    conda: "config/conda/pasta.yaml"
    shell:
        """
        mkdir -p {config[pastadir]} &&
        cd {config[pastadir]} &&
        [ -e .dopasta ] &&
        rm -f *postraxtree_tree.tre &&
        PASTA_TOOLS_DEVDIR=$CONDA_PREFIX/bin/ \
        run_pasta.py \
        -j ITS_LSU \
        --input {input.longasv} \
        --aligned \
        --datatype rna \
        --raxml-search-after \
        --num-cpus={threads} &&
        mv {config[pastadir]}/*postraxtree_tree.tre {output.tree} &&
        rm -f {config[pastadir]}/.dopasta &> {log} ||
        ( echo "alignment unchanged, skipping PASTA and using old tree.";
          touch {output.tree} )
        """

# do all remaining actions in the drake plan:  at the moment, this means making reports.
localrules: finish
rule finish:
    output: touch(".drake_finish")
    input:
        pasta     = rules.pasta.output.tree,
        drakedata = rules.drake_plan.output.drakedata,
        taxonomy  = taxon_outputs,
        script    = "{rdir}/drake-09-finish.R".format_map(config)
    conda: "config/conda/drake.yaml"
    threads: 1
    resources:
        walltime = 60
    log: "{logdir}/drake_finish.log".format_map(config)
    script: "{input.script}"
