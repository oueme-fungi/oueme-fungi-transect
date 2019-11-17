import pandas as pd
import os.path
from glob import glob
import re
import subprocess
from snakemake.utils import listfiles
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider

HTTP = HTTPRemoteProvider()

# For testing, parse the yaml file (this is automatically done by Snakemake)
#import yaml
#with open("config/config.yaml", 'r') as ymlfile: config = yaml.load(ymlfile)

configfile: "config/config.yaml"

# single part directory names are given in config.yaml
# this puts together the ones which are composed of previous values
config['moviedir']   = "{seqdir}/rawmovie".format_map(config)
config['ccsdir']     = "{seqdir}/ccs".format_map(config)
config['fastqdir']   = "{seqdir}/rawfastq".format_map(config)
config['demuxdir']   = "{seqdir}/demux".format_map(config)
config['trimdir']    = "{seqdir}/trim".format_map(config)
config['regiondir']  = "{seqdir}/regions".format_map(config)
config['filterdir']  = "{seqdir}/filter".format_map(config)
config['clusterdir'] = "{datadir}/clusters".format_map(config)
config['locarnadir'] = "{datadir}/mlocarna".format_map(config)
config['pastadir']   = "{datadir}/pasta".format_map(config)
config['plandir']    = "{datadir}/plan".format_map(config)
config['tagdir']     = "{labdir}/tags".format_map(config)

config['dataset']    = '{labdir}/datasets.csv'.format_map(config)
config['regions']    = '{labdir}/regions.csv'.format_map(config)
config['methods']    = '{labdir}/taxonomy_methods.csv'.format_map(config)
config['platemap']   = '{labdir}/Brendan_soil2.xlsx'.format_map(config)
config['gits7_tags'] = '{labdir}/Hectors_tag_primer_plates.xlsx'.format_map(config)
config['lr5_tags']   = '{labdir}/Brendan_soil2.xlsx'.format_map(config)

config['cm_5_8S'] = '{ref_root}/RF00002.cm'.format_map(config)
config['cm_32S']  = '{ref_root}/fungi_32S_LR5.cm'.format_map(config)

config['cmaln_long']   = "{locarnadir}/long_cmalign.aln".format_map(config),
config['guide_tree']   = "{locarnadir}/32S_guide.tree".format_map(config),
config['mlocarna_pp_dir'] = "{locarnadir}/consensus_pp".format_map(config)
config['mlocarna_dir'] = "{locarnadir}/consensus".format_map(config)
config['makelocarna']  = "{rdir}/snakemakelocarna.sh".format_map(config)
config['makelocarna_profile'] = "{configdir}/snakemakelocarnaUPPMAX".format_map(config)
config['makelocarna_conda'] = "{condadir}/snakemakelocarna.yaml".format_map(config)
config['mlocarna_aln'] = "{mlocarna_dir}/results/result.stk".format_map(config)
config['raxml_locarna_dir']    = "{datadir}/raxml/locarna".format_map(config)
config['raxml_decipher_dir']    = "{datadir}/raxml/decipher".format_map(config)

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
        ".drake_finish"
        #{outdir}/tech_compare.pdf".format_map(config)

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

# generate a fasta file of primer sequences for each plate in the format required by cutadapt
localrules: tagfiles
rule tagfiles:
    input:
        gits7_tags = config['gits7_tags'],
        lr5_tags = config['lr5_tags'],
        dataset = config['dataset'],
        script = "{config[rdir]}/tags.extract.R"
    output:
        expand("{tagdir}/{seqrun}.fasta", tagdir = config['tagdir'], seqrun = datasets['seq_run'])
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
           - 2>{log} is:open 
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


# put all the demultiplexed long PacBio reads in one file, with labeled samples
localrules: pacbio_singledemux
rule pacbio_singledemux:
    input:
        demux_find('pb_500_001'),
        demux_find('pb_500_002')
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
    # load the relevant platekey is:open 
    platekey = ion_platekey(wildcards.seqrun)
    # find the id which corresponds to the desired
    num = platekey.loc[lambda x: x['well'] == wildcards.well, 'id'].get_values()[0]
    # find the name of the dataset
    dset = datasets.loc[lambda x: x['seq_run'] == wildcards.seqrun, 'dataset'][0]
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
    wells = platekey.loc[ids, 'well'].get_values()
    return [("{trimdir}/{seqrun}_{plate}/{seqrun}_{plate}-{well}.trim.fastq.gz"
             .format(trimdir = config['trimdir'],
                     seqrun = seqrun,
                     plate = plate,
                     well = well))
            for well in wells]

#### Reference databases ####

# Download the Unite database
localrules: unite_download
rule unite_download:
    output: "{ref_root}/unite.raw.fasta.gz".format_map(config)
    input: 
      zip = HTTP.remote(config['unite_url'], allow_redirects = True, keep_local = True)
    shadow: "shallow"
    shell:
        """
        echo {config[unite_md5]} {input} |
        md5sum -c - &&
        mkdir -p $(dirname {output}) &&
        unzip {input} &&
        gzip -c {config[unite_filename]} > {output}
        """

# Download the RDP fungal LSU training set
localrules: rdp_download
rule rdp_download:
    output:
        fasta = "{ref_root}/rdp_train.raw.fasta.gz".format_map(config),
        taxa  = "{ref_root}/rdp_train.taxa.txt".format_map(config)
    input:
      zip = HTTP.remote(config['rdp_url'], allow_redirects = True, keep_local = True)
    shadow: "shallow"
    shell:
        """
        echo {config[rdp_md5]} {input.zip} |
        md5sum -c - &&
        mkdir -p $(dirname {output.fasta}) &&
        unzip {input.zip} &&
        gzip -c {config[rdp_filename]} > {output.fasta} &&
        mv {config[rdp_taxa]} {output.taxa}
        """

# Download the Warcup fungal ITS training set
localrules: warcup_download
rule warcup_download:
    output:
        fasta = "{ref_root}/warcup.fasta.gz".format_map(config),
        taxa  = "{ref_root}/warcup.taxa.txt".format_map(config)
    input:
      zip = HTTP.remote(config['warcup_url'], allow_redirects = True, keep_local = True)
    shadow: "shallow"
    shell:
        """
        echo {config[warcup_md5]} {input} |
        md5sum -c - &&
        mkdir -p $(dirname {output.fasta}) &&
        unzip {input} &&
        gzip -c {config[warcup_filename]} > {output.fasta} &&
        mv {config[warcup_taxa]} {output.taxa}
        """

# Process an ITS database (e.g., UNITE)
# Split into ITS1, ITS2 using ITSx.

rule itsx_reference:
    output:
        ITS2 = "{ref_root}/{{dbname}}.ITS2.fasta.gz".format_map(config),
        ITS1 = "{ref_root}/{{dbname}}.ITS1.fasta.gz".format_map(config)
    input: "{ref_root}/{{dbname}}.fasta.gz".format_map(config)
    threads: maxthreads
    shadow: "shallow"
    resources:
        walltime = 120
    params:
        shards = lambda wildcards, threads: max([threads // 4, 1]),
        cpu_per_shard = lambda wildcards, threads: max([threads // max([threads // 4, 1]), 1])
    conda: "config/conda/itsx.yaml"
    log: "{logdir}/{{dbname}}_ITSx.log".format_map(config)
    shell:
        """
        
        (   zcat {input} >temp.fasta &&
            fasta-splitter --n-parts {params.shards}\
                           temp.fasta &&
            for n in {{1..{params.shards}}}; do
             echo -o temp.part-$n -i temp.part-$n.fasta
            done |
                xargs -n 4 -P {params.shards} \
                    ITSx --complement F\
                 --cpu {params.cpu_per_shard}\
                 --summary F\
                 --graphical F\
                 --preserve T\
                 --positions F\
                 --not-found F\
                 --fasta F
            echo "xargs exited with value $?"
            echo "combining ITS1 files:" &&
            ls -l *.ITS1.fasta &&
            cat temp.part-{{1..{params.shards}}}.ITS1.fasta | gzip - > {output.ITS1} &&
            echo "combining ITS2 files:" &&
            ls -l *.ITS2.fasta &&
            cat temp.part-{{1..{params.shards}}}.ITS2.fasta | gzip - > {output.ITS2} ) &> {log}
        """

# Assume the entire database is ITS.
# UNITE also includes part of LSU in some sequences; this will be helpful when IDing long amplicons
localrules: its_reference
rule its_reference:
    output: "{ref_root}/{{dbname}}.ITS.fasta.gz".format_map(config)
    input: "{ref_root}/{{dbname}}.fasta.gz".format_map(config)
    params:
      fullinput = lambda wildcards, input: os.path.abspath(input[0]),
      fulloutput = lambda wildcards, output: os.path.abspath(output[0])
    threads: 1
    log: "{logdir}/{{dbname}}_ITS.log".format_map(config)
    shell: "rm -f {output} && cd $(dirname {input}) && ln -s $(basename {input}) $(basename {output})"

# Process an LSU reference database (e.g., RDP or Silva)
# At the moment this doesn't do anything
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
        ln -sr {input} {output}
        """

# Eukarya classification system proposed by Tedersoo
localrules: tedersoo_classification
rule tedersoo_classification:
    output: "{ref_root}/Tedersoo_Eukarya_classification.xlsx".format_map(config)
    input: HTTP.remote(config['tedersoo_url'], allow_redirects = True)
    shadow: "shallow"
    shell:
        """
        echo {config[tedersoo_md5]} {input} |
        md5sum -c - &&
        mv {input} {output}
        """

# generate taxonomy translation files
# These are files of sed commands that will translate the headers of the fasta
# files to the formats required by each database, as well as translating all
# the taxonomy to the Tedersoo system
localrules: translate_references
rule translate_references:
  output:
        expand("{ref_root}/{db}.{region}.{method}.fasta.gz",
               ref_root = config['ref_root'],
               db = ['warcup', 'unite'],
               region = ['ITS'],
               method = ['sintax', 'dada2']),
        expand("{ref_root}/{db}.{region}.{method}.fasta.gz",
               ref_root = config['ref_root'],
               db = ['rdp_train'],
               region = ['LSU'],
               method = ['sintax', 'dada2'])
  input:
    expand("{ref_root}/{dbname}.{type}",
           ref_root = config['ref_root'],
           dbname = ['rdp_train', 'warcup', 'unite'],
           type = ['fasta.gz', 'pre.sed']),
    "{rdir}/taxonomy.R".format_map(config),
    tedersoo = rules.tedersoo_classification.output,
    regions = config['regions'],
    script = "{rdir}/make_taxonomy.R".format_map(config)
  resources:
    walltime = 60
  conda: "config/conda/drake.yaml"
  log: "{logdir}/translate_references.log".format_map(config)
  script: "{rdir}/make_taxonomy.R".format_map(config)

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
        plan              = "{plandir}/plan.rds".format_map(config),
        itsx_meta         = "{plandir}/itsx_meta.rds".format_map(config),
        predada_meta      = "{plandir}/predada_meta.rds".format_map(config),
        dada_meta         = "{plandir}/dada_meta.rds".format_map(config),
        taxonomy_meta     = "{plandir}/taxonomy_meta.rds".format_map(config),
        drakedata         = "{plandir}/drake.Rdata".format_map(config),
        taxonomy_meta_csv = "{plandir}/taxonomy_meta.csv".format_map(config),
        tids              = "{plandir}/tids.txt".format_map(config)
    input:
        demux_find('pb_500_001'),
        demux_find('pb_500_002'),
        demux_find('pb_483_001'),
        demux_find('pb_483_002'),
        ion_find('is_057', '001'),
        "{rdir}/dada.R".format_map(config),
        "{rdir}/extract_regions.R".format_map(config),
        "{rdir}/inferrnal.R".format_map(config),
        "{rdir}/locarrna.R".format_map(config),
        "{rdir}/lsux.R".format_map(config),
        "{rdir}/mantel.R".format_map(config),
        "{rdir}/parallel_helpers.R".format_map(config),
        "{rdir}/plate_check.R".format_map(config),
        "{rdir}/raxml.R".format_map(config),
        "{rdir}/taxonomy.R".format_map(config),
        "{rdir}/utils.R".format_map(config),
        dataset  = config['dataset'],
        regions  = config['regions'],
        methods  = config['methods'],
        cm_5_8S  = config['cm_5_8S'],
        cm_32S   = config['cm_32S'],
        platemap = config['platemap'],
        script = "{rdir}/drake.R".format_map(config)
    conda: "config/conda/drake.yaml"
    threads: 1
    resources:
        walltime=1440
    log: "{logdir}/drakeplan.log".format_map(config)
    script: "{rdir}/drake.R".format_map(config)


# Dereplicate the sequences in the fastq.gz files and split them into equal-size groups for region detection using ITSx
rule preITSx:
    output:
        flag = touch(".preITSx")
    input:
        drakedata = rules.drake_plan.output.drakedata,
        script = "{rdir}/drake-01-preITSx.R".format_map(config)
    conda: "config/conda/drake.yaml"
    threads: 4
    resources:
        walltime=60
    log: "{logdir}/preITSx.log".format_map(config)
    script: "{rdir}/drake-01-preITSx.R".format_map(config)

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
    threads: maxthreads
    resources:
        walltime=360
    log: "{logdir}/ITSx.log".format_map(config)
    script: "{rdir}/drake-02-ITSx.R".format_map(config)

# Recombine the ITSx results, split the fastq files into different regions, and do quality filtering.
rule preDADA:
    output: touch(".preDADA")
    input:
        drakedata = rules.drake_plan.output.drakedata,
        ITSx      = rules.ITSx.output.flag,
        script    = "{rdir}/drake-03-preDADA.R".format_map(config)
    conda: "config/conda/drake.yaml"
    threads: maxthreads
    resources:
        walltime=120
    log: "{logdir}/preDADA.log".format_map(config)
    script: "{rdir}/drake-03-preDADA.R".format_map(config)

# Dereplicate, denoise, and remove chimeras for each region/plate combination
rule DADA:
    output: touch(".DADA")
    input:
        drakedata = rules.drake_plan.output.drakedata,
        preDADA   = ".preDADA",
        script    = "{rdir}/drake-04-DADA.R".format_map(config)
    conda: "config/conda/drake.yaml"
    threads: maxthreads
    resources:
        walltime = 120
    log: "{logdir}/DADA.log".format_map(config)
    script: "{rdir}/drake-04-DADA.R".format_map(config)


regions_table = datasets.assign(regions=datasets.regions.str.split(',')).explode('regions')

# Function to calculate which DADA results represent each region.
def region_inputs(wildcards):
    checkpoints.drake_plan.get()
    region_meta = regions_table.loc[regions_table.regions == wildcards.region]
    return expand('.nochim_{seq_run}_{region}', zip, seq_run = region_meta.seq_run, region = region_meta.regions)

# combine the dada results for each region
localrules: region_table
rule region_table:
    output:
        touch(".preconsensus")#,
        #bigfasta = "{datadir}/clusters/{{region}}.fasta.gz".format_map(config)
    input:
        drakedata = rules.drake_plan.output.drakedata,
        dada      = ".DADA",
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
        flag    = touch(".consensus")
    input:
        ".DADA",
        ".preconsensus",
        rules.translate_references.output,
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
        touch(".raxml")
    input:
        ".consensus",
        makelocarna = config['makelocarna'],
        cmaln_long = config['cmaln_long'],
        guide_tree = config['guide_tree'],
        drakedata = "{plandir}/drake.Rdata".format_map(config),
        script = "{rdir}/drake-09-raxml.R".format_map(config)
    log: "{logdir}/raxml.log".format_map(config)
    conda: "config/conda/drake.yaml"
    threads: maxthreads
    resources:
        walltime=60*96
    script: "{rdir}/drake-09-raxml.R".format_map(config)


# do all remaining actions in the drake plan:  at the moment, this means making reports.
localrules: finish
rule finish:
    output: touch(".drake_finish")
    input:
#        pasta     = rules.pasta.output.tree,
        drakedata = rules.drake_plan.output.drakedata,
        raxml  = ".raxml"
    conda: "config/conda/drake.yaml"
    threads: 1
    resources:
        walltime = 60
    log: "{logdir}/drake_finish.log".format_map(config)
    script: "{rdir}/drake-09-finish.R".format_map(config)
