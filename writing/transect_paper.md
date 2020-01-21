---
title: "A Tale of Two Transects"
subtitle: "Also Two Sequencing Technologies, and Two Amplicons"
author: 
  - Brendan Furneaux^1^,
    Roël Houdanon^2^,
    Mohammad Bahram^1^,
    Anna Rosling^3^,
    Nourou S. Yorou^2^,
    Martin Ryberg^1^
  - ^1^Program in Systematic Biology,
    Department of Organismal Biology,
    Uppsala University,
    Uppsala, Sweden
  - ^2^Research Unit in Tropical Mycology and Plant-Fungi Interactions,
    LEB,
    University of Parakou,
    Parakou, Benin
  - ^3^Program in Evolutionary Biology,
    Department of Ecology and Genetics,
    Uppsala University,
    Uppsala, Sweden
date: "December 05, 2019"
bibliography:
  - "all.bib"
  - "R.bib"
csl: fungal-ecology.csl
panflute-filters: figurefilter
panflute-path: '../../python'
output:
  bookdown::word_document2: 
    keep_md: yes
  bookdown::pdf_document2:
    fig_crop: no
    keep_tex: yes
    latex_engine: lualatex
    includes:
      in_header: preamble.tex
---


# Introduction

## Introduction to ECM communities

- Ectomycorrhizas (ECM) are a symbiosis between fungi and plant roots.
- x% of plants, y% of fungi
- The plants are often ecologically and economically important trees (Pinaceae, Fagaceae, Dipterocarpaceae, Amherstiae)
- The fungi provide nutrients and protection from pathogens, the plant provides sugar.
- Many well-known mushrooms are produced by ECM (Amanita, Cantharellus, Boletus) but many ECM produce inconspicuous (Tomentella) or no (Cenococcum) fruit bodies.
- Anyway the fruitbodies are ephemeral, so study of ECM communities requires looking at vegetative structures.
- Not easy to culture, so we need molecular methods.

## Turnover scale in ECM forests

- @lilleskov2004 - autocorrelation only <2.6 m at most sites.  Based on Sanger sequencing root tips in temperate forests.

- Meta analysis @Bahram2013
  - Includes data from Benin (same site? ask Nourou/Leho) but this is unpublished
  - conclude that there is greater range of distance decay in tropical forests; because of less host availability? older forests?  *Russulaceae*?
  - 65 m in Benin

- @pickles2012a
  - T-RFLP of root tips, temperate forest
  - Autocorrelation < 3.4 m

- Data from African wooded savannas - @tedersoo2011 but only southern band
  - Irregular sampling with no samples nearer than 8m
  - based on sanger sequencing root tips
  - Significant Moran's I only in <10m class

## Long reads vs short reads/phylogenetic methods

- Clustering vs. phylogenetic methods
    - Clustering-based OTU approach may both "clump" and "split" species.
    - Distance-based clustering conflates intra-species variation and sequencing error
    - Denoising (DADA2 [@callahan2017], Deblur [@amir2017], UNOISE2 [@edgar2016]) attempts to control for sequencing error while leaving intra-species variation.
    - Phylogenetic community distance measures like UNIFRAC [@wong2016, look for earlier reference] are relatively insensitive to species/OTU delimitation, but require a phylogenetic tree.

- Short reads
    - 454 [@buee2009], Ion Torrent, Illumina [@schmidt2013].
    - For fungi, frequently ITS1 or ITS2 [@schoch2012]
    - Dependent on matching to database for taxonomic placement.
    - Phylogenetic placement algorithms [@matsen2010; @berger2011] are not easy to apply because ITS is not in most LSU alignments (different situation than V4/V5 of 18S)
    - Possible to place on taxonomy-based tree [@tedersoo2018]

- Long reads
    - Can include both LSU/SSU and ITS in same amplicon for better taxonomic placement [@Tedersoo2018]
    - can build tree directly
    - short reads from same location can be mapped to tree.
    - Sequencing depth tradeoff [@Kennedy2018]

# Methods

## Sampling

Sampling was conducted at two sites approximately 30 km apart in the _Forêt Reservée de l'Ouémé Supérieur_ (Upper Ouémé Forest Reserve) in central Benin.
Both locations were located in woodlands dominated by the ECM host tree _Isoberlinia doka_ (Caesalpinioideae).
At each site, 25 soil samples were collected at intervals of 1 m in May 2015.
Additional samples were collected from one third of the sample locations (3-m spacing) one year later in June 2016.
For each sample, any coarse organic debris was removed from the soil surface and a sample of approximately 5cm×5cm×5cm was extracted with a knife blade.
Each sample was sealed in a plastic zipper bag and homogenized by shaking and manually breaking apart soil aggregations.
Approximately 50 mg total of soil was collected from two locations in the homogenized soil and placed in Xpedition lysis solution (Zymo Research) and lysed in the field using a handheld bead-beater (Terralyzer; Zymo Research).
The remainder of each sample was dried at 40°C with an electric dehydrator within 24-48 hours of collection for bulk chemical analysis.

## Soil Chemistry

  - Need to get methods from Nourou
  - Or just not use it?

## DNA extraction, amplification, and sequencing

After field lysis, DNA was extracted using the Zymo Research Soil/Fecal Prep kit **check name**.
DNA was quantified using PicoGreen *according to the manufacturer's protocol*?

Two different regions of the nuclear ribosomal DNA were amplified (figure \@ref(fig:rDNA)).
The short amplicon (approximately 300 bp) included the full ITS2 region as well as parts of the flanking 5.8S and large subunit (LSU) ribosomal RNA, using gITS7 [@ihrmark2012] as the forward primer and a mix of ITS4 [@white1990amplification] and ITS4A [@larena1999] (hereafter ITS4m) as the reverse primer.
The long amplicon (approximately 1550 bp) included approximately 12 bp at the 3' end of the ribosomal small subunit (SSU) RNA, the full ITS region including the 5.8S rRNA, and approximately 950 bp at the 5' end of the LSU, including the first three variable domains.

The gITS7 primers for the short amplicon were indexed for multiplexing (**supp info**).
Amplification was performed by polymerase chain reaction (PCR) in 20µl reactions containing 200 µM dNTP mix, 250 µM indexed gITS7 primer, 150 µM ITS4m, 2 mM McCl2, 0.1 U *Taq* polymerase (Dream *Taq*, Thermo Fisher) and 3--7 ng purified DNA in Dream *Taq* buffer.
The reaction conditions were 10 min at 95°, followed by 35 cycles of 60 s at 95°, 45 s at 56°, and 50 s at 72°, and finally 3 min at 72°.
Each reaction was conducted in three technical replicates to reduce the effect of PCR stochasticity (**ref**), which were pooled after amplification.

![(\#fig:rDNA)rDNA amplicons](transect_paper_files/figure-docx/rDNA-1.pdf)

Both primers for the long amplicon were indexed for multiplexing (**supp info**).
PCR was performed as for the short amplicons, but with 500 µM of each of the two primers.
Reaction conditions were 10 min at 95°, 30 cycles of 45 s at 95°, 45 s at 59°, and 90 s at 72°, and finally 10 min at 72°.
Each reaction was performed in three technical replicates as for short amplicons.

Amplicons were purified using SPRI beads [@vesterinen2016] and quantified fluorometrically using PicoGreen dsDNA (Quant-iT) fluorescent indicator dye on a Infinite F200 plate spectrofluorometer (Tecan) according to the manufacturer's protocol.
100 ng of DNA from each sample (or the total PCR product if less than 100 ng) was pooled into two libraries each for long and short amplicons.
Each library was sequenced using Single Molecule Real Time (SMRT) sequencing on a Pacific Biosciences RS II sequencer at the Uppsala Genome Center (UGS) at SciLifeLab in Uppsala, Sweden.
Short amplicon libraries were sequenced on two SMRT cells each, while long amplicon libraries were sequenced on four SMRT cells each.

Additionally the short amplicon libraries were combined and sequenced using an Ion S5  (Ion Torrent) sequencer using one 520 chip at UGS.
  
## Bioinformatics

  - Two workflows for short reads:
    - Mohammad's workflow in Pipecraft
    - Brendan's workflow using Snakemake/Drake
  
Circular consensus sequence (CCS) basecalls for SMRT sequences were made using `ccs` version 3.4 [@pacificbiosciences2019] using the default settings.
The resulting sequences were demultiplexed and sequencing primers were removed using `cutadapt` version 1.18 [@martin2011].
Sequencing primers were similarly removed from the Ion Torrent sequences, but interference between the tagged gITS7 primers and the Ion XPress tags used in library prep made full demultiplexing of the Ion Torrent sequences impossible.

Raw reads were divided into regions/domains by matching to covariance models (CM), which are similar to stochastic hidden markov models (HMM), but account for both nucleotide sequence and RNA secondary structure [@eddy1994].
First, the 5.8S rRNA was located in each read by searching for Rfam model RF0002 [@kalvari2018] using `cmsearch` from Infernal 1.1.2 [@nawrocki2013].
This approach was able to identify both the 5' and 3' ends of the 5.8S rRNA more consistently and quickly than ITSx [@bengtsson-palme2013, data not shown].
Each read was then aligned to a CM based on a modification of the fungal 28S RNA seed alignment from the Ribosomal Data Project (RDP) release 11.5 [@glockner2017; @cole2014] using `cmalign` from Infernal 1.1.2 [@nawrocki2013].
The reference line in the CM alignment for each read were then used to split the reads into alternating more-conserved and less-conserved regionss:
ITS1, 5.8S, ITS2, LSU1, D1, LSU2, D2, LSU3, D3, LSU4,
where the four "LSU" regions represent the conserved regions of the 28S rRNA flanking the "D" regions.
For short amplicons, only 5.8S, ITS2, and LSU1 were extracted.
No attempt was made to remove the approximately 12 bp fragment of the SSU from the 5' end of ITS1 in the long amplicons; it was too short to be reliably detected by a CM or the HMMs employed by ITSx.
The process of extracting the regions is automated in the new R package `LSUx`, available on github, CRAN, and bioconda (**but not yet**).

Each of the extracted regions was independently filtered for length (table **X**) and a maximum of three expected errors.
Sequences were then dereplicated and denoised into amplicon sequencing variants (ASVs) using DADA2 version 1.12.1.
The error model for DADA2 denoising was fit using the 5.8S RNA region for long amplicons, and using the entire read for short amplicons.
Independent error models were fit for each sequencing run (i.e., long *vs.* short amplicons, different sequencing technologies).
Chimeras within each region were removed using `removeBimeraDenovoTable` from DADA2.

For each ITS2 ASV from the long amplicon data set, the denoised sequences for the other regions corresponding to the same sequencing reads were concatenated to form a set of full-length reads.
For reads which were not assigned a denoised sequence for every region, the raw read for the region was used instead.
Because ITS2 is one of the most variable regions in the amplified regions, reads with identical ITS2 regions are expected have extremely similar sequences in the other regions, unless they are chimeric.
Thus, the concatenated regions for each ASV were aligned in R using the DECIPHER package [@wright2015], and outlier sequences (chimeras) were removed from each alignment using the odseq package [@jehl2015].
The consensus of the remaining aligned sequences was assigned as the full-length ASV sequence.
Full-length ASV sequences with more than three ambiguous bases (i.e., no nucleotide >50% at a given position) were removed.
A similar process was used to generate a consensus ITS (ITS1--5.8S--ITS2) and LSU (LSU1--D1--LSU2--D2--LSU3--D3--LSU4) sequence for each ASV.
The process of assigning consensus full-length ASVs was carried out using the new `tzara` package for R, available on github (https://github.com/brendanf/tzara), CRAN (not yet), and bioconda (not yet).

### Taxonomy assignment

Taxonomic annotations of the Ribosomal Data Project's LSU fungal training set (RDP) version 11.5 [@cole2014] and Warcup ITS training set [@wang2007] were mapped to the taxonomic classification system used in the Unite database version 8 [@nilsson2019a].
In particular, the classification for fungi was according to @tedersoo2018, and for non-fungal eukaryotes was according to the proposed system of @tedersoo2017c.
Although the latter system is not formally published, it seems to be in use for non-fungal eukaryotes in the Unite database.
Additionally, it is a system with only purportedly monophyletic taxa and a uniform set of taxon ranks, which make it more appropriate for sequence-based taxonomic assignment algorithms than more accepted classification systems such as @adl2019.
FASTA format files of the re-annotated RDP and Warcup training sets are available at (*somewhere*).

Taxonomic assignment was performed to genus level on the ITS and LSU regions using Unite/Warcup and RDP, respectively, as taxonomic references.
For each region/reference combination, taxonomy was assigned using three algorithms:
the RDP Naïve Bayesian Classifier (NBC) as implemented in DADA2;
SINTAX [@edgar2016a] as implemented in VSEARCH v2.9.1 [@rognes2016];
and IDTAXA [@murali2018].
Each ASV thus was thus given up to nine taxonomic assignments (three region/reference combinations $\times$ three algorithms).

## Alignment and phylogenetic inference

Unique LSU regions from long amplicon ASVs were aligned as RNA using DECIPHER [@wright2015] with up to 10 iterations of progressive alignment and conserved secondary structure calculation and 10 refinement iterations.
This alignment was truncated at a position after the D3 region corresponding to base 907 of the *Saccharomyces cerevisiae* S288C reference sequence for LSU, because several sequences had type 1 introns after this position [@holst-jensen1999].
Full length long amplicon ASVs (including ITS1, 5.8S, and ITS2 regions) were aligned and truncated in the same way.

An attempt was made to refine the long amplicon alignment using Pasta [@mirarab2014a], but the resulting maximum likelihood (ML) trees did not improve on the ML achieved using the original alignment from DECIPHER.
Since ML is the optimization criterion used by Pasta, the original DECIPHER alignment was retained.

ML trees were produced using RAxML version 8.2.12 [@stamatakis2014] using the GTR+GAMMA model and rapid bootstrapping with the MRE_IGN stopping criterion.
The LSU tree was not constrained topologically, but the long amplicon tree was constrained by the result for the LSU tree.

Short amplicon ASVs without a perfect match to one of the long amplicon ASVs were added to the long amplicon alignment using the `--add` option of MAFFT version 7.453 [@katoh2012] with a 10-mer based guide tree and two iterations.
The short ASVs were then placed on the long amplicon tree using EPA-NG version 0.3.6 [@barbera2019].
A "grafted" tree including the short amplicon ASVs in the most probable positions from the EPA-NG results was produced using GAPPA version 0.5.1 [@czech2019a].
When more than one short amplicon ASV was placed in the same position on the tree, GAPPA by default groups all of these sequences together in a polytomy.
The branch leading to each of these polytomies was deleted by setting its length to zero and then using the `d2multi` function from the `ape` package in R.
The resulting tree was used as a guide tree to construct a final ML tree in RAxML, using the GTR+GAMMA model, in order to resolve relationships between different sequences assigned to the same branch.
No bootstrapping was performed on the final tree.

The tree was rooted outside the fungi by using the most abundant ASV which was confidently assigned to the Chlorophyta by all 6 applicable taxonomic assignment methods.
Assignments based on Warcup were not used at this step because non-Fungi are not included in the dataset.
Inspection of the tree along with taxonomic assignments revealed that the genus *Tulasnella* (Basidiomycota), which is known to have accelerated evolution in the rDNA [@moncalvo2006], was nested within the Metazoa, so it was excluded from further analysis.
The kingdom Fungi was then identified as the minimal clade containing all remaining ASVs which were confidently identified (consensus of at least 6 of 9 assignments) to a fungal phylum.
ASVs falling outside this clade were excluded from further analysis.

<!--
          - preliminary alignment of 32S to RDP covariance model with `cmalign` from Infernal (aligns only conserved regions)
          - concatenate ITS1 to beginning of alignment
          - create guide tree using alignment of only conserved, gap-free columns in preliminary alignment, using UPGMA in DECIPHER [@wright2016]
          - realign nonconserved sites using MLocARNA (progressive simultaneous alignment and folding) from LocARNA 2.0.0RC8 [@will2007; @will2012; @will2013]. -->
          
## Statistics

Ecological community dissimilarity matrices were calculated using methods "bray" (both long and short amplicons) and "unifrac" (only long amplicons) in `phyloseq` version 1.26.0.
Each of these distance matrices was used to calculate a Mantel correlogram for distances of 0--12 m, 

# Results

# Discussion

# Bibliography {-}

<div id="refs"></div>