# schirmer-lab/metagear-pipeline: Citations

## [nf-core](https://pubmed.ncbi.nlm.nih.gov/32055031/)

> Ewels PA, Peltzer A, Fillinger S, Patel H, Alneberg J, Wilm A, Garcia MU, Di Tommaso P, Nahnsen S. The nf-core framework for community-curated bioinformatics pipelines. Nat Biotechnol. 2020 Mar;38(3):276-278. doi: 10.1038/s41587-020-0439-x. PubMed PMID: 32055031.

## [Nextflow](https://pubmed.ncbi.nlm.nih.gov/28398311/)

> Di Tommaso P, Chatzou M, Floden EW, Barja PP, Palumbo E, Notredame C. Nextflow enables reproducible computational workflows. Nat Biotechnol. 2017 Apr 11;35(4):316-319. doi: 10.1038/nbt.3820. PubMed PMID: 28398311.

## Pipeline tools

Tools are grouped by the entry-point workflow that runs them. A single run executes one workflow, so only that section applies — the MultiQC report's Methods Description lists the subset actually used for your run. This list is kept in sync with the citation registry in `subworkflows/local/utils_nfcore_metagear_pipeline/main.nf`.

- [MultiQC](https://pubmed.ncbi.nlm.nih.gov/27312411/) — all workflows

> Ewels P, Magnusson M, Lundin S, Käller M. MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics. 2016 Oct 1;32(19):3047-8. doi: 10.1093/bioinformatics/btw354. Epub 2016 Jun 16. PubMed PMID: 27312411; PubMed Central PMCID: PMC5039924.

### Quality control (`qc_dna`, `qc_rna`)

- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

> Andrews S. (2010) FastQC: a quality control tool for high throughput sequence data.

- [Trim Galore!](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)

> Krueger F. Trim Galore!: a wrapper around Cutadapt and FastQC.

- [Cutadapt](https://doi.org/10.14806/ej.17.1.200)

> Martin M. Cutadapt removes adapter sequences from high-throughput sequencing reads. EMBnet.journal. 2011;17(1):10-12. doi: 10.14806/ej.17.1.200.

- [KneadData](https://github.com/biobakery/kneaddata)

> McIver LJ, et al. KneadData: quality control for metagenomic sequencing data.

### Reference-based profiling (`microbial_profiles`)

- [MetaPhlAn 4](https://pubmed.ncbi.nlm.nih.gov/36823356/)

> Blanco-Míguez A, Beghini F, Cumbo F, McIver LJ, Thompson KN, Zolfo M, et al. Extending and improving metagenomic taxonomic profiling with uncharacterized species using MetaPhlAn 4. Nat Biotechnol. 2023;41:1633-1644. doi: 10.1038/s41587-023-01688-w.

- [HUMAnN 3](https://pubmed.ncbi.nlm.nih.gov/33944776/)

> Beghini F, McIver LJ, Blanco-Míguez A, Dubois L, Asnicar F, Maharjan S, et al. Integrating taxonomic, functional, and strain-level profiling of diverse microbial communities with bioBakery 3. eLife. 2021;10:e65088. doi: 10.7554/eLife.65088.

### Assembly, gene calling and catalogs (`genes`, `virus`, `classification`)

- [MEGAHIT](https://pubmed.ncbi.nlm.nih.gov/25609793/)

> Li D, Liu CM, Luo R, Sadakane K, Lam TW. MEGAHIT: an ultra-fast single-node solution for large and complex metagenomics assembly via succinct de Bruijn graph. Bioinformatics. 2015;31(10):1674-6. doi: 10.1093/bioinformatics/btv033.

- [Prodigal](https://pubmed.ncbi.nlm.nih.gov/20211023/)

> Hyatt D, Chen GL, LoCascio PF, Land ML, Larimer FW, Hauser LJ. Prodigal: prokaryotic gene recognition and translation initiation site identification. BMC Bioinformatics. 2010;11:119. doi: 10.1186/1471-2105-11-119.

- [MMseqs2](https://pubmed.ncbi.nlm.nih.gov/29035372/)

> Steinegger M, Söding J. MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets. Nat Biotechnol. 2017;35:1026-1028. doi: 10.1038/nbt.3988.

- [CD-HIT](https://pubmed.ncbi.nlm.nih.gov/23060610/)

> Fu L, Niu B, Zhu Z, Wu S, Li W. CD-HIT: accelerated for clustering the next-generation sequencing data. Bioinformatics. 2012;28(23):3150-2. doi: 10.1093/bioinformatics/bts565.

- [VAMB](https://pubmed.ncbi.nlm.nih.gov/33398153/)

> Nissen JN, Johansen J, Allesøe RL, Sønderby CK, Armenteros JJA, Grønbech CH, et al. Improved metagenome binning and assembly using deep variational autoencoders. Nat Biotechnol. 2021;39:555-560. doi: 10.1038/s41587-020-00777-4.

- [Biopython](https://pubmed.ncbi.nlm.nih.gov/19304878/)

> Cock PJA, Antao T, Chang JT, Chapman BA, Cox CJ, Dalke A, et al. Biopython: freely available Python tools for computational molecular biology and bioinformatics. Bioinformatics. 2009;25(11):1422-3. doi: 10.1093/bioinformatics/btp163.

### Read mapping and abundance (`genes`, `virus`, `mag`)

- [BWA-MEM2](https://doi.org/10.1109/IPDPS.2019.00041)

> Vasimuddin M, Misra S, Li H, Aluru S. Efficient Architecture-Aware Acceleration of BWA-MEM for Multicore Systems. IEEE IPDPS. 2019;314-324. doi: 10.1109/IPDPS.2019.00041.

- [SAMtools](https://pubmed.ncbi.nlm.nih.gov/33590861/)

> Danecek P, Bonfield JK, Liddle J, Marshall J, Ohan V, Pollard MO, et al. Twelve years of SAMtools and BCFtools. Gigascience. 2021;10(2):giab008. doi: 10.1093/gigascience/giab008.

- [CoverM](https://github.com/wwood/CoverM)

> Woodcroft BJ, et al. CoverM: read coverage calculator for metagenomics.

- [SeqKit](https://pubmed.ncbi.nlm.nih.gov/27706213/)

> Shen W, Le S, Li Y, Hu F. SeqKit: A Cross-Platform and Ultrafast Toolkit for FASTA/Q File Manipulation. PLoS One. 2016;11(10):e0163962. doi: 10.1371/journal.pone.0163962.

- [seqtk](https://github.com/lh3/seqtk)

> Li H. seqtk: a fast and lightweight tool for processing FASTA/Q sequences.

### Functional annotation (`genes`, `virus`)

- [InterProScan](https://pubmed.ncbi.nlm.nih.gov/24451626/)

> Jones P, Binns D, Chang HY, Fraser M, Li W, McAnulla C, et al. InterProScan 5: genome-scale protein function classification. Bioinformatics. 2014;30(9):1236-40. doi: 10.1093/bioinformatics/btu031.

- [AMRFinderPlus](https://pubmed.ncbi.nlm.nih.gov/34135355/)

> Feldgarden M, Brover V, Gonzalez-Escalona N, Frye JG, Haendiges J, Haft DH, et al. AMRFinderPlus and the Reference Gene Catalog facilitate examination of the genomic links among antimicrobial resistance, stress response, and virulence. Sci Rep. 2021;11:12728. doi: 10.1038/s41598-021-91456-0.

### Viral and plasmid analysis (`virus`, `classification`)

- [geNomad](https://pubmed.ncbi.nlm.nih.gov/37735266/)

> Camargo AP, Roux S, Schulz F, Babinski M, Xu Y, Hu B, et al. Identification of mobile genetic elements with geNomad. Nat Biotechnol. 2024;42:1303-1312. doi: 10.1038/s41587-023-01953-y.

- [CheckV](https://pubmed.ncbi.nlm.nih.gov/33349699/)

> Nayfach S, Camargo AP, Schulz F, Eloe-Fadrosh E, Roux S, Kyrpides NC. CheckV assesses the quality and completeness of metagenome-assembled viral genomes. Nat Biotechnol. 2021;39:578-585. doi: 10.1038/s41587-020-00774-7.

- [VirSorter2](https://pubmed.ncbi.nlm.nih.gov/33522966/)

> Guo J, Bolduc B, Zayed AA, Varsani A, Dominguez-Huerta G, Delmont TO, et al. VirSorter2: a multi-classifier, expert-guided approach to detect diverse DNA and RNA viruses. Microbiome. 2021;9:37. doi: 10.1186/s40168-020-00990-y.

- [DRAM / DRAM-v](https://pubmed.ncbi.nlm.nih.gov/32766782/)

> Shaffer M, Borton MA, McGivern BB, Zayed AA, La Rosa SL, Solden LM, et al. DRAM for distilling microbial metabolism to automate the curation of microbiome function. Nucleic Acids Res. 2020;48(16):8883-8900. doi: 10.1093/nar/gkaa621.

- [iPHoP](https://pubmed.ncbi.nlm.nih.gov/37083735/)

> Roux S, Camargo AP, Coutinho FH, Dabdoub SM, Dutilh BE, Nayfach S, Tritt A. iPHoP: An integrated machine learning framework to maximize host prediction for metagenome-derived viruses of archaea and bacteria. PLoS Biol. 2023;21(4):e3002083. doi: 10.1371/journal.pbio.3002083.

- [Pharokka](https://pubmed.ncbi.nlm.nih.gov/36453861/)

> Bouras G, Nepal R, Houtak G, Psaltis AJ, Wormald PJ, Vreugde S. Pharokka: a fast scalable bacteriophage annotation tool. Bioinformatics. 2023;39(1):btac776. doi: 10.1093/bioinformatics/btac776.

### Binning and MAGs (`classification`, `mag`)

- [SemiBin2](https://pubmed.ncbi.nlm.nih.gov/37387171/)

> Pan S, Zhao XM, Coelho LP. SemiBin2: self-supervised contrastive learning leads to better MAGs for short- and long-read sequencing. Bioinformatics. 2023;39(Suppl 1):i21-i29. doi: 10.1093/bioinformatics/btad209.

- [MetaBAT 2](https://pubmed.ncbi.nlm.nih.gov/31388474/)

> Kang DD, Li F, Kirton E, Thomas A, Egan R, An H, Wang Z. MetaBAT 2: an adaptive binning algorithm for robust and efficient genome reconstruction from metagenome assemblies. PeerJ. 2019;7:e7359. doi: 10.7717/peerj.7359.

- [Binette](https://doi.org/10.21105/joss.06782)

> Mainguy J, Hoede C. Binette: a fast and accurate bin refinement tool to construct high quality Metagenome Assembled Genomes. J Open Source Softw. 2024;9(102):6782. doi: 10.21105/joss.06782.

- [CheckM2](https://pubmed.ncbi.nlm.nih.gov/37500759/)

> Chklovski A, Parks DH, Woodcroft BJ, Tyson GW. CheckM2: a rapid, scalable and accurate tool for assessing microbial genome quality using machine learning. Nat Methods. 2023;20:1203-1212. doi: 10.1038/s41592-023-01940-w.

- [Tiara](https://pubmed.ncbi.nlm.nih.gov/34570171/)

> Karlicki M, Antonowicz S, Karnkowska A. Tiara: deep learning-based classification system for eukaryotic sequences. Bioinformatics. 2022;38(2):344-350. doi: 10.1093/bioinformatics/btab672.

- [dRep](https://pubmed.ncbi.nlm.nih.gov/28742071/)

> Olm MR, Brown CT, Brooks B, Banfield JF. dRep: a tool for fast and accurate genomic comparisons that enables improved genome recovery from metagenomes. ISME J. 2017;11:2864-2868. doi: 10.1038/ismej.2017.126.

- [skani](https://pubmed.ncbi.nlm.nih.gov/37735570/)

> Shaw J, Yu YW. Fast and robust metagenomic sequence comparison through sparse chaining with skani. Nat Methods. 2023;20:1661-1665. doi: 10.1038/s41592-023-02018-3.

- [GTDB-Tk](https://pubmed.ncbi.nlm.nih.gov/36218463/)

> Chaumeil PA, Mussig AJ, Hugenholtz P, Parks DH. GTDB-Tk v2: memory friendly classification with the genome taxonomy database. Bioinformatics. 2022;38(23):5315-5316. doi: 10.1093/bioinformatics/btac672.

- [GTDB](https://pubmed.ncbi.nlm.nih.gov/34520557/)

> Parks DH, Chuvochina M, Rinke C, Mussig AJ, Chaumeil PA, Hugenholtz P. GTDB: an ongoing census of bacterial and archaeal diversity through a phylogenetically consistent, rank normalized and complete genome-based taxonomy. Nucleic Acids Res. 2022;50(D1):D785-D794. doi: 10.1093/nar/gkab776.

### Metagenomic species pangenomes (`msp`)

- [MSPminer](https://doi.org/10.1093/bioinformatics/bty830)

> Plaza Oñate F, Le Chatelier E, Almeida M, Cervino ACL, Gauthier F, Magoulès F, Ehrlich SD, Pichaud M. MSPminer: abundance-based reconstitution of microbial pan-genomes from shotgun metagenomic data. Bioinformatics. 2019;35(9):1544-1552. doi: 10.1093/bioinformatics/bty830.

### Structural annotation (`structures`)

- [PHOLD](https://github.com/gbouras13/phold)

> Bouras G, et al. PHOLD: phage annotation using protein structures.

- [ProstT5](https://doi.org/10.1093/nargab/lqae150)

> Heinzinger M, Weissenow K, Sanchez JG, Henkel A, Mirdita M, Steinegger M, Rost B. Bilingual language model for protein sequence and structure. NAR Genom Bioinform. 2024;6(4):lqae150. doi: 10.1093/nargab/lqae150.

- [Foldseek](https://pubmed.ncbi.nlm.nih.gov/37156916/)

> van Kempen M, Kim SS, Tumescheit C, Mirdita M, Lee J, Gilchrist CLM, Söding J, Steinegger M. Fast and accurate protein structure search with Foldseek. Nat Biotechnol. 2024;42:243-246. doi: 10.1038/s41587-023-01773-0.

## Software packaging/containerisation tools

- [Anaconda](https://anaconda.com)

  > Anaconda Software Distribution. Computer software. Vers. 2-2.4.0. Anaconda, Nov. 2016. Web.

- [Bioconda](https://pubmed.ncbi.nlm.nih.gov/29967506/)

  > Grüning B, Dale R, Sjödin A, Chapman BA, Rowe J, Tomkins-Tinch CH, Valieris R, Köster J; Bioconda Team. Bioconda: sustainable and comprehensive software distribution for the life sciences. Nat Methods. 2018 Jul;15(7):475-476. doi: 10.1038/s41592-018-0046-7. PubMed PMID: 29967506.

- [BioContainers](https://pubmed.ncbi.nlm.nih.gov/28379341/)

  > da Veiga Leprevost F, Grüning B, Aflitos SA, Röst HL, Uszkoreit J, Barsnes H, Vaudel M, Moreno P, Gatto L, Weber J, Bai M, Jimenez RC, Sachsenberg T, Pfeuffer J, Alvarez RV, Griss J, Nesvizhskii AI, Perez-Riverol Y. BioContainers: an open-source and community-driven framework for software standardization. Bioinformatics. 2017 Aug 15;33(16):2580-2582. doi: 10.1093/bioinformatics/btx192. PubMed PMID: 28379341; PubMed Central PMCID: PMC5870671.

- [Docker](https://dl.acm.org/doi/10.5555/2600239.2600241)

  > Merkel, D. (2014). Docker: lightweight linux containers for consistent development and deployment. Linux Journal, 2014(239), 2. doi: 10.5555/2600239.2600241.

- [Singularity](https://pubmed.ncbi.nlm.nih.gov/28494014/)

  > Kurtzer GM, Sochat V, Bauer MW. Singularity: Scientific containers for mobility of compute. PLoS One. 2017 May 11;12(5):e0177459. doi: 10.1371/journal.pone.0177459. eCollection 2017. PubMed PMID: 28494014; PubMed Central PMCID: PMC5426675.
