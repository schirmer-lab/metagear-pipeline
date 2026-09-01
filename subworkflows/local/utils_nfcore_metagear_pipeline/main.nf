//
// Subworkflow with functionality specific to the schirmer-lab/metagear-pipeline pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFSCHEMA_PLUGIN     } from '../../nf-core/utils_nfschema_plugin'
include { paramsSummaryMap          } from 'plugin/nf-schema'
include { samplesheetToList         } from 'plugin/nf-schema'
include { paramsHelp                } from 'plugin/nf-schema'
include { completionEmail           } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'

include { INPUT_CHECK               } from "$projectDir/subworkflows/local/common/input_check"

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {

    take:
    version           // boolean: Display version and exit
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    monochrome_logs   // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where the results will be saved
    input             //  string: Path to input samplesheet
    help              // boolean: Display help message and exit
    help_full         // boolean: Show the full help message
    show_hidden       // boolean: Show hidden parameters in the help message

    main:

    ch_versions = channel.empty()

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE (
        version,
        false,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
    )

    // Dumped here rather than by UTILS_NEXTFLOW_PIPELINE: this file is written
    // straight to disk, so unlike every other pipeline_info output there is no
    // publishDir saveAs hook to prefix it with the workflow name.
    if (outdir) {
        dumpParametersToWorkflowJSON(outdir)
    }

    //
    // Validate parameters and generate parameter summary to stdout
    //

    def before_text = ""
    def after_text = ""
    if (monochrome_logs) {
        before_text = before_text.replaceAll(/\033\[[0-9;]*m/, '')
    }

    command = "nextflow run ${workflow.manifest.name} -profile <docker/singularity/.../institute> --input samplesheet.csv --outdir <OUTDIR>"

    UTILS_NFSCHEMA_PLUGIN (
        workflow,
        validate_params,
        null,
        help,
        help_full,
        show_hidden,
        before_text,
        after_text,
        command
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE (
        nextflow_cli_args
    )

    // NOTE: TEMPLATE constructs a samplesheet channel from `params.input` here
    // and emits it as `PIPELINE_INITIALISATION.out.samplesheet`. Our
    // SCHIRMERLAB workflow doesn't consume that channel (we have our own
    // INPUT_CHECK subworkflow that reads `params.input` directly) so we don't
    // build it. If a future refactor adopts the TEMPLATE-style flow, restore
    // the `channel.fromList(samplesheetToList(...))` block here and the
    // `validateInputSamplesheet` helper at the bottom of this file.

    emit:
    versions    = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_COMPLETION {

    take:
    email           //  string: email address
    email_on_fail   //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML
    outdir          //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    multiqc_report  //  string: Path to MultiQC report

    main:
    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def multiqc_reports = multiqc_report.toList()

    //
    // Completion email and summary
    //
    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(
                summary_params,
                email,
                email_on_fail,
                plaintext_email,
                outdir,
                monochrome_logs,
                multiqc_reports.getVal(),
            )
        }

        completionSummary(monochrome_logs)

    }

    workflow.onError {
        log.error "Pipeline failed. Please refer to troubleshooting docs for common issues: https://nf-co.re/docs/running/troubleshooting"
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// Generate methods description for MultiQC
//
//
// Citation registry. One entry per citable tool: `cite` is the in-text form,
// `ref` the bibliography entry. Defined inside a function rather than as a
// script-level `def` or @Field so it resolves from every calling function
// (see the BIOMES scope note in subworkflows/local/common/input_check.nf).
//
// Tools without a peer-reviewed publication are cited by their canonical URL.
// Keep this in sync with CITATIONS.md.
//
def dumpParametersToWorkflowJSON(outdir) {
    def timestamp = new java.util.Date().format('yyyy-MM-dd_HH-mm-ss')
    def temp_pf   = new File(workflow.launchDir.toString(), ".params_${timestamp}.json")
    temp_pf.text  = groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(params))

    nextflow.extension.FilesEx.copyTo(temp_pf.toPath(), "${outdir}/pipeline_info/${params.workflow}_params_${timestamp}.json")
    temp_pf.delete()
}


def metagearToolRegistry() {
    return [
        amrfinderplus: [cite: 'AMRFinderPlus (Feldgarden et al. 2021)', ref: '<li>Feldgarden, M., Brover, V., Gonzalez-Escalona, N., Frye, J. G., Haendiges, J., Haft, D. H., et al. (2021). AMRFinderPlus and the Reference Gene Catalog facilitate examination of the genomic links among antimicrobial resistance, stress response, and virulence. Scientific Reports, 11, 12728. doi: <a href="https://doi.org/10.1038/s41598-021-91456-0">10.1038/s41598-021-91456-0</a></li>'],
        binette: [cite: 'Binette (Mainguy &amp; Hoede 2024)', ref: '<li>Mainguy, J., &amp; Hoede, C. (2024). Binette: a fast and accurate bin refinement tool to construct high quality Metagenome Assembled Genomes. Journal of Open Source Software, 9(102), 6782. doi: <a href="https://doi.org/10.21105/joss.06782">10.21105/joss.06782</a></li>'],
        biopython: [cite: 'Biopython (Cock et al. 2009)', ref: '<li>Cock, P. J. A., Antao, T., Chang, J. T., Chapman, B. A., Cox, C. J., Dalke, A., et al. (2009). Biopython: freely available Python tools for computational molecular biology and bioinformatics. Bioinformatics, 25(11), 1422–1423. doi: <a href="https://doi.org/10.1093/bioinformatics/btp163">10.1093/bioinformatics/btp163</a></li>'],
        bwamem2: [cite: 'BWA-MEM2 (Vasimuddin et al. 2019)', ref: '<li>Vasimuddin, M., Misra, S., Li, H., &amp; Aluru, S. (2019). Efficient Architecture-Aware Acceleration of BWA-MEM for Multicore Systems. IEEE International Parallel and Distributed Processing Symposium (IPDPS), 314–324. doi: <a href="https://doi.org/10.1109/IPDPS.2019.00041">10.1109/IPDPS.2019.00041</a></li>'],
        cdhit: [cite: 'CD-HIT (Fu et al. 2012)', ref: '<li>Fu, L., Niu, B., Zhu, Z., Wu, S., &amp; Li, W. (2012). CD-HIT: accelerated for clustering the next-generation sequencing data. Bioinformatics, 28(23), 3150–3152. doi: <a href="https://doi.org/10.1093/bioinformatics/bts565">10.1093/bioinformatics/bts565</a></li>'],
        checkm2: [cite: 'CheckM2 (Chklovski et al. 2023)', ref: '<li>Chklovski, A., Parks, D. H., Woodcroft, B. J., &amp; Tyson, G. W. (2023). CheckM2: a rapid, scalable and accurate tool for assessing microbial genome quality using machine learning. Nature Methods, 20, 1203–1212. doi: <a href="https://doi.org/10.1038/s41592-023-01940-w">10.1038/s41592-023-01940-w</a></li>'],
        checkv: [cite: 'CheckV (Nayfach et al. 2021)', ref: '<li>Nayfach, S., Camargo, A. P., Schulz, F., Eloe-Fadrosh, E., Roux, S., &amp; Kyrpides, N. C. (2021). CheckV assesses the quality and completeness of metagenome-assembled viral genomes. Nature Biotechnology, 39, 578–585. doi: <a href="https://doi.org/10.1038/s41587-020-00774-7">10.1038/s41587-020-00774-7</a></li>'],
        coverm: [cite: 'CoverM (Woodcroft et al.)', ref: '<li>Woodcroft, B. J., et al. CoverM: read coverage calculator for metagenomics. <a href="https://github.com/wwood/CoverM">https://github.com/wwood/CoverM</a></li>'],
        cutadapt: [cite: 'Cutadapt (Martin 2011)', ref: '<li>Martin, M. (2011). Cutadapt removes adapter sequences from high-throughput sequencing reads. EMBnet.journal, 17(1), 10–12. doi: <a href="https://doi.org/10.14806/ej.17.1.200">10.14806/ej.17.1.200</a></li>'],
        dramv: [cite: 'DRAM-v (Shaffer et al. 2020)', ref: '<li>Shaffer, M., Borton, M. A., McGivern, B. B., Zayed, A. A., La Rosa, S. L., Solden, L. M., et al. (2020). DRAM for distilling microbial metabolism to automate the curation of microbiome function. Nucleic Acids Research, 48(16), 8883–8900. doi: <a href="https://doi.org/10.1093/nar/gkaa621">10.1093/nar/gkaa621</a></li>'],
        drep: [cite: 'dRep (Olm et al. 2017)', ref: '<li>Olm, M. R., Brown, C. T., Brooks, B., &amp; Banfield, J. F. (2017). dRep: a tool for fast and accurate genomic comparisons that enables improved genome recovery from metagenomes. The ISME Journal, 11, 2864–2868. doi: <a href="https://doi.org/10.1038/ismej.2017.126">10.1038/ismej.2017.126</a></li>'],
        fastqc: [cite: 'FastQC (Andrews 2010)', ref: '<li>Andrews, S. (2010). FastQC: a quality control tool for high throughput sequence data. <a href="https://www.bioinformatics.babraham.ac.uk/projects/fastqc/">https://www.bioinformatics.babraham.ac.uk/projects/fastqc/</a></li>'],
        foldseek: [cite: 'Foldseek (van Kempen et al. 2024)', ref: '<li>van Kempen, M., Kim, S. S., Tumescheit, C., Mirdita, M., Lee, J., Gilchrist, C. L. M., Söding, J., &amp; Steinegger, M. (2024). Fast and accurate protein structure search with Foldseek. Nature Biotechnology, 42, 243–246. doi: <a href="https://doi.org/10.1038/s41587-023-01773-0">10.1038/s41587-023-01773-0</a></li>'],
        genomad: [cite: 'geNomad (Camargo et al. 2024)', ref: '<li>Camargo, A. P., Roux, S., Schulz, F., Babinski, M., Xu, Y., Hu, B., et al. (2024). Identification of mobile genetic elements with geNomad. Nature Biotechnology, 42, 1303–1312. doi: <a href="https://doi.org/10.1038/s41587-023-01953-y">10.1038/s41587-023-01953-y</a></li>'],
        gtdb: [cite: 'GTDB (Parks et al. 2022)', ref: '<li>Parks, D. H., Chuvochina, M., Rinke, C., Mussig, A. J., Chaumeil, P.-A., &amp; Hugenholtz, P. (2022). GTDB: an ongoing census of bacterial and archaeal diversity through a phylogenetically consistent, rank normalized and complete genome-based taxonomy. Nucleic Acids Research, 50(D1), D785–D794. doi: <a href="https://doi.org/10.1093/nar/gkab776">10.1093/nar/gkab776</a></li>'],
        gtdbtk: [cite: 'GTDB-Tk (Chaumeil et al. 2022)', ref: '<li>Chaumeil, P.-A., Mussig, A. J., Hugenholtz, P., &amp; Parks, D. H. (2022). GTDB-Tk v2: memory friendly classification with the genome taxonomy database. Bioinformatics, 38(23), 5315–5316. doi: <a href="https://doi.org/10.1093/bioinformatics/btac672">10.1093/bioinformatics/btac672</a></li>'],
        humann: [cite: 'HUMAnN 3 (Beghini et al. 2021)', ref: '<li>Beghini, F., McIver, L. J., Blanco-Míguez, A., Dubois, L., Asnicar, F., Maharjan, S., et al. (2021). Integrating taxonomic, functional, and strain-level profiling of diverse microbial communities with bioBakery 3. eLife, 10, e65088. doi: <a href="https://doi.org/10.7554/eLife.65088">10.7554/eLife.65088</a></li>'],
        interproscan: [cite: 'InterProScan (Jones et al. 2014)', ref: '<li>Jones, P., Binns, D., Chang, H.-Y., Fraser, M., Li, W., McAnulla, C., et al. (2014). InterProScan 5: genome-scale protein function classification. Bioinformatics, 30(9), 1236–1240. doi: <a href="https://doi.org/10.1093/bioinformatics/btu031">10.1093/bioinformatics/btu031</a></li>'],
        iphop: [cite: 'iPHoP (Roux et al. 2023)', ref: '<li>Roux, S., Camargo, A. P., Coutinho, F. H., Dabdoub, S. M., Dutilh, B. E., Nayfach, S., &amp; Tritt, A. (2023). iPHoP: An integrated machine learning framework to maximize host prediction for metagenome-derived viruses of archaea and bacteria. PLoS Biology, 21(4), e3002083. doi: <a href="https://doi.org/10.1371/journal.pbio.3002083">10.1371/journal.pbio.3002083</a></li>'],
        kneaddata: [cite: 'KneadData (McIver et al.)', ref: '<li>McIver, L. J., et al. KneadData: quality control for metagenomic sequencing data. <a href="https://github.com/biobakery/kneaddata">https://github.com/biobakery/kneaddata</a></li>'],
        megahit: [cite: 'MEGAHIT (Li et al. 2015)', ref: '<li>Li, D., Liu, C.-M., Luo, R., Sadakane, K., &amp; Lam, T.-W. (2015). MEGAHIT: an ultra-fast single-node solution for large and complex metagenomics assembly via succinct de Bruijn graph. Bioinformatics, 31(10), 1674–1676. doi: <a href="https://doi.org/10.1093/bioinformatics/btv033">10.1093/bioinformatics/btv033</a></li>'],
        metabat2: [cite: 'MetaBAT 2 (Kang et al. 2019)', ref: '<li>Kang, D. D., Li, F., Kirton, E., Thomas, A., Egan, R., An, H., &amp; Wang, Z. (2019). MetaBAT 2: an adaptive binning algorithm for robust and efficient genome reconstruction from metagenome assemblies. PeerJ, 7, e7359. doi: <a href="https://doi.org/10.7717/peerj.7359">10.7717/peerj.7359</a></li>'],
        metaphlan: [cite: 'MetaPhlAn 4 (Blanco-Míguez et al. 2023)', ref: '<li>Blanco-Míguez, A., Beghini, F., Cumbo, F., McIver, L. J., Thompson, K. N., Zolfo, M., et al. (2023). Extending and improving metagenomic taxonomic profiling with uncharacterized species using MetaPhlAn 4. Nature Biotechnology, 41, 1633–1644. doi: <a href="https://doi.org/10.1038/s41587-023-01688-w">10.1038/s41587-023-01688-w</a></li>'],
        mmseqs2: [cite: 'MMseqs2 (Steinegger &amp; Söding 2017)', ref: '<li>Steinegger, M., &amp; Söding, J. (2017). MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets. Nature Biotechnology, 35, 1026–1028. doi: <a href="https://doi.org/10.1038/nbt.3988">10.1038/nbt.3988</a></li>'],
        mspminer: [cite: 'MSPminer (Plaza Oñate et al. 2019)', ref: '<li>Plaza Oñate, F., Le Chatelier, E., Almeida, M., Cervino, A. C. L., Gauthier, F., Magoulès, F., Ehrlich, S. D., &amp; Pichaud, M. (2019). MSPminer: abundance-based reconstitution of microbial pan-genomes from shotgun metagenomic data. Bioinformatics, 35(9), 1544–1552. doi: <a href="https://doi.org/10.1093/bioinformatics/bty830">10.1093/bioinformatics/bty830</a></li>'],
        multiqc: [cite: 'MultiQC (Ewels et al. 2016)', ref: '<li>Ewels, P., Magnusson, M., Lundin, S., &amp; Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics, 32(19), 3047–3048. doi: <a href="https://doi.org/10.1093/bioinformatics/btw354">10.1093/bioinformatics/btw354</a></li>'],
        pharokka: [cite: 'Pharokka (Bouras et al. 2023)', ref: '<li>Bouras, G., Nepal, R., Houtak, G., Psaltis, A. J., Wormald, P.-J., &amp; Vreugde, S. (2023). Pharokka: a fast scalable bacteriophage annotation tool. Bioinformatics, 39(1), btac776. doi: <a href="https://doi.org/10.1093/bioinformatics/btac776">10.1093/bioinformatics/btac776</a></li>'],
        phold: [cite: 'PHOLD (Bouras et al.)', ref: '<li>Bouras, G., et al. PHOLD: phage annotation using protein structures. <a href="https://github.com/gbouras13/phold">https://github.com/gbouras13/phold</a></li>'],
        prodigal: [cite: 'Prodigal (Hyatt et al. 2010)', ref: '<li>Hyatt, D., Chen, G.-L., LoCascio, P. F., Land, M. L., Larimer, F. W., &amp; Hauser, L. J. (2010). Prodigal: prokaryotic gene recognition and translation initiation site identification. BMC Bioinformatics, 11, 119. doi: <a href="https://doi.org/10.1186/1471-2105-11-119">10.1186/1471-2105-11-119</a></li>'],
        prostt5: [cite: 'ProstT5 (Heinzinger et al. 2024)', ref: '<li>Heinzinger, M., Weissenow, K., Sanchez, J. G., Henkel, A., Mirdita, M., Steinegger, M., &amp; Rost, B. (2024). Bilingual language model for protein sequence and structure. NAR Genomics and Bioinformatics, 6(4), lqae150. doi: <a href="https://doi.org/10.1093/nargab/lqae150">10.1093/nargab/lqae150</a></li>'],
        samtools: [cite: 'SAMtools (Danecek et al. 2021)', ref: '<li>Danecek, P., Bonfield, J. K., Liddle, J., Marshall, J., Ohan, V., Pollard, M. O., et al. (2021). Twelve years of SAMtools and BCFtools. GigaScience, 10(2), giab008. doi: <a href="https://doi.org/10.1093/gigascience/giab008">10.1093/gigascience/giab008</a></li>'],
        semibin2: [cite: 'SemiBin2 (Pan et al. 2023)', ref: '<li>Pan, S., Zhao, X.-M., &amp; Coelho, L. P. (2023). SemiBin2: self-supervised contrastive learning leads to better MAGs for short- and long-read sequencing. Bioinformatics, 39(Suppl 1), i21–i29. doi: <a href="https://doi.org/10.1093/bioinformatics/btad209">10.1093/bioinformatics/btad209</a></li>'],
        seqkit: [cite: 'SeqKit (Shen et al. 2016)', ref: '<li>Shen, W., Le, S., Li, Y., &amp; Hu, F. (2016). SeqKit: A Cross-Platform and Ultrafast Toolkit for FASTA/Q File Manipulation. PLoS ONE, 11(10), e0163962. doi: <a href="https://doi.org/10.1371/journal.pone.0163962">10.1371/journal.pone.0163962</a></li>'],
        seqtk: [cite: 'seqtk (Li)', ref: '<li>Li, H. seqtk: a fast and lightweight tool for processing FASTA/Q sequences. <a href="https://github.com/lh3/seqtk">https://github.com/lh3/seqtk</a></li>'],
        skani: [cite: 'skani (Shaw &amp; Yu 2023)', ref: '<li>Shaw, J., &amp; Yu, Y. W. (2023). Fast and robust metagenomic sequence comparison through sparse chaining with skani. Nature Methods, 20, 1661–1665. doi: <a href="https://doi.org/10.1038/s41592-023-02018-3">10.1038/s41592-023-02018-3</a></li>'],
        tiara: [cite: 'Tiara (Karlicki et al. 2022)', ref: '<li>Karlicki, M., Antonowicz, S., &amp; Karnkowska, A. (2022). Tiara: deep learning-based classification system for eukaryotic sequences. Bioinformatics, 38(2), 344–350. doi: <a href="https://doi.org/10.1093/bioinformatics/btab672">10.1093/bioinformatics/btab672</a></li>'],
        trimgalore: [cite: 'Trim Galore! (Krueger)', ref: '<li>Krueger, F. Trim Galore!: a wrapper around Cutadapt and FastQC. <a href="https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/">https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/</a></li>'],
        vamb: [cite: 'VAMB (Nissen et al. 2021)', ref: '<li>Nissen, J. N., Johansen, J., Allesøe, R. L., Sønderby, C. K., Armenteros, J. J. A., Grønbech, C. H., et al. (2021). Improved metagenome binning and assembly using deep variational autoencoders. Nature Biotechnology, 39, 555–560. doi: <a href="https://doi.org/10.1038/s41587-020-00777-4">10.1038/s41587-020-00777-4</a></li>'],
        virsorter2: [cite: 'VirSorter2 (Guo et al. 2021)', ref: '<li>Guo, J., Bolduc, B., Zayed, A. A., Varsani, A., Dominguez-Huerta, G., Delmont, T. O., et al. (2021). VirSorter2: a multi-classifier, expert-guided approach to detect diverse DNA and RNA viruses. Microbiome, 9, 37. doi: <a href="https://doi.org/10.1186/s40168-020-00990-y">10.1186/s40168-020-00990-y</a></li>'],
    ]
}

//
// Tools actually executed by each entry-point workflow. MetaGEAR dispatches to
// exactly one of these per run, so the methods text must describe only that
// workflow's tools — a qc_dna run should not cite iPHoP.
//
def metagearWorkflowTools() {
    return [
        download_databases: ['kneaddata', 'metaphlan', 'humann', 'gtdbtk', 'gtdb', 'genomad', 'checkv', 'virsorter2', 'dramv', 'iphop', 'pharokka', 'checkm2', 'mmseqs2', 'amrfinderplus', 'phold'],
        qc: ['fastqc', 'trimgalore', 'cutadapt', 'kneaddata', 'multiqc'],
        microbial_profiles: ['metaphlan', 'humann', 'multiqc'],
        genes: ['megahit', 'prodigal', 'vamb', 'mmseqs2', 'cdhit', 'biopython', 'bwamem2', 'samtools', 'coverm', 'seqkit', 'interproscan', 'amrfinderplus', 'multiqc'],
        virus: ['megahit', 'genomad', 'checkv', 'virsorter2', 'dramv', 'iphop', 'pharokka', 'prodigal', 'vamb', 'mmseqs2', 'cdhit', 'bwamem2', 'samtools', 'coverm', 'seqkit', 'seqtk', 'interproscan', 'amrfinderplus', 'multiqc'],
        classification: ['megahit', 'genomad', 'checkv', 'semibin2', 'metabat2', 'binette', 'checkm2', 'tiara', 'mmseqs2', 'vamb', 'bwamem2', 'samtools', 'seqtk', 'multiqc'],
        mag: ['drep', 'skani', 'gtdbtk', 'gtdb', 'coverm', 'bwamem2', 'samtools', 'multiqc'],
        msp: ['mspminer', 'gtdbtk', 'gtdb', 'metaphlan', 'multiqc'],
        structures: ['phold', 'prostt5', 'foldseek', 'seqkit', 'multiqc'],
    ]
}

//
// Resolve params.workflow to a key in metagearWorkflowTools(). Returns [] for
// an unrecognised value so the methods section degrades to the generic text
// rather than asserting tools that never ran.
//
def metagearActiveTools() {
    def selected = (params.workflow ?: '').trim()
    def key = selected.startsWith('qc_') ? 'qc' : selected
    return metagearWorkflowTools()[key] ?: []
}

def toolCitationText() {
    def registry = metagearToolRegistry()
    def tools = metagearActiveTools().findAll { registry.containsKey(it) }

    if (!tools) {
        return ''
    }

    def cites = tools.collect { registry[it].cite }
    return "Tools used in the workflow included: ${cites.join(', ')}."
}

def toolBibliographyText() {
    def registry = metagearToolRegistry()
    def tools = metagearActiveTools().findAll { registry.containsKey(it) }

    return tools.collect { registry[it].ref }.join(' ').trim()
}

def methodsDescriptionText(mqc_methods_yaml) {
    // Convert  to a named map so can be used as with familiar NXF ${workflow} variable syntax in the MultiQC YML file
    def meta = [:]
    meta.workflow = workflow.toMap()
    meta["manifest_map"] = workflow.manifest.toMap()

    // Pipeline DOI
    if (meta.manifest_map.doi) {
        // Using a loop to handle multiple DOIs
        // Removing `https://doi.org/` to handle pipelines using DOIs vs DOI resolvers
        // Removing ` ` since the manifest.doi is a string and not a proper list
        def temp_doi_ref = ""
        def manifest_doi = meta.manifest_map.doi.tokenize(",")
        manifest_doi.each { doi_ref ->
            temp_doi_ref += "(doi: <a href=\'https://doi.org/${doi_ref.replace("https://doi.org/", "").replace(" ", "")}\'>${doi_ref.replace("https://doi.org/", "").replace(" ", "")}</a>), "
        }
        meta["doi_text"] = temp_doi_ref.substring(0, temp_doi_ref.length() - 2)
    } else meta["doi_text"] = ""
    meta["nodoi_text"] = meta.manifest_map.doi ? "" : "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

    // Tool references, scoped to the tools the selected workflow actually runs.
    meta["tool_citations"] = toolCitationText().replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
    meta["tool_bibliography"] = toolBibliographyText()


    def methods_text = mqc_methods_yaml.text

    def engine =  new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}
