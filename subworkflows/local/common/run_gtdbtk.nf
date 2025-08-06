include { GTDBTK_CLASSIFYWF } from "$projectDir/modules/local/gtdbtk/classifywf"

workflow GTDBTK_ANN {

    main:

        meta = [id: "KCH_main_gtdbtk_ann"]

        genome_file = Channel.fromPath("/nfs/arxiv/shen/CLD_KCH_2025/analysis/gene_profile/results/mspminer/pangenome/*.pangenome.fasta")

        Channel
            .of(meta)
            .combine(genome_file)
            .set { genome_ch }


        db_path = "/nfs/data/database/GTDB_tk_database/release207_v2"
        db_ch = Channel.of(tuple(meta, db_path))

        GTDBTK_CLASSIFYWF ( genome_ch, db_ch, "/nfs/arxiv/shen/CLD_KCH_2025/analysis/gene_profile/results/mspminer/pangenome", false)

    emit:
        validated_input = GTDBTK_CLASSIFYWF.out.gtdb_outdir

}
