#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.expression_file = "$PWD/raw_data/*.bed"
params.genotype_file = "$PWD/raw_data/CEDAR_sc.vcf.gz"
params.metadata = "$PWD/raw_data/Covariates_biopsy.csv"

params.outdir = "$PWD/results"

params.run_trans_eqtl_analysis = false
params.run_trans_eqtl_permute = false
params.peer_approach = false //false - standard (iterational), true - use all peer factors

params.minor_AF = "0.05"
params.grouping = ""
params.cis_Noperm = "10000"
params.cis_window = "1000000"
params.trans_threshold = "0.0005"
params.trans_window = "5000000"
params.normalization = "median_Helene"
params.type_of_hidden = "PCA"
params.No_peer_factors = "auto"
params.built_in_normalization = "--normal"
params.number_of_peer_iterations = "1000"


params.script_merge_check = "$PWD/rscripts/merge_check.R"
params.script_normalization = "$PWD/rscripts/normalization.R"
params.script_peer = "$PWD/rscripts/peer_correction.R"
params.script_peer_2 = "$PWD/rscripts/peer_correction_2.R"
params.script_metaresid = "$PWD/rscripts/metadata_residualization.R"
params.script_cisresid = "$PWD/rscripts/cis_residualization.R"
params.script_cis_FDR = "$PWD/rscripts/qtltools_cis_fdr.R"
params.script_trans_FDR = "$PWD/rscripts/fdr_trans.R"
params.script_snps_extraction = "$PWD/rscripts/snps_extraction.R"
params.chunks_50 = "$PWD/rscripts/chunks_50.csv"
params.chunks_100 = "$PWD/rscripts/chunks_100.csv"
params.chunks_500 = "$PWD/rscripts/chunks_500.csv"
params.chunks_1000 = "$PWD/rscripts/chunks_1000.csv"


Ch_peer_test_50 = Channel
                         .fromPath(params.chunks_50)
                         .splitCsv(header:false,sep:',')
                         .map{it -> it[0]}
                              
Ch_trans_test_100 = Channel
                           .fromPath(params.chunks_100)
                           .splitCsv(header:false,sep:',')
                           .map{it -> it[0]}

Ch_trans_nom_test_500 = Channel
                               .fromPath(params.chunks_500)
                               .splitCsv(header:false,sep:',')
                               .map{it -> it[0]}
                        
Ch_trans_nom_test_1000 = Channel
                               .fromPath(params.chunks_1000)
                               .splitCsv(header:false,sep:',')
                               .map{it -> it[0]}


Ch_combined_seed_chunk_500 = Channel.from('1','2','3','4','5','6','7','8','9','10').combine(Ch_peer_test_50)
Ch_combined_chunk_seed_1000 = Channel
                                     .from('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','23')
                                     .combine(Ch_trans_test_100)


println """\
	    E Q T L - N F  P I P E L I N E
	    ==============================
	    transcriptome: ${params.expression_file}
	    genotype: ${params.genotype_file}
	    metadata: ${params.metadata}
	    normalization: ${params.normalization}
	    Hidden covariates: ${params.type_of_hidden}
	    Hidden covariates No: ${params.No_peer_factors}
	    outdir: ${params.outdir}
	    trans-eQTL analysis: ${params.run_trans_eqtl_analysis}
	"""
	.stripIndent()
	
	
Ch_unnorm_expression = Channel
                              .fromPath(params.expression_file)
                              .map { file -> tuple(file.baseName, file) }

process merge_check {
    
    label 'medium' 
        
    input:
    
    path genotype
    path metadata
    tuple val(sample_ID), path(sample_path)
    path script_merge_check
    val(MAF)
    
    
      
    
    output:
    
    tuple val(sample_ID), path("genotype_merged_${sample_ID}.vcf.gz"), path("genotype_merged_${sample_ID}.vcf.gz.tbi")
    tuple val(sample_ID), path("metadata_merged_${sample_ID}")
    
    
    script:
    
    """
    module load tabix
    module load bcftools
    module load R
    
    Rscript ${script_merge_check} ${sample_path} ${metadata} ${sample_ID}

    bcftools view --min-af ${MAF}:minor -Oz -S sample_list ${genotype} > genotype_merged_${sample_ID}.vcf.gz
    tabix -p vcf genotype_merged_${sample_ID}.vcf.gz
    
       
    """

}



process normalize_expression {

    label 'short'

    
    input:
    
    tuple val(sample_ID), path(sample_path)
    path script_normalization
    val norm
        
        
    publishDir "$params.outdir/normalized_expression_${sample_ID}", mode: 'copy'
    
    
    output:
    
    tuple val(sample_ID), path("expression_normalized_${sample_ID}.bed")
    
    
    script:
    
    """
    module load R
    Rscript ${script_normalization} ${sample_path} ${norm} expression_normalized_${sample_ID}.bed
    
    """
    
    
}


process peer_computation {
    
    label 'very_long'  
    container 'peer_singularity.sif'
    
    input:
    
    tuple val(sample_ID), path(expression), path(metadata)
    path script_peer    
    val(Hidden)
    
    publishDir "$params.outdir/peer_factors_${sample_ID}", mode: 'copy'
    
    output:
    
    tuple val(sample_ID), path("peer_factors_${sample_ID}") 
    
    
    script:
    """
    Rscript ${script_peer} expression_normalized_${sample_ID}.bed metadata_merged_${sample_ID} peer_factors_${sample_ID} peer_model ${Hidden}
    
    """
    
    
}



process peer_test_indexing {

    label 'medium'
         
    input:
    
    tuple val(sample_ID), path(peer_factors)
    
       
    output:
    
    tuple val(sample_ID), path("peer_factors_${sample_ID}.bed.gz"), path("peer_factors_${sample_ID}.bed.gz.tbi")
    
    
    script:
    
    """
    module load tabix
    
    mv ${peer_factors} peer_factors_${sample_ID}.bed
    bgzip -f peer_factors_${sample_ID}.bed
    tabix -f -p bed peer_factors_${sample_ID}.bed.gz
    
    """
    
    
}



process peer_test_run_nominal {
    
    label 'cpu_intensive'
    stageInMode = 'symlink'
    
    input:
    
    tuple val(sample_ID), path(peer_factors), path(peer_index), path(genotypes), path(genotypes_index)
    each chunk
    val(built_in_normalization)
       
    output:
    
    tuple val(sample_ID), path("*hits.txt.gz")
    
    
    script:
    
    """
    module load qtltools
    
    QTLtools trans --vcf genotype_merged_${sample_ID}.vcf.gz --bed peer_factors_${sample_ID}.bed.gz --nominal --chunk ${chunk} 50 --threshold 0.01 --window 1 --out trans.nominal.chunk_${chunk}_${sample_ID} ${built_in_normalization}
    
    """
    
    
}



process peer_test_run_nominal_aggregate{

    label 'medium'
    
    input:
    tuple val(sample_ID), path(peer_chunks)
    
    
    output:
    
    tuple val(sample_ID), path("trans.nominal_${sample_ID}_hits.txt.gz")
    
    
    script:
    
    """
    zcat trans.nominal.chunk_*_${sample_ID}* > trans.nominal_${sample_ID}_hits.txt
    gzip -c trans.nominal_${sample_ID}_hits.txt > trans.nominal_${sample_ID}_hits.txt.gz
        
    """


}



process peer_test_run_permutational {
    
    label 'cpu_intensive'
    stageInMode = 'symlink'
    
    input:
    
    tuple val(sample_ID), path(peer_factors), path(peer_index), path(genotypes), path(genotypes_index), val(seed), val(chunk)
    val(built_in_normalization)
    
    output:
    
    tuple val(sample_ID), path("*best.txt.gz")
    
    
    script:
    
    """
    module load qtltools
    
    
    QTLtools trans --vcf genotype_merged_${sample_ID}.vcf.gz --bed peer_factors_${sample_ID}.bed.gz --permute --chunk ${chunk} 50 --seed ${seed} --threshold 0.01 --window 1 --out trans.perm.chunk_${chunk}_${seed}_${sample_ID} ${built_in_normalization}
      
    """
    
    
}



process peer_test_run_permutational_aggregate{

    label 'medium'
    
    input:
    tuple val(sample_ID), path(peer_chunks)
    
    
    output:
    
    tuple val(sample_ID), path("trans.perm_${sample_ID}_best.txt.gz")
    
    
    script:
    
    """
    zcat trans.perm.chunk_*_${sample_ID}* > trans.perm_${sample_ID}_best.txt
    gzip -c trans.perm_${sample_ID}_best.txt > trans.perm_${sample_ID}_best.txt.gz
        
    """


}


process metadata_residualize {
    
    label 'short'
    
    input:
    
    tuple val(sample_ID), path(peer_test_norm), path(peer_test_perm), path(expression), path(metadata), path(peer_data)
    path(script_resid)
    val(No_peer_factors)
    
    publishDir "$params.outdir/residualised_expression_${sample_ID}_peer_${No_peer_factors}", mode: 'copy'
    
    output:
    
    tuple val(sample_ID), path("expression_residualised_${sample_ID}")
    tuple path("peer_factors_${sample_ID}_ass_w_genotype"), path("peer_factors_${sample_ID}_ass_names")
    
    script:
    
    """
    module load R
    
    Rscript ${script_resid} trans.nominal_${sample_ID}_hits.txt.gz trans.perm_${sample_ID}_best.txt.gz expression_normalized_${sample_ID}.bed metadata_merged_${sample_ID} peer_factors_${sample_ID} expression_residualised_${sample_ID} ${No_peer_factors}
    
    """
    
    
}


process peer_computation_2 {
    
    label 'very_long'  
    container 'peer_singularity.sif'
    
    input:
    
    tuple val(sample_ID), path(expression), path(metadata)
    path script_peer_2
    val(number_of_peer_iterations)
    
    publishDir "$params.outdir/residualised_expression_${sample_ID}_full_peer", mode: 'copy'
    
    output:
    
    tuple val(sample_ID), path("expression_residualised_${sample_ID}") 
    path("peer_model.pdf")
    
    script:
    """
    Rscript ${script_peer_2} expression_normalized_${sample_ID}.bed metadata_merged_${sample_ID} expression_residualised_${sample_ID} ${number_of_peer_iterations}
    
    """
    
    
}


process cis_eQTLs_indexing {
    
    label 'medium'
    
    input:
    
    tuple val(sample_ID), path(expression_residualised)
    
       
    output:
    
    tuple val(sample_ID), path("expression_residualised_${sample_ID}.bed.gz"), path("expression_residualised_${sample_ID}.bed.gz.tbi")
    
    
    script:
    
    """
    module load tabix
    
    mv expression_residualised_${sample_ID} expression_residualised_${sample_ID}.bed
    bgzip -f expression_residualised_${sample_ID}.bed
    tabix -f -p bed expression_residualised_${sample_ID}.bed.gz
    
    """
    
    
}


process cis_eQTLs_nominal {
    
    label 'cpu_intensive'
    stageInMode = 'symlink'
    
    input:
    
    tuple val(sample_ID), path(expression_file), path(expression_index), path(genotype_file), path(genotype_index)
    each chunk
    val(cis_window)
    val(built_in_normalization)
        
    output:
    
    tuple val(sample_ID), path("cis.nominal.chunk*")
    
    
    script:
    
    """
    module load qtltools
    
    QTLtools cis --vcf genotype_merged_${sample_ID}.vcf.gz --bed expression_residualised_${sample_ID}.bed.gz --nominal 1 --window ${cis_window} --std-err --chunk ${chunk} 500 --out cis.nominal.chunk_${sample_ID}_${chunk} ${built_in_normalization}
    
    """
   
    
}



process cis_eQTLs_nominal_agg {
    
    label 'medium'
    
    input:
    
    tuple val(sample_ID), path(results)
    val(No_peer_factors)
    
    publishDir "$params.outdir/cis_eQTLs_nominal_${sample_ID}_peer_${No_peer_factors}", mode: 'copy'
    
    
    output:
    
    tuple val(sample_ID), path("cis.nominal_${sample_ID}.txt.gz")
    
    
    script:
    
    """
    module load R
    
    cat cis.nominal.chunk_${sample_ID}* > cis.nominal_${sample_ID}.txt
    gzip -c cis.nominal_${sample_ID}.txt > cis.nominal_${sample_ID}.txt.gz

        
    """
    
    
}


process cis_eQTLs_permutational {
    
    label 'cpu_intensive'
    stageInMode = 'symlink'
    
    input:
    
    tuple val(sample_ID), path(expression_file), path(expression_index), path(genotype_file), path(genotype_index)
    each chunk
    val(grouping)
    val(Noperm)
    val(cis_window)
    val(built_in_normalization)
        
    output:
    
    tuple val(sample_ID), path("cis.permutational.chunk*")
    
    
    script:
    
    """
    module load qtltools
    
    QTLtools cis --vcf genotype_merged_${sample_ID}.vcf.gz --bed expression_residualised_${sample_ID}.bed.gz --permute ${Noperm} --window ${cis_window} --std-err --chunk ${chunk} 500 --out cis.permutational.chunk_${chunk}_${sample_ID} ${grouping} ${built_in_normalization}
    
    """
   
    
}



process cis_eQTLs_permutational_agg {
    
    label 'medium'
    
    input:
    
    tuple val(sample_ID), path(results)
    path(cis_FDR)
    val(No_peer_factors)
    
    publishDir "$params.outdir/cis_eQTLs_permutational_${sample_ID}_peer_${No_peer_factors}", mode: 'copy'
    
    
    output:
    
    tuple val(sample_ID), path("cis.permutational_${sample_ID}_FDR.thresholds.txt")
    tuple val(sample_ID), path("cis.permutational_${sample_ID}_FDR.significant.txt")
    tuple val(sample_ID), path("cis.permutational_${sample_ID}.txt.gz")
    
    
    script:
    
    """
    module load R
    
    cat cis.permutational.chunk_*_${sample_ID}* > cis.permutational_${sample_ID}.txt
    gzip -c cis.permutational_${sample_ID}.txt > cis.permutational_${sample_ID}.txt.gz
    
    Rscript ${cis_FDR} cis.permutational_${sample_ID}.txt.gz 0.05 cis.permutational_${sample_ID}_FDR
        
    """
    
    
}



process cis_eQTLs_conditional {
    
    label 'mem_cons'
    
    input:
    tuple val(sample_ID), path(expression_file), path(expression_index), path(genotype_file), path(genotype_index), path(cis_thresholds)
    val(grouping)
    val(No_peer_factors)
    val(built_in_normalization)

    publishDir "$params.outdir/cis_eQTLs_conditional_${sample_ID}_peer_${No_peer_factors}", mode: 'copy'
        
    output:
    tuple val(sample_ID), path("cis.conditional_${sample_ID}")
        
    script:
    
    """
    module load qtltools
    
    QTLtools cis --vcf genotype_merged_${sample_ID}.vcf.gz --bed expression_residualised_${sample_ID}.bed.gz --mapping cis.permutational_${sample_ID}_FDR.thresholds.txt --out cis.conditional_${sample_ID} ${grouping} ${built_in_normalization}
    
    """
    
    
}



process cis_eQTLs_for_correct {
    
    label 'medium'
    
    input:
    tuple val(sample_ID), path(permutational_results), path(genotype_file), path(genotype_index)
    path(script_snps_extraction)
         
    output:
    tuple val(sample_ID), path("to_correct_expr_${sample_ID}.vcf.gz")
    
    when:
    params.run_trans_eqtl_analysis
    
    script:
    
    """
    module load R
    module load bcftools
    
    Rscript ${script_snps_extraction} cis.permutational_${sample_ID}_FDR.significant.txt sniplist_${sample_ID}
    bcftools view -Oz -T sniplist_${sample_ID} genotype_merged_${sample_ID}.vcf.gz > to_correct_expr_${sample_ID}.vcf.gz
        
    """
    
    
}



process cis_residualise {
    
    label 'mem_cons'
    
    input:
    
    tuple val(sample_ID), path(residualised_expression), path(vcf_to_correct), path(permutational_results) 
    path(script_cis_resid)
    val(No_peer_factors)
    
    publishDir "$params.outdir/expression_cis_residualized_${sample_ID}_peer_${No_peer_factors}", mode: 'copy'
    
    output:
    tuple val(sample_ID), path("expression_cis_residualised_${sample_ID}")

    when:
    params.run_trans_eqtl_analysis 
    
    script:
    
    """
    module load R
    
    Rscript ${script_cis_resid} cis.permutational_${sample_ID}_FDR.significant.txt expression_residualised_${sample_ID} to_correct_expr_${sample_ID}.vcf.gz expression_cis_residualised_${sample_ID}
    
    """
    
    
}


process trans_eQTLs_indexing {
    
    label 'medium'
    
    input:
    tuple val(sample_ID), path(residualised_expression)
    
       
    output:
    tuple val(sample_ID), path("expression_cis_residualised_${sample_ID}.bed.gz"), path("expression_cis_residualised_${sample_ID}.bed.gz.tbi")

    when:
    params.run_trans_eqtl_analysis 
    
    script:
    
    """
    module load tabix
    
    mv expression_cis_residualised_${sample_ID} expression_cis_residualised_${sample_ID}.bed
    bgzip -f expression_cis_residualised_${sample_ID}.bed
    tabix -f -p bed expression_cis_residualised_${sample_ID}.bed.gz
    
    """
    
    
}


process trans_eQTLs_nominal {
    
    label 'cpu_intensive'
    stageInMode = 'symlink'
    
    input:
    
    tuple val(sample_ID), path(expression_file), path(expression_index), path(genotype_file), path(genotype_index)
    each chunk
    val(trans_threshold)
    val(trans_window)
    val(built_in_normalization)
            
    output:
    tuple val(sample_ID), path("trans.nominal.chunk*hits.txt.gz")
    tuple val(sample_ID), path("trans.nominal.chunk*bins.txt.gz")
    tuple val(sample_ID), path("trans.nominal.chunk*best.txt.gz")

    when:
    params.run_trans_eqtl_analysis 
    
    script:
    
    """
    module load qtltools
    
    QTLtools trans --vcf genotype_merged_${sample_ID}.vcf.gz --bed expression_cis_residualised_${sample_ID}.bed.gz --nominal --chunk ${chunk} 1000 --threshold ${trans_threshold} --window ${trans_window} --out trans.nominal.chunk_${sample_ID}_${chunk} ${built_in_normalization}
        
    """
   
    
}


process trans_eQTLs_nominal_agg_hits {
    
    label 'medium'
    
    input:
    
    tuple val(sample_ID), path(chunks_trans_hits)
    val(No_peer_factors)
        
    publishDir "$params.outdir/trans_eQTLs_nominal_results_${sample_ID}_peer_${No_peer_factors}", mode: 'copy'
        
    output:
    tuple val(sample_ID), path("trans.nominal_${sample_ID}_hits.txt.gz")

    when:
    params.run_trans_eqtl_analysis    
    
    script:
    
    """
    cat trans.nominal.chunk_${sample_ID}*.hits.txt.gz > trans.nominal_${sample_ID}_hits.txt.gz
        
    """
    
    
}


process trans_eQTLs_nominal_agg_fdr {
    
    label 'medium'
    
    input:
    
    tuple val(sample_ID), path(trans_results), path(genotype_file), path(genotype_index), path(expression_normalized)
    path(script_fdr)
    val(No_peer_factors)
        
    publishDir "$params.outdir/trans_eQTLs_nominal_results_${sample_ID}_peer_${No_peer_factors}", mode: 'copy'
        
    output:
    tuple val(sample_ID), path("trans.nominal_${sample_ID}_fdr.tsv")

    when:
    params.run_trans_eqtl_analysis    
    
    script:
    
    """
    module load R
    module load bcftools
    
    n_SNPs="\$(bcftools query -f '%POS\n' genotype_merged_${sample_ID}.vcf.gz | wc -l)"
    Rscript ${script_fdr} trans.nominal_${sample_ID}_hits.txt.gz expression_normalized_${sample_ID}.bed \${n_SNPs} trans.nominal_${sample_ID}_fdr.tsv

    
    """
    
    
}

process trans_eQTLs_nominal_agg_bins {
    
    label 'medium'
    
    input:
    
    tuple val(sample_ID), path(chunks_trans_hits)
    val(No_peer_factors)
        
    publishDir "$params.outdir/trans_eQTLs_nominal_results_${sample_ID}_peer_${No_peer_factors}", mode: 'copy'
        
    output:
    tuple val(sample_ID), path("trans.nominal_${sample_ID}_bins.txt.gz")

    when:
    params.run_trans_eqtl_analysis    
    
    script:
    
    """
    cat trans.nominal.chunk_${sample_ID}*.bins.txt.gz > trans.nominal_${sample_ID}_bins.txt.gz
       
    """
    
    
}



process trans_eQTLs_nominal_agg_best {
    
    label 'medium'
    
    input:
    
    tuple val(sample_ID), path(chunks_trans_hits)
    val(No_peer_factors)
        
    publishDir "$params.outdir/trans_eQTLs_nominal_results_${sample_ID}_peer_${No_peer_factors}", mode: 'copy'
        
    output:
    tuple val(sample_ID), path("trans.nominal_${sample_ID}_best.txt.gz")

    when:
    params.run_trans_eqtl_analysis    
    
    script:
    
    """
    cat trans.nominal.chunk_${sample_ID}*.best.txt.gz > trans.nominal_${sample_ID}_best.txt.gz
        
    """
    
    
}




process trans_eQTLs_permutational {
    
    label 'cpu_intensive'
    stageInMode = 'symlink'
    
    input:
    
    tuple val(sample_ID), path(expression_file), path(expression_index), path(genotype_file), path(genotype_index), val(chunk), val(seed)
    val(trans_threshold)
    val(trans_window)
    val(built_in_normalization)
        
        
    output:
    
    tuple val(sample_ID), path("trans_${sample_ID}_perm.chunk*hits.txt.gz")
    tuple val(sample_ID), path("trans_${sample_ID}_perm.chunk*bins.txt.gz")
    tuple val(sample_ID), path("trans_${sample_ID}_perm.chunk*best.txt.gz")
    
    when:
    params.run_trans_eqtl_analysis && params.run_trans_eqtl_permute
    
    script:
    
    """
    module load qtltools
    
    
    QTLtools trans --vcf genotype_merged_${sample_ID}.vcf.gz --bed expression_cis_residualised_${sample_ID}.bed.gz --permute --chunk ${chunk} 23 --seed ${seed} --threshold ${trans_threshold} --window ${trans_window} --out trans_${sample_ID}_perm.chunk_${seed}_${chunk} ${built_in_normalization}
        
    """
   
    
}




process trans_eQTLs_permutational_agg_hits {
    
    label 'medium'
    
    input:
    tuple val(sample_ID), path(chunks_trans_hits)
    val(No_peer_factors)
        
    publishDir "$params.outdir/trans_eQTLs_permutational_results_${sample_ID}_peer_${No_peer_factors}", mode: 'copy'
        
    output:
    tuple val(sample_ID), path("trans_perm_${sample_ID}_hits.tar.gz")

    when:
    params.run_trans_eqtl_analysis && params.run_trans_eqtl_permute
    
    script:
    
    """
    mkdir trans_perm_${sample_ID}_hits
    
    for i in {1..100}

    do

    cat trans_${sample_ID}_perm.chunk_\${i}_*.hits.txt.gz > trans_perm_${sample_ID}_hits/trans_${sample_ID}_perm.chunk_\${i}_hits.txt.gz
    

    done
    
    tar -zcvf trans_perm_${sample_ID}_hits.tar.gz trans_perm_${sample_ID}_hits  
        
    """
    
    
}


process trans_eQTLs_permutational_agg_bins {
    
    label 'medium'
    
    input:
    tuple val(sample_ID), path(chunks_trans_hits)
    val(No_peer_factors)
        
    publishDir "$params.outdir/trans_eQTLs_permutational_results_${sample_ID}_peer_${No_peer_factors}", mode: 'copy'
        
    output:
    tuple val(sample_ID), path("trans_perm_${sample_ID}_bins.tar.gz")

    when:
    params.run_trans_eqtl_analysis && params.run_trans_eqtl_permute
    
    script:
    
    """
    mkdir trans_perm_${sample_ID}_bins
    
    for i in {1..100}

    do

    cat trans_${sample_ID}_perm.chunk_\${i}_*.bins.txt.gz > trans_perm_${sample_ID}_bins/trans_${sample_ID}_perm.chunk_\${i}_bins.txt.gz
    

    done
    
    tar -zcvf trans_perm_${sample_ID}_bins.tar.gz trans_perm_${sample_ID}_bins   
        
    """
    
    
}


process trans_eQTLs_permutational_agg_best {
    
    label 'medium'
    
    input:
    tuple val(sample_ID), path(chunks_trans_hits)
    val(No_peer_factors)
        
    publishDir "$params.outdir/trans_eQTLs_permutational_results_${sample_ID}_peer_${No_peer_factors}", mode: 'copy'
        
    output:
    tuple val(sample_ID), path("trans_perm_${sample_ID}_best.tar.gz")

    when:
    params.run_trans_eqtl_analysis && params.run_trans_eqtl_permute
    
    script:
    
    """
    mkdir trans_perm_${sample_ID}_best
    
    for i in {1..100}

    do

    cat trans_${sample_ID}_perm.chunk_\${i}_*.best.txt.gz > trans_perm_${sample_ID}_best/trans_${sample_ID}_perm.chunk_\${i}_best.txt.gz
    

    done
    
    tar -zcvf trans_perm_${sample_ID}_best.tar.gz trans_perm_${sample_ID}_best  
        
    """
    
    
}




workflow {

    normalize_expression(Ch_unnorm_expression, params.script_normalization, params.normalization)
    merge_check(params.genotype_file, params.metadata, normalize_expression.out, params.script_merge_check, params.minor_AF)
    
    if ( params.peer_approach ) {
    
    peer_computation_2(normalize_expression.out.join(merge_check.out[1], by: 0), params.script_peer_2, params.number_of_peer_iterations)
    cis_eQTLs_indexing(peer_computation_2.out[0])
        
    cis_eQTLs_nominal(cis_eQTLs_indexing.out.join(merge_check.out[0], by: 0), Ch_trans_nom_test_500, params.cis_window, params.built_in_normalization)
    cis_eQTLs_nominal_agg(cis_eQTLs_nominal.out.groupTuple(), params.No_peer_factors)
        
    cis_eQTLs_permutational(cis_eQTLs_indexing.out.join(merge_check.out[0], by: 0), Ch_trans_nom_test_500, params.grouping, params.cis_Noperm, params.cis_window, params.built_in_normalization)
    cis_eQTLs_permutational_agg(cis_eQTLs_permutational.out.groupTuple(), params.script_cis_FDR, params.No_peer_factors)
    
    cis_eQTLs_conditional(cis_eQTLs_indexing.out.join(merge_check.out[0].join(cis_eQTLs_permutational_agg.out[0], by: 0), by: 0), params.grouping, params.No_peer_factors, params.built_in_normalization)
    
        
    } else {
    
    peer_computation(normalize_expression.out.join(merge_check.out[1], by: 0), params.script_peer, params.type_of_hidden)
    peer_test_indexing(peer_computation.out)
    
    peer_test_run_nominal(peer_test_indexing.out.join(merge_check.out[0], by: 0), Ch_peer_test_50, params.built_in_normalization)
    peer_test_run_nominal_aggregate(peer_test_run_nominal.out.groupTuple())
    
    peer_test_run_permutational(peer_test_indexing.out.join(merge_check.out[0], by: 0).combine(Ch_combined_seed_chunk_500), params.built_in_normalization)
    peer_test_run_permutational_aggregate(peer_test_run_permutational.out.groupTuple())
        
    metadata_residualize(peer_test_run_nominal_aggregate.out
                                                            .join(peer_test_run_permutational_aggregate.out
                                                            .join(normalize_expression.out
                                                            .join(merge_check.out[1].join(peer_computation.out, by: 0), by: 0), by: 0), by: 0), params.script_metaresid, params.No_peer_factors)
    
    cis_eQTLs_indexing(metadata_residualize.out[0])
    
    cis_eQTLs_nominal(cis_eQTLs_indexing.out.join(merge_check.out[0], by: 0), Ch_trans_nom_test_500, params.cis_window, params.built_in_normalization)
    cis_eQTLs_nominal_agg(cis_eQTLs_nominal.out.groupTuple(), params.No_peer_factors)
        
    cis_eQTLs_permutational(cis_eQTLs_indexing.out.join(merge_check.out[0], by: 0), Ch_trans_nom_test_500, params.grouping, params.cis_Noperm, params.cis_window, params.built_in_normalization)
    cis_eQTLs_permutational_agg(cis_eQTLs_permutational.out.groupTuple(), params.script_cis_FDR, params.No_peer_factors)
    cis_eQTLs_conditional(cis_eQTLs_indexing.out.join(merge_check.out[0].join(cis_eQTLs_permutational_agg.out[0], by: 0), by: 0), params.grouping, params.No_peer_factors, params.built_in_normalization)
    
    cis_eQTLs_for_correct(cis_eQTLs_permutational_agg.out[1].join(merge_check.out[0], by: 0), params.script_snps_extraction)
    cis_residualise(metadata_residualize.out[0].join(cis_eQTLs_for_correct.out.join(cis_eQTLs_permutational_agg.out[1], by: 0), by: 0), params.script_cisresid, params.No_peer_factors)
    trans_eQTLs_indexing(cis_residualise.out)
    
    trans_eQTLs_nominal(trans_eQTLs_indexing.out.join(merge_check.out[0], by: 0), Ch_trans_nom_test_1000, params.trans_threshold, params.trans_window, params.built_in_normalization)
    trans_eQTLs_nominal_agg_hits(trans_eQTLs_nominal.out[0].groupTuple(), params.No_peer_factors)
    trans_eQTLs_nominal_agg_fdr(trans_eQTLs_nominal_agg_hits.out.join(merge_check.out[0].join(normalize_expression.out, by: 0), by: 0), params.script_trans_FDR, params.No_peer_factors)
    trans_eQTLs_nominal_agg_bins(trans_eQTLs_nominal.out[1].groupTuple(), params.No_peer_factors)
    trans_eQTLs_nominal_agg_best(trans_eQTLs_nominal.out[2].groupTuple(), params.No_peer_factors)
    
    trans_eQTLs_permutational(trans_eQTLs_indexing.out.join(merge_check.out[0], by: 0).combine(Ch_combined_chunk_seed_1000), params.trans_threshold, params.trans_window, params.built_in_normalization)
    trans_eQTLs_permutational_agg_hits(trans_eQTLs_permutational.out[0].groupTuple(), params.No_peer_factors)
    trans_eQTLs_permutational_agg_bins(trans_eQTLs_permutational.out[1].groupTuple(), params.No_peer_factors)
    trans_eQTLs_permutational_agg_best(trans_eQTLs_permutational.out[2].groupTuple(), params.No_peer_factors)
    
    
    }
      
    
}






