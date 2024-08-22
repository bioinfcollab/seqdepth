#!/usr/bin/env nextflow

median_per_chr           = params.median_per_chr
chrom_prefix             = params.chrom_prefix
nproc                    = params.nproc
medianNproc              = 2*nproc

/*
Check that chr prefix matches

so far it is assumed that chrom_prefixs in bam files either start with chr
or are just plain numbers with X and Y for 23 and 24
see chr_ids channel

*/
process bamFileChromPrefix_test {
    label 'low_res'
    input:
    errorStrategy 'terminate'
    tuple val(sampleID), file(sampleFile)
    val(chromosome)

    output:
    val(sampleID)

    script:
    """
    ${projectDir}/scripts/check_chr.sh ${sampleID} ${chrom_prefix}${chromosome}
    """

}
/* Convert bam files to bed files compressed with bgzip
Bgzip  compresses files in a similar manner to, and compatible with, gzip(1).
The file is compressed into a series of small (less than 64K) 'BGZF' blocks.
This allows indexes to be built against the compressed file and used to retrieve
portions of the data without having to decompress the entire file.
Build tabix indexes (.tbi)
See bam header with for example samtools view -H <path_to_bam_file>
*/
process bam2bedGenomecov {
    label 'bam2bed'
    input:
    tuple val(sampleID), file(sampleFile)
    val(chr_test_out)

    output:
    path("${sampleID}.*bed.gz"), emit: bedfiles_ch
    path("${sampleID}.*tbi"), emit: bedfiles_tbi

    script:
    if (! chrom_prefix )
    """
    bedtools genomecov -ibam ${sampleID}.bam -dz | bgzip > ${sampleID}.bed.gz
    tabix -b2 -e2 ${sampleID}.bed.gz
    """
    else if ( chrom_prefix  == 'chr' )
    """
    bedtools genomecov -ibam ${sampleID}.bam -dz | cut -b4- | bgzip > ${sampleID}.bed.gz
    tabix -b2 -e2 ${sampleID}.bed.gz
    """
    else
    error "Not recognized chromosome prefix"
    /*
    """
    bedtools genomecov -ibam ${sampleID}.bam -dz | sed s"^${chrom_prefix}//" | bgzip > ${sampleID}.bed.gz
    tabix -b2 -e2 ${sampleID}.bed.gz
    """
    */

}

/* Split the bed files per chromosome and save them in the directory structure allowing
furher processing with SeqDepth scripts
Split chromosomes to subsets each of $step
*/
process extractChromosomes {
    label 'extractChromosomes'
    cache false
    input:
    file(bedfiles)
    file(tbis)
    val(step)
    each chr_id
    each subid

    output:
    path("${chr_id}/*"), emit: chr_by_sample
    tuple val(chr_id),val(subid), emit: chr_ids_processed

    script:
    """
    step=${step}
    from="\$(( ${subid}*\${step}+1 ))"
    to="\$(( ${subid}*\${step}+\${step} ))"
    mkdir "${chr_id}/${subid}" -p;
    ls -1 *.bed.gz | parallel -P $nproc "tabix --verbosity=2 -b2 -e2 {} ${chr_id}:\${from}-\${to} | bgzip > ${chr_id}/${subid}/{/}"
    """
}

/*Build table of reads per position over all samples, omit zero values to save space.*/
process buildVectorPerPosition {
    label 'mid_res'
    // iterate over chr_ids to split bed files per chromosome with tabix
    input:
    tuple val (chr_id),val(subid)
    path(chr_by_sample)

    output:
    tuple val(chr_id), val(subid), path("${chr_id}/${subid}/result_table.txt.gz")

    script:
    """
    mkdir ${chr_id}/${subid} -p && ${projectDir}/scripts/samples_table.pl -n $nproc -v -s ./${chr_by_sample} -n ${nproc} -t ${params.temp} -d ./${chr_id}/${subid}
    """
}

/* Calculate median for the tables generated in buildVectorPerPositon process */
process medianPerVector {
    label 'mid_res'
    input:
    tuple val(chr_id), val(subid), path(table_file)
    val(fnumber)

    publishDir median_per_chr, mode: params.mode, overwrite: true

    output:
    tuple val(chr_id), path("${chr_id}/result_median*.txt.gz")

    script:
    """
    mkdir -p ${chr_id}/${subid}
    ${projectDir}/scripts/median.pl -i "${table_file}" -N ${fnumber} -n ${medianNproc} -d ${chr_id}/${subid}
    mv ${chr_id}/${subid}/result_median.txt.gz ${chr_id}/result_median_${subid}.txt.gz
    rm -vrf ${chr_id}/${subid}
    """
}

process joinResultMedianFiles{
    memory '8 GB'
    cpus 4
    input:
    tuple val(chr_id), path(result)

    output:
    tuple val(chr_id), path("${chr_id}_median.txt.gz")

    publishDir median_per_chr, mode: params.mode, overwrite: true

    //${result} contains mulitple files! if the number of files is too long, cat may not work
    script:
    """
    cat ${result} >> merged.txt.gz
    zcat merged.txt.gz | sort -nk1 | pigz -p 4 >> ${chr_id}_median.txt.gz
    """
}

/* Run hypergeometric.pl
Output:  Chromosome, Position, Pvalue? Expect_depth Prob_expect_depth Prob_number_tock
*/
process hyperGeometric{
    memory '2 GB'
    cpus 4
    input:
    tuple val(chr_id), path(median_file)

    output:
    tuple val(chr_id), path("${chr_id}")

    script:
    """
    mkdir -p ${chr_id}; ${projectDir}/scripts/hypergeometric.pl -s "${median_file}" -C ${chr_id} -n $nproc -c ${params.tocheck} \
    --minap ${params.minAcceptedProb} --ctd ${params.currentTotalDepth} --mindepth ${params.minDepth} -r ${chr_id}/hypergeometric_unsorted.txt
    """
}

process hyperGeometricSort{
    memory '3 GB'
    cpus 4
    input:
    tuple val(chr_id), path(my_chr_id)
    val(hypergeom_publish_dir)

    output:
    path("${chr_id}")

    publishDir "${hypergeom_publish_dir}", mode: params.mode, overwrite: true

    script:
    """
    zcat ${chr_id}/hypergeometric_unsorted.txt.gz | sort -nk 2 | pigz -p 4 > ${chr_id}/hypergeometric.txt.gz
    #rm -v ${chr_id}/hypergeometric_unsorted.txt.gz
    """
}
