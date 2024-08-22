#!/usr/bin/env nextflow
nextflow.enable.dsl = 2
//vim: set filetype=groovy

if (params.help) {
    log.info ''
    log.info 'S E Q D E P T H'
    log.info '---------------------------------'
    log.info 'Calucalte median and hypergeometric probabylity over all reads.'
    log.info ''
    log.info 'Usage: '
    log.info '   nextflow run seqdepth.nf [OPTIONS]'
    log.info ''
    log.info 'Options: '
    log.info '    --help                                    Show this message and exit.'
    log.info '    --datadir   <dir>                         Directory with bam files, possibly spread over multiple subdirectories.'
    log.info '    --outdir    <dir|workflow_results>        Optional result directory, by default $baseDir/workflow_results'
    log.info '    --temp      <dir|/tmp>                    Temporary directory, by default /tmp, ignores TMPDIR env'
    log.info '    --container <path|$baseDir/seqdepth>      Singularity container location, by default $baseDir/seqdepth'
    log.info '    --chromosome chromosome|range             Specify a chromosome or chromosome range, see channel definition'
    log.info '    --step  <10|25|50>                        Split every chromosome in a $step of subranges and process them as different tasks'
    log.info '    --hpcproject <same as -A for srun>        Specify project when runnin on hpc, the default is hardcoded in nextflow.config'
    log.info '    --task  <hypergeometric or median>         Hypergeomtric or median only, hypergeometric includes median which is the default'
    log.info '    --hypergeom parameters                    TODO'
    log.info '      --minDepth                              <todo>   '
    log.info '      --currentTotalDepth                     <todo>   '
    log.info '      --minAcceptedProb                       <todo>   '
    log.info '      --tocheck                               <todo>   '

    log.info ''
    log.info 'Nextflow options: '
    log.info '     -profile <local|slurm>       Local for the local execution, SLURM for ZIH/Biocluster4'
    log.info '     -resume                      Resume failed/interruped run, use cached results'
    log.info ''
    exit 1
}

/* some rudimental param checking
*/
if (! params.datadir) {
    exit 1, "Missing --datadir parameter"
}
if (! file(params.datadir,type: 'dir').isDirectory()) {
    exit 1, "Directory ${params.datadir} doesn't exist"
}

// Input validation

tempdir = file(params.temp, type: 'dir')
if (! tempdir.isDirectory()) {
    exit 1, "Directory ${tempdir} doesn't exist"
}

container = file(params.container)
if (! container.isFile()) {
    exit 1, "Can't find the container: ${container}, check the path or specify --container option"
}

// Params

task                                    = params.task
outdir                                  = file(params.outdir, type: 'dir')
datadir                                 = file(params.datadir)
input_bam_files                         = params.datadir + '/**.bam'
input_bam_index_files                   = params.datadir + '/**.bai'
bam_files_list                          = file(input_bam_files)
hypergeom_publish_dir                   = file(params.outdir + '/hypergeometric', type: 'dir')

input_bed_files                         = params.datadir + '/**.bed.gz'
input_bed_index_files                   = params.datadir + '/**.bed.gz.tbi'
bed_files_list                          = file(input_bed_files)
tabix_index_files_list                  = file(input_bed_index_files)

nproc                   = params.nproc
//result_bed_per_chr      = file(outdir + '/bedfiles_per_chr_temp', type: 'dir')
//table_dir               = file(outdir + '/table', type: 'dir')
split                   = params.split
step                    = splitStep(split)
chrom_prefix            = params.chrom_prefix;

if (bed_files_list) {
    fnumber     = bed_files_list.size()
}
else if (bam_files_list) {
    fnumber     =  bam_files_list.size()
}

def splitStep(split) {
        maxposition = 250000000
        defined_split = split.toInteger()
        if (maxposition.mod(defined_split) != 0) {
                exit 1, "${maxposition} can't be divided on the split value: ${split}"
        }
        result_s = maxposition/defined_split
        return result_s
}

def checkInputFiles(index_files,data_files,datadir) {
        if (! file(index_files) ) {
                exit 1, "Can't find index files ${index_files}"
        }
        if ( file(index_files).size() != file(data_files).size() ) {
                exit 1, "\n\033[1mNumber of the index files don't the match the number of data files in the datadir: ${datadir}\033[0m\n"
        }
}

include { bamFileChromPrefix_test; bam2bedGenomecov; extractChromosomes; buildVectorPerPosition; medianPerVector;joinResultMedianFiles; hyperGeometric; hyperGeometricSort} from './modules/seqdepth.nf'

workflow {
    outdir.mkdir()

    def chrom_prefix_msg = params.chrom_prefix ?: 'Default, no prefix'
    chromosome          = chrom_prefix ? chrom_prefix+params.chromosome : params.chromosome

    /*  ======================================================================================================
    *  RUN INFO
    *  ======================================================================================================
    */
    log.info "===================================="
    log.info 'S E Q D E P T H pipeline         '
    log.info "===================================="
    log.info "Container                 : ${params.container}"
    log.info "Input data                : ${params.datadir}"
    log.info "Output path               : ${outdir}"
    log.info "Temporary dir             : ${params.temp}"
    log.info "Nproc                     : ${nproc}"
    log.info "Fnumber                   : ${fnumber}"
    log.info "Split                     : ${split}"
    log.info "Chromosome                : ${chromosome}"
    log.info "HPC Project               : ${params.hpcproject}"
    log.info "Chrom prefix:             : ${chrom_prefix_msg}"
    log.info "Task                      : ${task}"
    if ( task == 'hypergeometric') {
            log.info "minDepth:                 : ${params.minDepth}"
            log.info "currentTotalDepth:        : ${params.currentTotalDepth}"
            log.info "minAcceptedProb:          : ${params.minAcceptedProb}"
            log.info "tocheck:                  : ${params.tocheck}"
    }
    log.info "====================================\n"


    // Chromosome ids, list from 1 to 22 plus X and Y
    if (! params.chromosome) {
        chr_ids = channel.of(1..22,'X','Y')
        //        if (chrom_prefix) {
        //            chr_ids = chr_ids.map{chrom_prefix + it}
        //        }
    }
    else {
        chr_ids = channel.of(params.chromosome)
    }
    chr_subids = channel.of(0 .. split - 1)

    initial_bed_files_list = bed_files_list
    if (bam_files_list) {
        checkInputFiles(input_bam_index_files,input_bam_files,datadir)
        //bamfiles = Channel.fromPath(input_bam_files).map { file_path -> tuple(file_path.baseName, file_path, file(params.datadir).relativize(file_path).getParent()) }
        // check https://nextflow-io.github.io/patterns/index.html example for bam & bai
        bamfiles = Channel.fromFilePairs(params.datadir + '/*.{bam,bai}'){ file -> file.name.replaceAll(/.bam|.bai$/,'') }
        //bamfiles.ifEmpty("\n\033[1mCan't find corresponding index (*.bai) files in [ ${params.datadir} ], exitting\033[0m\n").view()
        //first_chr_id = chrom_prefix ? chrom_prefix+chr_ids.first() : chr_ids.first()
        //first_chr_id.view()
        bamFileChromPrefix_test(bamfiles.first(),chr_ids.first())
        bam2bedGenomecov(bamfiles,bamFileChromPrefix_test.out)
        bed_files_list = bam2bedGenomecov.out.bedfiles_ch.collect()
        tabix_index_files_list = bam2bedGenomecov.out.bedfiles_tbi.collect()
    }
    if (initial_bed_files_list) {
        checkInputFiles(input_bed_index_files,input_bed_files,datadir)
    }
    extractChromosomes(bed_files_list,tabix_index_files_list,step,chr_ids,chr_subids)
    chr_by_sample = extractChromosomes.out.chr_by_sample
    chr_ids_processed = extractChromosomes.out.chr_ids_processed
    buildVectorPerPosition(chr_ids_processed,chr_by_sample)
    medianPerVector(buildVectorPerPosition.out,fnumber)
    collected_ch = medianPerVector.out.groupTuple()
    joinResultMedianFiles(collected_ch)
    if (task == 'hypergeometric') {
        hyperGeometric(joinResultMedianFiles.out)
        hyperGeometricSort(hyperGeometric.out,hypergeom_publish_dir)
    }

}

workflow.onComplete = {
    println "Pipeline complete"
    println "Command line: $workflow.commandLine"
}
