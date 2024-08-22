1.  [Runninng demo](#org4259da5)
2.  [Plotting](#org3d2ce7b)
3.  [params.yaml](#orge628284)


<a id="org4259da5"></a>

# Runninng demo

The code in this repository is designed to run on a Linux-based machine with up-to-date versions of Nextflow and Singularity-CE.
For installation instructions, please refer to the relevant links:
 - https://github.com/nextflow-io/nextflow
 - https://sylabs.io/singularity/

The container with all the neccessary software can be build if the host machine has enought privileges with `make build_container`
command.
Alternatively the container can be downloaded with `make download_container` command.

Running a demo (assuming that singularity and nextflow with java, wget, make, git and other tools are present):
```
mkdir /tmp/demo && cd /tmp/demo
git clone <https://github.com/bioinfcollab/seqdepth>
cd seqdepth
make download_container
make test
```

<a id="org3d2ce7b"></a>

# Plotting

```
singularity exec -B $PWD containers/seqdepth.sif Rscript scripts/plotting.r --help
Uage: scripts/plotting.r [options]

Options:
        -c CRM, --crm=CRM
                Chromosome

        -p POSITION, --position=POSITION
                Chromosome position

        -m DATADIR, --datadir=DATADIR
                Directory with median files per chromosome

        -o OUTPUTDIR, --outputdir=OUTPUTDIR
                Directory for the result plot

        -t TYPE, --type=TYPE
                DataType - median or hypergeom

        -h, --help
                Show this help message and exit
```
Example:

```
mkdir /tmp/plot
singularity exec -B $PWD containers/seqdepth.sif Rscript scripts/plotting.r --crm 19 --position 51417359 --datadir "results/hypergeometric/" -t hypergeom --outputdir /tmp/plot
```


<a id="orge628284"></a>

# params.yaml

```
# yaml file, use a text editor!
# validate: perl -MYAML -e 'use YAML;YAML::LoadFile("./params.yaml")'

# datadir: path to bam or bed files, in case when data files are distributed
# over directory structure, use symlinks to create a directory pointing to
# specific set of files from a dataset
datadir:                /path_to/bamfiles
# outdir: directory with the  results
outdir:                 /path_to_results
# seqdepth singularity container with scripts and tools
container:              ./containers/seqdepth.sif
# singularity container options, for example bind mounts '-B /lustre'
containerOptions:
# temporary directory, must be big enough for the results, can be also set
# via TMPDIR env variable
temp:                   '/tmp/test'
# Defines on how many parts a single bed file will be split,
# affects pipeline concurency, since every part will be processed
# as a separate task
# the size of the part is determined as 250000000/split
# thus the maxpossiblepostitions must be divisible on the split number
split:                  5
# number of processes for per task
nproc:                  8
# If defined only a specific chromosome will be extracted and processed
# if not set, the range 1..22, X, Y is assumed
chromosome:             6
# some bam files come with chr prefix, other without,
# check with samtools idxstat what is the correct prefix for the
# current dataset
# current pipeline code is aware only of 'chr' prefix
chrom_prefix:
# optinal, needed for a slurm cluster
hpcproject:	<hpcProj>
# by default only median, second option hypergeometric
task:
## hypergeometric calcucation parameters
#minDepth:
#cureentTotalDepth:
#minAcceptedProb:
#tocheck:
```
