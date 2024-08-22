#! /bin/bash
set -e

sample=$1
chrom_prefix=$2


chr_prefix=$(head -1 <(samtools idxstat ${sample}.bam ) | cut -f1 |cut -b -3)
# silly check that both prefixes are 1,2 char size - i.e digits
# more precise checking shall be done in a proper language like perl/python
#if [[ $(expr length "${chr_prefix}${chrom_prefix}") -le 3 ]]
#then
#	exit 0
#fi
# undefined --chrom_prefix is equal to 'null' in this case
if [[ "$chrom_prefix" == 'null'*  ]]
then
	if [[ "$chr_prefix"  == 'chr' ]]; then
		echo -e "\033[1m----->>>>>    Emtpy --chrom_prefix while detected chr prefix is equal to 'chr'    <<<<<-----\033[0m"
		exit 2
	fi
exit 0
fi
if [[ ! -z $chrom_prefix  ]]
then
	if [[ "$chrom_prefix"  != "$chr_prefix"* ]]; then
		echo -e "\033[1m----->>>>>    Detected (via samtools idxstat) chr_prefix($chr_prefix) and provided --chrom_prefix($chrom_prefix) don't match    <<<<<-----\033[0m"
		exit 3
	fi
fi
