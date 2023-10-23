function samtoolsStat(){
    echo "Usage: samtoolsStat *.bam"
    local bamfile=$1
    local thread=${2:-10}
    if [ -z $bamfile ];then
        echo "Please provide bamfile,this function include sort,index,flagstats command."
        echo "Usage: samtoolsStat <*.bam>"
    else

        local sortbamOpt="$bamfile.sort"
        local bamStatOpt="$bamfile.idxstatstat"
        local bamflagStat="$bamfile.flagstat"
        local bamcoverage="$bamfile.coverage"
        local bamdepth="$bamfile.depth"
        #bowtie2 -x $ref -1 $seq1 -2 $seq2 --fast-local -p 10| samtools view -bS - >$bamOpt
        samtools sort -@ $thread -o $sortbamOpt $bamfile
        samtools index $sortbamOpt
        echo -ne "reference_sequence_name\tsequence_length\t#mapped_read-segments\t#unmapped_read-segments\n" >$bamStatOpt
        samtools idxstats $sortbamOpt >>$bamStatOpt
        samtools flagstat $sortbamOpt >$bamflagStat
        #samtools coverage $sortbamOpt >$bamcoverage
        #samtools bedcov $sortbamOpt >$bamcoverage
        samtools depth -a $sortbamOpt > $bamdepth
    fi
}

samtoolsStat $1
