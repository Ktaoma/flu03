shopt -s expand_aliases
alias samtools=$(find "tools/samtools-1.22" -name "samtools")
makeblastdb -in db/sequences_DNA.fasta -dbtype nucl

#flu A primer
## fp="AGTAGAAACAAGG", "AGCAAAAGCAGG"
#flu B primer
## fp="AGCAGAAGCAGAGC", "AGTAGTAACAAGAGC"


fl=$(ls sample/barcode*)
type=both
leng=350
fp=$1
rp=$2

for i in $fl;
do
    #make sure that your fastq is .fastq.gz
    sn=$(echo $i | cut -d '/' -f2 | sed -e 's/.fastq.gz//g')
    echo $sn

    path=result/$sn/${type}
    mkdir -p $path

    echo 'Get summary statistic of uncurated file'
    seqkit stat $i -Ta | sed -e 's/,//g' >> ${path}/read_count.tsv

    echo 'Get primer pattern'
    pp=$path/pattern/$type
    mkdir -p $pp
    seqkit locate -m 2 -p ${fp} $i | awk -v ss=$sn '{print $0,ss}' | gzip -9 > ${pp%/type}/reverse.tsv.gz
    seqkit locate -m 2 -p ${rp}  $i | awk -v ss=$sn '{print $0,ss}' | gzip -9 > ${pp%/type}/forward.tsv.gz

    echo 'Extract primer'
    mkdir -p ${pp}
    Rscript ./02_parse_primer_pattern.R $pp $type

    echo 'Export as Fastq'
    seqkit grep -f $pp/QC_id.tsv $i | gzip -9 > $path/filtered.fq.gz

    echo 'Trim primer with error rate at 0.25 and quality 10'
    cutadapt --revcomp -a ${fp}...${rp} -e 0.25 -m $leng -q 10 -o $path/sample_final_qc.fq.gz $path/filtered.fq.gz
    seqkit fq2fa $path/sample_final_qc.fq.gz | gzip -9 > $path/sample_final.fa.gz

    echo 'Sort reads into bucket'
    zcat $path/sample_final.fa.gz | blastn -query - -db db/sequences_DNA.fasta -outfmt '6 qseqid sseqid pident length qlen' > $path/blast_hits.tsv
    cat ${path}/blast_hits.tsv | awk '{prop=100*($4/$5);print $0 "\t" prop}'| awk '$3 >= 90 && $6 >= 90' > ${path}/blast_hits_filtered.tsv

    echo 'Get summary statistics of QC fasq'
    seqkit stat $path/sample_final.fa.gz -Ta | sed -e 's/,//g' | awk 'NR!=1 {print $0}' >> ${path}/read_count.tsv

    for k in 1 2 3 4 5 6 7 8;
    do
    
        echo 'Get read of ${k} segment'
        cat db/sequences_DNA.fasta | grep ">" | grep "segment $k"  |\
            awk '{print $1}' | sed 's/>//g' | grep -f - $path/blast_hits_filtered.tsv | cut -f1 | sort -u |\
                seqkit grep -f - $path/sample_final.fa.gz | gzip -9 > $path/sample_final_segment_${k}.fa.gz

        echo 'Get summary statistics of k segment'
        seqkit stat $path/sample_final_segment_${k}.fa.gz -Ta | sed -e 's/,//g' | awk 'NR!=1 {print $0}' >> ${path}/read_count.tsv

        echo 'Get the name of best reference of ${k} segment'
        name=$(cat db/sequences_DNA.fasta | grep ">" | grep "segment $k"  | awk '{print $1}' |\
               sed 's/>//g' | grep -f - ${path}/blast_hits_filtered.tsv | cut -f2 | sort | uniq -c | sort -k1,1nr | awk 'NR==1{print $2}')

        echo 'Remove supplmentary reads'
        seqkit fx2tab db/sequences_DNA.fasta | grep $name | seqkit tab2fx |\
            minimap2 -a - $path/sample_final_segment_${k}.fa.gz |\
                samtools sort -o $path/sample_final_segment_${k}_raw.bam && samtools index $path/sample_final_segment_${k}_raw.bam

        echo 'Select mapped read and remove supplementary read'
        samtools view -b -F 4 $path/sample_final_segment_${k}_raw.bam > $path/sample_final_segment_${k}_tmp.bam
        samtools view -f 2048  $path/sample_final_segment_${k}_tmp.bam | cut -f1 | sort -u > $path/supplement_${k}.txt
        samtools view -h $path/sample_final_segment_${k}_tmp.bam | grep -vf $path/supplement_${k}.txt | samtools view -b -o $path/sample_final_segment_${k}.bam -

        echo 'Remove read with excessive INDEL and softclip'
        samtools view -h $path/sample_final_segment_${k}.bam | awk '{print $1,$6}' | grep -v ^@ > ${path}/ReadCIGAR.tsv

        Rscript ./03_parse_CIGAR.R $path
        mv ${path}/ReadCIGAR.tsv ${path}/ReadCIGAR_${k}.tsv
        mv ${path}/CIGAR_Tab.tsv ${path}/CIGAR_Tab_${k}.tsv

        cat ${path}/CIGAR_Tab_${k}.tsv  | awk '$10 != "keep" {print $1}' > $path/sample_final_segment_${k}_ID_filtered.txt
        samtools view -h $path/sample_final_segment_${k}.bam | grep -vf $path/sample_final_segment_${k}_ID_filtered.txt | samtools view -b -o $path/sample_final_segment_${k}_filtered.bam -
        samtools sort $path/sample_final_segment_${k}_filtered.bam -o $path/sample_final_segment_${k}_filtered_sorted.bam
        samtools index $path/sample_final_segment_${k}_filtered_sorted.bam

        echo 'Reference based assembly'
        samtools consensus -m simple -d 30 -A $path/sample_final_segment_${k}_filtered_sorted.bam |\
               seqkit fx2tab | awk -v sample=$sn -v seg=$k 'BEGIN{OFS="\t"}{print sample "_" seg,$2}' |\
                        seqkit tab2fx > $path/draft_segment_${k}.fasta

        echo 'Calculate depth'
        samtools depth $path/sample_final_segment_${k}_filtered_sorted.bam |\
               awk -v sample=$sn -v seg=$k '{print sample "_" seg,$2,$3}' > $path/depth_segment_${k}.tsv

        echo 'Self alignment'
        samtools fastq $path/sample_final_segment_${k}_filtered_sorted.bam | gzip -9 > $path/sample_final_segment_${k}_filtered_sorted.fq.gz
        minimap2 -a $path/draft_segment_${k}.fasta $path/sample_final_segment_${k}_filtered_sorted.fq.gz > $path/self_aln_${k}.sam
        samtools sort $path/self_aln_${k}.sam -o $path/self_aln_${k}.bam


    done
done
