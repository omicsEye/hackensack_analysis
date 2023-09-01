mkdir results_"$filename"
mkdir results_"$filename"/pathomap-Pathoscope
mkdir results_"$filename"/pathomap-Pathoscope/bacterialRefBowtie2
mkdir results_"$filename"/pathomap-Pathoscope/humanRef
mkdir results_"$filename"/pathomap-Pathoscope/viralRef
cutadapt -a AGATCGGAAGAG -A AGATCGGAAGAG -o tr_"$filename"_R1.fastq -p tr_"$filename"_R2.fastq AllFastqFiles/"$filename"_R1_cat.fastq AllFastqFiles/"$filename"_R2_cat.fastq
pathoscope MAP -1 tr_"$filename"_R1.fastq -2 tr_"$filename"_R2.fastq \
-numThreads 40 \
-indexDir . \
-targetIndexPrefixes Pathoscope/bacterialRefBowtie2/allBacGenomes \
-filterIndexPrefixes Pathoscope/humanRef/GCF_000001405.39_GRCh38.p13_genomic.fna_bt2,Pathoscope/viralRef/allViralRefSeqGenomes.fna_bt2 \
-outDir results_"$filename" \
-outAlign "$filename".sam
rm results_"$filename"/pathomap-Pathoscope/bacterialRefBowtie2/*.sam
rm results_"$filename"/pathomap-appendAlign.sam
pathoscope ID -alignFile results_"$filename"/"$filename".sam -fileType sam -outDir results_"$filename" -expTag "$filename"
rm results_"$filename"/"$filename".sam
