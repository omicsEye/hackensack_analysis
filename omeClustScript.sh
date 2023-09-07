#!/usr/bin/env bash
SAVEIFS=$IFS
IFS="
"
list=("2'-O-ribose methyltransferase"
"3'-to-5' exonuclease"
"3C-like proteinase"
"E"
"Genome"
"M"
"N"
"ORF10"
"ORF1a"
"ORF1ab"
"ORF1b"
"ORF3a"
"ORF6"
"ORF7a"
"ORF7b"
"ORF8"
"RNA-dependent RNA polymerase"
"S"
"endoRNAse"
"helicase"
"leader protein"
"nsp10"
"nsp11"
"nsp2"
"nsp3"
"nsp4"
"nsp6"
"nsp7"
"nsp8"
"nsp9")
for item in ${list[@]}
  do
    echo $item
    omeClust -i ../omeClust/$item/adist.txt --metadata cleanMetadataNonDiscretized.tsv -o ./"$item"_output --plot
  done
IFS=$SAVEIFS
