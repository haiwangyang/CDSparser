# CDSparser
fetch CDS from genome and gff

## get CDS in exons
mkdir data/<br>
cat UCSC1.final.noSEQ.gff | awk -v OFS='\t' -v FS='\t' '$2=="maker"&&$3=="CDS" {print $1, $4-1, $5, $1 ":" $4 "-" $5 ":" $9, "1", $7}' | bedtools getfasta -s -name -fi UCSC1_CLC_de_novo_rmhost_mod.fa -bed - -fo data/iCDS.fa<br>

## concatenate to make full CDS
python3 concatenate_iCDS.py

## translate
python3 translate_CDS.py



