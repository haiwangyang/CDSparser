# CDSparser
fetch CDS from genome and gff

## get genome (fasta)
put data in CDSparser/data/<br>

## get annotation (gff; and remove sequence at the end if neccessary)
put data in CDSparser/data/<br>

## get CDS in exons
mkdir data/<br>
cat UCSC1.final.noSEQ.gff | awk -v OFS='\t' -v FS='\t' '$2=="maker"&&$3=="CDS" {print $1, $4-1, $5, $1 ":" $4 "-" $5 ";" $9, "1", $7}' | bedtools getfasta -s -name -fi UCSC1_CLC_de_novo_rmhost_mod.fa -bed - -fo data/iCDS.fa<br>

## concatenate to make full CDS
python3 concatenate_iCDS.py<br><br>

(1) it conconcate iCDS to make CDS<br>
(2) it checks if CDS ends with stop codon (TAA, TGA, TAG)<br>
(3) if check if CDS is multiple of three<br>
(4) if translate the CDS into PEP<br>
(5) if output all the above results in CDSparser/output/ folder<br><br>


