# Usage

Installation

```
chmod u+x install_requirement.sh
./install_requirement.sh
```

After installing all required packages, you can run the script with following command



```
#1. make nucleotide database from BLAST
makeblastdb -in db/sequences_DNA.fasta -dbtype nucl

#2. make sure that you have put your fastq in sample folder

#3. running assembly from main script
#influenza A
./01_run.sh both 300 AGTAGAAACAAGG AGCAAAAGCAGG
#influenza B
./01_run.sh both 300 AGCAGAAGCAGAGC AGTAGTAACAAGAGC
```

