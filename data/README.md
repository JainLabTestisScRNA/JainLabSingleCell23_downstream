## Upstream results

From preprocessing workflow

- data/cellranger/
- data/dfam_sequences.fasta*
- genes.gtf.gz

## DFAM classification

Obtained from upstream results via:

```
~/miniconda3/envs/rpm/share/RepeatMasker/famdb.py -i ~/amarel-matt/Jain/JainLabSingleCell23/results/dfam/Dfam.h5 families -f summary --class SINE --curated -ad 'mus musculus' | grep DF | awk '{print $0, "SINE"}' > data/dfam_classification.txt

~/miniconda3/envs/rpm/share/RepeatMasker/famdb.py -i ~/amarel-matt/Jain/JainLabSingleCell23/results/dfam/Dfam.h5 families -f summary --class LINE --curated -ad 'mus musculus' | grep DF | awk '{print $0, "LINE"}' >> data/dfam_classification.txt

~/miniconda3/envs/rpm/share/RepeatMasker/famdb.py -i ~/amarel-matt/Jain/JainLabSingleCell23/results/dfam/Dfam.h5 families -f summary --class LTR --curated -ad 'mus musculus' | grep DF | awk '{print $0, "LTR"}' >> data/dfam_classification.txt

~/miniconda3/envs/rpm/share/RepeatMasker/famdb.py -i ~/amarel-matt/Jain/JainLabSingleCell23/results/dfam/Dfam.h5 families -f summary --class DNA --curated -ad 'mus musculus' | grep DF | awk '{print $0, "DNA"}' >> data/dfam_classification.txt

~/miniconda3/envs/rpm/share/RepeatMasker/famdb.py -i ~/amarel-matt/Jain/JainLabSingleCell23/results/dfam/Dfam.h5 families -f summary --class RC --curated -ad 'mus musculus' | grep DF | awk '{print $0, "RC"}' >> data/dfam_classification.txt
```

This is then parsed with `workflow/scripts/resource_prep/parse_te_classification.R` and produced data/dfam_classification.parsed.txt.

## Green et al. 2018 data

green2018_supp{2,3}.xlsx from that manuscript's supplementary data.
GSE112393_MergedAdultMouseST25_12GermCellClusters_AllGeneExp.xlsx from GEO accession page processed data section.