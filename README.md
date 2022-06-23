## Team mitten_TDC19
[DREAM Tumor Deconvolution Challenge](https://www.synapse.org/#!Synapse:syn15589870/wiki/582446), team mitten_TDC19 from Michigan State University, East Lansing, Michigan: 

+ Dr. Sergii Domanskyi (domansk6@msu.edu)
+ Thomas Bertus (bertusth@msu.edu)
+ Prof. Carlo Piermarocchi (piermaro@msu.edu)

## Method description

This repository contains source code from [synapse.org: mitten_TDC19](https://www.synapse.org/#!Synapse:syn22175932) and documentation on how to run the scripts.

The method writeup can be found on synapse [Wiki page](https://www.synapse.org/#!Synapse:syn20330693/wiki/604067) of our team.

Here we provide pre-trained models in form of the signature matrices for Submissions 1, 2, 3 of the Final Round of the challenge for both coarse-grained and fine-grained subchallenges. The spreadsheet files are in `signature` directory of this repository. The name of each of the 6 signature matrices contains subchallenge identifier and final round number:

+ File `coarse_final_round_2.xlsx` contains all markers from `coarse_final_round_1.xlsx` plus 63 additional markers for CD4.T.cells and 5 additional markers for CD8.T.cells. File `coarse_final_round_3.xlsx` is identical to `coarse_final_round_1.xlsx`.

+ File `fine_final_round_2.xlsx` contains all markers from `fine_final_round_1.xlsx` plus 150 additional markers for myeloid.dendritic.cells and 9 additional markers for memory.CD4.T.cells. File `fine_final_round_3.xlsx` contains all markers from `fine_final_round_1.xlsx` plus 5 additional markers for myeloid.dendritic.cells and 9 additional markers for memory.CD4.T.cells. 


## Installation

Pre-Installation Requirement: to use this software on any platform you need Python version 3.7 or higher. Install packages `numpy`, `pandas`, `scipy`, `xlrd`, `matplotlib`, `mygene`, `sklearn` used in the scripts.

To clone this repository run:

```
git clone https://github.com/sdomanskyi/mitten_TDC19.git
```

## Typical install time
Installation of the packages `numpy`, `pandas`, `scipy`, `xlrd`, `matplotlib`, `mygene`, `sklearn` on a normal desktop computer takes about 5 minutes.
Installation of this git repository takes a few seconds.

## Usage

Run ```python run_model.py --help``` to see description of all command line options.

> Note: This repository contains source code that was used in Submissions 1, 2, 3 of the Final Round of the challenge, as well as the Leaderboard Rounds. The parameters exposed in the script `run_model.py` constitute the minimal set to reproduce the Final Rounds submissions and process any new user-provided data.

Input file `input.csv` should contain at least 3 columns as described below. Every row in the input file corresponds to a dataset:

| Column         | Description  |
|----------------|--------------|
| dataset.name   | Unique identifier of the dataset |
| scale | Descriptor of how the gene expression data was scaled, e.g. linear, log2 |
| hugo.expr.file  | Name of the file with normalized gene expression and gene symbol gene identifiers  |


Coarse-grained model run for Final Round 1:
```
python run_model.py --input=input/input.csv --model-level=coarse --signature-matrix=signature/coarse_final_round_1.xlsx
```

Fine-grained model run for Final Round 1:
```
python run_model.py --input=input/input.csv --model-level=fine --signature-matrix=signature/fine_final_round_1.xlsx
```

Coarse-grained model run for Final Round 2:
```
python run_model.py --input=input/input.csv --model-level=coarse --signature-matrix=signature/coarse_final_round_2.xlsx
```

Fine-grained model run for Final Round 2:
```
python run_model.py --input=input/input.csv --model-level=fine --signature-matrix=signature/fine_final_round_2.xlsx
```

Coarse-grained model run for Final Round 3:
```
python run_model.py --number-hvg=15 --input=input/input.csv --model-level=coarse --signature-matrix=signature/coarse_final_round_3.xlsx
```

Fine-grained model run for Final Round 3:
```
python run_model.py --number-hvg=15 --input=input/input.csv --model-level=fine --signature-matrix=signature/fine_final_round_3.xlsx
```

> Note: In Final Round 3 we limited the number of highly variable genes to 15 by passing parameter `--number-hvg=15`.


Output file(s) is generated in the directory `output`. Below is the output file from the coarse-grained model run. Prediction is relative abundance score of a given sample for a given cell type.

|dataset.name|sample.id|cell.type        |prediction |
|------------|---------|-----------------|-----------|
|FIAS5       |S1       |B.cells          |2.225265108|
|FIAS5       |S1       |NK.cells         |1.465895084|
|FIAS5       |S1       |neutrophils      |2.135799171|
|FIAS5       |S1       |monocytic.lineage|2.449364958|
|FIAS5       |S1       |CD4.T.cells      |0.606665336|
|FIAS5       |S1       |CD8.T.cells      |0.305858142|
|FIAS5       |S1       |fibroblasts      |0.657552838|
|FIAS5       |S1       |endothelial.cells|0.80516959 |
|FIAS5       |S2       |B.cells          |2.200313867|
| ~trimmed~  |

## Expected run time
It took approximately 1 minute to run each of the 6 scenarios detailed above.

## Licensing
This software is released under an MIT License. Please also consult the file LICENSE in this repository regarding Licensing information for use of external associated content.

## Funding

This work was supported by the National Institutes of Health, Grant No. R01GM122085. 

## Acknowledgements

We acknowledge discussions with Prof. Giovanni Paternostro and Alex Hakansson. 

