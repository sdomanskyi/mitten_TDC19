## Team mitten_TDC19
[DREAM Tumor Deconvolution Challenge](https://www.synapse.org/#!Synapse:syn15589870/wiki/582446), team mitten_TDC19 from Michigan State University, East Lansing, Michigan: 

+ Dr. Sergii Domanskyi (domansk6@msu.edu)
+ Thomas Bertus (bertusth@msu.edu)
+ Prof. Carlo Piermarocchi (piermaro@msu.edu)


This repository contains source code from [synapse.org: mitten_TDC19](https://www.synapse.org/#!Synapse:syn22175932) and documentation on how to run the scripts.

The method writeup can be found on synapse [Wiki page](https://www.synapse.org/#!Synapse:syn20330693/wiki/604067) of our team.

## Installation

Pre-Installation Requirement: to use this software on any platform you need Python version 3.7 or higher. Also install packages `numpy`, `pandas`, `scipy`, `xlrd`, `matplotlib`, `mygene`, `sklearn` used in the scripts.

To clone this repository run:

```
git clone https://github.com/sdomanskyi/mitten_TDC19.git
```

## Usage

Run ```python run_model.py --help``` to see description of all command line options.

> Note: This repository contains source code that was used in Rounds 1, 2, 3 and Final Round of the chalenge. The parameters exposed in the script `run_model.py` constitute the minimal set to process new data. The details of each run can be found at [synapse.org: mitten_TDC19](https://www.synapse.org/#!Synapse:syn22175932).

Input file `input.csv` should contain at least 3 columns as described below. Every row in the input file corresponds to a dataset:

| Column         | Description  |
|----------------|--------------|
| dataset.name   | Unique identifier of the dataset |
| scale | Descriptor of how the gene expression data was scaled, e.g. linear, log2 |
| hugo.expr.file  | Name of the file with normalized gene expression and gene symbol gene identifiers  |

Example of coarse-grained model run:

```
python run_model.py --input=input/input.csv --model-level=coarse --signature-matrix=signature/signature_matrix_combo_coarse.xlsx
```

Example of fine-grained model run:

```
python run_model.py --input=input/input.csv --model-level=fine --signature-matrix=signature/signature_matrix_combo_fine.xlsx
```

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


## Licensing
This software is released under an MIT License. Please also consult the file LICENSE in this repository regarding Licensing information for use of external associated content.

## Funding

This work was supported by the National Institutes of Health, Grant No. R01GM122085. 

## Acknowledgements

We acknowledge discussions with Prof. Giovanni Paternostro and Alex Hakansson. 

