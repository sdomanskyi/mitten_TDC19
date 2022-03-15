FROM ubuntu:18.04
FROM python:3.7.3-stretch

COPY output/ /output/
COPY tools/ /tools/
COPY requirements.txt /requirements.txt
COPY Sub_Challenge_2/signature_matrix_combo_fine.xlsx /signature_matrix_combo_fine.xlsx
COPY Sub_Challenge_1/signature_matrix_combo_coarse.xlsx /signature_matrix_combo_coarse.xlsx

COPY gencode.v31.annotation.gtf.len.ensembl.hugo.pklz /gencode.v31.annotation.gtf.len.ensembl.hugo.pklz

COPY run_model.py /run_model.py
COPY tdc_models.py /tdc_models.py
COPY tdc_metrics.py /tdc_metrics.py

RUN chmod a+x run_model.py

RUN pip install --upgrade pip
RUN pip install -r requirements.txt

ENTRYPOINT ["python","run_model.py"]
