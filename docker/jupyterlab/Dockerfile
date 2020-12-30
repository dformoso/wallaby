FROM jupyter/r-notebook:latest

COPY install_libraries.R .
RUN Rscript install_libraries.R

RUN python3 -m pip install papermill