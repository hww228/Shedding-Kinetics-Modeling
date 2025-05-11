FROM rocker/rstudio:4.4.3

# Update package list and install dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    gcc \
    g++ \
    gfortran

# Install rstan package
RUN R -e "install.packages('rstan')"
RUN R -e "install.packages('readxl')"
RUN R -e "install.packages('dplyr')"

COPY . /home/rstudio/
WORKDIR /home/rstudio/