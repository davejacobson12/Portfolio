FROM ubuntu:22.04
# FROM rocker/r-base

ENV DEBIAN_FRONTEND=noninteractive

# RUN apt-get update && apt-get install -y --no-install-recommends build-essential r-base python3.10 python3-pip python3-setuptools python3-dev python3-pandas python3-openpyxl

RUN apt-get update && apt-get install -y --no-install-recommends  unzip cmake build-essential  curl pandoc pandoc-citeproc \
        autoconf \
        automake \
        bzip2 \
        lbzip2 \
        libbz2-dev \
        libcurl4-openssl-dev \
        liblzma-dev \
        libncurses5-dev \
        libssl-dev \
        libxml2-dev \
        libfontconfig1-dev \
        libmagick++-dev \
        pbzip2 \
        pigz \
        zlib1g-dev && \
        rm -rf /var/lib/apt/lists/* && \
        apt-get autoclean
        # gcc-multilib


RUN apt-get update && apt-get install -y r-base-core 
RUN apt-get clean

RUN R -e "install.packages('xml2',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('systemfonts',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('svglite',dependencies=TRUE, repos='http://cran.rstudio.com/')"


RUN R -e "install.packages('kableExtra',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('knitr',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('dplyr',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('stringr',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('lubridate',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('rmarkdown',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('tibble',dependencies=TRUE, repos='http://cran.rstudio.com/')"

# RUN install2.r -e kableExtra knitr dplyr stringr lubridate rmarkdown tibble ggplot2
# RUN install2.r -e xml2
# RUN install2.r -e systemfonts
# RUN install2.r -e svglite
# # RUN install2.r -e kableExtra

# RUN pip3 install dash --break-system-packages
# RUN pip3 install urllib3==1.26.6 --break-system-packages
# RUN pip3 install dash-bootstrap-components --break-system-packages
# # RUN pip3 install gunicorn --break-system-packages




