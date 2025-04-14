## Background

This folder hosts input files and code used to make rmarkdown html reports and dashboards. The input files, codes, and output were made to mimic what I generated at CDC, as part of the Domestic Parastic Surveillance team. The html reports would be shared with epidemiologists at CDC and State Public Health Labs. The purpose of the html reports was to share updates on the status of the specimens that we have analyzed and provide information on limited genotyping results.

## Data

I am not able to share any of the data or reports that I used/created at CDC. Therefore, I made fake data to mimic the structure of the data I used, and shows the extracting/transforming of data necessary to create the reports. 

## Reports

The rmarkdown reports are html documents and are not rendered in GitHub. You will need to download the gzipped example repeorts, unzip them, and then open the html files in a file browser to view the full report

The dashboards I created at CDC were primarily in PowerBI (with some Rshiny as well). I have decided to recreate these dashboards using plotly and dash in python. The idea is to showcase what was included in the various reports.

## Running code

I have included a docker file in the different folder that has the environment used to create the rmarkdown and dashboards.