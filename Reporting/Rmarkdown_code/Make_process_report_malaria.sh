#! /usr/bin/bash

currentPath="/Users/davejacobson/Desktop/DataScience_Trainings/Portfolio/Reporting/Rmarkdown_code"
echo $currentPath
docker run --mount type=bind,src=/Users/davejacobson/Desktop/,dst=/Users/davejacobson/Desktop/ reporting_image Rscript $currentPath/DMS_report_portfolio.R