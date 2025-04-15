## Background

This folder hosts demo reports/dashboards that I created while working at CDC. The primary purpose of these reports was to share updates on the status of the specimens that CDC analyzed for parasite domestic surveillance and provide limited genotyping reports. The reports' audience were collaborators at State Public Health Labs (both epidemiologists and laboratorians)

The reports shared with SPHLs were in html format, which cannot be rendered on GitHub. You will need to navigate to the HTML_process_reports folder to download and unzip the reports, which will then open in a web browser.

## Data / Code

I am not able to share any of the data or reports that I used/created at CDC. Therefore, I made fake data to mimic the structure that were used to generate the reports and dashboards. The code and input files used to generate the HTML reports are in the Rmarkdown_code directory.

Please note that the mutations reported in the *P.falciparum* section of the *Plasmodium* report are not actual true mutations that were detected, they are randomly generated mutations and thus do not represent any resistance profile that was ever detected.

## Reports

The dashboards I created at CDC were primarily in PowerBI (with some Rshiny as well). I have decided to recreate these dashboards using plotly and dash in python. The idea is to showcase what was included in the various reports.

## Running code

I have included a docker file in the different folder that has the environment used to create the rmarkdown and dashboards.