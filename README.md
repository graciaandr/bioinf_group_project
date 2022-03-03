# MSc Bioinformatics BIO727P Software Group Project
### SNP Data Portal for Population Genomics  
##### Gracia A., Celine L. and Amanah L.W.

The aim of the website is to provide researchers with a  tool to investigate the effect of different evolutionary processes such as genetic drift and natural selection on genetic variation within population using Shannon Diversity index, Tajima’s D and heterozygosity as well as measure the genetic differentiation between populations indicated by the FST statistic.

### Running the software locally:

In order to deploy the website locally, download the team_celine_webapp} directory from the GitHub page (see https://github.com/graciaandr/bioinf_group_project) and the database file in \texttt{.db} format from Google Drive using the following link:
https://drive.google.com/drive/folders/1g5-2Yr40tCxM30NLZQ0zmyzaBFjKPJIq?usp=sharing

Install the required packages to the local machine by running:

````
$python -m pip install -upgrade pip
$pip install -r requirements.txt
````
Run the software on http://localhost:5000/ using:

```
$git clone https://github.com/graciaandr/bioinf_group_project.git
$cd webapp
$FLASK_APP=webapp.py
$flask run
```

### Required packages:
````
click==8.0.4
cloudpickle==2.0.0
cycler==0.11.0
dask==2022.2.1
dominate==2.6.0
Flask==2.0.3
Flask-Bootstrap==3.3.7.1
fonttools==4.29.1
fsspec==2022.2.0
itsdangerous==2.1.0
Jinja2==3.0.3
kiwisolver==1.3.2
locket==0.2.1
MarkupSafe==2.1.0
matplotlib==3.5.1
numpy==1.22.2
packaging==21.3
pandas==1.4.1
partd==1.2.0
Pillow==9.0.1
pyparsing==3.0.7
python-dateutil==2.8.2
pytz==2021.3
PyYAML==6.0
scikit-allel==1.3.5
scipy==1.8.0
seaborn==0.11.2
six==1.16.0
toolz==0.11.2
visitor==0.1.3
Werkzeug==2.0.3
````
#### Note:
If an error such as 

```
ERROR: Could not find a version that satisfies the requirement pandas==1.4.1 
(from versions: 0.1, 0.2, 0.3.0, 0.4.0, 0.4.1, 0.4.2, 0.4.3, 0.5.0, 0.6.0, 0.6.1, 0.7.0, 0.7.1, 0.7.2, 0.7.3, 0.8.0, 0.8.1, 0.9.0, 0.9.1, 0.10.0, 0.10.1, 0.11.0, 0.12.0, 0.13.0, 0.13.1, 0.14.0, 0.14.1, 0.15.0, 0.15.1, 0.15.2, 0.16.0, 0.16.1, 0.16.2, 0.17.0, 0.17.1, 0.18.0, 0.18.1, 0.19.0, 0.19.1, 0.19.2, 0.20.0, 0.20.1, 0.20.2, 0.20.3, 0.21.0, 0.21.1, 0.22.0, 0.23.0, 0.23.1, 0.23.2, 0.23.3, 0.23.4, 0.24.0, 0.24.1, 0.24.2, 0.25.0, 0.25.1, 0.25.2, 0.25.3, 1.0.0, 1.0.1, 1.0.2, 1.0.3, 1.0.4, 1.0.5, 1.1.0, 1.1.1, 1.1.2, 1.1.3, 1.1.4, 1.1.5, 1.2.0, 1.2.1, 1.2.2, 1.2.3, 1.2.4, 1.2.5, 1.3.0, 1.3.1, 1.3.2, 1.3.3, 1.3.4, 1.3.5)
``` 
occurs during installation of the requirements, please edit the ```requirements.txt``` such that the last recognised version in your device (here: ``1.3.5``) is in the file.

### Features of Software

* SNP search
* Browse for SNPs via rs ID, Gene Name (or Aliases), or Genomic Coordinates
* Summary statistics selection​
* Select between Shannon Diversity, Expected Heterozygosity and Tajima’s D
* Select FST Analysis
* Population selection​
* Choose one or more of the 5 provided populations 
* Download stats as TXT
* Download the FST or other stats​ data frame as a txt file
* Visualise summary stats
