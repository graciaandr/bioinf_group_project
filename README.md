# MSc Bioinformatics BIO727P Software Group Project
### SNP Data Base for Population Genomics  
##### Gracia A., Celine L. and Amanah L.W.

The aim of the website is to provide researchers with a  tool to investigate the effect of different evolutionary processes such as genetic drift and natural selection on genetic variation within population using Shannon Diversity index, Tajima’s D and heterozygosity as well as measure the genetic differentiation between populations indicated by the FST statistic.

### Running the software locally:

In order to deploy the website locally, download the team_celine_webapp} directory from the GitHub page (see https://github.com/graciaandr/bioinf_group_project) and the database file in \texttt{.db} format from Google Drive using the following link:

-- drivelink missing --

Install the required packages to the local machine by running:

````
$python -m pip install -upgrade pip
$pip install -r requirements.txt
\end{lstlisting}
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
