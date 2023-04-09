## Pipeline

### Before getting started

Make sure that your SLKB pipeline package is appropriately installed. In order to run GEMINI Score and MAGeCK Score, you need to follow two additional steps:

1. GEMINI Score: Make sure that you have an R environment with GEMINI installed (version >= 2.1.1).
2. MAGeCK Score: Make sure that you have MAGECK installed into your system path. 

If you cannot change your system path and/or need to load R environment seperately (e.g. working on HPC systems), this will be covered later.

You can check whether you can access your R environment and MAGeCK on your computer.

```
import shutil
shutil.which('R') ## should yield accessed R environment location
shutil.which('mageck') ## should yield MAGeCK location
```

<hr>

### Starting with a local database (SLKB Schema)

First, we start with creating a local database to store the CRISPR synthetic lethality data at hand. We start by creating a database named **myCDKO_db** at the desired location.

```
create_SLKB(location = os.getcwd(), name = 'myCDKO_db')
```

### Preparing Data for Insert

#### sgRNA sequences

Sequences reference file consists of 3 columns:

1. sgRNA_guide_name: Name of your sgRNA guide
2. sgRNA_guide_sequence: Sequence for your sgRNA guide
3. sgRNA_guide_name: Target gene of the sgRNA guide

Example:
<center>

|sgRNA_guide_name|sgRNA_guide_sequence|sgRNA_guide_target|
|:-:|:-:|:-:|
|Gene1_sg1| CGCGC | Gene1 |
|Gene1_sg2| CGTGC  |  Gene1 |
|Gene2_sg1| CAAGC  |  Gene2 |

</center>

#### Dual sgRNA counts

Counts file contains X columns. Make sure that your counts annotation matches with the prior sequence information.

1. sgRNA_guide_1_name: Name of your sgRNA guide at first location
2. sgRNA_guide_1_target: Target of your sgRNA guide at first location
3. sgRNA_guide_2_name: Name of your sgRNA guide at second location
4. sgRNA_guide_2_target: Target of your sgRNA guide at second location
5. study_origin: Study identifier of the added data (i.e. PubmedID, MYCDKO)
6. cell_line_origin: Cell line identifier of the added data (e.g. 22Rv1)
7. replicate_names: Replicate names combined, separated by ';'
8. count_replicates: sgRNA counts from each replicate combined, separated by ';'


Example:
<!-- <center> -->

|sgRNA_guide_1_name|sgRNA_guide_1_target|sgRNA_guide_2_name|sgRNA_guide_2_target|study_origin|cell_line_origin|replicate_names|count_replicates|
|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|
|   |   |   |   |   |   |   |   |
|   |   |   |   |   |   |   |   |
|   |   |   |   |   |   |   |   |

#### Calculated SL Scores

If you have calculated the gene level SL scores for your data already, you can add them here with the following columns. Otherwise, leave empty. 

1. Gene_A
2. Gene_B
3. study_origin
4. cell_line_origin
5. SL_score
6. SL_score_cutoff
7. Stat_score
8. stat_score_cutoff

Example:

|sgRNA_guide_name|sgRNA_guide_sequence|sgRNA_guide_target|
|:-:|:-:|:-:|
|Gene1_sg1| CGCGC | Gene1 |
|Gene1_sg2| CGTGC  |  Gene1 |
|Gene2_sg1| CAAGC  |  Gene2 |

<hr>

### Inserting Data to Database

Following data preperation, the data is now ready to be processed. We will need to declare two additonal variables before calling the processing function:

1. control_list: List of control targets. These need to be included at counts reference; **sgRNA_guide_1/2_name**
2. timepoint_list: A list of two elements; T0 replicates and TEnd replicates. Make sure that the replicate names align.

Example:
```
sequence_ref = ...
counts_ref = ...
scores_ref = ...
control_list = ['CONTROL']
timepoint_list = [['T0_rep1', 'T0_rep2', 'T0_rep3'],
                    ['TEnd_rep1', 'TEnd_rep2', 'TEnd_rep3']]

db_inserts = (sequence_ref = sequence_ref, 
                counts_ref = counts_ref,
                scores_ref = scores_ref,
                control_list = control_list,
                timepoint_list = timepoint_list)
```

If passed checks successfully, you will notice that db_inserts contains the 3 items: sequence, counts, and score reference. In the event no scores reference was given, a dummy score of 0 was given to each possible gene pair. This is done in order to make sure that gene pairs are unique to each study and cell line.

If SL scores are supplied, by default, SL scores and statistical scores below the specified threshold are deemed as SL (column SL_or_not). Otherwise, they are not SL. You can customize this behavior by accessing the yielding db_inserts['scores_ref'].

By default, control genes supplied in scores file are removed. 


Finally, data can be inserted to the database.

```
db_inserts
```

### Accessing Database

Database contents can be accessed externally. In this pipeline, sqlalchemy will be used to load in the database we have just created. The following codes need to be ran to ensure score calculation is correct. 

### Calculating SL Scores and Inserting to Database

Score calculation methods are independent of each other. They can be ran in any order. The details of each scoring method are located in the original paper. Each score is accompanied with two helper functions; checking whether scores have been added to the database and inserting scores to the database.

1. check
2. insert

#### Initial steps

SL scores for the gene pairs are calculated for each cell line individually under each study. First, we filter the counts to obtain the study counts, followed by the cell line counts.

```
study_counts = ...
cl_counts = study_counts[...]

```

For MAGeCK Score, GEMINI Score, and Horlbeck Score, files will be created in process. You can specify the location to save your files (default: current working directory). This is done in order to enable quick loading to database for repeated analyses. GEMINI Score and MAGeCK score require file generation in order to run. In the event of updated counts file (e.g., adding additional counts), setting the parameter ```restart_analysis=TRUE``` will restart the analysis from scratch. 


#### Median-B/NB Score

```
if ~check:
    run_median
    insert
```

#### sgRNA-Derived-B/NB Score

```
if ~check:
    run_median
    insert
```

#### MAGeCK Score

In MAGeCK score, files will be created in process. You can specify the location to save your files (default: current working directory). 

```
if ~check:
    run_median
    insert
```

#### Horlbeck Score

In Horlbeck score, files will be created in process. You can specify the location to save your files (default: current working directory). 

```
if ~check:
    run_median
    insert
```

#### GEMINI Score

In GEMINI score, files will be created in process. You can specify the location to save your files (default: current working directory). Scores will be stored following GEMINI analysis for quick inserts to the database. 

```
if ~check:
    run_median
    insert
```

### Query Results

Following the score calculations, the query is relatively easy. In this snippet of code, we will be iterating over the all available data in the database and creating a score table that can be exported

```
loop and export
```

### Further Analyses

SLKB web application is available for download to help analyze your generated data. You can access the website at the following link [], and it's code at the link []. You will need to enter the database into the KB/ folder and calculated scores into the www/ folder with their appropriate names (db: SLKB.sqlite3, www: SLKB_calculated_scores.csv)