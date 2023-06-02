## Pipeline

### Before getting started

Make sure that your SLKB pipeline package is appropriately installed. In order to run GEMINI Score and MAGeCK Score, you need to follow two additional steps:

1. GEMINI Score: Make sure that you have an R environment with GEMINI installed (version >= 2.1.1), and ggplot2 (will be installed alongside GEMINI).
2. MAGeCK Score: Make sure that you have MAGECK installed into your system path. 

If you cannot change your system path and/or need to load R environment seperately (e.g. working on HPC systems), this will be covered later.

You can check whether you can access your R environment and MAGeCK on your computer.

```
import shutil
shutil.which('R') ## should yield accessed R environment location
shutil.which('mageck') ## should yield MAGeCK location
```

<hr>

### SLKB Pipeline Template

A template has been prepared that follows through the guide's steps. Feel free to install it from its GitHub [link](https://github.com/BirkanGokbag/SLKB-Analysis-Pipeline/blob/main/SLKB/files/pipeline.ipynb) following SLKB package's installation. 
### Starting with a local database (SLKB Schema)

First, we start with creating a local database to store the CRISPR synthetic lethality data at hand. SLKB supports mysql and sqlite3 at this time. A URL object is passed that will then create the database via stored schemas.

```
url_object = sqlalchemy.URL.create(
    "mysql+mysqlconnector",
    username="root",
    password="password",  # plain (unescaped) text
    host="localhost",
    port = '3306',
    database="SLKB_mysql",
)
url_object = 'sqlite:///SLKB_sqlite3'

SLKB_engine = sqlalchemy.create_engine(url_object, pool_size = 0)

# create the database at the url_object
SLKB.create_SLKB(engine = SLKB_engine, db_type = 'sqlite3') # or mysql
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
|0Safe-safe-ACOC-204550.4522| GTGTATTTGGCTTCCAAAA | control |
|0Safe-safe-ACOC-204550.4525| GCATGGCCTCCACTTGCAA  |  control |
|AKT3-1| GTAAGGTAAATCCACATCTTG  |  AKT3 |

</center>

#### Dual sgRNA counts

Counts file contains X columns. Make sure that your counts annotation matches with the prior sequence information.

1. guide_1: Name of your sgRNA guide at first location
2. gene_1: Target of your sgRNA guide at first location (CONTROL if targeting control)
3. guide_2: Name of your sgRNA guide at second location
4. gene_2: Target of your sgRNA guide at second location (CONTROL if targeting control)
5. study_origin: Study identifier of the added data (i.e. PubmedID, MYCDKO)
6. cell_line_origin: Cell line identifier of the added data (e.g. 22Rv1)
7. study_conditions: Replicate names combined, separated by ';'
8. count_replicates: sgRNA counts from each replicate combined, separated by ';'


Example:

|guide_1|guide_2|gene_1|gene_2|count_replicates|cell_line_origin|study_conditions|study_origin|
|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|
|0Safe-safe-ACOC-204550.4522|0Safe-safe-ACOC-204550.4522|0Safe-safe-ACOC|0Safe-safe-ACOC|22;17;36;43|22Rv1|T0_1;T0_2;T12_1;T12_2|36060092|
|0Safe-safe-ACOC-204550.4522|0Safe-safe-ACOC-204550.4525|0Safe-safe-ACOC|0Safe-safe-ACOC|57;73;107;98|22Rv1|T0_1;T0_2;T12_1;T12_2|36060092|
|0Safe-safe-ACOC-204550.4522|AKT3-1|0Safe-safe-ACOC|AKT3|47;45;68;85|22Rv1|T0_1;T0_2;T12_1;T12_2|36060092|

#### Calculated SL Scores

If you have calculated the gene level SL scores for your data already, you can add them here with the following columns. Otherwise, leave empty. 

1. gene_1
2. gene_2
3. study_origin
4. cell_line_origin
5. SL_score
6. SL_score_cutoff (If NaN, no cutoff applied)
7. statistical_score
8. statistical_score_cutoff (If NaN, no cutoff applied)

Example:

|gene_1|gene_2|study_origin|cell_line_origin|SL_score|SL_score_cutoff|statistical_score|statistical_score_cutoff|
|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|
|AKT3|AR|36060092|22Rv1|-0.977141|-0.361|NaN|NaN|

<hr>

### Accessing Database

Database contents can be accessed externally, and also to insert records. In this pipeline, sqlalchemy will be used to load in the database we have just created. The following codes need to be ran to ensure score calculation is correct. 

```
SLKB_engine = sqlalchemy.create_engine(db_engine)
```

### Inserting Data to Database

Following data preperation, the data is now ready to be processed. We will need to declare two additonal variables before calling the processing function:

1. control_list: List of control targets. These need to be included at counts reference; **gene_1, gene_2 names to use as controls **
2. timepoint_list: A list of two elements; T0 replicates and TEnd replicates. Make sure that the replicate names align.

Example:
```
sequence_ref = ...
counts_ref = ...
scores_ref = ...
study_controls = ['CONTROL']
study_conditions = [['T0_rep1', 'T0_rep2', 'T0_rep3'],
                    ['TEnd_rep1', 'TEnd_rep2', 'TEnd_rep3']]

db_inserts = SLKB.prepare_study_for_export(sequence_ref = sequence_ref, 
                                                     counts_ref = counts_ref,
                                                     scores_ref = scores_ref,
                                                     study_controls = study_controls,
                                                     study_conditions = study_conditions)
```

If passed checks successfully, you will notice that db_inserts contains the 3 items: sequence, counts, and score reference. In the event no scores reference was given, a dummy score of 0 was given to each possible gene pair. This is done in order to make sure that gene pairs are unique to each study and cell line.

If SL scores are supplied, by default, SL scores and statistical scores below the specified threshold are deemed as SL (column SL_or_not). Otherwise, they are not SL. You can customize this behavior by accessing the yielding db_inserts['scores_ref'].

By default, control genes supplied in scores file are removed. 

Finally, data can be inserted to the database.

```
SLKB.insert_study_to_db(SLKB_engine, db_inserts)
```

### Calculating SL Scores and Inserting to Database

Score calculation methods are independent of each other. They can be ran in any order. The details of each scoring method are located in the original paper. Each score is accompanied with two helper functions; checking whether scores have been added to the database and inserting scores to the database.

1. check
2. insert

#### Initial steps

SL scores for the gene pairs are calculated for each cell line individually under each study and cell line. First, we filter the counts to obtain the study counts, followed by the cell line counts.

```

# read the data

# experiment design
experiment_design = pd.read_sql_query(con=SLKB_engine.connect(), 
                              sql=sqlalchemy.text('SELECT * from CDKO_EXPERIMENT_DESIGN'.lower()), index_col = 'sgRNA_id')
experiment_design.reset_index(drop = True, inplace = True)
experiment_design.index.rename('sgRNA_id', inplace = True)

# counts
counts = pd.read_sql_query(con=SLKB_engine.connect(), 
                              sql=sqlalchemy.text('SELECT * from joined_counts'.lower()), index_col = 'sgRNA_pair_id')

# scores
scores = pd.read_sql_query(con=SLKB_engine.connect(), 
                              sql=sqlalchemy.text('SELECT * from CDKO_ORIGINAL_SL_RESULTS'.lower()), index_col = 'id')
scores.reset_index(drop = True, inplace = True)
scores.index.rename('gene_pair_id', inplace = True)

```

For all scores, files will be created in the process. You can specify the location to save your files (default: current working directory). This is done in order to enable quick loading to database for repeated analyses. GEMINI Score and MAGeCK score require file generation in order to run. In the event of updated counts file (e.g., adding additional counts), setting the parameter ```re_run=TRUE``` will restart the analysis from scratch. 


#### Median-B/NB Score

```
if SLKB.check_if_added_to_table(curr_counts.copy(), 'median_nb_score', SLKB_engine):
    median_res = SLKB.run_median_scores(curr_counts.copy(), curr_study = curr_study, curr_cl = curr_cl, store_loc = os.getcwd(), save_dir = 'MEDIAN_Files')
    SLKB.add_table_to_db(curr_counts.copy(), median_res['MEDIAN_NB_SCORE'], 'median_nb_score', SLKB_engine)
    if median_res['MEDIAN_B_SCORE'] is not None:
        SLKB.add_table_to_db(curr_counts.copy(), median_res['MEDIAN_B_SCORE'], 'median_b_score', SLKB_engine)
```

#### sgRNA-Derived-B/NB Score

```
if not SLKB.check_if_added_to_table(curr_counts.copy(), 'sgrna_derived_nb_score', SLKB_engine):
    sgRNA_res = SLKB.run_sgrna_scores(curr_counts.copy(), curr_study = curr_study, curr_cl = curr_cl, store_loc = os.getcwd(), save_dir = 'sgRNA-DERIVED_Files')
    SLKB.add_table_to_db(curr_counts.copy(), sgRNA_res['SGRNA_DERIVED_NB_SCORE'], 'sgrna_derived_nb_score', SLKB_engine)
    if sgRNA_res['SGRNA_DERIVED_B_SCORE'] is not None:
        SLKB.add_table_to_db(curr_counts.copy(), sgRNA_res['SGRNA_DERIVED_B_SCORE'], 'sgrna_derived_b_score', SLKB_engine)
```

#### MAGeCK Score

In MAGeCK score, files will be created in process. You can specify the location to save your files (default: current working directory). If you wish to re-run to store new results in its stead, set ```re_run``` to True.

MAGeCK is run through a script file at the designated location. If you need to load in any packages or set path to your mageck installation, please supply them to cmd_params.

```
cmd_params = []
if not SLKB.check_if_added_to_table(curr_counts.copy(), 'mageck_score', SLKB_engine):
    mageck_res = SLKB.run_mageck_score(curr_counts.copy(), curr_study = curr_study, curr_cl = curr_cl, store_loc = os.getcwd(), save_dir = 'MAGECK_Files', command_line_params = cmd_params,re_run = False)
    SLKB.add_table_to_db(curr_counts.copy(), mageck_res['MAGECK_SCORE'], 'mageck_score', SLKB_engine)
        
```

#### Horlbeck Score

In Horlbeck score, files will be created in process. You can specify the location to save your files (default: current working directory). If you wish to re-run to store new results in its stead, set ```re_run``` to True.

```
if not SLKB.check_if_added_to_table(curr_counts.copy(), 'horlbeck_score', SLKB_engine):
    horlbeck_res = SLKB.run_horlbeck_score(curr_counts.copy(), curr_study = curr_study, curr_cl = curr_cl, store_loc = os.getcwd(), save_dir = 'HORLBECK_Files', do_preprocessing = True, re_run = False)
    SLKB.add_table_to_db(curr_counts.copy(), horlbeck_res['HORLBECK_SCORE'], 'horlbeck_score', SLKB_engine)
```

#### GEMINI Score

In GEMINI score, files will be created in process. You can specify the location to save your files (default: current working directory). Scores will be stored following GEMINI analysis for quick inserts to the database. If you wish to re-run to store new results in its stead, set ```re_run``` to True.

Similarly to MAGeCK, GEMINI is run through a script file at the designated location. If you need to load in any packages or set path to your mageck installation, please supply them to cmd_params.

```
cmd_params = ['module load R/4.1.0']
if not SLKB.check_if_added_to_table(curr_counts.copy(), 'gemini_score', SLKB_engine):
    gemini_res = SLKB.run_gemini_score(curr_counts.copy(), curr_study = curr_study, curr_cl = curr_cl, store_loc = os.getcwd(), save_dir = 'GEMINI_Files', command_line_params = cmd_params, re_run = False)
    SLKB.add_table_to_db(curr_counts.copy(), gemini_res['GEMINI_SCORE'], 'gemini_score', SLKB_engine)
```

### Query Results (For one table)

Following the score calculations, the query is relatively easy. In this snippet of code, we will access the scores for one of the tables.

```
score = SLKB.query_result_table(curr_counts.copy(), 'median_b_score', curr_study, curr_cl, SLKB_engine)
```

### Query Results (For all tables)

Here, we will query all scores using the view within the database.

```
all_scores = pd.read_sql_query(con=SLKB_engine.connect(), 
                              sql=sqlalchemy.text('SELECT * from calculated_sl_table'))
```

### Further Analyses

SLKB web application is available for download to help analyze your generated data. You can access the website at the following [link](https://slkb.osubmi.org/), and it's code at the [link](https://github.com/BirkanGokbag/SLKB-Analysis-Pipeline/blob/main/SLKB/files/SLKB_webapp.zip). To use user generated data, check ```server.r``` within the web app.