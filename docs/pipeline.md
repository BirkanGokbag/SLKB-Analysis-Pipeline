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

### Starting with a local database (SLKB Schema)

First, we start with creating a local database to store the CRISPR synthetic lethality data at hand. We start by creating a database named **myCDKO_db** at the desired location.

Note: The pipeline will add in the scores one at a time **after** the counts and sequence files will be deposited to the database. Thus, foreign key (FK) constrait will fail. For that reason, in the pipeline, the starting database starts with foreign keys disabled by default. 

Alternatively, you could enable them, run every scoring method, and insert into the tables in one transaction. However, due to the runtimes in some of the scores being very high, this is not recommended. 

```
db_engine = create_SLKB(location = os.getcwd(), name = 'myCDKO_db', disable_foreign_keys = True)
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

### Accessing Database

Database contents can be accessed externally, and also to insert records. In this pipeline, sqlalchemy will be used to load in the database we have just created. The following codes need to be ran to ensure score calculation is correct. 

```
SLKB_engine = sqlalchemy.create_engine(db_engine)
```

### Inserting Data to Database

Following data preperation, the data is now ready to be processed. We will need to declare two additonal variables before calling the processing function:

1. control_list: List of control targets. These need to be included at counts reference; **sgRNA_guide_1/2_name**
2. timepoint_list: A list of two elements; T0 replicates and TEnd replicates. Make sure that the replicate names align.

Example:
```
sequence_ref = ...
counts_ref = ...
scores_ref = ...
study_controls = ['CONTROL']
study_conditions = [['T0_rep1', 'T0_rep2', 'T0_rep3'],
                    ['TEnd_rep1', 'TEnd_rep2', 'TEnd_rep3']]

db_inserts = prepare_study_for_export(sequence_ref = sequence_ref, 
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
insert_study_to_db(dSLKB_engine, db_inserts)
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
experiment_design = pd.read_sql_table('CDKO_EXPERIMENT_DESIGN', SLKB_engine, index_col = 'sgRNA_id')
experiment_design.reset_index(drop = True, inplace = True)
experiment_design.index.rename('sgRNA_id', inplace = True)

# counts
counts = pd.read_sql_table('CDKO_SGRNA_COUNTS', SLKB_engine, index_col = 'sgRNA_pair_id')
counts.reset_index(drop = True, inplace = True)
counts.index.rename('sgRNA_pair_id', inplace = True)

# scores
scores = pd.read_sql_table('CDKO_ORIGINAL_SL_RESULTS', SLKB_engine, index_col = 'id')
scores.reset_index(drop = True, inplace = True)
scores.index.rename('gene_pair_id', inplace = True)

# join the tables together
counts = counts.merge(experiment_design, how = 'left', left_on = 'guide_1_id', right_index = True, suffixes = ('', '_g1'))
counts = counts.merge(experiment_design, how = 'left', left_on = 'guide_2_id', right_index = True, suffixes = ('', '_g2'))
# rename
counts = counts.rename({'sgRNA_guide_name': 'sgRNA_guide_name_g1',
                        'sgRNA_guide_seq': 'sgRNA_guide_seq_g1',
                        'sgRNA_target_name': 'sgRNA_target_name_g1',
                        'study_origin_x': 'study_origin',
                        'cell_line_origin_x': 'cell_line_origin'}, axis = 1)


curr_counts = counts[(counts['study_origin'] == 'myStudy') & (counts['cell_line_origin'] == 'myCL')]

```

For MAGeCK Score, GEMINI Score, and Horlbeck Score, files will be created in process. You can specify the location to save your files (default: current working directory). This is done in order to enable quick loading to database for repeated analyses. GEMINI Score and MAGeCK score require file generation in order to run. In the event of updated counts file (e.g., adding additional counts), setting the parameter ```restart_analysis=TRUE``` will restart the analysis from scratch. 


#### Median-B/NB Score

```
if check_if_added_to_table(curr_counts.copy(), 'MEDIAN_NB_SCORE', SLKB_engine):
    median_res = run_median_scores(curr_counts.copy())
    add_table_to_db(curr_counts.copy(), median_res['MEDIAN_NB_SCORE'], 'MEDIAN_NB_SCORE', SLKB_engine)
    if median_res['MEDIAN_B_SCORE'] is not None:
        add_table_to_db(curr_counts.copy(), median_res['MEDIAN_B_SCORE'], 'MEDIAN_B_SCORE', SLKB_engine)
```

#### sgRNA-Derived-B/NB Score

```
if not check_if_added_to_table(curr_counts.copy(), 'SGRA_DERIVED_NB_SCORE', SLKB_engine):
    sgRNA_res = run_sgrna_scores(curr_counts.copy())
    add_table_to_db(curr_counts.copy(), sgRNA_res['SGRA_DERIVED_NB_SCORE'], 'SGRA_DERIVED_NB_SCORE', SLKB_engine)
    if sgRNA_res['SGRA_DERIVED_B_SCORE'] is not None:
        add_table_to_db(curr_counts.copy(), sgRNA_res['SGRA_DERIVED_B_SCORE'], 'SGRA_DERIVED_B_SCORE', SLKB_engine)
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