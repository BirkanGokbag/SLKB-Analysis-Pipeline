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

A template has been prepared that follows through the guide's steps. Feel free to install it from its GitHub [link](https://www.google.com) following SLKB package's installation. 
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
insert_study_to_db(SLKB_engine, db_inserts)
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

In MAGeCK score, files will be created in process. You can specify the location to save your files (default: current working directory). If you wish to re-run to store new results in its stead, set ```re_run``` to True.

MAGeCK is run through a script file at the designated location. If you need to load in any packages or set path to your mageck installation, please supply them to cmd_params.

```
cmd_params = []
if not check_if_added_to_table(curr_counts.copy(), 'MAGECK_SCORE', SLKB_engine):
    mageck_res = run_mageck_score(curr_counts.copy(), curr_study = curr_study, curr_cl = curr_cl, store_loc = os.getcwd(), save_dir = 'MAGECK_Files', command_line_params = cmd_params,re_run = False)
    add_table_to_db(curr_counts.copy(), mageck_res['MAGECK_SCORE'], 'MAGECK_SCORE', SLKB_engine)
        
```

#### Horlbeck Score

In Horlbeck score, files will be created in process. You can specify the location to save your files (default: current working directory). If you wish to re-run to store new results in its stead, set ```re_run``` to True.

```
if not check_if_added_to_table(curr_counts.copy(), 'HORLBECK_SCORE', SLKB_engine):
    horlbeck_res = run_horlbeck_score(curr_counts.copy(), curr_study = curr_study, curr_cl = curr_cl, store_loc = os.getcwd(), save_dir = 'HORLBECK_Files', do_preprocessing = True, re_run = False)
    add_table_to_db(curr_counts.copy(), horlbeck_res['HORLBECK_SCORE'], 'HORLBECK_SCORE', SLKB_engine)
```

#### GEMINI Score

In GEMINI score, files will be created in process. You can specify the location to save your files (default: current working directory). Scores will be stored following GEMINI analysis for quick inserts to the database. If you wish to re-run to store new results in its stead, set ```re_run``` to True.

Similarly to MAGeCK, GEMINI is run through a script file at the designated location. If you need to load in any packages or set path to your mageck installation, please supply them to cmd_params.

```
cmd_params = ['module load R/4.1.0']
if not check_if_added_to_table(curr_counts.copy(), 'GEMINI_SCORE', SLKB_engine):
    gemini_res = run_gemini_score(curr_counts.copy(), curr_study = curr_study, curr_cl = curr_cl, store_loc = os.getcwd(), save_dir = 'GEMINI_Files', command_line_params = cmd_params, re_run = False)
    add_table_to_db(curr_counts.copy(), gemini_res['GEMINI_SCORE'], 'GEMINI_SCORE', SLKB_engine)
```

### Query Results

Following the score calculations, the query is relatively easy. In this snippet of code, we will be iterating over the all available data in the database and creating a score table that can be exported.

```

# read the data

# experiment design
experiment_design = pd.read_sql_table('CDKO_EXPERIMENT_DESIGN', SLKB_engine, index_col = 'sgRNA_id')
experiment_design.drop(['study_origin'], axis = 1, inplace = True)
experiment_design.reset_index(drop = True, inplace = True)
experiment_design.index.rename('sgRNA_id', inplace = True)

# counts
counts = pd.read_sql_table('CDKO_SGRNA_COUNTS', SLKB_engine, index_col = 'sgRNA_pair_id')
counts.reset_index(drop = True, inplace = True)
counts.index.rename('sgRNA_pair_id', inplace = True)

# scores
scores = pd.read_sql_table('CDKO_ORIGINAL_SL_RESULTS', SLKB_engine, index_col = 'gene_pair_id')

# join the tables together
counts = counts.merge(scores, how = 'left', left_on = 'gene_pair_id', right_index = True)
counts = counts.merge(experiment_design, how = 'left', left_on = 'guide_1_id', right_index = True, suffixes = ('', '_g1'))
counts = counts.merge(experiment_design, how = 'left', left_on = 'guide_2_id', right_index = True, suffixes = ('', '_g2'))
# rename
counts = counts.rename({'sgRNA_guide_name': 'sgRNA_guide_name_g1',
                        'sgRNA_guide_seq': 'sgRNA_guide_seq_g1',
                        'sgRNA_target_name': 'sgRNA_target_name_g1',
                        'study_origin_x': 'study_origin',
                        'cell_line_origin_x': 'cell_line_origin'}, axis = 1)

experiment_design = pd.read_sql_table('CDKO_EXPERIMENT_DESIGN', SLKB_engine, index_col = 'sgRNA_id')
experiment_design.reset_index(drop = True, inplace = True)
experiment_design.index.rename('sgRNA_id', inplace = True)

# tables to obtain the data from
all_results_tables = ['HORLBECK_SCORE', 
                      'MAGECK_SCORE', 
                      'MEDIAN_NB_SCORE', 
                      'MEDIAN_B_SCORE', 
                      'SGRA_DERIVED_NB_SCORE', 
                      'SGRA_DERIVED_B_SCORE', 
                      'GEMINI_SCORE']

#################

all_scores = []
for curr_study in available_studies:
    print('Working on study: ' + curr_study)

    # get study counts and seq
    study_counts = counts.loc[counts['study_origin'] == study_name_to_pubmed_id[curr_study]].copy()

    curr_seq_ids = np.array(sorted(list(set(study_counts['guide_1_id'].tolist() + study_counts['guide_2_id'].tolist()))))
    study_sequences = experiment_design.loc[curr_seq_ids]

    # the analysis runs for each individual cell line
    available_cell_lines = set(study_counts['cell_line_origin'])


    for curr_cl in available_cell_lines:
        # store results here
        study_scores = []
    
        print('Working on cell line: ' + curr_cl)
        curr_counts = study_counts.loc[study_counts['cell_line_origin'] == curr_cl].copy()
        
        for table_name in all_results_tables:
            # add the result of the table to the list
            study_scores.append(query_result_table(curr_counts.copy(), table_name, curr_study, curr_cl, SLKB_engine))
    
        # remove duplicate annotation columns
        study_scores = pd.concat(study_scores, axis = 1, ignore_index = False)
        study_scores = study_scores.loc[:,~study_scores.columns.duplicated(keep = 'last')].copy()
        
        # make sure the annotations are all filled
        study_scores['gene_pair'] = study_scores.index
        study_scores['study_origin'] = curr_study
        study_scores['cell_line_origin'] = curr_cl
        
        # reset the index, gene_pair -> id
        study_scores.reset_index(drop = True, inplace = True)

        # add to big table 
        all_scores.append(study_scores)
    
    print('-----')
    
print('Done getting all data!')
    
# combine the scores at the end
all_scores = pd.concat(all_scores, axis = 0, ignore_index = True)

# add individual genes
all_scores['gene_1'] = [i.split('|')[0] for i in all_scores['gene_pair']]
all_scores['gene_2'] = [i.split('|')[1] for i in all_scores['gene_pair']]

# sort such that all annotations are at the front
all_columns = sorted(list(all_scores.columns))
annotation_columns = ['gene_pair', 'gene_1', 'gene_2', 'study_origin', 'cell_line_origin']

# get the final scores
all_scores = all_scores.loc[:, annotation_columns + [i for i in all_columns if i not in annotation_columns]]
```

### Further Analyses

SLKB web application is available for download to help analyze your generated data. You can access the website at the following [link](https://www.google.com), and it's code at the link [link](https://www.google.com). You will need to enter the database into the KB/ folder and calculated scores into the www/ folder with their appropriate names (db: SLKB.sqlite3, www: SLKB_calculated_scores.csv)

Alternatively, you can access its helper functions. Discussed in the [API](API.md).