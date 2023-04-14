# API

## File Formats

For more details, follow the pipeline.

### Sequence File

### Counts File

### Score File

## Functions For Data Insertion to KB

### generate_database

Creates a local sqlite3 database, using SLKB schema.

```
db_engine = create_SLKB(location = os.getcwd(), name = 'myCDKO_db')
```

**Params**:

* location: Location to store the database
* name: Name of the database

**Returns**:

* db_engine: Database engine link with sqlite3.

<hr>

### prepare_study_for_export

Prepares the counts, scores, and sequences files for insertion into the DB.

```
db_inserts = prepare_study_for_export(score_ref, sequence_ref = None, counts_ref = None, study_controls = None, study_conditions = None, can_control_be_substring = True, remove_unrelated_counts = False)
```

**Params**:

* score_ref: A pandas table that adheres to the scores table template. 
* sequence_ref: A pandas table that adheres to the sequence table template (default: None). 
* counts_ref: A pandas table that adheres to the counts table template (default: None). 
* study_controls: A list of control targets of the sgRNAs (default: None).
* study_conditions: A list of two lists; first list contains the replicate names of initial time point, and second list contains the same for final time point (default: None).
* can_control_be_substring: Can the controls be a substring of gene targets (in case of possible name conventions: default: True)
* remove_unrelated_counts = Remove dual counts with targets that are outside of supplied scores targets? (default: False)

**Returns**:

* A dictionary of three items:
    * scores_ref: Contains the procesed scores table (if supplied)
    * sequences_ref: Contains the procesed sequences table (if supplied)
    * counts_ref: Contains the procesed counts table (if supplied)

### insert_study_to_db

Inserts the counts to the designated DB.

```
insert_study_to_db(SLKB_engine, db_inserts)
```

**Params**:

* SLKB_engine: SQLAlchemy engine link
* db_inserts: Processed data, obtained via ```prepare_study_for_export```

**Returns**:

* None

<hr>

## Functions for Scoring

command_line_params = []

### ...

### check_if_added_to_table

If running the scoring methods multiple times, the method may be useful in skipping over the computation if there are gene pair records already in the database.

```
check_if_added_to_table(curr_counts, score_name, SLKB_engine)
```

*Params**:

* curr_counts: Counts to calculate the scores to.
* score_name: Table to insert the scores to. Must be any of the 7 scoring table names:
    * HORLBECK_SCORE
    * MEDIAN_B_SCORE
    * MEDIAN_NB_SCORE
    * GEMINI_SCORE
    * MAGECK_SCORE
    * SGRA_DERIVED_B_SCORE
    * SGRA_DERIVED_NB_SCORE

**Returns**:

* Boolean. True if records are inserted into the DB, False otherwise.

**Helper Functions**

Helper functions are used across different functions within SLKB. You may use them as is or modify them to suit your needs.

## Helper Python Functions

###

## Helper R Functions

SLKB's helper functions are also located within the R environment. Functions for majority vote calculation, network visualization, and venn diagram creatios are located under the website's GitHub link. They are not included within the python package. 

To have access to those methods, you need to load in the functions to your R environment either through the URL or download them seperately.

```
source('R_loc')
```

### create_venn_diagram

Returns the venn diagram for the yielded results.

```
create_venn_diagram(results, top_p = '10percent', output_dir = 'venn_results')
```

**Params**:

* results: Obtain as a result of running the ```SLKB_query``` on Python.
* top_p: Top percentage of pairs to consider as text (Default: '10percent')
    * If entered a number and without the keyword 'percent', then that many top pairs are considered.
* output_dir: Directory to produce results (Default: 'venn_results' folder in current working directory)

**Returns**:

* None

### create_SL_network

Create the SL network based on the pair file provided via VisNetwork. Normally, this is either the scores file or results file that is filtered by the user to include only the specified pairs. The network's edges are color coded based on ```color_base```, and the interactions between SL pairs are highlighted if ```target_genes``` is provided.

```
network <- create_SL_network(pair_file, color_base = None, target_genes = None, output_dir = 'network_results')
```

**Params**:

* pair_file: DataFrame of score file (via dplyr and dbPool) and/or result file from ```SLKB_query```.
* color_base: Criteria to color the edges of the SL network. Options: None (Default), 'cell_line_origin', 'study_origin'
* target_genes: A list of genes that are targeted for color coding. (Default: None.)
* output_dir: Directory to produce results (Default: 'network_results' folder in current working directory)

**Returns**:

* A visNetwork object


**Additional Functions**

The following two visualization functions are included in the package to allow the user to visualize their results quickly on the spot.

### Dim_Red_Plot

Creates a 2D plot of the data, displays it and returns it.

```
figure = DGCyTOF.Dim_Red_Plot(name, data, labels, no_classes, class_names)
```

**Params**:

* name: Name of the plot 
* data: 2 dimensional embedding of data, for instance UMAP transformed data
* no_classes: List of classes
* class_names: List of class names

**Returns**:

* figure: Plt plot

<hr>

### Dim_Red_Plot_3d

Creates a 3D plot of the data, displays it and returns it.

```
figure = DGCyTOF.Dim_Red_Plot_3d(data, labels, all_celltypes)
```

**Params**:

* data: Low dimension embedding of data, for instance UMAP transformed data
* labels: List of labels corresponding to each cell
* all_celltypes: List of celltype names, where each celltype is in its respective label position

**Returns**:

* figure: Plt plot

<hr>

© Copyright 2023, The Ohio State University, College of Medicine, Biomedical Informatics