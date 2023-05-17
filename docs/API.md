# API

## File Formats

A demo data is available for loading. Additional details can be found in the [pipeline](pipeline.md).

```
demo_data = SLKB.load_demo_data()
```
**Params**:

* None.

**Returns**:

* demo_data. A list of 3 items: sequence file, counts, fle, and score file.

## Functions For Data Insertion to KB

### generate_database

Creates a sqlite3 or mysql database, using SLKB schema.

```
SLKB.create_SLKB(engine = 'sqlite:///SLKB_sqlite3', db_type = 'sqlite3')
```

**Params**:

* engine: sqlalchemy url object. (Default: sqlite:///SLKB_sqlite3)
* db_type: Type of database to use schema for, currently available in mysql and sqlite3. (Default: sqlite3)

**Returns**:

* None.

<hr>

## extract_SLKB_webapp

Extracts the SLKB webapp to the specified location.

```
SLKB.extract_SLKB_webapp(location = os.getcwd())
```

**Params**:

* location: Location to extract SLKB files. (Default: Current working directory)

**Returns**:

* None.

<hr>

### prepare_study_for_export

Prepares the counts, scores, and sequences files for insertion into the DB.

```
db_inserts = SLKB.prepare_study_for_export(score_ref, sequence_ref = None, counts_ref = None, study_controls = None, study_conditions = None, can_control_be_substring = True, remove_unrelated_counts = False)
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
SLKB.insert_study_to_db(SLKB_engine, db_inserts)
```

**Params**:

* SLKB_engine: SQLAlchemy engine link
* db_inserts: Processed data, obtained via ```prepare_study_for_export```

**Returns**:

* None

<hr>

### Scoring Functions

#### Median-B/NB Score

Calculates Median B/NB Scores.

```
median_res = SLKB.run_median_scores(curr_counts, curr_study, curr_cl, full_normalization = False, re_run = False, store_loc = os.getcwd(), save_dir = 'MEDIAN_Files')
```

**Params**:

* curr_counts: Counts to calculate scores to.)
* curr_study: String, name of study to analyze data for.
* curr_cl: String, name of cell line to analyze data for.
* full_normalization: Whether to normalize counts across the whole sample or according to target type (Default: False)
* re_run: Boolean. Recreate and rerun the results instead of loading for subsequent analyses (Default: False)
* store_loc: String: Directory to store the MAGeCK files to. (Default: current working directory)
* save_dir: String: Folder name to store the MAGeCK files to. (Default: 'MEDIAN_Files')

**Returns**:

* median_res: A dictionary of two pandas dataframes: Median-B and Median-NB.

#### sgRNA-Derived-B/NB Score

Calculates sgRNA Derived N/NB scores.

sgRNA_res = SLKB.run_sgrna_scores(curr_counts, curr_study, curr_cl, full_normalization = False, re_run = False, store_loc = os.getcwd(), save_dir = 'sgRNA-DERIVED_Files')

**Params**:

* curr_counts: Counts to calculate scores to.)
* curr_study: String, name of study to analyze data for.
* curr_cl: String, name of cell line to analyze data for.
* full_normalization: Whether to normalize counts across the whole sample or according to target type (Default: False)
* re_run: Boolean. Recreate and rerun the results instead of loading for subsequent analyses (Default: False)
* store_loc: String: Directory to store the MAGeCK files to. (Default: current working directory)
* save_dir: String: Folder name to store the MAGeCK files to. (Default: 'sgRNA-DERIVED_Files')

**Returns**:

* sgRNA_res: A dictionary of two pandas dataframes: sgRNA_derived_B and sgRNA_derived_NB. 

#### MAGeCK Score

Calculates MAGeCK Score. Score files will created at the designated store location and save directory. 

```
mageck_res = SLKB.run_mageck_score(curr_counts.copy(), curr_study, curr_cl, store_loc = os.getcwd(), save_dir = 'MAGECK_Files', command_line_params = [],re_run = False)   
```

**Params**:
* curr_counts: Counts to calculate scores to.)
* curr_study: String, name of study to analyze data for.
* curr_cl: String, name of cell line to analyze data for.
* store_loc: String: Directory to store the MAGeCK files to. (Default: current working directory)
* save_dir: String: Folder name to store the MAGeCK files to. (Default: 'MAGECK_Files')
* command_line_params: Optional list to load programming environment(s) to be able to run mageck tool (i.e. loading path, activating python environment). 
* re_run: Boolean. Recreate and rerun the results instead of loading for subsequent analyses (Default: False)


**Returns**:

* mageck_res: A dict that contains a pandas dataframe for MAGeCK Score.

#### Horlbeck Score

Calculates Horlbeck score. Score files will created at the designated store location and save directory. 
```
horlbeck_res = SLKB.run_horlbeck_score(curr_counts.copy(), curr_study = curr_study, curr_cl = curr_cl, store_loc = os.getcwd(), save_dir = 'HORLBECK_Files', do_preprocessing = True, re_run = False)
```

**Params**:
* curr_counts: Counts to calculate scores to.)
* curr_study: String, name of study to analyze data for.
* curr_cl: String, name of cell line to analyze data for.
* store_loc: String: Directory to store the MAGeCK files to. (Default: current working directory)
* save_dir: String: Folder name to store the MAGeCK files to. (Default: 'Horlbeck_Files')
* do_preprocessing: Boolean. Run Horlbeck preprocessing (Default: True)
* re_run: Boolean. Recreate and rerun the results instead of loading for subsequent analyses (Default: False)

**Returns**:

* horlbeck_res: A dict that contains a pandas dataframe for Horlbeck Score.

#### GEMINI Score

Calculates GEMINI Score. Score files will created at the designated store location and save directory. 

```
gemini_res = run_gemini_score(curr_counts.copy(), curr_study = curr_study, curr_cl = curr_cl, store_loc = os.getcwd(), save_dir = 'GEMINI_Files', command_line_params = cmd_params, re_run = False)
```

**Params**:
* curr_counts: Counts to calculate scores to.)
* curr_study: String, name of study to analyze data for.
* curr_cl: String, name of cell line to analyze data for.
* store_loc: String: Directory to store the MAGeCK files to. (Default: current working directory)
* save_dir: String: Folder name to store the MAGeCK files to. (Default: 'GEMINI_Files')
* command_line_params: Optional list to load programming environment(s) to be able to run GEMINI through R (i.e. loading path, activating R environment). 
* re_run: Boolean. Recreate and rerun the results instead of loading for subsequent analyses (Default: False)

**Returns**:

* gemini_res: A dict that contains a pandas dataframe for GEMINI Score.


### check_if_added_to_table

If running the scoring methods multiple times, the method may be useful in skipping over the computation if there are gene pair records already in the database.

```
SLKB.check_if_added_to_table(curr_counts, score_name, SLKB_engine)
```

**Params**:

* curr_counts: Counts to calculate the scores to.
* score_name: Table to insert the scores to. Must be any of the 7 scoring table names:
    * HORLBECK_SCORE
    * MEDIAN_B_SCORE
    * MEDIAN_NB_SCORE
    * GEMINI_SCORE
    * MAGECK_SCORE
    * SGRNA_DERIVED_B_SCORE
    * SGRNA_DERIVED_NB_SCORE

**Returns**:

* Boolean. True if records are inserted into the DB, False otherwise.

### query_results_table

Obtain SL Scores from the specified scoring table.

```
result = SLKB.query_result_table(curr_counts, table_name, curr_study, curr_cl, engine_link)
```

**Params**:

* curr_counts: Counts to obtain the scores from.
* table_name: Must be any of the 7 scoring table names:
    * HORLBECK_SCORE
    * MEDIAN_B_SCORE
    * MEDIAN_NB_SCORE
    * GEMINI_SCORE
    * MAGECK_SCORE
    * SGRNA_DERIVED_B_SCORE
    * SGRNA_DERIVED_NB_SCORE
* curr_study: String, name of the study to obtain the results for.
* curr_cl: String, name of the cell line to obtain the results for.
* engine_link: SQLAlchemy connection for the database.

**Returns**:

* result: A pandas dataframe of the inserted results. Includes annotations for gene pair, study origin, and cell line origin.

© Copyright 2023, The Ohio State University, College of Medicine, Department of Biomedical Informatics