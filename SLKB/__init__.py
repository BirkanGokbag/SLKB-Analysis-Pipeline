# imports
import shutil
import numpy as np
import pandas as pd
import os
import math
import pickle
import sqlalchemy
from sqlalchemy.orm import sessionmaker
from scipy import optimize
from scipy.stats import sem
import subprocess

import pkg_resources
PACKAGE_PATH = pkg_resources.resource_filename('SLKB', '/')

def load_demo_data():
    '''
    A demo data is available for loading. Additional details can be found in the [pipeline](pipeline.md).

    **Params**:

    * None.

    **Returns**:

    * demo_data. A list of 3 items: sequence file, counts, fle, and score file.
    '''
    # load up the demo data and return
    data_loc = os.path.join(PACKAGE_PATH, 'files', 'demo_data.pickle')
    with open(data_loc, 'rb') as handle:
        toy_data = pickle.load(handle)
    return(toy_data)
    
    
def check_repeated_constructs(x, index_loc):
    '''
    Helper function, Returns location of counts with respect to study conditions/replicate names.
    '''
    if len(x) < max(index_loc):
        sub = index_loc[index_loc < len(x)]
        return(x[sub])
    else:
        return(x[index_loc])
    
def create_SLKB(engine = 'sqlite:///SLKB_sqlite3', db_type = 'sqlite3'):
    '''
    Creates a sqlite3 or mysql database, using SLKB schema.

    **Params**:

    * engine: sqlalchemy url object. (Default: sqlite:///SLKB_sqlite3)
    * db_type: Type of database to use schema for, currently available in mysql and sqlite3. (Default: sqlite3)

    **Returns**:

    * None.
    '''
    schema_loc = os.path.join(PACKAGE_PATH, 'files')
    if db_type == 'sqlite3':
        schema_loc = os.path.join(schema_loc, 'SLKB_sqlite3_schema.sql')
    elif db_type == 'mysql':
        schema_loc = os.path.join(schema_loc, 'SLKB_mysql_schema.sql')
    else:
        print('Unavailable. Please choose either sqlite3 or mysql.')

    # read the schema
    with open(schema_loc) as f:
        command = f.read()

    # execute
    with engine.begin() as transaction:
        for com in command.split(';'):
            transaction.execute(sqlalchemy.text(com)) 


def extract_SLKB_webapp(location = os.getcwd()):
    '''
    Extracts a SLKB webapp to the specified location.

    **Params**:

    * location: Location to extract SLKB files. (Default: Current working directory)

    **Returns**:

    * None.
    '''
    webapp_loc = os.path.join(PACKAGE_PATH, 'files', 'SLKB_webapp.zip')
    print('Extracting to location: ' + location)
    shutil.unpack_archive(webapp_loc, location)
    print('Done!')

###### Data Preperation Functions

def create_placeholder_scores(curr_counts, sequence_ref):
    # we should add genes to the KB that can later be modified following scoring
    
    # first, set the controls so they can be removed
    curr_counts.loc[curr_counts['guide_1'].isin(sequence_ref.loc[sequence_ref['sgRNA_target_name'] == 'control', 'sgRNA_guide_name'].values), 'Gene 1'] = 'CONTROL'
    curr_counts.loc[curr_counts['guide_2'].isin(sequence_ref.loc[sequence_ref['sgRNA_target_name'] == 'control', 'sgRNA_guide_name'].values), 'Gene 2'] = 'CONTROL'

    idx = (curr_counts['gene_1'] == 'CONTROL') | (curr_counts['gene_2'] == 'CONTROL')
    # remove them
    curr_counts = curr_counts[~idx]

    # add sorted genes so they can be removed
    curr_counts['sorted_genes'] = ['|'.join(sorted([curr_counts['gene_1'].iloc[i], curr_counts['gene_2'].iloc[i]])) for i in range(curr_counts.shape[0])]
    curr_counts.drop_duplicates(subset = ['sorted_genes', 'cell_line_origin'], keep = 'first', inplace = True)

    # drop the same genes as well
    curr_counts = curr_counts.loc[curr_counts['gene_1'] != curr_counts['gene_2']]
    curr_counts.reset_index(drop = True, inplace = True)

    # proceed to create the GI and return
    curr_GI = pd.DataFrame(columns = ["gene_1", "gene_2", "study_origin", "cell_line_origin", "SL_score", "SL_score_cutoff", "statistical_score", "statistical_score_cutoff"])
    curr_GI["gene_1"] = curr_counts['gene_1'].values
    curr_GI["gene_2"] = curr_counts['gene_2'].values
    curr_GI["study_origin"] = [curr_counts['study_origin'].iloc[0]] * curr_GI.shape[0]
    curr_GI["cell_line_origin"] = curr_counts['cell_line_origin'].values
    curr_GI = curr_GI.fillna(0)
    
    return(curr_GI)

def prepare_study_for_export(sequence_ref, counts_ref, score_ref, study_controls = None, study_conditions = None, can_control_be_substring = True, remove_unrelated_counts = False):
    '''
        
    Prepares the counts, scores, and sequences files for insertion into the DB.

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
    '''
    ## make sure the columns are within each table, if not return error
    sequence_ref_needed_columns = {'sgRNA_guide_name', 'sgRNA_guide_seq', 'sgRNA_target_name'}
    
    if sequence_ref is not None:
        if len(sequence_ref_needed_columns.difference(sequence_ref.columns)) > 0:
            print('Error')
            print('sequence_ref')
            return
        # reset index by default
        sequence_ref.sort_values('sgRNA_target_name', ignore_index = True, inplace = True)
        sequence_ref.reset_index(drop = True, inplace = True)
    
    counts_ref_needed_columns = {'guide_1', 'guide_2', 'gene_1', 'gene_2', 'count_replicates', 'cell_line_origin', 'study_origin', 'study_conditions'}
    if counts_ref is not None:
        if len(counts_ref_needed_columns.difference(counts_ref.columns)) > 0:
            print('Error')
            print('counts_ref')
            return
        # reset index by default
        counts_ref.reset_index(drop = True, inplace = True)
    
    score_ref_needed_columns = {'gene_1', 'gene_2', 'study_origin', 'cell_line_origin', 'SL_score', 'SL_score_cutoff', 'statistical_score', 'statistical_score_cutoff'}
    if (score_ref is None) and (counts_ref is not None):
        print('There are no scores, but there are counts...Generating Placeholder...')
        score_ref = create_placeholder_scores(counts_ref.copy(), sequence_ref.copy())
    if len(score_ref_needed_columns.difference(score_ref.columns)) > 0:
        print('Error')
        print('score_ref')
        return
    # reset index by default
    score_ref.reset_index(drop = True, inplace = True)
    
    if study_controls is not None:
        study_controls = [i.upper() for i in study_controls]
    
    
    ## prepare each table to be inserted to their respective tables
    
    print("Starting processing...")
    
    ################################# first, handle the scores ref
    print('Score reference...')
    
    # fill NA
    score_ref = score_ref.fillna(0)
    
    for col in ['gene_1', 'gene_2', 'cell_line_origin']:
        score_ref[col] = [i.upper() for i in score_ref[col]]
    
    genes_A = score_ref['gene_1'].values
    genes_B = score_ref['gene_2'].values
    sorted_genes = []
    for i in range(score_ref.shape[0]):
        sorted_genes.append('_'.join(sorted([genes_A[i], genes_B[i]])))
    score_ref["gene_pair"] = sorted_genes

    # remove same ones
    score_ref = score_ref.loc[~(score_ref["gene_1"].values == score_ref["gene_2"].values)]

    # remove controls from SL scores
    if study_controls is not None:
        control_idx = np.array([False] * score_ref.shape[0])

        for curr_control in study_controls:
            
            if can_control_be_substring:
                control_idx = control_idx | np.array([True if curr_control in i else False for i in score_ref["gene_1"]]) | np.array([True if curr_control in i else False for i in score_ref["gene_2"]])
                
            control_idx = control_idx | np.array([True if i in curr_control else False for i in score_ref["gene_1"]]) | np.array([True if i in curr_control else False for i in score_ref["gene_2"]])

        print('Controls within SL score that are removed: ')
        print(control_idx.sum())
        print('---')

        score_ref = score_ref.loc[~control_idx]
        
    if (score_ref['statistical_score_cutoff'].iloc[0] != 0) and (score_ref['SL_score_cutoff'].iloc[0] != 0):
        print('Both GI and Stat cutoffs are present...')
        score_ref['SL_or_not'] = (score_ref['SL_score'] <= (score_ref['SL_score_cutoff'].iloc[0])) & (score_ref['statistical_score'] <= (score_ref['statistical_score_cutoff'].iloc[0]))
    elif score_ref['SL_score_cutoff'].iloc[0] != 0:
        print('Only GI cutoff is present...')
        score_ref['SL_or_not'] = score_ref['SL_score'] <= score_ref['SL_score_cutoff'].iloc[0]
    elif score_ref['statistical_score'].iloc[0] != 0:
        print('Only Stat cutoff is present...')
        score_ref['SL_or_not'] = score_ref['statistical_score'] <= score_ref['statistical_score_cutoff'].iloc[0]
    else:
        print('No scores/stats cutoffs are available, possibly generated. Setting all to be NOT SL')
        score_ref['SL_or_not'] = [False] * score_ref.shape[0]
        
    score_ref.loc[score_ref['SL_or_not'], 'SL_or_not'] = 'SL'
    score_ref.loc[score_ref['SL_or_not'] != 'SL', 'SL_or_not'] = 'Not SL'
    
    ################################# score ref - DONE
    
    print('Counts reference...')
    
    if counts_ref is not None:

        for col in ['guide_1', 'guide_2', 'gene_1', 'gene_2', 'cell_line_origin']:
            counts_ref[col] = [i.upper() for i in counts_ref[col]]


        # label whether single, double, or control
        sgRNA_true_pair_index = np.array([i for i in range(counts_ref.shape[0]) if (str(counts_ref["gene_1"].iloc[i]) not in study_controls) and (str(counts_ref["gene_2"].iloc[i]) not in study_controls) and (str(counts_ref["gene_1"].iloc[i]) != str(counts_ref["gene_2"].iloc[i]))])
        print(' '.join(["Number of double pairs:", str(len(sgRNA_true_pair_index))]))

        sgRNA_control_pair_index = np.array([i for i in range(counts_ref.shape[0]) if (str(counts_ref["gene_1"].iloc[i]) in study_controls) and (str(counts_ref["gene_2"].iloc[i]) in study_controls)])
        print(' '.join(["Number of controls:", str(len(sgRNA_control_pair_index))]))

        sgRNA_single_gene_index = np.array(sorted(list(set(range(counts_ref.shape[0])).difference(set(np.concatenate((sgRNA_true_pair_index, sgRNA_control_pair_index)))))))
        print(' '.join(["Number of singles:", str(len(sgRNA_single_gene_index))]))

        if (len(sgRNA_single_gene_index) + len(sgRNA_control_pair_index) + len(sgRNA_true_pair_index)) != counts_ref.shape[0]:
            print('Missing annotation')
            print(counts_ref.shape[0] - (len(sgRNA_single_gene_index) + len(sgRNA_control_pair_index) + len(sgRNA_true_pair_index)))

        counts_ref['target_type'] = 'N/A'
        counts_ref['target_type'].iloc[sgRNA_true_pair_index] = 'Dual'
        counts_ref['target_type'].iloc[sgRNA_control_pair_index] = 'Control'
        counts_ref['target_type'].iloc[sgRNA_single_gene_index] = 'Single'


        if 'Type' in counts_ref.columns:
            counts_ref = counts_ref.drop(columns = ['Type'])
        if 'Sequencing' in counts_ref.columns:
            counts_ref = counts_ref.drop(columns = ['Sequencing'])


        ## seperate the replicate counts across T0 and TEnd
        counts_ref['T0_counts'] = ""
        counts_ref['T0_replicate_names'] = ""
        counts_ref['TEnd_counts'] = ""
        counts_ref['TEnd_replicate_names'] = ""

        if isinstance(study_conditions, dict):
            # for different cell_line_origins within a study

            for cell_line_origin in study_conditions:
                curr_conditions = study_conditions[cell_line_origin]

                # access the cell_line_origin counts
                access_level = counts_ref.loc[counts_ref['cell_line_origin'] == cell_line_origin].copy()

                ## get all conditions
                condition = access_level['study_conditions'].value_counts().index.tolist()
                condition = condition[0].split(';')

                # time point T_0
                t_0_index = np.array([i for i in range(len(condition)) if condition[i] in curr_conditions[0]])
                # time point T_end
                t_end_index = np.array([i for i in range(len(condition)) if condition[i] in curr_conditions[1]])

                # get counts
                replicate_sep =  access_level["count_replicates"].apply(    
                    lambda x: np.array(x.split(";"), dtype = np.float64)
                )

                # get time point 
                t_0_comb = replicate_sep.apply(
                    lambda x: check_repeated_constructs(x, t_0_index)
                ).apply(lambda x: ';'.join(x.astype(np.str_)))

                t_end_comb = replicate_sep.apply(
                    lambda x: check_repeated_constructs(x, t_end_index)#np.median(x[t_end_index])
                ).apply(lambda x: ';'.join(x.astype(np.str_)))

                access_level['T0_counts'] = t_0_comb
                access_level['T0_replicate_names'] = ';'.join(curr_conditions[0])
                access_level['TEnd_counts'] = t_end_comb
                access_level['TEnd_replicate_names'] = ';'.join(curr_conditions[1])

                counts_ref.loc[counts_ref['cell_line_origin'] == cell_line_origin] = access_level
        else:
            # for only one cell_line_origin
            curr_conditions = study_conditions

            ## get all conditions
            condition = counts_ref['study_conditions'].value_counts().index.tolist()
            condition = condition[0].split(';')

            # time point T_0
            t_0_index = np.array([i for i in range(len(condition)) if condition[i] in curr_conditions[0]])
            # time point T_end
            t_end_index = np.array([i for i in range(len(condition)) if condition[i] in curr_conditions[1]])

            # get counts
            replicate_sep =  counts_ref["count_replicates"].apply(    
                lambda x: np.array(x.split(";"), dtype = np.float64)
            )

            # get time point 
            t_0_comb = replicate_sep.apply(
                lambda x: check_repeated_constructs(x, t_0_index)
            ).apply(lambda x: ';'.join(x.astype(np.str_)))

            t_end_comb = replicate_sep.apply(
                lambda x: check_repeated_constructs(x, t_end_index)#np.median(x[t_end_index])
            ).apply(lambda x: ';'.join(x.astype(np.str_)))

            counts_ref['T0_counts'] = t_0_comb
            counts_ref['T0_replicate_names'] = ';'.join(curr_conditions[0])
            counts_ref['TEnd_counts'] = t_end_comb
            counts_ref['TEnd_replicate_names'] = ';'.join(curr_conditions[1])


        # proceed to add the orientation

        unsorted_orientations = np.array(['|'.join([counts_ref["gene_1"].iloc[i], counts_ref["gene_2"].iloc[i]]) for i in range(counts_ref.shape[0])])
        sorted_orientations = np.array(['|'.join(sorted([counts_ref["gene_1"].iloc[i], counts_ref["gene_2"].iloc[i]])) for i in range(counts_ref.shape[0])])

        counts_ref['gene_pair'] = sorted_orientations
        counts_ref['gene_pair_orientation'] = 'A_B'
        counts_ref['gene_pair_orientation'].loc[sorted_orientations != unsorted_orientations] = 'B_A'
        
        counts_ref_list = []
        # applied for HORLBECK
        if remove_unrelated_counts:
            print('remove_unrelated_counts = TRUE')
            counts_ref_list = []
            for cl in sorted(list(set(counts_ref['cell_line_origin_origin']))):
                print('For cl = ' + cl)
                temp_counts = counts_ref.loc[counts_ref['cell_line_origin_origin'] == cl]
                temp_scores = score_ref.loc[score_ref['cell_line_origin_origin'] == cl]
                
                # get score genes + control
                score_genes = list(set(temp_scores['gene_1'].tolist() + temp_scores['gene_2'].tolist())) + study_controls
                if study_controls is not None:
                    score_genes = score_genes + study_controls
                # filter down
                temp_counts_filt = temp_counts.loc[(temp_counts['gene_1'].isin(score_genes) & temp_counts['gene_2'].isin(score_genes)) | (temp_counts['target_type'] == 'Control')]
                
                removed = temp_counts.shape[0] - temp_counts_filt.shape[0]
                if removed != 0:
                    print('Removed a total of {rem} sgRNAs...'.format(rem = removed))
                else:
                    print('No unrelated counts found!')
                    
                counts_ref_list.append(temp_counts_filt)
                
            counts_ref = pd.concat(counts_ref_list, axis = 0)
            counts_ref.reset_index(drop = True, inplace = True)

    ################################# counts ref - DONE
    
    print('Sequence reference...')
    if sequence_ref is not None:
        for col in ['sgRNA_guide_name', 'sgRNA_guide_seq', 'sgRNA_target_name']:
            sequence_ref[col] = [i.upper() for i in sequence_ref[col]]
            
        # set the target names to control
        control_idx = np.array([True if str(sequence_ref['sgRNA_target_name'].iloc[i]) in study_controls else False for i in range(sequence_ref.shape[0])]) | np.array([True if str(sequence_ref['sgRNA_guide_name'].iloc[i]) in study_controls else False for i in range(sequence_ref.shape[0])])
        sequence_ref.loc[control_idx, 'sgRNA_target_name'] = 'CONTROL'    
        # add study origin as well
        sequence_ref['study_origin'] = [counts_ref['study_origin'].iloc[0]] * sequence_ref.shape[0]
        
        sequence_ref.reset_index(drop=True, inplace = True)
    
    ################################# sequence ref - DONE
    
    print('Done! Returning...')
    return({'sequence_ref': sequence_ref,
            'counts_ref': counts_ref,
            'score_ref': score_ref})


def insert_study_to_db(engine_link, db_inserts):
    '''
    Inserts the counts to the designated DB.

    **Params**:

    * SLKB_engine: SQLAlchemy engine link
    * db_inserts: Processed data, obtained via ```prepare_study_for_export```

    **Returns**:

    * None

    '''
    # first, get the metadata
    db_metadata = sqlalchemy.MetaData()
    db_metadata.reflect(bind=engine_link)

    # access the tables
    sequence_table = db_metadata.tables['cdko_experiment_design']
    counts_table = db_metadata.tables['cdko_sgrna_counts']
    scores_table = db_metadata.tables['cdko_original_sl_results']

    # then, start the session
    engine_session = sessionmaker(bind=engine_link)
    curr_session = engine_session()

    # get number of records in each table
    sequence_records_num = curr_session.query(sequence_table).count()
    counts_records_num = curr_session.query(counts_table).count()
    scores_records_num = curr_session.query(scores_table).count()

    # get available gene pairs
    available_gene_pairs = curr_session.query(sqlalchemy.func.max(counts_table.c.gene_pair_id)).first()[0]

    # for the first instance
    if available_gene_pairs is None:
        available_gene_pairs = -1

    # in the case of studies that doesn't have counts, but were added
    temp = curr_session.query(sqlalchemy.func.max(scores_table.c.gene_pair_id)).first()[0]

    if temp is not None:
        available_gene_pairs = max(available_gene_pairs, temp)

    # proceed to add the IDs to each table and reindex
    if db_inserts['sequence_ref'] is not None:
        sequence_insert = db_inserts['sequence_ref'].reset_index(drop=True)
    else:
        sequence_insert = None

    if db_inserts['counts_ref'] is not None:
        counts_insert = db_inserts['counts_ref'].reset_index(drop=True)
    else:
        counts_insert = None
    score_insert = db_inserts['score_ref'].reset_index(drop=True)

    # add the current number of records for proper insertion
    if sequence_insert is not None:
        sequence_insert.index += sequence_records_num
    if counts_insert is not None:
        counts_insert.index += counts_records_num
    score_insert.index += scores_records_num

    # set IDs
    if sequence_insert is not None:
        sequence_insert['sgRNA_id'] = sequence_insert.index
    if counts_insert is not None:
        counts_insert['sgRNA_pair_id'] = counts_insert.index
    score_insert['id'] = score_insert.index
    score_insert['gene_pair_id'] = score_insert['id'].copy()
    
    # update the gene pairs
    print('Updating gene pairs with seperator |...')
    if counts_insert is not None:
        counts_insert['gene_pair'] = np.array(['|'.join(sorted([counts_insert["gene_1"].iloc[i], counts_insert["gene_2"].iloc[i]])) for i in range(counts_insert.shape[0])])
    score_insert['gene_pair'] = np.array(['|'.join(sorted([score_insert["gene_1"].iloc[i], score_insert["gene_2"].iloc[i]])) for i in range(score_insert.shape[0])])

    if (sequence_insert is not None) or (counts_insert is not None):
        for_merging = sequence_insert.copy()
        for_merging['ref_id'] = for_merging.index

        # # add the foreign keys
        counts_insert['FK_guide_1_id'] = counts_insert.merge(for_merging, how = 'left', left_on = 'guide_1', right_on = 'sgRNA_guide_name')['ref_id'].values
        counts_insert['FK_guide_2_id'] = counts_insert.merge(for_merging, how = 'left', left_on = 'guide_2', right_on = 'sgRNA_guide_name')['ref_id'].values

    if counts_insert is not None:
        for_merging = score_insert.copy()

        # in the case of multiple cell lines
        for_merging['gene_pair+cell_line+study_origin'] = for_merging['gene_pair'] + '+' + for_merging['cell_line_origin'] + '+' + for_merging['study_origin']
        counts_insert['gene_pair+cell_line+study_origin'] = counts_insert['gene_pair'] + '+' + counts_insert['cell_line_origin'] + '+' + counts_insert['study_origin']

        # set the gene pair id
        #counts_insert['gene_pair_id_all'] = np.NaN
        # dual have the gene pair id
        #counts_insert.loc[counts_insert['target_type'] == 'Dual', 'gene_pair_id_all'] = counts_insert.loc[counts_insert['target_type'] == 'Dual'].groupby(['gene_pair+cell_line+study_origin']).ngroup() + (available_gene_pairs + 1)
        counts_insert['gene_pair_id_all'] = counts_insert.groupby(['gene_pair+cell_line+study_origin']).ngroup() + (available_gene_pairs + 1)

        score_insert['gene_pair_id'] = for_merging.merge(counts_insert.drop_duplicates(subset = 'gene_pair+cell_line+study_origin'), how = 'left', left_on = 'gene_pair+cell_line+study_origin', right_on = 'gene_pair+cell_line+study_origin')['gene_pair_id_all'].values

    ## check if there is any NA in the references
    if (sequence_insert is not None) or (counts_insert is not None):
        for col in ['FK_guide_1_id', 'FK_guide_2_id']:
            if counts_insert[col].isna().sum() > 0:
                print('NA in foreign keys: ' + col)
    else:
        print('No counts and sequences together')

    print('Final QC...')
    # Final quality control
    if sequence_insert is not None:
        sequence_insert = sequence_insert.applymap(lambda x: x.strip() if isinstance(x, str) else x, na_action='ignore')
    if counts_insert is not None:
        counts_insert = counts_insert.applymap(lambda x: x.strip() if isinstance(x, str) else x, na_action='ignore')
    score_insert = score_insert.applymap(lambda x: x.strip() if isinstance(x, str) else x, na_action='ignore')

    # proceed to insert to the database

    # start the transaction, # insert only the columns we need
    with engine_link.begin() as transaction:
        # insert sequence
        print('Beginning transaction...')

        if sequence_insert is not None:
            sequence_insert = sequence_insert.loc[:,['sgRNA_guide_name', 'sgRNA_guide_seq', 'sgRNA_target_name', 'study_origin', 'sgRNA_id']]
            sequence_insert.to_sql(name = 'cdko_experiment_design', con = transaction, if_exists = 'append', index = False, index_label = 'sgRNA_id')

            print('Done sequence')

        # insert CDKO counts
        if counts_insert is not None:
            counts_insert = counts_insert.loc[:,['sgRNA_pair_id', 'FK_guide_1_id', 'FK_guide_2_id', 'gene_pair_id_all', 'gene_pair_orientation', 'T0_counts', 'T0_replicate_names', 'TEnd_counts', 'TEnd_replicate_names', 'target_type', 'study_origin', 'cell_line_origin']]

            counts_insert = counts_insert.rename(columns = {'FK_guide_1_id': 'guide_1_id',
                                                                   'FK_guide_2_id': 'guide_2_id',
                                                                   'gene_pair_id_all': 'gene_pair_id'})

            counts_insert.to_sql(name = 'cdko_sgrna_counts', con = transaction, if_exists = 'append', index = False, index_label = 'sgRNA_pair_id')

            print('Done counts')

        # finally, insert scores
        score_insert = score_insert.loc[:, ['gene_1', 'gene_2', 'study_origin', 'cell_line_origin', 'SL_score', 'SL_score_cutoff', 'statistical_score', 'statistical_score_cutoff', 'gene_pair', 'SL_or_not', 'gene_pair_id', 'id']]
        score_insert.to_sql(name = 'cdko_original_sl_results', con = transaction, if_exists = 'append', index = False, index_label = 'id')

        print('Done score')

        print('Successfully inserted!')

        print('Added Record stats...')
        if sequence_insert is not None:
            print(' '.join(['Sequence insert:', str(sequence_insert.shape[0])]))
        if counts_insert is not None:
            print(' '.join(['Counts insert:', str(counts_insert.shape[0])]))
        print(' '.join(['Score insert:', str(score_insert.shape[0])]))

    print('Done!')

###### Score Analysis Functions

def get_raw_counts(curr_counts):
    '''
    Helper function, gets the raw counts based on the T0 and TEnd annotations of the sample names
    '''
    print('Getting raw counts...')
    # get counts
    T0_counts = curr_counts['T0_counts'].apply(    
        lambda x: np.array(x.split(";"), dtype = np.float64)
    )

    T0_counts = pd.DataFrame(data = T0_counts.tolist(),
                   index = T0_counts.index, columns = curr_counts['T0_replicate_names'].iloc[0].split(';'))

    TEnd_counts = curr_counts['TEnd_counts'].apply(    
        lambda x: np.array(x.split(";"), dtype = np.float64)
    )

    TEnd_counts = pd.DataFrame(data = TEnd_counts.tolist(),
                   index = TEnd_counts.index, columns = curr_counts['TEnd_replicate_names'].iloc[0].split(';'))
    
    # make sure no columns are filled with NAs completely (in case of additional annotations)
    NA_replicate = T0_counts.isna().sum()
    if (NA_replicate == T0_counts.shape[0]).sum() > 0:
        print('Removing NA replicate from T0...')
        T0_counts.drop(NA_replicate.index[NA_replicate == T0_counts.shape[0]], axis = 1, inplace = True)
    
    NA_replicate = TEnd_counts.isna().sum()
    if (NA_replicate == TEnd_counts.shape[0]).sum() > 0:
        print('Removing NA replicate from TEnd...')
        TEnd_counts.drop(NA_replicate.index[NA_replicate == TEnd_counts.shape[0]], axis = 1, inplace = True)

    T0_counts = T0_counts.fillna(0)
    TEnd_counts = TEnd_counts.fillna(0)
    
    return((T0_counts, TEnd_counts))

def filter_counts(curr_counts, filtering_counts = 35):
    '''
    Helper function, filters sgRNAs with counts less than the threshold
    '''
    print(' '.join(["Filtering enabled... Condition:", str(filtering_counts), "counts"]))
    
    curr_counts[curr_counts < filtering_counts] = np.nan
    
    # drop the entire sgRNAs
    curr_counts = curr_counts.dropna()    
    
    return(curr_counts)

def normalize_counts(curr_counts, set_normalization = 1e6):
    
    print("Normalization enabled...")
    
    if set_normalization is not None:
        print("Current counts:")
        print(curr_counts.sum(axis = 0))
        
        norm_value = set_normalization
        print(' '.join(["Normalize based on a specific value...", str(set_normalization), "counts"]))
    else:
        print("Normalize based on sample counts... Current counts:")
        print(curr_counts.sum(axis = 0))
        
        norm_value = np.median(curr_counts.sum(axis = 0))
        print(' '.join(["Normalize value...", str(norm_value), "counts"]))    
    
    filt_locations = curr_counts.isna()
    #print(curr_counts.isna().sum())
    
    curr_counts = (curr_counts * norm_value) / curr_counts.sum(axis = 0)
    curr_counts[filt_locations] = np.nan
    
    return(curr_counts)


def sort_pairs_and_guides(curr_counts):
    # sort the genes and guides based on gene ordering
    print('Sorting gene pairs and guides based on ordering gene ordering...')
    gene_pairs = []
    gene_pair_guides = []
    for i in range(curr_counts.shape[0]):

        guide_1 = curr_counts['sgRNA_guide_name_g1'].iloc[i]
        guide_2 = curr_counts['sgRNA_guide_name_g2'].iloc[i]

        gene_1 = curr_counts['sgRNA_target_name_g1'].iloc[i]
        gene_2 = curr_counts['sgRNA_target_name_g2'].iloc[i]

        t_gene_1, t_gene_2 = sorted([gene_1, gene_2])


        if (t_gene_1 == gene_1) and (t_gene_2 == gene_2):
            gene_1 = t_gene_1
            gene_2 = t_gene_2
        else:
            gene_1 = t_gene_1
            gene_2 = t_gene_2

            # swap the guides accordingly
            temp = guide_1
            guide_1 = guide_2
            guide_2 = temp

        gene_pairs.append('|'.join([gene_1, gene_2]))
        gene_pair_guides.append('|'.join([guide_1, guide_2]))

    return(gene_pairs, gene_pair_guides)

# taken from Horlbeck et al., https://github.com/mhorlbeck/GImap_tools/blob/601cd22126432edadb30202e952859195c73a841/GImap_analysis.py
def quadFitForceIntercept(xdata, ydata, bdata):
    m1 = optimize.fmin(lambda m, x, y: ((m[0]*(x**2) + m[1]*x + bdata - y)**2).sum(), x0=[0.1,0.1], args=(xdata, ydata), disp=0)
    
    return lambda x1: m1[0]*(np.array(x1)**2) + m1[1]*np.array(x1) + bdata


def run_horlbeck_preprocessing(curr_counts, filterThreshold = 35, pseudocount = 10):
        
    T0_counts, TEnd_counts = get_raw_counts(curr_counts.copy())
    
    # horlbeck uses single x single as double, proceed to move them to dual instead
    replace_idx = (curr_counts['target_type'] == 'Single') & (curr_counts['sgRNA_target_name_g1'] == curr_counts['sgRNA_target_name_g2'])
    curr_counts.loc[replace_idx, 'target_type'] = 'Dual'

    if T0_counts.shape[1] != TEnd_counts.shape[1]:
        print("Mismatch times, averaging...")

        T0_counts = pd.DataFrame(data = T0_counts.apply(lambda x: np.mean(x), axis = 1).values,
                             index = T0_counts.index)

        TEnd_counts = pd.DataFrame(data = TEnd_counts.apply(lambda x: np.mean(x), axis = 1).values,
                 index = TEnd_counts.index)

    T0_counts = pd.concat([T0_counts, curr_counts['sgRNA_guide_name_g1'], curr_counts['sgRNA_guide_name_g2']], axis = 1)
    TEnd_counts = pd.concat([TEnd_counts, curr_counts['sgRNA_guide_name_g1'], curr_counts['sgRNA_guide_name_g2']], axis = 1)
    all_sgRNAs = set(TEnd_counts['sgRNA_guide_name_g1']).union(set(TEnd_counts['sgRNA_guide_name_g2']))

    # add sorted targets
    sorted_gene_pairs, sorted_gene_guides = sort_pairs_and_guides(curr_counts.copy())
    curr_counts['sgRNA_pair'] = sorted_gene_guides
    curr_counts['gene_pair'] = sorted_gene_pairs
    
    replicate_list = []
    for replicate_i in range(len(T0_counts.columns)-2):
        print("For replicate " + str(replicate_i + 1))
        meanCounts = pd.concat((TEnd_counts.iloc[:,replicate_i].groupby(TEnd_counts['sgRNA_guide_name_g1']).agg(np.median),TEnd_counts.iloc[:,replicate_i].groupby(TEnd_counts['sgRNA_guide_name_g2']).agg(np.median)),axis=1, keys=['sgRNA_guide_name_g1', 'sgRNA_guide_name_g2'])
        sgsToFilter = set(meanCounts.loc[meanCounts.loc[:,'sgRNA_guide_name_g1'] < filterThreshold].index).union(set(meanCounts.loc[meanCounts.loc[:,'sgRNA_guide_name_g2'] < filterThreshold].index))
        print(" ".join(["Total of", str(len(sgsToFilter)), 'sgRNAs were filtered out of', str(len(all_sgRNAs))]))

        chosen_idx = np.array([True if i not in sgsToFilter else False for i in TEnd_counts['sgRNA_guide_name_g1']]) & np.array([True if i not in sgsToFilter else False for i in TEnd_counts['sgRNA_guide_name_g2']])
        TEnd_counts_curr = TEnd_counts.iloc[chosen_idx, replicate_i]
        T0_counts_curr = T0_counts.iloc[chosen_idx, replicate_i]

        counts_ratio = ((T0_counts_curr + pseudocount).sum()*1.0)/(TEnd_counts_curr + pseudocount).sum()

        # calculate FC like in horlbeck
        replicate_FC = np.log2((TEnd_counts_curr + pseudocount)/(T0_counts_curr + pseudocount)/counts_ratio)
        replicate_FC.columns = ['Replicate_' + str(replicate_i+1) + "_FC"]
        replicate_FC.name = 'Replicate_' + str(replicate_i+1) + "_FC"

        # get control
        control_effect = 0
        if 'Control' in set(curr_counts['target_type']):
            control_index = curr_counts['target_type'] == 'Control'
            if control_index.sum() != 0:
                control_effect = replicate_FC.loc[control_index].median()

        replicate_FC -= control_effect

        # doubling differences, taken from original code
        replicate_FC /= 6.3

        curr_counts = curr_counts.join(replicate_FC)

        replicate_list.append(replicate_FC)

    # save the results to original data
    replicate_list = pd.concat(replicate_list, axis = 1)
    replicate_list = replicate_list.dropna()

    replicate_list = replicate_list.mean(axis = 1)

    replicate_list.columns = ['FC_Averaged']
    replicate_list.name = 'FC_Averaged'

    curr_counts = curr_counts.join(replicate_list)

    average_of_transpose = curr_counts.groupby('sgRNA_pair')['FC_Averaged'].apply(np.nanmean)
    curr_counts = curr_counts.join(average_of_transpose,
                             on = 'sgRNA_pair',
                             rsuffix = "_abbaAveraged")
    
    return(curr_counts)

def run_horlbeck_score(curr_counts, curr_study, curr_cl, do_preprocessing = True, store_loc = os.getcwd(), save_dir = 'HORLBECK_Files', re_run = False):
    '''
    
    Calculates Horlbeck score. Score files will created at the designated store location and save directory. 

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

    '''
    print('Running horlbeck score...')
    
    ######### preprocessing
    
    print('Running preprocessing...')

    if do_preprocessing:
        curr_counts = run_horlbeck_preprocessing(curr_counts)


    #########/ preprocessing
    
    # get save location
    save_loc = os.path.join(store_loc, save_dir, curr_study, curr_cl)
    os.makedirs(save_loc, exist_ok = True)
    
    ######### original horlbeck scoring

    # first, drop the rows with nan replicateFCname
    curr_counts.dropna(subset = ['FC_Averaged_abbaAveraged'], inplace = True)

    # get ab/ba
    a_average, b_average = curr_counts.loc[curr_counts['target_type'] != 'Dual'].copy(), curr_counts.loc[curr_counts['target_type'] != 'Dual'].copy()
    #curr_counts.loc[curr_counts['target_type'] == 'Single'].copy(), curr_counts.loc[curr_counts['target_type'] == 'Single'].copy()
    #
    a_average = a_average[a_average['sgRNA_target_name_g2'] == "CONTROL"]
    b_average = b_average[b_average['sgRNA_target_name_g1'] == "CONTROL"]

    a_average = a_average.groupby('sgRNA_guide_name_g1')['FC_Averaged_abbaAveraged'].apply(np.mean)
    b_average = b_average.groupby('sgRNA_guide_name_g2')['FC_Averaged_abbaAveraged'].apply(np.mean)

    # single, control, and dual phenotypes are used in calculation
    all_pairs = set(curr_counts['sgRNA_guide_name_g1']).union(set(curr_counts['sgRNA_guide_name_g2']))
    curr_counts['GI_Averaged'] = 0

    # for missing pairs, update a_average, b_average
    a_average_0s = list(all_pairs.difference(set(a_average.index)))
    a_average = pd.concat([a_average, pd.Series(data = np.zeros(len(a_average_0s)), index = a_average_0s)])

    b_average_0s = list(all_pairs.difference(set(b_average.index)))
    b_average = pd.concat([b_average, pd.Series(data = np.zeros(len(b_average_0s)), index = b_average_0s)])

    # store in a matrix
    GI_Score_1 = pd.DataFrame(0, index = sorted(list(all_pairs)), columns = sorted(list(all_pairs)))
    GI_Score_2 = pd.DataFrame(0, index = sorted(list(all_pairs)), columns = sorted(list(all_pairs)))
    
    
    # scores have already been computed
    if os.path.exists(os.path.join(save_loc, "GI_Score_1.gzip")) and (not re_run):
        print('Scores exist For GI_Score_1! Loading...')
        GI_Score_1 = pd.read_pickle(os.path.join(save_loc, "GI_Score_1.gzip"))
    else:
        print('Calculating GI_Score_1...')
        
        ## A orientation ()

        ## go through all query sgRNAs
        for query_sgRNA in all_pairs:

            ## get all the pairs with the given query
            idx_loc = (curr_counts['sgRNA_guide_name_g2'] == query_sgRNA)

            if len(idx_loc) == 0:
                continue

            ## all pairs 
            curr_filtered_pairs = curr_counts.loc[idx_loc, :]

            ## get sgRNAs assayed together with the query sgRNA
            selected_sgRNAs = curr_filtered_pairs['sgRNA_guide_name_g1'].values

            if 'Control' in set(curr_counts['target_type']):
                control_sgRNAs = np.where(curr_filtered_pairs['sgRNA_target_name_g1'] == "CONTROL")[0]

            # Fit to a quadratic formula, where the x is the single phenotypes and y is the pair phenotypes

            xs = a_average.loc[selected_sgRNAs].values # a -> b
            ys = curr_filtered_pairs['FC_Averaged_abbaAveraged'].values
            bs = b_average.loc[query_sgRNA] # b -> a

            res_fn = quadFitForceIntercept(xs, ys, bs)

            # get expected
            expected_phenotype = res_fn(xs)

            # the difference is the GI score
            GI_Score = ys - expected_phenotype

            if ('Control' in set(curr_counts['target_type'])) and len(control_sgRNAs) > 0:
                if GI_Score[control_sgRNAs].std() != 0:
                    GI_Score /= GI_Score[control_sgRNAs].std()

            GI_Score_1.loc[query_sgRNA, selected_sgRNAs] = GI_Score
            
        # save scores for future loading
        GI_Score_1.to_pickle(os.path.join(save_loc, "GI_Score_1.gzip"))
    
    if os.path.exists(os.path.join(save_loc, "GI_Score_2.gzip")) and (not re_run):
        print('Scores exist For GI_Score_2! Loading...')
        GI_Score_2 = pd.read_pickle(os.path.join(save_loc, "GI_Score_2.gzip"))
    else:
        print('Calculating GI_Score_2...')

        ## B orientation ()
        ## go through all query sgRNAs
        for query_sgRNA in all_pairs:

            ## get all the pairs with the given query
            idx_loc = (curr_counts['sgRNA_guide_name_g1'] == query_sgRNA)

            if len(idx_loc) == 0:
                continue

            ## all pairs 
            curr_filtered_pairs = curr_counts.loc[idx_loc, :]

            ## get sgRNAs assayed together with the query sgRNA
            selected_sgRNAs = curr_filtered_pairs['sgRNA_guide_name_g2'].values

            if 'Control' in set(curr_counts['target_type']):
                control_sgRNAs = np.where(curr_filtered_pairs['sgRNA_target_name_g2'] == "CONTROL")[0]

            # Fit to a quadratic formula, where the x is the single phenotypes and y is the pair phenotypes

            xs = b_average.loc[selected_sgRNAs].values # b -> a
            ys = curr_filtered_pairs['FC_Averaged_abbaAveraged'].values
            bs = a_average.loc[query_sgRNA] # a -> b

            res_fn = quadFitForceIntercept(xs, ys, bs)

            # get expected
            expected_phenotype = res_fn(xs)

            # the difference is the GI score
            GI_Score = ys - expected_phenotype

            if ('Control' in set(curr_counts['target_type'])) and len(control_sgRNAs) > 0:
                if GI_Score[control_sgRNAs].std() != 0:
                    GI_Score /= GI_Score[control_sgRNAs].std()

            # set the 
            #curr_res['sgRNA_level']['dual'].loc[idx_loc, replicate_GI_name] += GI_Score
            GI_Score_2.loc[query_sgRNA, selected_sgRNAs] = GI_Score


        # save scores for future loading
        GI_Score_2.to_pickle(os.path.join(save_loc, "GI_Score_2.gzip"))
    
    
    # average between A and B orientations
    #curr_res['sgRNA_level']['dual'][replicate_GI_name] /= 2
    GI_Score_avg = (GI_Score_1 + GI_Score_2)/2
    GI_Score_avg = (GI_Score_avg + GI_Score_avg.T)/2

    for i in range(len(curr_counts['GI_Averaged'])):
        guide_1 = curr_counts['sgRNA_guide_name_g1'].iloc[i]
        guide_2 = curr_counts['sgRNA_guide_name_g2'].iloc[i]

        curr_counts['GI_Averaged'].iloc[i] = GI_Score_avg.loc[guide_1, guide_2]

    
    ######### /original horlbeck scoring
    
    
    # store results
    SL_score = curr_counts.groupby('gene_pair')['GI_Averaged'].apply(lambda x: np.mean(x))
    SE = curr_counts.groupby('gene_pair')['GI_Averaged'].apply(lambda x: sem(x, ddof=1))

    genes_1 = [i.split('|')[0] for i in SL_score.index]
    genes_2 = [i.split('|')[1] for i in SL_score.index]
    
    horlbeck_results = pd.DataFrame(data = {'SL_score' : SL_score.values,
                                             'standard_error' : SE.values,
                                             'Gene 1' : genes_1,
                                             'Gene 2' : genes_2}, index = SL_score.index)
    
    
    # remove possible controls
    control_idx = np.array([True if 'CONTROL' in i else False for i in horlbeck_results.index])
    horlbeck_results = horlbeck_results.loc[~control_idx]
    
    results = {}
    results['HORLBECK_SCORE'] = horlbeck_results

    return(results)

def run_median_scores(curr_counts, curr_study, curr_cl, full_normalization = False, re_run = False, store_loc = os.getcwd(), save_dir = 'MEDIAN_Files'):
    '''
    Calculates Median B/NB Scores.

    **Params**:

    * curr_counts: Counts to calculate scores to.)
    * curr_study: String, name of study to analyze data for.
    * curr_cl: String, name of cell line to analyze data for.
    * full_normalization: Whether to normalize counts across the whole sample or according to target type (Default: False)
    * re_run: Boolean. Recreate and rerun the results instead of loading for subsequent analyses (Default: False)
    * store_loc: String: Directory to store the Median files to. (Default: current working directory)
    * save_dir: String: Folder name to store the Median files to. (Default: 'MEDIAN_Files')

    **Returns**:

    * median_res: A dictionary of two pandas dataframes: Median-B and Median-NB.
    '''

    # for standard error
    median_SE_constant = 1.25
    
    print('Running median scores...')
    
    # get save location
    save_loc = os.path.join(store_loc, save_dir, curr_study, curr_cl)
    os.makedirs(save_loc, exist_ok = True)

    
    if os.path.exists(os.path.join(save_loc, "median_results.p")) and (not re_run):
        print('Loading final results!')
        #results =  pd.read_pickle(os.path.join(save_loc, "median_results.p"))
        with open(os.path.join(save_loc, "median_results.p"), 'rb') as handle:
            results = pickle.load(handle)
    else:
    
        ######### preprocessing
        t_0_comb, t_end_comb = get_raw_counts(curr_counts)

        # filter counts, only at T0
        t_0_comb = filter_counts(t_0_comb, filtering_counts = 35)
        print(' '.join(['Filtered a total of', str(t_end_comb.shape[0] - t_0_comb.shape[0]), "out of", str(t_end_comb.shape[0]), "sgRNAs."]))
        print("\n---\n")

        # add pseudocount of 10 after filtering
        t_0_comb = t_0_comb + 10
        t_end_comb = t_end_comb + 10

        # some sgRNAs were filtered out
        overlapping_sgRNAs = sorted(list(set(t_0_comb.index).intersection(set(t_end_comb.index))))

        t_0_comb = t_0_comb.loc[overlapping_sgRNAs,:]
        t_end_comb = t_end_comb.loc[overlapping_sgRNAs,:]
        curr_counts = curr_counts.loc[overlapping_sgRNAs,:]

        # normalize to the median of the all time points
        if full_normalization:
            print('Full normalization...')
            normalization_value = np.median(pd.concat([t_0_comb, t_end_comb], axis = 1).sum(axis = 0))

            t_0_comb = normalize_counts(t_0_comb, set_normalization = normalization_value)
            t_end_comb = normalize_counts(t_end_comb, set_normalization = normalization_value)

        else:
            print('Not full normalization...')
            for subset in set(curr_counts['target_type']):
                idx = curr_counts.loc[curr_counts['target_type'] == subset,:].index

                # normalize to the median of the all time points
                normalization_value = np.median(pd.concat([t_0_comb.loc[idx,:], t_end_comb.loc[idx,:]], axis = 1).sum(axis = 0))

                t_0_comb.loc[idx,:] = normalize_counts(t_0_comb.loc[idx,:], set_normalization = normalization_value)
                t_end_comb.loc[idx,:] = normalize_counts(t_end_comb.loc[idx,:], set_normalization = normalization_value)



        # get median of counts 
        t_0_comb = t_0_comb.apply(lambda x: np.median(x), axis = 1)
        t_end_comb = t_end_comb.apply(lambda x: np.median(x), axis = 1)

        # get LFC
        FC = np.log2(t_end_comb) - np.log2(t_0_comb)

        # set FC
        curr_counts['FC'] = FC

        # add sorted targets
        sorted_gene_pairs, sorted_gene_guides = sort_pairs_and_guides(curr_counts.copy())
        curr_counts['sgRNA_pair'] = sorted_gene_guides
        curr_counts['gene_pair'] = sorted_gene_pairs

        ######### /preprocessing

        # store results
        results = {}
        results['MEDIAN_B_SCORE'] = None
        results['MEDIAN_NB_SCORE'] = None

        ######### scoring

        # get the three target categories
        single = curr_counts.loc[curr_counts['target_type'] == 'Single']
        dual = curr_counts.loc[curr_counts['target_type'] == 'Dual']
        control = curr_counts.loc[curr_counts['target_type'] == 'Control']

        print('Available singles: ' + str(single.shape[0]))
        print('Available duals: ' + str(dual.shape[0]))
        print('Available control: ' + str(control.shape[0]))

        temp_repeat = single.copy()
        temp_repeat['sgRNA_guide_name_g1'] = single["sgRNA_guide_name_g2"]
        temp_repeat['sgRNA_target_name_g1'] = single["sgRNA_target_name_g2"]
        temp_repeat['sgRNA_guide_name_g2'] = single["sgRNA_guide_name_g1"]
        temp_repeat['sgRNA_target_name_g2'] = single["sgRNA_target_name_g1"]

        single_repeat = pd.concat([single, temp_repeat])

        # get single sgRNA impact
        EC_single = single_repeat.groupby("sgRNA_guide_name_g1")['FC'].apply(
                lambda x: np.median(x))

        # get control sgRNA impact
        EC_control = None
        if control.shape[0] != 0:

            temp_repeat = control.copy()
            temp_repeat['sgRNA_guide_name_g1'] = control["sgRNA_guide_name_g2"]
            temp_repeat['sgRNA_guide_name_g2'] = control["sgRNA_guide_name_g1"]

            EC_control = pd.concat([control, temp_repeat]).groupby("sgRNA_guide_name_g1")['FC'].apply(
                lambda x: np.median(x)
            )

            EC_single = EC_single.drop(set(EC_control.index).intersection(set(EC_single.index)))

        # all available dual sgRNAs
        all_pairs = set(dual['sgRNA_guide_name_g1']).union(set(dual['sgRNA_guide_name_g2']))

        # fill for empty
        missing_pairs = np.array(list(all_pairs.difference(set(EC_single.index))))

        print(' '.join(["Filtered single sgRNA count:", str(len(set(missing_pairs)))]))

        # add them as 0s
        EC_single = pd.concat([EC_single, pd.Series(index = missing_pairs, data = np.zeros(len(missing_pairs)))])

        # get EC for each
        EC_1 = EC_single[dual['sgRNA_guide_name_g1']]
        EC_2 = EC_single[dual['sgRNA_guide_name_g2']]

        # calculate Impact Scores (sgRNA level)

        dual['Median-NB-dual-IS'] = dual['FC'].values
        dual['Median-NB-single-IS-Guide-1'] = EC_1.values
        dual['Median-NB-single-IS-Guide-2'] = EC_2.values
        dual['Median-NB-dual-SL-sgRNA'] = dual['FC'].values - EC_1.values - EC_2.values

        ## calculate SL scores (sgRNA)
        gene_pair_SL = dual.groupby('gene_pair')['Median-NB-dual-IS'].apply(lambda x: np.median(x))
        gene_pair_SE = dual.groupby('gene_pair')['Median-NB-dual-IS'].apply(lambda x: np.var(x) / np.size(x))

        ## calculate SL scores (gene)
        gene_SL = single_repeat.groupby("sgRNA_target_name_g1")['FC'].apply(
        lambda x: np.median(x))
        gene_SE = single_repeat.groupby("sgRNA_target_name_g1")['FC'].apply(
            lambda x: np.var(x) / np.size(x))

        genes_1 = np.array([i.split('|')[0] for i in gene_pair_SL.index])
        genes_2 = np.array([i.split('|')[1] for i in gene_pair_SL.index])

        all_genes = set(genes_1).union(set(genes_2))
        missing_genes = all_genes.difference(set(gene_SL.index))
        print(' '.join(["Filtered gene count:", str(len(set(missing_genes)))]))

        # add them as 0s
        gene_SL = pd.concat([gene_SL, pd.Series(index = missing_genes, data = np.zeros(len(missing_genes)))])
        gene_SE = pd.concat([gene_SE, pd.Series(index = missing_genes, data = np.zeros(len(missing_genes)))])

        median_nb_SL = gene_pair_SL.values - gene_SL[genes_1].values - gene_SL[genes_2].values
        median_nb_SE = np.sqrt(gene_pair_SE.values + gene_SE[genes_1].values + gene_SE[genes_2].values) * median_SE_constant
        median_nb_Z = median_nb_SL/median_nb_SE

        median_nb_results = pd.DataFrame(data = {'SL_score' : median_nb_SL,
                                                 'standard_error' : median_nb_SE,
                                                 'Z_SL_score' : median_nb_Z,
                                                 'Gene 1' : genes_1,
                                                 'Gene 2' : genes_2}, index = gene_pair_SL.index)

        results['MEDIAN_NB_SCORE'] = median_nb_results

        if EC_control is not None:
            control_median = np.median(EC_control)

            dual['Median-B-dual-IS'] = dual['FC'].values - control_median
            dual['Median-B-single-IS-Guide-1'] = EC_1.values - control_median
            dual['Median-B-single-IS-Guide-2'] = EC_2.values - control_median
            dual['Median-B-dual-SL-sgRNA'] = (dual['FC'].values - control_median) - (EC_1.values - control_median) - (EC_2.values - control_median)

            ## calculate SL scores (sgRNA)
            gene_pair_SL = dual.groupby('gene_pair')['Median-B-dual-IS'].apply(lambda x: np.median(x))
            gene_pair_SE = dual.groupby('gene_pair')['Median-B-dual-IS'].apply(lambda x: np.var(x) / np.size(x))

            # remove controls first
            single_repeat['FC'] = single_repeat['FC'] - control_median
            ## calculate SL scores (gene)
            gene_SL = single_repeat.groupby("sgRNA_target_name_g1")['FC'].apply(
            lambda x: np.median(x))
            gene_SE = single_repeat.groupby("sgRNA_target_name_g1")['FC'].apply(
                lambda x: np.var(x) / np.size(x))

            genes_1 = np.array([i.split('|')[0] for i in gene_pair_SL.index])
            genes_2 = np.array([i.split('|')[1] for i in gene_pair_SL.index])

            all_genes = set(genes_1).union(set(genes_2))
            missing_genes = all_genes.difference(set(gene_SL.index))
            print(' '.join(["Filtered gene count:", str(len(set(missing_genes)))]))

            # add them as 0s
            gene_SL = pd.concat([gene_SL, pd.Series(index = missing_genes, data = np.zeros(len(missing_genes)))])
            gene_SE = pd.concat([gene_SE, pd.Series(index = missing_genes, data = np.zeros(len(missing_genes)))])

            median_b_SL = gene_pair_SL.values - gene_SL[genes_1].values - gene_SL[genes_2].values
            median_b_SE = np.sqrt(gene_pair_SE.values + gene_SE[genes_1].values + gene_SE[genes_2].values) * median_SE_constant
            median_b_Z = median_b_SL/median_b_SE

            median_b_results = pd.DataFrame(data = {'SL_score' : median_b_SL,
                                                     'standard_error' : median_b_SE,
                                                     'Z_SL_score' : median_b_Z,
                                                     'Gene 1' : genes_1,
                                                     'Gene 2' : genes_2}, index = gene_pair_SL.index)

            results['MEDIAN_B_SCORE'] = median_b_results
            
        # save for easy loading
        with open(os.path.join(save_loc, "median_results.p"), 'wb') as handle:
            pickle.dump(results, handle, protocol=pickle.HIGHEST_PROTOCOL)
    
    ######### /scoring
    
    
    # return computed scores
    return(results)

def run_sgrna_scores(curr_counts, curr_study, curr_cl, full_normalization = False, re_run = False, store_loc = os.getcwd(), save_dir = 'sgRNA-DERIVED_Files'):
    '''
    Calculates sgRNA Derived N/NB scores.

    **Params**:

    * curr_counts: Counts to calculate scores to.)
    * curr_study: String, name of study to analyze data for.
    * curr_cl: String, name of cell line to analyze data for.
    * full_normalization: Whether to normalize counts across the whole sample or according to target type (Default: False)
    * re_run: Boolean. Recreate and rerun the results instead of loading for subsequent analyses (Default: False)
    * store_loc: String: Directory to store the sgRNA-Derived files to. (Default: current working directory)
    * save_dir: String: Folder name to store the sgRNA-Derived files to. (Default: 'sgRNA-DERIVED_Files')


    **Returns**:

    * sgRNA_res: A dictionary of two pandas dataframes: sgRNA_derived_B and sgRNA_derived_NB. 
    '''
    # for standard error
    median_SE_constant = 1.25

    print('Running sgrna derived score...')
    
    # get save location
    save_loc = os.path.join(store_loc, save_dir, curr_study, curr_cl)
    os.makedirs(save_loc, exist_ok = True)

    
    if os.path.exists(os.path.join(save_loc, "sgRNA_results.p")) and (not re_run):
        print('Loading final results!')
        #results =  pd.read_pickle(os.path.join(save_loc, "sgRNA_results.gzip"))
        with open(os.path.join(save_loc, "sgRNA_results.p"), 'rb') as handle:
            results = pickle.load(handle)
    else:

        ######### preprocessing
        t_0_comb, t_end_comb = get_raw_counts(curr_counts)

        # filter counts, only at T0
        t_0_comb = filter_counts(t_0_comb, filtering_counts = 35)
        print(' '.join(['Filtered a total of', str(t_end_comb.shape[0] - t_0_comb.shape[0]), "out of", str(t_end_comb.shape[0]), "sgRNAs."]))
        print("\n---\n")

        # add pseudocount of 10 after filtering
        t_0_comb = t_0_comb + 10
        t_end_comb = t_end_comb + 10

        # some sgRNAs were filtered out
        overlapping_sgRNAs = sorted(list(set(t_0_comb.index).intersection(set(t_end_comb.index))))

        t_0_comb = t_0_comb.loc[overlapping_sgRNAs,:]
        t_end_comb = t_end_comb.loc[overlapping_sgRNAs,:]
        curr_counts = curr_counts.loc[overlapping_sgRNAs,:]

        if full_normalization:
            print('Full normalization...')

            # normalize to the median of the all time points
            normalization_value = np.median(pd.concat([t_0_comb, t_end_comb], axis = 1).sum(axis = 0))

            t_0_comb = normalize_counts(t_0_comb, set_normalization = normalization_value)
            t_end_comb = normalize_counts(t_end_comb, set_normalization = normalization_value)

        else:
            print('Not full normalization...')

            for subset in set(curr_counts['target_type']):
                idx = curr_counts.loc[curr_counts['target_type'] == subset,:].index

                # normalize to the median of the all time points
                normalization_value = np.median(pd.concat([t_0_comb.loc[idx,:], t_end_comb.loc[idx,:]], axis = 1).sum(axis = 0))

                t_0_comb.loc[idx,:] = normalize_counts(t_0_comb.loc[idx,:], set_normalization = normalization_value)
                t_end_comb.loc[idx,:] = normalize_counts(t_end_comb.loc[idx,:], set_normalization = normalization_value)

        # if mismatch, average
        if t_0_comb.shape[1] != t_end_comb.shape[1]:
            print("Mismatch times, averaging...")
            t_0_comb = pd.DataFrame(data = t_0_comb.apply(lambda x: np.mean(x), axis = 1).values,
                         index = t_0_comb.index)
            t_end_comb = pd.DataFrame(data = t_end_comb.apply(lambda x: np.mean(x), axis = 1).values,
                 index = t_end_comb.index)

        # set FC
        curr_counts['FC'] = 0

        # add NOT sorted targets
    #     curr_counts['sgRNA_pair'] = ['|'.join(sorted([curr_counts['sgRNA_guide_name_g1'].iloc[i], curr_counts['sgRNA_guide_name_g2'].iloc[i]])) for i in range(curr_counts.shape[0])]
    #     curr_counts['gene_pair'] = ['|'.join(sorted([curr_counts['sgRNA_target_name_g1'].iloc[i], curr_counts['sgRNA_target_name_g2'].iloc[i]])) for i in range(curr_counts.shape[0])]

        # add sorted targets
        sorted_gene_pairs, sorted_gene_guides = sort_pairs_and_guides(curr_counts.copy())
        curr_counts['sgRNA_pair'] = sorted_gene_guides
        curr_counts['gene_pair'] = sorted_gene_pairs

        # get count annotations
        count_annotations = curr_counts.loc[:,['sgRNA_guide_name_g1', 'sgRNA_guide_name_g2', 'sgRNA_target_name_g1', 'sgRNA_target_name_g2', 'target_type', 'sgRNA_pair', 'gene_pair']].copy()


        ######### /preprocessing

        print('Starting scoring..')

        replicate_results = []
        for i in range(t_0_comb.shape[1]):
            print('calculating for replicate ' + str(i))

            replicate_fc = pd.DataFrame(data = np.log2(t_end_comb.iloc[:, i]/t_0_comb.iloc[:, i]).values,
                                        index = t_end_comb.index,
                                        columns = ['FC'])

            # merge
            replicate_fc = replicate_fc.merge(count_annotations, left_index = True, right_index = True)

            # get the three target categories
            single = replicate_fc.loc[curr_counts['target_type'] == 'Single']
            dual = replicate_fc.loc[curr_counts['target_type'] == 'Dual']
            control = replicate_fc.loc[curr_counts['target_type'] == 'Control']

            ## proceed with GI calculation
            temp_repeat = single.copy()
            temp_repeat['sgRNA_guide_name_g1'] = single["sgRNA_guide_name_g2"]
            temp_repeat['sgRNA_target_name_g1'] = single["sgRNA_target_name_g2"]
            temp_repeat['sgRNA_guide_name_g2'] = single["sgRNA_guide_name_g1"]
            temp_repeat['sgRNA_target_name_g2'] = single["sgRNA_target_name_g1"]

            single_repeat = pd.concat([single, temp_repeat])

            # get single sgRNA impact
            EC_single = single_repeat.groupby("sgRNA_guide_name_g1")['FC'].apply(
                    lambda x: np.median(x))
            sgRNA_SE = single_repeat.groupby("sgRNA_guide_name_g1")['FC'].apply(
                lambda x: median_SE_constant * np.sqrt(np.var(x) / np.size(x)))


            EC_control = None
            if control.shape[0] != 0:# and (study != 'parrish_data')

                temp_repeat = control.copy()
                temp_repeat['sgRNA_guide_name_g1'] = control["sgRNA_guide_name_g2"]
                temp_repeat['sgRNA_guide_name_g2'] = control["sgRNA_guide_name_g1"]

                EC_control = np.median(pd.concat([control, temp_repeat])['FC'])

            ## get all pairs
            all_pairs = set(dual['sgRNA_guide_name_g1']).union(set(dual['sgRNA_guide_name_g2']))

            missing_pairs = np.array(list(all_pairs.difference(set(EC_single.index))))

            print(' '.join(["Filtered single sgRNA count:", str(len(set(missing_pairs)))]))

            # add them as 0s

            EC_single = pd.concat([EC_single, pd.Series(index = missing_pairs, data = np.zeros(len(missing_pairs)))])
            sgRNA_SE = pd.concat([sgRNA_SE, pd.Series(index = missing_pairs, data = np.zeros(len(missing_pairs)))])

            sgRNA_level_scores = dual.groupby(['gene_pair', 'sgRNA_pair'], as_index = False)['FC'].apply(lambda x: np.mean(x))
            sgRNA_level_SE = dual.groupby(['gene_pair', 'sgRNA_pair'], as_index = False)['FC'].apply(lambda x: np.sqrt(np.var(x) / np.size(x)))

            guide_1 = np.array([i.split('|')[0] for i in sgRNA_level_scores['sgRNA_pair']])
            guide_2 = np.array([i.split('|')[1] for i in sgRNA_level_scores['sgRNA_pair']])
            EC_1 = EC_single[guide_1]
            EC_2 = EC_single[guide_2]

            SE_1 = sgRNA_SE[guide_1]
            SE_2 = sgRNA_SE[guide_2]

            sgRNA_level_scores['SL'] = sgRNA_level_scores['FC'].values - EC_1.values - EC_2.values
            sgRNA_level_scores['SE'] = np.sqrt(np.square(sgRNA_level_SE['FC'].values) + np.square(SE_1.values) + np.square(SE_2.values))
            sgRNA_level_scores['SE'].loc[sgRNA_level_scores['SE'].isna()] = 1
            sgRNA_level_scores['SE'].loc[sgRNA_level_scores['SE'] == 0] = 1
            sgRNA_level_scores['Z-Score'] = sgRNA_level_scores['SL'].values/sgRNA_level_scores['SE'].values

            gene_SL_scores_nobackground = sgRNA_level_scores.groupby('gene_pair')['Z-Score'].apply(lambda x: np.median(x))
            gene_SL_scores_SE = sgRNA_level_scores.groupby('gene_pair')['Z-Score'].apply(lambda x:  median_SE_constant * np.sqrt(np.var(x) / np.size(x)))
            gene_SL_scores_SE.loc[gene_SL_scores_SE.isna()] = 1
            gene_SL_scores_SE.loc[gene_SL_scores_SE == 0] = 1
            gene_SL_scores_nobackground_Z = gene_SL_scores_nobackground/gene_SL_scores_SE

            results_nb = pd.concat([gene_SL_scores_nobackground, gene_SL_scores_SE, gene_SL_scores_nobackground_Z], axis = 1)
            results_nb.columns = ['sgRNA-Score-NB_' + str(i), 'sgRNA-Score-NB SE_' + str(i), 'sgRNA-Score-NB SL_' + str(i)]

            replicate_results.append(results_nb)

            if EC_control is not None:

                single['FC'] = single['FC'] - EC_control
                dual['FC'] = dual['FC'] - EC_control

                ## proceed with GI calculation
                temp_repeat = single.copy()
                temp_repeat['sgRNA_guide_name_g1'] = single["sgRNA_guide_name_g2"]
                temp_repeat['sgRNA_target_name_g1'] = single["sgRNA_target_name_g2"]
                temp_repeat['sgRNA_guide_name_g2'] = single["sgRNA_guide_name_g1"]
                temp_repeat['sgRNA_target_name_g2'] = single["sgRNA_target_name_g1"]

                single_repeat = pd.concat([single, temp_repeat])

                # get single sgRNA impact
                EC_single = single_repeat.groupby("sgRNA_guide_name_g1")['FC'].apply(
                        lambda x: np.median(x))
                sgRNA_SE = single_repeat.groupby("sgRNA_guide_name_g1")['FC'].apply(
                    lambda x: median_SE_constant * np.sqrt(np.var(x) / np.size(x)))


                EC_control = None
                if control.shape[0] != 0:# and (study != 'parrish_data')

                    temp_repeat = control.copy()
                    temp_repeat['sgRNA_guide_name_g1'] = control["sgRNA_guide_name_g2"]
                    temp_repeat['sgRNA_guide_name_g2'] = control["sgRNA_guide_name_g1"]

                    EC_control = np.median(pd.concat([control, temp_repeat])['FC'])

                ## get all pairs
                all_pairs = set(dual['sgRNA_guide_name_g1']).union(set(dual['sgRNA_guide_name_g2']))

                missing_pairs = np.array(list(all_pairs.difference(set(EC_single.index))))

                print(' '.join(["Filtered single sgRNA count:", str(len(set(missing_pairs)))]))

                # add them as 0s

                EC_single = pd.concat([EC_single, pd.Series(index = missing_pairs, data = np.zeros(len(missing_pairs)))])
                sgRNA_SE = pd.concat([sgRNA_SE, pd.Series(index = missing_pairs, data = np.zeros(len(missing_pairs)))])

                sgRNA_level_scores = dual.groupby(['gene_pair', 'sgRNA_pair'], as_index = False)['FC'].apply(lambda x: np.mean(x))
                sgRNA_level_SE = dual.groupby(['gene_pair', 'sgRNA_pair'], as_index = False)['FC'].apply(lambda x: np.sqrt(np.var(x) / np.size(x)))

                guide_1 = np.array([i.split('|')[0] for i in sgRNA_level_scores['sgRNA_pair']])
                guide_2 = np.array([i.split('|')[1] for i in sgRNA_level_scores['sgRNA_pair']])
                EC_1 = EC_single[guide_1]
                EC_2 = EC_single[guide_2]

                SE_1 = sgRNA_SE[guide_1]
                SE_2 = sgRNA_SE[guide_2]

                sgRNA_level_scores['SL'] = sgRNA_level_scores['FC'].values - EC_1.values - EC_2.values
                sgRNA_level_scores['SE'] = np.sqrt(np.square(sgRNA_level_SE['FC'].values) + np.square(SE_1.values) + np.square(SE_2.values))
                sgRNA_level_scores['SE'].loc[sgRNA_level_scores['SE'].isna()] = 1
                sgRNA_level_scores['SE'].loc[sgRNA_level_scores['SE'] == 0] = 1
                sgRNA_level_scores['Z-Score'] = sgRNA_level_scores['SL'].values/sgRNA_level_scores['SE'].values

                gene_SL_scores_w_background = sgRNA_level_scores.groupby('gene_pair')['Z-Score'].apply(lambda x: np.median(x))
                gene_SL_scores_SE = sgRNA_level_scores.groupby('gene_pair')['Z-Score'].apply(lambda x:  median_SE_constant * np.sqrt(np.var(x) / np.size(x)))
                gene_SL_scores_SE.loc[gene_SL_scores_SE.isna()] = 1
                gene_SL_scores_SE.loc[gene_SL_scores_SE == 0] = 1
                gene_SL_scores_w_background_Z = gene_SL_scores_w_background/gene_SL_scores_SE


                results_b = pd.concat([gene_SL_scores_w_background, gene_SL_scores_SE, gene_SL_scores_w_background_Z], axis = 1)
                results_b.columns = ['sgRNA-Score-B_' + str(i), 'sgRNA-Score-B SE_' + str(i), 'sgRNA-Score-B SL_' + str(i)]

                replicate_results.append(results_b)

        # save results
        results = {}
        results['SGRNA_DERIVED_NB_SCORE'] = None
        results['SGRNA_DERIVED_B_SCORE'] = None

        merged = pd.concat(replicate_results, axis = 1)

        # sort the names
        merged.index = ['|'.join(sorted(i.split('|'))) for i in merged .index]

        merged['sgRNA-Score_Average_NB'] = merged.loc[:,['sgRNA-Score-NB SL_' + str(i) for i in range(t_end_comb.shape[1])]].mean(axis = 1)
        results['SGRNA_DERIVED_NB_SCORE'] = pd.DataFrame(merged['sgRNA-Score_Average_NB'])
        results['SGRNA_DERIVED_NB_SCORE'].columns = ['SL_score']
        results['SGRNA_DERIVED_NB_SCORE']['Gene 1'] = [i.split('|')[0] for i in results['SGRNA_DERIVED_NB_SCORE'].index]
        results['SGRNA_DERIVED_NB_SCORE']['Gene 2'] = [i.split('|')[1] for i in results['SGRNA_DERIVED_NB_SCORE'].index]

        if 'sgRNA-Score-B_0' in merged.columns:
            merged['sgRNA-Score_Average_B'] = merged.loc[:,['sgRNA-Score-B SL_' + str(i) for i in range(t_end_comb.shape[1])]].mean(axis = 1)
            results['SGRNA_DERIVED_B_SCORE'] = pd.DataFrame(merged['sgRNA-Score_Average_B'])
            results['SGRNA_DERIVED_B_SCORE'].columns = ['SL_score']
            results['SGRNA_DERIVED_B_SCORE']['Gene 1'] = [i.split('|')[0] for i in results['SGRNA_DERIVED_B_SCORE'].index]
            results['SGRNA_DERIVED_B_SCORE']['Gene 2'] = [i.split('|')[1] for i in results['SGRNA_DERIVED_B_SCORE'].index]
            
        # save for easy loading
        with open(os.path.join(save_loc, "sgRNA_results.p"), 'wb') as handle:
            pickle.dump(results, handle, protocol=pickle.HIGHEST_PROTOCOL)

    return(results)




def run_mageck_score(curr_counts, curr_study, curr_cl, store_loc = os.getcwd(), save_dir = 'MAGECK_Files', command_line_params = [], re_run = False):
    '''

    Calculates MAGeCK Score. Score files will created at the designated store location and save directory. 

    **Params**:
    * curr_counts: Counts to calculate scores to.)
    * curr_study: String, name of study to analyze data for.
    * curr_cl: String, name of cell line to analyze data for.
    * store_loc: String: Directory to store the MAGeCK files to. (Default: current working directory)
    * save_dir: String: Folder name to store the MAGeCK files to. (Default: 'MAGECK_Files')
    * command_line_params: Optional list to load programming environment(s) to be able to run mageck tool (i.e. loading path, activating python environment). 
    * re_run: Boolean. Recreate and rerun the results instead of loading for subsequent analyses (Default: False)


    **Returns**:

    * mageck_res: A dict that contains a pandas dataframe for MAGeck Score.
    '''
    print('Running mageck score...')

    # !no preprocessing!
    T0_counts, TEnd_counts = get_raw_counts(curr_counts)

    # due to mageck, don't have any comma on columns
    T0_counts.columns = ['T0_' + str(i) for i in range(T0_counts.shape[1])]
    TEnd_counts.columns = ['TEnd_' + str(i) for i in range(TEnd_counts.shape[1])]

    # get the annotations
    curr_counts['sgRNA_pair'] = ['|'.join([curr_counts['sgRNA_guide_name_g1'].iloc[i], curr_counts['sgRNA_guide_name_g2'].iloc[i]]) for i in range(curr_counts.shape[0])]
    curr_counts['gene_pair'] = ['|'.join([curr_counts['sgRNA_target_name_g1'].iloc[i], curr_counts['sgRNA_target_name_g2'].iloc[i]]) for i in range(curr_counts.shape[0])]

    curr_counts['sgRNA_pair_mageck_id'] = curr_counts['sgRNA_pair'].values + "|" + np.array(range(curr_counts.shape[0]), dtype = str)

    # combine them
    comb = pd.concat([curr_counts.loc[:, ['sgRNA_pair_mageck_id', 'gene_pair']], T0_counts, TEnd_counts], axis = 1)
    comb = comb.fillna(0)
    
    # based on T0 and TEnd counts, we can run them paired
    paired = True if len(T0_counts.columns) == len(TEnd_counts.columns) else False
    print('Paired Status = ' + str(paired))

    ######### save

    # get save location 
    save_loc = os.path.join(store_loc, save_dir, curr_study, curr_cl)
    os.makedirs(save_loc, exist_ok = True)

    # save the counts
    comb.to_csv(os.path.join(save_loc, "counts.csv"), sep = ',', index = False)

    ######### /save

    ######### create script and run

    file_loc = os.path.join(save_loc, 'MAGECK_commands.sh')

    fp = open(file_loc, '+w')
    fp.write("#!/bin/sh\n")

    for line in command_line_params:
        fp.write(line + '\n')
    # get index of last time point columns

    t_end_col_locs = []
    for i in range(comb.shape[1]):
        if comb.columns[i] in TEnd_counts.columns:
            t_end_col_locs.append(str(i-2))

    fp.write("cd \"" + os.path.join(os.getcwd(), save_loc) + "\"\n")

    mageck_control_loc = os.path.join(PACKAGE_PATH, 'files', 'mageck_control.txt')
    if 'Control' in set(curr_counts['target_type']):
        if paired:
            command = "mageck test -k counts.csv -t \"" + ','.join(t_end_col_locs) + "\" --paired --norm-method control --control-gene " + mageck_control_loc + " --normcounts-to-file -n out --pdf-report"
        else:
            command = "mageck test -k counts.csv -t \"" + ','.join(t_end_col_locs) + "\" --norm-method control --control-gene " + mageck_control_loc + " --normcounts-to-file -n out --pdf-report"

    fp.write(command)
    fp.close()

    # set chmod
    os.chmod(file_loc, 0o0777)
    
    # scores have already been computed
    if os.path.exists(os.path.join(save_loc, "out.sgrna_summary.txt")) and (not re_run):
        print('Scores exist!')
    else:
        print("Running mageck...")
        process = subprocess.run([file_loc], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if process.returncode == 1:
            print("Error in mageck!!!")
            return(process.stdout.splitlines())
        else:
            print("Finished running mageck!")


    ######### load results

    print('Loading computed results...')

    res = pd.read_csv(os.path.join(save_loc, "out.sgrna_summary.txt"), index_col = 0, sep = "\t")
    res.index = ['|'.join(i.split('|')[:len(i.split('|'))-1]) for i in res.index]

    gene1 = np.array([genes.split('|')[0] for genes in res['Gene']])
    gene2 = np.array([genes.split('|')[1] for genes in res['Gene']])

    res_df = pd.DataFrame(data = {'Gene Pair': res['Gene'],
                                  'Gene 1' : gene1,
                                  'Gene 2' : gene2,
                                  'MAGECK-FC': res['LFC']})

    # add a column for gene pairs, but sorted
    sorted_pairs = []
    sorted_gene_1 = []
    sorted_gene_2 = []
    for i in range(res_df.shape[0]):

        gene_1 = str(res_df["Gene 1"].iloc[i])
        gene_2 = str(res_df["Gene 2"].iloc[i])

        t_gene_1, t_gene_2 = sorted([gene_1, gene_2])


        if (t_gene_1 == gene_1) and (t_gene_2 == gene_2):
            gene_1 = t_gene_1
            gene_2 = t_gene_2
        else:
            gene_1 = t_gene_1
            gene_2 = t_gene_2

        pair = '|'.join([gene_1, gene_2])

        sorted_gene_1.append(gene_1)
        sorted_gene_2.append(gene_2)

        sorted_pairs.append(pair)
    res_df['Gene Pair'] = sorted_pairs
    res_df['Gene 1'] = sorted_gene_1
    res_df['Gene 2'] = sorted_gene_2

    # get dual controls and remove them
    dual_controls_idx = np.array([True if i == 'CONTROL|CONTROL' else False for i in res_df['Gene Pair']])
    dual_controls = res_df.loc[dual_controls_idx]
    res_df = res_df.loc[~dual_controls_idx]

    # get single targets and remove them
    singles_idx = np.array([True if i == 'CONTROL' else False for i in res_df['Gene 1']]) | np.array([True if i == 'CONTROL' else False for i in res_df['Gene 2']]) | (res_df['Gene 1'] == res_df['Gene 2'])
    singles = res_df.loc[singles_idx]
    res_df = res_df.loc[~singles_idx]

    temp_repeat = singles.copy()
    temp_repeat['Gene 1'] = singles["Gene 2"]
    temp_repeat['Gene 2'] = singles["Gene 1"]

    single_repeat = pd.concat([singles, temp_repeat])

    # now res_df contains strictly dual targets
    ## calculate SL scores
    gene_pair_SL = res_df.groupby('Gene Pair')['MAGECK-FC'].apply(lambda x: np.median(x))
    gene_pair_SE = res_df.groupby('Gene Pair')['MAGECK-FC'].apply(lambda x: np.var(x) / np.size(x))

    gene_SL = single_repeat.groupby("Gene 1")['MAGECK-FC'].apply(
        lambda x: np.median(x))
    gene_SE = single_repeat.groupby("Gene 1")['MAGECK-FC'].apply(
        lambda x: np.var(x) / np.size(x))

    genes_1 = np.array([i.split('|')[0] for i in gene_pair_SL.index])
    genes_2 = np.array([i.split('|')[1] for i in gene_pair_SL.index])

    all_genes = set(genes_1).union(set(genes_2))
    missing_genes = all_genes.difference(set(gene_SL.index))
    print(' '.join(["Filtered gene count:", str(len(set(missing_genes)))]))

    # add them as 0s
    gene_SL = pd.concat([gene_SL, pd.Series(index = missing_genes, data = np.zeros(len(missing_genes)))])
    gene_SE = pd.concat([gene_SE, pd.Series(index = missing_genes, data = np.zeros(len(missing_genes)))])

    mageck_SL = gene_pair_SL.values - gene_SL[genes_1].values - gene_SL[genes_2].values
    mageck_SE = np.sqrt(gene_pair_SE.values + gene_SE[genes_1].values + gene_SE[genes_2].values) * math.sqrt(2)
    mageck_Z = mageck_SL/mageck_SE

    mageck_results = pd.DataFrame(data = {'SL_score' : mageck_SL,
                                             'standard_error' : mageck_SE,
                                             'Z_SL_score' : mageck_Z,
                                             'Gene 1' : genes_1,
                                             'Gene 2' : genes_2}, index = ['|'.join(sorted(i.split('|'))) for i in gene_pair_SL.index])

    ######### /load results

    results = {}
    results['MAGECK_SCORE'] = mageck_results

    return(results)

def run_gemini_score(curr_counts, curr_study, curr_cl, store_loc = os.getcwd(), save_dir = 'GEMINI_Files', command_line_params = [], re_run = False):
    '''
    Calculates GEMINI Score. Score files will created at the designated store location and save directory. 

    **Params**:
    * curr_counts: Counts to calculate scores to.)
    * curr_study: String, name of study to analyze data for.
    * curr_cl: String, name of cell line to analyze data for.
    * store_loc: String: Directory to store the GEMINI files to. (Default: current working directory)
    * save_dir: String: Folder name to store the GEMINI files to. (Default: 'GEMINI_Files')
    * command_line_params: Optional list to load programming environment(s) to be able to run GEMINI through R (i.e. loading path, activating R environment). 
    * re_run: Boolean. Recreate and rerun the results instead of loading for subsequent analyses (Default: False)

    **Returns**:

    * gemini_res: A dict that contains a pandas dataframe for GEMINI Score.
    '''
    print('Running gemini score...')
    
    # !no preprocessing!
    T0_counts, TEnd_counts = get_raw_counts(curr_counts)

    T0_counts.columns = ['T0_' + str(i) for i in range(T0_counts.shape[1])]
    TEnd_counts.columns = ['TEnd_' + str(i) for i in range(TEnd_counts.shape[1])]
    
    # get save location 
    save_loc = os.path.join(store_loc, save_dir, curr_study, curr_cl)
    os.makedirs(save_loc, exist_ok = True)
    
    # save the sequences
    study_sequences = pd.DataFrame({'Guide_ID' : curr_counts['sgRNA_guide_name_g1'].tolist() + curr_counts['sgRNA_guide_name_g2'].tolist(),
                                    'Sequence' : curr_counts['sgRNA_guide_seq_g1'].tolist() + curr_counts['sgRNA_guide_seq_g2'].tolist()})

    study_sequences.drop_duplicates(subset = 'Sequence', inplace = True)
    study_sequences.index = study_sequences['Guide_ID']
    study_sequences.to_csv(os.path.join(save_loc, "sequences.csv"), sep = ',', index = False)

    # save guidexgene annotation
    guide_gene_annotation = pd.DataFrame({'Sequences Comb': curr_counts[['sgRNA_guide_seq_g1', 'sgRNA_guide_seq_g2']].agg(';'.join, axis = 1),
                                          'Gene 1': curr_counts['sgRNA_target_name_g1'],
                                          'Gene 2': curr_counts['sgRNA_target_name_g2']})

    gemini_counts = pd.DataFrame(guide_gene_annotation['Sequences Comb'].copy())
    gemini_counts = gemini_counts.merge(T0_counts, how = 'left', right_index = True, left_index = True)
    gemini_counts = gemini_counts.merge(TEnd_counts, how = 'left', right_index = True, left_index = True)

    guide_gene_annotation.reset_index(drop = True, inplace = True)
    guide_gene_annotation.to_csv(os.path.join(save_loc, "guide_gene_annotation.csv"), sep = ',', index = False)

    # save counts
    gemini_counts.reset_index(drop = True, inplace = True)
    gemini_counts.to_csv(os.path.join(save_loc, "counts.csv"), sep = ',', index = False)
    
    # write a gemini bash file
    file_loc = os.path.join(save_loc, 'GEMINI_commands.sh')
    fp = open(file_loc, '+w')
    fp.write("#!/bin/sh\n")
    for line in command_line_params:
        fp.write(line + '\n')
    fp.write('Rscript --vanilla ' + os.path.join(PACKAGE_PATH, 'files', 'GEMINI.R') + ' --args ' + os.path.join(save_dir, curr_study, curr_cl) + '\n')
    fp.close()

    # set chmod
    os.chmod(file_loc, 0o0777)

    
    #### scoring
    
    # scores have already been computed
    if os.path.exists(os.path.join(save_loc, 'GEMINI_Scores.csv')) and (not re_run):
        print('Scores exist!')
    else:
        print("Running GEMINI...")
        process = subprocess.run([file_loc], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if process.returncode == 1:
            print("Error in GEMINI!!!")
            return(process.stdout.splitlines())
        else:
            print("Finished running GEMINI!")

        
    ######### load results
    
    res = pd.read_csv(os.path.join(save_loc, 'GEMINI_Scores.csv'), index_col = 0)

    res.index = ['|'.join(sorted(i.split(';'))) for i in res.index]
    res.columns = ['GEMINI Score Strong', 'GEMINI Sensitive Lethality', 'GEMINI Sensitive Recovery']

    # # get only dual SL
    only_dual_idx = [False if 'CONTROL' in i else True for i in res.index]
    res = res.loc[only_dual_idx]

    # set results
    gemini_results = pd.DataFrame(data = {'SL_score_Strong' : res['GEMINI Score Strong'].values,
                                          'SL_score_SensitiveLethality' : res['GEMINI Sensitive Lethality'].values,
                                          'SL_score_SensitiveRecovery' : res['GEMINI Sensitive Recovery'].values,
                                             'Gene 1' : [i.split('|')[0] for i in res.index],
                                             'Gene 2' : [i.split('|')[1] for i in res.index]}, 
                                  index = ['|'.join(sorted(i.split('|'))) for i in res.index])

    
    
    ######### /load results

    results = {}
    results['GEMINI_SCORE'] = gemini_results
    
    return(results)

###### Adding Scores to Database Functions

def add_table_to_db(curr_counts, curr_results, table_name, engine_link):
    
    print('---------ADDING-TO-DB---------')
    
    # print table
    print('Processing table for: ' + table_name)
    
    # add sorted targets
    # add a sorted gene pair column
    curr_counts['gene_pair'] = ['|'.join(sorted([curr_counts['sgRNA_target_name_g1'].iloc[i], curr_counts['sgRNA_target_name_g2'].iloc[i]])) for i in range(curr_counts.shape[0])]

    # remove the same ones
    curr_results = curr_results.loc[curr_results['Gene 1'] != curr_results['Gene 2'],:]

    # keep only score columns
    curr_results.drop(['Gene 1', 'Gene 2'], axis = 1, inplace = True, errors = 'ignore')

    # make the index
    # curr_counts['gene_pair'] = ['|'.join(sorted([curr_counts['sgRNA_target_name_g1'].iloc[i], curr_counts['sgRNA_target_name_g2'].iloc[i]])) for i in range(curr_counts.shape[0])]

    # merge and get final table
    curr_results = curr_results.merge(curr_counts.drop_duplicates(subset = 'gene_pair'), how = 'left', left_index = True, right_on ='gene_pair').loc[:, ['gene_pair_id'] + list(curr_results.columns)]

    # first, get the metadata
    db_metadata = sqlalchemy.MetaData()
    db_metadata.reflect(bind=engine_link)

    # access the tables
    curr_table = db_metadata.tables[table_name]

    # then, start the session
    engine_session = sessionmaker(bind=engine_link)
    curr_session = engine_session()

    # get number of records in each table
    curr_records_num = curr_session.query(curr_table).count()

    curr_results.reset_index(drop = True, inplace = True)
    curr_results.index += curr_records_num
    # set index
    curr_results['id'] = curr_results.index
    
    if curr_results['gene_pair_id'].isna().sum() > 0:
        print('NA found')
        return()
        
    with engine_link.begin() as transaction:
        # insert sequence
        print('Beginning transaction...')

        # insert scores
        curr_results.to_sql(name = table_name, con = transaction, if_exists = 'append', index = False, index_label = 'id')

        print('Successfully inserted!')

        print('Added Record stats...')
        print(' '.join(['Score insert:', str(curr_results.shape[0])]))



def check_if_added_to_table(curr_counts, table_name, engine_link):
    '''
        
    If running the scoring methods multiple times, the method may be useful in skipping over the computation if there are gene pair records already in the database.

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
    '''
    print('Checking if score already computed: ' + table_name)
    
    # get available results
    res = pd.read_sql_query(con=engine_link.connect(), 
                              sql=sqlalchemy.text('SELECT * from ' + table_name.lower()), index_col = 'id')
    
    if res.shape[0] == 0:
        # none added, so proceed
        return(False)
    else:
        # get all IDs
        available_ids = set(res['gene_pair_id']).difference(set(curr_counts['gene_pair_id']))
        
        # already added
        if len(available_ids) == 0:
            print('Scores already in database!')
            print('Inserted scores: ' + str(len(set(curr_counts['gene_pair_id']))))
            print('---------NOT-TO-DB---------')
            return(True)
        elif len(available_ids) == len(set(res['gene_pair_id'])):
            # none added, so proceed
            return(False)
        else:
            print('Scores already in database.')
            print('Inserted scores: ' + str(res.shape[0]))
            print('---------NOT-TO-DB---------')
            return(True)

###### Score Query Functions

def query_result_table(curr_counts, table_name, curr_study, curr_cl, engine_link):
    '''
    
    Obtain SL Scores from the specified scoring table.

    **Params**:

    * curr_counts: Counts to obtain the scores from.
    * table_name: Must be any of the 7 scoring table names, unless customly added:
        * horlbeck_score
        * median_b_score
        * median_nb_score
        * gemini_score
        * mageck_score
        * sgrna_derived_b_score
        * sgrna_derived_nb_score
    * curr_study: String, name of the study to obtain the results for.
    * curr_cl: String, name of the cell line to obtain the results for.
    * engine_link: SQLAlchemy connection for the database.

    **Returns**:

    * result: A pandas dataframe of the inserted results. Includes annotations for gene pair, study origin, and cell line origin.

    '''
    print('Accessing table: ' + table_name)
    
    # get available results
    res = pd.read_sql_query(con=engine_link.connect(), 
                              sql=sqlalchemy.text('SELECT * from ' + table_name.lower()), index_col = 'id')
    
    # possible gene pairs
    curr_counts['gene_pair'] = ['|'.join(sorted([curr_counts['sgRNA_target_name_g1'].iloc[i], curr_counts['sgRNA_target_name_g2'].iloc[i]])) for i in range(curr_counts.shape[0])]

    # get results
    query_res = curr_counts.loc[curr_counts['target_type'] == 'Dual', ['gene_pair', 'gene_pair_id']].drop_duplicates(subset = ['gene_pair_id'])
    query_res = query_res.merge(res, left_on = 'gene_pair_id', right_on = 'gene_pair_id').drop('gene_pair_id', axis = 1)
    
    # add column names to the front
    names_dict = {i: table_name + '_' + i for i in query_res.columns[1:]}
    query_res.rename(columns = names_dict, inplace = True)
    
    print('Available gene pairs: ' + str(query_res.shape[0]))
    
    # add name of study
    query_res['study_origin'] = curr_study
    
    # add name of cell line
    query_res['cell_line_origin'] = curr_cl
    
    return(query_res)