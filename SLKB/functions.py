# imports
import numpy as np
import pandas as pd
import os
import math
from itertools import chain
import pickle
import sqlalchemy
from sqlalchemy.orm import sessionmaker
from scipy import optimize
from scipy.stats import sem
import subprocess
import shlex
import sqlite3
import sqlalchemy

from sqlalchemy.engine import Engine
from sqlalchemy import event
# to enable foreign keys constraint
@event.listens_for(Engine, "connect")
def set_sqlite_pragma(dbapi_connection, connection_record):
    cursor = dbapi_connection.cursor()
    cursor.execute("PRAGMA foreign_keys=ON")
    cursor.close()

def create_SLKB(location = os.getcwd(), name = 'myCDKO_db', disable_foreign_keys = True):

    # create connection
    conn = sqlite3.connect('myCDKO_db') 

    # get cursor
    c = conn.cursor()

    # turn on foreign keys
    c.execute("PRAGMA foreign_keys=ON")

    ## 1. CDKO Experiment Design
    c.execute('''
            CREATE TABLE IF NOT EXISTS CDKO_EXPERIMENT_DESIGN
            ([sgRNA_id] INTEGER, 
            [sgRNA_guide_name] TEXT NOT NULL,
            [sgRNA_guide_seq] TEXT NOT NULL,
            [sgRNA_target_name] TEXT NOT NULL,
            [study_origin] TEXT NOT NULL,
            PRIMARY KEY (sgRNA_id)
            )
            ''')

    conn.commit()

    ## 2. CDKO sgRNA counts
    c.execute('''
          CREATE TABLE IF NOT EXISTS CDKO_SGRNA_COUNTS
          ([sgRNA_pair_id] INTEGER, 
          [guide_1_id] INTEGER,
          [guide_2_id] INTEGER,
          [gene_pair_id] INTEGER,
          [gene_pair_orientation] TEXT,
          [T0_counts] TEXT,
          [T0_replicate_names] TEXT,
          [TEnd_counts] TEXT,
          [TEnd_replicate_names] TEXT,
          [target_type] TEXT,
          [study_origin] TEXT NOT NULL,
          [cell_line_origin] TEXT NOT NULL,
          PRIMARY KEY (sgRNA_pair_id),
          FOREIGN KEY(guide_1_id) REFERENCES CDKO_EXPERIMENT_DESIGN(sgRNA_id),
          FOREIGN KEY(guide_2_id) REFERENCES CDKO_EXPERIMENT_DESIGN(sgRNA_id)
          
          )
          ''')

    conn.commit()

    ## 3. CDKO Original SL Results
    c.execute('''
          CREATE TABLE IF NOT EXISTS CDKO_ORIGINAL_SL_RESULTS
          ([id] INTEGER,
          [gene_pair_id] INTEGER, 
          [gene_pair] TEXT NOT NULL,
          [study_origin] TEXT NOT NULL,
          [cell_line_origin] TEXT NOT NULL,
          [gene_1] TEXT NOT NULL,
          [gene_2] TEXT NOT NULL,
          [SL_or_not] TEXT NOT NULL,
          [SL_score] REAL,
          [statistical_score] REAL,
          [SL_score_cutoff] REAL,
          [statistical_score_cutoff] REAL,
          PRIMARY KEY (id)
          )
          ''')

    conn.commit()

    ## Horlbeck Score
    c.execute('''
          CREATE TABLE IF NOT EXISTS HORLBECK_SCORE
          ([id] INTEGER,
          [gene_pair_id] INTEGER, 
          [SL_score] REAL,
          [standard_error] REAL,
          PRIMARY KEY (id)
          )
          ''')
    conn.commit()

    ## Median-B Score
    c.execute('''
          CREATE TABLE IF NOT EXISTS MEDIAN_B_SCORE
          ([id] INTEGER,
          [gene_pair_id] INTEGER, 
          [SL_score] REAL,
          [standard_error] REAL,
          [Z_SL_score] REAL,
          PRIMARY KEY (id)
          )
          ''')

    conn.commit()

    ## Median-NB Score
    c.execute('''
          CREATE TABLE IF NOT EXISTS MEDIAN_NB_SCORE
          ([id] INTEGER,
          [gene_pair_id] INTEGER, 
          [SL_score] REAL,
          [standard_error] REAL,
          [Z_SL_score] REAL,
          PRIMARY KEY (id)
          )
          ''')

    conn.commit()

    ## GEMINI Score
    c.execute('''
          CREATE TABLE IF NOT EXISTS GEMINI_SCORE
          ([id] INTEGER,
          [gene_pair_id] INTEGER, 
          [SL_score_Strong] REAL,
          [SL_score_SensitiveLethality] REAL,
          [SL_score_SensitiveRecovery] REAL,
          PRIMARY KEY (id)
          )
          ''')
    conn.commit()

    ## MAGeCK Score
    c.execute('''
          CREATE TABLE IF NOT EXISTS MAGECK_SCORE
          ([id] INTEGER,
          [gene_pair_id] INTEGER, 
          [SL_score] REAL,
          [standard_error] REAL,
          [Z_SL_score] REAL,
          PRIMARY KEY (id)
          )
          ''')
    conn.commit()

    ## sgRNA-Derived B Score
    c.execute('''
          CREATE TABLE IF NOT EXISTS SGRA_DERIVED_B_SCORE
          ([id] INTEGER,
          [gene_pair_id] INTEGER, 
          [SL_score] REAL,
          PRIMARY KEY (id)
          )
          ''')
    conn.commit()

    ## sgRNA-Derived NB Score
    c.execute('''
          CREATE TABLE IF NOT EXISTS SGRA_DERIVED_NB_SCORE
          ([id] INTEGER,
          [gene_pair_id] INTEGER, 
          [SL_score] REAL,
          PRIMARY KEY (id)
          )
          ''')
    conn.commit()

    # add the foreign keys
    if not disable_foreign_keys:
        pass

    ## Save changes
    conn.close()







