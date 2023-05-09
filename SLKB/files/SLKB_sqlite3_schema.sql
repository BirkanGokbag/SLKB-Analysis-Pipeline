DROP TABLE IF EXISTS cdko_experiment_design;
CREATE TABLE cdko_experiment_design
          ([sgRNA_id] INTEGER, 
          [sgRNA_guide_name] TEXT NOT NULL,
          [sgRNA_guide_seq] TEXT NOT NULL,
          [sgRNA_target_name] TEXT NOT NULL,
          [study_origin] TEXT NOT NULL,
          PRIMARY KEY (sgRNA_id)
          );
DROP TABLE IF EXISTS cdko_sgrna_counts;
CREATE TABLE cdko_sgrna_counts
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
          FOREIGN KEY(guide_1_id) REFERENCES cdko_experiment_design(sgRNA_id),
          FOREIGN KEY(guide_2_id) REFERENCES cdko_experiment_design(sgRNA_id)
          
          );
DROP TABLE IF EXISTS cdko_original_sl_results;
CREATE TABLE cdko_original_sl_results
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
          PRIMARY KEY (id),
          FOREIGN KEY(gene_pair_id) REFERENCES cdko_sgrna_counts(gene_pair_id)
          );
DROP TABLE IF EXISTS horlbeck_score;
CREATE TABLE horlbeck_score
          ([id] INTEGER,
          [gene_pair_id] INTEGER, 
          [SL_score] REAL,
          [standard_error] REAL,
          PRIMARY KEY (id),
          FOREIGN KEY(gene_pair_id) REFERENCES cdko_sgrna_counts(gene_pair_id)
          );
DROP TABLE IF EXISTS median_b_score;
CREATE TABLE median_b_score
          ([id] INTEGER,
          [gene_pair_id] INTEGER, 
          [SL_score] REAL,
          [standard_error] REAL,
          [Z_SL_score] REAL,
          PRIMARY KEY (id),
          FOREIGN KEY(gene_pair_id) REFERENCES cdko_sgrna_counts(gene_pair_id)
          );
DROP TABLE IF EXISTS median_nb_score;
CREATE TABLE median_nb_score
          ([id] INTEGER,
          [gene_pair_id] INTEGER, 
          [SL_score] REAL,
          [standard_error] REAL,
          [Z_SL_score] REAL,
          PRIMARY KEY (id),
          FOREIGN KEY(gene_pair_id) REFERENCES cdko_sgrna_counts(gene_pair_id)
          );
DROP TABLE IF EXISTS gemini_score;
CREATE TABLE gemini_score
          ([id] INTEGER,
          [gene_pair_id] INTEGER, 
          [SL_score_Strong] REAL,
          [SL_score_SensitiveLethality] REAL,
          [SL_score_SensitiveRecovery] REAL,
          PRIMARY KEY (id),
          FOREIGN KEY(gene_pair_id) REFERENCES cdko_sgrna_counts(gene_pair_id)
          );
DROP TABLE IF EXISTS mageck_score;
CREATE TABLE mageck_score
          ([id] INTEGER,
          [gene_pair_id] INTEGER, 
          [SL_score] REAL,
          [standard_error] REAL,
          [Z_SL_score] REAL,
          PRIMARY KEY (id),
          FOREIGN KEY(gene_pair_id) REFERENCES cdko_sgrna_counts(gene_pair_id)
          );
DROP TABLE IF EXISTS sgrna_derived_b_score;
CREATE TABLE sgrna_derived_b_score
          ([id] INTEGER,
          [gene_pair_id] INTEGER, 
          [SL_score] REAL,
          PRIMARY KEY (id),
          FOREIGN KEY(gene_pair_id) REFERENCES cdko_sgrna_counts(gene_pair_id)
          );
DROP TABLE IF EXISTS sgrna_derived_nb_score;
CREATE TABLE sgrna_derived_nb_score
          ([id] INTEGER,
          [gene_pair_id] INTEGER, 
          [SL_score] REAL,
          PRIMARY KEY (id),
          FOREIGN KEY(gene_pair_id) REFERENCES cdko_sgrna_counts(gene_pair_id)
          );

DROP VIEW IF EXISTS joined_counts;

CREATE VIEW joined_counts         AS            SELECT d.sgRNA_guide_name as sgRNA_guide_name_g1, d.sgRNA_guide_seq as sgRNA_guide_seq_g1, d.sgRNA_target_name as sgRNA_target_name_g1,                       e.sgRNA_guide_name as sgRNA_guide_name_g2, e.sgRNA_guide_seq as sgRNA_guide_seq_g2, e.sgRNA_target_name as sgRNA_target_name_g2,                                         c.* FROM cdko_sgrna_counts c                                        LEFT JOIN cdko_experiment_design d                                        ON c.guide_1_id = d.sgRNA_id                                         LEFT JOIN cdko_experiment_design e                                        ON c.guide_2_id = e.sgRNA_id
/* joined_counts(sgRNA_guide_name_g1,sgRNA_guide_seq_g1,sgRNA_target_name_g1,sgRNA_guide_name_g2,sgRNA_guide_seq_g2,sgRNA_target_name_g2,sgRNA_pair_id,guide_1_id,guide_2_id,gene_pair_id,gene_pair_orientation,T0_counts,T0_replicate_names,TEnd_counts,TEnd_replicate_names,target_type,study_origin,cell_line_origin) */;
DROP VIEW IF EXISTS calculated_sl_table;

CREATE VIEW calculated_sl_table          AS             SELECT joined_counts.sgRNA_target_name_g1 gene_1,                                   joined_counts.sgRNA_target_name_g2 gene_2,                                   joined_counts.study_origin study_origin,                                   joined_counts.cell_line_origin cell_line_origin,                                   joined_counts.gene_pair_id gene_pair_id,                                   median_nb_score.SL_score median_nb_score_SL_score,                                   median_nb_score.standard_error median_nb_score_standard_error,                                   median_nb_score.Z_SL_score median_nb_score_Z_SL_score,                                   median_b_score.SL_score median_b_score_SL_score,                                   median_b_score.standard_error median_b_score_standard_error,                                   median_b_score.Z_SL_score median_b_score_Z_SL_score,                                   sgrna_derived_b_score.SL_score sgrna_derived_b_score_SL_score,                                   sgrna_derived_nb_score.SL_score sgrna_derived_nb_score_SL_score,                                   horlbeck_score.SL_score horlbeck_score_SL_score,                                   horlbeck_score.standard_error horlbeck_score_standard_error,                                   mageck_score.SL_score mageck_score_SL_score,                                   mageck_score.standard_error mageck_score_standard_error,                                   mageck_score.Z_SL_score mageck_score_Z_SL_score,                                   gemini_score.SL_score_Strong gemini_score_SL_score_Strong,                                   gemini_score.SL_score_SensitiveLethality gemini_score_SL_score_SensitiveLethality,                                   gemini_score.SL_score_SensitiveRecovery gemini_score_SL_score_SensitiveRecovery                                             FROM gemini_score                                             LEFT JOIN joined_counts                                             ON gemini_score.gene_pair_id = joined_counts.gene_pair_id                                             LEFT JOIN median_b_score                                             ON gemini_score.gene_pair_id = median_b_score.gene_pair_id                                             LEFT JOIN sgrna_derived_b_score                                             ON gemini_score.gene_pair_id = sgrna_derived_b_score.gene_pair_id                                             LEFT JOIN sgrna_derived_nb_score                                             ON gemini_score.gene_pair_id = sgrna_derived_nb_score.gene_pair_id                                             LEFT JOIN horlbeck_score                                             ON gemini_score.gene_pair_id = horlbeck_score.gene_pair_id                                             LEFT JOIN mageck_score                                             ON gemini_score.gene_pair_id = mageck_score.gene_pair_id                                             LEFT JOIN median_nb_score                                             ON gemini_score.gene_pair_id = median_nb_score.gene_pair_id                                                 GROUP BY gemini_score.gene_pair_id
/* calculated_sl_table(gene_1,gene_2,study_origin,cell_line_origin,gene_pair_id,median_nb_score_SL_score,median_nb_score_standard_error,median_nb_score_Z_SL_score,median_b_score_SL_score,median_b_score_standard_error,median_b_score_Z_SL_score,sgrna_derived_b_score_SL_score,sgrna_derived_nb_score_SL_score,horlbeck_score_SL_score,horlbeck_score_standard_error,mageck_score_SL_score,mageck_score_standard_error,mageck_score_Z_SL_score,gemini_score_SL_score_Strong,gemini_score_SL_score_SensitiveLethality,gemini_score_SL_score_SensitiveRecovery) */;
