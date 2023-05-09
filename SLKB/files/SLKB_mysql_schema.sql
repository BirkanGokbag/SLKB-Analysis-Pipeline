-- MySQL dump 10.13  Distrib 8.0.33, for macos13 (arm64)
--
-- Host: localhost    Database: SLKB_mysql
-- ------------------------------------------------------
-- Server version	8.0.33

/*!40101 SET @OLD_CHARACTER_SET_CLIENT=@@CHARACTER_SET_CLIENT */;
/*!40101 SET @OLD_CHARACTER_SET_RESULTS=@@CHARACTER_SET_RESULTS */;
/*!40101 SET @OLD_COLLATION_CONNECTION=@@COLLATION_CONNECTION */;
/*!50503 SET NAMES utf8mb4 */;
/*!40103 SET @OLD_TIME_ZONE=@@TIME_ZONE */;
/*!40103 SET TIME_ZONE='+00:00' */;
/*!40014 SET @OLD_UNIQUE_CHECKS=@@UNIQUE_CHECKS, UNIQUE_CHECKS=0 */;
/*!40014 SET @OLD_FOREIGN_KEY_CHECKS=@@FOREIGN_KEY_CHECKS, FOREIGN_KEY_CHECKS=0 */;
/*!40101 SET @OLD_SQL_MODE=@@SQL_MODE, SQL_MODE='NO_AUTO_VALUE_ON_ZERO' */;
/*!40111 SET @OLD_SQL_NOTES=@@SQL_NOTES, SQL_NOTES=0 */;

--
-- Current Database: `SLKB_mysql_live`
--

CREATE DATABASE /*!32312 IF NOT EXISTS*/ `SLKB_mysql_live` /*!40100 DEFAULT CHARACTER SET utf8mb4 COLLATE utf8mb4_0900_ai_ci */ /*!80016 DEFAULT ENCRYPTION='N' */;

USE `SLKB_mysql_live`;

--
-- Table structure for table `cdko_experiment_design`
--

DROP TABLE IF EXISTS `cdko_experiment_design`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `cdko_experiment_design` (
  `sgRNA_id` int NOT NULL AUTO_INCREMENT,
  `sgRNA_guide_name` text COLLATE utf8mb4_general_ci NOT NULL,
  `sgRNA_guide_seq` text COLLATE utf8mb4_general_ci NOT NULL,
  `sgRNA_target_name` text COLLATE utf8mb4_general_ci NOT NULL,
  `study_origin` text COLLATE utf8mb4_general_ci NOT NULL,
  PRIMARY KEY (`sgRNA_id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_general_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `cdko_original_sl_results`
--

DROP TABLE IF EXISTS `cdko_original_sl_results`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `cdko_original_sl_results` (
  `id` int NOT NULL AUTO_INCREMENT,
  `gene_pair_id` int DEFAULT NULL,
  `gene_pair` text COLLATE utf8mb4_general_ci NOT NULL,
  `study_origin` text COLLATE utf8mb4_general_ci NOT NULL,
  `cell_line_origin` text COLLATE utf8mb4_general_ci NOT NULL,
  `gene_1` text COLLATE utf8mb4_general_ci NOT NULL,
  `gene_2` text COLLATE utf8mb4_general_ci NOT NULL,
  `SL_or_not` text COLLATE utf8mb4_general_ci NOT NULL,
  `SL_score` double DEFAULT NULL,
  `statistical_score` double DEFAULT NULL,
  `SL_score_cutoff` double DEFAULT NULL,
  `statistical_score_cutoff` double DEFAULT NULL,
  PRIMARY KEY (`id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_general_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `cdko_sgrna_counts`
--

DROP TABLE IF EXISTS `cdko_sgrna_counts`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `cdko_sgrna_counts` (
  `sgRNA_pair_id` int NOT NULL AUTO_INCREMENT,
  `guide_1_id` int DEFAULT NULL,
  `guide_2_id` int DEFAULT NULL,
  `gene_pair_id` int DEFAULT NULL,
  `gene_pair_orientation` text COLLATE utf8mb4_general_ci,
  `T0_counts` text COLLATE utf8mb4_general_ci,
  `T0_replicate_names` text COLLATE utf8mb4_general_ci,
  `TEnd_counts` text COLLATE utf8mb4_general_ci,
  `TEnd_replicate_names` text COLLATE utf8mb4_general_ci,
  `target_type` text COLLATE utf8mb4_general_ci,
  `study_origin` text COLLATE utf8mb4_general_ci NOT NULL,
  `cell_line_origin` text COLLATE utf8mb4_general_ci NOT NULL,
  PRIMARY KEY (`sgRNA_pair_id`),
  INDEX (`gene_pair_id`),
  CONSTRAINT FOREIGN KEY (`guide_2_id`) REFERENCES `CDKO_EXPERIMENT_DESIGN` (`sgRNA_id`),
  CONSTRAINT FOREIGN KEY (`guide_1_id`) REFERENCES `CDKO_EXPERIMENT_DESIGN` (`sgRNA_id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_general_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `gemini_score`
--

DROP TABLE IF EXISTS `gemini_score`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `gemini_score` (
  `id` int NOT NULL AUTO_INCREMENT,
  `gene_pair_id` int DEFAULT NULL,
  `SL_score_Strong` double DEFAULT NULL,
  `SL_score_SensitiveLethality` double DEFAULT NULL,
  `SL_score_SensitiveRecovery` double DEFAULT NULL,
  PRIMARY KEY (`id`),
  FOREIGN KEY (`gene_pair_id`) REFERENCES `cdko_sgrna_counts` (`gene_pair_id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_general_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `horlbeck_score`
--

DROP TABLE IF EXISTS `horlbeck_score`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `horlbeck_score` (
  `id` int NOT NULL AUTO_INCREMENT,
  `gene_pair_id` int DEFAULT NULL,
  `SL_score` double DEFAULT NULL,
  `standard_error` double DEFAULT NULL,
  PRIMARY KEY (`id`),
  FOREIGN KEY (`gene_pair_id`) REFERENCES `cdko_sgrna_counts` (`gene_pair_id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_general_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `mageck_score`
--

DROP TABLE IF EXISTS `mageck_score`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `mageck_score` (
  `id` int NOT NULL AUTO_INCREMENT,
  `gene_pair_id` int DEFAULT NULL,
  `SL_score` double DEFAULT NULL,
  `standard_error` double DEFAULT NULL,
  `Z_SL_score` double DEFAULT NULL,
  PRIMARY KEY (`id`),
  FOREIGN KEY (`gene_pair_id`) REFERENCES `cdko_sgrna_counts` (`gene_pair_id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_general_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `median_b_score`
--

DROP TABLE IF EXISTS `median_b_score`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `median_b_score` (
  `id` int NOT NULL AUTO_INCREMENT,
  `gene_pair_id` int DEFAULT NULL,
  `SL_score` double DEFAULT NULL,
  `standard_error` double DEFAULT NULL,
  `Z_SL_score` double DEFAULT NULL,
  PRIMARY KEY (`id`),
  FOREIGN KEY (`gene_pair_id`) REFERENCES `cdko_sgrna_counts` (`gene_pair_id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_general_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `median_nb_score`
--

DROP TABLE IF EXISTS `median_nb_score`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `median_nb_score` (
  `id` int NOT NULL AUTO_INCREMENT,
  `gene_pair_id` int DEFAULT NULL,
  `SL_score` double DEFAULT NULL,
  `standard_error` double DEFAULT NULL,
  `Z_SL_score` double DEFAULT NULL,
  PRIMARY KEY (`id`),
  FOREIGN KEY (`gene_pair_id`) REFERENCES `cdko_sgrna_counts` (`gene_pair_id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_general_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `sgrna_derived_b_score`
--

DROP TABLE IF EXISTS `sgrna_derived_b_score`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `sgrna_derived_b_score` (
  `id` int NOT NULL AUTO_INCREMENT,
  `gene_pair_id` int DEFAULT NULL,
  `SL_score` double DEFAULT NULL,
  PRIMARY KEY (`id`),
  FOREIGN KEY (`gene_pair_id`) REFERENCES `cdko_sgrna_counts` (`gene_pair_id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_general_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `sgrna_derived_nb_score`
--

DROP TABLE IF EXISTS `sgrna_derived_nb_score`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `sgrna_derived_nb_score` (
  `id` int NOT NULL AUTO_INCREMENT,
  `gene_pair_id` int DEFAULT NULL,
  `SL_score` double DEFAULT NULL,
  PRIMARY KEY (`id`),
  FOREIGN KEY (`gene_pair_id`) REFERENCES `cdko_sgrna_counts` (`gene_pair_id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_general_ci;

--
-- Create view
--

DROP VIEW IF EXISTS `joined_counts`;

CREATE VIEW joined_counts         AS            SELECT d.sgRNA_guide_name as sgRNA_guide_name_g1, d.sgRNA_guide_seq as sgRNA_guide_seq_g1, d.sgRNA_target_name as sgRNA_target_name_g1,                       e.sgRNA_guide_name as sgRNA_guide_name_g2, e.sgRNA_guide_seq as sgRNA_guide_seq_g2, e.sgRNA_target_name as sgRNA_target_name_g2,                                         c.* FROM cdko_sgrna_counts c                                        LEFT JOIN cdko_experiment_design d                                        ON c.guide_1_id = d.sgRNA_id                                         LEFT JOIN cdko_experiment_design e                                        ON c.guide_2_id = e.sgRNA_id
/* joined_counts(sgRNA_guide_name_g1,sgRNA_guide_seq_g1,sgRNA_target_name_g1,sgRNA_guide_name_g2,sgRNA_guide_seq_g2,sgRNA_target_name_g2,sgRNA_pair_id,guide_1_id,guide_2_id,gene_pair_id,gene_pair_orientation,T0_counts,T0_replicate_names,TEnd_counts,TEnd_replicate_names,target_type,study_origin,cell_line_origin) */;

DROP VIEW IF EXISTS `calculated_sl_table`;

CREATE VIEW calculated_sl_table          AS             SELECT joined_counts.sgRNA_target_name_g1 gene_1,                                   joined_counts.sgRNA_target_name_g2 gene_2,                                   joined_counts.study_origin study_origin,                                   joined_counts.cell_line_origin cell_line_origin,                                   joined_counts.gene_pair_id gene_pair_id,                                   median_nb_score.SL_score median_nb_score_SL_score,                                   median_nb_score.standard_error median_nb_score_standard_error,                                   median_nb_score.Z_SL_score median_nb_score_Z_SL_score,                                   median_b_score.SL_score median_b_score_SL_score,                                   median_b_score.standard_error median_b_score_standard_error,                                   median_b_score.Z_SL_score median_b_score_Z_SL_score,                                   sgrna_derived_b_score.SL_score sgrna_derived_b_score_SL_score,                                   sgrna_derived_nb_score.SL_score sgrna_derived_nb_score_SL_score,                                   horlbeck_score.SL_score horlbeck_score_SL_score,                                   horlbeck_score.standard_error horlbeck_score_standard_error,                                   mageck_score.SL_score mageck_score_SL_score,                                   mageck_score.standard_error mageck_score_standard_error,                                   mageck_score.Z_SL_score mageck_score_Z_SL_score,                                   gemini_score.SL_score_Strong gemini_score_SL_score_Strong,                                   gemini_score.SL_score_SensitiveLethality gemini_score_SL_score_SensitiveLethality,                                   gemini_score.SL_score_SensitiveRecovery gemini_score_SL_score_SensitiveRecovery                                             FROM gemini_score                                             LEFT JOIN joined_counts                                             ON gemini_score.gene_pair_id = joined_counts.gene_pair_id                                             LEFT JOIN median_b_score                                             ON gemini_score.gene_pair_id = median_b_score.gene_pair_id                                             LEFT JOIN sgrna_derived_b_score                                             ON gemini_score.gene_pair_id = sgrna_derived_b_score.gene_pair_id                                             LEFT JOIN sgrna_derived_nb_score                                             ON gemini_score.gene_pair_id = sgrna_derived_nb_score.gene_pair_id                                             LEFT JOIN horlbeck_score                                             ON gemini_score.gene_pair_id = horlbeck_score.gene_pair_id                                             LEFT JOIN mageck_score                                             ON gemini_score.gene_pair_id = mageck_score.gene_pair_id                                             LEFT JOIN median_nb_score                                             ON gemini_score.gene_pair_id = median_nb_score.gene_pair_id                                                 GROUP BY gemini_score.gene_pair_id
/* calculated_sl_table(gene_1,gene_2,study_origin,cell_line_origin,gene_pair_id,median_nb_score_SL_score,median_nb_score_standard_error,median_nb_score_Z_SL_score,median_b_score_SL_score,median_b_score_standard_error,median_b_score_Z_SL_score,sgrna_derived_b_score_SL_score,sgrna_derived_nb_score_SL_score,horlbeck_score_SL_score,horlbeck_score_standard_error,mageck_score_SL_score,mageck_score_standard_error,mageck_score_Z_SL_score,gemini_score_SL_score_Strong,gemini_score_SL_score_SensitiveLethality,gemini_score_SL_score_SensitiveRecovery) */;

/*!40101 SET character_set_client = @saved_cs_client */;
/*!40103 SET TIME_ZONE=@OLD_TIME_ZONE */;

/*!40101 SET SQL_MODE=@OLD_SQL_MODE */;
/*!40014 SET FOREIGN_KEY_CHECKS=@OLD_FOREIGN_KEY_CHECKS */;
/*!40014 SET UNIQUE_CHECKS=@OLD_UNIQUE_CHECKS */;
/*!40101 SET CHARACTER_SET_CLIENT=@OLD_CHARACTER_SET_CLIENT */;
/*!40101 SET CHARACTER_SET_RESULTS=@OLD_CHARACTER_SET_RESULTS */;
/*!40101 SET COLLATION_CONNECTION=@OLD_COLLATION_CONNECTION */;
/*!40111 SET SQL_NOTES=@OLD_SQL_NOTES */;

-- Dump completed on 2023-04-26 15:29:37
