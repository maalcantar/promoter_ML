# Promoter discovery using natural language processing 

This repo contains code for promoter discovery and biodiversity mining using natural language processing. This effort is being pursued through three specific aims:

1. Aim 1: Develop natural processing-based model for promoter identification
2. Aim 2: Extend model to identify inducible promoter sequences
3. Aim 3: Experimentally validate promoter predictions 

Promoter sequences were collected from three main databases: EPDnew, RegulonDB, DBTBS. 

# Directory structure

### data

All data files are found in and/or will be written to <code>data/</code>

* <code>data/DBTBS/</code>
  * Contains raw data from DBTBS: *Bacillus subtilis* promoter database 
* <code>data/EPDnew/</code>
  * Contains raw data from EPDnew: *Eukaryote* promoter database (promoter data for 15 different organisms
* <code>data/RegulonDB/</code>
  * Contains raw data from RegulonDB: *Escherichia coli* promoter database
* <code>data/parsed_promoter_data/</code>
  * Promoter data parsed from each database
* <code>data/20191114promoter_identification_ML_curation</code>
  * Manually curated information on other state-of-the-art ML models for promoter prediction

### src

All code are found in and/or will be written to <code>src/</code> in either notebook or script form

**Notebooks**:
* <code>src/notebooks/20191125_promoter_database_parsing.ipynb</code>
  * Notebook containing code for parsing promoter data
  
### figs

All raw and edits figures will be writted to <code>figs/</code>


