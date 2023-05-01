# STEMI: MOFA analysis

This repository contains the code to generate the results and figures from

**Multi-Omic Factor Analysis uncovers immunological signatures with pathophysiologic and clinical implications in coronary syndromes**

Kami Pekayvaz*, Corinna Losert*, Viktoria Knottenberg*, Irene V. van Blokland, Roy Oelen, Hilde E. Groot, Jan Walter Benjamins, Sophia Brambs, Rainer Kaiser, Luke Eivers, Vivien Polewka, Raphael Escaig, Markus Joppich, Aleksandar Janjic, Oliver Popp, Tobias Petzold, Ralf Zimmer, Wolfgang Enard, Kathrin Saar, Philipp Mertins, Norbert Huebner, Pim van der Harst, Lude H. Franke, Monique G. P. van der Wijst, Steffen Massberg, Matthias Heinig#, Leo Nicolai# and Konstantin Stark#

\* These authors contributed equally<br>
# These authors contributed equally<br>

Preprint: TBD


## Overview

The code and analysis is structured in different modules indicated by the letter the script name starts with (the sequence in which scripts are executed is indicated by the number behind the letter):

* 00: Loads and prepares all sample-metadata + proteomic, cytokine, neutrophil and clinical data for later usage
* A: Pre-processing of the single-cell RNA-seq data per library
* B: Integration of the different libraries, clustering and annotation
* C: Aggregation of sc data to pseudobulk level
* E: MOFA Analysis and Interpretation
* F: Ligand-Target Analysis
* G: Replication of findings on Groningen cohort

A more detailled documentation can be found in the excel file in the documentation folder.
