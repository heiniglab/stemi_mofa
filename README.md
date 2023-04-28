# STEMI: MOFA analysis

This repository contains the code to generate the results and figures from

**Multi-Omic Factor Analysis uncovers immunological signatures with pathophysiologic and clinical implications in coronary syndromes **

Kami Pekayvaz*,#, Corinna Losert*, Viktoria Knottenberg*, Sophia Brambs, Rainer Kaiser, Luke Eivers, Vivien Polewka, Raphael Escaig, Markus Joppich, Aleksandar Janjic, Oliver Popp, Tobias Petzold, Ralf Zimmer, Wolfgang Enard, Kathrin Saar, Philipp Mertins, Norbert Huebner, Steffen Massberg, Matthias Heinig#,§, Leo Nicolai#,§ and Konstantin Stark#,§ 

\* These authors contributed equally<br>
# These authors contributed equally<br>

Preprint: TBD


## Overview

The code and analysis is structured in different modules indicated by the letter the script name starts with (the sequence in which scripts are executed is indicated by the number behind the letter):

**00: Loads and prepares all sample-metadata + proteomic, cytokine, neutrophil and clinical data for later usage
**A: Pre-processing of the single-cell RNA-seq data per library
**B: Integration of the different libraries, clustering and annotation
**C: Aggregation of sc data to pseudobulk level
**E: MOFA Analysis and Interpretation
**F: Ligand-Target Analysis
**G: Replication of findings on Groningen cohort