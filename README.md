# global_epistasis_review

Code to generate the panels in the figures of the paper:

Juan Diaz-Colunga, Abigail Skwara, Karna Gowda, Ramon Diaz-Uriarte, Mikhail Tikhonov, Djordje Bajic & Alvaro Sanchez (2022). Global epistasis patterns on fitness landscapes.

The ./scripts directory contains R scripts and data files necessary to produce the plots. Plots are saved under the ./plots directory. Files description:

**eFunctions.R**   Contains auxiliary functions for data preprocessing (described in [[1](https://doi.org/10.1101/2022.06.21.496987)]).

**landscape_regions.R**   Produces the panels for Figs. 1 and 2, plus some additional plots for slightly more complex landscape topologies (not discussed in the paper).

**khan_fij.R**   Generates the panels of Fig. 3. Imports data from the genotype_fitness.txt file (from [[2](https://doi.org/10.1126/science.1203801)]).

**landscape_nonlineartrans.R**   Produces the panels for Fig. 4.



[1] Juan Diaz-Colunga, Abigail Skwara, Jean C. C. Vila, Djordje Bajic, Alvaro Sanchez (2022). Emergent ecosystem functions follow simple quantitative rules. *biorXiv*. DOI: [10.1101/2022.06.21.496987](https://doi.org/10.1101/2022.06.21.496987)

[2] Aisha I. Khan, Duy M. Dinh, Dominique Schneider, Richard E. Lenski, Tim F. Cooper (2011). Negative epistasis between beneficial mutations in an evolving bacterial population. *Science* **332(6034)**:1193–1196. DOI: [10.1126/science.1203801](https://doi.org/10.1126/science.1203801)