# A Comparison of Causal Forests and the DR-Learner for Estimating Conditional Average Treatment Effects 
This repository shows the code, truth data and output for the paper "A Comparison of Causal Forests and the DR-Learner for Estimating Conditional Average Treatment Effects"

## Abstract
Conditional average treatment effects (CATEs) hold great promise for precision medicine, particularly in settings where effect modification is likely. Theoretical work has developed methods to estimate CATEs, including the double-robust (DR) learner and the causal forest algorithm. Here, we conduct a simulation study comparing the finite sample properties of the DR learner and the causal forest algorithm. We explore performance in a range of scenarios with a binary effect modifier and when a set of conditioning variables are included with varying degrees of effect modifiers present. Scenarios explored different effect parametrizations, sample sizes, proportions of modifying to non-modifying variables, and number of confounding variables. For all analyses, we used 10-fold cross fitting, and linear projection approach to identify pre-specified modifiers. Our results suggest that both the causal forest and the DR learner have good 95% confidence interval coverage in most settings. However, the DR learner outperformed the causal forest in coverage (93% vs. 88%) under the scenarios of strong treatment effect but low heterogeneity. We also found that the best linear projections may not always reliably identify pre-specified effect modifiers when either method is used, especially in small sample sizes (with a successful identification of 53% at best scenario).  This will provide practical insights to guide method selection for estimating CATEs in empirical research. 

## Data simulation parameters
To get the parameters for the data generating process that targeted different risk differences across scenarios, run `code
/dgm_true_RD.R` to produce the truth file "truth_data_simulation_all_update.csv". 

## Simulation and analysis
Run `code/run_cluster_sim.R` to intiate the data simulation, estimating CATE using DR Learner and causal forest, evaluting finite sample property and identifying the effect modifier with best linear projection.

## Evaluation and figures
* Best linear projection: `code/blp.R`
* Finite sample property: `code/figure.R`
* Comparison of performance of causal forest with vs. without the Super Learner: `code/figure_cf.R`


