Simulation study and empirical data application for JSPCA

## Authors
- Tra T. Le | Tilburg University
- Dr. Katrijn Van Deun (Supervisor)

## Description
These are the R codes used for the simulation study and empirical data application of Joint Sparse Principal Component Analysis. The analyses were done as one of the four traineeships of the Research Master in Social andn Behavioral Science at Tilburg University. The main method (joint sparse PCA) was proposed and developed by Dr. Katrijn Van Deun, who supervised this traineeship. 

The repository consists of the codes for the **Simulation Study** and **Empirical Data Analysis** using the Genetic Dataset (provided by Dr. Van Deun). 

## Simulation Study
Conditions
- Sample sizes: 50, 200 600, 1000
- Type of loading differences: Primary Loading shift, Cross Loading of 0.4 and 0.2, Primary Loading decrease of 0.2 and 0.4
- Size of loading differences: 4, 16
- Number of groups: 2, 4
- Number of components: 2, 4
### Joint Sparse PCA 
Performance measures: Tucker's congruence and zero/non-zero recovery rate. The complete simulation study can be done using R. Here, [2groups_2components.R](2groups_2components.R) is an example for the combination of conditions: 2 groups, 2 components, Primary loading shift, 4 differences with all sample sizes. 

### Multigroup factor rotation
These are preliminary codes that were translated from the authors' original Matlab codes, and should be used with great caution! Here, [2G2R_PLshift_4diff.R](2G2R_PLshift.R) is an example for the combination of condition: 2 groups, 2 components, Primary loading shift, 4 differences with sample size N = 50 (the codes can easily be tweaked for all sample sizes). For each combination of condition mentioned above, there are **10 combinations of rotation criteria** (2 types of rotation with 5 weights). 

The analysis is done in LatentGOLD6.0 following the steps:
1. Step 1: Run [NoRotation_N50sim1_G2R2PLshift_4diff.lgs](NoRotation_N50sim1_G2R2PLshift_4diff.lgs) to get the parameters, which will serve as starting values for the next step. The starting values are saved in [startingvalues_N50sim1G2R2PLshift_4diff.txt](startingvalues_N50sim1G2R2PLshift_4diff.txt).
2. Step 2: Use [data_N50sim1G2R2PLshift_4diff.sav](data_N50sim1G2R2PLshift_4diff.sav) as infile and the starting values above, run [Rotated_OA_N50OA0.01sim1_G2R2PLshift_4diff.lgs](Rotated_OA_N50OA0.01sim1_G2R2PLshift_4diff.lgs) for Oblimin Alignment and [Rotated_OP_N50OP0.01sim1_G2R2PLshift_4diff.lgs](Rotated_OP_N50OP0.01sim1_G2R2PLshift_4diff.lgs) for Oblimin Proscuste to get the final loadings.

Performance measures: Only Tucker's congruence! In the latest version of LatentGOLD, the rotation criteria can be specified using the command 
```
<simple structure criterion> <agreement criterion><w> sets = 1;
```
Here in this example, we used 
```
oblimin alignment = 0.01 sets=1;
```
However, LatentGOLD can be called from R, so running the R script [2G2R_PLshift_4diff.R](2G2R_PLshift.R) alone is sufficient to generate data, execute LatentGOLD, and read the results. It is not advised to parallel process LatentGOLD since the license can get overwhelmed.
