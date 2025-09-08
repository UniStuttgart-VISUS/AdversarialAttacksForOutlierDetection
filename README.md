# Supplemental Material for the paper ["Robustness and accuracy of mean opinion scores with hard and soft outlier detection"](https://ieeexplore.ieee.org/Xplore/home.jsp)
Dietmar Saupe, Tim Bleile

---

This supplemental material provides the source code for all methods evaluated in our study, enabling reproduction of the main experiments described in the paper. 


### Reference
The reference that should be cited for their usage is:
- **Saupe, D., Bleile, T., Robustness and accuracy of mean opinion scores with hard and soft outlier detection, 17th International Conference on Quality of Multimedia Experience (QoMEX), Sept./Oct. 2025, Madrid, Spain.**

Included are:

### Code

- **main.m** — The main script to run accuracy tests for selected outlier detection methods. Users can choose which methods to evaluate, set the number of iterations, and configure the dataset dimensions, including the number of items, subjects, and attacker subjects. Upon completion, the script generates a LaTeX-formatted table summarizing the evaluation metrics, similar to Table 2 in the paper.
- **geneticAlgorithm.m** — The genetic algorithm implementation used for optimization of the attacks.
- **perform_CB.m** through **perform_ZREC.m** — MATLAB implementations of all outlier detection methods evaluated in this work. Among these, **perform_SUREAL.m**, **perform_ESQR.m**, and **perform_ZREC.m** are adapted from the implementation provided by Altieri *et al.* (IEEE Trans. Multimedia, 2024), available [here](https://doi.org/10.1109/TMM.2024.3390113).
- **simulation.m** — Simulates a subjective rating matrix with *I* subjects and *J* items based on the SUREAL model, incorporating subject bias and inconsistency. This code uses the `brcw.mat` file containing bias and inconsistency values derived from the KonIQ-10k dataset.


### Helper Functions

- **calculateMaximalDeviation.m** — Computes the maximal deviation metric used in attack evaluation.
- **generateAttackSet.m** — Generates set of attackers, with random ratings.
- **entropyCalculation.m** — Performs entropy calculations needed in the **perform_HB.m** method.
- **fisher_z_transformation.m** — Supporting function used by the **perform_ESQR.m** implementation.


