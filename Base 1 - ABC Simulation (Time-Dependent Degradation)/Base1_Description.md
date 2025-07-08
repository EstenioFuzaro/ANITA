## Base 1 – Synthetic Data from Reduced-Order Model with ABC Calibration

This dataset simulates two concentric shell-and-tube heat exchangers with different geometries, both using water as working fluid. The thermal degradation over time is modeled through a reduced-order set of ODEs, where the global heat transfer coefficient $\mathbb{U}$ decays exponentially to emulate fouling. The parameters are calibrated using Approximate Bayesian Computation (ABC), generating temperature curves in both clean and fouled conditions.

### Execution Order

To reproduce the full pipeline for Base 1, execute the following scripts in order:

1. **`trocador_1.m`**  
   Generates the simulated temperature curves for the first (source) heat exchanger using a reduced-order model and stochastic $\mathbb{U}(t)$ decay. Also includes ABC-based calibration.

2. **`trocador_2.m`**  
   Repeats the simulation for a second (target) heat exchanger with different geometric and operational parameters.

3. **`transfer_learning_trocadores.m`**  
   Applies Support Vector Machine (SVM) classification and Joint Distribution Adaptation (JDA) to evaluate fouling detection across domains using Transfer Learning.

The resulting temperature datasets and diagnostic performance metrics are saved automatically at the end of each script.
