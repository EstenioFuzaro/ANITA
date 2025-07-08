## Base 2 – Synthetic Data Using AR Modeling and Monte Carlo Simulation

This dataset simulates the behavior of two heat exchangers with different fluid configurations. The source exchanger (TC1) operates with water–water, while the target (TC2) uses oil–water. Temperature curves are generated synthetically via Monte Carlo simulation by sampling the global heat transfer coefficient $\mathbb{U}$ from a Gamma distribution. Each curve is modeled using AutoRegressive (AR) models to extract key statistical features for fouling detection.

### Execution Order

To reproduce the full analysis for Base 2, follow this order:

1. **`Simula_TC1.m`**  
   Simulates clean and fouled curves for the source exchanger (TC1), saving them to `.mat` files.

2. **`Simula_TC2.m`**  
   Simulates clean and fouled curves for the target exchanger (TC2) under a different fluid regime.

3. **`Identifica_AR_TC1.m`**  
   Identifies an AR model for the source domain, extracts features from all curves, and computes AIC.

4. **`Identifica_AR_TC2.m`**  
   Identifies an AR model for the target domain and extracts features using the reference model.

5. **`Transfer_Learning_TC1_TC2.m`**  
   Performs classification with Support Vector Machines (SVM) and applies Joint Distribution Adaptation (JDA) to align the domains and evaluate transfer learning effectiveness.

All data, models, and features are saved in `.mat` files for transparency and reproducibility.
