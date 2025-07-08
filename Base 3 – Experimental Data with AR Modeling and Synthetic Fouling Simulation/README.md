## Base 3 – Experimental Data from Real Heat Exchangers

This dataset is based on experimental measurements collected from two real heat exchangers with distinct geometries: a shell-and-tube (CascoTubo) and a brazed-plate (PlacaPlana) unit. Both devices operate with water and were instrumented under clean conditions. The fouling condition was simulated by altering the thermal response parameters identified from experimental curves. Synthetic curves were generated using Monte Carlo simulation with Gamma-distributed parameters.

### Execution Order

To reproduce the analysis for Base 3, follow this order:

1. **`MonteCarloCurves_CascoTubo.m`**  
   Generates synthetic clean and fouled curves for the shell-and-tube exchanger based on fitted dynamic parameters.

2. **`MonteCarloCurves_PlacaPlana.m`**  
   Generates synthetic clean and fouled curves for the plate exchanger using the same Monte Carlo methodology.

3. **`IdentificaAR_CascoTubo.m`**  
   Identifies AR models from the CascoTubo data and extracts statistical features.

4. **`IdentificaAR_PlacaPlana.m`**  
   Identifies AR models from the PlacaPlana data and extracts features using the reference model.

5. **`Transfer_Learning_Trocadores.m`**  
   Applies SVM classification and evaluates Transfer Learning performance via Joint Distribution Adaptation (JDA).

All experimental data, model parameters, and extracted features are stored in `.mat` files to ensure reproducibility.
