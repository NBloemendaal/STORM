# STORM
This is the STORM (Synthetic Tropical cyclOne geneRation Model) Model. For details, see https://doi.org/10.1038/s41597-020-0381-2

**IMPORTANT: Please be aware that these scripts are not maintained and NO support is provided!!**

These python programs are all modules:
1. SELECT_BASIN.py <--- generates the basin boundaries, no of cyclones, month of occurrence
2. SAMPLE_STARTING_POINT.py <--- generates the genesis locations
3. SAMPLE_TC_MOVEMENT.py <--- creates the tropical cyclone track
4. SAMPLE_TC_PRESSURE.py <--- creates the pressure and wind speed changes along track
5. SAMPLE_RMAX.py <--- samples the radius to maximum wind along the track

MASTER.py is the master-program, which is the one you should run to generate the synthetic tracks.
Input data for each of the modules is generated using the pre-processing scripts, see respective repository for those scripts.
