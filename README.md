# SurfCo - Protein Corona Prediction Framework

A modular kinetic Monte Carlo (KMC) simulation framework for predicting nanoparticle-protein interactions and corona formation.

<img width="2950" height="2616" alt="final_surface" src="https://github.com/user-attachments/assets/3e4402cb-ac94-437e-85e9-3a874f2ff036" />

## Installation

### Prerequisites
- Python 3.7+
- C++ compiler (g++)
- Boost C++ libraries
- Make build system

### Quick Install
```bash
git clone https://github.com/DanaResearchGroup/SurfCo.git
cd SurfCo
chmod +x install.sh
./install.sh
```

The installer will:
1. Set up Python virtual environment (optional)
2. Install Python dependencies
3. Install system dependencies
4. Compile UnitedAtom component
5. Verify installation

## Running SurfCo

### Basic Usage
```bash
# Activate virtual environment (if created)
source surfco_env/bin/activate

# Run simulation
python run_pipeline.py <project_name>
```

### Command Options
```bash
python run_pipeline.py --help                    # Show help
python run_pipeline.py --skip-energy <project>   # Skip energy calculations  
python run_pipeline.py --skip-projection <project> # Skip projections
python run_pipeline.py --skip-simulation <project> # Skip simulation
python run_pipeline.py --show-c                    # Show citations
python run_pipeline.py --show-w                    # Show warranty
```

## Required Inputs

### Directory Structure
```
input/
└── <project_name>/
    ├── config.yaml          # Simulation parameters
    ├── protein1.pdb         # Protein structure files
    ├── protein2.pdb
    └── ...
```

### Configuration File (config.yaml)
```yaml
# SurfCo Configuration File - Probability-Based Penalties

# Nanoparticle properties
nanoparticle:
  radius: 35.0                    # nm - nanoparticle radius
  material: "carbonblack"         # must match MaterialSet.csv
  zeta_potential: -0.015          # V - surface potential 

# Simulation parameters
simulation:
  resolution: 0.1                 # nm² per grid cell
  iteration_limit: 200000          # maximum iterations
  equilibrium_window: 190000        # iterations for equilibrium check
  
  # Energy normalization - maps raw energies to probability range
  min_adsorption_probability: 0.70  # weakest binder
  max_adsorption_probability: 0.95  # strongest binder
  
  # Event selection ratio - what fraction of events are adsorption attempts
  adsorption_event_probability: 0.99  
  
  # Vroman effect control - coverage-dependent displacement
  max_displacement_probability: 0.8   # Maximum displacement chance at high coverage (0-1)
  coverage_threshold_for_displacement: 0.38  # Coverage % where displacement starts (0-1)
  
  # Probability-based penalty system 
  probability_penalty_per_extra_protein: 0.2  # % probability reduction per extra protein displaced
  min_displacement_probability: 0.05          # Never go below % displacement chance
  base_displacement_penalty: 0.1             # % reduction for any displacement attempt

# Protein concentrations (mg/mL)
proteins:
  HSA: 45.0      # Human Serum Albumin
  FIB-B: 3.0     # Fibrinogen B chain
  HDL: 15.0      # High-density lipoprotein
  FIB-A: 3.0     # Fibrinogen A chain
```

### Required Data Files
- `data/MaterialSet.csv` - Material properties database
- `data/PMF/` - Potential of mean force files  
- `data/Hamaker/` - Hamaker constant files
- Protein PDB structures in `input/<project_name>/`

## Generated Outputs

### Directory Structure
```
projects/<project_name>/
├── logs/
│   ├── simulation_log.csv      # Complete iteration log
│   ├── summary.json            # Final results summary
│   └── timing_summary.json     # Performance metrics
├── visualizations/
│   ├── coverage/
│   │   ├── coverage_complete_smooth.png
│   │   ├── final_composition.png
│   │   └── protein_accumulation_smooth.png
│   ├── masks/
│   │   └── *_mask.png          # Protein projection masks
│   └── surface/
│       └── final_surface.png   # Final surface state
└── working_files/
    ├── energies/               # UnitedAtom outputs
    ├── projections/            # Protein masks and spheres
    └── simulation/             # Final surface matrix
```

### Key Output Files
- **simulation_log.csv**: Complete iteration history with protein counts, coverage, and events
- **summary.json**: Final composition, equilibrium status, and statistics  
- **final_surface.png**: Visual representation of final protein corona
- **coverage plots**: Time evolution of protein coverage dynamics

## Mathematical Framework

### Energy Normalization
Raw binding energies are mapped to probability ranges:
```
E_target_min = ln((1/P_min) - 1)    # For weakest binder
E_target_max = ln((1/P_max) - 1)    # For strongest binder

E_norm = E_target_min + (E_target_max - E_target_min) * fraction
```
where `fraction = (E_raw - E_raw_weakest)/(E_raw_strongest - E_raw_weakest)`

### Adsorption Probability
Probability of successful binding upon surface encounter:
```
P_ads(E_norm) = 1/(1 + exp(E_norm))
```

### Desorption Probability  
Probability of spontaneous protein release:
```
P_des(E_norm) = 1/(1 + exp(-E_norm))
```

### Coverage-Dependent Displacement (Vroman Effect)
```
P_displacement(θ) = P_max * (θ - θ_threshold)/(1 - θ_threshold)    if θ > θ_threshold
                  = 0                                              if θ ≤ θ_threshold
```
where θ is surface coverage fraction.

### Probability-Based Displacement Penalties
For multi-protein displacement events:
```
P_replace = P_base * (1 - penalty_base) * (1 - n_extra * penalty_per_protein)
P_final = max(P_replace, P_min)
```

### Protein Selection
Concentration-weighted random selection:
```
P_select(protein_i) = C_i / Σ(C_j)
```
where C_i is the concentration of protein i.
