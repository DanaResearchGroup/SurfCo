# SurfCo - Protein Corona Prediction Framework

A modular kinetic Monte Carlo (KMC) simulation framework for predicting nanoparticle-protein interactions and corona formation.

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
nanoparticle:
  radius: 35.0              # nm
  material: "carbonblack"   # Must match MaterialSet.csv
  zeta_potential: -0.015    # V

simulation:
  resolution: 0.1           # nm²/cell
  iteration_limit: 200000
  equilibrium_window: 190000
  min_adsorption_probability: 0.70
  max_adsorption_probability: 0.95
  adsorption_event_probability: 0.99
  max_displacement_probability: 0.8
  coverage_threshold_for_displacement: 0.38
  probability_penalty_per_extra_protein: 0.2
  min_displacement_probability: 0.05
  base_displacement_penalty: 0.1

proteins:
  HSA: 45.0                 # mg/mL
  FIB-B: 3.0
  HDL: 15.0
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

## Parameter Descriptions

### Nanoparticle Parameters
- **radius**: Nanoparticle radius in nanometers
- **material**: Material identifier (must exist in MaterialSet.csv)
- **zeta_potential**: Surface potential in volts (used by UnitedAtom)

### Simulation Parameters  
- **resolution**: Grid cell area in nm²/cell (controls spatial discretization)
- **iteration_limit**: Maximum number of Monte Carlo steps
- **equilibrium_window**: Number of iterations with unchanged state to declare equilibrium
- **min/max_adsorption_probability**: Target probability range for weakest/strongest binders
- **adsorption_event_probability**: Fraction of events that are adsorption attempts (vs desorption)

### Vroman Effect Parameters
- **max_displacement_probability**: Maximum displacement chance at full coverage
- **coverage_threshold_for_displacement**: Coverage fraction where displacement begins
- **probability_penalty_per_extra_protein**: Probability reduction per additional displaced protein
- **min_displacement_probability**: Minimum displacement probability floor
- **base_displacement_penalty**: Base probability reduction for any displacement

### Protein Parameters
- **concentrations**: Protein concentrations in mg/mL (determines selection probability)