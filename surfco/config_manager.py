"""
SurfCo - Protein Corona Prediction Framework
Copyright (C) 2024 Dana Research Group

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""

import yaml
import numpy as np
import pandas as pd
from pathlib import Path


class ConfigManager:
    def __init__(self, project_name):
        self.project_name = project_name
        self.base_dir = Path.cwd()
        
        self._setup_directories()
        self._load_config()
        self._detect_proteins()
        self._load_material_data()
        self._calculate_surface_properties()
    
    def _setup_directories(self):
        self.input_dir = self.base_dir / "input" / self.project_name
        if not self.input_dir.exists():
            raise ValueError(f"Project folder not found: {self.input_dir}")
        
        self.project_dir = self.base_dir / "projects" / self.project_name
        self.working_dir = self.project_dir / "working_files"
        self.viz_dir = self.project_dir / "visualizations"
        self.log_dir = self.project_dir / "logs"
        
        self.energy_dir = self.working_dir / "energies"
        self.proj_dir = self.working_dir / "projections"
        self.sim_dir = self.working_dir / "simulation"
        
        self.mask_viz_dir = self.viz_dir / "masks"
        self.coverage_viz_dir = self.viz_dir / "coverage"
        self.surface_viz_dir = self.viz_dir / "surface"
        
        for directory in [self.project_dir, self.working_dir, self.viz_dir, self.log_dir,
                         self.energy_dir, self.proj_dir, self.sim_dir,
                         self.mask_viz_dir, self.coverage_viz_dir, self.surface_viz_dir]:
            directory.mkdir(parents=True, exist_ok=True)
    
    def _load_config(self):
        config_file = self.input_dir / "config.yaml"
        if not config_file.exists():
            raise ValueError(f"Config file not found: {config_file}")
        
        with open(config_file, 'r') as f:
            config = yaml.safe_load(f)
        
        try:
            np_config = config['nanoparticle']
            self.np_radius = np_config['radius']  # nm
            self.material = np_config['material']
            self.zeta_potential_volts = np_config['zeta_potential']  # V (for UA)
            self.zeta_potential_mv = self.zeta_potential_volts * 1000  # mV (for display)
        except KeyError as e:
            raise ValueError(f"Missing required nanoparticle parameter: {e}")
        
        try:
            sim_config = config['simulation']
            self.resolution = sim_config['resolution']  # nm² per cell
            self.iteration_limit = sim_config['iteration_limit']
            self.equilibrium_window = sim_config['equilibrium_window']
            
            self.min_adsorption_probability = sim_config['min_adsorption_probability']
            self.max_adsorption_probability = sim_config['max_adsorption_probability']
            
            self.adsorption_event_probability = sim_config['adsorption_event_probability']
            
            self.max_displacement_probability = sim_config['max_displacement_probability']
            self.coverage_threshold_for_displacement = sim_config['coverage_threshold_for_displacement']
            
            self.probability_penalty_per_extra_protein = sim_config['probability_penalty_per_extra_protein']
            self.min_displacement_probability = sim_config['min_displacement_probability']
            self.base_displacement_penalty = sim_config['base_displacement_penalty']
            
        except KeyError as e:
            raise ValueError(f"Missing required simulation parameter: {e}. "
                           f"Please ensure your config.yaml file contains all required parameters.")
        
        try:
            self.protein_concentrations = config['proteins']  # mg/mL
        except KeyError:
            raise ValueError("Missing required 'proteins' section in config file")
        
        self._validate_parameters()
    
    def _validate_parameters(self):
        """Validate all configuration parameters with clear error messages."""
        
        if self.resolution <= 0:
            raise ValueError("resolution must be positive")
        if self.iteration_limit <= 0:
            raise ValueError("iteration_limit must be positive")
        if self.equilibrium_window <= 0:
            raise ValueError("equilibrium_window must be positive")
        if self.equilibrium_window >= self.iteration_limit:
            raise ValueError("equilibrium_window must be less than iteration_limit")
        
        if not (0.01 <= self.min_adsorption_probability <= 0.99):
            raise ValueError("min_adsorption_probability must be between 0.01 and 0.99")
        if not (0.01 <= self.max_adsorption_probability <= 0.99):
            raise ValueError("max_adsorption_probability must be between 0.01 and 0.99")
        if self.min_adsorption_probability >= self.max_adsorption_probability:
            raise ValueError("min_adsorption_probability must be less than max_adsorption_probability")
        
        if not (0.01 <= self.adsorption_event_probability <= 0.99):
            raise ValueError("adsorption_event_probability must be between 0.01 and 0.99")
        
        if not (0.0 <= self.max_displacement_probability <= 1.0):
            raise ValueError("max_displacement_probability must be between 0.0 and 1.0")
        if not (0.0 <= self.coverage_threshold_for_displacement <= 1.0):
            raise ValueError("coverage_threshold_for_displacement must be between 0.0 and 1.0")
        
        if not (0.0 <= self.probability_penalty_per_extra_protein <= 1.0):
            raise ValueError("probability_penalty_per_extra_protein must be between 0.0 and 1.0")
        if not (0.0 <= self.min_displacement_probability <= 1.0):
            raise ValueError("min_displacement_probability must be between 0.0 and 1.0")
        if not (0.0 <= self.base_displacement_penalty <= 1.0):
            raise ValueError("base_displacement_penalty must be between 0.0 and 1.0")
        
        if self.np_radius <= 0:
            raise ValueError("nanoparticle radius must be positive")
        if not isinstance(self.material, str) or len(self.material) == 0:
            raise ValueError("material must be a non-empty string")
        
        if not isinstance(self.protein_concentrations, dict):
            raise ValueError("proteins section must be a dictionary")
        if len(self.protein_concentrations) == 0:
            raise ValueError("at least one protein must be specified")
        
        for protein_name, concentration in self.protein_concentrations.items():
            if not isinstance(concentration, (int, float)):
                raise ValueError(f"concentration for {protein_name} must be a number")
            if concentration <= 0:
                raise ValueError(f"concentration for {protein_name} must be positive")
    
    def _detect_proteins(self):
        pdb_files = list(self.input_dir.glob("*.pdb"))
        if not pdb_files:
            raise ValueError(f"No PDB files found in {self.input_dir}")
        
        self.proteins = {}
        for pdb_file in pdb_files:
            protein_name = pdb_file.stem
            
            if protein_name not in self.protein_concentrations:
                raise ValueError(f"Protein {protein_name} found in PDB files but not specified in config proteins section")
            
            self.proteins[protein_name] = {
                'pdb_path': pdb_file,
                'concentration': self.protein_concentrations[protein_name],
                'raw_energy': None,
                'normalized_energy': None,
                'mask': None,
                'area': None,
                'spheres': None
            }
        
        for protein_name in self.protein_concentrations:
            if protein_name not in self.proteins:
                raise ValueError(f"Protein {protein_name} specified in config but no corresponding PDB file found")
    
    def _load_material_data(self):
        materials_csv = self.base_dir / "data" / "MaterialSet.csv"
        if not materials_csv.exists():
            raise ValueError(f"MaterialSet.csv not found at {materials_csv}")
        
        df = pd.read_csv(materials_csv, comment='#')
        df.columns = df.columns.str.strip()
        
        material_row = df[df.iloc[:, 0] == self.material]
        if material_row.empty:
            available_materials = df.iloc[:, 0].tolist()
            raise ValueError(f"Material '{self.material}' not found in MaterialSet.csv. "
                           f"Available materials: {available_materials}")
        
        row = material_row.iloc[0]
        
        pmf_path = str(row.iloc[1]).strip()
        hamaker_path = str(row.iloc[2]).strip()
        
        if pmf_path.startswith("data/"):
            self.pmf_folder = self.base_dir / pmf_path
        else:
            self.pmf_folder = self.base_dir / "data" / pmf_path
            
        if hamaker_path.startswith("data/"):
            self.hamaker_file = self.base_dir / hamaker_path
        else:
            self.hamaker_file = self.base_dir / "data" / hamaker_path
        
        self.shape = int(row.iloc[3])
        self.pmf_cutoff = float(row.iloc[4])
        
        print(f"  PMF folder: {self.pmf_folder}")
        print(f"  Hamaker file: {self.hamaker_file}")
        
        if not self.pmf_folder.exists():
            raise ValueError(f"PMF folder not found: {self.pmf_folder}")
        if not self.hamaker_file.exists():
            raise ValueError(f"Hamaker file not found: {self.hamaker_file}")
    
    def _calculate_surface_properties(self):
        self.surface_area = 4 * np.pi * (self.np_radius ** 2)  # nm²
        self.n_cells = int(self.surface_area / self.resolution)
        self.grid_size = int(np.ceil(np.sqrt(self.n_cells)))
        self.cell_size = np.sqrt(self.resolution)  # nm
    
    @property 
    def zeta_potential(self):
        """Return zeta potential in volts (for UA compatibility)"""
        return self.zeta_potential_volts