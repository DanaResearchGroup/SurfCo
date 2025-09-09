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

import subprocess
import numpy as np
from pathlib import Path
import os
import time
import threading
import sys


class Spinner:
    """Simple spinner implementation"""
    def __init__(self, message="Working"):
        self.message = message
        self.spinning = False
        self.spinner_chars = "|/-\\"
        self.idx = 0
        
    def start(self):
        self.spinning = True
        self.thread = threading.Thread(target=self._spin)
        self.thread.daemon = True
        self.thread.start()
        
    def _spin(self):
        while self.spinning:
            sys.stdout.write(f'\r      {self.message}... {self.spinner_chars[self.idx]}')
            sys.stdout.flush()
            self.idx = (self.idx + 1) % len(self.spinner_chars)
            time.sleep(0.2)
            
    def stop(self):
        self.spinning = False
        if hasattr(self, 'thread'):
            self.thread.join(timeout=0.5)
        sys.stdout.write(f'\r      {self.message}... Done\n')
        sys.stdout.flush()


class EnergyCalculator:
    def __init__(self, config):
        self.config = config
        self.ua_binary = self.find_ua_binary()
    
    def find_ua_binary(self):
        """Find the UnitedAtom binary in various possible locations"""
        possible_paths = [
            Path("./UnitedAtom"),                    
            Path("./ua/build/UnitedAtom"),          
            Path("./ua/UnitedAtom"),                
        ]
        
        for path in possible_paths:
            if path.exists() and path.is_file():
                print(f"  Found UA binary at: {path}")
                return path.absolute()
        
        raise FileNotFoundError(
            "UnitedAtom binary not found. Please ensure UA is properly installed.\n"
            "Try running: bash install_ua.sh\n"
            f"Searched locations: {[str(p) for p in possible_paths]}"
        )
    
    def generate_ua_config(self, protein_name):
        config_path = self.config.energy_dir / f"{protein_name}_ua.config"
        pdb_path = self.config.proteins[protein_name]['pdb_path']
        
        T = 310.15  # K
        T_C = T - 273.15
        I = 0.15  # M
        
        eps_rel = 87.740 - 0.4008 * T_C + 9.398e-4 * T_C**2 - 1.410e-6 * T_C**3
        eps0 = 8.854e-12  # F/m
        n_avo = 6.022e23
        kb = 1.38e-23  # J/K
        e = 1.6e-19  # C
        
        debye_length = 1e9 * np.sqrt(eps_rel * eps0 * kb * T / (2 * e**2 * I * n_avo * 1000))
        bjerrum_length = 1e9 * e**2 / (4 * np.pi * eps0 * eps_rel * kb * T)
        
        output_dir = "."
        
        with open(config_path, 'w') as f:
            f.write(f"output-directory = {output_dir}\n")
            f.write(f"pdb-target = {pdb_path}\n")
            f.write(f"nanoparticle-radius = [{self.config.np_radius}]\n")
            f.write(f"np-type = {self.config.shape}\n")
            f.write(f"zeta-potential = [{self.config.zeta_potential}]\n")
            f.write(f"pmf-directory = {self.config.pmf_folder}\n")
            f.write(f"hamaker-file = {self.config.hamaker_file}\n")
            f.write("enable-surface\n")
            f.write("enable-core\n")
            f.write("enable-electrostatic\n")
            f.write("simulation-steps = 2000\n")
            f.write("angle-delta = 5.0\n")
            f.write("potential-cutoff = 5.0\n")
            f.write("potential-size = 1000\n")
            f.write(f"bjerum-length = {bjerrum_length:.3f}\n")
            f.write(f"debye-length = {debye_length:.3f}\n")
            f.write(f"temperature = {T}\n")
            f.write(f"pmf-cutoff = {self.config.pmf_cutoff}\n")
            f.write("amino-acids = [ALA, ARG, ASN, ASP, CYS, GLN, GLU, GLY, HIS, ILE, ")
            f.write("LEU, LYS, MET, PHE, PRO, SER, THR, TRP, TYR, VAL]\n")
            f.write("amino-acid-charges = [0.0, 1.0, 0.0, -1.0, 0.0, 0.0, -1.0, 0.0, 0.5, ")
            f.write("0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]\n")
            f.write("amino-acid-radii = [0.323, 0.429, 0.362, 0.356, 0.352, 0.386, 0.376, ")
            f.write("0.285, 0.302, 0.401, 0.401, 0.405, 0.402, 0.421, 0.362, 0.328, ")
            f.write("0.357, 0.449, 0.425, 0.377]\n")
        
        return config_path
    
    def run_ua(self, protein_name):
        config_file = self.generate_ua_config(protein_name)
        
        original_cwd = Path.cwd()
        energy_dir = self.config.energy_dir
        
        try:
            os.chdir(energy_dir)
            
            cmd = [str(self.ua_binary), f"--config-file={config_file.name}"]
            
            spinner = Spinner(f"Scanning orientations for {protein_name} in UA")
            spinner.start()
            
            try:
                result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
            finally:
                spinner.stop()
            
            if result.returncode != 0:
                print(f"      UA failed with return code: {result.returncode}")
                if result.stderr:
                    print(f"      Error: {result.stderr.strip()}")
                raise RuntimeError(
                    f"UnitedAtom failed for protein {protein_name}\n"
                    f"Return code: {result.returncode}\n"
                    f"STDERR: {result.stderr}"
                )
                
        except subprocess.TimeoutExpired:
            spinner.stop()
            raise RuntimeError(f"UnitedAtom timed out for protein {protein_name}")
        except Exception as e:
            if 'spinner' in locals():
                spinner.stop()
            raise RuntimeError(f"Failed to run UnitedAtom for protein {protein_name}: {str(e)}")
        finally:
            os.chdir(original_cwd)
        
        all_uam_files = []
        for pattern in ["*.uam", "**/*.uam"]:
            all_uam_files.extend(list(energy_dir.glob(pattern)))
        
        uam_files = []
        for uam_file in all_uam_files:
            if protein_name.lower() in uam_file.name.lower():
                uam_files.append(uam_file)
        
        if not uam_files and all_uam_files:
            uam_files = [all_uam_files[0]]
        
        if not uam_files:
            raise RuntimeError(
                f"UnitedAtom completed but no .uam file found for protein {protein_name}\n"
                f"Expected .uam files in: {energy_dir}"
            )
        
        return self.parse_uam_file(uam_files[0])
    
    def parse_uam_file(self, uam_file):
        try:
            data = np.loadtxt(uam_file)
            if data.ndim == 1:
                data = data.reshape(1, -1)
            
            if data.shape[1] < 3:
                raise ValueError(f"Invalid .uam file format: {uam_file}")
            
            min_idx = np.argmin(data[:, 2])
            
            return {
                'energy': data[min_idx, 2],
                'phi': data[min_idx, 0],
                'theta': data[min_idx, 1]
            }
        except Exception as e:
            raise RuntimeError(f"Failed to parse .uam file {uam_file}: {str(e)}")
    
    def calculate_all_energies(self):
        print("  Calculating binding energies...")
        
        for protein_name in self.config.proteins:
            print(f"    Processing {protein_name}...")
            result = self.run_ua(protein_name)
            self.config.proteins[protein_name]['raw_energy'] = result['energy']
            self.config.proteins[protein_name]['optimal_phi'] = result['phi']
            self.config.proteins[protein_name]['optimal_theta'] = result['theta']
            
            print(f"    {protein_name}: E = {result['energy']:.2f} kT")
        
        self._normalize_energies()
    
    def _normalize_energies(self):
        """Probability range normalization - FIXED VERSION"""
        p_min = self.config.min_adsorption_probability
        p_max = self.config.max_adsorption_probability
        
        E_target_for_min_prob = np.log((1/p_min) - 1)  
        E_target_for_max_prob = np.log((1/p_max) - 1)  
        
        energies = []
        for protein in self.config.proteins.values():
            if protein['raw_energy'] is not None:
                energies.append(protein['raw_energy'])
        
        if not energies:
            raise RuntimeError("No valid energies calculated for any proteins")
        
        energies = np.array(energies)
        
        E_raw_weakest = np.max(energies)    
        E_raw_strongest = np.min(energies)  
        
        print(f"\n  Probability-based normalization:")
        print(f"    Target probability range: {p_min:.0%} to {p_max:.0%}")
        print(f"    Raw energy range: [{E_raw_strongest:.1f}, {E_raw_weakest:.1f}] kT")
        print(f"    Target energy range: [{E_target_for_max_prob:.2f}, {E_target_for_min_prob:.2f}]")
        
        if abs(E_raw_strongest - E_raw_weakest) < 0.01:
            p_middle = (p_min + p_max) / 2
            E_middle = np.log((1/p_middle) - 1)
            
            for protein_name, protein in self.config.proteins.items():
                if protein['raw_energy'] is not None:
                    protein['normalized_energy'] = E_middle
                    print(f"    {protein_name}: E_norm = {E_middle:.2f}, P_ads(θ=0) = {p_middle:.1%}")
            return
        
        for protein_name, protein in self.config.proteins.items():
            if protein['raw_energy'] is not None:
                E_raw = protein['raw_energy']
                
                fraction = (E_raw - E_raw_weakest) / (E_raw_strongest - E_raw_weakest)
                
                E_norm = E_target_for_min_prob + (E_target_for_max_prob - E_target_for_min_prob) * fraction
                
                protein['normalized_energy'] = E_norm
                
                p_ads_actual = 1.0 / (1 + np.exp(E_norm))
                p_des_actual = 1.0 / (1 + np.exp(-E_norm))
                
                print(f"    {protein_name}: E_raw = {E_raw:.1f} → E_norm = {E_norm:.2f}")
                print(f"      P_ads(θ=0) = {p_ads_actual:.1%}, P_des = {p_des_actual:.1%}")