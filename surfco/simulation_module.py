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

import numpy as np
import pandas as pd
import random
import json
import time
from pathlib import Path

from .surface import Surface
from .proteins import ProteinProbabilities


class SimulationModule:
    def __init__(self, config):
        self.config = config
        self.surface = Surface(config)
        self.probabilities = ProteinProbabilities(config)
        
        self.iteration = 0
        self.equilibrium_counter = 0
        self.previous_counts = {}
        
        self.log_data = []
        
        self.event_stats = {
            'simple_adsorption': 0,
            'adsorption_with_replacement': 0,
            'failed_adsorption_energy': 0,
            'failed_adsorption_space': 0,
            'failed_displacement_vroman': 0,
            'failed_displacement_probability': 0,
            'successful_desorption': 0,
            'failed_desorption_energy': 0,
            'total_adsorption_attempts': 0,
            'total_desorption_attempts': 0
        }
        
        self.timing_data = {
            'initialization_time': 0.0,
            'main_loop_time': 0.0,
            'equilibrium_checks_time': 0.0,
            'logging_time': 0.0,
            'save_results_time': 0.0,
            'iterations_per_second': 0.0
        }
        
        self.intermediate_log_interval = 10000
        self.last_intermediate_save = 0
    
    def find_optimal_orientation(self, protein_name, center_row, center_col):
        """Find orientation(s) with minimal collisions using fixed-size masks."""
        rotated_masks = self.config.proteins[protein_name]['rotated_masks']
        
        min_collisions = float('inf')
        optimal_orientations = []
        
        for i, mask in enumerate(rotated_masks):
            colliding = self.surface.check_collision(mask, center_row, center_col)
            n_collisions = len(colliding)
            
            if n_collisions < min_collisions:
                min_collisions = n_collisions
                optimal_orientations = [(i, mask, colliding)]
            elif n_collisions == min_collisions:
                optimal_orientations.append((i, mask, colliding))
        
        if optimal_orientations:
            selected_orientation = random.choice(optimal_orientations)
            angle_idx, best_mask, colliding_proteins = selected_orientation
            return angle_idx, best_mask, colliding_proteins
        else:
            return 0, rotated_masks[0], set()
    
    def attempt_adsorption(self, protein_name):
        """Enhanced adsorption with probability-based penalties and fixed-size masks."""
        self.event_stats['total_adsorption_attempts'] += 1
        
        if 'rotated_masks' not in self.config.proteins[protein_name]:
            return False, False, False, "no_rotated_masks"
        
        coverage = self.surface.get_coverage_fraction()
        p_binding = self.probabilities.calculate_adsorption_probability(protein_name, coverage)
        
        if random.random() > p_binding:
            self.event_stats['failed_adsorption_energy'] += 1
            return False, False, False, "binding_energy_failed"
        
        for attempt in range(20):
            center_row = random.randint(0, self.config.grid_size - 1)
            center_col = random.randint(0, self.config.grid_size - 1)
            
            angle_idx, best_mask, colliding = self.find_optimal_orientation(
                protein_name, center_row, center_col
            )
            
            if not colliding:
                success = self.surface.place_protein(protein_name, best_mask, center_row, center_col)
                if success:
                    self.event_stats['simple_adsorption'] += 1
                    return True, False, False, "simple_adsorption"
                else:
                    continue
            else:
                p_replace = self.probabilities.calculate_replacement_probability_with_probability_penalties(
                    protein_name, colliding, coverage
                )
                
                if p_replace == 0.0:
                    self.event_stats['failed_displacement_vroman'] += 1
                    continue
                
                if random.random() < p_replace:
                    self.surface.remove_proteins(colliding)
                    success = self.surface.place_protein(protein_name, best_mask, center_row, center_col)
                    if success:
                        self.event_stats['adsorption_with_replacement'] += 1
                        return True, True, True, "adsorption_with_replacement"
                    else:
                        continue
                else:
                    self.event_stats['failed_displacement_probability'] += 1
                    continue
        
        self.event_stats['failed_adsorption_space'] += 1
        return False, True, False, "no_space_available"
    
    def attempt_desorption(self, protein_name):
        """Attempt protein desorption."""
        self.event_stats['total_desorption_attempts'] += 1
        
        if self.surface.protein_counts.get(protein_name, 0) == 0:
            return False, False, False, "not_on_surface"
        
        p_des = self.probabilities.calculate_desorption_probability(protein_name)
        
        if random.random() > p_des:
            self.event_stats['failed_desorption_energy'] += 1
            return False, False, False, "desorption_energy_failed"
        
        success = self.surface.remove_random_instance(protein_name)
        
        if success:
            self.event_stats['successful_desorption'] += 1
            return True, False, False, "successful_desorption"
        else:
            return False, False, False, "removal_failed"
    
    def log_iteration(self, protein_name, event_type, success, collision, exchange, reason=""):
        """Log iteration results."""
        log_start = time.time()
        
        entry = {
            'iteration_number': self.iteration,
            'selected_protein': protein_name,
            'selected_event': event_type,
            'event_success': success,
            'collision_occurred': collision,
            'replacement_occurred': exchange,
            'event_outcome': reason,
            'surface_coverage_fraction': self.surface.get_coverage_fraction(),
            'total_surface_area': self.config.surface_area
        }
        
        for prot_name in self.config.proteins:
            entry[f'{prot_name}_count'] = self.surface.protein_counts.get(prot_name, 0)
        
        coverage_fractions = self.surface.get_protein_coverage_fractions()
        for prot_name in self.config.proteins:
            entry[f'{prot_name}_coverage'] = coverage_fractions.get(prot_name, 0.0)
        
        self.log_data.append(entry)
        self.timing_data['logging_time'] += time.time() - log_start
    
    def save_intermediate_log(self):
        """Save intermediate log data."""
        if len(self.log_data) == 0:
            return
            
        log_df = pd.DataFrame(self.log_data)
        intermediate_path = self.config.log_dir / f"simulation_log_intermediate_{self.iteration}.csv"
        log_df.to_csv(intermediate_path, index=False)
        
        intermediate_summary = {
            'current_iteration': self.iteration,
            'current_coverage': self.surface.get_coverage_fraction(),
            'current_protein_counts': dict(self.surface.protein_counts),
            'event_statistics': self.event_stats.copy(),
            'timestamp': time.time()
        }
        
        summary_path = self.config.log_dir / f"intermediate_summary_{self.iteration}.json"
        with open(summary_path, 'w') as f:
            json.dump(intermediate_summary, f, indent=2)
    
    def check_equilibrium(self):
        """Check for equilibrium."""
        eq_start = time.time()
        
        current_counts = dict(self.surface.protein_counts)
        
        if current_counts == self.previous_counts:
            self.equilibrium_counter += 1
        else:
            self.equilibrium_counter = 0
            self.previous_counts = current_counts
        
        self.timing_data['equilibrium_checks_time'] += time.time() - eq_start
        return self.equilibrium_counter >= self.config.equilibrium_window
    
    def print_event_statistics(self):
        """Print detailed event statistics."""
        total_attempts = self.event_stats['total_adsorption_attempts'] + self.event_stats['total_desorption_attempts']
        
        if total_attempts == 0:
            return
        
        print(f"\n  Event Statistics (after {self.iteration} iterations):")
        print(f"    Total event attempts: {total_attempts}")
        
        print(f"\n    Adsorption Events:")
        print(f"      Total attempts: {self.event_stats['total_adsorption_attempts']}")
        print(f"      Simple adsorption: {self.event_stats['simple_adsorption']}")
        print(f"      Adsorption with replacement: {self.event_stats['adsorption_with_replacement']}")
        print(f"      Failed (energy): {self.event_stats['failed_adsorption_energy']}")
        print(f"      Failed (coverage threshold): {self.event_stats['failed_displacement_vroman']}")
        print(f"      Failed (probability penalties): {self.event_stats['failed_displacement_probability']}")
        print(f"      Failed (no space): {self.event_stats['failed_adsorption_space']}")
        
        print(f"\n    Desorption Events:")
        print(f"      Total attempts: {self.event_stats['total_desorption_attempts']}")
        print(f"      Successful: {self.event_stats['successful_desorption']}")
        print(f"      Failed (energy): {self.event_stats['failed_desorption_energy']}")
        
        if self.event_stats['total_adsorption_attempts'] > 0:
            ads_success_rate = (self.event_stats['simple_adsorption'] + 
                               self.event_stats['adsorption_with_replacement']) / self.event_stats['total_adsorption_attempts']
            print(f"    Adsorption success rate: {ads_success_rate:.1%}")
        
        if self.event_stats['total_desorption_attempts'] > 0:
            des_success_rate = self.event_stats['successful_desorption'] / self.event_stats['total_desorption_attempts']
            print(f"    Desorption success rate: {des_success_rate:.1%}")
        
        if self.timing_data['iterations_per_second'] > 0:
            print(f"\n    Performance: {self.timing_data['iterations_per_second']:.0f} iterations/second")
    
    def run(self):
        """Run simulation with probability-based penalties and fixed-size masks."""
        init_start = time.time()
        
        print(f"  Starting simulation...")
        print(f"    Surface: {self.config.surface_area:.1f} nm²")
        print(f"    Grid: {self.config.grid_size}×{self.config.grid_size}")
        print(f"    Proteins: {len(self.config.proteins)}")
        print(f"    Rotational masks: 72 per protein (5° increments)")
        
        self.probabilities.print_concentration_analysis()
        print()
        self.probabilities.print_event_preferences()
        print()
        self.probabilities.print_displacement_settings()
        print()
        
        self.timing_data['initialization_time'] = time.time() - init_start
        
        main_loop_start = time.time()
        loop_start_time = time.time()
        
        while self.iteration < self.config.iteration_limit:
            self.iteration += 1
            
            protein_name, event_type = self.probabilities.select_protein_and_event(
                self.surface.get_coverage_fraction(), self.surface.protein_counts
            )
            
            if protein_name is None:
                self.log_iteration('None', 'none', False, False, False, "no_protein_selected")
                continue
            
            if event_type == 'adsorption':
                success, collision, exchange, reason = self.attempt_adsorption(protein_name)
            else:
                success, collision, exchange, reason = self.attempt_desorption(protein_name)
            
            self.log_iteration(protein_name, event_type, success, collision, exchange, reason)
            
            if self.iteration % 1000 == 0:
                coverage_pct = self.surface.get_coverage_fraction() * 100
                n_proteins = sum(self.surface.protein_counts.values())
                
                current_time = time.time()
                time_for_1000 = current_time - loop_start_time
                if time_for_1000 > 0:
                    current_rate = 1000 / time_for_1000
                    print(f"    Iteration {self.iteration}: {coverage_pct:.1f}% coverage, {n_proteins} proteins ({current_rate:.0f} it/s)")
                else:
                    print(f"    Iteration {self.iteration}: {coverage_pct:.1f}% coverage, {n_proteins} proteins")
                
                loop_start_time = current_time
            
            if (self.iteration - self.last_intermediate_save) >= self.intermediate_log_interval:
                self.save_intermediate_log()
                self.last_intermediate_save = self.iteration
            
            if self.check_equilibrium():
                print(f"  Equilibrium reached at iteration {self.iteration}")
                break
        
        main_loop_time = time.time() - main_loop_start
        self.timing_data['main_loop_time'] = main_loop_time
        
        if main_loop_time > 0:
            self.timing_data['iterations_per_second'] = self.iteration / main_loop_time
        
        print(f"  Simulation complete: {self.iteration} iterations")
        self.print_event_statistics()
        self._save_results()
        self._print_final_summary()
    
    def _print_final_summary(self):
        """Print final simulation state."""
        print(f"\n  Final Results:")
        print(f"    Total coverage: {self.surface.get_coverage_fraction()*100:.1f}%")
        print(f"    Total proteins: {sum(self.surface.protein_counts.values())}")
        
        coverage_fractions = self.surface.get_protein_coverage_fractions()
        for protein_name in self.config.proteins:
            count = self.surface.protein_counts.get(protein_name, 0)
            coverage = coverage_fractions.get(protein_name, 0.0) * 100
            E_norm = self.config.proteins[protein_name].get('normalized_energy', 0.0)
            print(f"    {protein_name}: {count} proteins ({coverage:.2f}% coverage, E_norm={E_norm:.2f})")
    
    def _save_results(self):
        """Save simulation results."""
        save_start = time.time()
        
        log_df = pd.DataFrame(self.log_data)
        log_path = self.config.log_dir / "simulation_log.csv"
        log_df.to_csv(log_path, index=False)
        
        matrix, protein_ids = self.surface.export_matrix()
        matrix_path = self.config.sim_dir / "final_surface.npy"
        np.save(matrix_path, matrix)
        
        id_df = pd.DataFrame(list(protein_ids.items()), columns=['protein', 'id'])
        id_path = self.config.sim_dir / "protein_ids.csv"
        id_df.to_csv(id_path, index=False)
        
        summary = {
            'total_iterations': self.iteration,
            'final_coverage': self.surface.get_coverage_fraction(),
            'total_proteins': sum(self.surface.protein_counts.values()),
            'protein_counts': dict(self.surface.protein_counts),
            'protein_coverage_fractions': self.surface.get_protein_coverage_fractions(),
            'equilibrium_reached': self.equilibrium_counter >= self.config.equilibrium_window,
            'concentration_weights': self.probabilities.calculate_concentration_weights(),
            'event_statistics': self.event_stats,
            'timing_data': self.timing_data
        }
        
        summary_path = self.config.log_dir / "summary.json"
        with open(summary_path, 'w') as f:
            json.dump(summary, f, indent=2)
        
        self.timing_data['save_results_time'] = time.time() - save_start
        
        for intermediate_file in self.config.log_dir.glob("*intermediate*"):
            intermediate_file.unlink()