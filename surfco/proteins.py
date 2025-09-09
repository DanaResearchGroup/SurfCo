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


class ProteinProbabilities:
    def __init__(self, config):
        self.config = config
        self.proteins = config.proteins
    
    def calculate_adsorption_probability(self, protein_name, coverage_fraction):
        """Calculate probability that a protein successfully binds upon encounter - ENERGY ONLY."""
        if protein_name not in self.proteins:
            return 0.0
        
        E_norm = self.proteins[protein_name].get('normalized_energy', 0.0)
        
        energy_term = 1 + np.exp(E_norm)
        probability = 1.0 / energy_term
        
        return np.clip(probability, 0.001, 0.999)
    
    def calculate_desorption_probability(self, protein_name):
        """Calculate probability that a protein desorbs from the surface."""
        if protein_name not in self.proteins:
            return 0.0
        
        E_norm = self.proteins[protein_name].get('normalized_energy', 0.0)
        
        probability = 1.0 / (1 + np.exp(-E_norm))
        
        return np.clip(probability, 0.001, 0.999)
    
    def calculate_coverage_dependent_displacement_probability(self, coverage_fraction):
        """Calculate coverage-dependent displacement probability (Vroman effect)."""
        threshold = getattr(self.config, 'coverage_threshold_for_displacement', 0.15)
        max_prob = getattr(self.config, 'max_displacement_probability', 0.8)
        
        if coverage_fraction <= threshold:
            return 0.0
        
        normalized_coverage = (coverage_fraction - threshold) / (1.0 - threshold)
        displacement_prob = max_prob * normalized_coverage
        
        return np.clip(displacement_prob, 0.0, max_prob)
    
    def calculate_replacement_probability_with_probability_penalties(self, incoming_name, existing_names, coverage_fraction):
        """Enhanced replacement probability using probability-based penalties."""
        if incoming_name not in self.proteins:
            return 0.0
        
        if not existing_names:
            return 0.0
        
        vroman_prob = self.calculate_coverage_dependent_displacement_probability(coverage_fraction)
        if vroman_prob == 0.0:
            return 0.0
        
        E_new = self.proteins[incoming_name].get('normalized_energy', 0.0)
        E_old_total = sum(self.proteins[name].get('normalized_energy', 0.0) 
                         for name in existing_names if name in self.proteins)
        
        delta_E = E_new - E_old_total
        
        if delta_E >= 0:
            base_probability = np.exp(-delta_E)
        else:
            base_probability = 1.0 / (1.0 + np.exp(delta_E))
        
        n_displaced = len(existing_names)
        
        penalty_per_protein = getattr(self.config, 'probability_penalty_per_extra_protein', 0.25)
        min_prob = getattr(self.config, 'min_displacement_probability', 0.05)
        base_penalty = getattr(self.config, 'base_displacement_penalty', 0.15)
        
        penalized_probability = base_probability * (1.0 - base_penalty)
        
        extra_proteins = max(0, n_displaced - 1)
        total_extra_penalty = penalty_per_protein * extra_proteins
        
        final_probability = penalized_probability * (1.0 - total_extra_penalty)
        final_probability = max(final_probability, min_prob)
        final_probability = min(final_probability, 0.3)
        
        return final_probability
    
    def select_protein_based_on_concentration(self):
        """Select protein based on concentration."""
        protein_names = list(self.proteins.keys())
        concentrations = [self.proteins[name]['concentration'] for name in protein_names]
        
        total_concentration = sum(concentrations)
        if total_concentration <= 0:
            return None
        
        probabilities = [c / total_concentration for c in concentrations]
        selected_idx = np.random.choice(len(protein_names), p=probabilities)
        return protein_names[selected_idx]
    
    def select_event_type(self, protein_name, surface_counts):
        """Select event type using user-configured probability."""
        proteins_on_surface = surface_counts.get(protein_name, 0)
        
        if proteins_on_surface == 0:
            return 'adsorption'
        
        p_adsorption = self.config.adsorption_event_probability
        return 'adsorption' if np.random.random() < p_adsorption else 'desorption'
    
    def select_protein_and_event(self, coverage_fraction, surface_counts):
        """Two-step selection process."""
        selected_protein = self.select_protein_based_on_concentration()
        if selected_protein is None:
            return None, None
        
        event_type = self.select_event_type(selected_protein, surface_counts)
        return selected_protein, event_type
    
    def calculate_event_probabilities(self, protein_name, surface_counts):
        """Calculate event selection probabilities."""
        proteins_on_surface = surface_counts.get(protein_name, 0)
        
        if proteins_on_surface == 0:
            return 1.0, 0.0
        
        p_adsorption = self.config.adsorption_event_probability
        return p_adsorption, 1.0 - p_adsorption
    
    def calculate_concentration_weights(self):
        """Calculate relative concentration weights."""
        weights = {}
        concentrations = []
        
        for protein_name, protein_data in self.proteins.items():
            concentrations.append(protein_data['concentration'])
        
        total_concentration = sum(concentrations)
        
        for protein_name, protein_data in self.proteins.items():
            concentration = protein_data['concentration']
            relative_weight = concentration / total_concentration if total_concentration > 0 else 0
            
            weights[protein_name] = {
                'concentration': concentration,
                'selection_probability': relative_weight,
                'relative_frequency': concentration / min(concentrations) if min(concentrations) > 0 else 1
            }
        
        return weights
    
    def print_concentration_analysis(self):
        """Print concentration analysis."""
        weights = self.calculate_concentration_weights()
        
        print("  Protein selection probabilities:")
        total_conc = sum(w['concentration'] for w in weights.values())
        
        for protein_name, data in weights.items():
            conc = data['concentration']
            prob = data['selection_probability']
            freq = data['relative_frequency']
            
            print(f"    {protein_name}: {conc:.1f} mg/mL ({prob:.1%} selection chance, {freq:.1f}x relative frequency)")
        
        print(f"    Total concentration: {total_conc:.1f} mg/mL")
    
    def print_event_preferences(self):
        """Print event selection analysis."""
        baseline_prob = self.config.adsorption_event_probability
        print(f"  Event selection probabilities:")
        print(f"    Adsorption probability: {baseline_prob:.1%} (user-defined)")
        print(f"    Desorption probability: {1.0-baseline_prob:.1%} (user-defined)")
        
        surface_counts = {name: 5 for name in self.proteins.keys()}
        
        for protein_name in self.proteins.keys():
            E_norm = self.proteins[protein_name].get('normalized_energy', 0.0)
            print(f"    {protein_name}: E_norm = {E_norm:.2f}")
    
    def print_displacement_settings(self):
        """Print displacement control settings."""
        penalty_per_protein = getattr(self.config, 'probability_penalty_per_extra_protein', 0.25)
        min_prob = getattr(self.config, 'min_displacement_probability', 0.05)
        base_penalty = getattr(self.config, 'base_displacement_penalty', 0.15)
        threshold = getattr(self.config, 'coverage_threshold_for_displacement', 0.15)
        max_prob = getattr(self.config, 'max_displacement_probability', 0.8)
        
        print(f"  Displacement control parameters:")
        print(f"    Coverage threshold: {threshold:.1%} (user-defined)")
        print(f"    Maximum displacement probability: {max_prob:.1%} (user-defined)")
        print(f"    Base displacement penalty: {base_penalty:.1%} (user-defined)")
        print(f"    Penalty per extra protein: {penalty_per_protein:.1%} (user-defined)")
        print(f"    Minimum displacement probability: {min_prob:.1%} (user-defined)")