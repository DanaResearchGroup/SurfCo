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
from collections import defaultdict


class Surface:
    def __init__(self, config):
        self.config = config
        self.grid_size = config.grid_size
        self.resolution = config.resolution
        self.total_area = config.surface_area
        
        self.occupied_cells = {}
        self.protein_instances = {}
        self.next_instance_id = 1
        
        self.protein_counts = defaultdict(int)
        self.protein_areas = defaultdict(float)
    
    def check_collision(self, mask, center_row, center_col):
        h, w = mask.shape
        colliding_proteins = set()
        
        start_row = center_row - h // 2
        start_col = center_col - w // 2
        
        for dy in range(h):
            for dx in range(w):
                if not mask[dy, dx]:
                    continue
                
                row = (start_row + dy) % self.grid_size
                col = (start_col + dx) % self.grid_size
                
                if (row, col) in self.occupied_cells:
                    protein_name, _ = self.occupied_cells[(row, col)]
                    colliding_proteins.add(protein_name)
        
        return colliding_proteins
    
    def place_protein(self, protein_name, mask, center_row, center_col):
        h, w = mask.shape
        instance_id = self.next_instance_id
        self.next_instance_id += 1
        
        cells = []
        start_row = center_row - h // 2
        start_col = center_col - w // 2
        
        for dy in range(h):
            for dx in range(w):
                if not mask[dy, dx]:
                    continue
                
                row = (start_row + dy) % self.grid_size
                col = (start_col + dx) % self.grid_size
                
                self.occupied_cells[(row, col)] = (protein_name, instance_id)
                cells.append((row, col))
        
        self.protein_instances[instance_id] = (protein_name, cells)
        
        self.protein_counts[protein_name] += 1
        self.protein_areas[protein_name] += len(cells) * self.resolution
        
        return instance_id
    
    def remove_proteins(self, protein_names):
        instances_to_remove = []
        
        for inst_id, (prot_name, cells) in self.protein_instances.items():
            if prot_name in protein_names:
                instances_to_remove.append(inst_id)
                
                for cell in cells:
                    if cell in self.occupied_cells:
                        del self.occupied_cells[cell]
                
                self.protein_counts[prot_name] -= 1
                self.protein_areas[prot_name] -= len(cells) * self.resolution
        
        for inst_id in instances_to_remove:
            del self.protein_instances[inst_id]
        
        for name in list(self.protein_counts.keys()):
            if self.protein_counts[name] <= 0:
                self.protein_counts[name] = 0
                self.protein_areas[name] = 0.0
    
    def remove_random_instance(self, protein_name):
        instances = [
            inst_id for inst_id, (prot_name, _) in self.protein_instances.items()
            if prot_name == protein_name
        ]
        
        if not instances:
            return False
        
        inst_id = np.random.choice(instances)
        _, cells = self.protein_instances[inst_id]
        
        for cell in cells:
            if cell in self.occupied_cells:
                del self.occupied_cells[cell]
        
        self.protein_counts[protein_name] -= 1
        self.protein_areas[protein_name] -= len(cells) * self.resolution
        
        del self.protein_instances[inst_id]
        
        return True
    
    def get_coverage_fraction(self):
        occupied_area = len(self.occupied_cells) * self.resolution
        return occupied_area / self.total_area
    
    def get_protein_coverage_fractions(self):
        fractions = {}
        for protein_name, area in self.protein_areas.items():
            fractions[protein_name] = area / self.total_area
        return fractions
    
    def export_matrix(self):
        matrix = np.zeros((self.grid_size, self.grid_size), dtype=int)
        
        protein_ids = {}
        current_id = 1
        
        for (row, col), (protein_name, _) in self.occupied_cells.items():
            if protein_name not in protein_ids:
                protein_ids[protein_name] = current_id
                current_id += 1
            
            matrix[row, col] = protein_ids[protein_name]
        
        return matrix, protein_ids