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
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from pathlib import Path
from scipy.ndimage import rotate


class ProjectionCalculator:
    def __init__(self, config):
        self.config = config
    
    def rotate_pdb(self, protein_name):
        protein = self.config.proteins[protein_name]
        phi = protein.get('optimal_phi', 0)
        theta = protein.get('optimal_theta', 0)
        
        input_pdb = protein['pdb_path']
        output_pdb = self.config.proj_dir / f"{protein_name}_rotated.pdb"
        
        with open(input_pdb, 'r') as f:
            lines = f.readlines()
        
        atom_lines = [l for l in lines if l.startswith("ATOM")]
        coords = []
        for line in atom_lines:
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            coords.append([x, y, z])
        
        coords = np.array(coords)
        coords -= np.mean(coords, axis=0)
        
        phi_rad = -np.deg2rad(phi)
        theta_rad = np.pi - np.deg2rad(theta)
        
        x_phi = coords[:, 0] * np.cos(phi_rad) - coords[:, 1] * np.sin(phi_rad)
        y_phi = coords[:, 0] * np.sin(phi_rad) + coords[:, 1] * np.cos(phi_rad)
        z_phi = coords[:, 2]
        
        coords_phi = np.column_stack([x_phi, y_phi, z_phi])
        
        x_theta = coords_phi[:, 0] * np.cos(theta_rad) + coords_phi[:, 2] * np.sin(theta_rad)
        y_theta = coords_phi[:, 1]
        z_theta = -coords_phi[:, 0] * np.sin(theta_rad) + coords_phi[:, 2] * np.cos(theta_rad)
        
        rotated_coords = np.column_stack([x_theta, y_theta, z_theta])
        
        with open(output_pdb, 'w') as f:
            for i, line in enumerate(atom_lines):
                x, y, z = rotated_coords[i]
                new_line = f"{line[:30]}{x:8.3f}{y:8.3f}{z:8.3f}{line[54:]}"
                f.write(new_line)
        
        return output_pdb
    
    def parse_pdb_to_residues(self, pdb_path):
        atoms = []
        with open(pdb_path, 'r') as f:
            for line in f:
                if not line.startswith("ATOM"):
                    continue
                if len(line) > 16 and line[16] != ' ':
                    continue
                try:
                    atom = {
                        'resid': int(line[22:26].strip()),
                        'resname': line[17:20].strip(),
                        'x': float(line[30:38].strip()) / 10.0,
                        'y': float(line[38:46].strip()) / 10.0,
                        'z': float(line[46:54].strip()) / 10.0
                    }
                    atoms.append(atom)
                except ValueError:
                    continue
        
        residues = {}
        for atom in atoms:
            if atom['resid'] not in residues:
                residues[atom['resid']] = []
            residues[atom['resid']].append(atom)
        
        return residues
    
    def coarse_grain(self, residues):
        spheres = []
        
        for resid, atom_list in residues.items():
            coords = np.array([[a['x'], a['y'], a['z']] for a in atom_list])
            center = coords.mean(axis=0)
            
            distances = np.linalg.norm(coords - center, axis=1)
            radius = np.max(distances) if len(distances) > 0 else 0.3
            
            spheres.append({
                'resid': resid,
                'resname': atom_list[0]['resname'],
                'x': center[0],
                'y': center[1],
                'z': center[2],
                'radius': radius
            })
        
        return spheres
    
    def generate_mask(self, spheres):
        if not spheres:
            return np.zeros((10, 10), dtype=bool)
        
        x_coords = [s['x'] for s in spheres]
        y_coords = [s['y'] for s in spheres]
        radii = [s['radius'] for s in spheres]
        
        padding = 2.0
        x_min = min(x_coords) - max(radii) - padding
        x_max = max(x_coords) + max(radii) + padding
        y_min = min(y_coords) - max(radii) - padding
        y_max = max(y_coords) + max(radii) + padding
        
        nx = int(np.ceil((x_max - x_min) / self.config.cell_size))
        ny = int(np.ceil((y_max - y_min) / self.config.cell_size))
        
        mask = np.zeros((ny, nx), dtype=bool)
        
        for sphere in spheres:
            cx = int((sphere['x'] - x_min) / self.config.cell_size)
            cy = int((sphere['y'] - y_min) / self.config.cell_size)
            cr = int(np.ceil(sphere['radius'] / self.config.cell_size))
            
            for dy in range(-cr, cr + 1):
                for dx in range(-cr, cr + 1):
                    gx, gy = cx + dx, cy + dy
                    
                    if 0 <= gx < nx and 0 <= gy < ny:
                        dist_sq = (dx * self.config.cell_size)**2 + (dy * self.config.cell_size)**2
                        if dist_sq <= sphere['radius']**2:
                            mask[gy, gx] = True
        
        return mask
    
    def generate_rotational_masks_fixed_size(self, base_mask, n_angles=72):
        """
        Generate rotational masks with CONSISTENT sizing to prevent collision detection issues.
        All masks will have the same canvas size and approximately the same number of True cells.
        """
        rotated_masks = []
        
        h, w = base_mask.shape
        diagonal = int(np.ceil(np.sqrt(h**2 + w**2)))
        
        if diagonal % 2 == 0:
            diagonal += 1
        
        padded_size = diagonal
        padded_mask = np.zeros((padded_size, padded_size), dtype=bool)
        
        start_row = (padded_size - h) // 2
        start_col = (padded_size - w) // 2
        padded_mask[start_row:start_row+h, start_col:start_col+w] = base_mask
        
        original_cell_count = np.sum(base_mask)
        
        for i in range(n_angles):
            angle_deg = i * (360 / n_angles)  # 5° increments
            
            rotated = rotate(padded_mask.astype(float), angle_deg, 
                           reshape=False, order=0, prefilter=False)
            
            rotated_bool = rotated > 0.5
            
            rotated_cell_count = np.sum(rotated_bool)
            cell_count_diff = abs(rotated_cell_count - original_cell_count)
            
            if cell_count_diff > original_cell_count * 0.1:  
                rotated_bool = rotated > 0.7
                rotated_cell_count = np.sum(rotated_bool)
            
            rotated_masks.append(rotated_bool)
        
        print(f"        Original mask: {original_cell_count} cells")
        cell_counts = [np.sum(mask) for mask in rotated_masks]
        print(f"        Rotated masks: {min(cell_counts)}-{max(cell_counts)} cells (variation: {max(cell_counts)-min(cell_counts)})")
        
        return rotated_masks
    
    def visualize_mask(self, mask, protein_name):
        fig, ax = plt.subplots(figsize=(8, 8))
        
        ax.imshow(mask, cmap='binary', origin='lower')
        ax.set_title(f'{protein_name} Projection Mask\nResolution: {self.config.resolution} nm²/cell')
        ax.set_xlabel('Grid X (cells)')
        ax.set_ylabel('Grid Y (cells)')
        
        border_touch = (
            np.any(mask[0, :]) or np.any(mask[-1, :]) or
            np.any(mask[:, 0]) or np.any(mask[:, -1])
        )
        
        status_text = 'WARNING: Mask touches border!' if border_touch else 'Mask OK: Clear borders'
        color = 'red' if border_touch else 'green'
        ax.text(0.5, 0.95, status_text, transform=ax.transAxes, 
               ha='center', color=color, fontweight='bold', fontsize=12)
        
        plt.tight_layout()
        save_path = self.config.mask_viz_dir / f"{protein_name}_mask.png"
        plt.savefig(save_path, dpi=150)
        plt.close()
        
        return not border_touch
    
    def calculate_projection_area(self, spheres):
        if not spheres:
            return 0.0
        
        total_area = sum(np.pi * s['radius']**2 for s in spheres)
        overlap_factor = 0.85
        
        return total_area * overlap_factor
    
    def save_projection_data(self, protein_name, spheres):
        df = pd.DataFrame(spheres)
        csv_path = self.config.proj_dir / f"{protein_name}_spheres.csv"
        df.to_csv(csv_path, index=False)
        
        mask = self.generate_mask(spheres)
        mask_path = self.config.proj_dir / f"{protein_name}_mask.npy"
        np.save(mask_path, mask)
        
        return mask
    
    def process_protein(self, protein_name):
        rotated_pdb = self.rotate_pdb(protein_name)
        
        residues = self.parse_pdb_to_residues(rotated_pdb)
        spheres = self.coarse_grain(residues)
        
        base_mask = self.generate_mask(spheres)
        
        print(f"      Generating 72 rotational masks for {protein_name} (fixed size)...")
        rotated_masks = self.generate_rotational_masks_fixed_size(base_mask, n_angles=72)
        
        mask_ok = self.visualize_mask(base_mask, protein_name)
        
        area = self.calculate_projection_area(spheres)
        
        self.save_projection_data(protein_name, spheres)
        
        self.config.proteins[protein_name]['mask'] = base_mask
        self.config.proteins[protein_name]['rotated_masks'] = rotated_masks
        self.config.proteins[protein_name]['spheres'] = spheres
        self.config.proteins[protein_name]['area'] = area
        
        return area, mask_ok
    
    def calculate_all_projections(self):
        print("  Calculating 2D projections...")
        
        for protein_name in self.config.proteins:
            area, mask_ok = self.process_protein(protein_name)
            status = "✓" if mask_ok else "⚠ Border contact"
            print(f"    {protein_name}: Area = {area:.2f} nm² [{status}]")