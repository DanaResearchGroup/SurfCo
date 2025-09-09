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
import json
from pathlib import Path


class Visualizer:
    def __init__(self, config):
        self.config = config
        
        plt.rcParams['figure.dpi'] = 100
        plt.rcParams['savefig.dpi'] = 300
        plt.rcParams['font.size'] = 11
        plt.rcParams['axes.labelsize'] = 12
        plt.rcParams['axes.titlesize'] = 14
        plt.rcParams['legend.fontsize'] = 10
    
    def smooth_data(self, df, column, window_fraction=0.02):
        """Apply smoothing to data for cleaner visualization"""
        if len(df) < 10:
            return df[column]
        
        window_size = max(int(len(df) * window_fraction), 5)
        window_size = min(window_size, len(df) // 4)
        
        return df[column].rolling(window=window_size, center=True, min_periods=1).mean()
    
    def subsample_data(self, df, max_points=2000):
        """Subsample data for plotting while preserving trends"""
        if len(df) <= max_points:
            return df
        
        step = len(df) // max_points
        return df.iloc[::step].copy()
    
    def create_coverage_plots(self):
        log_path = self.config.log_dir / "simulation_log.csv"
        if not log_path.exists():
            print("    No simulation log found")
            return
        
        df = pd.read_csv(log_path)
        protein_names = list(self.config.proteins.keys())
        
        colors = plt.cm.Set2(np.linspace(0, 1, len(protein_names)))
        
        self._plot_early_dynamics(df, protein_names, colors)
        self._plot_late_equilibrium(df, protein_names, colors)
        self._plot_complete_dynamics(df, protein_names, colors)
        
        print("    Coverage plots created")
    
    def _plot_early_dynamics(self, df, proteins, colors):
        cutoff = min(len(df) // 10, 10000)
        early_df = df.iloc[:cutoff].copy()
        
        plot_df = self.subsample_data(early_df, max_points=1000)
        
        fig, ax = plt.subplots(figsize=(10, 6))
        
        for i, protein in enumerate(proteins):
            col = f'{protein}_coverage'
            if col in plot_df.columns:
                smooth_coverage = self.smooth_data(plot_df, col, window_fraction=0.05)
                
                ax.plot(plot_df['iteration_number'], smooth_coverage, 
                       label=protein, color=colors[i], linewidth=2.5, alpha=0.9)
        
        ax.set_xlabel('Iteration Number')
        ax.set_ylabel('Surface Coverage Fraction')
        ax.set_title('Protein Corona Formation - Early Dynamics (Smoothed)')
        ax.legend(loc='best', frameon=True, fancybox=True, shadow=True)
        ax.grid(True, alpha=0.3, linestyle='--')
        ax.set_ylim(bottom=0)
        
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        
        plt.tight_layout()
        save_path = self.config.coverage_viz_dir / "coverage_early_smooth.png"
        plt.savefig(save_path, bbox_inches='tight')
        plt.close()
    
    def _plot_late_equilibrium(self, df, proteins, colors):
        cutoff = len(df) // 2
        late_df = df.iloc[cutoff:].copy()
        
        plot_df = self.subsample_data(late_df, max_points=1000)
        
        fig, ax = plt.subplots(figsize=(10, 6))
        
        for i, protein in enumerate(proteins):
            col = f'{protein}_coverage'
            if col in plot_df.columns:
                smooth_coverage = self.smooth_data(plot_df, col, window_fraction=0.03)
                
                ax.plot(plot_df['iteration_number'], smooth_coverage, 
                       label=protein, color=colors[i], linewidth=2.5, alpha=0.9)
        
        ax.set_xlabel('Iteration Number')
        ax.set_ylabel('Surface Coverage Fraction')
        ax.set_title('Protein Corona Formation - Equilibrium Phase (Smoothed)')
        ax.legend(loc='best', frameon=True, fancybox=True, shadow=True)
        ax.grid(True, alpha=0.3, linestyle='--')
        
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        
        plt.tight_layout()
        save_path = self.config.coverage_viz_dir / "coverage_equilibrium_smooth.png"
        plt.savefig(save_path, bbox_inches='tight')
        plt.close()
    
    def _plot_complete_dynamics(self, df, proteins, colors):
        plot_df = self.subsample_data(df, max_points=2500)
        
        fig, ax = plt.subplots(figsize=(12, 7))
        
        for i, protein in enumerate(proteins):
            col = f'{protein}_coverage'
            if col in plot_df.columns:
                smooth_coverage = self.smooth_data(plot_df, col, window_fraction=0.02)
                
                ax.plot(plot_df['iteration_number'], smooth_coverage, 
                       label=protein, color=colors[i], linewidth=2.5, alpha=0.9)
        
        ax.set_xlabel('Iteration Number')
        ax.set_ylabel('Surface Coverage Fraction')
        ax.set_title('Protein Corona Formation - Complete Dynamics (Smoothed)')
        ax.legend(loc='best', frameon=True, fancybox=True, shadow=True)
        ax.grid(True, alpha=0.3, linestyle='--')
        ax.set_ylim(bottom=0)
        
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        
        final_coverage = plot_df['surface_coverage_fraction'].iloc[-1]
        total_iterations = plot_df['iteration_number'].iloc[-1]
        
        stats_text = f'Final Coverage: {final_coverage:.1%}\nIterations: {total_iterations:,}'
        ax.text(0.02, 0.98, stats_text, transform=ax.transAxes, 
               verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
        
        plt.tight_layout()
        save_path = self.config.coverage_viz_dir / "coverage_complete_smooth.png"
        plt.savefig(save_path, bbox_inches='tight')
        plt.close()
    
    def create_cumulative_coverage_plot(self, df, proteins, colors):
        """Create a cumulative coverage plot showing protein accumulation over time"""
        plot_df = self.subsample_data(df, max_points=2000)
        
        fig, ax = plt.subplots(figsize=(12, 7))
        
        for i, protein in enumerate(proteins):
            count_col = f'{protein}_count'
            if count_col in plot_df.columns:
                smooth_counts = self.smooth_data(plot_df, count_col, window_fraction=0.03)
                
                ax.plot(plot_df['iteration_number'], smooth_counts, 
                       label=f'{protein} (count)', color=colors[i], linewidth=2.5, alpha=0.9)
        
        ax.set_xlabel('Iteration Number')
        ax.set_ylabel('Number of Proteins on Surface')
        ax.set_title('Protein Accumulation Over Time (Smoothed)')
        ax.legend(loc='best', frameon=True, fancybox=True, shadow=True)
        ax.grid(True, alpha=0.3, linestyle='--')
        ax.set_ylim(bottom=0)
        
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        
        plt.tight_layout()
        save_path = self.config.coverage_viz_dir / "protein_accumulation_smooth.png"
        plt.savefig(save_path, bbox_inches='tight')
        plt.close()
        
        print("    Protein accumulation plot created")
    
    def create_final_composition_plot(self):
        summary_path = self.config.log_dir / "summary.json"
        if not summary_path.exists():
            return
        
        with open(summary_path, 'r') as f:
            summary = json.load(f)
        
        counts = summary.get('protein_counts', {})
        if not counts or sum(counts.values()) == 0:
            return
        
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
        
        proteins = list(counts.keys())
        values = list(counts.values())
        colors = plt.cm.Set2(np.linspace(0, 1, len(proteins)))
        
        bars = ax1.bar(proteins, values, color=colors, edgecolor='black', linewidth=1.5, alpha=0.8)
        ax1.set_xlabel('Protein')
        ax1.set_ylabel('Number of Molecules')
        ax1.set_title('Final Corona Composition - Molecular Count')
        ax1.grid(True, alpha=0.3, axis='y', linestyle='--')
        
        for bar, val in zip(bars, values):
            height = bar.get_height()
            ax1.text(bar.get_x() + bar.get_width()/2., height + 0.5,
                    f'{int(val)}', ha='center', va='bottom', fontweight='bold', fontsize=11)
        
        wedges, texts, autotexts = ax2.pie(values, labels=proteins, colors=colors,
                                            autopct='%1.1f%%', startangle=90,
                                            textprops={'fontsize': 11}, 
                                            wedgeprops=dict(edgecolor='white', linewidth=2))
        ax2.set_title('Final Corona Composition - Relative Abundance')
        
        for autotext in autotexts:
            autotext.set_fontweight('bold')
            autotext.set_color('white')
        
        plt.tight_layout()
        save_path = self.config.coverage_viz_dir / "final_composition.png"
        plt.savefig(save_path, bbox_inches='tight')
        plt.close()
        
        print("    Final composition plot created")
    
    def create_surface_visualization(self):
        matrix_path = self.config.sim_dir / "final_surface.npy"
        if not matrix_path.exists():
            return
        
        matrix = np.load(matrix_path)
        
        id_path = self.config.sim_dir / "protein_ids.csv"
        if id_path.exists():
            id_df = pd.read_csv(id_path)
            id_map = dict(zip(id_df['id'], id_df['protein']))
        else:
            id_map = {}
        
        fig, ax = plt.subplots(figsize=(10, 10))
        
        from matplotlib.colors import ListedColormap
        n_proteins = len(np.unique(matrix[matrix > 0]))
        base_colors = plt.cm.Set2(np.linspace(0, 1, max(n_proteins, 1)))
        colors = ['white'] + list(base_colors[:n_proteins])
        cmap = ListedColormap(colors)
        
        im = ax.imshow(matrix, cmap=cmap, interpolation='nearest')
        
        ax.set_title(f'Final Surface State\nNP: R={self.config.np_radius}nm, {self.config.material}', fontsize=14)
        ax.set_xlabel('Grid X (cells)')
        ax.set_ylabel('Grid Y (cells)')
        
        cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        cbar.set_label('Protein Type', fontsize=12)
        
        if id_map:
            cbar.set_ticks(list(id_map.keys()))
            cbar.set_ticklabels(list(id_map.values()))
        
        plt.tight_layout()
        save_path = self.config.surface_viz_dir / "final_surface.png"
        plt.savefig(save_path, bbox_inches='tight')
        plt.close()
        
        print("    Surface visualization created")
    
    def create_all_visualizations(self):
        print("  Creating visualizations...")
        
        log_path = self.config.log_dir / "simulation_log.csv"
        if log_path.exists():
            df = pd.read_csv(log_path)
            protein_names = list(self.config.proteins.keys())
            colors = plt.cm.Set2(np.linspace(0, 1, len(protein_names)))
            
            self.create_cumulative_coverage_plot(df, protein_names, colors)
        
        self.create_coverage_plots()
        self.create_final_composition_plot()
        self.create_surface_visualization()
        print("  Visualizations complete")