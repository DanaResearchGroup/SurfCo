#!/usr/bin/env python3

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

import sys
import argparse
import time
import json
import logging
from pathlib import Path
from datetime import datetime

sys.path.append(str(Path(__file__).parent))

from surfco.config_manager import ConfigManager
from surfco.energy_calculator import EnergyCalculator
from surfco.projection_calculator import ProjectionCalculator
from surfco.simulation_module import SimulationModule
from surfco.visualizer import Visualizer


class TeeLogger:
    """Custom class to redirect output to both console and file."""
    def __init__(self, log_file):
        self.terminal = sys.stdout
        self.log = open(log_file, 'w')

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)
        self.log.flush()

    def flush(self):
        self.terminal.flush()
        self.log.flush()

    def close(self):
        self.log.close()


def format_time(seconds):
    """Format time in a human-readable way."""
    if seconds < 60:
        return f"{seconds:.1f}s"
    elif seconds < 3600:
        minutes = seconds / 60
        return f"{minutes:.1f}min"
    else:
        hours = seconds / 3600
        return f"{hours:.1f}h"


def show_warranty():
    """Display warranty disclaimer as per GPL v3 requirements."""
    warranty_text = """
THERE IS NO WARRANTY FOR THE PROGRAM, TO THE EXTENT PERMITTED BY
APPLICABLE LAW. EXCEPT WHEN OTHERWISE STATED IN WRITING THE COPYRIGHT
HOLDERS AND/OR OTHER PARTIES PROVIDE THE PROGRAM "AS IS" WITHOUT WARRANTY
OF ANY KIND, EITHER EXPRESSED OR IMPLIED, INCLUDING, BUT NOT LIMITED TO,
THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
PURPOSE. THE ENTIRE RISK AS TO THE QUALITY AND PERFORMANCE OF THE PROGRAM
IS WITH YOU. SHOULD THE PROGRAM PROVE DEFECTIVE, YOU ASSUME THE COST OF
ALL NECESSARY SERVICING, REPAIR OR CORRECTION.
"""
    print(warranty_text)
    sys.exit(0)


def show_conditions():
    """Display license conditions as per GPL v3 requirements."""
    conditions_text = """
SurfCo - Protein Corona Prediction Framework
Copyright (C) 2024 Dana Research Group

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

For the complete license terms, see the LICENSE file distributed with
this program or visit: https://www.gnu.org/licenses/gpl-3.0.txt

CITATION REQUIREMENTS:
When publishing results using SurfCo, you must cite:

1. Power et al. (2019) MSMSE 27(8) - UnitedAtom algorithm
2. Lopez & Lobaskin (2015) JCP 143(24) - Coarse-graining methodology

Additional material-specific citations may be required.
See LICENSE file for complete citation requirements.
"""
    print(conditions_text)
    sys.exit(0)


def print_startup_notice():
    """Print GPL v3 startup notice."""
    notice = """SurfCo  Copyright (C) 2024 Dana Research Group
This program comes with ABSOLUTELY NO WARRANTY; for details use '--show-w'.
This is free software, and you are welcome to redistribute it
under certain conditions; use '--show-c' for details.
"""
    print(notice)


def main():
    parser = argparse.ArgumentParser(
        description="SurfCo - Protein Corona Prediction Framework",
        epilog="This program is distributed under GPL v3. Use --show-c for details."
    )
    
    parser.add_argument('--show-w', '--show-warranty', action='store_true',
                       help='Show warranty disclaimer')
    parser.add_argument('--show-c', '--show-conditions', action='store_true',
                       help='Show license conditions and citation requirements')
    
    parser.add_argument('project', type=str, nargs='?',
                       help='Project name (folder in input/)')
    parser.add_argument('--skip-energy', action='store_true',
                       help='Skip energy calculations')
    parser.add_argument('--skip-projection', action='store_true',
                       help='Skip projection calculations')
    parser.add_argument('--skip-simulation', action='store_true',
                       help='Skip simulation')
    parser.add_argument('--quiet', action='store_true',
                       help='Suppress startup notice')
    
    args = parser.parse_args()
    
    if args.show_w:
        show_warranty()
    
    if args.show_c:
        show_conditions()
    
    if not args.project:
        parser.print_help()
        sys.exit(1)
    
    if not args.quiet:
        print_startup_notice()
    
    total_start_time = time.time()
    timing_data = {}
    
    print("\n" + "="*70)
    print("  SurfCo - Protein Corona Prediction Framework")
    print("  Modular KMC Simulation for Nanoparticle-Protein Interactions")
    print("="*70 + "\n")
    
    try:
        config_start = time.time()
        config = ConfigManager(args.project)
        config_time = time.time() - config_start
        timing_data['configuration'] = config_time
        
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        log_file = config.log_dir / f"terminal_output_{timestamp}.log"
        tee_logger = TeeLogger(log_file)
        sys.stdout = tee_logger
        
        print("\n" + "="*70)
        print("  SurfCo - Protein Corona Prediction Framework")
        print("  Modular KMC Simulation for Nanoparticle-Protein Interactions")
        print("="*70 + "\n")
        print(f"  Run started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        print(f"  Terminal output logged to: {log_file}\n")
        
        print(f"Project: {config.project_name}")
        print(f"NP: R={config.np_radius}nm, Material={config.material}, ζ={config.zeta_potential_mv:.1f}mV")
        print(f"Proteins: {', '.join(config.proteins.keys())}")
        print(f"Grid: {config.grid_size}×{config.grid_size} ({config.resolution} nm²/cell)")
        print(f"Configuration loaded in {format_time(config_time)}\n")
        
        if not args.skip_energy:
            print("-"*50)
            print("  Module 1: Energy Calculator")
            print("-"*50)
            energy_start = time.time()
            energy_calc = EnergyCalculator(config)
            energy_calc.calculate_all_energies()
            energy_time = time.time() - energy_start
            timing_data['energy_calculation'] = energy_time
            print(f"  Energy calculation completed in {format_time(energy_time)}")
            print()
        else:
            timing_data['energy_calculation'] = 0.0
            print("  Skipping energy calculations\n")
        
        if not args.skip_projection:
            print("-"*50)
            print("  Module 2: Projection Calculator")
            print("-"*50)
            projection_start = time.time()
            proj_calc = ProjectionCalculator(config)
            proj_calc.calculate_all_projections()
            projection_time = time.time() - projection_start
            timing_data['projection_calculation'] = projection_time
            print(f"  Projection calculation completed in {format_time(projection_time)}")
            print()
        else:
            timing_data['projection_calculation'] = 0.0
            print("  Skipping projection calculations\n")
        
        if not args.skip_simulation:
            print("-"*50)
            print("  Module 3: KMC Simulation")
            print("-"*50)
            simulation_start = time.time()
            simulation = SimulationModule(config)
            simulation.run()
            simulation_time = time.time() - simulation_start
            timing_data['simulation'] = simulation_time
            if hasattr(simulation, 'timing_data'):
                timing_data['simulation_details'] = simulation.timing_data
            print(f"  Simulation completed in {format_time(simulation_time)}")
            print()
            
            print("-"*50)
            print("  Module 4: Analysis & Visualization")
            print("-"*50)
            viz_start = time.time()
            visualizer = Visualizer(config)
            visualizer.create_all_visualizations()
            viz_time = time.time() - viz_start
            timing_data['visualization'] = viz_time
            print(f"  Visualization completed in {format_time(viz_time)}")
        else:
            timing_data['simulation'] = 0.0
            timing_data['visualization'] = 0.0
            print("  Skipping simulation and visualization\n")
        
        total_time = time.time() - total_start_time
        timing_data['total_runtime'] = total_time
        
        timing_path = config.log_dir / "timing_summary.json"
        with open(timing_path, 'w') as f:
            json.dump(timing_data, f, indent=2)
        
        print("\n" + "="*70)
        print(f"  Complete! Results saved to: projects/{config.project_name}/")
        print(f"  Total runtime: {format_time(total_time)}")
        print("="*70)
        
        print("\n  Timing Breakdown:")
        print(f"    Configuration: {format_time(timing_data['configuration'])}")
        print(f"    Energy calculation: {format_time(timing_data['energy_calculation'])}")
        print(f"    Projection calculation: {format_time(timing_data['projection_calculation'])}")
        print(f"    Simulation: {format_time(timing_data['simulation'])}")
        print(f"    Visualization: {format_time(timing_data['visualization'])}")
        print(f"    Total: {format_time(timing_data['total_runtime'])}")
        
        print(f"\n  Run completed: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        print(f"  Terminal output saved to: {log_file}\n")
        
        sys.stdout = tee_logger.terminal
        tee_logger.close()
        
    except Exception as e:
        print(f"\nError: {e}")
        import traceback
        traceback.print_exc()
        
        if 'tee_logger' in locals():
            sys.stdout = tee_logger.terminal
            tee_logger.close()
        
        sys.exit(1)


if __name__ == "__main__":
    main()