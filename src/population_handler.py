import os
import pandas as pd
import numpy as np
from typing import List, Dict, Optional
from .genome_build import GenomeBuildHandler

class PopulationHandler:
    """Module for managing population-specific analyses."""
    
    def __init__(self, temp_dir="./temp", chain_dir="./data/chain_files"):
        """Initialize with directory for temporary files."""
        self.temp_dir = temp_dir
        os.makedirs(temp_dir, exist_ok=True)
        self.genome_build_handler = GenomeBuildHandler(chain_dir)
    
    def load_population_data(self, population_file):
        """Load population mapping data from file."""
        try:
            # Load population file
            pop_df = pd.read_csv(population_file, sep=None, engine='python')
            required_cols = ['Individual', 'Population']
            if not all(col in pop_df.columns for col in required_cols):
                raise ValueError(f"Population file must contain columns: {required_cols}")
            return pop_df
        except Exception as e:
            raise Exception(f"Error loading population data: {str(e)}")
    
    def create_population_keep_file(self, population_data, population, output_file=None):
        """Create a PLINK-compatible keep file for a specific population."""
        if output_file is None:
            output_file = os.path.join(self.temp_dir, f"keep_{population}.txt")
        
        # Filter by population
        pop_subset = population_data[population_data['Population'] == population]
        
        # Write keep file with FID and IID columns (PLINK format)
        with open(output_file, 'w') as f:
            for _, row in pop_subset.iterrows():
                # If FID is available, use it, otherwise use Individual as both FID and IID
                if 'FID' in population_data.columns:
                    f.write(f"{row['FID']}\t{row['Individual']}\n")
                else:
                    f.write(f"{row['Individual']}\t{row['Individual']}\n")
        
        print(f"Created keep file for {population} with {len(pop_subset)} individuals")
        return output_file
    
    def summarize_population_data(self, population_data):
        """Generate summary statistics for available populations."""
        # Count individuals per population
        pop_counts = population_data['Population'].value_counts().reset_index()
        pop_counts.columns = ['Population', 'Count']
        
        print("Population summary:")
        for _, row in pop_counts.iterrows():
            print(f"  {row['Population']}: {row['Count']} individuals")
        
        return pop_counts
    
    def calculate_population_statistics(self, profile_file, output_file=None):
        """
        Calculate statistics for a population PRS distribution.
        
        Parameters:
        -----------
        profile_file : str
            Path to PLINK profile file for population
        output_file : str, optional
            Path to save statistics
            
        Returns:
        --------
        Dict with population statistics
        """
        # Load profile data
        pop_df = pd.read_csv(profile_file, delim_whitespace=True)
        
        # Calculate distribution statistics
        stats = {
            'mean': pop_df['SCORE'].mean(),
            'std': pop_df['SCORE'].std(),
            'median': pop_df['SCORE'].median(),
            'percentiles': {
                p: np.percentile(pop_df['SCORE'], p) 
                for p in [1, 5, 10, 20, 25, 50, 75, 80, 90, 95, 99]
            },
            'n_samples': len(pop_df)
        }
        
        # Save to file if specified
        if output_file:
            with open(output_file, 'w') as f:
                f.write(f"Population Statistics\n")
                f.write(f"-------------------\n")
                f.write(f"Number of samples: {stats['n_samples']}\n")
                f.write(f"Mean: {stats['mean']:.4f}\n")
                f.write(f"Standard deviation: {stats['std']:.4f}\n")
                f.write(f"Median: {stats['median']:.4f}\n")
                f.write(f"\nPercentiles:\n")
                for p, val in stats['percentiles'].items():
                    f.write(f"{p}th: {val:.4f}\n")
        
        return stats
