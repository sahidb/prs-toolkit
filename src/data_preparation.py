import os
import pandas as pd
import numpy as np
from typing import Optional, List
from .genome_build import GenomeBuildHandler

class DataPreparation:
    """Module for preparing data files for PLINK-based PRS calculation."""
    
    def __init__(self, temp_dir="./temp", chain_dir="./data/chain_files"):
        """Initialize with directory for temporary files."""
        self.temp_dir = temp_dir
        os.makedirs(temp_dir, exist_ok=True)
        self.genome_build_handler = GenomeBuildHandler(chain_dir)
    
    def prepare_score_file(self, gwas_file, output_file=None, 
                          snp_col="SNP", effect_col="BETA", 
                          allele_col="A1", pval_col="P", 
                          pval_threshold=5e-8, source_build=None, 
                          target_build=None):
        """Create PLINK-compatible score file from GWAS summary statistics."""
        if output_file is None:
            output_file = os.path.join(self.temp_dir, f"score_p{pval_threshold}.txt")
        
        # Check if build conversion is needed
        if source_build and target_build and source_build != target_build:
            # Create temporary file for build conversion
            temp_gwas_file = os.path.join(self.temp_dir, os.path.basename(gwas_file) + ".temp")
            
            # Convert build
            self.genome_build_handler.convert_variant_build(
                gwas_file, temp_gwas_file,
                source_build=source_build,
                target_build=target_build,
                chr_col="CHR", pos_col="BP"
            )
            
            # Use converted file for further processing
            gwas_file = temp_gwas_file
        
        # Load GWAS data
        gwas_df = pd.read_csv(gwas_file, sep=None, engine='python')
        
        # Filter by p-value
        filtered_df = gwas_df[gwas_df[pval_col] <= pval_threshold].copy()
        
        # Select required columns and rename
        score_df = filtered_df[[snp_col, allele_col, effect_col]]
        
        # Write to score file (no header, tab-separated)
        score_df.to_csv(output_file, sep='\t', header=False, index=False)
        
        print(f"Created score file with {len(score_df)} SNPs at p-value <= {pval_threshold}")
        return output_file
    
    def create_range_file(self, thresholds=None, output_file=None):
        """Create a range file for PLINK's --q-score-range option."""
        if thresholds is None:
            thresholds = [5e-8, 1e-6, 1e-4, 0.001, 0.01, 0.05, 0.1, 0.5, 1.0]
        
        if output_file is None:
            output_file = os.path.join(self.temp_dir, "pvalue_ranges.txt")
        
        with open(output_file, 'w') as f:
            for i, threshold in enumerate(thresholds):
                f.write(f"S{i+1}\t0\t{threshold}\n")
        
        print(f"Created range file with {len(thresholds)} thresholds")
        return output_file
    
    def prepare_pvalue_file(self, gwas_file, output_file=None,
                           snp_col="SNP", pval_col="P", source_build=None, 
                           target_build=None):
        """Create p-value file for PLINK's --q-score-file option."""
        if output_file is None:
            output_file = os.path.join(self.temp_dir, "pval_file.txt")
        
        # Check if build conversion is needed
        if source_build and target_build and source_build != target_build:
            # Create temporary file for build conversion
            temp_gwas_file = os.path.join(self.temp_dir, os.path.basename(gwas_file) + ".temp")
            
            # Convert build
            self.genome_build_handler.convert_variant_build(
                gwas_file, temp_gwas_file,
                source_build=source_build,
                target_build=target_build,
                chr_col="CHR", pos_col="BP"
            )
            
            # Use converted file for further processing
            gwas_file = temp_gwas_file
        
        # Load GWAS data
        gwas_df = pd.read_csv(gwas_file, sep=None, engine='python')
        
        # Select required columns
        pval_df = gwas_df[[snp_col, pval_col]]
        
        # Write to p-value file (no header, tab-separated)
        pval_df.to_csv(output_file, sep='\t', header=False, index=False)
        
        print(f"Created p-value file with {len(pval_df)} SNPs")
        return output_file
    
    def extract_clumped_snps(self, clumped_file, output_file=None):
        """Extract SNP IDs from PLINK clumped file for use with --extract."""
        if output_file is None:
            output_file = os.path.join(self.temp_dir, "extract_snps.txt")
        
        # Read clumped file
        clumped_df = pd.read_csv(clumped_file, delim_whitespace=True)
        
        # Extract SNP column and write to file
        with open(output_file, 'w') as f:
            for snp in clumped_df["SNP"]:
                f.write(f"{snp}\n")
        
        print(f"Extracted {len(clumped_df)} independent SNPs from clumped file")
        return output_file
