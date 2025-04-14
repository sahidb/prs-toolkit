import os
import subprocess
import shutil
from typing import List, Optional, Dict
from .genome_build import GenomeBuildHandler

class PlinkInterface:
    """Interface for executing PLINK commands for PRS calculation."""
    
    def __init__(self, plink_path="plink", chain_dir="./data/chain_files"):
        """Initialize with path to PLINK executable."""
        self.plink_path = plink_path
        self._verify_plink()
        self.genome_build_handler = GenomeBuildHandler(chain_dir)
    
    def _verify_plink(self):
        """Verify PLINK executable exists and is accessible."""
        if shutil.which(self.plink_path) is None and not os.path.exists(self.plink_path):
            raise FileNotFoundError(f"PLINK executable not found at: {self.plink_path}")
    
    def run_command(self, command):
        """Run a PLINK command and capture output."""
        try:
            result = subprocess.run(
                command,
                capture_output=True,
                text=True,
                check=True
            )
            return result.stdout
        except subprocess.CalledProcessError as e:
            print(f"PLINK error: {e.stderr}")
            raise
    
    def clump(self, bfile, clump_file, output_prefix,
             p1=5e-8, p2=1.0, r2=0.1, kb=250, source_build=None, target_build=None):
        """
        Perform SNP clumping to identify independent signals.
        
        Parameters:
        -----------
        bfile : str
            Prefix for PLINK binary files
        clump_file : str
            File containing SNP IDs and p-values
        output_prefix : str
            Prefix for output files
        p1, p2 : float
            P-value thresholds for index and clumped SNPs
        r2 : float
            LD r-squared threshold
        kb : int
            Physical distance threshold in kb
        source_build, target_build : str
            Genome builds for input files, if conversion is needed
        """
        # Check if build conversion is needed
        if source_build and target_build and source_build != target_build:
            print(f"Converting clump file from {source_build} to {target_build}")
            # Convert clump file to match the bfile build
            temp_clump_file = f"{output_prefix}_temp_clump.txt"
            self.genome_build_handler.convert_variant_build(
                clump_file, temp_clump_file,
                source_build=source_build,
                target_build=target_build
            )
            clump_file = temp_clump_file
        
        command = [
            self.plink_path,
            "--bfile", bfile,
            "--clump", clump_file,
            "--clump-p1", str(p1),
            "--clump-p2", str(p2),
            "--clump-r2", str(r2),
            "--clump-kb", str(kb),
            "--out", output_prefix
        ]
        
        self.run_command(command)
        return f"{output_prefix}.clumped"
    
    def calculate_prs(self, bfile, score_file, output_prefix,
                     score_cols="1 2 3", no_mean_imputation=False,
                     extract_snps=None, keep_samples=None,
                     source_build=None, target_build=None):
        """
        Calculate Polygenic Risk Scores.
        
        Parameters:
        -----------
        bfile : str
            Prefix for PLINK binary files
        score_file : str
            Path to score file with SNP IDs, alleles, and weights
        output_prefix : str
            Prefix for output files
        score_cols : str
            Column positions for SNP ID, allele, score in score file
        no_mean_imputation : bool
            Whether to disable mean imputation for missing genotypes
        extract_snps : str
            File containing SNPs to extract
        keep_samples : str
            File containing samples to keep
        source_build, target_build : str
            Genome builds for input files, if conversion is needed
        """
        # Check if build conversion is needed for the bfile
        temp_bfile = None
        if source_build and target_build and source_build != target_build:
            print(f"Converting PLINK files from {source_build} to {target_build}")
            temp_bfile = f"{output_prefix}_temp_bfile"
            self.genome_build_handler.convert_plink_build(
                bfile, temp_bfile,
                source_build=source_build,
                target_build=target_build
            )
            bfile = temp_bfile
        
        command = [
            self.plink_path,
            "--bfile", bfile,
            "--score", score_file
        ]
        
        # Add score column specification
        command.extend(score_cols.split())
        
        # Add option to skip mean imputation for missing genotypes
        if no_mean_imputation:
            command.append("--score-no-mean-imputation")
        
        # Add extract option if provided
        if extract_snps:
            command.extend(["--extract", extract_snps])
        
        # Add keep option for population-specific analysis
        if keep_samples:
            command.extend(["--keep", keep_samples])
        
        # Add output prefix
        command.extend(["--out", output_prefix])
        
        self.run_command(command)
        
        # Clean up temporary files if created
        if temp_bfile:
            for ext in ['.bed', '.bim', '.fam']:
                temp_file = f"{temp_bfile}{ext}"
                if os.path.exists(temp_file):
                    os.unlink(temp_file)
        
        return f"{output_prefix}.profile"
    
    def calculate_multiple_prs(self, bfile, score_file, pvalue_file, 
                              range_file, output_prefix, score_cols="1 2 3",
                              extract_snps=None, keep_samples=None,
                              source_build=None, target_build=None):
        """
        Calculate multiple PRS for different p-value thresholds.
        
        Parameters:
        -----------
        bfile : str
            Prefix for PLINK binary files
        score_file : str
            Path to score file with SNP IDs, alleles, and weights
        pvalue_file : str
            File with SNP IDs and p-values
        range_file : str
            File defining score ranges and thresholds
        output_prefix : str
            Prefix for output files
        score_cols : str
            Column positions for SNP ID, allele, score in score file
        extract_snps : str
            File containing SNPs to extract
        keep_samples : str
            File containing samples to keep
        source_build, target_build : str
            Genome builds for input files, if conversion is needed
        """
        # Check if build conversion is needed for the bfile
        temp_bfile = None
        if source_build and target_build and source_build != target_build:
            print(f"Converting PLINK files from {source_build} to {target_build}")
            temp_bfile = f"{output_prefix}_temp_bfile"
            self.genome_build_handler.convert_plink_build(
                bfile, temp_bfile,
                source_build=source_build,
                target_build=target_build
            )
            bfile = temp_bfile
        
        command = [
            self.plink_path,
            "--bfile", bfile,
            "--score", score_file
        ]
        
        # Add score column specification
        command.extend(score_cols.split())
        
        # Add q-score options
        command.extend([
            "--q-score-file", pvalue_file,
            "--q-score-range", range_file
        ])
        
        # Add extract option if provided
        if extract_snps:
            command.extend(["--extract", extract_snps])
        
        # Add keep option for population-specific analysis
        if keep_samples:
            command.extend(["--keep", keep_samples])
        
        # Add output prefix
        command.extend(["--out", output_prefix])
        
        self.run_command(command)
        
        # Clean up temporary files if created
        if temp_bfile:
            for ext in ['.bed', '.bim', '.fam']:
                temp_file = f"{temp_bfile}{ext}"
                if os.path.exists(temp_file):
                    os.unlink(temp_file)
        
        # Get threshold labels from range file
        threshold_labels = []
        with open(range_file, 'r') as f:
            for line in f:
                parts = line.strip().split()
                if parts:
                    threshold_labels.append(parts[0])
        
        # Return dict mapping labels to profile files
        return {label: f"{output_prefix}.{label}.profile" for label in threshold_labels}
