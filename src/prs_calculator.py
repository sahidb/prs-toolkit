import os
from typing import Dict, List, Optional, Any
from .genome_build import GenomeBuildHandler

class PRSCalculator:
    """Core calculator class for PRS computation workflows."""
    
    def __init__(self, data_prep, plink, output_dir="./output", chain_dir="./data/chain_files"):
        """Initialize with data preparation and PLINK interface objects."""
        self.data_prep = data_prep
        self.plink = plink
        self.output_dir = output_dir
        os.makedirs(output_dir, exist_ok=True)
        self.genome_build_handler = GenomeBuildHandler(chain_dir)
    
    def basic_prs(self, gwas_file, bfile, output_name,
                 snp_col="SNP", effect_col="BETA", allele_col="A1", pval_col="P",
                 pval_threshold=5e-8, keep_samples=None, 
                 source_build=None, target_build=None):
        """
        Calculate basic PRS using GWAS summary statistics.
        
        Parameters:
        -----------
        gwas_file : str
            Path to GWAS summary statistics file
        bfile : str
            Prefix for PLINK binary files
        output_name : str
            Prefix for output files
        snp_col, effect_col, allele_col, pval_col : str
            Column names in GWAS file
        pval_threshold : float
            P-value threshold for SNP selection
        keep_samples : str
            File containing samples to keep
        source_build, target_build : str
            Genome builds for GWAS and genotype data
            
        Returns:
        --------
        str : Path to the profile file
        """
        # Create output prefix
        output_prefix = os.path.join(self.output_dir, output_name)
        
        # Detect or validate builds if needed
        gwas_build = source_build
        genotype_build = target_build
        
        if gwas_build is None:
            gwas_build = self.genome_build_handler.detect_genome_build(
                gwas_file, chr_col="CHR", pos_col="BP"
            )
            print(f"Detected GWAS build: {gwas_build}")
        
        if genotype_build is None:
            # Try to detect from .bim file
            bim_file = f"{bfile}.bim"
            if os.path.exists(bim_file):
                genotype_build = self.genome_build_handler.detect_genome_build(
                    bim_file, chr_col=0, pos_col=3
                )
                print(f"Detected genotype build: {genotype_build}")
        
        # Prepare score file with build conversion if needed
        score_file = self.data_prep.prepare_score_file(
            gwas_file=gwas_file,
            snp_col=snp_col,
            effect_col=effect_col,
            allele_col=allele_col,
            pval_col=pval_col,
            pval_threshold=pval_threshold,
            source_build=gwas_build,
            target_build=genotype_build
        )
        
        # Calculate PRS
        profile_file = self.plink.calculate_prs(
            bfile=bfile,
            score_file=score_file,
            output_prefix=output_prefix,
            keep_samples=keep_samples,
            source_build=gwas_build,
            target_build=genotype_build
        )
        
        return profile_file
    
    def clumped_prs(self, gwas_file, bfile, output_name,
                   clump_p1=1e-5, clump_r2=0.1, clump_kb=250,
                   snp_col="SNP", effect_col="BETA", allele_col="A1", pval_col="P",
                   pval_threshold=5e-8, keep_samples=None,
                   source_build=None, target_build=None):
        """
        Calculate PRS using clumped SNPs to reduce LD.
        
        Parameters as in basic_prs, with additional clumping parameters.
        """
        # Create base output prefix
        base_prefix = os.path.join(self.output_dir, output_name)
        clump_prefix = f"{base_prefix}_clumped"
        
        # Detect or validate builds if needed
        gwas_build = source_build
        genotype_build = target_build
        
        if gwas_build is None:
            gwas_build = self.genome_build_handler.detect_genome_build(
                gwas_file, chr_col="CHR", pos_col="BP"
            )
            print(f"Detected GWAS build: {gwas_build}")
        
        if genotype_build is None:
            # Try to detect from .bim file
            bim_file = f"{bfile}.bim"
            if os.path.exists(bim_file):
                genotype_build = self.genome_build_handler.detect_genome_build(
                    bim_file, chr_col=0, pos_col=3
                )
                print(f"Detected genotype build: {genotype_build}")
        
        # Prepare p-value file for clumping with build conversion if needed
        pval_file = self.data_prep.prepare_pvalue_file(
            gwas_file=gwas_file,
            snp_col=snp_col,
            pval_col=pval_col,
            source_build=gwas_build,
            target_build=genotype_build
        )
        
        # Perform clumping
        clumped_file = self.plink.clump(
            bfile=bfile,
            clump_file=pval_file,
            output_prefix=clump_prefix,
            p1=clump_p1,
            r2=clump_r2,
            kb=clump_kb,
            source_build=gwas_build,
            target_build=genotype_build
        )
        
        # Extract clumped SNPs
        extract_file = self.data_prep.extract_clumped_snps(clumped_file)
        
        # Prepare score file with build conversion if needed
        score_file = self.data_prep.prepare_score_file(
            gwas_file=gwas_file,
            snp_col=snp_col,
            effect_col=effect_col,
            allele_col=allele_col,
            pval_col=pval_col,
            pval_threshold=pval_threshold,
            source_build=gwas_build,
            target_build=genotype_build
        )
        
        # Calculate PRS using clumped SNPs
        profile_file = self.plink.calculate_prs(
            bfile=bfile,
            score_file=score_file,
            output_prefix=f"{base_prefix}_final",
            extract_snps=extract_file,
            keep_samples=keep_samples,
            source_build=gwas_build,
            target_build=genotype_build
        )
        
        return profile_file
    
    def multi_threshold_prs(self, gwas_file, bfile, output_name,
                           thresholds=None,
                           snp_col="SNP", effect_col="BETA", allele_col="A1", pval_col="P",
                           keep_samples=None, use_clumping=True,
                           clump_p1=1e-5, clump_r2=0.1, clump_kb=250,
                           source_build=None, target_build=None):
        """
        Calculate PRS at multiple p-value thresholds.
        
        Parameters as in clumped_prs, with additional parameters for thresholds.
        """
        # Create base output prefix
        base_prefix = os.path.join(self.output_dir, output_name)
        
        # Set default thresholds if none provided
        if thresholds is None:
            thresholds = [5e-8, 1e-6, 1e-4, 0.001, 0.01, 0.05, 0.1, 0.5, 1.0]
        
        # Detect or validate builds if needed
        gwas_build = source_build
        genotype_build = target_build
        
        if gwas_build is None:
            gwas_build = self.genome_build_handler.detect_genome_build(
                gwas_file, chr_col="CHR", pos_col="BP"
            )
            print(f"Detected GWAS build: {gwas_build}")
        
        if genotype_build is None:
            # Try to detect from .bim file
            bim_file = f"{bfile}.bim"
            if os.path.exists(bim_file):
                genotype_build = self.genome_build_handler.detect_genome_build(
                    bim_file, chr_col=0, pos_col=3
                )
                print(f"Detected genotype build: {genotype_build}")
        
        # Prepare p-value range file
        range_file = self.data_prep.create_range_file(thresholds)
        
        # Prepare p-value file for q-score with build conversion if needed
        pval_file = self.data_prep.prepare_pvalue_file(
            gwas_file=gwas_file,
            snp_col=snp_col,
            pval_col=pval_col,
            source_build=gwas_build,
            target_build=genotype_build
        )
        
        # Prepare score file with build conversion if needed
        score_file = self.data_prep.prepare_score_file(
            gwas_file=gwas_file,
            snp_col=snp_col,
            effect_col=effect_col,
            allele_col=allele_col,
            pval_col=pval_col,
            pval_threshold=1.0,  # Use max p-value to include all SNPs
            source_build=gwas_build,
            target_build=genotype_build
        )
        
        extract_file = None
        if use_clumping:
            # Perform clumping if requested
            clump_prefix = f"{base_prefix}_clumped"
            clumped_file = self.plink.clump(
                bfile=bfile,
                clump_file=pval_file,
                output_prefix=clump_prefix,
                p1=clump_p1,
                r2=clump_r2,
                kb=clump_kb,
                source_build=gwas_build,
                target_build=genotype_build
            )
            extract_file = self.data_prep.extract_clumped_snps(clumped_file)
        
        # Calculate multiple PRS
        profile_files = self.plink.calculate_multiple_prs(
            bfile=bfile,
            score_file=score_file,
            pvalue_file=pval_file,
            range_file=range_file,
            output_prefix=f"{base_prefix}_multi",
            extract_snps=extract_file,
            keep_samples=keep_samples,
            source_build=gwas_build,
            target_build=genotype_build
        )
        
        return profile_files
