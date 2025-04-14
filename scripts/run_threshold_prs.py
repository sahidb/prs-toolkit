#!/usr/bin/env python3
"""
PRS calculation at multiple P-value thresholds.

Example usage:
python run_threshold_prs.py --config config/myconfig.ini --gwas data/gwas/cad_gwas.csv 
                            --bfile data/genotypes/1000g --thresholds 5e-8,1e-6,1e-4,0.001,0.01,0.05,0.1,0.5,1.0
                            --use_clumping --gwas_build GRCh37 --genotype_build GRCh38
                            --out cad_threshold_prs
"""

import os
import sys
import argparse
import logging

# Add parent directory to path to import modules
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from src.config import Config
from src.data_preparation import DataPreparation
from src.plink_interface import PlinkInterface
from src.prs_calculator import PRSCalculator
from src.utils import Utils
from src.genome_build import GenomeBuildHandler

def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Calculate PRS at multiple P-value thresholds.")
    
    parser.add_argument("--config", help="Path to configuration file")
    parser.add_argument("--gwas", required=True, help="Path to GWAS summary statistics file")
    parser.add_argument("--bfile", required=True, help="Prefix for PLINK binary files")
    parser.add_argument("--thresholds", help="Comma-separated list of P-value thresholds")
    parser.add_argument("--use_clumping", action="store_true", help="Perform clumping to select independent SNPs")
    parser.add_argument("--out", required=True, help="Output name prefix")
    
    # Population option
    parser.add_argument("--keep", help="File containing samples to keep")
    
    # Clumping parameters
    parser.add_argument("--clump_p1", type=float, default=1e-5, help="Significance threshold for index SNPs")
    parser.add_argument("--clump_r2", type=float, default=0.1, help="LD r² threshold for clumping")
    parser.add_argument("--clump_kb", type=int, default=250, help="Physical distance threshold in kb")
    
    # GWAS columns
    parser.add_argument("--snp_col", default="SNP", help="SNP column name in GWAS file")
    parser.add_argument("--effect_col", default="BETA", help="Effect size column name in GWAS file")
    parser.add_argument("--allele_col", default="A1", help="Effect allele column name in GWAS file")
    parser.add_argument("--pval_col", default="P", help="P-value column name in GWAS file")
    
    # Genome build options
    parser.add_argument("--gwas_build", choices=["GRCh37", "GRCh38", "auto"],
                      help="Genome build of GWAS data (default: auto-detect)")
    parser.add_argument("--genotype_build", choices=["GRCh37", "GRCh38", "auto"],
                      help="Genome build of genotype data (default: auto-detect)")
    
    return parser.parse_args()

def main():
    """Main function to run multi-threshold PRS calculation."""
    # Parse arguments
    args = parse_args()
    
    # Load configuration
    config = Config(args.config)
    
    # Set up logging
    logger = Utils.setup_logging(os.path.join(config.get_path("output_dir"), "threshold_prs.log"))
    logger.info("Starting multi-threshold PRS calculation")
    
    try:
        # Initialize components
        chain_dir = config.get_path("chain_dir")
        data_prep = DataPreparation(config.get_path("temp_dir"), chain_dir)
        plink = PlinkInterface(config.get_path("plink_path"), chain_dir)
        calculator = PRSCalculator(data_prep, plink, config.get_path("output_dir"), chain_dir)
        
        # Initialize genome build handler
        genome_build_handler = GenomeBuildHandler(chain_dir)
        
        # Determine GWAS build
        gwas_build = args.gwas_build
        if gwas_build == "auto" or gwas_build is None:
            if config.get_build_setting("auto_detect"):
                logger.info("Auto-detecting GWAS genome build")
                gwas_build = genome_build_handler.detect_genome_build(
                    args.gwas, chr_col="CHR", pos_col="BP"
                )
                logger.info(f"Detected GWAS build: {gwas_build}")
            else:
                gwas_build = config.get_build_setting("default_source_build")
                logger.info(f"Using default GWAS build from config: {gwas_build}")
        
        # Determine genotype build
        genotype_build = args.genotype_build
        if genotype_build == "auto" or genotype_build is None:
            if config.get_build_setting("auto_detect"):
                logger.info("Auto-detecting genotype data genome build")
                bim_file = f"{args.bfile}.bim"
                if os.path.exists(bim_file):
                    genotype_build = genome_build_handler.detect_genome_build(
                        bim_file, chr_col=0, pos_col=3
                    )
                    logger.info(f"Detected genotype build: {genotype_build}")
                else:
                    logger.warning(f"Cannot auto-detect genotype build: BIM file not found")
                    genotype_build = config.get_build_setting("default_target_build")
            else:
                genotype_build = config.get_build_setting("default_target_build")
                logger.info(f"Using default genotype build from config: {genotype_build}")
        
        # Parse thresholds
        if args.thresholds:
            thresholds = [float(t) for t in args.thresholds.split(",")]
        else:
            thresholds = None  # Use default thresholds in calculator
        
        logger.info(f"Calculating PRS at multiple thresholds")
        logger.info(f"GWAS build: {gwas_build}, Genotype build: {genotype_build}")
        logger.info(f"Using clumping: {args.use_clumping}")
        
        # Calculate multi-threshold PRS
        profile_files = calculator.multi_threshold_prs(
            gwas_file=args.gwas,
            bfile=args.bfile,
            output_name=args.out,
            thresholds=thresholds,
            snp_col=args.snp_col,
            effect_col=args.effect_col,
            allele_col=args.allele_col,
            pval_col=args.pval_col,
            keep_samples=args.keep,
            use_clumping=args.use_clumping,
            clump_p1=args.clump_p1,
            clump_r2=args.clump_r2,
            clump_kb=args.clump_kb,
            source_build=gwas_build,
            target_build=genotype_build
        )
        
        # Log results
        logger.info(f"PRS calculation complete. Results:")
        for threshold, profile_file in profile_files.items():
            logger.info(f"  Threshold {threshold}: {profile_file}")
        
        # Optional: clean temporary files
        if config.get_parameter("clean_temp", as_type=bool):
            Utils.clean_temp_files(config.get_path("temp_dir"))
            
        return 0
        
    except Exception as e:
        logger.error(f"Error in multi-threshold PRS calculation: {str(e)}", exc_info=True)
        return 1

if __name__ == "__main__":
    sys.exit(main())
