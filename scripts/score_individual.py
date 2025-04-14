#!/usr/bin/env python3
"""
Individual PRS scoring script.

This script calculates PRS for an individual and compares against EAS and SAS population distributions.

Example usage:
python score_individual.py --config config/myconfig.ini --gwas data/gwas/cad_gwas.csv 
                          --individual data/genotypes/individual_sample --eas_pop data/population/eas_scores.profile
                          --sas_pop data/population/sas_scores.profile --gwas_build GRCh37 --genotype_build GRCh38
                          --out individual_prs
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
from src.individual_scorer import IndividualScorer
from src.utils import Utils
from src.genome_build import GenomeBuildHandler

def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Calculate and compare individual PRS to EAS/SAS populations.")
    
    parser.add_argument("--config", help="Path to configuration file")
    parser.add_argument("--gwas", required=True, help="Path to GWAS summary statistics file")
    parser.add_argument("--individual", required=True, help="Prefix for individual's PLINK files")
    parser.add_argument("--p_threshold", type=float, default=5e-8, help="P-value threshold for SNP selection")
    parser.add_argument("--use_clumping", action="store_true", help="Perform clumping to select independent SNPs")
    
    # Population references
    parser.add_argument("--eas_pop", required=True, help="Path to EAS population PRS profile file")
    parser.add_argument("--sas_pop", required=True, help="Path to SAS population PRS profile file")
    parser.add_argument("--standardize", action="store_true", help="Standardize scores for comparison")
    
    # Genome build options
    parser.add_argument("--gwas_build", choices=["GRCh37", "GRCh38", "auto"],
                      help="Genome build of GWAS data (default: auto-detect)")
    parser.add_argument("--genotype_build", choices=["GRCh37", "GRCh38", "auto"],
                      help="Genome build of genotype data (default: auto-detect)")
    
    parser.add_argument("--out", required=True, help="Output name prefix")
    
    return parser.parse_args()

def main():
    """Main function to run individual PRS scoring."""
    # Parse arguments
    args = parse_args()
    
    # Load configuration
    config = Config(args.config)
    
    # Set up logging
    logger = Utils.setup_logging(os.path.join(config.get_path("output_dir"), "individual_scoring.log"))
    logger.info("Starting individual PRS scoring")
    
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
                bim_file = f"{args.individual}.bim"
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
        
        # Create individual scorer
        scorer = IndividualScorer(calculator, plink, config.get_path("output_dir"))
        
        # Calculate individual PRS
        logger.info(f"Calculating PRS for individual sample")
        logger.info(f"GWAS build: {gwas_build}, Genotype build: {genotype_build}")
        
        individual_score, profile_file = scorer.calculate_individual_prs(
            gwas_file=args.gwas,
            individual_file=args.individual,
            output_name=args.out,
            pval_threshold=args.p_threshold,
            use_clumping=args.use_clumping,
            source_build=gwas_build,
            target_build=genotype_build
        )
        
        logger.info(f"Individual PRS score: {individual_score}")
        
        # Load population distributions
        logger.info("Loading population distributions")
        population_files = {
            "EAS": args.eas_pop,
            "SAS": args.sas_pop
        }
        
        scorer.load_population_distributions(population_files)
        
        # Compare to populations
        logger.info("Comparing to population distributions")
        comparison_results = scorer.compare_to_populations(
            individual_score=individual_score,
            standardize=args.standardize
        )
        
        # Generate visualization
        logger.info("Generating visualization")
        viz_file = scorer.visualize_comparison(
            individual_score=individual_score,
            output_file=os.path.join(config.get_path("output_dir"), f"{args.out}_comparison.png"),
            standardize=args.standardize
        )
        
        # Generate report
        logger.info("Generating report")
        report_file = scorer.generate_report(
            individual_score=individual_score,
            comparison_results=comparison_results,
            output_file=os.path.join(config.get_path("output_dir"), f"{args.out}_report.html")
        )
        
        logger.info(f"Individual scoring complete. Results saved to:")
        logger.info(f"  - PRS profile: {profile_file}")
        logger.info(f"  - Visualization: {viz_file}")
        logger.info(f"  - Report: {report_file}")
        
        # Log the assessment against each population
        for pop, results in comparison_results.items():
            logger.info(f"  - {pop} assessment: {results['percentile']:.1f}th percentile, {results['risk_category']} risk")
        
        # Optional: clean temporary files
        if config.get_parameter("clean_temp", as_type=bool):
            Utils.clean_temp_files(config.get_path("temp_dir"))
        
        return 0
        
    except Exception as e:
        logger.error(f"Error in individual PRS scoring: {str(e)}", exc_info=True)
        return 1

if __name__ == "__main__":
    sys.exit(main())
