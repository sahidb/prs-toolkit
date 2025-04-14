#!/usr/bin/env python3
"""
Genome build conversion script.

This script converts genomic coordinates between GRCh37 (hg19) and GRCh38 (hg38).

Example usage:
python convert_build.py --input data/gwas/cad_gwas.csv --output data/gwas/cad_gwas_hg38.csv 
                        --source GRCh37 --target GRCh38
"""

import os
import sys
import argparse
import logging

# Add parent directory to path to import modules
# sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from src.genome_build import GenomeBuildHandler
from src.config import Config
from src.utils import Utils

def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Convert genomic coordinates between builds.")
    
    parser.add_argument("--config", help="Path to configuration file")
    parser.add_argument("--input", required=True, help="Path to input file")
    parser.add_argument("--output", required=True, help="Path for output file")
    parser.add_argument("--source", choices=["GRCh37", "GRCh38", "auto"], default="auto",
                      help="Source genome build (default: auto-detect)")
    parser.add_argument("--target", choices=["GRCh37", "GRCh38"], required=True,
                      help="Target genome build")
    parser.add_argument("--file_type", choices=["gwas", "plink"], default="gwas",
                      help="File type to convert (gwas or plink)")
    parser.add_argument("--chr_col", default="CHR", help="Chromosome column name (for GWAS files)")
    parser.add_argument("--pos_col", default="BP", help="Position column name (for GWAS files)")
    
    return parser.parse_args()

def main():
    """Main function to run genome build conversion."""
    # Parse arguments
    args = parse_args()
    
    # Load configuration
    config = Config(args.config)
    
    # Set up logging
    logger = Utils.setup_logging(os.path.join(config.get_path("output_dir"), "build_conversion.log"))
    logger.info("Starting genome build conversion")
    
    try:
        # Initialize genome build handler
        chain_dir = config.get_path("chain_dir")
        genome_build_handler = GenomeBuildHandler(chain_dir)
        
        # Determine source build if set to auto
        source_build = args.source
        if source_build == "auto":
            logger.info("Auto-detecting source genome build")
            if args.file_type == "gwas":
                source_build = genome_build_handler.detect_genome_build(
                    args.input, chr_col=args.chr_col, pos_col=args.pos_col
                )
            else:  # plink
                bim_file = f"{args.input}.bim"
                if os.path.exists(bim_file):
                    source_build = genome_build_handler.detect_genome_build(
                        bim_file, chr_col=0, pos_col=3
                    )
                else:
                    logger.error(f"BIM file not found: {bim_file}")
                    return 1
            
            logger.info(f"Detected source build: {source_build}")
        
        # Validate source build
        if source_build not in ["GRCh37", "GRCh38"]:
            logger.error(f"Invalid or undetectable source build: {source_build}")
            return 1
        
        # Check if conversion is necessary
        if source_build == args.target:
            logger.info(f"Source and target builds are the same ({source_build}), no conversion needed.")
            return 0
        
        # Convert coordinates
        logger.info(f"Converting from {source_build} to {args.target}")
        
        if args.file_type == "gwas":
            output_file = genome_build_handler.convert_variant_build(
                args.input, args.output, 
                source_build=source_build, 
                target_build=args.target,
                chr_col=args.chr_col, 
                pos_col=args.pos_col
            )
        else:  # plink
            output_prefix = genome_build_handler.convert_plink_build(
                args.input, args.output, 
                source_build=source_build, 
                target_build=args.target
            )
            output_file = f"{output_prefix}.bim"
        
        logger.info(f"Conversion completed successfully. Output saved to: {args.output}")
        return 0
        
    except Exception as e:
        logger.error(f"Error in genome build conversion: {str(e)}", exc_info=True)
        return 1

if __name__ == "__main__":
    sys.exit(main())
