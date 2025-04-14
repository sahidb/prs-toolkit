import os
import logging
import shutil
import pandas as pd
from typing import List, Dict, Any, Optional

class Utils:
    """Utility functions for the PRS pipeline."""
    
    @staticmethod
    def setup_logging(log_file=None, level=logging.INFO):
        """Set up logging configuration."""
        # Create directory for log file if needed
        if log_file and not os.path.exists(os.path.dirname(log_file)):
            os.makedirs(os.path.dirname(log_file), exist_ok=True)
            
        # Configure basic logging format
        logging.basicConfig(
            level=level,
            format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )
        
        # Create logger
        logger = logging.getLogger('prs_toolkit')
        
        # Add file handler if log file specified
        if log_file:
            file_handler = logging.FileHandler(log_file)
            file_formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
            file_handler.setFormatter(file_formatter)
            logger.addHandler(file_handler)
        
        return logger
    
    @staticmethod
    def get_plink_version(plink_path):
        """Get the version of the PLINK executable."""
        try:
            import subprocess
            result = subprocess.run(
                [plink_path, "--version"],
                capture_output=True,
                text=True
            )
            # Extract version number from output
            lines = result.stdout.strip().split('\n')
            if lines:
                return lines[0]
            return "Unknown"
        except Exception:
            return "Error detecting PLINK version"
    
    @staticmethod
    def clean_temp_files(temp_dir, keep_log=True):
        """Clean temporary files but optionally preserve logs."""
        if not os.path.exists(temp_dir):
            return
            
        for item in os.listdir(temp_dir):
            item_path = os.path.join(temp_dir, item)
            # Keep log files if requested
            if keep_log and item.endswith('.log'):
                continue
            
            if os.path.isfile(item_path):
                os.unlink(item_path)
            elif os.path.isdir(item_path):
                shutil.rmtree(item_path)
    
    @staticmethod
    def validate_gwas_file(gwas_file, required_cols=None):
        """Validate that a GWAS file contains required columns."""
        if required_cols is None:
            required_cols = ['SNP', 'A1', 'BETA', 'P']
        
        try:
            # Try to read the first few rows
            df = pd.read_csv(gwas_file, sep=None, engine='python', nrows=5)
            
            # Check for required columns
            missing_cols = [col for col in required_cols if col not in df.columns]
            
            if missing_cols:
                return False, f"Missing required columns: {', '.join(missing_cols)}"
            
            return True, "GWAS file is valid"
        except Exception as e:
            return False, f"Error validating GWAS file: {str(e)}"
