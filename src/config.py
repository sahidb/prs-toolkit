import os
import configparser
from typing import List, Dict

class Config:
    """Configuration manager for the PRS pipeline."""
    
    def __init__(self, config_file=None):
        """Initialize configuration with default or custom settings."""
        self.config = configparser.ConfigParser()
        self._set_defaults()
        
        if config_file and os.path.exists(config_file):
            self.load_config(config_file)
    
    def _set_defaults(self):
        """Set default configuration values."""
        self.config["PATHS"] = {
            "plink_path": "plink",
            "output_dir": "./output",
            "temp_dir": "./temp",
            "gwas_dir": "./data/gwas",
            "genotype_dir": "./data/genotypes",
            "population_dir": "./data/population",
            "chain_dir": "./data/chain_files"
        }
        
        self.config["PARAMETERS"] = {
            "p_value_threshold": "5e-8",
            "clump_p1": "1e-5",
            "clump_r2": "0.1",
            "clump_kb": "250",
            "score_no_mean_imputation": "False",
            "clean_temp": "True"
        }
        
        self.config["POPULATIONS"] = {
            "target_populations": "EAS,SAS"
        }
        
        self.config["GWAS_COLUMNS"] = {
            "snp_col": "SNP",
            "effect_col": "BETA",
            "allele_col": "A1",
            "pval_col": "P",
            "chr_col": "CHR",
            "pos_col": "BP"
        }
        
        self.config["GENOME_BUILD"] = {
            "default_source_build": "GRCh37",
            "default_target_build": "GRCh37",
            "auto_detect": "True",
            "validate_build": "True"
        }
    
    def load_config(self, config_file):
        """Load configuration from file."""
        self.config.read(config_file)
    
    def get_path(self, key):
        """Get a path from configuration."""
        path = self.config["PATHS"].get(key, "")
        # Create directory if it doesn't exist (except for executable paths)
        if key not in ["plink_path"] and path and not os.path.exists(path):
            os.makedirs(path, exist_ok=True)
        return path
    
    def get_parameter(self, key, as_type=str):
        """Get a parameter with type conversion."""
        value = self.config["PARAMETERS"].get(key, "")
        if as_type == bool:
            return value.lower() in ("yes", "true", "t", "1")
        return as_type(value)
    
    def get_populations(self):
        """Get list of target populations."""
        pops_str = self.config["POPULATIONS"].get("target_populations", "")
        return [p.strip() for p in pops_str.split(",") if p.strip()]
    
    def get_gwas_column(self, key):
        """Get column name from GWAS_COLUMNS section."""
        return self.config["GWAS_COLUMNS"].get(key, "")
    
    def get_build_setting(self, key):
        """Get setting from GENOME_BUILD section."""
        if key == "auto_detect" or key == "validate_build":
            value = self.config["GENOME_BUILD"].get(key, "True")
            return value.lower() in ("yes", "true", "t", "1")
        return self.config["GENOME_BUILD"].get(key, "")
