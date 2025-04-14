import os
import subprocess
import pandas as pd
import numpy as np
import gzip
import tempfile
from typing import List, Dict, Optional, Union, Tuple

class GenomeBuildHandler:
    """
    Class for handling genomic coordinates across different genome builds.
    Supports conversion between GRCh37 (hg19) and GRCh38 (hg38).
    """
    
    # Define mappings between UCSC and Ensembl chromosome naming
    CHR_MAPPING = {
        # UCSC -> Ensembl
        'chr1': '1', 'chr2': '2', 'chr3': '3', 'chr4': '4', 'chr5': '5',
        'chr6': '6', 'chr7': '7', 'chr8': '8', 'chr9': '9', 'chr10': '10',
        'chr11': '11', 'chr12': '12', 'chr13': '13', 'chr14': '14', 'chr15': '15',
        'chr16': '16', 'chr17': '17', 'chr18': '18', 'chr19': '19', 'chr20': '20',
        'chr21': '21', 'chr22': '22', 'chrX': 'X', 'chrY': 'Y', 'chrM': 'MT',
        # Ensembl -> UCSC
        '1': 'chr1', '2': 'chr2', '3': 'chr3', '4': 'chr4', '5': 'chr5',
        '6': 'chr6', '7': 'chr7', '8': 'chr8', '9': 'chr9', '10': 'chr10',
        '11': 'chr11', '12': 'chr12', '13': 'chr13', '14': 'chr14', '15': 'chr15',
        '16': 'chr16', '17': 'chr17', '18': 'chr18', '19': 'chr19', '20': 'chr20',
        '21': 'chr21', '22': 'chr22', 'X': 'chrX', 'Y': 'chrY', 'MT': 'chrM'
    }
    
    # Build-specific SNPs to aid in detection
    BUILD_SPECIFIC_MARKERS = {
        'GRCh37': [
            {'chr': '6', 'pos': 31033797, 'ref': 'T', 'alt': 'G'},
            {'chr': '12', 'pos': 112843254, 'ref': 'C', 'alt': 'G'},
            {'chr': '17', 'pos': 7531642, 'ref': 'G', 'alt': 'T'}
        ],
        'GRCh38': [
            {'chr': '6', 'pos': 31353872, 'ref': 'T', 'alt': 'G'},
            {'chr': '12', 'pos': 111351574, 'ref': 'C', 'alt': 'G'},
            {'chr': '17', 'pos': 7668402, 'ref': 'G', 'alt': 'T'}
        ]
    }
    
    def __init__(self, chain_dir="./data/chain_files"):
        """
        Initialize GenomeBuildHandler with directory containing chain files.
        
        Parameters:
        -----------
        chain_dir : str
            Directory containing chain files for liftOver operations
        """
        self.chain_dir = chain_dir
        os.makedirs(chain_dir, exist_ok=True)
        
        self.chain_files = {
            '37to38': os.path.join(chain_dir, "hg19ToHg38.over.chain.gz"),
            '38to37': os.path.join(chain_dir, "hg38ToHg19.over.chain.gz")
        }
        
        self._check_chain_files()
        self._check_crossmap()
    
    def _check_chain_files(self):
        """Check if necessary chain files exist and download if missing."""
        for purpose, file_path in self.chain_files.items():
            if not os.path.exists(file_path):
                print(f"Chain file not found: {file_path}")
                self._download_chain_file(purpose, file_path)
    
    def _download_chain_file(self, purpose, file_path):
        """Download chain file if missing."""
        try:
            import urllib.request
            
            # Define URLs for chain files
            chain_urls = {
                '37to38': "http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz",
                '38to37': "http://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz"
            }
            
            url = chain_urls.get(purpose)
            if not url:
                raise ValueError(f"No URL defined for chain file purpose: {purpose}")
                
            print(f"Downloading chain file from {url}...")
            urllib.request.urlretrieve(url, file_path)
            print(f"Downloaded chain file to {file_path}")
        except Exception as e:
            print(f"Error downloading chain file: {str(e)}")
            print("Please download the chain files manually and place them in the chain_dir.")
    
    def _check_crossmap(self):
        """Verify that CrossMap is installed."""
        try:
            result = subprocess.run(
                ["CrossMap.py", "--help"],
                capture_output=True,
                text=True
            )
            return True
        except FileNotFoundError:
            print("Warning: CrossMap.py not found in PATH. Build conversion functionality will be limited.")
            print("Please install CrossMap: pip install CrossMap")
            return False
    
    def detect_genome_build(self, variant_file, 
                           chr_col="CHR", pos_col="BP", 
                           ref_col=None, alt_col=None):
        """
        Detect the genome build of a variant file.
        
        Parameters:
        -----------
        variant_file : str
            Path to the variant file
        chr_col, pos_col, ref_col, alt_col : str
            Column names for chromosome, position, reference and alternate alleles
            
        Returns:
        --------
        str: 'GRCh37', 'GRCh38', or 'Unknown'
        """
        # Load the variant data
        try:
            # First determine if the file is PLINK format or CSV
            if variant_file.endswith('.bim'):
                # Handle PLINK BIM file format (special handling)
                variants = pd.read_csv(variant_file, sep='\t', header=None)
                chr_col, pos_col = 0, 3  # BIM format has chromosome in column 0, position in column 3
            else:
                # Try to auto-detect separator for standard files
                variants = pd.read_csv(variant_file, sep=None, engine='python')
        except Exception as e:
            print(f"Error reading variant file: {str(e)}")
            return "Unknown"
        
        # Handle PLINK BIM file format - adjust column indices to names
        if variant_file.endswith('.bim'):
            # Create temporary column names for reference
            if isinstance(chr_col, int) and isinstance(pos_col, int):
                chr_col_data, pos_col_data = chr_col, pos_col
                # Convert to column names for consistent code below
                variants.columns = ["CHR", "SNP", "CM", "BP", "A1", "A2"]
                chr_col, pos_col = "CHR", "BP"
                ref_col, alt_col = "A1", "A2"
        
        # Check format of chromosome column to determine if format conversion is needed
        if variants[chr_col].iloc[0].startswith('chr'):
            # Convert UCSC format to Ensembl format for checking
            variants[chr_col] = variants[chr_col].map(
                lambda x: self.CHR_MAPPING.get(x, x) if x.startswith('chr') else x
            )
        
        # Check for build-specific markers
        confidence_scores = {'GRCh37': 0, 'GRCh38': 0}
        
        for build, markers in self.BUILD_SPECIFIC_MARKERS.items():
            for marker in markers:
                # Check if marker exists in variant file
                marker_variants = variants[
                    (variants[chr_col] == marker['chr']) & 
                    (variants[pos_col] == marker['pos'])
                ]
                
                if not marker_variants.empty:
                    # Check reference and alternate alleles if available
                    if ref_col in variants.columns and alt_col in variants.columns:
                        matching_variants = marker_variants[
                            (marker_variants[ref_col] == marker['ref']) &
                            (marker_variants[alt_col] == marker['alt'])
                        ]
                        if not matching_variants.empty:
                            confidence_scores[build] += 2
                    else:
                        confidence_scores[build] += 1
        
        # Check for alt contigs which are specific to GRCh38
        if any('_alt' in str(chrom) for chrom in variants[chr_col].unique()):
            confidence_scores['GRCh38'] += 3
        
        # Make a determination based on confidence scores
        if confidence_scores['GRCh37'] > confidence_scores['GRCh38']:
            return 'GRCh37'
        elif confidence_scores['GRCh38'] > confidence_scores['GRCh37']:
            return 'GRCh38'
        else:
            # If ties or no markers found, make a best guess based on position patterns
            # GRCh38 has different coordinates for certain chromosomes
            
            # Check some known regions with large coordinate differences
            if variants[chr_col].isin(['6', '12', '17']).any():
                chr6_variants = variants[variants[chr_col] == '6']
                chr12_variants = variants[variants[chr_col] == '12']
                chr17_variants = variants[variants[chr_col] == '17']
                
                # Check position patterns in MHC region on chromosome 6
                if not chr6_variants.empty:
                    mean_pos = chr6_variants[(chr6_variants[pos_col] > 28000000) & 
                                          (chr6_variants[pos_col] < 34000000)][pos_col].mean()
                    if not np.isnan(mean_pos):
                        if mean_pos > 32000000:  # GRCh38 tends to have higher coords in MHC
                            return 'GRCh38'
                
                # Check chr12 patterns - higher positions in GRCh37 for q arm
                if not chr12_variants.empty:
                    high_pos_variants = chr12_variants[chr12_variants[pos_col] > 110000000]
                    if not high_pos_variants.empty:
                        if high_pos_variants[pos_col].mean() > 112000000:
                            return 'GRCh37'
                
                # Check chr17 p arm around TP53
                if not chr17_variants.empty:
                    tp53_region = chr17_variants[(chr17_variants[pos_col] > 7500000) & 
                                             (chr17_variants[pos_col] < 7800000)]
                    if not tp53_region.empty:
                        if tp53_region[pos_col].mean() < 7600000:
                            return 'GRCh37'
                        else:
                            return 'GRCh38'
            
            # Default to GRCh37 if we can't determine (since it's more common in older datasets)
            return 'GRCh37'
    
    def convert_variant_build(self, input_file, output_file, 
                             source_build=None, target_build=None,
                             chr_col="CHR", pos_col="BP"):
        """
        Convert variants from one genome build to another.
        
        Parameters:
        -----------
        input_file : str
            Path to input variant file
        output_file : str
            Path for output converted file
        source_build : str
            Source genome build ('GRCh37' or 'GRCh38'), auto-detected if None
        target_build : str
            Target genome build ('GRCh37' or 'GRCh38')
        chr_col, pos_col : str
            Column names for chromosome and position
            
        Returns:
        --------
        str: Path to the converted output file
        """
        if not os.path.exists(input_file):
            raise FileNotFoundError(f"Input file not found: {input_file}")
        
        # Auto-detect source build if not specified
        if source_build is None:
            source_build = self.detect_genome_build(input_file, chr_col, pos_col)
            print(f"Detected source build: {source_build}")
        
        # Validate builds
        if source_build not in ['GRCh37', 'GRCh38']:
            raise ValueError(f"Invalid source build: {source_build}. Must be 'GRCh37' or 'GRCh38'.")
        
        if target_build not in ['GRCh37', 'GRCh38']:
            raise ValueError(f"Invalid target build: {target_build}. Must be 'GRCh37' or 'GRCh38'.")
        
        # If builds are the same, just copy the file
        if source_build == target_build:
            print(f"Source and target builds are the same ({source_build}), no conversion needed.")
            import shutil
            shutil.copy2(input_file, output_file)
            return output_file
        
        # Determine conversion direction and chain file
        conversion_key = '37to38' if source_build == 'GRCh37' else '38to37'
        chain_file = self.chain_files[conversion_key]
        
        # Check if CrossMap is available
        have_crossmap = self._check_crossmap()
        if not have_crossmap:
            print("CrossMap not available. Using built-in coordinate conversion...")
            return self._convert_without_crossmap(input_file, output_file, source_build, target_build, chr_col, pos_col)
        
        # Load variant data
        variants = pd.read_csv(input_file, sep=None, engine='python')
        
        # Format for CrossMap (requires BED format)
        # Create a temporary BED file with columns: chr, start, end, name
        with tempfile.NamedTemporaryFile(mode='w+t', suffix='.bed', delete=False) as temp_bed:
            temp_bed_path = temp_bed.name
            # Convert chromosome format if needed (CrossMap expects 'chr' prefix)
            for i, row in variants.iterrows():
                chrom = row[chr_col]
                pos = row[pos_col]
                
                # Ensure chromosome has 'chr' prefix
                if not str(chrom).startswith('chr'):
                    chrom = self.CHR_MAPPING.get(str(chrom), f"chr{chrom}")
                
                # BED format: 0-based start, 1-based end
                temp_bed.write(f"{chrom}\t{int(pos)-1}\t{int(pos)}\t{i}\n")
        
        # Create temporary output file for CrossMap
        with tempfile.NamedTemporaryFile(mode='w+t', suffix='.bed', delete=False) as temp_output:
            temp_output_path = temp_output.name
        
        # Run CrossMap
        try:
            cmd = [
                "CrossMap.py", "bed", chain_file, temp_bed_path, temp_output_path
            ]
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            
            # Load CrossMap results
            crossmap_results = {}
            with open(temp_output_path, 'r') as f:
                for line in f:
                    if line.strip():
                        parts = line.strip().split('\t')
                        # CrossMap output: chr, start, end, name
                        if len(parts) >= 4:
                            # Extract the variant index and new position
                            # Remember BED is 0-based, so convert back to 1-based
                            variant_idx = int(parts[3])
                            new_chrom = parts[0]
                            new_pos = int(parts[1]) + 1  # Convert 0-based to 1-based
                            
                            # Convert chromosome back to original format if needed
                            if chr_col in variants and not str(variants[chr_col].iloc[0]).startswith('chr'):
                                new_chrom = self.CHR_MAPPING.get(new_chrom, new_chrom.replace('chr', ''))
                            
                            crossmap_results[variant_idx] = {'chrom': new_chrom, 'pos': new_pos}
            
            # Update original dataframe with new coordinates
            updated_variants = variants.copy()
            for idx, mapping in crossmap_results.items():
                updated_variants.loc[idx, chr_col] = mapping['chrom']
                updated_variants.loc[idx, pos_col] = mapping['pos']
            
            # Save updated variants
            updated_variants.to_csv(output_file, index=False)
            
            # Log conversion stats
            total = len(variants)
            converted = len(crossmap_results)
            print(f"Conversion complete: {converted}/{total} variants successfully converted ({converted/total:.1%})")
            
            # Clean up temporary files
            os.unlink(temp_bed_path)
            os.unlink(temp_output_path)
            
            return output_file
            
        except subprocess.CalledProcessError as e:
            print(f"Error running CrossMap: {e.stderr}")
            # Clean up temporary files
            os.unlink(temp_bed_path)
            if os.path.exists(temp_output_path):
                os.unlink(temp_output_path)
            raise
    
    def _convert_without_crossmap(self, input_file, output_file, source_build, target_build, chr_col, pos_col):
        """Simplified coordinate conversion for common regions when CrossMap is not available."""
        # Load variant data
        variants = pd.read_csv(input_file, sep=None, engine='python')
        
        # Conversion offsets for key chromosomal regions
        # These are simplified approximations for common regions
        conversion_offsets = {
            # Direction: GRCh37 -> GRCh38
            'GRCh37_to_GRCh38': {
                '1': 150000,  # Chr1 positions generally increased
                '6': 320000,  # MHC region
                '12': -1490000,  # Chr12q (subtract)
                '17': 136760,  # Around TP53
                '19': 280000,  # Chr19 (gene-dense)
            },
            # Direction: GRCh38 -> GRCh37
            'GRCh38_to_GRCh37': {
                '1': -150000,  # Reverse of above
                '6': -320000,
                '12': 1490000,
                '17': -136760,
                '19': -280000,
            }
        }
        
        # Determine conversion direction
        direction = 'GRCh37_to_GRCh38' if source_build == 'GRCh37' else 'GRCh38_to_GRCh37'
        offsets = conversion_offsets[direction]
        
        # Create a copy of the dataframe for conversion
        converted = variants.copy()
        
        # Apply offsets based on chromosome
        for chrom, offset in offsets.items():
            mask = (converted[chr_col] == chrom)
            converted.loc[mask, pos_col] = converted.loc[mask, pos_col] + offset
        
        # Save converted file
        converted.to_csv(output_file, index=False)
        
        print(f"Completed simplified conversion from {source_build} to {target_build}.")
        print("Note: This is an approximate conversion. For precise results, install CrossMap.")
        
        return output_file
    
    def convert_plink_build(self, input_prefix, output_prefix, 
                           source_build=None, target_build=None):
        """
        Convert PLINK binary files from one genome build to another.
        
        Parameters:
        -----------
        input_prefix : str
            Prefix for input PLINK files (.bed, .bim, .fam)
        output_prefix : str
            Prefix for output PLINK files
        source_build : str
            Source genome build ('GRCh37' or 'GRCh38'), auto-detected if None
        target_build : str
            Target genome build ('GRCh37' or 'GRCh38')
            
        Returns:
        --------
        str: Prefix for the converted PLINK files
        """
        # Check input files
        bim_file = f"{input_prefix}.bim"
        if not os.path.exists(bim_file):
            raise FileNotFoundError(f"BIM file not found: {bim_file}")
        
        # Auto-detect source build if not specified
        if source_build is None:
            source_build = self.detect_genome_build(bim_file, chr_col=0, pos_col=3)
            print(f"Detected source build: {source_build}")
        
        # Validate builds
        if source_build not in ['GRCh37', 'GRCh38']:
            raise ValueError(f"Invalid source build: {source_build}. Must be 'GRCh37' or 'GRCh38'.")
        
        if target_build not in ['GRCh37', 'GRCh38']:
            raise ValueError(f"Invalid target build: {target_build}. Must be 'GRCh37' or 'GRCh38'.")
        
        # If builds are the same, just copy the files
        if source_build == target_build:
            print(f"Source and target builds are the same ({source_build}), no conversion needed.")
            for ext in ['.bed', '.bim', '.fam']:
                import shutil
                source = f"{input_prefix}{ext}"
                dest = f"{output_prefix}{ext}"
                if os.path.exists(source):
                    shutil.copy2(source, dest)
            return output_prefix
        
        # First, convert the BIM file coordinates
        temp_output_bim = f"{output_prefix}_temp.bim"
        self.convert_variant_build(
            bim_file, temp_output_bim, source_build, target_build, 
            chr_col=0, pos_col=3
        )
        
        # Now use PLINK to update the binary files
        try:
            import subprocess
            cmd = [
                "plink", "--bfile", input_prefix,
                "--update-map", temp_output_bim, "--update-chr", temp_output_bim, "1", "2",
                "--make-bed", "--out", output_prefix
            ]
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            
            # Clean up temporary file
            os.unlink(temp_output_bim)
            
            print(f"Successfully converted PLINK files from {source_build} to {target_build}")
            return output_prefix
            
        except subprocess.CalledProcessError as e:
            print(f"Error updating PLINK files: {e.stderr}")
            if os.path.exists(temp_output_bim):
                os.unlink(temp_output_bim)
            raise
