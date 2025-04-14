import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from typing import Dict, List, Optional, Union
from sklearn.metrics import roc_curve, auc

class PRSAnalysis:
    """Module for analyzing PRS results and generating visualizations."""
    
    def __init__(self, output_dir="./output"):
        """Initialize with directory for output files."""
        self.output_dir = output_dir
        os.makedirs(output_dir, exist_ok=True)
    
    def load_profile(self, profile_file):
        """Load a PLINK profile file containing PRS results."""
        return pd.read_csv(profile_file, delim_whitespace=True)
    
    def load_multiple_profiles(self, profile_files):
        """Load multiple PLINK profile files."""
        profiles = {}
        for label, file_path in profile_files.items():
            profiles[label] = self.load_profile(file_path)
        return profiles
    
    def merge_with_phenotype(self, prs_data, phenotype_file, id_col='IID', pheno_col='PHENO'):
        """Merge PRS results with phenotype data for evaluation."""
        # Load phenotype data
        pheno_df = pd.read_csv(phenotype_file, sep=None, engine='python')
        
        # Ensure ID column exists in both datasets
        if id_col not in prs_data.columns or id_col not in pheno_df.columns:
            raise ValueError(f"ID column '{id_col}' not found in both datasets")
        
        # Merge on the ID column
        merged = pd.merge(prs_data, pheno_df[[id_col, pheno_col]], on=id_col, how='inner')
        
        print(f"Merged PRS data with phenotypes: {len(merged)} individuals matched")
        return merged
    
    def evaluate_prs(self, merged_data, score_col='SCORE', pheno_col='PHENO', 
                    continuous=False, output_prefix=None):
        """Evaluate PRS performance against phenotype."""
        results = {}
        
        if continuous:
            # For quantitative traits: calculate variance explained (R²)
            from sklearn.linear_model import LinearRegression
            from sklearn.metrics import r2_score
            
            X = merged_data[[score_col]]
            y = merged_data[pheno_col]
            
            model = LinearRegression().fit(X, y)
            y_pred = model.predict(X)
            r2 = r2_score(y, y_pred)
            
            results['r2'] = r2
            results['beta'] = model.coef_[0]
            results['intercept'] = model.intercept_
            
            print(f"R² = {r2:.4f}, Beta = {results['beta']:.4f}")
            
            # Optional: Create scatter plot
            if output_prefix:
                plt.figure(figsize=(8, 6))
                plt.scatter(merged_data[score_col], merged_data[pheno_col], alpha=0.5)
                plt.plot(merged_data[score_col], y_pred, color='red', linewidth=2)
                plt.xlabel('Polygenic Risk Score')
                plt.ylabel('Phenotype')
                plt.title(f'PRS vs. Phenotype (R² = {r2:.4f})')
                plt.savefig(f"{output_prefix}_scatter.png", dpi=300, bbox_inches='tight')
                plt.close()
        
        else:
            # For binary traits: calculate AUC
            from sklearn.metrics import roc_curve, roc_auc_score
            
            if len(merged_data[pheno_col].unique()) != 2:
                raise ValueError("Phenotype must be binary for AUC calculation")
            
            # Calculate ROC and AUC
            fpr, tpr, _ = roc_curve(merged_data[pheno_col], merged_data[score_col])
            auc_value = auc(fpr, tpr)
            
            results['auc'] = auc_value
            results['fpr'] = fpr
            results['tpr'] = tpr
            
            print(f"AUC = {auc_value:.4f}")
            
            # Optional: Create ROC curve
            if output_prefix:
                plt.figure(figsize=(8, 6))
                plt.plot(fpr, tpr, color='blue', lw=2, 
                         label=f'ROC curve (AUC = {auc_value:.4f})')
                plt.plot([0, 1], [0, 1], color='gray', linestyle='--')
                plt.xlim([0.0, 1.0])
                plt.ylim([0.0, 1.05])
                plt.xlabel('False Positive Rate')
                plt.ylabel('True Positive Rate')
                plt.title('Receiver Operating Characteristic (ROC) Curve')
                plt.legend(loc="lower right")
                plt.savefig(f"{output_prefix}_roc.png", dpi=300, bbox_inches='tight')
                plt.close()
        
        return results
    
    def compare_thresholds(self, profiles_dict, phenotype_file, metric='auc',
                          id_col='IID', pheno_col='PHENO', continuous=False,
                          output_prefix=None):
        """Compare PRS performance across different p-value thresholds."""
        results = {}
        labels = []
        scores = []
        
        # Evaluate each threshold
        for label, profile_df in profiles_dict.items():
            # Merge with phenotype
            merged = self.merge_with_phenotype(
                profile_df, phenotype_file, id_col, pheno_col
            )
            
            # Evaluate PRS
            eval_results = self.evaluate_prs(
                merged, 'SCORE', pheno_col, continuous
            )
            
            # Store results
            score_value = eval_results['r2'] if continuous else eval_results['auc']
            results[label] = score_value
            labels.append(label)
            scores.append(score_value)
        
        # Create comparison plot
        if output_prefix:
            plt.figure(figsize=(10, 6))
            plt.bar(labels, scores, color='skyblue')
            plt.xlabel('P-value Threshold')
            plt.ylabel('R²' if continuous else 'AUC')
            plt.title('PRS Performance Across P-value Thresholds')
            plt.xticks(rotation=45)
            plt.tight_layout()
            plt.savefig(f"{output_prefix}_threshold_comparison.png", dpi=300)
            plt.close()
        
        # Identify best threshold
        best_threshold = max(results, key=results.get)
        print(f"Best threshold: {best_threshold} with {'R²' if continuous else 'AUC'} = {results[best_threshold]:.4f}")
        
        return results, best_threshold
    
    def compare_populations(self, pop_profiles, pheno_files=None, 
                           output_prefix=None, standardize=True):
        """Compare PRS distributions across different populations."""
        pop_dfs = {}
        
        # Load profile data for each population
        for pop, profile_file in pop_profiles.items():
            pop_dfs[pop] = self.load_profile(profile_file)
            
            # If phenotype file provided, merge and evaluate
            if pheno_files and pop in pheno_files:
                merged = self.merge_with_phenotype(pop_dfs[pop], pheno_files[pop])
                # Evaluation can be added here if needed
        
        # Standardize scores for comparison (if requested)
        if standardize:
            for pop in pop_dfs:
                pop_dfs[pop]['STD_SCORE'] = (pop_dfs[pop]['SCORE'] - pop_dfs[pop]['SCORE'].mean()) / pop_dfs[pop]['SCORE'].std()
        
        # Create density plot for visual comparison
        if output_prefix:
            plt.figure(figsize=(10, 6))
            
            for pop, df in pop_dfs.items():
                score_col = 'STD_SCORE' if standardize else 'SCORE'
                df[score_col].plot.density(label=pop)
            
            plt.xlabel('Polygenic Risk Score' + (' (Standardized)' if standardize else ''))
            plt.ylabel('Density')
            plt.title('PRS Distribution Across Populations')
            plt.legend()
            plt.savefig(f"{output_prefix}_population_comparison.png", dpi=300)
            plt.close()
        
        # Calculate summary statistics
        summary_stats = {}
        for pop, df in pop_dfs.items():
            score_col = 'STD_SCORE' if standardize else 'SCORE'
            summary_stats[pop] = {
                'mean': df[score_col].mean(),
                'std': df[score_col].std(),
                'min': df[score_col].min(),
                'max': df[score_col].max(),
                'n': len(df)
            }
        
        return pop_dfs, summary_stats
