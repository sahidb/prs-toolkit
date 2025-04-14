import os
import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
from typing import Dict, List, Optional, Union, Tuple

class IndividualScorer:
    """Module for scoring individual samples against population distributions."""
    
    def __init__(self, prs_calculator, plink_interface, output_dir="./output"):
        """Initialize with calculator and interface objects."""
        self.prs_calculator = prs_calculator
        self.plink_interface = plink_interface
        self.output_dir = output_dir
        os.makedirs(output_dir, exist_ok=True)
        
        # Population statistics cache
        self.pop_stats = {}
    
    def calculate_individual_prs(self, gwas_file, individual_file, output_name,
                                snp_col="SNP", effect_col="BETA", allele_col="A1", pval_col="P",
                                pval_threshold=5e-8, use_clumping=True,
                                source_build=None, target_build=None):
        """Calculate PRS for an individual sample."""
        # Prepare temporary output path
        temp_output = os.path.join(self.output_dir, f"{output_name}_temp")
        
        # Calculate PRS for the individual sample
        if use_clumping:
            profile_file = self.prs_calculator.clumped_prs(
                gwas_file=gwas_file,
                bfile=individual_file,
                output_name=temp_output,
                pval_threshold=pval_threshold,
                source_build=source_build,
                target_build=target_build
            )
        else:
            profile_file = self.prs_calculator.basic_prs(
                gwas_file=gwas_file,
                bfile=individual_file,
                output_name=temp_output,
                pval_threshold=pval_threshold,
                source_build=source_build,
                target_build=target_build
            )
        
        # Load resulting PRS
        prs_data = pd.read_csv(profile_file, delim_whitespace=True)
        
        # Extract the individual's PRS score
        if len(prs_data) == 1:
            prs_score = prs_data['SCORE'].iloc[0]
        else:
            # In case there are multiple rows, take the average
            prs_score = prs_data['SCORE'].mean()
            
        return prs_score, profile_file
    
    def load_population_distributions(self, pop_files: Dict[str, str]):
        """
        Load PRS distributions for reference populations.
        
        Parameters:
        -----------
        pop_files : Dict[str, str]
            Dictionary mapping population names to profile files
        """
        population_data = {}
        
        for pop, file_path in pop_files.items():
            # Load population PRS data
            pop_df = pd.read_csv(file_path, delim_whitespace=True)
            
            # Calculate distribution statistics
            stats = {
                'mean': pop_df['SCORE'].mean(),
                'std': pop_df['SCORE'].std(),
                'percentiles': {
                    p: np.percentile(pop_df['SCORE'], p) 
                    for p in [1, 5, 10, 20, 25, 50, 75, 80, 90, 95, 99]
                },
                'raw_scores': pop_df['SCORE'].values
            }
            
            population_data[pop] = stats
            
        self.pop_stats = population_data
        return population_data
    
    def compare_to_populations(self, individual_score: float, 
                              populations: Optional[List[str]] = None,
                              standardize: bool = True) -> Dict[str, Dict]:
        """
        Compare individual PRS to population distributions.
        
        Parameters:
        -----------
        individual_score : float
            The individual's PRS score
        populations : Optional[List[str]]
            List of population names to compare against (default: all loaded populations)
        standardize : bool
            Whether to standardize scores for comparison (recommended across populations)
            
        Returns:
        --------
        Dict mapping populations to comparison results
        """
        if not self.pop_stats:
            raise ValueError("No population distributions loaded. Call load_population_distributions first.")
        
        if populations is None:
            populations = list(self.pop_stats.keys())
            
        comparison_results = {}
        
        for pop in populations:
            if pop not in self.pop_stats:
                continue
                
            # Get population statistics
            pop_mean = self.pop_stats[pop]['mean']
            pop_std = self.pop_stats[pop]['std']
            
            # Standardize scores if requested
            if standardize:
                std_individual_score = (individual_score - pop_mean) / pop_std
                
                # Calculate percentile using standardized score
                percentile = stats.norm.cdf(std_individual_score) * 100
                
                # Determine risk category based on standardized score
                if std_individual_score <= -1.5:
                    risk_category = "Very Low"
                elif std_individual_score <= -0.5:
                    risk_category = "Low"
                elif std_individual_score <= 0.5:
                    risk_category = "Average"
                elif std_individual_score <= 1.5:
                    risk_category = "High"
                else:
                    risk_category = "Very High"
            else:
                # Calculate percentile using raw score
                percentile = stats.percentileofscore(self.pop_stats[pop]['raw_scores'], individual_score)
                std_individual_score = (individual_score - pop_mean) / pop_std
                
                # Determine risk category based on percentiles
                if percentile <= 10:
                    risk_category = "Very Low"
                elif percentile <= 30:
                    risk_category = "Low"
                elif percentile <= 70:
                    risk_category = "Average"
                elif percentile <= 90:
                    risk_category = "High"
                else:
                    risk_category = "Very High"
            
            # Store results
            comparison_results[pop] = {
                'percentile': percentile,
                'standardized_score': std_individual_score,
                'risk_category': risk_category,
                'population_mean': pop_mean,
                'population_std': pop_std
            }
            
        return comparison_results
    
    def visualize_comparison(self, individual_score: float, 
                            target_populations: List[str] = ["EAS", "SAS"],
                            output_file: Optional[str] = None,
                            standardize: bool = True) -> str:
        """
        Create visualization comparing individual PRS to target populations.
        
        Parameters:
        -----------
        individual_score : float
            The individual's PRS score
        target_populations : List[str]
            List of target populations to compare against
        output_file : Optional[str]
            Path to save visualization
        standardize : bool
            Whether to standardize scores for visualization
            
        Returns:
        --------
        Path to saved visualization
        """
        if not output_file:
            output_file = os.path.join(self.output_dir, "prs_comparison.png")
            
        # Ensure we have the target populations
        available_pops = [pop for pop in target_populations if pop in self.pop_stats]
        
        if not available_pops:
            raise ValueError(f"None of the target populations {target_populations} are available")
            
        # Set up the plot
        fig, axes = plt.subplots(len(available_pops), 1, figsize=(10, 4*len(available_pops)))
        
        # Handle case with only one population
        if len(available_pops) == 1:
            axes = [axes]
            
        for i, pop in enumerate(available_pops):
            # Extract population data
            pop_mean = self.pop_stats[pop]['mean']
            pop_std = self.pop_stats[pop]['std']
            
            # Generate distribution curve
            x = np.linspace(pop_mean - 4*pop_std, pop_mean + 4*pop_std, 1000)
            y = stats.norm.pdf(x, pop_mean, pop_std)
            
            # Plot distribution
            axes[i].plot(x, y, label=f"{pop} Distribution")
            
            # Standardize or use raw scores
            if standardize:
                # In standardized mode, we convert the individual score to population-specific z-score
                std_individual_score = (individual_score - pop_mean) / pop_std
                
                # Draw individual score as vertical line at standardized position
                x_loc = pop_mean + std_individual_score * pop_std
                
                # Get percentile for standardized score
                percentile = stats.norm.cdf(std_individual_score) * 100
            else:
                # In raw mode, we use the actual score
                x_loc = individual_score
                
                # Calculate percentile from empirical distribution
                percentile = stats.percentileofscore(self.pop_stats[pop]['raw_scores'], individual_score)
            
            # Plot individual score
            axes[i].axvline(x=x_loc, color='red', linestyle='--', 
                         label=f"Individual (Percentile: {percentile:.1f}%)")
            
            # Add percentile lines for reference
            for p, color in [(10, 'blue'), (50, 'green'), (90, 'purple')]:
                if standardize:
                    # For standardized display, convert percentile to z-score
                    z = stats.norm.ppf(p/100)
                    perc_loc = pop_mean + z * pop_std
                else:
                    # For raw display, get actual percentile value
                    perc_loc = np.percentile(self.pop_stats[pop]['raw_scores'], p)
                
                axes[i].axvline(x=perc_loc, color=color, alpha=0.5, linestyle=':',
                             label=f"{p}th Percentile")
            
            # Configure plot
            axes[i].set_title(f"PRS Comparison to {pop} Population")
            axes[i].set_xlabel("Polygenic Risk Score" + (" (Raw)" if not standardize else ""))
            axes[i].set_ylabel("Density")
            axes[i].legend()
            
        # Adjust layout and save
        plt.tight_layout()
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        return output_file
    
    def generate_report(self, individual_score: float,
                       comparison_results: Dict[str, Dict],
                       output_file: Optional[str] = None) -> str:
        """
        Generate a detailed report of PRS comparison.
        
        Parameters:
        -----------
        individual_score : float
            The individual's PRS score
        comparison_results : Dict[str, Dict]
            Results from compare_to_populations
        output_file : Optional[str]
            Path to save the report
            
        Returns:
        --------
        Path to saved report
        """
        if not output_file:
            output_file = os.path.join(self.output_dir, "prs_report.html")
            
        # Generate HTML report
        html_content = f"""
        <!DOCTYPE html>
        <html>
        <head>
            <title>Polygenic Risk Score (PRS) Report</title>
            <style>
                body {{ font-family: Arial, sans-serif; margin: 40px; }}
                .header {{ text-align: center; margin-bottom: 30px; }}
                .section {{ margin-bottom: 20px; }}
                table {{ border-collapse: collapse; width: 100%; }}
                th, td {{ padding: 8px; border: 1px solid #ddd; text-align: left; }}
                th {{ background-color: #f2f2f2; }}
                tr:nth-child(even) {{ background-color: #f9f9f9; }}
                .very-high {{ background-color: #ff8c8c; }}
                .high {{ background-color: #ffcc80; }}
                .average {{ background-color: #ffffb3; }}
                .low {{ background-color: #b3ffb3; }}
                .very-low {{ background-color: #b3d9ff; }}
            </style>
        </head>
        <body>
            <div class="header">
                <h1>Polygenic Risk Score (PRS) Report</h1>
                <p>Generated on {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
            </div>
            
            <div class="section">
                <h2>Individual PRS Information</h2>
                <p>Raw PRS Score: {individual_score:.4f}</p>
            </div>
            
            <div class="section">
                <h2>Population Comparisons</h2>
                <table>
                    <tr>
                        <th>Population</th>
                        <th>Percentile</th>
                        <th>Standardized Score</th>
                        <th>Risk Category</th>
                    </tr>
        """
        
        # Add rows for each population
        for pop, results in comparison_results.items():
            risk_class = results['risk_category'].lower().replace(' ', '-')
            html_content += f"""
                    <tr>
                        <td>{pop}</td>
                        <td>{results['percentile']:.2f}%</td>
                        <td>{results['standardized_score']:.2f}</td>
                        <td class="{risk_class}">{results['risk_category']}</td>
                    </tr>
            """
            
        html_content += """
                </table>
            </div>
            
            <div class="section">
                <h2>Interpretation</h2>
                <p>
                    This report compares an individual's polygenic risk score (PRS) against 
                    reference distributions from different populations. The PRS represents 
                    the cumulative effect of many genetic variants associated with a specific trait or disease.
                </p>
                <p>
                    <strong>Percentile:</strong> Indicates where the individual's score falls within the 
                    population distribution. A percentile of 90% means the score is higher than 90% of people 
                    in that population.
                </p>
                <p>
                    <strong>Standardized Score:</strong> Represents how many standard deviations the individual's 
                    score is from the population mean. Positive values are above average, negative are below average.
                </p>
                <p>
                    <strong>Risk Categories:</strong>
                    <ul>
                        <li><strong>Very High:</strong> Score is in the top 10% of the population (>90th percentile)</li>
                        <li><strong>High:</strong> Score is between the 70th and 90th percentile</li>
                        <li><strong>Average:</strong> Score is between the 30th and 70th percentile</li>
                        <li><strong>Low:</strong> Score is between the 10th and 30th percentile</li>
                        <li><strong>Very Low:</strong> Score is in the bottom 10% of the population (<10th percentile)</li>
                    </ul>
                </p>
                <p>
                    <strong>Important Note:</strong> PRS comparisons should be made against ancestry-matched populations 
                    for the most accurate risk assessment. PRS distributions vary significantly between populations 
                    due to differences in genetic background.
                </p>
            </div>
            
            <div class="section">
                <h2>Limitations</h2>
                <p>
                    Polygenic risk scores have several important limitations:
                    <ul>
                        <li>PRS capture only the genetic component of risk and do not account for environmental, 
                            lifestyle, or other non-genetic factors.</li>
                        <li>The predictive ability of PRS varies by trait and is generally modest for complex diseases.</li>
                        <li>PRS performance may be reduced for individuals whose genetic ancestry differs from 
                            the populations used to develop the score.</li>
                        <li>A high PRS does not guarantee that someone will develop a condition, nor does a low 
                            PRS guarantee that they will not.</li>
                    </ul>
                </p>
            </div>
        </body>
        </html>
        """
        
        # Write to file
        with open(output_file, 'w') as f:
            f.write(html_content)
            
        return output_file
