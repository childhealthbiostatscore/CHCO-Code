import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import sys
import logging
from scipy import stats # For Mann-Whitney U test and Kruskal-Wallis

# =========================================================================
# --- USER CONFIGURATION START ---
# Define which groups to load and combine from the results directory.
# OPTIONS:
# 1. To plot ONLY Dapagliflozin: TARGET_GROUPS = ['Dapagliflozin']
# 2. To plot ONLY Control:      TARGET_GROUPS = ['Control']
# 3. To plot BOTH groups:       TARGET_GROUPS = ['Control', 'Dapagliflozin']
TARGET_GROUPS = ['Control'] # <-- Set to plot ONLY Control group
# Input/Output paths
INPUT_DIR = '/gscratch/togo/rameshsh/sceptic_results/' 
OUTPUT_DIR = '/gscratch/togo/rameshsh/sceptic_results/' 
# --- END CONFIGURATION ---
# =========================================================================

# Define the consistent order for plotting the combined labels
COMBINED_LABEL_NAMES = [
    'PRE Control', 
    'POST Control', 
    'PRE Dapagliflozin', 
    'POST Dapagliflozin'
]

# Set up logging
logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')


# Helper function to get significance markers
def get_sig_marker(p_value):
    """Returns significance stars based on common thresholds."""
    if p_value < 0.001:
        return '***'
    elif p_value < 0.01:
        return '**'
    elif p_value < 0.05:
        return '*'
    else:
        return 'ns' # not significant


logger.info("STEP 1: Loading and combining results data...")
# (Loading and combining logic remains the same)
all_dfs = []
for group in TARGET_GROUPS:
    csv_file_path = os.path.join(INPUT_DIR, f'sceptic_results_{group}.csv')
    if not os.path.exists(csv_file_path):
        logger.warning(f"Input file not found for {group} group: {csv_file_path}. Skipping.")
    else:
        try:
            df = pd.read_csv(csv_file_path)
            all_dfs.append(df)
            logger.info(f"  ✓ Loaded data for {group} group ({len(df):,} cells).")
        except Exception as e:
            logger.error(f"✗ Failed to read CSV for {group}: {e}")

if not all_dfs:
    logger.error("✗ Failed to load any group data. Please ensure the analysis script has run for the specified groups.")
    sys.exit(1)
    
results_df = pd.concat(all_dfs, ignore_index=True)
logger.info(f"✓ Data combined. Total cells for plotting: {len(results_df):,}")

# --- Data Cleaning and Preparation ---
results_df['combined_label'] = results_df['combined_label'].str.strip()
present_labels = results_df['combined_label'].unique().tolist()
ACTUAL_PLOTTING_ORDER = [label for label in COMBINED_LABEL_NAMES if label in present_labels]

if not ACTUAL_PLOTTING_ORDER:
    logger.error("✗ After cleaning, no valid combined labels remain in the data for plotting.")
    sys.exit(1)

results_df['combined_label'] = pd.Categorical(
    results_df['combined_label'], 
    categories=COMBINED_LABEL_NAMES, 
    ordered=True
)
results_df.dropna(subset=['combined_label'], inplace=True)


# Set seaborn style and palette
sns.set_style("whitegrid")
sns.set_context("poster", font_scale=1.0) 
color_palette = {
    'PRE Control': '#d62728',         # Red
    'POST Control': '#1f77b4',        # Blue
    'PRE Dapagliflozin': '#2ca02c',   # Green
    'POST Dapagliflozin': '#9467bd'   # Purple
}


logger.info("\nSTEP 2: Generating 4-Group Pseudotime Plot...")

try:
    # --- STEP 2a: Run Omnibus Kruskal-Wallis H-test, calculate medians, and prepare result string ---
    data_dict = {}
    median_strings = []

    for label in COMBINED_LABEL_NAMES:
        data = results_df[results_df['combined_label'] == label]['pseudotime'].values
        if len(data) > 0:
            data_dict[label] = data
            median_val = np.median(data)
            median_strings.append(f"{label}: {median_val:.3f}")

    valid_labels_to_test = list(data_dict.keys())
    valid_groups_to_test = [data_dict[label] for label in valid_labels_to_test]
    
    kw_result_string = "Kruskal-Wallis Test: Not Applicable (Need >= 2 Groups)"
    kw_h_stat = None
    kw_p_value = None
    
    if len(valid_groups_to_test) >= 2:
        kw_h_stat, kw_p_value = stats.kruskal(*valid_groups_to_test)
        
        # Format P-value for clear display on the plot
        if kw_p_value < 0.0001:
            kw_p_str = " < 0.0001"
        else:
            kw_p_str = f" = {kw_p_value:.4f}"
            
        kw_result_string = f"Kruskal-Wallis Test\nH={kw_h_stat:.3f}, p{kw_p_str}"
        logger.info(f"  KW Omnibus Result: H={kw_h_stat:.3f}, p{kw_p_str}")
        
    median_result_string = "Group Medians:\n" + "\n".join(median_strings)
    
    # Combine KW and Median into a single box string for the top-right annotation
    main_stat_string = f"{kw_result_string}\n\n{median_result_string}"


    # (Heatmap plotting logic remains the same)
    logger.info("  Generating 4x2 Classification Summary (True Combined Label vs. Predicted Timepoint)...")
    code_to_timepoint = {0: 'PRE', 1: 'POST'}
    results_df['predicted_timepoint_name'] = results_df['predicted_label_code'].map(code_to_timepoint)
    raw_matrix = pd.crosstab(
        results_df['combined_label'], 
        results_df['predicted_timepoint_name'],
        dropna=False
    ).reindex(index=ACTUAL_PLOTTING_ORDER, columns=['PRE', 'POST'], fill_value=0) 
    norm_matrix = raw_matrix.div(raw_matrix.sum(axis=1), axis=0)

    plt.figure(figsize=(8, 7)) 
    ax_cm = sns.heatmap(
        norm_matrix, annot=True, cmap="YlGnBu", fmt='.2f', linewidths=.5, linecolor='black',
        cbar_kws={'label': 'Proportion Predicted'}, annot_kws={"size": 14}
    )
    for i in range(norm_matrix.shape[0]):
        for j in range(norm_matrix.shape[1]):
            count = raw_matrix.iloc[i, j]
            ax_cm.text(j + 0.5, i + 0.5, f'\n({count:,})', ha='center', va='top', color='black', fontsize=10)
    ax_cm.set_title("Classification Summary: True Group vs. Predicted Timepoint", fontsize=16)
    ax_cm.set_ylabel("True Combined Label", fontsize=14)
    ax_cm.set_xlabel("Sceptic's Predicted Timepoint (2 Classes)", fontsize=14)
    ax_cm.set_yticklabels(ax_cm.get_yticklabels(), rotation=0)

    output_filename_cm = os.path.join(OUTPUT_DIR, 'CM_Control.png')
    plt.tight_layout()
    plt.savefig(output_filename_cm, dpi=300, bbox_inches='tight')
    plt.close() 
    logger.info(f"✓ Classification summary plot saved to: {output_filename_cm}")
    
    
    # --- Pseudotime Distribution Plot with Significance Annotation ---
    fig, axes = plt.subplots(
        2, 1, sharex=True, figsize=(12, 12), 
        gridspec_kw={'height_ratios': [1, 3]}
    )
    plt.subplots_adjust(hspace=0.4) 
    
    # --- Box Plot (Top Subplot) ---
    ax_box = axes[0]
    sns.boxplot(
        x="pseudotime", y="combined_label", data=results_df, order=ACTUAL_PLOTTING_ORDER, 
        palette=color_palette, ax=ax_box, orient='h', notch=True, width=0.8, fliersize=3
    )
    
    # Get y-positions for annotations based on the plotting order (not needed for drawing lines now)
    y_map = {label: i for i, label in enumerate(ACTUAL_PLOTTING_ORDER)}
    
    # --- Perform Mann-Whitney U test and prepare pairwise result string ---
    
    # Define all key comparisons
    ALL_KEY_COMPARISONS = [
        ('PRE Control', 'POST Control'),
        ('PRE Dapagliflozin', 'POST Dapagliflozin'),
        ('POST Control', 'POST Dapagliflozin'), 
        ('PRE Control', 'PRE Dapagliflozin'),
    ]
    
    logger.info("\nRunning Pairwise Mann-Whitney U tests for key comparisons:")
    
    pairwise_stat_lines = ["Pairwise Mann-Whitney U Tests:"]
    comparison_results = {} # Map to store full results for logging (not plotting)

    for label1, label2 in ALL_KEY_COMPARISONS:
        
        if label1 in y_map and label2 in y_map:
            
            data1 = results_df[results_df['combined_label'] == label1]['pseudotime'].values
            data2 = results_df[results_df['combined_label'] == label2]['pseudotime'].values
            
            if len(data1) > 0 and len(data2) > 0:
                # Mann-Whitney U Test
                _, p_value = stats.mannwhitneyu(data1, data2, alternative='two-sided')
                marker = get_sig_marker(p_value)
                
                # Log the result for the user in the console
                logger.info(f"  {label1} vs {label2} (MWU): P={p_value:.8f}, Marker: {marker}")
                
                # Format for the plot box (Shortening names to fit nicely)
                short_l1 = label1.replace('Dapagliflozin', 'Dapa').replace('Control', 'Ctrl')
                short_l2 = label2.replace('Dapagliflozin', 'Dapa').replace('Control', 'Ctrl')
                
                # Append to the pairwise stats box string
                pairwise_stat_lines.append(f"  {short_l1} vs {short_l2}: p={p_value:.3g} ({marker})")
                
                # Store for full results logging
                comparison_results[(label1, label2)] = {'p_value': p_value, 'marker': marker}
            else:
                logger.warning(f"  Cannot run MWU: Missing data for {label1} or {label2}.")
        else:
            logger.warning(f"  Skipping test: Missing labels {label1} or {label2} in data.")
            
    pairwise_result_string = "\n".join(pairwise_stat_lines)
    
    # --- Annotation Drawing (No lines, just cleaning up the plot) ---
    
    # Clean up the box plot
    ax_box.set_ylabel("Time and Group", fontsize=16) 
    ax_box.set_xlabel("") 
    ax_box.tick_params(axis='x', which='both', bottom=True, top=False, labelbottom=True)
    ax_box.set_xlim(0, 1)
    ax_box.grid(axis='y', linestyle='--', alpha=0.7)
    
    
    # --- Density Plot (Bottom Subplot) ---
    ax_dens = axes[1]
    for label_name in ACTUAL_PLOTTING_ORDER:
        group_data = results_df[results_df['combined_label'] == label_name]
        if not group_data.empty:
            sns.kdeplot(
                x=group_data['pseudotime'], ax=ax_dens, fill=True, alpha=0.6, 
                linewidth=1.5, color=color_palette[label_name], label=label_name, 
                cut=0, clip=(0, 1)
            )

    ax_dens.set_xlim(0, 1) 
    ax_dens.set_xlabel("Sceptic's Pseudotime (Weighted Average)", fontsize=16)
    ax_dens.set_ylabel("Probability Density", fontsize=16)
    ax_dens.legend(title="Group", loc='upper right', frameon=True, fontsize=12)
    ax_dens.grid(axis='x', linestyle='--', alpha=0.7)

    # --- Add KW/Median Result Annotation (Top Right) ---
    fig.suptitle(
        "Pseudotime Distribution by Timepoint and Group", 
        fontsize=20, 
        fontweight='bold',
        y=0.97 # Positioned high above both boxes
    )
    
    # 1. Place the main statistics (KW + Medians) box (Top Right)
    plt.figtext(
        0.15, 0.55, # Adjusted Y position to avoid title
        main_stat_string,
        fontsize=12,
        fontweight='bold',
        ha='right',
        bbox=dict(boxstyle="round,pad=0.5", fc="#e0e0e0", alpha=0.9, ec="gray") 
    )

    # 2. Place the pairwise comparison results box (Top Left)
    plt.figtext(
        0.65, 0.85, 
        pairwise_result_string,
        fontsize=12,
        fontweight='bold',
        ha='left',
        bbox=dict(boxstyle="round,pad=0.5", fc="#e0e0e0", alpha=0.9, ec="gray") 
    )
    
    # Save the figure
    output_filename = os.path.join(OUTPUT_DIR, 'pseudotime_distribution_annotated_Control.png')
    plt.tight_layout(rect=[0, 0.03, 1, 0.9])
    plt.savefig(output_filename, dpi=300, bbox_inches='tight')
    plt.close(fig)
    logger.info(f"✓ Plot saved successfully to: {output_filename}")


    # --- STEP 3: Omnibus Kruskal-Wallis H-test (Logging Only) ---
    logger.info("\nSTEP 3: Omnibus Kruskal-Wallis H-test Summary (from plot annotation data):")

    if kw_h_stat is not None:
        logger.info("-" * 60)
        logger.info(f"Omnibus Kruskal-Wallis H-test ({len(valid_groups_to_test)} Groups):")
        logger.info(f"  Groups Tested: {', '.join(valid_labels_to_test)}")
        logger.info(f"  H Statistic: {kw_h_stat:.4f}")
        logger.info(f"  Omnibus P-value: {kw_p_value:.10f}")
        
        if kw_p_value < 0.05:
            logger.info("  CONCLUSION: The overall pseudotime medians are significantly different (p < 0.05) across the tested groups.")
        else:
            logger.info("  CONCLUSION: The overall pseudotime medians are not significantly different (p >= 0.05) across the tested groups.")
        logger.info("-" * 60)
            
    else:
        logger.error("✗ Kruskal-Wallis test was skipped due to insufficient groups (less than 2).")


except Exception as e:
    logger.error(f"✗ Failed to generate plot or run analysis: {e}", exc_info=True)
    sys.exit(1)

logger.info("\n" + "="*80)
logger.info("PLOTTING AND ANALYSIS COMPLETE!")
logger.info("="*80)
