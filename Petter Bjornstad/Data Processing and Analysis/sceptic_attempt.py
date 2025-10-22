import numpy as np
import pandas as pd
import boto3
import json
from datetime import datetime
import os
import tempfile
import logging
import sys
from sklearn import preprocessing
from sklearn.model_selection import train_test_split 


# =========================================================================
# To switch the analysis group, change the value below:
TARGET_GROUP = 'Dapagliflozin' 
# =========================================================================

CELL_DOWNSAMPLE_N = None 
GROUP_LABELS_S3_KEY = 'cleaned_data/Attempt_sceptic_labels_group_PT.txt' 

# Map the actual numeric timepoint labels found in the file to their names.
TIMEPOINT_MAPPING = {
    '0': 'PRE',
    '1': 'POST'
}

# Map the actual group codes found in the file ('B') to their group names.
GROUP_CODE_MAPPING = {
    'A': 'Control',
    'B': 'Dapagliflozin' 
    # If there are other codes, they should be added here!
}

# Define the consistent order for plotting the combined labels
COMBINED_LABEL_NAMES = [
    'PRE Control', 
    'POST Control', 
    'PRE Dapagliflozin', 
    'POST Dapagliflozin'
]
# --- END CONFIGURATION ---

# Set up logging
log_dir = '/gscratch/togo/rameshsh/sceptic_logs/'
os.makedirs(log_dir, exist_ok=True)

log_file = os.path.join(log_dir, f'sceptic_run_{TARGET_GROUP}_{datetime.now().strftime("%Y%m%d_%H%M%S")}.log')

# Configure logging to write to both file and console
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler(log_file),
        logging.StreamHandler(sys.stdout)
    ]
)

logger = logging.getLogger(__name__)

logger.info("="*80)
logger.info(f"STARTING SCEPTIC PSEUDOTIME ANALYSIS (GROUP: {TARGET_GROUP} ONLY)")
logger.info(f"Target cell count: {CELL_DOWNSAMPLE_N if CELL_DOWNSAMPLE_N else 'All'}")
logger.info("="*80)

# Set up Kopah S3 connection
logger.info("STEP 1: Connecting to Kopah S3...")

try:
    with open('/mmfs1/home/rameshsh/keys.json', 'r') as f:
        keys = json.load(f)

    session = boto3.session.Session(
        aws_access_key_id=keys['MY_ACCESS_KEY'],
        aws_secret_access_key=keys['MY_SECRET_KEY']
    )
    s3 = session.client('s3', endpoint_url='https://s3.kopah.uw.edu')
    logger.info("✓ Successfully connected to Kopah S3")
except Exception as e:
    logger.error(f"✗ Failed to connect to Kopah S3: {e}", exc_info=True)
    sys.exit(1)


# Download files
logger.info("\nSTEP 2: Downloading data from Kopah S3...")

temp_dir = tempfile.mkdtemp()
data_file = os.path.join(temp_dir, 'Attempt_sceptic_PT.txt')
labels_file = os.path.join(temp_dir, 'Attempt_sceptic_labels_PT.txt')
group_labels_file = os.path.join(temp_dir, 'Attempt_sceptic_labels_group_PT.txt')


logger.info(f"  Downloading data and labels...")
try:
    s3.download_file('attempt', 'cleaned_data/Attempt_sceptic_PT.txt', data_file)
    s3.download_file('attempt', 'cleaned_data/Attempt_sceptic_labels_PT.txt', labels_file)
    s3.download_file('attempt', GROUP_LABELS_S3_KEY, group_labels_file)
    
    logger.info(f"  ✓ Downloads complete. Saved to temporary directory.")
except Exception as e:
    logger.error(f"✗ Failed to download data from Kopah S3: {e}", exc_info=True)
    import shutil
    shutil.rmtree(temp_dir)
    sys.exit(1)


# Load and preprocess data
logger.info("\nSTEP 3: Loading and processing data (CRITICAL LABEL CLEANING)...")
try:
    # Load feature matrix
    data_concat = np.loadtxt(data_file)
    total_cells, total_genes = data_concat.shape
    logger.info(f"  ✓ Raw Data shape: {total_cells:,} cells × {total_genes:,} genes.")

    # --- Timepoint Label Loading and Cleaning ---
    y_raw_full = np.loadtxt(labels_file, dtype=str) 
    # Strip whitespace
    y_raw_stripped = np.array([s.strip() for s in y_raw_full])
    logger.info(f"  [DEBUG] Raw Timepoint labels (first 10): {y_raw_full[:10]}")
    
    # APPLY TIMEPOINT MAPPING
    y_timepoint_name = np.array([TIMEPOINT_MAPPING.get(code, f'Unknown Timepoint {code}') for code in y_raw_stripped])
    logger.info(f"  [DEBUG] Mapped Timepoint names (first 10): {y_timepoint_name[:10]}")
    
    # The `label` array for SCEPTIC must still be the numeric encoding (0, 1), 
    # so we use the stripped numeric codes for `lab_encoder`.
    y_raw_numeric = y_raw_stripped 
    
    # --- Group Code Loading and Cleaning ---
    y_group_codes_full = np.loadtxt(group_labels_file, dtype=str) 
    # Strip whitespace
    y_group_codes_stripped = np.array([s.strip() for s in y_group_codes_full])
    logger.info(f"  [DEBUG] Raw Group codes (first 10): {y_group_codes_full[:10]}")
    
    # Apply mapping to get human-readable group names
    y_group_name = np.array([GROUP_CODE_MAPPING.get(code, f'Unknown Group {code}') for code in y_group_codes_stripped])
    logger.info(f"  [DEBUG] Mapped Group names (first 10): {y_group_name[:10]}")

    # Consistency Check
    if len(y_raw_numeric) != len(y_group_codes_stripped):
         logger.error(f"Timepoint labels ({len(y_raw_numeric):,}) and Group labels ({len(y_group_codes_stripped):,}) are inconsistent in size. Exiting.")
         sys.exit(1)
    
    # --- FILTERING FOR TARGET GROUP ---
    logger.info(f"  Filtering data to include only the '{TARGET_GROUP}' group...")

    filter_mask = (y_group_name == TARGET_GROUP)
    num_filtered = np.sum(filter_mask)

    if num_filtered == 0:
        logger.error(f"✗ Failed to find any cells belonging to the '{TARGET_GROUP}' group. Exiting.")
        sys.exit(1)

    # Apply filter to all relevant arrays
    data_concat = data_concat[filter_mask]
    y_raw_numeric = y_raw_numeric[filter_mask]
    y_timepoint_name = y_timepoint_name[filter_mask]
    y_group_name = y_group_name[filter_mask]

    # Recalculate combined labels as the arrays are now filtered
    combined_labels = [f'{time} {group}' for time, group in zip(y_timepoint_name, y_group_name)]

    total_cells, total_genes = data_concat.shape
    logger.info(f"  ✓ Filtering complete. New Data shape: {total_cells:,} cells × {total_genes:,} genes.")
    # --- END FILTERING ---
    
    # Re-encode the numeric timepoint codes for SCEPTIC training (0, 1) using the filtered data
    lab_encoder = preprocessing.LabelEncoder()
    # Use the cleaned numeric timepoint codes for encoding
    label = lab_encoder.fit_transform(y_raw_numeric.flatten())
    label = label.astype(np.int64)
    
    label_list = np.unique(label)
    
    logger.info(f"  ✓ Labels shape: {label.shape}. Encoded label list: {label_list}")
    
    # progression: 0 -> PRE, 1 -> POST (Based on TIMEPOINT_MAPPING)
    class_names = ['PRE', 'POST']
    ENCODED_MAPPING = {0: 'PRE', 1: 'POST'}

    
    #  Show the exact combined labels being generated
    logger.info(f"  [DEBUG] Final Combined Labels (first 10): {combined_labels[:10]}")
    
    # Downsampling (Logic remains the same, but only applied to the filtered data)
    if CELL_DOWNSAMPLE_N and total_cells > CELL_DOWNSAMPLE_N:
        logger.info(f"  Applying stratified downsampling to {CELL_DOWNSAMPLE_N:,} cells...")
        
        indices = np.arange(total_cells)

        sampled_indices, _, _, _ = train_test_split(
            indices, 
            label, 
            train_size=CELL_DOWNSAMPLE_N, 
            random_state=42, 
            stratify=label
        )
        
        data_concat = data_concat[sampled_indices]
        label = label[sampled_indices]
        
        # Must sample the raw and combined labels for accurate plotting
        y_timepoint_name = y_timepoint_name[sampled_indices]
        y_group_name = y_group_name[sampled_indices] 
        combined_labels = np.array(combined_labels)[sampled_indices].tolist() 
        
        logger.info(f"  ✓ Downsampling complete. New Data shape: {data_concat.shape}")
    else:
        if CELL_DOWNSAMPLE_N and total_cells <= CELL_DOWNSAMPLE_N:
             logger.info(f"  Downsampling skipped. Using all cells.")
        else:
             logger.info(f"  Downsampling skipped. Using all {total_cells:,} cells.")


    label_counts = np.bincount(label)
    logger.info(f"  Final Label distribution in sampled data:")
    for i, count in enumerate(label_counts):
        class_name = class_names[i]
        logger.info(f"    - Class {i} ({class_name}): {count:,} cells ({count/len(label)*100:.1f}%)")    
    
except Exception as e:
    logger.error(f"✗ Failed during data loading/preprocessing: {e}", exc_info=True)
    sys.exit(1)

# Run SCEPTIC
logger.info("\nSTEP 4: Running SCEPTIC analysis with XGBoost...")
logger.info("="*80)
try:
    from sklearn.model_selection import KFold, GridSearchCV, StratifiedKFold 
    # Assuming sceptic_fixed.py is in this path
    sys.path.insert(0, '/gscratch/togo/rameshsh') 
    import sceptic_fixed as sceptic
    logger.info("  ✓ SCEPTIC module imported successfully")
    
    parameters = {
        "max_depth": [3, 5],
        "learning_rate": [0.1, 0.3],
        "n_estimators": [100],
        "subsample": [0.8]
    }
    
    logger.info("  Configuration:")
    logger.info(f"    - Method: XGBoost")
    logger.info(f"    - External CV folds (eFold): 3")
    logger.info(f"    - Internal CV folds (iFold): 4")
    
    start_time = datetime.now()
    
    cm, label_predicted, pseudotime, sceptic_prob = sceptic.run_sceptic_and_evaluate(
        data=data_concat, 
        labels=label, 
        label_list=label_list, 
        parameters=parameters,
        method="xgboost",
        use_gpu=False
    )
    
    end_time = datetime.now()
    duration = end_time - start_time
    
    logger.info("  " + "-"*76)
    logger.info(f"  ✓ Analysis completed. Total runtime: {duration}")
    
except Exception as e:
    logger.error(f"✗ SCEPTIC analysis failed: {e}", exc_info=True)
    sys.exit(1)

# Analyze results
logger.info("\nSTEP 5: Analyzing results...")
try:
    class_labels = class_names
    
    logger.info(f"  Confusion matrix shape: {cm.shape}")
    logger.info(f"  Confusion matrix:\n{cm}")
    
    if np.sum(cm) > 0:
        accuracy = np.trace(cm) / np.sum(cm)
        logger.info(f"\n  Overall accuracy (CV): {accuracy:.3f} ({accuracy*100:.1f}%)")
    else:
        logger.warning("\n  Confusion matrix sum is zero. Accuracy could not be calculated.")
    
    logger.info(f"\n  Pseudotime statistics (Overall):")
    logger.info(f"    - Range: [{pseudotime.min():.4f}, {pseudotime.max():.4f}]")
    logger.info(f"    - Mean: {pseudotime.mean():.4f}")
    
   
    for i, class_name in enumerate(class_labels):
        class_pseudotime = pseudotime[label == i]
        logger.info(f"\n  Pseudotime for {class_name} (Class {i}):")
        logger.info(f"    - Mean: {class_pseudotime.mean():.4f}")
        logger.info(f"    - Std: {class_pseudotime.std():.4f}")
        logger.info(f"    - Median: {np.median(class_pseudotime):.4f}")
    
except Exception as e:
    logger.error(f"✗ Failed to analyze results: {e}", exc_info=True)
    pass

# Save outputs 
logger.info("\nSTEP 6: Saving numerical outputs...")
OUTPUT_DIR = '/gscratch/togo/rameshsh/sceptic_results/' 
os.makedirs(OUTPUT_DIR, exist_ok=True)
logger.info(f"  Output directory: {OUTPUT_DIR}")
try:
    
    # Save files needed by the plotting script
    # --- MODIFICATION: Group-specific filenames added ---
    np.savetxt(os.path.join(OUTPUT_DIR, f'confusion_matrix_{TARGET_GROUP}.txt'), cm, fmt='%d')
    np.savetxt(os.path.join(OUTPUT_DIR, f'pseudotime_{TARGET_GROUP}.txt'), pseudotime, fmt='%.6f')
    
    # Create and save results DataFrame
    code_to_name = ENCODED_MAPPING
    results_df = pd.DataFrame({
        'true_timepoint_code': label, 
        'true_timepoint_name': y_timepoint_name, # Use the mapped name
        'true_group_name': y_group_name, # Use the mapped name
        'combined_label': combined_labels,
        'predicted_label_code': label_predicted.astype(int),
        'predicted_label_name': pd.Series(label_predicted.astype(int)).map(code_to_name).to_numpy(),
        'pseudotime': pseudotime,
    })
    
    for i, name in enumerate(class_names):
        results_df[f'prob_{name.lower()}'] = sceptic_prob[:, i]
        
    results_df.to_csv(os.path.join(OUTPUT_DIR, f'sceptic_results_{TARGET_GROUP}.csv'), index=False)
    logger.info(f"  ✓ Saved core output files for {TARGET_GROUP} group.")
    logger.info(f"    - confusion_matrix_{TARGET_GROUP}.txt")
    logger.info(f"    - pseudotime_{TARGET_GROUP}.txt")
    logger.info(f"    - sceptic_results_{TARGET_GROUP}.csv")
    # --- END MODIFICATION ---
    
except Exception as e:
    logger.error(f"✗ Failed to save outputs: {e}", exc_info=True)
    pass

# Clean up
logger.info("\nSTEP 7: Cleaning up temporary files...")
try:
    import shutil
    shutil.rmtree(temp_dir)
    logger.info("  ✓ Temporary files removed")
except Exception as e:
    logger.warning(f"⚠ Failed to remove temporary files: {e}")

# Final summary
logger.info("\n" + "="*80)
logger.info(f"ANALYSIS COMPLETE FOR GROUP: {TARGET_GROUP}!")
logger.info("Now run sceptic_plots.py to generate visualizations.")
logger.info("="*80)
logger.info(f"\nAll outputs saved to: {OUTPUT_DIR}")
logger.info(f"Log file: {log_file}")
logger.info("="*80)
sys.exit(0)
