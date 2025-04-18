import sys
import tensorflow as tf
import logging
import os
import platform

# Args passed from Snakemake
log_file = snakemake.log[0]
output_flag_file = snakemake.output.gpu_check_flag # To signal success

# Setup logging
logging.basicConfig(
    filename=log_file,
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)

logging.info("--- Starting TensorFlow GPU Check ---")
logging.info(f"Python Version: {sys.version}")
logging.info(f"TensorFlow Version: {tf.__version__}")
logging.info(f"Platform: {platform.platform()}")
logging.info(f"Hostname: {os.uname().nodename}")

physical_gpus = []
all_ops_successful = False # Flag to track if operations worked on at least one GPU

try:
    physical_gpus = tf.config.list_physical_devices('GPU')
    logging.info(f"Physical GPUs found: {len(physical_gpus)}")

    if physical_gpus:
        all_ops_successful = True # Assume success initially
        for i, gpu in enumerate(physical_gpus):
            logging.info(f"--- GPU {i} ---")
            logging.info(f"Name (from list_physical_devices): {gpu.name}")

            # --- Construct device string for tf.device ---
            # Extract type ('GPU') and index ('0')
            try:
                # A common way, though might vary slightly across TF versions
                device_type = gpu.device_type
                # Get index - name might be '/physical_device:GPU:0'
                device_index_str = gpu.name.split(':')[-1]
                device_index = int(device_index_str)
                tf_device_string = f"/{device_type}:{device_index}" # e.g., /GPU:0
                logging.info(f"Constructed tf.device string: {tf_device_string}")
            except Exception as e_parse:
                logging.error(f"Could not parse device name '{gpu.name}' to create tf.device string: {e_parse}")
                all_ops_successful = False
                continue # Skip trying operations on this device

            # --- Get Details (Optional) ---
            try:
                details = tf.config.experimental.get_device_details(gpu)
                logging.info(f"Details: {details}")
            except Exception as detail_err:
                logging.warning(f"Could not get experimental device details for GPU {i}: {detail_err}")

            # --- Try Operation ---
            try:
                logging.info(f"Attempting operation on {tf_device_string}...")
                # Use the constructed string with tf.device
                with tf.device(tf_device_string):
                    a = tf.constant([[1.0, 2.0], [3.0, 4.0]], dtype=tf.float32)
                    b = tf.matmul(a, a)
                logging.info(f"Successfully performed matmul on {tf_device_string}. Result device: {b.device}")
                del a, b
            except Exception as op_err:
                 logging.error(f"Failed to perform operation on {tf_device_string}: {op_err}")
                 all_ops_successful = False # Mark failure if any device op fails
                 # Optionally break here if one failure means the whole test fails

    else:
        logging.warning("TensorFlow did not find any Physical GPUs configured.")
        all_ops_successful = False

except Exception as e_list:
    logging.exception("Error occurred while listing physical GPU devices:")
    all_ops_successful = False

# --- Signal Success or Failure ---
if physical_gpus and all_ops_successful:
    logging.info(f"--- TensorFlow GPU Check Successful ({len(physical_gpus)} device(s) found and usable) ---")
    with open(output_flag_file, 'w') as f:
        f.write("TensorFlow GPU check successful.\n")
        f.write(f"Found {len(physical_gpus)} physical device(s).\n")
else:
    logging.error("--- TensorFlow GPU Check Failed (No physical GPUs found or operation failed) ---")
    sys.exit(1) # Fail the rule