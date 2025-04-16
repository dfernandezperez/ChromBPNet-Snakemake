import sys
import pandas as pd
import numpy as np # For infinity if needed
import logging

# Args passed from Snakemake
input_peaks_file = snakemake.input.peaks
output_peaks_file = snakemake.output.top_peaks
log_file = snakemake.log[0]
max_peaks = int(snakemake.params.top_n_peaks) # Ensure it's an integer

# Setup logging
logging.basicConfig(filename=log_file, level=logging.INFO,
                    format='%(asctime)s %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
logging.info(f"Starting peak filtering.")
logging.info(f"Input file: {input_peaks_file}")
logging.info(f"Output file: {output_peaks_file}")
logging.info(f"Max peaks requested (N): {max_peaks}")

try:
    # --- Step 1: Determine the score threshold ---
    logging.info("Reading scores from column 7...")
    # Use iterator and chunking for potentially large files
    scores = []
    try:
        # Read only the 7th column (index 6), assuming tab-separated, no header
        # Use float64 for scores to handle scientific notation properly
        score_iterator = pd.read_csv(input_peaks_file, sep='\t', header=None, usecols=[6],
                                     chunksize=100000, iterator=True, dtype={6: np.float64})
        for chunk in score_iterator:
            scores.extend(chunk[6].tolist())
    except pd.errors.EmptyDataError:
        logging.warning(f"Input file {input_peaks_file} is empty or contains no valid score data.")
        scores = []
    except Exception as e:
         logging.error(f"Error reading scores from {input_peaks_file}: {e}")
         sys.exit(1)


    score_threshold = -np.inf # Default if file is empty or has < N peaks

    if scores:
        scores.sort(reverse=True) # Sort scores descending (highest first)
        num_actual_peaks = len(scores)
        logging.info(f"Found {num_actual_peaks} peaks.")

        if num_actual_peaks >= max_peaks:
            # Get the score of the Nth peak (0-based index is N-1)
            score_threshold = scores[max_peaks - 1]
            logging.info(f"Score threshold for top {max_peaks} peaks: {score_threshold}")
        else:
            # Fewer than N peaks, keep them all (threshold is the lowest score)
            score_threshold = scores[-1]
            logging.info(f"Fewer than {max_peaks} peaks found. Using minimum score as threshold: {score_threshold}")
    else:
         logging.warning("No scores found. Output file will be empty.")
         # Create empty output file and exit cleanly
         with open(output_peaks_file, 'w') as f_out:
             pass # Creates an empty file
         logging.info("Created empty output file.")
         sys.exit(0)


    # --- Step 2: Filter the original file using the threshold ---
    logging.info(f"Filtering original file using threshold >= {score_threshold}...")
    written_count = 0
    # Read input and write output line by line to preserve order and save memory
    with open(input_peaks_file, 'r') as f_in, open(output_peaks_file, 'w') as f_out:
        for line in f_in:
            try:
                fields = line.strip().split('\t')
                # Check if line has enough fields before accessing index 6
                if len(fields) > 6:
                    current_score = float(fields[6]) # 7th column is index 6
                    if current_score >= score_threshold:
                        f_out.write(line)
                        written_count += 1
                else:
                    logging.warning(f"Skipping malformed line (fewer than 7 columns): {line.strip()}")
            except ValueError:
                # Handle cases where score is not a valid number
                 logging.warning(f"Skipping line with non-numeric score in column 7: {line.strip()}")
            except Exception as e:
                 logging.error(f"Error processing line: {line.strip()} - {e}")
                 # Decide whether to continue or exit based on error type if desired

    logging.info(f"Finished filtering. Wrote {written_count} peaks to {output_peaks_file}.")

except Exception as e:
    logging.exception("An error occurred during peak filtering:") # Logs traceback
    sys.exit(1)