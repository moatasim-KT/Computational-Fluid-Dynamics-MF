#!/bin/bash

# Set error handling
set -e

# Function to log messages
log_message() {
  echo "[$(date)] $1" | tee -a conversion_log.txt
}

log_message "Starting MATLAB to Python conversion process"

# Install SMOP if not already installed
if ! pip list | grep -q smop; then
  log_message "Installing SMOP (MATLAB to Python converter)..."
  pip install smop
  log_message "SMOP installation complete"
else
  log_message "SMOP is already installed"
fi

# Create a directory for Python output
log_message "Creating Python output directory..."
mkdir -p python_output

# Finding all MATLAB files
log_message "Locating MATLAB files..."
matlab_files=$(find . -name "*.m" -o -name "*.mlx")

# Check if we found any files
if [ -z "$matlab_files" ]; then
  log_message "No MATLAB files found in the repository."
  exit 1
fi

log_message "Found MATLAB files, beginning conversion..."

# Convert each file
for file in $matlab_files; do
  log_message "Converting $file to Python..."
  python_file="python_output/$(basename "${file%.*}.py")"
  
  # Make sure the parent directory exists
  mkdir -p "$(dirname "$python_file")"
  
  # Try to convert the file using smop command directly
  # Method 1: Try using smop directly as a command
  smop "$file" > "$python_file" 2>> conversion_log.txt || {
    log_message "Direct smop command failed, trying alternative method..."
    
    # Method 2: Try using python -c to invoke SMOP programmatically
    python -c "from smop.main import main; main(['$file'])" > "$python_file" 2>> conversion_log.txt || {
      log_message "Failed to convert $file. See conversion_log.txt for details."
      continue
    }
  }
  
  log_message "Created $python_file"
done

log_message "Conversion process complete."
log_message "Python files are stored in the python_output directory."
log_message "Check conversion_log.txt for the detailed log."