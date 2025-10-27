#!/bin/bash
set -e

# FULL_IMAGE is exported from the calling .bat script
if [ -z "$FULL_IMAGE" ]; then
    echo "Error: FULL_IMAGE environment variable is not set." >&2
    exit 1
fi

echo "--- FEniCSx Docker Image Verification ---"

# --- 1. Check DolfinX Version (No Change Needed) ---
echo "Checking dolfinx version inside the container: $FULL_IMAGE"
echo "------------------------------------------------"
# This command already outputs the version clearly.
if ! docker run --rm "$FULL_IMAGE" python3 -c "import dolfinx; print('dolfinx version: ' + dolfinx.__version__)" ; then
    echo " ERROR: Failed to retrieve dolfinx version."
    exit 1
fi
echo "------------------------------------------------"


# --- 2. Check mpirun Version ---
echo "Checking mpirun (Open MPI) version inside the container (Printing first 3 lines of output):"
echo "------------------------------------------------"
# New Command: Capture the first 3 lines of mpirun --version output, redirecting stderr to stdout.
if ! docker run --rm "$FULL_IMAGE" mpirun --version 2>&1 | head -n 3 ; then
    echo " ERROR: Failed to retrieve mpirun version."
    exit 1
fi
echo "------------------------------------------------"

echo " All version checks passed successfully."