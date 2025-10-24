#!/bin/bash
#
# Simplified Docker Build Script for FEniCSx Environment
# This script is executed inside the WSL2 distribution.

# Variables are passed from the Windows .bat file via 'wsl -e'
# CUSTOM_PACKAGE_FLAG is the key differentiator.
IMAGE_NAME=$IMAGE_NAME
TAG=$TAG
CUSTOM_PACKAGE_FLAG=$CUSTOM_PACKAGE_FLAG # 1 if pyproject.toml found, 0 otherwise
FULL_IMAGE="${IMAGE_NAME}:${TAG}"
DOCKERFILE_PATH="./Dockerfile.fenicsx"

# The build context is the current directory (set by the .bat file's 'cd')
TARGET_DIR="." 

echo "Starting Docker image build for ${FULL_IMAGE}..."

# --- 1. Define the parameterized Dockerfile content using a single HERE-document ---
# We use 'EODF' with quotes to prevent shell variable substitution inside the Dockerfile block
# (like $TAG), then use it selectively outside to pass the FEniCSx version (dolfinx/dolfinx:${TAG})

cat <<'EODF' > "${DOCKERFILE_PATH}"
# Base image from FEniCSx Project, parameterized by the script's TAG
FROM dolfinx/dolfinx:TAG_PLACEHOLDER

# Install system dependencies (for plotting/GUI/utility)
# RUN apt-get update && \
#     apt-get install -y xvfb libgl1 libglu1-mesa mesa-utils && \
#     python3 -m pip install --upgrade pip
RUN apt-get update && \
    apt-get install -y xvfb libgl1-mesa-dri libglx-mesa0 libglu1-mesa mesa-utils && \
    python3 -m pip install --upgrade pip

# Other standard Python packages used by common demos/post-processing
RUN pip3 install pandas imageio pyvista meshio h5py && \
    rm -rf /var/lib/apt/lists/*
EODF

# --- 2. Inject custom package installation if the flag is set ---
# Since you stated "Custom package always exists" in your prompt, 
# this block *should* be set to always run or the .bat file should be adjusted.
# We will honor the existing conditional logic but simplify the steps.
if [ "$CUSTOM_PACKAGE_FLAG" -eq "1" ]; then
    echo "Custom package flag is set. Injecting custom package install steps."
    
    cat <<'EOC' >> "${DOCKERFILE_PATH}"

# --- Custom Package Installation ---
# Copy the custom package source folder relative to the build context (TARGET_DIR)
COPY src_install/custom_packages /tmp/custom_packages

# Install the custom Python package from the copied directory
RUN pip3 install /tmp/custom_packages

# --- End Custom Package Installation ---
EOC
else
    echo "No custom package flag set (CUSTOM_PACKAGE_FLAG=0). Skipping custom install."
fi

# --- 3. Add Final Configuration and Replace Placeholder ---
cat <<'EODF_FINAL' >> "${DOCKERFILE_PATH}"

# Final image configuration
WORKDIR /home/fenics
EODF_FINAL

# Replace the placeholder with the actual TAG variable (FENICSx base image version)
sed -i "s|TAG_PLACEHOLDER|${TAG}|g" "${DOCKERFILE_PATH}"

echo "Final Dockerfile written to ${DOCKERFILE_PATH}"

# --- 4. Build the Docker image ---
echo "Building Docker image ${FULL_IMAGE} from ${DOCKERFILE_PATH}"

# Build command: -f specifies the Dockerfile, . specifies the build context
if docker build -f "${DOCKERFILE_PATH}" -t "${FULL_IMAGE}" .; then
    echo "Docker build successful."
else
    echo "ERROR: Docker build failed. Review the output above." >&2
    # Clean up Dockerfile on failure
    rm -f "${DOCKERFILE_PATH}"
    exit 1
fi

# --- 5. Cleanup ---
rm -f "${DOCKERFILE_PATH}"

echo "Docker image build script finished."