#!/bin/bash
set -e

DOCKERFILE=Dockerfile

# Create Dockerfile
cat > "$DOCKERFILE" <<EOF
FROM dolfinx/dolfinx:v0.9.0
RUN apt-get update && \
    apt-get install -y xvfb libgl1 libglu1-mesa mesa-utils && \
    python3 -m pip install --upgrade pip && \
    pip3 install pandas imageio pyvista && \
    rm -rf /var/lib/apt/lists/*
EOF

#     In case you have custom packages, add to the Dockerfile (inside the .bat)
#     python3 -m pip install ./src/pyproject.toml && \
#     apt update 

# Build image
if docker image inspect fenicsx:v0.9.0 >/dev/null 2>&1; then
    echo "Docker image fenicsx:v0.9.0 exists. Deleting..."
    docker rmi -f fenicsx:v0.9.0
fi

echo "Building Docker image..."
docker build -f "$DOCKERFILE" -t fenicsx:v0.9.0 .

rm -f "$DOCKERFILE"

# Test Docker installation
docker run --rm hello-world
