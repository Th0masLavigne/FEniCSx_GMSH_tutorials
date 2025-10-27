#!/bin/bash
set -e

# --- 1. Robust Check for Docker Desktop Integration (Preferred) ---
if command -v docker >/dev/null 2>&1; then
    echo "Docker command found. Checking if Docker daemon is accessible..."
    
    # Check if the Docker daemon is running/reachable (likely via Docker Desktop integration).
    if timeout 10 docker info >/dev/null 2>&1; then
        echo "Docker daemon is reachable (likely via Docker Desktop integration). Skipping local installation and service start."
        
        # Ensure the current user is in the 'docker' group.
        sudo groupadd docker || true 
        sudo usermod -aG docker "$USER"
        echo "User added to 'docker' group."
        
        exit 0
    else
        echo "Docker command found but daemon is unreachable (or timing out). Proceeding with local installation..."
    fi
else
    echo "Docker command not found. Proceeding with installation..."
fi

# --- 2. Local Docker Engine Installation (Fallback for non-Docker Desktop users) ---

# Enable Systemd in WSL for proper service management (Ubuntu 24.04 standard)
echo "Ensuring Systemd is enabled for proper service start..."
if [ ! -f /etc/wsl.conf ] || ! grep -q "systemd=true" /etc/wsl.conf; then
    echo "[boot]" | sudo tee -a /etc/wsl.conf > /dev/null
    echo "systemd=true" | sudo tee -a /etc/wsl.conf > /dev/null
    echo "NOTE: Systemd configuration added to /etc/wsl.conf. WSL Distribution restart may be required for full stability."
fi

# Remove old Docker packages
for pkg in docker.io docker-doc docker-compose docker-compose-v2 podman-docker containerd runc; do sudo apt-get -y remove $pkg; done

# Install prerequisites
sudo apt-get update
sudo apt-get install -y ca-certificates curl

# Add Docker's official GPG key:
sudo install -m 0755 -d /etc/apt/keyrings
sudo curl -fsSL https://download.docker.com/linux/ubuntu/gpg -o /etc/apt/keyrings/docker.asc
sudo chmod a+r /etc/apt/keyrings/docker.asc

# Add the repository to Apt sources:
echo \
  "deb [arch=$(dpkg --print-architecture) signed-by=/etc/apt/keyrings/docker.asc] https://download.docker.com/linux/ubuntu \
  $(. /etc/os-release && echo "$VERSION_CODENAME") stable" | \
  sudo tee /etc/apt/sources.list.d/docker.list > /dev/null
sudo apt-get update

# Install Docker packages
sudo apt-get install -y docker-ce docker-ce-cli containerd.io docker-buildx-plugin docker-compose-plugin

# Post-installation steps
sudo groupadd docker || true # Add group if it doesn't exist
sudo usermod -aG docker "$USER"
echo "User added to 'docker' group."

# --- 3. Robust Service Start ---

# Attempt to start the service using systemctl (if Systemd is active)
echo "Attempting to start Docker service via systemctl..."
if sudo systemctl start docker.service 2>/dev/null; then
    echo "Docker service started via systemctl."
    exit 0
fi

# Fallback: Attempt to start the service using the traditional 'service' command
echo "Systemctl failed or Systemd not running. Attempting to start Docker service via 'service' command..."
if sudo service docker start; then
    echo "Docker service started via service command."
    exit 0
else
    echo "ERROR: Failed to start Docker service via both systemctl and service commands." >&2
    echo "Please ensure Systemd is enabled in /etc/wsl.conf and restart the distribution." >&2
    exit 1
fi