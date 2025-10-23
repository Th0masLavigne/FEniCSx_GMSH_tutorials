# Installation Guide for FEniCSx v0.9.0

This guide provides detailed instructions to install **FEniCSx v0.9.0** and its dependencies on **Linux**, **Windows (via WSL 2)**, and **macOS**. You can choose between **local installation** or **Docker-based installation** depending on your system and requirements.

## Linux Installation

### Local Installation of FEniCSx v0.9.0

The recommended way to install [FEniCSx](https://fenicsproject.org/download/) on Debian/Ubuntu is via the official APT repository.

> **Note:** These instructions are tested for **Ubuntu 24.04 LTS**. Versions ≥22.04 may work. For older distributions, consider using Docker.

#### Step 1 — Add the FEniCSx Repository
First you need to access the package and update your shell environment:
```sh
sudo add-apt-repository ppa:fenics-packages/fenics
sudo apt update
```

#### Step 2 — Install FEniCSx
Install FEniCSx, with the desired version:
```sh
sudo apt install -y 'fenicsx=2:0.9.0.1*'
```

#### Step 3 — Install Additional Packages
To run the codes you further require to load the following packages:
```sh
sudo apt install python3-pip -y &&\
sudo apt install python3-venv -y &&\
sudo apt update &&\
sudo apt install python3-pandas
```
And 
```sh
sudo apt install xvfb libgl1 libglu1-mesa mesa-utils -y && \
sudo apt install gmsh -y  && \
sudo apt update
```


#### Step 4 — Check Versions
Python module version:
```sh
pip freeze | grep <module-name>
```
APT package version:
```sh
apt show <package-name>
```

Python installed module:
```sh
pip list
```
APT installed module:
```sh
apt list --installed
```

#### Step 5 — Setup Virtual Environment (Recommended for Ubuntu ≥24.04)

Starting from Ubuntu 24.04, it is recommended to work in a **Python virtual environment** for externally managed packages, allowing you to safely use `pip` without affecting system packages.

```sh
cd <your_install_directory>
python3 -m venv fenicsx-env --system-site-packages
```

>**Note:** The --system-site-packages option allows the virtual environment to access system-wide packages, including FEniCSx.

Activate the environment:
```sh
source <your_install_directory>/fenicsx-env/bin/activate
```

Upgrade pip and install additional Python packages:
```sh
pip install --upgrade pip
pip install matplotlib
```

Deactivate the environment:
```sh
deactivate
```

Remove the environment if needed:
```sh
# Delete the current venv
rm -rf fenicsx-env
```

*Remark*: In case of a virtual environment, a quick check to verify that the MPI consistency is:
```sh
mpirun -n 2 python3 -c "from mpi4py import MPI; print(MPI.COMM_WORLD.Get_rank())"
```
You should see ranks 0 and 1 printed — confirming MPI is consistent.


#### Step 6 — Prevent Automatic Package Upgrades

Hold packages:
```sh
sudo apt-mark hold <package-name>
```

List held packages:
```sh
sudo apt-mark showhold
```


Remove hold:
```sh
sudo apt-mark unhold <package-name>
```

Hold all FEniCSx v0.9.0 dependencies:
```sh
sudo apt-mark hold python-ufl-doc python3-ufl python3-pusimp python3-basix \
libdolfinx-real0.9 python3-ffcx libfmt-dev libspdlog-dev fenicsx dolfinx-doc basix-doc \
gcc-13 g++-13 gfortran-13 gcc-13-x86-64-linux-gnu g++-13-x86-64-linux-gnu gfortran-13-x86-64-linux-gnu \
cmake gmsh libdolfinx-dev libdolfinx-real-dev libdolfinx-real0.9 libbasix-dev libbasix0.9 \
libadios2-common-c++11-dev libadios2-common-core-dev libadios2-mpi-auxiliary-2.10 \
libadios2-mpi-auxiliary-dev libadios2-mpi-c++11-2.10 libadios2-mpi-c++11-dev \
libadios2-mpi-core-2.10 libadios2-mpi-core-dev libadios2-mpi-plugins \
libeigen3-dev libpetsc-real-dev libpetsc-real3.19-dev libpetsc-real3.19t64 \
libpetsc3.19-dev-common libpetsc3.19-dev-examples libhypre-2.28.0 libhypre-dev \
libmumps-5.6t64 libmumps-dev libmumps-headers-dev libgraphblas-dev libgraphblas7 \
libparpack2-dev libparpack2t64 libmetis5 python3.12 python3.12-dev python3.12-minimal \
python3.12-venv libpython3.12-dev libpython3.12-stdlib python3-dolfinx python3-dolfinx-real \
trilinos-dev
```

### Docker Engine (in Linux)

If you prefer not to alter system packages, Docker is a robust alternative.
The Directives to install it are provided on the [official website](https://docs.docker.com/engine/install/). 

#### Step 1 — Remove Conflicting Packages
```sh
for pkg in docker.io docker-doc docker-compose docker-compose-v2 podman-docker containerd runc; do sudo apt-get remove $pkg; done
```

#### Step 2 — Add Docker Repository
Set up Docker's apt repository:
```sh
# Add Docker's official GPG key:
sudo apt-get update
sudo apt-get install ca-certificates curl
sudo install -m 0755 -d /etc/apt/keyrings
sudo curl -fsSL https://download.docker.com/linux/ubuntu/gpg -o /etc/apt/keyrings/docker.asc
sudo chmod a+r /etc/apt/keyrings/docker.asc

# Add the repository to Apt sources:
echo \
  "deb [arch=$(dpkg --print-architecture) signed-by=/etc/apt/keyrings/docker.asc] https://download.docker.com/linux/ubuntu \
  $(. /etc/os-release && echo "$VERSION_CODENAME") stable" | \
  sudo tee /etc/apt/sources.list.d/docker.list > /dev/null
sudo apt-get update
```

>**Note**: According to the official website, if you use an Ubuntu derivative distro, such as Linux Mint, you may need to use UBUNTU_CODENAME instead of VERSION_CODENAME.

#### Step 3 — Install Docker Engine
Install the Docker Engine:
```sh
sudo apt-get install docker-ce docker-ce-cli containerd.io docker-buildx-plugin docker-compose-plugin
```

Verify that the Docker Engine installation is successful by running the hello-world image.
```sh
sudo docker run hello-world
```

**Remark**: On linux, you might further need to define a rootless context (if files are locked/secured when the codes are computed with sudo). Please refer to the guide of the [official website](https://docs.docker.com/engine/security/rootless/).

#### Manage Docker as a non-root user:
Create the docker group.
```sh
sudo groupadd docker
```

Add your user to the docker group.
```sh
sudo usermod -aG docker $USER
```

Activate the changes to groups:
```sh
newgrp docker
```

Verify that you can run docker commands without sudo.
```sh
docker run hello-world
```

If you initially ran Docker CLI commands using sudo before adding your user to the docker group, you may see the following error:
```sh
WARNING: Error loading config file: /home/user/.docker/config.json -
stat /home/user/.docker/config.json: permission denied
```

This error indicates that the permission settings for the ~/.docker/ directory are incorrect, due to having used the sudo command earlier.

To fix it:
```sh
sudo chown "$USER":"$USER" /home/"$USER"/.docker -R
sudo chmod g+rwx "$HOME/.docker" -R
```

#### Autostart Docker

From the [official website](https://docs.docker.com/engine/install/linux-postinstall/), to automatically start Docker and containerd on boot for other Linux distributions using systemd, run the following commands:
```sh
sudo systemctl enable docker.service
sudo systemctl enable containerd.service
```

To stop it :
```sh
sudo systemctl disable docker.service
sudo systemctl disable containerd.service
```

In case systemctl is not working with WSL, configure the /etc/wsl.conf file as described in the following [link](https://learn.microsoft.com/fr-fr/windows/wsl/systemd). 

#### Command to start and stop docker manually

If the docker daemon does not start automatically, run:
```sh
sudo systemctl start docker
```

To stop it, run:
```sh
sudo systemctl stop docker
```

## Windows

The use of FEniCSx and Docker is here based on the Windows Subsystem Linux ([WSL 2](https://learn.microsoft.com/fr-fr/windows/wsl/install)). WSL 2 allows to get a Ubuntu distribution within windows (10 or 11)

### Check WSL status and install Ubuntu 24.04 (if not existing)
You can check if WSL is already installed with:
```sh
wsl --status
```

If not installed, after having enabled Microsoft-Windows-Subsystem-Linux and VirtualMachinePlatform features and completed a reboot, launch Windows Powershell or Terminal (admin):
```sh
wsl --install
```

Then ensure using WSL2 with:
```sh
wsl --set-default-version 2
```

Update the version:
```sh
wsl --update
```

You can then list the installed versions distribution using:
```sh
wsl.exe -l -v
```

Available versions can be seen using:
```sh
wsl.exe --list --online
```


If no distribution is currently installed use:
```sh
wsl.exe --install Ubuntu --name <Whatever the name you want to give>
```

To remove a distribution from WSL and delete all of the data associated with that Linux distribution, run:
```sh 
wsl --unregister <distroName>
``` 
where `sh <distroName>` is the name of your Linux distro, which can be seen from the list in the `sh wsl --list --online` command. In the present case, it is Ubuntu-24.04 .

### Install FEniCSx inside WSL
Once Ubuntu is installed you can follow the <ins>here-above installation procedures for Linux</ins> (either local or with Docker) in the WSL 2 terminal. Prefer to use [Docker Desktop](https://desktop.docker.com/win/main/amd64/165256/Docker%20Desktop%20Installer.exe) v4.34.0 or later, please refer to the installation .exe on the official website (please use the following link: [Docker Desktop](https://desktop.docker.com/win/main/amd64/165256/Docker%20Desktop%20Installer.exe) to get the right .exe) and **link it with WSL 2**.



>**Note**: If you have more than a single *default* WSL distro (let's say FEniCSxenv), you need to tell Docker Desktop to integrate with all distros. Step-by-step:

1. Open Docker Desktop (on Windows).

2. o to Settings → Resources → WSL Integration.

3. Under “Enable integration with additional distros,” check the box next to FEniCSxenv

### macOS

For macOS consider using a **conda** environment or install Docker Desktop (with [Apple chip](https://desktop.docker.com/mac/main/arm64/165256/Docker.dmg) or [Intel chip](https://desktop.docker.com/mac/main/amd64/165256/Docker.dmg)). In case of use of docker desktop please refer to the official website for guidelines.

To create the conda environment, run:
```sh
conda create -n fenicsx-env
conda activate fenicsx-env
conda install -c conda-forge fenics-dolfinx mpich pyvista
```
