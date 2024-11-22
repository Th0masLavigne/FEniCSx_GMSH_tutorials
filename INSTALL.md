## Installation of the prerequisites

### Linux

#### Local Version of FEniCSx v0.8.0 or v0.9.0

From the [official website](https://fenicsproject.org/download/), the easiest way to install FEniCSx on Debian or Ubuntu Linux is via apt.

**Remark**: FEniCSx v0.8.0 requires at least **Ubuntu 22.04**. If you are using **Ubuntu 20.04**, consider using a docker. Previous distro versions might require an upgrade.

#### Install FEniCSx
First you need to access the package:
```sh
sudo add-apt-repository ppa:fenics-packages/fenics
```

Update the environment:
```sh
sudo apt update
```

Install FEniCSx:
```sh
sudo apt install fenicsx
```

#### Additional packages
To run the codes you further require to load the following packages:
```sh
pip3 install h5py \
             imageio \
             pyvista \
             gmsh
```

```sh
sudo apt install libgl1-mesa-glx xvfb -y && \
sudo apt install gmsh -y  && \
sudo apt update
```

If one further need to check a python module version, one can use:
```sh
pip3 freeze | grep module
```


#### Avoid an upgrade of the packages

To avoid an upgrade of the packages you can set them 'on hold'.
To hold a package:
```sh
sudo apt-mark hold <package-name>
```

To show all packages on hold:
```sh
sudo apt-mark showhold
```


To remove the hold:
```sh
sudo apt-mark unhold <package-name>
```

For FEniCSx 0.8.0:
```sh

sudo apt-mark hold  basix-doc \
                    dolfinx-doc \
                    fenicsx \
                    libbasix-dev \
                    libbasix0.8 \
                    libdolfinx-dev \
                    libdolfinx-real-dev \
                    libdolfinx-real0.8 \
                    python-petsc4py-doc \
                    python-ufl-doc \
                    python3-basix \
                    python3-dolfinx \
                    python3-dolfinx-real \
                    python3-ffcx \
                    python3-petsc4py \
                    python3-petsc4py-real \
                    python3-petsc4py-real3.15 \
                    python3-ufl
```

For FEniCSx 0.9.0 :
```sh

sudo apt-mark hold basix-doc \
                   libbasix0.9 \
                   libfmt8 \
                   python-ufl-doc \
                   libspdlog1 \
                   libbasix-dev \
                   python3-pusimp \
                   python3-basix \
                   libdolfinx-real0.9 \
                   python3-dolfin \
                   python3-ffcx \
                   libdolfinx-real-dev \
                   libfmt-dev \
                   dolfin-bin \
                   libspdlog-dev \
                   libdolfinx-dev \
                   python3-dolfinx-real \
                   python3-dolfinx fenicsx
```

### Docker Engine

If you already have a version of FEniCSx you prefer not to upgrade, it is possible to use Docker (alternatively singularity for cluster computations).
The Directives to install it are provided on the [official website](https://docs.docker.com/engine/install/). 

First, uninstall potential conflicting packages:
```sh
for pkg in docker.io docker-doc docker-compose docker-compose-v2 podman-docker containerd runc; do sudo apt-get remove $pkg; done
```

#### Install using the apt repository. 
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

**Remark**: According to the official website, if you use an Ubuntu derivative distro, such as Linux Mint, you may need to use UBUNTU_CODENAME instead of VERSION_CODENAME.

Install the Docker packages:
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

#### Configure Docker to start on boot with systemd

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

#### Command to start and stop docker

If the docker daemon does not start automatically, run:
```sh
sudo systemctl start docker
```

To stop it, run:
```sh
sudo systemctl stop docker
```

### Windows

The use of FEniCSx and Docker is here based on the Windows Subsystem Linux ([WSL 2](https://learn.microsoft.com/fr-fr/windows/wsl/install)). WSL 2 allows to get a Ubuntu distribution within windows (10 or 11)

#### Check WSL status and install Ubuntu 22.04 (if not existing)
In Windows Powershell, run:
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

You can then list the available distribution using:
```sh
wsl.exe -l -v
```

If no distribution is currently installed use:
```sh
wsl.exe --install Ubuntu
```

To remove a distribution from WSL and delete all of the data associated with that Linux distribution, run:
```sh 
wsl --unregister <distroName>
``` 
where `sh <distroName>` is the name of your Linux distro, which can be seen from the list in the `sh wsl -l` command.

#### Install FEniCSx
Once Ubuntu is installed you can follow the <ins>here-above installation procedures for Linux</ins> (either local or with Docker) in the WSL 2 terminal. The non-root user procedure in Docker does not look necessary within WSL 2. Ensure adding sudo at the beginning of your commands/aliases in that case. 

If one prefer to use [Docker Desktop](https://desktop.docker.com/win/main/amd64/165256/Docker%20Desktop%20Installer.exe) v4.34.0 or later, please refer to the installation .exe on the official website (please use the following link: [Docker Desktop](https://desktop.docker.com/win/main/amd64/165256/Docker%20Desktop%20Installer.exe) to get the right .exe) and link it with WSL 2.

### macOS

For macOS consider using a **conda** environment or install Docker Desktop (with [Apple chip](https://desktop.docker.com/mac/main/arm64/165256/Docker.dmg) or [Intel chip](https://desktop.docker.com/mac/main/amd64/165256/Docker.dmg)). In case of use of docker desktop please refer to the official website for guidelines.

To create the conda environment, run:
```sh
conda create -n fenicsx-env
conda activate fenicsx-env
conda install -c conda-forge fenics-dolfinx mpich pyvista
```
