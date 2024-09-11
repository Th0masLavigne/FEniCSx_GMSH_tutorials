# FEniCSx_GMSH_tutorials

This repository holds all the documents related to the workshop conducted at I2M Bordeaux in September 2024. The objective of the workshop is to introduce open-source softwares for finite element modelling. More specifically, it focuses on the use of FEniCSx (version 0.8.0) and GMSH (version >4.11). Their documentation as well as other softwares are available at the end of this document. 

The following elements are required to be able to run the examples:
- A local installation of [FEniCSx v0.8.0](https://fenicsproject.org/download/) or [Docker Desktop](https://desktop.docker.com/win/main/amd64/165256/Docker%20Desktop%20Installer.exe) v4.34.0 or later *(macOS or linked with WSL2 on Windows)* **or** [Docker Engine](https://docs.docker.com/engine/install/ubuntu/) *(macOS or Linux or Ubuntu WSL2)*,
- [GMSH](https://gmsh.info/#Download) software,
- [Paraview](https://www.paraview.org/download/) software.


This README is organised as follows:
1. Description of the Workshop contents
2. A described description of the installation procedure is recalled for FEniCSx/Docker (in the INSTALL.md file). To install [GMSH](https://gmsh.info/#Download) and [Paraview](https://www.paraview.org/download/), please refer to their official website,
3. The instructions to use a virtual workspace is described (docker commands and aliases definition),
4. All links used as references are provided.

## Contents of the workshop

The workshop indrocudes the following items:
- Brief reminder about the finite element method,
- Creation of a mesh using GMSH:
  * From a sketch,
  * From elementary geometries,
  * Using symmetries,
  * Using boolean operations,
  * From a STL file,
  * with local refinement procedures (from *[Lavigne et al., 2023](https://doi.org/10.1016/j.jmbbm.2023.105902)*),
  * Boundary and domain tagging,
  * Export compatibility: FEniCSx vs Fenics Legacy. 
- Finite Element computation using FEniCSx:
  * Stationary and transient thermal problems,
  * Solid Continous Mechanics problem (elastic, hyper-elastic, penalty contact, updated lagrangian formulation and evaluation of a quantity),
  * Stokes Equation solving in 2D and 3D,
  * Linear and Non-Linear resolutions,
  * Updated mesh resolution,
  * Terzaghi poromechanical model (from *[Lavigne et al., 2023](https://doi.org/10.1016/j.jmbbm.2023.105902)*).

One can also refer to the tutorial from *[Lavigne et al., 2023](https://doi.org/10.1016/j.jmbbm.2023.105902)* ([Github repository](https://github.com/Th0masLavigne/Dolfinx_Porous_Media.git)). Please cite this work if you use codes from this workshop that can be related to this tutorial. One can also refer to the work presented in *[Lavigne et al., 2024]()* for an example of an <ins>**updated Lagrangian**</ins> (*i.e.* [mesh update](https://fenicsproject.discourse.group/t/how-to-do-updated-lagrangian-when-the-displacement-lives-in-a-different-space-to-the-mesh-geometry/10760/2)) poro-elastic model with imposed displacement: [Github repository](https://github.com/Th0masLavigne/Skin_porous_modelling.git). The reaction force is evaluated and volume tags from gmsh are used to map the material parameters with a test case. It completes the multimaterial codes proposed in this workshop.
The Corresponding **Jupyter lab** notebooks are available for an interactive use.

## Installation of the prerequisites

### Linux

#### Local Version of FEniCSx v0.8.0

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

For FEniCSx:
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

## Creating a virtual workspace

### If Docker as a non-root user is configured:
To set an interactive working directory (in an ubuntu environment), respectively using Docker and FEniCSx, the following commands can be used:
```sh
docker run -ti -v $(pwd):/home/fenicsx/shared -w /home/fenicsx/shared th0maslavigne/dolfinx:v0.8.0
```

**Remark**: On a cluster, singularity is often installed instead of docker. In such cases, based on the image.sif file, the command is similar:
```sh
singularity exec /modules/containers/images/dolfinx/dolfinx-0.8.0.sif python3 file.py
```

To create a jupyter container, compute:
```sh
docker run --init -p 8888:8888 -v "$(pwd)":/root/shared --name=jupyter_dolfinx dolfinx/lab:v0.8.0
```

Then to use it, consider using:
```sh
docker container start -i jupyter_dolfinx
```


**Remark:** Docker can store the images and therefore fill a huge amount of space which you can purge with:
```sh
docker stop `docker ps -qa` > /dev/null 2>&1; docker system prune --volumes --all;
```

Sometimes a docker image is missing some python library one want to use. A new image based on an existing image can be created including additionnal packages. For instance, the image used for this workshop have been created using the following Dockerfile:
```sh
FROM dolfinx/dolfinx:v0.8.0
RUN apt update && \
    apt upgrade -y && \
    apt update && \
    apt install libgl1-mesa-glx xvfb -y && \
    python3 -m pip install --upgrade pip && \
    apt update
   
RUN pip3 install pandas \
         imageio \
         pyvista
```

Then to build the image, run in the folder where the 'Dockerfile' is present: 
```sh
docker build .
```
**Remark**: Be careful, it is sensitive to the case so ensure your file is named 'Dockerfile'.

You can list your local images using :
```sh
docker images
```

You can tag the images based on their ID:
```sh
docker tag ImageID meaningful_name
```

To save an image:
```sh
docker save ImageTag > name.tar
```
or
```sh
docker save -o name.tar ImageTag
```

Further commands are described in the cheat sheets available in the resources section.

### Else, Docker as a non-root user is not configured:
All the above commands are working but one need to put `sudo` before.

To set an interactive working directory (in an ubuntu environment), respectively using Docker and FEniCSx, the following commands can be used:
```sh
sudo docker run -ti -v $(pwd):/home/fenicsx/shared -w /home/fenicsx/shared th0maslavigne/dolfinx:v0.8.0
```

To create a jupyter container, compute:
```sh
sudo docker run --init -p 8888:8888 -v "$(pwd)":/root/shared --name=jupyter_dolfinx dolfinx/lab:v0.8.0
```

Then to use it, consider using:
```sh
sudo docker container start -i jupyter_dolfinx
```

And so on. **Remark**: if Docker as a non-root user is not configured, it is required to put sudo as part of the commands for the aliases defined here-after too.

## Create an Alias
The repeated use of a command can be reduced by the use of aliases (see *[create an alias fot linux](https://www.malekal.com/comment-creer-un-alias-linux/)*). 
The idea of the alias is to execute commands from a meaningful name by introducing them in a `sh ~/.bash_aliases`. 

To access this file, run:
```sh
gedit ~/.bash_aliases
```

It will open a window in which you can write your aliases on the following basis:
```sh
alias <Meaningful_name>='<the command>'
```

For example, in the present workshop, one could create:
```sh
alias fenicsx_v0_8_0='docker run -ti -v $(pwd):/home/fenicsx/shared -w /home/fenicsx/shared th0maslavigne/dolfinx:v0.8.0'
```

This is of interest and will allow you to have different versions of a same environment without conflicting package:
```sh
alias fenicsx_v0_6_0='docker run -ti -v $(pwd):/home/fenicsx/shared -w /home/fenicsx/shared dolfinx/dolfinx:v0.6.0'
alias fenicsx_v0_7_3='docker run -ti -v $(pwd):/home/fenicsx/shared -w /home/fenicsx/shared dolfinx/dolfinx:v0.7.3'
alias fenicsx_v0_8_0='docker run -ti -v $(pwd):/home/fenicsx/shared -w /home/fenicsx/shared dolfinx/dolfinx:v0.8.0'
```

If you prefer using the jupyter lab notebooks, the following commands can be created:
Run once:
```sh
alias create_fenicsx_v0_8_0_jupyter='docker run --init -p 8888:8888 -v "$(pwd)":/root/shared --name=jupyter_dolfinx dolfinx/lab:v0.8.0'
```

Then use:
```sh
alias run_fenicsx_v0_8_0_jupyter='docker container start -i jupyter_dolfinx'
alias stop_fenicsx_v0_8_0_jupyter='docker stop jupyter_dolfinx'
alias logs_fenicsx_v0_8_0_jupyter='docker logs jupyter_dolfinx'
```

**Remark:** Including the bash term at the end allows to exit the python environnment to the linux command. Here are two examples respectively for fenics legacy and pymesh:
```sh
alias fenics2019='docker run -ti -v $(pwd):/home/fenics/shared -w /home/fenics/shared pymor/fenics_py3.9 bash'
alias pymesh='docker run -ti -v $(pwd):/home/pymesh/shared -w /home/pymesh/shared pymesh/pymesh bash'
```

## Acknowledgments 
This repository is inspired from the work of [Jørgen S. Dokken](https://jsdokken.com/tutorials.html) and [Christophe Geuzaine](https://gitlab.onelab.info/gmsh/gmsh/tree/master/tutorials). The author also wants to thank F. Daghia and Ludovic Chamoin for the quality of the finite element courses at the \'Ecole Normale Supérieure Paris-Saclay, Jack Hale and Stéphane Urcun for their help in debugging the codes throughout the PhD work and Giuseppe Sciumè for the invite.

This activity is part of Thomas Lavigne PhD work. This research was funded in whole, or in part, by the Luxembourg National Research Fund (FNR), grant reference No. 17013182. For the purpose of open access, the author has applied a Creative Commons Attribution 4.0 International (CC BY 4.0) license to any Author Accepted Manuscript version arising from this submission.

## Resources

**Remark:** If you do not find the documentation for a specific item, you can use ipython3 and the help() command. For example if you have a class object element that we call mesh, executing "mesh." + "tab" you will be able to navigate the attributes to the object and then apply the help command.

### Good Practice for coding conferences:
- *[Clean Code - Uncle Bob / Lesson 1](https://www.youtube.com/watch?v=7EmboKQH8lM)*
- *[Clean Code - Conference](https://www.youtube.com/watch?v=7EmboKQH8lM&list=PLmmYSbUCWJ4x1GO839azG_BBw8rkh-zOj&index=1)*

### Docker (alternatively singularity)
- *[Docker website](https://www.docker.com/products/docker-desktop/)*
- *[Docker Installation Ressource](https://docs.docker.com/desktop/install/windows-install/)*
- *[Docker hub](https://hub.docker.com/)*
- *[Docker cheat sheet](https://docs.docker.com/get-started/docker_cheatsheet.pdf)*
- *[Docker cheat sheet 2](https://dockerlabs.collabnix.com/docker/cheatsheet/)*
- *[Windows docker cheat sheet](https://gist.github.com/danijeljw/a7a2553bd06742648172363ce3983a9a)*
- *[Singularity](https://docs.sylabs.io/guides/3.5/user-guide/introduction.html)*

### Paraview
- *[Paraview Download](https://www.paraview.org/download/)*

### GMSH 
- *[GMSH Download](https://gmsh.info/)*
- *[GMSH Introductive Presentation](https://gmsh.info/doc/course/general_overview.pdf)*
- *[GMSH Manual](https://gmsh.info/doc/texinfo/gmsh.html)*
- *[GMSH GitLab](https://gitlab.onelab.info/gmsh/gmsh)*
- *[GMSH API tutoraials](https://bthierry.pages.math.cnrs.fr/tutorial/gmsh/api/)*
- *[OpenCascade Commands](https://koehlerson.github.io/gmsh.jl/dev/occ/occ/)*

### FEniCSx
- *[FEniCS Website](https://fenicsproject.org/)*
- *[FEniCSx Tutorial](https://jsdokken.com/dolfinx-tutorial/)*
- *[FEniCSx Changelog](https://github.com/FEniCS/dolfinx/releases)*
- *[FEniCSx Discourse Forum](https://fenicsproject.discourse.group/)*
- *[FEniCSx Github](https://github.com/orgs/FEniCS/repositories)*
- *[Dolfinx Documentation](https://docs.fenicsproject.org/dolfinx/v0.8.0/python/)*
- *[Basix Documentation](https://docs.fenicsproject.org/basix/v0.8.0/python/)*
- *[UFL Documentation](https://fenics.readthedocs.io/projects/ufl/en/latest/)*
- *[FFCx Documentation](https://docs.fenicsproject.org/ffcx/main/)*
- *[FEniCS book](https://launchpadlibrarian.net/83776282/fenics-book-2011-10-27-final.pdf)*
- *[FEniCS Legacy Documentation](https://fenicsproject.org/olddocs/)*
- *[Mesh Update](https://fenicsproject.discourse.group/t/how-to-do-updated-lagrangian-when-the-displacement-lives-in-a-different-space-to-the-mesh-geometry/10760/2)*

**Remark:** On ubuntu jammy, the FEniCSx version 0.8.0 had a conflict on my laptop requiring to remove python3-numba:
```sh
sudo apt remove python3-numba
```
### GIT
- *[GIT Docs](https://docs.github.com/en/get-started)*
- *[GIT Reference](https://git-scm.com/docs)*
- *[GIT interactive tutorial](https://learngitbranching.js.org/?locale=fr_FR)*

### Others
#### Deal.ii
- *[Deal.ii website](https://www.dealii.org/)*
- *[Deal.ii tutorials](https://www.dealii.org/current/doxygen/deal.II/Tutorial.html)*
- *[Deal.ii documentation](https://www.dealii.org/current/index.html)*
- *[Deal.ii library](https://www.dealii.org/current/doxygen/deal.II/index.html)*

#### Doxygen
- *[Doxygen](https://www.doxygen.nl/index.html)*

#### Neper
- *[Neper](https://neper.info/)*

#### Pygmsh
- *[Pygmsh](https://pypi.org/project/pygmsh/)*
- *[Pygmsh github](https://github.com/nschloe/pygmsh)*

#### Meshio
- *[Meshio](https://pypi.org/project/meshio/)*
- *[Meshio github](https://github.com/nschloe/meshio)*

#### Pymesh
- *[Pymesh](https://pymesh.readthedocs.io/en/latest/)*
- *[Pymesh github](https://github.com/PyMesh/PyMesh)*

#### Netgen/NGSolve
- *[Website](https://ngsolve.org/)*
- *[Tutorials](https://jschoeberl.github.io/iFEM/intro.html)*

#### Boolean Operations for STL Files
OpenCascade Kernel does not support stls for boolean operations.
There a few libraries out there that support boolean operations for meshes, you might want to give them a try. Here's a personal list I have of mostly C++ and Python packages from a comment of [@gnikit](https://gitlab.onelab.info/gmsh/gmsh/-/issues/2493):
- [CGAL](https://github.com/CGAL/cgal)
- [VTK-based improved booleans](https://github.com/zippy84/vtkbool)
- [MeshLib](https://github.com/MeshInspector/MeshLib)
- [mcut](https://github.com/cutdigital/mcut)
- [trimesh](https://github.com/mikedh/trimesh)
- [vedo](https://github.com/marcomusy/vedo/blob/master/examples/basic/boolean.py)
- [PyMeshLab](https://github.com/cnr-isti-vclab/PyMeshLab) (and its full application [MeshLab](https://www.meshlab.net/))
- [PyMesh](https://github.com/PyMesh/PyMesh)
- [Pycork](https://pypi.org/project/pycork/) (and it's C base library [cork](https://github.com/gilbo/cork))
- [Blender](https://www.blender.org/) (also has a python module)


#### Tomography to conform mutlipart meshes
- *[Tomo2FE github](https://github.com/ANR-MultiFIRE/TomoToFE/blob/main/workflow2/Workflow2-Python.ipynb)*
- *[Tomo2FE article](https://letters.rilem.net/index.php/rilem/article/view/184)*

#### Function to export lists in a csv file
```python
    def export_to_csv(data, filename, header=None):
        import csv
        try:
            with open(filename, 'w', newline='') as file:
                writer = csv.writer(file)
                if header:
                    writer.writerow(header)
                writer.writerows(data)
            print(f"Data exported to {filename} successfully")
        except Exception as e:
            print(f"An error occurred while exporting data to {filename}: {e}")
```