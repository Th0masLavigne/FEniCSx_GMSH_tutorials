# FEniCSx_GMSH_tutorials

This repository holds all the documents related to the workshop conducted at I2M Bordeaux in September 2024. The objective of the workshop is to introduce open-source softwares for finite element modelling. More specifically, it focuses on the use of FEniCSx (version 0.8.0) and GMSH (version >4.11). Their documentation as well as other softwares are available at the end of this document. 

The following elements are required to be able to run the examples:
- A **local installation** of FEniCSx v0.8.0 or [Docker Desktop](https://docs.docker.com/desktop/install/windows-install/) v4.34.0 or later *(macOS or linked with WSL2 on Windows)* **or** [Docker Engine](https://docs.docker.com/engine/install/ubuntu/) *(macOS or Linux or Ubuntu WSL2)*,
- GMSH software,
- Paraview software.

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
The Corresponding **jupyter** notebooks are available for an interactive use.

## Installation of the prerequisites

### Linux

#### Local Version of FEniCSx v0.8.0

From the [official website](https://fenicsproject.org/download/), the easiest way to install FEniCSx on Debian or Ubuntu Linux is via apt.

**Remark**: FEniCSx v0.8.0 requires at least **Ubuntu 22.04**. If you are using **Ubuntu 20.04**, consider using a docker. Previous distro versions might require an upgrade.

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


### Docker Engine

If you already have a version of FEniCSx you prefer not to upgrade, it is possible to use Docker (alternatively singularity for cluster computations).
The Directives to install it are provided on the (official website)[https://docs.docker.com/engine/install/]. 

First, uninstall conflicting packages:
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

**Remark**: On linux, you might further need to define a rootless context (if files are locked/secured when the codes are computed with sudo). Please refer to the guide of the (official website)[https://docs.docker.com/engine/security/rootless/].

#### Manage Docker as a non-root user:
Create the docker group.
```sh
sudo groupadd docker


Add your user to the docker group.
```sh
sudo usermod -aG docker $USER
```

Activate the changes to groups:
```sh
newgrp docker

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

From the (official website)[https://docs.docker.com/engine/install/linux-postinstall/], to automatically start Docker and containerd on boot for other Linux distributions using systemd, run the following commands:
```sh
sudo systemctl enable docker.service
sudo systemctl enable containerd.service
```

To stop it :
```sh
sudo systemctl disable docker.service
sudo systemctl disable containerd.service
```

### Windows












### macOS

Use of conda or binaries (docker desktop / docker engine)
















## Creating a virtual workspace

To set an interactive working directory (in an ubuntu environment), respectively using Docker and FEniCSx, the following commands can be used:
```sh
docker run -ti -v $(pwd):/home/fenicsx/shared -w /home/fenicsx/shared dolfinx/dolfinx:v0.8.0
```

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

The repeated use of a command can be reduced by the use of aliases (see *[create an alias fot linux](https://www.malekal.com/comment-creer-un-alias-linux/)*). Several containers can be considered based on the version you need:

```sh
alias fenicsx_v0_5_2='docker run -ti -v $(pwd):/home/fenicsx/shared -w /home/fenicsx/shared th0maslavigne/dolfinx:v0.5.2'
alias fenicsx_v0_6_0='docker run -ti -v $(pwd):/home/fenicsx/shared -w /home/fenicsx/shared dolfinx/dolfinx:v0.6.0'
alias fenicsx_v0_7_3='docker run -ti -v $(pwd):/home/fenicsx/shared -w /home/fenicsx/shared dolfinx/dolfinx:v0.7.3'
alias fenicsx_v0_8_0='docker run -ti -v $(pwd):/home/fenicsx/shared -w /home/fenicsx/shared dolfinx/dolfinx:v0.8.0'
alias create_fenicsx_v0_8_0_jupyter='docker run --init -p 8888:8888 -v "$(pwd)":/root/shared --name=jupyter_dolfinx dolfinx/lab:v0.8.0'
alias run_fenicsx_v0_8_0_jupyter='docker container start -i jupyter_dolfinx'
alias stop_fenicsx_v0_8_0_jupyter='docker stop jupyter_dolfinx'
alias logs_fenicsx_v0_8_0_jupyter='docker logs jupyter_dolfinx'
alias fenics2019='docker run -ti -v $(pwd):/home/fenics/shared -w /home/fenics/shared pymor/fenics_py3.9 bash'
alias pymesh='docker run -ti -v $(pwd):/home/pymesh/shared -w /home/pymesh/shared pymesh/pymesh bash'
```
**Remark:** Including the bash term at the end allows to exit the python environnment to the linux command.

**Remark:** Docker can store the images and therefore fill a huge amount of space which you can purge with:
```sh
alias dockerRemoveAll="docker stop `docker ps -qa` > /dev/null 2>&1; docker system prune --volumes --all;"
```

Sametimes a docker image is missing some python library you'd need. You can create a new image (with a dockerfile that you will build) based on an existing image. For instance, you want the dolfinx image with pandas library. Your Dockerfile will contain:
```sh
FROM dolfinx/dolfinx:v0.5.2
RUN pip3 install pandas
```
Then to build the image, run in the folder: 
```sh
docker build .
```

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

**Remark:** On ubuntu jammy, the FEniCSx version 0.8.0 had a conflict requiring to remove python3-numba:
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

### Boolean Operations for STL Files
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


### Tomography to conform mutlipart meshes
- *[Tomo2FE github](https://github.com/ANR-MultiFIRE/TomoToFE/blob/main/workflow2/Workflow2-Python.ipynb)*
- *[Tomo2FE article](https://letters.rilem.net/index.php/rilem/article/view/184)*

### Function to export lists in a csv file
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

