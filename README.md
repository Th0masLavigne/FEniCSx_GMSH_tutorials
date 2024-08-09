# FEniCSx_GMSH_tutorials
This repository holds all the documents related to the workshop conducted at I2M Bordeaux in September 2024. The objective of the workshop is to introduce open-source softwares for finite element modelling. More specifically, it focuses on the use of FEniCSx and GMSH. Their documentation as well as other softwares are available at the end of this document. 

The following elements are required to be able to run the examples:
- Docker or Singularity with super-user rights (or a local installation of the softwares),
- GMSH software,
- Paraview software.

To set an interactive working directory, respectively using Docker and FEniCSx, the following commands can be used:
```cmd
docker run -ti -v $(pwd):/home/fenicsx/shared -w /home/fenicsx/shared dolfinx/dolfinx:v0.8.0
```

```cmd
singularity exec /modules/containers/images/dolfinx/dolfinx-0.8.0.sif python3 file.py
```

The repeated use of a command can be reduced by the use of aliases (see *[create an alias fot linux](https://www.malekal.com/comment-creer-un-alias-linux/)*). Several containers can be considered based on the version you need:

```cmd
alias fenicsx_v0_5_2='docker run -ti -v $(pwd):/home/fenicsx/shared -w /home/fenicsx/shared th0maslavigne/dolfinx:v0.5.2'
alias fenicsx_v0_6_0='docker run -ti -v $(pwd):/home/fenicsx/shared -w /home/fenicsx/shared dolfinx/dolfinx:v0.6.0'
alias fenicsx_v0_7_3='docker run -ti -v $(pwd):/home/fenicsx/shared -w /home/fenicsx/shared dolfinx/dolfinx:v0.7.3'
alias fenicsx_v0_8_0='docker run -ti -v $(pwd):/home/fenicsx/shared -w /home/fenicsx/shared dolfinx/dolfinx:v0.8.0'
alias fenics2019='docker run -ti -v $(pwd):/home/fenics/shared -w /home/fenics/shared pymor/fenics_py3.9 bash'
alias pymesh='docker run -ti -v $(pwd):/home/pymesh/shared -w /home/pymesh/shared pymesh/pymesh bash'
```
**Remark:** Including the bash term at the end allows to exit the python environnment to the linux command.
**Remark:** Docker can store the images and therefore fill a huge amount of space which you can purge with:
```cmd
alias dockerRemoveAll="docker stop `docker ps -qa` > /dev/null 2>&1; docker system prune --volumes --all;"
```

Sametimes a docker image is missing some python library you'd need. You can create a new image (with a dockerfile that you will build) based on an existing image. For instance, you want the dolfinx image with pandas library. Your dockerfile will contain:
```cmd
FROM dolfinx/dolfinx:v0.5.2
RUN pip3 install pandas
```


** hyper el poutre dans local finger tutorial. Ajout avec conditionnel contact**



1 mot sur 

Reprendre les commandes tuto docker et git sur 


Cette formation a comme objectif d’introduire des outils « open-source » pour la modélisation élément finis. L’ensemble des codes et supports seront laissés en libre accès.


Essayer docker en windows pour équivalent.

**Reprendre dans ipad mon tuto docker et git**


Penser à montrer comment trouver les fonctions d’une classe avec ipython3 quand on a du mal à trouver la documentation (exemple booléen GMSH).

# Ressources

## Docker (alternatively singularity)
- *[Docker website](https://www.docker.com/products/docker-desktop/)**
- *[Docker hub](https://hub.docker.com/)*
- *[Docker cheat sheet](https://docs.docker.com/get-started/docker_cheatsheet.pdf)
- *[Docker cheat sheet 2](https://dockerlabs.collabnix.com/docker/cheatsheet/)*
- *[Windows docker cheat sheet](https://gist.github.com/danijeljw/a7a2553bd06742648172363ce3983a9a)*

## Paraview
- *[Paraview Download](https://www.paraview.org/download/)*

## GMSH 
- *[GMSH Download](https://gmsh.info/)*
- *[GMSH Introductive Presentation]([https://duckduckgo.com](https://gmsh.info/doc/course/general_overview.pdf))*
- *[GMSH Manual](https://gmsh.info/doc/texinfo/gmsh.html)*
- *[GMSH GitLab](https://gitlab.onelab.info/gmsh/gmsh)*
- *[GMSH API tutoraials](https://bthierry.pages.math.cnrs.fr/tutorial/gmsh/api/)*
- *[OpenCascade Commands](https://koehlerson.github.io/gmsh.jl/dev/occ/occ/)*

## FEniCSx
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

**Remark:** On ubuntu jammy, the FEniCSx version 0.8.0 had a conflict requiring to remove python3-numba:
```cmd
sudo apt remove python3-numba
```
## GIT
- *[GIT Reference](https://git-scm.com/docs)*
- *[GIT interactive tutorial](https://learngitbranching.js.org/?locale=fr_FR)*

## Others
### Deal.ii
- *[Deal.ii website](https://www.dealii.org/)*
- *[Deal.ii tutorials](https://www.dealii.org/current/doxygen/deal.II/Tutorial.html)*
- *[Deal.ii documentation](https://www.dealii.org/current/index.html)*
- *[Deal.ii library](https://www.dealii.org/current/doxygen/deal.II/index.html)*

### Neper
- *[Neper](https://neper.info/)*

### Pygmsh
- *[Pygmsh](https://pypi.org/project/pygmsh/)*
- *[Pygmsh github](https://github.com/nschloe/pygmsh)*

### Meshio
- *[Meshio](https://pypi.org/project/meshio/)*
- *[Meshio github](https://github.com/nschloe/meshio)*

### Pymesh
- *[Pymesh](https://pymesh.readthedocs.io/en/latest/)*
- *[Pymesh githbu](https://github.com/PyMesh/PyMesh)*
