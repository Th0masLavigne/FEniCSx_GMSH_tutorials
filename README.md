# FEniCSx_GMSH_tutorials
This repository holds all the documents related to the workshop conducted at I2M Bordeaux in September 2024.

1 mot sur *[créer un alias sur linux](https://www.malekal.com/comment-creer-un-alias-linux/)*




Cette formation a comme objectif d’introduire des outils « open-source » pour la modélisation élément finis. L’ensemble des codes et supports seront laissés en libre accès.


Essayer docker en windows pour équivalent.




Penser à montrer comment trouver les fonctions d’une classe avec ipython3 quand on a du mal à trouver la documentation (exemple booléen GMSH).

# Ressources

## Docker (alternatively singularity)
- *[Docker website](https://www.docker.com/products/docker-desktop/)*

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



# Mail envoyé :
Prérequis :

    Ordinateur (Windows, Linux, Mac) avec Docker installé (si Linux, vérifier les droits super utilisateur) (Docker Desktop: The #1 Containerization Tool for Developers | Docker)
    Installer le logiciel ParaView (Download ParaView)
    Installer le logiciel GMSH (Gmsh: a three-dimensional finite element mesh generator with built-in pre- and post-processing facilities)

Cette formation a pour objectif d’introduire des outils « open-source » pour la modélisation par éléments finis. L’ensemble des codes et supports seront laissés en libre accès.



09/09/2024 matin – Séminaire d’introduction (ouvert à tous)

    Rappel rapide sur les éléments finis.
    Introduction à GMSH/FenicsX pour la modélisation par éléments finis.
    Cas d’application en thermique.
    Cas d’application en mécanique des fluides.
    Cas d’application en mécanique du solide (élasticité, hyperélasticité et poroélasticité).

9 après-midi - 11 septembre – Sessions de formation (inscription nécessaire)

    Génération d’un maillage dans GMSH (à partir d’un schéma, par symétrie de révolution, par opérations booléennes, par importation d’un fichier .stl) et identification des frontières.
    Génération d’un maillage multi-matériaux (tags des volumes).
    Génération d’un maillage à partir d’un fichier .geo.
    Importation de maillage et des tags des frontières dans Fenics.
    Modélisation d’un problème thermique multi-matériaux.

        Thermique stationnaire.
        Thermique transitoire.

    Équation de Stokes.
    Modélisation d’un problème mécanique (élastique (poutre) et hyper-élastique (sein)).
    Modélisation d’un problème de filtration (Darcy).
    Poro-mécanique (Terzaghi axi-symétrique).
