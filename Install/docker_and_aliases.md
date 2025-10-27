## Creating a virtual workspace

### If Docker as a non-root user is configured:
To set an interactive working directory (in an ubuntu environment), respectively using Docker and FEniCSx, the following commands can be used:
```sh
docker run -ti -v $(pwd):/home/fenics/shared -w /home/fenics/shared th0maslavigne/dolfinx:v0.9.0
```

>**Note**: If you have several WSL Distros, enable them all in Docker Desktop (if used, in Windows). See Install.md

>**Remark**: On a cluster, singularity is often installed instead of docker. In such cases, based on the image.sif file, the command is similar:
```sh
singularity exec /modules/containers/images/dolfinx/dolfinx-0.9.0.sif python3 file.py
```

To create a jupyter container, compute:
```sh
docker run --init -p 8888:8888 -v "$(pwd)":/root/shared --name=jupyter_dolfinx dolfinx/lab:v0.9.0
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
FROM dolfinx/dolfinx:v0.9.0
RUN apt update && \
    apt upgrade -y && \
    apt update && \
    apt install xvfb libgl1-mesa-dri libglx-mesa0 libglu1-mesa mesa-utils -y && \
    python3 -m pip install --upgrade pip && \
    apt update 
   
RUN pip3 install pandas \
         imageio \
         pyvista
```

Then to build the image, run in the folder where the 'Dockerfile' is present: 
```sh
docker build -f Dockerfile -t FEniCSx:v0.9.0 .
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
sudo docker run -ti -v $(pwd):/home/fenicsx/shared -w /home/fenicsx/shared th0maslavigne/dolfinx:v0.9.0
```

To create a jupyter container, compute:
```sh
sudo docker run --init -p 8888:8888 -v "$(pwd)":/root/shared --name=jupyter_dolfinx dolfinx/lab:v0.9.0
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
alias fenicsx_v0_9_0='docker run -ti -v $(pwd):/home/fenicsx/shared -w /home/fenicsx/shared th0maslavigne/dolfinx:v0.9.0'
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
alias create_fenicsx_v0_8_0_jupyter='docker run --init -p 8888:8888 -v "$(pwd)":/root/shared --name=jupyter_dolfinx08 dolfinx/lab:v0.8.0'
```
or
```sh
alias create_fenicsx_v0_9_0_jupyter='docker run --init -p 8888:8888 -v "$(pwd)":/root/shared --name=jupyter_dolfinx09 dolfinx/lab:v0.9.0'
```

Then use:
```sh
alias run_fenicsx_v0_8_0_jupyter='docker container start -i jupyter_dolfinx08'
alias stop_fenicsx_v0_8_0_jupyter='docker stop jupyter_dolfinx08'
alias logs_fenicsx_v0_8_0_jupyter='docker logs jupyter_dolfinx08'
```
or
```sh
alias run_fenicsx_v0_9_0_jupyter='docker container start -i jupyter_dolfinx09'
alias stop_fenicsx_v0_9_0_jupyter='docker stop jupyter_dolfinx09'
alias logs_fenicsx_v0_9_0_jupyter='docker logs jupyter_dolfinx09'
```

**Remark:** Including the bash term at the end allows to exit the python environnment to the linux command. Here are two examples respectively for fenics legacy and pymesh:
```sh
alias fenics2019='docker run -ti -v $(pwd):/home/fenics/shared -w /home/fenics/shared pymor/fenics_py3.9 bash'
alias pymesh='docker run -ti -v $(pwd):/home/pymesh/shared -w /home/pymesh/shared pymesh/pymesh bash'
```