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