FROM dolfinx/dolfinx:v0.8.0
RUN apt update && \
    apt upgrade -y && \
    apt update && \
    apt install libgl1-mesa-glx xvfb && \
    python3 -m pip install --upgrade pip && \
    apt update
   
RUN pip3 install pandas \
		 imageio \
		 pyvista
