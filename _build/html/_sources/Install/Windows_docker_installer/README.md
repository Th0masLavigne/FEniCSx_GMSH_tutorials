# FEniCSx Environment Setup with Docker Desktop (Windows Direct)

This repository offers a robust and portable solution for running **FEniCSx** simulations directly on Windows using a **Docker container**. This approach ensures a clean, isolated, and consistent environment, eliminating dependency conflicts with the host operating system.

>**IMPORTANT NOTE:** **All `.bat` files in this repository are designed to be run exclusively from the Windows Command Prompt or PowerShell. This setup relies on a correctly installed and running instance of Docker Desktop (using either the Windows native backend or WSL2 integration) to function.**


## Performance and Isolation

This configuration focuses on providing a **reproducible and isolated execution environment**. By leveraging Docker, FEniCSx and all its dependencies are contained within the image, guaranteeing that simulations run consistently regardless of local system changes.

-----

## Directory Structure

The following tree represents the structure of the project folder. The presence of the `src` directory allows for the integration of custom Python libraries into the Docker image build process.

```
Windows_docker_installer/
├── FEniCSx_AGGREGATE.bat               # Run Test.py (or Main.py) if image exists.
├── FEniCSx_docker_-0-9-0_delete_image.bat # Uninstall/delete the custom FEniCSx image.
├── FEniCSx_docker_0-9-0.bat            # Build and verify the FEniCSx image.
├── Interactive_FEniCSx_0-9-0.bat       # Launch an interactive shell in the container.
├── Test.py                             # Example script to run.
└── src
    ├── pyproject.toml                  # Configuration for custom Python packages.
    └── src
        └── MLacour_fenicsx             # Custom FEniCSx library package source.
            ├── Mesh_Bio.py
            ├── Operators.py
            ├── Retrieve_and_plot.py
            └── __init__.py
```

-----

## Key Features

  * **Custom Image Build:** The setup automatically builds a custom Docker image (`fenicsx:v0.9.0`) based on the official FEniCSx base image.
  * **Custom Package Integration:** The contents of the **`src`** folder, including the custom Python library (`MLacour_fenicsx`) defined in **`src/pyproject.toml`**, are automatically integrated and installed into the Docker image during the build process, ensuring all project dependencies are met.
  * **Simplified Launchers:** Provides dedicated batch files for common tasks: building the image, launching an interactive shell, running a specific script, and cleaning up the environment.
  * **Direct Windows Interaction:** Utilizes Docker Desktop's native integration, allowing scripts to be run and file systems to be mounted directly from Windows file paths (`C:\...`).

-----

## Installation

The installation process primarily involves building the custom Docker image.

1.  **Preparation:** Ensure **Docker Desktop** is installed and running on your Windows machine.
2.  **Build Image:** Execute the build script from the repository's root folder:
    ```bash
    FEniCSx_docker_0-9-0.bat
    ```
    This script will:
      * Start Docker Desktop if it is not running.
      * Generate a Dockerfile (if necessary) and build the custom image (`fenicsx:v0.9.0`).
      * Test the image by checking the dolfinx version.
3.  **Completion:** Upon a successful build, the image is ready for use. The script may optionally run the `Test.py` script if it is present.

-----

## Usage

All usage scenarios require **Docker Desktop to be running**. The batch files automatically mount your current working directory (the directory you run the batch file from) to `/home/fenicsx/shared` inside the container.

### 1\. Interactive Shell

Use this to open a bash prompt directly inside the FEniCSx container. This is useful for debugging, installing temporary packages, or using other Linux tools.

```bash
Interactive_FEniCSx_0-9-0.bat
```

### 2\. Run Single Script (Build/Verify)

The build script can also be used to run a test script (`Test.py` or `Main.py`) in the current directory after building the image.

```bash
FEniCSx_docker_0-9-0.bat
```

### 3\. Aggregate Run (Non-Interactive Execution)

Use the aggregate launcher to quickly execute a primary script (e.g., `Main.py` or `Test.py`) non-interactively, assuming the Docker image already exists. This is ideal for production runs.

```bash
FEniCSx_AGGREGATE.bat
```

-----

## Uninstallation

To remove the FEniCSx Docker image and free up disk space:

1.  **Run Uninstaller:** Execute the deletion batch file:
    ```bash
    FEniCSx_docker_-0-9-0_delete_image.bat
    ```
    This script will stop and remove any containers currently running from the image before deleting the image itself (`fenicsx:v0.9.0`).