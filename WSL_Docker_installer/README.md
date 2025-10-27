-----

>**IMPORTANT NOTE FOR INSTALLATION:** **The `install_wsl2_fenicsx_docker.bat` file, along with its dependency folder (`src_install/`) and all other files from the root of this repository (the entire WSL\_Docker\_installer folder), must be run from a dedicated Windows folder containing ONLY the necessary installation files. This is crucial because the entire contents of this folder will be copied to your Linux home directory within WSL2 (`~/Softwares/WSL2_Docker_FEniCSx`) as the installation source. After a successful installation, you can safely delete this original Windows folder, as a safe version is stored within the WSL distribution. The other launcher batch files (`FEniCSx_interactive.bat`, `FEniCSx_command.bat`, `uninstall_wsl2_fenicsx_docker.bat`) can be placed anywhere else on your Windows system for convenience.**

-----

# 🐳 FEniCSx Environment Setup with Docker and WSL2

This repository provides an automated, professional, and highly efficient setup for running **FEniCSx** simulations. It utilizes a **Docker container** orchestrated within the **Windows Subsystem for Linux 2 (WSL2)**, delivering a consistent, isolated, and high-performance simulation environment on Windows systems.

## Performance and Efficiency Statement

For local development and execution, the use of **Docker within WSL2 (Windows 11)** is often sufficient and demonstrably efficient. This methodology has specifically proven to be more performant than running Docker directly on Windows (via Docker Desktop's native integration) or managing a local installation with a Python 3 virtual environment (`venv`) of FEniCSx. This efficiency is achieved through WSL2's optimized I/O performance and direct integration with the Linux kernel.

-----

## Directory Structure

The following tree represents the complete structure of the installation folder (`WSL2_Docker_FEniCSx`) that is copied into the WSL distribution's home directory during the installation process:

```
WSL2_Docker_FEniCSx/
├── FEniCSx_command.bat
├── FEniCSx_interactive.bat
├── README.md
├── install_wsl2_fenicsx_docker.bat
├── report.log
├── report_uninstall.log
├── src_install
│   ├── apt_environment_aliases
│   │   ├── fenicsx_aliases.sh
│   │   └── verification.sh
│   ├── custom_packages
│   │   ├── package.txt
│   │   ├── pyproject.toml
│   │   └── src
│   │       └── MLacour_fenicsx
│   │           ├── Mesh_Bio.py
│   │           ├── Operators.py
│   │           ├── Retrieve_and_plot.py
│   │           └── __init__.py
│   └── docker
│       ├── docker_fenicsx_env.sh
│       └── docker_install_check.sh
├── test_codes
│   └── Test.py
└── uninstall_wsl2_fenicsx_docker.bat
```

-----

## Key Features

  * **Isolated and Consistent:** FEniCSx runs inside a secure Docker container, ensuring environmental stability and reproducibility.
  * **Optimized Performance:** Leverages the performance benefits and native Linux compatibility of **WSL2**.
  * **Simple Image Update:** The Docker image can be updated or rebuilt (e.g., after custom modifications) by simply **re-running `install_wsl2_fenicsx_docker.bat`** (respecting the insitial Note about file dependencies). The installer will update the image without requiring a full reinstallation of the WSL distribution or re-engaging the interactive user prompts.
  * **Extensible Environment:** Other Python `pip` packages can be added to the Docker image by modifying the Docker build logic within the **`src_install/docker/docker_fenicsx_env.sh`** file and placing any necessary configuration (like `pyproject.toml`) inside the **`src_install/custom_packages`** directory.
  * **Windows-Native Launchers:** All **`.bat` files must be run from Windows** (Command Prompt or PowerShell) to interact with and launch the environment inside WSL2.

-----

## Installation

1.  **Execution:** Run the main installer script from the dedicated Windows folder by double clicking on it:
    ```bash
    install_wsl2_fenicsx_docker.bat
    ```
2.  **Process:** The script will guide you through:
      * Naming your WSL distribution.
      * Creating a WSL **username and password**. You will need to type `exit` to pursue the installation once this step is complete.
      * The necessary steps to install Docker and build the FEniCSx image.
3.  **Completion:** Upon success, a `report.log` file is generated, and the FEniCSx environment is ready for use.

-----

## Usage

### 1\. Launching from Windows (Using .bat files)

All launcher batch files (`.bat`) are designed to be run directly from **Windows**.

  * **Interactive Shell:** Use `FEniCSx_interactive.bat` to open a command line session **inside the FEniCSx Docker container**. You will be prompted for the project folder path to mount.
  * **Run Script/Command:** Use `FEniCSx_command.bat` to execute a single Python script or command **non-interactively** within the container. (You must edit the batch file to specify the target script path).

### 2\. Direct Use in WSL2 Command Line (Using Aliases)

Once your WSL distribution is running, you can open a shell directly and utilize the custom bash aliases created during installation.

1.  **Access WSL:** Open the WSL distribution command line (e.g., by typing `wsl` in Windows search or opening Windows Terminal).

2.  **Navigate:** Change directory to the folder containing your FEniCSx code (e.g., `cd ./Documents/MyProject`).

3.  **Run Interactive Session:** To open an interactive shell inside the FEniCSx Docker environment:

    ```bash
    fenicsx
    ```

    *This command mounts your current WSL directory into the container at `/home/fenicsx/shared`.*

4.  **Run Non-Interactive Script:** To execute a script using Python or MPI without entering the interactive shell:

    ```bash
    # Run a simple Python script
    fenicsx-run "python3 filename.py"

    # Run a script using mpirun for parallel execution
    fenicsx-run "mpirun -np 4 python3 filename.py"
    ```

    *The `fenicsx-run` alias passes the entire quoted string to a non-interactive bash session within the container.*

-----

## Uninstallation

To completely remove the installed FEniCSx environment and all associated files:

1.  **Run Uninstaller:** Execute the uninstallation batch file from Windows:
    ```bash
    uninstall_wsl2_fenicsx_docker.bat
    ```
2.  **Warning:** This is an **irreversible** action and will delete the entire `FEniCSxenv` WSL distribution, including all user files and the copied installation source code.