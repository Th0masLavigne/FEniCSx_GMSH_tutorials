# FEniCSx Local Environment Setup with WSL2

This repository provides a simplified, automated setup for a **native FEniCSx** installation using a **Python Virtual Environment (`venv`)** directly within a dedicated **Windows Subsystem for Linux 2 (WSL2)** distribution. This approach offers a clean, isolated, and locally controlled environment for simulation development.

>**IMPORTANT NOTE:** **The `install_local_fenicsx.bat` file must be executed from a dedicated Windows folder as well as `uninstall_local_fenicsx.bat`).**



## Directory Structure

The following tree represents the required structure of the installer folder. All files shown here are necessary for the installation process.

```
WSL_Local_installer/
├── install_local_fenicsx.bat
├── uninstall_local_fenicsx.bat
├── apt-installed-list.txt
└── pip-installed-list.txt
```

-----

## Key Features

  * **Native WSL2 Installation:** Installs FEniCSx directly onto the Linux filesystem within a dedicated Ubuntu 24.04 distribution (`FEniCSxenv`), providing near-native Linux performance.
  * **Isolated Virtual Environment:** Uses a **Python Virtual Environment (`venv`)** named `fenicsx-env` to manage FEniCSx and its dependencies, ensuring system-wide Python installations remain clean.
  * **Pre-defined Dependencies:** System (via `apt`) and Python (via `pip`) packages are automatically installed based on the content of the provided `.txt` list files, ensuring a consistent and reproducible setup.
  * **Automatic Activation:** The virtual environment is automatically sourced upon launching the WSL distribution, meaning you can start running FEniCSx commands immediately without manual activation.
  * **Windows-Native Launchers:** The installation and uninstallation processes are entirely driven by Windows `.bat` files.

-----

## Installation

1.  **Preparation:** Ensure **WSL2** is enabled and configured on your Windows machine.
2.  **Execution:** From the dedicated Windows folder containing all installer files, run the main installer script:
    ```bash
    install_local_fenicsx.bat
    ```
3.  **Prompts:** The script will first install the Ubuntu 24.04 distribution (`FEniCSxenv`). You will be prompted to:
      * Set a WSL **username and password** in the new Ubuntu window.
      * After the Ubuntu setup is complete, you must re-enter the **WSL username** in the installer's command prompt to allow the script to proceed with the FEniCSx installation, virtual environment setup, and auto-activation.
4.  **Completion:** The script will automatically verify the installation by checking the dolfinx version. A success or failure message will be displayed upon completion.

-----

## Usage

Since the environment is configured for **automatic activation**, usage is simple:

### 1\. Launching the FEniCSx Environment

Open your WSL distribution directly from the Windows search bar or by running this command in Windows Command Prompt/PowerShell:

```bash
wsl -d FEniCSxenv
```

You will be logged in, and the `fenicsx-env` virtual environment will be automatically activated. You are now ready to run your FEniCSx Python scripts (e.g., `python3 myscript.py` or using `mpirun`).

### 2\. Running Scripts

You can execute your FEniCSx scripts directly from the WSL command line:

```bash
# Run a single Python script
python3 filename.py

# Run a parallel simulation using mpirun
mpirun -np 4 python3 parallel_script.py
```

-----

## Uninstallation

To completely remove the installed FEniCSx environment:

1.  **Run Uninstaller:** Execute the uninstallation batch file from Windows:
    ```bash
    uninstall_local_fenicsx.bat
    ```
2.  **Warning:** This uses the `wsl --unregister` command, which is an **irreversible** action. It will delete the entire `FEniCSxenv` WSL distribution, including all data, the Python virtual environment, and any user files stored inside it.