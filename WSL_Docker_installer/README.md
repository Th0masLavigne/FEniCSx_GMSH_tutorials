# 🐳 FEniCSx Docker Environment Installer for WSL2

This repository contains the scripts necessary to automate the installation of a **FEniCSx** simulation environment running inside a **Docker container** on **Windows Subsystem for Linux 2 (WSL2)**.

This setup prioritizes a clean, isolated, and easily reproducible environment, avoiding the complexities and potential instability of native Python virtual environments or Docker Desktop's direct Windows integration for FEniCSx installations.

## 🚀 Key Features

* **Isolated Environment:** FEniCSx and its dependencies run in a secure Docker container, ensuring a consistent environment regardless of your host system's state.
* **WSL2 Integration:** Leverages the performance and Linux compatibility of WSL2.
* **Idempotent Installation:** The installer can be run multiple times; it checks for existing installations and offers a full reinstall option to clean up and start fresh.
* **Self-Contained:** All installation files are copied into the WSL distribution (`~/Softwares/WSL2_Docker_FEniCSx`), allowing the original Windows folder to be deleted after a successful install.
* **Custom Packages:** The setup automatically checks for a `pyproject.toml` file in `src_install/custom_packages` and integrates it into the Docker image build.
* **Simplified Usage:** Custom bash aliases (`fenicsx` and `fenicsx-run`) are automatically created for quick access to the container's shell and for running non-interactive scripts.

## 📦 Prerequisites

Before running the installer, ensure you have:

1.  **Windows 10/11** (with WSL2 support enabled).
2.  **Administrator privileges** for installing/enabling WSL components and Docker.

The installer will handle the installation of WSL features and the required Ubuntu distribution if they are missing.

## 📋 Installation Guide

1.  **Download/Clone:** Download or clone this repository to a local folder (e.g., `C:\Users\YourUser\WSL2_Docker_FEniCSx`).
2.  **Add Custom Packages (Optional):** If you have custom Python packages, place your build configuration (e.g., `pyproject.toml`) inside the `src_install/custom_packages/` folder. The installer will automatically detect and integrate it.
3.  **Run Installer:** Execute the main installation batch file:
    ```bash
    install_wsl2_fenicsx_docker.bat
    ```
4.  **Follow Prompts:**
    * The script will guide you through enabling WSL features (may require a reboot).
    * It will install the **Ubuntu 24.04** distribution as `FEniCSxenv` and prompt you to create a **username and password**.
    * **Crucially:** You will be prompted to enter your **WSL password** in a secure, hidden window to allow non-interactive `sudo` commands for the rest of the installation (e.g., Docker installation).
    * If the environment already exists, you will be asked if you wish to proceed with a **full reinstall**.
5.  **Completion:** Upon successful completion, a `report.log` file will contain the installation details. The original Windows folder can be deleted.

## 🛠 Usage

1.  **Launch WSL Environment:**
    Open a command prompt or PowerShell and launch the specific WSL distribution:
    ```bash
    wsl -d FEniCSxenv
    ```

2.  **Interactive Shell (`fenicsx`):**
    To start an interactive bash shell inside the FEniCSx Docker container (mounting your current directory as `/home/fenicsx/shared`):
    ```bash
    fenicsx
    ```

3.  **Run Script/Command (`fenicsx-run`):**
    To execute a command or a Python script non-interactively within the container:
    ```bash
    # Example: Run a Python script named Main.py in the current directory
    fenicsx-run "python3 Main.py"
    ```

## 🧹 Uninstallation

To completely remove the installed FEniCSx environment and its files:

1.  **Run Uninstaller:** Execute the uninstallation batch file:
    ```bash
    uninstall_wsl2_fenicsx_docker.bat
    ```
    This script unregisters the `FEniCSxenv` WSL distribution, deleting all its contents.