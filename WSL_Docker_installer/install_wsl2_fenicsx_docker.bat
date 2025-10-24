@echo off
setlocal enabledelayedexpansion

REM =====================================================
REM Idempotent WSL2 + Ubuntu 24.04 + FEniCSx + Docker Installer
REM =====================================================

set "DIST_NAME=FEniCSxenv"
set "UBUNTU_DIST=Ubuntu-24.04"
set "IMAGE_NAME=fenicsx"
set "TAG=v0.9.0"
set "FULL_IMAGE=%IMAGE_NAME%:%TAG%"

REM -------------------------
REM Step 1: Ensure WSL installed and default version 2
REM -------------------------
echo Checking WSL installation...
wsl.exe --status >nul 2>&1
if %errorlevel% neq 0 (
    echo WSL not found. Enabling WSL and Virtual Machine Platform...
    dism.exe /online /enable-feature /featurename:Microsoft-Windows-Subsystem-Linux /all /norestart >nul
    dism.exe /online /enable-feature /featurename:VirtualMachinePlatform /all /norestart >nul
    echo WSL features enabled. Please reboot your system and rerun this installer.
    pause
    exit /b
)

wsl --set-default-version 2

REM -------------------------
REM Step 2: Check if WSL distro exists
REM -------------------------
echo Installing Ubuntu 24.04 as %DIST_NAME% (if not already installed)...
echo Press any key to continue - you will create your WSL username and password. [press any key]
echo -------------------------
echo If first install create username and password
echo The password keeps invisible when typed. Remember it.
echo Type exit when the credentials are created to continue
echo -------------------------
echo When message has been read [press any key]
pause >nul
REM New installation triggered
powershell -Command "Add-Type -AssemblyName PresentationFramework; [System.Windows.MessageBox]::Show('During the first install, username and password have to be created. The password will be invisible. Once done, type exit to return here.');"

wsl --install -d %UBUNTU_DIST% --name %DIST_NAME% 


REM -------------------------
REM Step 3: Get WSL default username
REM -------------------------
for /f "delims=" %%U in ('wsl -d %DIST_NAME% -e whoami') do set "USER_NAME=%%U"
echo Detected WSL username: %USER_NAME%

REM -------------------------
REM Step 4: Update Ubuntu system
REM -------------------------
echo Updating WSL system...
wsl -d %DIST_NAME% -u %USER_NAME% -- bash -c "sudo apt update && sudo apt -y upgrade && sudo apt -y autoremove"

echo Installation and update complete.

REM -------------------------
REM Step 5: Check if Docker is installed inside WSL and ensure daemon is running
REM -------------------------

echo Installing Docker Engine...
wsl -d %DIST_NAME% -u %USER_NAME% -- bash docker_install.sh


echo Docker installation and daemon setup complete.


REM -------------------------
REM Step 6: Create Dockerfile and build it
REM -------------------------
echo Building FEniCSx Docker image...
wsl -d %DIST_NAME% -u %USER_NAME% -- bash docker_fenicsx_env.sh

REM -------------------------
REM Step 7: Create aliases if missing
REM -------------------------
wsl -d %DIST_NAME% -u %USER_NAME% -- bash -c "ALIAS_FILE=~/.bash_aliases && touch \$ALIAS_FILE && grep -qxF 'alias fenicsx=\"docker run -it --rm -v $(pwd):/home/fenicsx/shared -w /home/fenicsx/shared %FULL_IMAGE% bash\"' \$ALIAS_FILE || echo 'alias fenicsx=\"docker run -it --rm -v $(pwd):/home/fenicsx/shared -w /home/fenicsx/shared %FULL_IMAGE% bash\"' >> \$ALIAS_FILE && grep -qxF 'alias fenicsx-run=\"docker run --rm -v $(pwd):/home/fenicsx/shared -w /home/fenicsx/shared %FULL_IMAGE% bash -c\"' \$ALIAS_FILE || echo 'alias fenicsx-run=\"docker run --rm -v $(pwd):/home/fenicsx/shared -w /home/fenicsx/shared %FULL_IMAGE% bash -c\"' >> \$ALIAS_FILE && echo 'source ~/.bash_aliases' >> ~/.bashrc"


REM -------------------------
REM Step 8: Verification: Docker image, DolfinX version, MPI
REM -------------------------
echo =====================================================
echo Running verification inside WSL2...
wsl -d %DIST_NAME% -u %USER_NAME% -- bash -c "docker run --rm -v $(pwd):/home/fenicsx/shared -w /home/fenicsx/shared fenicsx:v0.9.0 bash -c \"python3 -c 'import dolfinx; print(\\\"DolfinX version:\\\", dolfinx.__version__)'\" "

wsl -d %DIST_NAME% -u %USER_NAME% -- bash -c "docker run --rm -v $(pwd):/home/fenicsx/shared -w /home/fenicsx/shared fenicsx:v0.9.0 bash -c \"mpirun -n 2 python3 -c 'from mpi4py import MPI; print(\\\"MPI ranks:\\\", MPI.COMM_WORLD.Get_rank())'\" "

echo =====================================================
echo Installation and verification complete!
echo Launch your environment with:
echo     wsl -d %DIST_NAME%
echo Use:
echo     fenicsx       - interactive FEniCSx shell
echo     fenicsx-run   - execute a bash command in FEniCSx container
echo =====================================================


pause
endlocal
