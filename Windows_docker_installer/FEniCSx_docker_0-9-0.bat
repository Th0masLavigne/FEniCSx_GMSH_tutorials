@echo off
setlocal enabledelayedexpansion

REM ============================================================
REM  Automatic FEniCSx Docker Setup, Verification, and Runner
REM ============================================================

set "IMAGE_NAME=fenicsx"
set "TAG=v0.9.0"
set "FULL_IMAGE=%IMAGE_NAME%:%TAG%"
set "DOCKERFILE=Dockerfile"
set "MAIN_SCRIPT=Test.py"

REM ============================================================
REM Step 1: Ensure WSL2 is enabled (Docker backend requirement)
REM ============================================================
echo Checking WSL installation...
wsl --status >nul 2>&1
if %errorlevel% neq 0 (
    echo Enabling WSL and Virtual Machine Platform features...
    dism.exe /online /enable-feature /featurename:Microsoft-Windows-Subsystem-Linux /all /norestart >nul
    dism.exe /online /enable-feature /featurename:VirtualMachinePlatform /all /norestart >nul
    echo Please reboot your system, then rerun this installer.
    pause
    exit /b
)
wsl --set-default-version 2 >nul

REM ============================================================
REM Step 2: Check Docker Desktop installation
REM ============================================================
echo Checking Docker Desktop installation...
where docker >nul 2>&1
if %errorlevel% neq 0 (
    echo Docker not found. Downloading installer...
    powershell -Command "Invoke-WebRequest -Uri 'https://desktop.docker.com/win/stable/Docker%20Desktop%20Installer.exe' -OutFile '%TEMP%\DockerDesktopInstaller.exe'"
    echo Installing Docker Desktop silently...
    "%TEMP%\DockerDesktopInstaller.exe" install --quiet
    echo Please log out and back in after installation to complete setup.
    pause
)

REM ============================================================
REM Step 3: Start Docker Desktop
REM ============================================================
echo Starting Docker Desktop if not already running...
start "" "C:\Program Files\Docker\Docker\Docker Desktop.exe"
timeout /t 15 >nul

echo Waiting for Docker to become available...
:wait_docker
docker info >nul 2>&1
if %errorlevel% neq 0 (
    echo Docker not ready yet... waiting 5 seconds.
    timeout /t 5 >nul
    goto wait_docker
)

REM ============================================================
REM Step 4: Create Dockerfile dynamically
REM ============================================================
echo Creating Dockerfile...
(
    echo FROM dolfinx/dolfinx:%TAG%
    echo RUN apt-get update ^&^& \
    echo     apt-get install -y xvfb libgl1 libglu1-mesa mesa-utils ^&^& \
    echo     python3 -m pip install --upgrade pip ^&^& \
    echo     pip3 install pandas imageio pyvista ^&^& \
    echo     rm -rf /var/lib/apt/lists/*
) > "%DOCKERFILE%"
echo Dockerfile created successfully.

REM     In case you have custom packages, add to the Dockerfile (inside the .bat)
REM     python3 -m pip install ./src/pyproject.toml && \
REM     apt update 

REM ============================================================
REM Step 5: Build Docker image locally
REM ============================================================
echo Building Docker image %FULL_IMAGE% ...
docker build -f "%DOCKERFILE%" -t "%FULL_IMAGE%" .
if %errorlevel% neq 0 (
    echo ERROR: Docker image build failed.
    powershell -Command "Add-Type -AssemblyName PresentationFramework; [System.Windows.MessageBox]::Show('Docker image build FAILED!','FEniCSx Docker Installer',[System.Windows.MessageBoxButton]::OK,[System.Windows.MessageBoxImage]::Error)"
    pause
    exit /b
)
del "%DOCKERFILE%" >nul 2>&1
echo Docker image built successfully.

REM ============================================================
REM Step 6: Verify FEniCSx functionality with MPI parallel test
REM ============================================================
echo Running verification test inside container...

set "VERSION_LOG=version.log"
set "MPI_LOG=mpi.log"

echo Running DolfinX version check inside container...
docker run --rm -v "%cd%":/home/fenicsx/shared -w /home/fenicsx/shared %FULL_IMAGE% bash -c "python3 -c \"import dolfinx; print(dolfinx.__version__)\" > /home/fenicsx/shared/%VERSION_LOG% 2>&1"

echo Running MPI test inside container...
docker run --rm -v "%cd%":/home/fenicsx/shared -w /home/fenicsx/shared %FULL_IMAGE% bash -c "mpirun -n 2 python3 -c \"from mpi4py import MPI; print(MPI.COMM_WORLD.Get_rank())\" > /home/fenicsx/shared/%MPI_LOG% 2>&1"

REM Wait briefly to ensure file sync
timeout /t 2 >nul

REM ============================================================
REM Check logs, display, then delete
REM ============================================================
if exist "%VERSION_LOG%" (
    echo DolfinX version:
    type "%VERSION_LOG%"
    del "%VERSION_LOG%" >nul 2>&1
) else (
    echo ERROR: DolfinX version log missing. Docker test failed.
)
timeout /t 5 >nul

if exist "%MPI_LOG%" (
    echo MPI ranks (should print 0 1)
    type "%MPI_LOG%"
    del "%MPI_LOG%" 2>nul
) else (
    echo ERROR: MPI log missing. MPI test failed.
)


echo ============================================================
echo Verification complete.
echo ============================================================

REM ============================================================
REM Step 7: Optional — Run user Main.py if present
REM ============================================================
if exist "%MAIN_SCRIPT%" (
    echo Running user script: %MAIN_SCRIPT%
    docker run -it --rm ^
        -v "%cd%":/home/fenicsx/shared ^
        -w /home/fenicsx/shared ^
        %FULL_IMAGE% ^
        python3 "%MAIN_SCRIPT%"
) else (
    echo No Main.py found — skipping execution step.
)

echo ============================================================
echo Installation and verification complete!
echo Docker image: %FULL_IMAGE%
echo ============================================================
pause
endlocal