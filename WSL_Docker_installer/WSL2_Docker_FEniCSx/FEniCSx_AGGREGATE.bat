@echo off
setlocal enabledelayedexpansion


set "IMAGE_NAME=fenicsx"
set "TAG=v0.9.0"
set "FULL_IMAGE=%IMAGE_NAME%:%TAG%"
set "DOCKERFILE=Dockerfile"
set "MAIN_SCRIPT=Main.py"


REM Step 0: Start Docker Desktop
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

REM Step 1: Check if Docker image exists
REM ============================================================
docker run --rm %FULL_IMAGE% echo "test" >nul 2>&1
if %errorlevel% neq 0 (
    echo ============================================================
    echo ERROR: Docker image %FULL_IMAGE% not found or unusable.
    echo Please run FEniCSx_docker_0-9-0 to build it first.
    echo ============================================================
    pause
    exit /b 1
)

REM ============================================================
REM Run user Main.py if present
REM ============================================================
if exist "%MAIN_SCRIPT%" (
    echo Running user script: %MAIN_SCRIPT%
    docker run -it --rm -v "%cd%":/home/fenicsx/shared -w /home/fenicsx/shared %FULL_IMAGE% python3 "%MAIN_SCRIPT%"
) else (
    echo No Main.py found — skipping execution step.
)

echo ============================================================
echo Docker image: %FULL_IMAGE%
echo ============================================================
pause
endlocal