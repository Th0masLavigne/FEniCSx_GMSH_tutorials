@echo off
setlocal enabledelayedexpansion

set "IMAGE_NAME=fenicsx"
set "TAG=v0.9.0"
set "FULL_IMAGE=%IMAGE_NAME%:%TAG%"
set "DOCKERFILE=Dockerfile"
set "MAIN_SCRIPT=Main.py"

REM Step 0: Start Docker Desktop
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
docker run --rm %FULL_IMAGE% echo "test" >nul 2>&1
if %errorlevel% neq 0 (
    echo ============================================================
    echo ERROR: Docker image %FULL_IMAGE% not found or unusable.
    echo Please run FEniCSx_docker_0-9-0 to build it first.
    echo ============================================================
    pause
    exit /b 1
)

:LOOP
REM ============================================================
REM Determine which folder to mount
REM ============================================================
if "%~1"=="" (
    REM Interactive prompt
    set "TARGET_DIR="
    set /p TARGET_DIR="Enter folder to mount (leave empty for current directory): "
    if "!TARGET_DIR!"=="" set "TARGET_DIR=%cd%"
) else (
    REM Use path from argument
    set "TARGET_DIR=%~1"
)

REM Convert relative paths to absolute
pushd "%TARGET_DIR%" >nul 2>&1
if errorlevel 1 (
    echo Invalid path: %TARGET_DIR%
    pause
    exit /b 1
)
set "TARGET_DIR=%cd%"
popd >nul

echo Using folder: %TARGET_DIR%

REM ============================================================
REM Open Interactive Docker
REM ============================================================
docker run -it --rm -v "%TARGET_DIR%":/home/fenicsx/shared -w /home/fenicsx/shared %FULL_IMAGE%

echo ============================================================
echo Docker image: %FULL_IMAGE%
echo ============================================================

REM ============================================================
REM Ask user if they want to continue
REM ============================================================
:ASK_CONTINUE
set /p CHOICE="Do you want to run again? (Y/N): "
if /i "!CHOICE!"=="Y" goto LOOP
if /i "!CHOICE!"=="N" goto END
echo Please answer Y or N.
goto ASK_CONTINUE

:END
echo Exiting script.
pause
endlocal
