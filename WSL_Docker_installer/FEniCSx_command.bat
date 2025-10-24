@echo off
setlocal enabledelayedexpansion

REM =====================================================
REM FEniCSx Command Runner (DIRECT DOCKER EXECUTION)
REM FIX: Bypass the unreliable alias sourcing by executing
REM the full Docker command directly.
REM =====================================================

set "DIST_NAME=FEniCSxenv"
set "FULL_IMAGE=fenicsx:v0.9.0"
set "DEFAULT_USER_NAME="
set "KNOWN_USER=tlavigne002"

echo Retrieving active WSL username for %DIST_NAME%...
for /f "delims=" %%U in ('wsl -d %DIST_NAME% -e whoami 2^>nul') do set "DEFAULT_USER_NAME=%%U"

if not defined DEFAULT_USER_NAME (
    set "DEFAULT_USER_NAME=!KNOWN_USER!"
    echo WARNING: Could not dynamically determine username, using fallback: !DEFAULT_USER_NAME!
)

REM --- The Path to the Script (Must be absolute for Docker) ---
REM Need the full Linux path to the script: /home/tlavigne002/Softwares/WSL2_Docker_FEniCSx/test_codes/Test.py
set "WSL_SCRIPT_PATH=./test_codes/Test.py"
set "WSL_PROJECT_DIR=/home/!DEFAULT_USER_NAME!/Softwares/WSL2_Docker_FEniCSx"

REM --- The actual Docker command that 'fenicsx-run' was meant to execute ---
REM 1. Get into the project directory
REM 2. Run Docker, mounting the project directory, and executing the Python script inside the container.
set "INNER_COMMAND=cd !WSL_PROJECT_DIR! && docker run --rm -v !WSL_PROJECT_DIR!:/home/fenics/shared -w /home/fenics/shared %FULL_IMAGE% python3 !WSL_SCRIPT_PATH!"

echo Detected User: !DEFAULT_USER_NAME!
echo Executing Docker Command Directly:
echo !INNER_COMMAND!
echo.

echo =================================================================
echo STARTING WSL EXECUTION
echo =================================================================

REM Execute the command directly as a single string via 'bash -c'
wsl -d %DIST_NAME% -u !DEFAULT_USER_NAME! -- bash -c "!INNER_COMMAND!"

set "EXIT_CODE=%errorlevel%"

echo.
echo =================================================================
echo WSL EXECUTION FINISHED
echo =================================================================
echo.

if %EXIT_CODE% neq 0 (
    echo !!! CRITICAL FAILURE !!!
    echo The script exited with code: %EXIT_CODE%.
    echo Please ensure the Docker daemon is running inside WSL.
) else (
    echo Command executed successfully.
)

pause
endlocal