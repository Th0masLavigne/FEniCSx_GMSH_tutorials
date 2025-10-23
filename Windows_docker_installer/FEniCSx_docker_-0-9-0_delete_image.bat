@echo off
setlocal enabledelayedexpansion

REM ============================================================
REM Automatic FEniCSx Docker Uninstaller
REM ============================================================

set "IMAGE_NAME=fenicsx"
set "TAG=v0.9.0"
set "FULL_IMAGE=%IMAGE_NAME%:%TAG%"

echo ============================================================
echo Stopping any running containers using image %FULL_IMAGE%...
echo ============================================================

REM List running containers using this image
for /f "tokens=*" %%i in ('docker ps -q --filter "ancestor=%FULL_IMAGE%"') do (
    echo Stopping container %%i ...
    docker stop %%i >nul 2>&1
    docker rm %%i >nul 2>&1
)

echo ============================================================
echo Removing Docker image %FULL_IMAGE%...
echo ============================================================

REM Check if image exists before removing
docker images -q %FULL_IMAGE% >nul 2>&1
if %errorlevel%==0 (
    docker rmi -f %FULL_IMAGE%
    if %errorlevel%==0 (
        echo Docker image %FULL_IMAGE% removed successfully.
    ) else (
        echo ERROR: Failed to remove Docker image %FULL_IMAGE%.
    )
) else (
    echo Docker image %FULL_IMAGE% not found.
)

echo ============================================================
echo FEniCSx Docker uninstallation complete.
echo ============================================================
pause
endlocal
