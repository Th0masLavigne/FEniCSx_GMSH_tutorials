@echo off
setlocal enabledelayedexpansion

REM =====================================================
REM FEniCSx Interactive Environment Launcher (WSL Version)
REM FIX: Prioritizes WSL UNC path parsing over Windows 'if exist' check.
REM =====================================================

set "DIST_NAME=FEniCSxenv"
set "IMAGE_NAME=fenicsx"
set "TAG=v0.9.0"
set "FULL_IMAGE=%IMAGE_NAME%:%TAG%"
set "DEFAULT_USER_NAME="
set "KNOWN_USER=tlavigne002"
set "DOCKER_MOUNT_POINT=/home/fenicsx/shared"

REM -----------------------------------------------------------------
REM Step 1: Get WSL Username
REM -----------------------------------------------------------------
echo Retrieving active WSL username for %DIST_NAME%...
for /f "delims=" %%U in ('wsl -d %DIST_NAME% -e whoami 2^>nul') do set "DEFAULT_USER_NAME=%%U"

if not defined DEFAULT_USER_NAME (
    set "DEFAULT_USER_NAME=!KNOWN_USER!"
    echo WARNING: Could not determine username, using fallback: !DEFAULT_USER_NAME!
)

REM -----------------------------------------------------------------
REM Step 2: Determine Mount Folder (Windows Side)
REM -----------------------------------------------------------------
:GET_TARGET_DIR
set "TARGET_DIR_INPUT="
if "%~1"=="" (
    set /p TARGET_DIR_INPUT="Enter folder to mount (leave empty for current directory): "
    if "!TARGET_DIR_INPUT!"=="" set "TARGET_DIR_INPUT=%cd%"
) else (
    set "TARGET_DIR_INPUT=%~1"
)

REM --- NEW ROBUST PATH PROCESSING ---
set "TARGET_DIR=!TARGET_DIR_INPUT!"

echo Checking path type for: !TARGET_DIR!

REM 1. HIGH PRIORITY: Check for WSL UNC path first
if "!TARGET_DIR:wsl.localhost=!" neq "!TARGET_DIR!" (
    
    echo Path is a WSL UNC path. Converting to native Linux host path...

    set "LINUX_ONLY_PATH=!TARGET_DIR: =!"
    
    REM Find the position of 'home' in the string (e.g., /home/tlavigne002/...)
    set "UNC_BASE_POS=0"
    for /l %%i in (1,1,256) do (
        if "!TARGET_DIR:~%%i,4!"=="home" (
            set /a UNC_BASE_POS=%%i
            goto FOUND_HOME_POS
        )
    )

    :FOUND_HOME_POS
    if !UNC_BASE_POS! equ 0 (
        echo ERROR: Could not parse /home path from WSL UNC.
        pause
        exit /b 1
    )

    REM Extract the Linux path component starting from /home/
    set "LINUX_ONLY_PATH=!TARGET_DIR: =!"
    set "LINUX_ONLY_PATH=!LINUX_ONLY_PATH:~%UNC_BASE_POS%!"
    
    REM Fix the backslashes to forward slashes, and prepend /
    set "WSL_MOUNT_DIR=^/!LINUX_ONLY_PATH:\=/!"
    
    goto PATH_CONVERSION_DONE

) else if exist "!TARGET_DIR!" (
    REM 2. SECOND PRIORITY: Path exists on Windows (e.g., C:\Data)

    echo Path is a standard Windows path. Converting to /mnt/c/ style...
    
    pushd "!TARGET_DIR!" >nul 2>&1
    if errorlevel 1 (
        echo Invalid path: !TARGET_DIR!
        pause
        exit /b 1
    )
    set "WINDOWS_TARGET_DIR=!cd!"
    popd >nul

    set "DRIVE_LETTER=!WINDOWS_TARGET_DIR:~0,1!"
    set "PATH_WITHOUT_DRIVE=!WINDOWS_TARGET_DIR:~3!"
    set "DRIVE_LOWER=!DRIVE_LETTER:A=a!"
    set "DRIVE_LOWER=!DRIVE_LOWER:B=b!"
    set "DRIVE_LOWER=!DRIVE_LOWER:C=c!"
    set "DRIVE_LOWER=!DRIVE_LOWER:D=d!"
    set "LINUX_PATH_PART=!PATH_WITHOUT_DRIVE:\=/!"
    set "WSL_MOUNT_DIR=/mnt/!DRIVE_LOWER!/!LINUX_PATH_PART!"
    
    goto PATH_CONVERSION_DONE

) else (
    REM 3. FAILURE: Path is neither a recognized WSL path nor an existing Windows path.
    echo ERROR: Path not found on Windows or recognized as WSL path.
    pause
    exit /b 1
)

:PATH_CONVERSION_DONE
echo Converted to FINAL WSL path: !WSL_MOUNT_DIR!

REM -----------------------------------------------------------------
REM Step 3: Execute Interactive Docker Session (inside WSL)
REM -----------------------------------------------------------------

REM Execute the command interactively. No changes needed here, as the problem was the path.
set "DOCKER_CMD=docker run -it --rm -v ""!WSL_MOUNT_DIR!""^:!DOCKER_MOUNT_POINT! -w !DOCKER_MOUNT_POINT! %FULL_IMAGE% bash"

echo.
echo Launching FEniCSx environment in WSL... (Ctrl+D or 'exit' to return)
echo Command: !DOCKER_CMD!

wsl -d %DIST_NAME% -u !DEFAULT_USER_NAME! -- bash -c "!DOCKER_CMD!"

echo.
echo ============================================================
echo FEniCSx interactive session ended.
echo ============================================================

REM -----------------------------------------------------------------
REM Step 5: Ask user if they want to run again
REM -----------------------------------------------------------------
:ASK_CONTINUE
set /p CHOICE="Do you want to run again? (Y/N): "
if /i "!CHOICE!"=="Y" (
    set "%~1="
    goto GET_TARGET_DIR
)
if /i "!CHOICE!"=="N" goto END
echo Please answer Y or N.
goto ASK_CONTINUE

:END
echo Exiting script.
pause
endlocal