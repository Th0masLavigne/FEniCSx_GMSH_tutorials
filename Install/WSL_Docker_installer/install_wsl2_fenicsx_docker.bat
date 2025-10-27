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
set "LOG_FILE=report.log"
set "WINDOWS_INSTALL_DIR=%cd%"
set "LOCAL_FOLDER_NAME=WSL2_Docker_FEniCSx"
set "CUSTOM_PACKAGE_FLAG=0" 
set "USER_NAME=" REM Initialize USER_NAME as empty

REM -------------------------
REM Log File Initialization
REM -------------------------
echo Installation started on %date% at %time% > "%LOG_FILE%"
echo ---------------------------------------------------------------->> "%LOG_FILE%"

REM -------------------------
REM Custom Package Check
REM -------------------------
call :check_custom_package


REM -------------------------
REM NEW: Step 0.5: Get Custom Distro Name from User
REM -------------------------
echo Prompting user for custom distribution name...
echo Prompting user for custom distribution name... >> "%LOG_FILE%"

set "DEFAULT_DIST_NAME=FEniCSxenv"
set "INPUT_FILE=DistroNameInput.txt"

powershell -Command "Add-Type -AssemblyName PresentationFramework; Add-Type -AssemblyName Microsoft.VisualBasic; [Microsoft.VisualBasic.Interaction]::InputBox('Enter the custom name for your WSL2 FEniCSx environment (e.g., FEniCSxenv, MySimEnv).\n\nIf the distro already exists, the installer will check for it in the next step.','Custom WSL Distro Name','!DEFAULT_DIST_NAME!') | Out-File -FilePath '!INPUT_FILE!' -Encoding ASCII"

set "USER_DIST_NAME="
for /f "usebackq delims=" %%i in ("!INPUT_FILE!") do set "USER_DIST_NAME=%%i"

del "!INPUT_FILE!"

REM Remove carriage returns and newlines that PowerShell sometimes adds
set "USER_DIST_NAME=!USER_DIST_NAME:\r=!"
set "USER_DIST_NAME=!USER_DIST_NAME:\n=!"

REM Use the user input if it's not empty, otherwise default
if not defined USER_DIST_NAME (
    set "DIST_NAME=!DEFAULT_DIST_NAME!"
) else (
    set "DIST_NAME=!USER_DIST_NAME!"
)

echo Using distribution name: %DIST_NAME%
echo Using distribution name: %DIST_NAME% >> "%LOG_FILE%"

REM -------------------------
REM Step 1: Ensure WSL installed and default version 2
REM -------------------------
echo Checking WSL installation...
echo Checking WSL installation... >> "%LOG_FILE%"
wsl.exe --status >nul 2>&1
if %errorlevel% neq 0 goto :wsl_not_found

:wsl_check_done

REM WSL is installed, set default version 2 and continue
REM FIX: Discard output (>nul 2>&1) to avoid logging messy localized text.
wsl --set-default-version 2 >nul 2>&1
if %errorlevel% neq 0 goto :wsl_fail

REM Fall through to Step 2

REM -------------------------
REM Step 2: Check if WSL distro exists and handle reinstall
REM -------------------------
echo Checking if WSL distribution %DIST_NAME% exists using direct access...
echo Checking if WSL distribution %DIST_NAME% exists using direct access... >> "%LOG_FILE%"

REM Robustly check for distribution existence
wsl -d %DIST_NAME% -- exit 2>nul 
if %errorlevel% equ 0 goto :distro_exists

REM Distro does not exist (errorlevel was non-zero)
goto :new_installation_flow

:distro_exists
echo WARNING: The environment %DIST_NAME% already exists.
echo WARNING: The environment %DIST_NAME% already exists. >> "%LOG_FILE%"
    
powershell -Command "Add-Type -AssemblyName PresentationFramework; $choice = [System.Windows.MessageBox]::Show('The environment %DIST_NAME% exists. Do you want a full reinstall (Yes), continue configuration (No), or exit (Cancel)?','WSL Environment Exists','YesNoCancel','Warning'); [IO.File]::WriteAllText('ReinstallChoice.txt', $choice.ToString().Trim(), [System.Text.Encoding]::ASCII)"

set "REINSTALL_CHOICE="
for /f "usebackq tokens=1" %%i in ("ReinstallChoice.txt") do set "REINSTALL_CHOICE=%%i"

del ReinstallChoice.txt

set "REINSTALL_CHOICE=!REINSTALL_CHOICE: =!"
set "REINSTALL_CHOICE=!REINSTALL_CHOICE:\r=!"
set "REINSTALL_CHOICE=!REINSTALL_CHOICE:\n=!"

if /i "!REINSTALL_CHOICE!"=="Yes" goto :confirm_loss
if /i "!REINSTALL_CHOICE!"=="No" goto :skip_reinstall
goto :exit_flow
    
:skip_reinstall
echo Skipping full reinstall.
echo Skipping full reinstall. Continuing with existing distro configuration (Docker, Aliases). >> "%LOG_FILE%"
        
for /f "delims=" %%U in ('wsl -d %DIST_NAME% -e whoami') do set "USER_NAME=%%U"
echo Detected existing WSL username: !USER_NAME!
echo Detected existing WSL username: !USER_NAME! >> "%LOG_FILE%"
        
set "WSL_TARGET_DIR=/home/!USER_NAME!/Softwares/WSL2_Docker_FEniCSx"
echo WSL target directory set to: !WSL_TARGET_DIR!
echo WSL target directory set to: !WSL_TARGET_DIR! >> "%LOG_FILE%"

goto SETUP_SUDO_ACCESS


:confirm_loss
powershell -Command "Add-Type -AssemblyName PresentationFramework; $choice = [System.Windows.MessageBox]::Show('Do you understand that all files inside the WSL2 %DIST_NAME% will be lost? Do you still want a full reinstall (Yes/No)?','Confirm Reinstall','YesNo','Warning'); [IO.File]::WriteAllText('ConfirmLossChoice.txt', $choice.ToString().Trim(), [System.Text.Encoding]::ASCII)"
set "CONFIRM_LOSS="
for /f "usebackq tokens=1" %%i in ("ConfirmLossChoice.txt") do set "CONFIRM_LOSS=%%i"
del ConfirmLossChoice.txt

if /i "!CONFIRM_LOSS!"=="Yes" goto :unregister_distro

echo Reinstallation aborted.
echo Reinstallation aborted. Exiting. >> "%LOG_FILE%"
goto :exit_flow


:unregister_distro
echo Unregistering existing %DIST_NAME%...
echo Unregistering existing %DIST_NAME%... >> "%LOG_FILE%"
wsl --unregister %DIST_NAME% >> "%LOG_FILE%" 2>&1
if %errorlevel% neq 0 goto :wsl_fail
REM Fall through to :new_installation_flow


:new_installation_flow
echo Installing Ubuntu 24.04 as %DIST_NAME%...
echo Installing Ubuntu 24.04 as %DIST_NAME%... >> "%LOG_FILE%"
powershell -Command "Add-Type -AssemblyName PresentationFramework; [System.Windows.MessageBox]::Show('During the first install, a username and password need to be created in the new terminal window. The password will be invisible when typed. Once done, type exit to return here.','First-Time WSL Setup','OK','Information');"
wsl --install -d %UBUNTU_DIST% --name %DIST_NAME% 
if %errorlevel% neq 0 goto :wsl_fail

REM -------------------------
REM Step 3: Get WSL default username and set up path variables
REM -------------------------
echo Retrieving WSL username and setting path variables...
echo Retrieving WSL username and setting path variables... >> "%LOG_FILE%"
for /f "delims=" %%U in ('wsl -d %DIST_NAME% -e whoami') do set "USER_NAME=%%U"
echo Detected WSL username: !USER_NAME!
echo Detected WSL username: !USER_NAME! >> "%LOG_FILE%"

set "WSL_TARGET_DIR=/home/!USER_NAME!/Softwares/WSL2_Docker_FEniCSx"
echo WSL target directory set to: !WSL_TARGET_DIR!
echo WSL target directory set to: !WSL_TARGET_DIR! >> "%LOG_FILE%"

REM Fall through to SETUP_SUDO_ACCESS

REM -------------------------
REM NEW CENTRALIZED BLOCK: Setup Sudo Access
REM -------------------------
:SETUP_SUDO_ACCESS
echo Setting up sudo access for non-interactive commands...
echo Setting up sudo access for non-interactive commands... >> "%LOG_FILE%"

REM The password collection block is removed because the subsequent command
REM runs as the 'root' user, which does not require the user's password to execute.
wsl -d %DIST_NAME% -u root -- bash -c "echo '!USER_NAME! ALL=(ALL) NOPASSWD:ALL' | tee /etc/sudoers.d/99-nopasswd-fenicsx" >> "%LOG_FILE%" 2>&1
if %errorlevel% neq 0 goto :wsl_fail

:CONFIGURE_FS
REM -------------------------
REM Step 4: Configure WSL Filesystem (create directories and copy files)
REM -------------------------
echo Configuring WSL filesystem and copying installation files...
echo Configuring WSL filesystem and copying installation files... >> "%LOG_FILE%"
wsl -d %DIST_NAME% -u !USER_NAME! -- bash -c "mkdir -p ~/Softwares && mkdir -p ~/Documents" >> "%LOG_FILE%" 2>&1
if %errorlevel% neq 0 goto :wsl_fail
echo Created ~/Softwares and ~/Documents.
echo Created ~/Softwares and ~/Documents. >> "%LOG_FILE%"

echo Copying %LOCAL_FOLDER_NAME% from %WINDOWS_INSTALL_DIR% to !WSL_TARGET_DIR!...
echo Copying %LOCAL_FOLDER_NAME% from %WINDOWS_INSTALL_DIR% to !WSL_TARGET_DIR!... >> "%LOG_FILE%"

set "WIN_INSTALL_DRIVE=!WINDOWS_INSTALL_DIR:~0,1!"
set "WIN_INSTALL_PATH_REL=!WINDOWS_INSTALL_DIR:~3!"

set "DRIVE_LOWER=!WIN_INSTALL_DRIVE:A=a!"
set "DRIVE_LOWER=!DRIVE_LOWER:B=b!"
set "DRIVE_LOWER=!DRIVE_LOWER:C=c!"
set "DRIVE_LOWER=!DRIVE_LOWER:D=d!"
if /i "!DRIVE_LOWER!" equ "!WIN_INSTALL_DRIVE!" set "DRIVE_LOWER=!DRIVE_LOWER!"
set "WIN_INSTALL_PATH_REL=!WIN_INSTALL_PATH_REL:\=/!"

set "LINUX_SOURCE_DIR=/mnt/!DRIVE_LOWER!/!WIN_INSTALL_PATH_REL!"

REM 1. Remove old copy in case we skipped reinstall
wsl -d %DIST_NAME% -u !USER_NAME! -- bash -c "rm -rf !WSL_TARGET_DIR!" >> "%LOG_FILE%" 2>&1
if %errorlevel% neq 0 goto :wsl_fail
REM 2. Create the target directory
wsl -d %DIST_NAME% -u !USER_NAME! -- bash -c "mkdir -p !WSL_TARGET_DIR!" >> "%LOG_FILE%" 2>&1
if %errorlevel% neq 0 goto :wsl_fail

REM 3. Copy contents from Windows mount to WSL.
echo Linux source path: !LINUX_SOURCE_DIR!
echo Linux source path: !LINUX_SOURCE_DIR! >> "%LOG_FILE%"
wsl -d %DIST_NAME% -u !USER_NAME! -- bash -c "cp -r '!LINUX_SOURCE_DIR!/.' '!WSL_TARGET_DIR!'" >> "%LOG_FILE%" 2>&1
if %errorlevel% neq 0 goto :wsl_fail

echo Files copied and permissions set.
echo Files copied and permissions set. >> "%LOG_FILE%"

REM Set execution permissions
wsl -d %DIST_NAME% -u !USER_NAME! -- bash -c "chmod +x !WSL_TARGET_DIR!/src_install/docker/*.sh !WSL_TARGET_DIR!/src_install/apt_environment_aliases/*.sh" >> "%LOG_FILE%" 2>&1
if %errorlevel% neq 0 goto :wsl_fail


REM -------------------------
REM Step 5: Update Ubuntu system
REM -------------------------
echo Updating WSL system... This might take a while...
echo Updating WSL system... This might take a while... >> "%LOG_FILE%"
wsl -d %DIST_NAME% -u !USER_NAME! -- bash -c "sudo apt update && sudo apt -y upgrade && sudo apt -y autoremove" >> "%LOG_FILE%" 2>&1
if %errorlevel% neq 0 goto :wsl_fail
echo WSL system updated.
echo WSL system updated. >> "%LOG_FILE%"

REM -------------------------
REM Step 6: Install Docker Engine
REM -------------------------
echo Installing Docker Engine...
echo Installing Docker Engine... >> "%LOG_FILE%"
wsl -d %DIST_NAME% -u !USER_NAME! -- bash !WSL_TARGET_DIR!/src_install/docker/docker_install_check.sh >> "%LOG_FILE%" 2>&1
if %errorlevel% neq 0 goto :wsl_fail
echo Docker installation and daemon setup complete.
echo Docker installation and daemon setup complete. >> "%LOG_FILE%"

REM -------------------------
REM Step 7: Create Dockerfile and build it (FIXED CONTEXT)
REM -------------------------
echo Building FEniCSx Docker image...
echo Building FEniCSx Docker image... >> "%LOG_FILE%"
wsl -d %DIST_NAME% -u !USER_NAME! -- bash -c "export IMAGE_NAME=%IMAGE_NAME%; export TAG=%TAG%; export FULL_IMAGE=%FULL_IMAGE%; export CUSTOM_PACKAGE_FLAG=%CUSTOM_PACKAGE_FLAG%; cd %WSL_TARGET_DIR% && ./src_install/docker/docker_fenicsx_env.sh" >> "%LOG_FILE%" 2>&1
if %errorlevel% neq 0 goto :wsl_fail
echo FEniCSx Docker image built.
echo FEniCSx Docker image built. >> "%LOG_FILE%"

REM -------------------------
REM Step 8: Create aliases (FIXED CONTEXT)
REM -------------------------
echo Creating aliases...
echo Creating aliases... >> "%LOG_FILE%"
wsl -d %DIST_NAME% -u !USER_NAME! -- bash -c "export FULL_IMAGE=%FULL_IMAGE%; cd %WSL_TARGET_DIR% && ./src_install/apt_environment_aliases/fenicsx_aliases.sh" >> "%LOG_FILE%" 2>&1
if %errorlevel% neq 0 goto :wsl_fail
echo Aliases created.
echo Aliases created. >> "%LOG_FILE%"

REM -------------------------
REM Step 9: Verification (FIXED CONTEXT)
REM -------------------------
echo Running verification inside WSL2...
echo Running verification inside WSL2... >> "%LOG_FILE%"
wsl -d %DIST_NAME% -u !USER_NAME! -- bash -c "export FULL_IMAGE=%FULL_IMAGE%; cd %WSL_TARGET_DIR% && ./src_install/apt_environment_aliases/verification.sh" >> "%LOG_FILE%" 2>&1
if %errorlevel% neq 0 goto :wsl_fail
echo Verification complete.
echo Verification complete. >> "%LOG_FILE%"

REM Cleanup: Remove NOPASSWD setting
wsl -d %DIST_NAME% -u root -- bash -c "rm -f /etc/sudoers.d/99-nopasswd-fenicsx" >> "%LOG_FILE%" 2>&1


echo =====================================================
echo ===================================================== >> "%LOG_FILE%"
echo Installation and verification complete!
echo Installation and verification complete! >> "%LOG_FILE%"

REM -------------------------
REM NEW: Final Success Pop-up (English, using temporary file for robust display)
REM -------------------------

set "DIST_NAME_ESC=!DIST_NAME!"
set "FULL_IMAGE_ESC=!FULL_IMAGE!"
set "WINDOWS_INSTALL_DIR_ESC=!WINDOWS_INSTALL_DIR!"
set "WSL_TARGET_DIR_ESC=!WSL_TARGET_DIR!"

REM 1. Create a clean, multi-line message file (%TEMP%\success_msg.txt)

(
echo Installation of WSL2 Distribution %DIST_NAME_ESC% is successful!
echo.
echo Docker image %FULL_IMAGE_ESC% is ready.
echo Verification checks passed successfully.
echo.
echo Aliases created:
echo   - fenicsx: Interactive shell in the container.
echo   - fenicsx-run: Execute a bash command [Ex: fenicsx-run "mpirun -n 4 python3 my_script.py"].
echo.
echo To run, double-click on 'execute.bat' or 'execute_interactive.bat'.
echo To uninstall, double-click on 'uninstall.bat'.
echo.
echo NOTE: These .bat files are in your Windows file explorer at:
echo %WINDOWS_INSTALL_DIR_ESC%
echo.
echo To launch the environment manually: wsl -d %DIST_NAME_ESC%
echo report.log copied to WSL location: !WSL_TARGET_DIR_ESC!
echo You may safely **DELETE** the Windows folder: %WINDOWS_INSTALL_DIR_ESC%
) > "%TEMP%\success_msg.txt"


set "SUCCESS_TITLE=FEniCSx WSL2 Installation Complete"

REM 2. Use PowerShell to read the file, replace Windows newlines with PowerShell newlines (`n), and display the message.
powershell -Command "$msg = [System.IO.File]::ReadAllText('%TEMP%\success_msg.txt'); $msg = $msg -replace \"`r`n\", \"\`n\"; Add-Type -AssemblyName PresentationFramework; [System.Windows.MessageBox]::Show($msg,'%SUCCESS_TITLE%','OK','Information')"

REM 3. Clean up the temporary file
del "%TEMP%\success_msg.txt"

echo Pop-up displayed and awaiting user acknowledgement. >> "%LOG_FILE%"


echo =====================================================
echo ===================================================== >> "%LOG_FILE%"
echo Installation and verification complete!
echo Installation and verification complete! >> "%LOG_FILE%"
echo Launch your environment with:
echo Launch your environment with: >> "%LOG_FILE%"
echo     wsl -d %DIST_NAME%
echo     wsl -d %DIST_NAME% >> "%LOG_FILE%"
echo Use:
echo Use: >> "%LOG_FILE%"
echo     fenicsx       - interactive FEniCSx shell
echo     fenicsx       - interactive FEniCSx shell >> "%LOG_FILE%"
echo     fenicsx-run   - execute a bash command in FEniCSx container
echo     fenicsx-run   - execute a bash command in FEniCSx container >> "%LOG_FILE%"
echo =====================================================
echo ===================================================== >> "%LOG_FILE%"
echo.
echo NOTE: Since all necessary files are now copied inside the WSL distribution at:
echo NOTE: Since all necessary files are now copied inside the WSL distribution at: >> "%LOG_FILE%"
echo       \\wsl.localhost\%DIST_NAME%\home\!USER_NAME!\Softwares\%LOCAL_FOLDER_NAME%
echo       \\wsl.localhost\%DIST_NAME%\home\!USER_NAME!\Softwares\%LOCAL_FOLDER_NAME% >> "%LOG_FILE%"
echo       You may safely **DELETE** the Windows folder: %WINDOWS_INSTALL_DIR%
echo       You may safely **DELETE** the Windows folder: %WINDOWS_INSTALL_DIR% >> "%LOG_FILE%"
echo.

call :CONVERT_LOG_TO_ASCII

pause
goto :eof


:wsl_not_found
echo WSL not found.
echo WSL not found. Enabling WSL and Virtual Machine Platform... >> "%LOG_FILE%"
dism.exe /online /enable-feature /featurename:Microsoft-Windows-Subsystem-Linux /all /norestart >nul
dism.exe /online /enable-feature /featurename:VirtualMachinePlatform /all /norestart >nul
echo WSL features enabled.
echo WSL features enabled. Please reboot your system and rerun this installer. >> "%LOG_FILE%"
pause
exit /b

:exit_flow
echo Exiting installer as requested by user.
echo Exiting installer as requested by user. >> "%LOG_FILE%"

call :CONVERT_LOG_TO_ASCII

pause
endlocal
exit /b 0

:wsl_fail
REM -----------------------------------------------------------------
REM CENTRALIZED FAILURE EXIT BLOCK
REM -----------------------------------------------------------------
echo.
echo =====================================================
echo !!! CRITICAL FAILURE !!!
echo =====================================================
echo A critical WSL command failed with exit code %errorlevel%.
echo Please review "%LOG_FILE%" for the complete error output.
echo =====================================================
echo.
echo A critical WSL command failed with exit code %errorlevel%. >> "%LOG_FILE%"
echo Please review "%LOG_FILE%" for the complete error output. >> "%LOG_FILE%"
echo.

call :CONVERT_LOG_TO_ASCII

pause
endlocal
exit /b 1

REM =====================================================
REM SUBROUTINES
REM =====================================================

:CONVERT_LOG_TO_ASCII
REM -----------------------------------------------------------------
REM --- Subroutine to convert the entire log file to ASCII safely ---
REM --- This includes copying the cleaned version back to WSL.    ---
REM -----------------------------------------------------------------
set "TEMP_LOG=report_temp_conversion.log"

echo.
echo =====================================================
echo Cleaning up and converting report.log to standard ASCII...
echo =====================================================
echo.

REM 1. Convert content to a temporary file (OEM to ASCII fix)
powershell -Command "Get-Content -Path '%LOG_FILE%' -Encoding OEM | Out-File -FilePath '%TEMP_LOG%' -Encoding ASCII"

REM 2. Check if the temporary file was created successfully
if not exist "%TEMP_LOG%" (
    echo WARNING: Failed to create temporary log file. Skipping conversion.
    goto :eof
)

REM 3. Delete the original (broken) log file
del "%LOG_FILE%"

REM 4. Rename the clean, converted temporary file to the final log file name
rename "%TEMP_LOG%" "%LOG_FILE%"

REM 5. Copy the clean log file back to the WSL installation folder.
if defined USER_NAME (
    echo Copying clean report.log to WSL...
    
    REM Re-calculate LINUX_SOURCE_FILE for the Windows log file location
    set "WIN_INSTALL_DRIVE=!WINDOWS_INSTALL_DIR:~0,1!"
    set "WIN_INSTALL_PATH_REL=!WINDOWS_INSTALL_DIR:~3!"
    set "DRIVE_LOWER=!WIN_INSTALL_DRIVE:A=a!"
    set "DRIVE_LOWER=!DRIVE_LOWER:B=b!"
    set "DRIVE_LOWER=!DRIVE_LOWER:C=c!"
    set "DRIVE_LOWER=!DRIVE_LOWER:D=d!"
    if /i "!DRIVE_LOWER!" equ "!WIN_INSTALL_DRIVE!" set "DRIVE_LOWER=!DRIVE_LOWER!"
    set "WIN_INSTALL_PATH_REL=!WIN_INSTALL_PATH_REL:\=/!"
    set "LINUX_SOURCE_FILE=/mnt/!DRIVE_LOWER!/!WIN_INSTALL_PATH_REL!/%LOG_FILE%"

    REM Execute Linux copy command using the existing WSL paths
    wsl -d %DIST_NAME% -u !USER_NAME! -- bash -c "cp -f \"!LINUX_SOURCE_FILE!\" \"!WSL_TARGET_DIR!/report.log\"" >> nul 2>&1
    
    if %errorlevel% equ 0 (
        echo report.log copied to WSL location: !WSL_TARGET_DIR!
    ) else (
        echo WARNING: Failed to copy report.log to WSL. The log file remains in your Windows folder.
    )
)

echo Report.log conversion and sync complete.
echo.
goto :eof


:check_custom_package
if not exist "%WINDOWS_INSTALL_DIR%\src_install\custom_packages\pyproject.toml" goto :no_custom_package

REM Custom package found
set "CUSTOM_PACKAGE_FLAG=1"
echo Custom package (pyproject.toml) detected.
echo Custom package (pyproject.toml) detected. It will be installed in the Docker image. >> "%LOG_FILE%"
goto :eof

:no_custom_package
REM No custom package found
set "CUSTOM_PACKAGE_FLAG=0"
echo No custom package (pyproject.toml) detected in src_install\custom_packages.
echo No custom package (pyproject.toml) detected in src_install\custom_packages. >> "%LOG_FILE%"
goto :eof