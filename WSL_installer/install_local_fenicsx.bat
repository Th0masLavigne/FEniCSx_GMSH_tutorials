@echo off
REM ============================================
REM Semi-interactive Ubuntu 24.04 + FEniCSx 0.9.0 Installer for WSL
REM Shows progress bar during installation
REM ============================================

setlocal
set "DIST_NAME=FEniCSxenv"
set "UBUNTU_DIST=Ubuntu-24.04"

REM Ask user for WSL username
set /p USER_NAME="Enter the WSL username you created during first-time setup: "
if "%USER_NAME%"=="" (
    echo ERROR: You must provide the username created in WSL.
    exit /b
)

REM ============================================
REM Step 1: Ensure WSL installed and WSL 2 default
REM ============================================
echo Checking WSL installation...
wsl.exe --status 
if %errorlevel% neq 0 (
    echo WSL not found. Enabling WSL and Virtual Machine Platform...
    dism.exe /online /enable-feature /featurename:Microsoft-Windows-Subsystem-Linux /all /norestart >nul
    dism.exe /online /enable-feature /featurename:VirtualMachinePlatform /all /norestart >nul
    echo WSL features enabled. Reboot required.
    pause
    exit /b
)
wsl --set-default-version 2 

REM ============================================
REM Step 2: Install Ubuntu 24.04 as FEniCSxenv
REM ============================================
echo Installing Ubuntu 24.04...
echo Follow instructions in the Ubuntu window to create your username and password.
echo Once done, type "exit" to return here.
wsl --install -d %UBUNTU_DIST% --name %DIST_NAME%
if %errorlevel% neq 0 (
    echo Ubuntu installation may already exist. Continuing...
)

REM ============================================
REM Step 3: Update Ubuntu system
REM ============================================
echo Updating WSL system...
wsl -d %DIST_NAME% -u %USER_NAME% -- bash -c "sudo apt update && sudo apt -y upgrade && sudo apt -y autoremove"

REM ============================================
REM Step 4: Install FEniCSx 0.9.0
REM ============================================
echo Installing FEniCSx 0.9.0...
wsl -d %DIST_NAME% -u %USER_NAME% -- bash -c "dpkg -l | grep -q fenicsx || (sudo apt install -y software-properties-common && sudo add-apt-repository -y ppa:fenics-packages/fenics && sudo apt update && sudo apt install -y 'fenicsx=2:0.9.0.1*' && sudo apt-mark hold fenicsx python3-dolfinx python3-ufl)"

REM ============================================
REM Step 5: Verify installation
REM ============================================
echo make the test
for /f %%v in ('wsl -d %DIST_NAME% -u %USER_NAME% -- python3 -c "import dolfinx; print(dolfinx.__version__)"') do set "FENICSX_VER=%%v"

if "%FENICSX_VER%"=="" (
    powershell -Command "Add-Type -AssemblyName PresentationFramework; [System.Windows.MessageBox]::Show('FEniCSx verification FAILED!','FEniCSx Installer',[System.Windows.MessageBoxButton]::OK,[System.Windows.MessageBoxImage]::Error)"
) else (
    powershell -Command "Add-Type -AssemblyName PresentationFramework; [System.Windows.MessageBox]::Show('FEniCSx successfully installed! Version: %FENICSX_VER%','FEniCSx Installer',[System.Windows.MessageBoxButton]::OK,[System.Windows.MessageBoxImage]::Information)"
)

echo ============================================
echo Installation/resume complete!
echo Launch environment using:
echo     wsl -d %DIST_NAME%
echo You will be automatically logged in as: %USER_NAME%
echo ============================================
pause
endlocal
