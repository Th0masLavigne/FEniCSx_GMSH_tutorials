@echo off
setlocal enabledelayedexpansion

REM ====================================================================
REM DANGER ZONE: Uninstall FEniCSxenv WSL distribution
REM WARNING: This action is irreversible and will delete ALL data
REM          inside the FEniCSxenv distribution, including user files
REM          and the copied installation source code.
REM ====================================================================

set "DIST_NAME=FEniCSxenv"
set "LOG_FILE=report_uninstall.log"

echo Uninstallation process started on %date% at %time% >> "%LOG_FILE%"
echo -------------------------------------------------------------------- >> "%LOG_FILE%"

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

echo ====================================================================
echo WARNING: Preparing to UNREGISTER and DELETE the %DIST_NAME% environment.
echo ====================================================================
echo.
echo THIS ACTION IS IRREVERSIBLE. ALL FILES, including your documents,
echo settings, and the installed Docker image, will be permanently lost.
echo.

:: === MODIFICATION: Enhanced Pre-Unregister Pop-up Security Question ===
powershell -Command "Add-Type -AssemblyName PresentationFramework; $choice = [System.Windows.MessageBox]::Show('FINAL WARNING: You are about to permanently delete the WSL distribution %DIST_NAME%. ALL DATA (user files, Docker image, copied source) inside this environment will be lost. Click YES to proceed to the final console confirmation, NO to abort.','IRREVERSIBLE WSL DELETION WARNING','YesNo','Stop'); [IO.File]::WriteAllText('ConfirmWSLChoice.txt', $choice.ToString().Trim(), [System.Text.Encoding]::ASCII)"

set "WSL_CONFIRM="
for /f "usebackq tokens=1" %%i in ("ConfirmWSLChoice.txt") do set "WSL_CONFIRM=%%i"
del ConfirmWSLChoice.txt

if /i "!WSL_CONFIRM!"=="Yes" goto :CONFIRM_DELETION_CONSOLE
echo Graphical uninstallation aborted by user.
echo ============================================================>&2
echo UNINSTALLATION ABORTED.
echo ============================================================>&2
pause
exit /b 1


:CONFIRM_DELETION_CONSOLE
set /p USER_CONFIRM="Type 'YES' (in capital letters) to proceed with irreversible WSL deletion: "

if /i "%USER_CONFIRM%"=="YES" (
    echo Confirmed by user. Proceeding with unregistration... >> "%LOG_FILE%"
    goto :UNREGISTER_DISTRIBUTION
) else (
    goto :UNINSTALL_ABORTED
)

:: --- Subroutine for Unregistering the Distribution ---
:UNREGISTER_DISTRIBUTION
    REM Unregister the WSL distribution
    wsl --unregister %DIST_NAME%

    if %errorlevel% equ 0 (
        goto :UNREGISTER_SUCCESS
    ) else (
        goto :UNREGISTER_FAILURE
    )

:: --- Subroutine for Successful Unregistration ---
:UNREGISTER_SUCCESS
    echo %DIST_NAME% unregistered successfully. >> "%LOG_FILE%"

    :: === MODIFICATION: Post-Unregister Windows File Deletion Prompt ===
    echo ============================================================
    echo WSL UNREGISTRATION (unregistration) COMPLETE.
    goto :END_SCRIPT :: Redirect to the end label instead of :EOF

:: --- Subroutine for Failed Unregistration ---
:UNREGISTER_FAILURE
    echo ERROR: Failed to unregister %DIST_NAME%. It may not exist. >> "%LOG_FILE%"
    echo ============================================================
    echo UNINSTALLATION ERROR.
    echo ============================================================
    goto :END_SCRIPT :: Redirect to the end label instead of :EOF

:: --- Subroutine for Aborted Uninstallation ---
:UNINSTALL_ABORTED
    echo Uninstallation aborted by user. No changes were made to %DIST_NAME%. >> "%LOG_FILE%"
    echo ============================================================
    echo UNINSTALLATION ABORTED.
    echo ============================================================
    goto :END_SCRIPT :: Redirect to the end label instead of :EOF

:: --- End of Script Execution ---
:END_SCRIPT
pause
endlocal