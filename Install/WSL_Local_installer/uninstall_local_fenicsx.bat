@echo off
REM ============================================
REM Uninstall FEniCSxenv WSL distribution and Desktop shortcut
REM Fully consistent with the install_fenicsx.bat modifications
REM ============================================

setlocal
set "DIST_NAME=FEniCSxenv"
wsl --unregister %DIST_NAME%
echo %DIST_NAME% unregistered successfully.


echo ============================================
echo Uninstallation complete.
echo All FEniCSxenv files and settings removed.
echo ============================================
pause
endlocal
