@echo off
setlocal

REM ============================================================================
REM FinancePy API Documentation Generator
REM ============================================================================

REM This batch file lives in:
REM     financepy-git\scripts\
set "SCRIPT_DIR=%~dp0"

REM Repository root:
REM     financepy-git\
for %%I in ("%SCRIPT_DIR%..") do set "PROJECT_ROOT=%%~fI"

REM Generated documentation:
REM     financepy-git\docs\index.html
set "INDEX_FILE=%PROJECT_ROOT%\html\index.html"

echo.
echo ============================================================================
echo Generating FinancePy API Documentation
echo ============================================================================
echo.
echo Repository:
echo     %PROJECT_ROOT%
echo.

cd /d "%PROJECT_ROOT%"

python "%PROJECT_ROOT%\scripts\generate_api_docs.py"

if errorlevel 1 (
    echo.
    echo ============================================================================
    echo ERROR: Documentation generation failed.
    echo ============================================================================
    echo.
    pause
    exit /b 1
)

echo.
echo ============================================================================
echo Documentation generated successfully.
echo ============================================================================
echo.

if exist "%INDEX_FILE%" (
    echo Opening:
    echo     %INDEX_FILE%
    echo.
    start "" "%INDEX_FILE%"
) else (
    echo ERROR: Generated index.html was not found.
    echo.
    echo Expected:
    echo     %INDEX_FILE%
    echo.
    pause
    exit /b 1
)

endlocal