@echo off
setlocal enabledelayedexpansion

REM Путь к первому дереву, BEAST (!)
set "TREE1=C:\Users\gdvov\OneDrive\Desktop\TREES_RENAMED_SHORTEST\BAYESIAN_REGIONAL_NO_LOGS\all_serotypes_3368-4090_strict_expop_100_combined.tree"

REM Путь к папке с .treefile
set "TREE2_DIR=C:\Users\gdvov\OneDrive\Desktop\TREES_RENAMED_SHORTEST\BAYESIAN_REGIONAL_NO_LOGS"

REM Перебрать все .treefile в папке + python-скрипт
for %%f in (%TREE2_DIR%\*.tree) do (
    echo Запуск для файла %%f
    python .\get_RF_halflife.py -tree1 %TREE1% -tree2 "%%f" -method all -thr 0.9 -thr2 0.9
)

pause