WINDOWS
1. to run: double-click on run.bat
2. to install required python libraries if required (listed in requirements.txt) - Do this in CMD/PowerShell (not WSL)
cd "C:\Users\adribot\Documents\Adri_study\UOE\Classes\Sem_2\MSc_Dissertation\project_alignate"
call pypackages_windows\Scripts\activate.bat
pip install -r requirements.txt
# then double-click on run.bat

LINUX
1. to run: ./run.sh
2. to install required python libraries if required (listed in requirements.txt)
source pypackages/bin/activate
pip install PySide6
pip install biopython
pip install matplotlib
# then run ./run.sh
