========================================================================
Alignate V1.0 was completed on 01/08/2025 by Alignate Team.
========================================================================


INTRODUCTION
This tool aligns amino acid and codon sequences and generates consensus, conservation percentage and secondary structure prediction. It provides a better group comparison with multiple group-specific details, allowing users to select a reference group to base the comparisons on.


INSTALLATION
1. Install python packages. While in the project alignate folder, please execute as below.
% python3 -m venv pypackages
% source pypackages/bin/activate
% pip install -r requirements.txt

2. Install tcsh shell. Please run,
% sudo apt install tcsh


3. Install external tools in external_folder
1. psipred v + BLAST+
2. Place the folder in external_tools/ and rename the psipred folder to psipred


TO RUN,
./run.sh
















# Packaging sofware
pyinstaller alignate.py --add-data "external_tools/psipred:psipred"
# Have to run pyinstaller separately for each platform
## Eg. linux:
( please confirm if alignate.py or run.sh/run.bat ) pyinstaller alignate.py --onefile \
  --add-data "external_tools/clustalo_linux/clustalo:external_tools/clustalo" \
  --add-data "external_tools/mafft_linux/mafft:external_tools/mafft" \
  --add-data "external_tools/psipred:external_tools/psipred"


LINUX
1. to run: ./run.sh
2. to install required python libraries if required (listed in requirements.txt)
python -m venv pypackages
source pypackages/bin/activate
pip install PySide6
pip install biopython
pip install matplotlib
# then run ./run.sh

WINDOWS
1. to run: double-click on run.bat
2. to install required python libraries if required (listed in requirements.txt) - Do this in CMD/PowerShell (not WSL)
cd "C:\Users\adribot\Documents\Adri_study\UOE\Classes\Sem_2\MSc Dissertation\project_alignate"
call pypackages_windows\Scripts\activate.bat
pip install -r requirements.txt
# then double-click on run.bat
