========================================================================
Alignate V1.0 was completed on 01/08/2025 by Alignate Team.
========================================================================


INTRODUCTION
This tool aligns amino acid and codon sequences and generates consensus, conservation percentage and secondary structure prediction. It provides a better group comparison with multiple group-specific details, allowing users to select a reference group to base the comparisons on.


------
To run this Alignate v1.0, please type ./run.sh. This tool requires some additional installation, please see INSTALLATION for details.
------

RUN PROGRAM
./run.sh


INSTALLATION
1. Install python packages. While in the project_alignate folder, please execute as below.
% python3 -m venv pypackages
% source pypackages/bin/activate
% pip install -r requirements.txt

2. Install tcsh shell (Only if not exist). Please run,
% sudo apt install tcsh

3. Install external tools in external_folder
1. psipred v + BLAST+
2. Place the folder in external_tools/ and rename the psipred folder to psipred














# Packaging sofware
pyinstaller alignate.py --add-data "external_tools/psipred:psipred"
# Have to run pyinstaller separately for each platform
## Eg. linux:
( please confirm if alignate.py or run.sh/run.bat ) pyinstaller alignate.py --onefile \
  --add-data "external_tools/clustalo_linux/clustalo:external_tools/clustalo" \
  --add-data "external_tools/mafft_linux/mafft:external_tools/mafft" \
  --add-data "external_tools/psipred:external_tools/psipred"
