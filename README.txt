========================================================================
Alignate V1.0 was completed on 15/08/2025 by Alignate Team for MSc Dissertation.
========================================================================


INTRODUCTION
Optimized for use on Ubuntu/Debian-based Linux system.
This tool aligns amino acid and codon sequences and generates consensus, conservation percentage and secondary structure prediction. It provides a better group comparison with multiple group-specific details, allowing users to select a reference group to base the comparisons on.


------
To run this Alignate v1.0, please type ./run.sh. This tool requires some additional installation, please see INSTALLATION for details. Only 1 & 4 are required.
------



RUN PROGRAM
chmod -R +x .
./run.sh




### [NOT REQUIRED] ADDITIONAL NOTES: TO DOWNLOAD/INSTALL INDEPENDENTLY, FOLLOW BELOW.
1. Install python packages. While in the project_alignate folder, please execute as below.
% python3 -m venv pypackages
% source pypackages/bin/activate
% pip install -r requirements.txt

2. [OPTIONAL: Only to use PSIPRED locally. To use PSIPRED website, skip this step.] Install tcsh shell (Only if not exist). Please run,
% sudo apt install tcsh

3. [OPTIONAL: Only to use PSIPRED locally. To use PSIPRED website, skip this step. Please select accordingly when running the program (after pressing button Align/Align All)] Install BLAST+ protein database in psipred folder (see https://ftp.ncbi.nlm.nih.gov/blast/db/ for ref)
% mkdir external_tools/psipred/BLAST+/blastdb
% cd external_tools/psipred/BLAST+/blastdb
% wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/swissprot.tar.gz		# note: example database. please install accordingly.
% tar -xzf swissprot.tar.gz						# note: change according to your database
% chmod -R +x .





### IMPORTANT!
- If you move this folder, please remove and re-install python packages (1) if software failed to work.
- If you are in trouble, leave me an email at adriana971026@gmail.com
