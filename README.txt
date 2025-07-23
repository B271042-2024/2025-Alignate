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

3. [OPTIONAL: Don't need to install to use website-based pripred] Install BLAST+ protein database in psipred folder (see https://ftp.ncbi.nlm.nih.gov/blast/db/ for ref)
% mkdir external_tools/psipred/BLAST+/blastdb
% cd external_tools/psipred/BLAST+/blastdb
% wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/swissprot.tar.gz		# note: example database. please install accordingly.
% tar -xzf swissprot.tar.gz						# note: change according to your database

4. Allow executions in all files
% chmod -R +x .

5. If you move this folder, please remove and re-install python packages if software failed to work.
