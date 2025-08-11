#!/bin/bash
set -euo pipefail	# e: exit immediately, u: treat unset var as error, o pipefail: fail if any command fails, not just the last one

# Folder pypackages present?
if [ ! -d pypackages ]; then
	echo "First time opening file. Initial setup may take around 5 minutes..."
	read -p "Would you like to install local PSIPRED? (y/n) " ans

	# if user select y
	if [ "$ans" = "y" ] || [ "$ans" = "Y" ]; then

		#DOWNLOAD NCBI DB
		echo "NCBI Protein Databases: https://ftp.ncbi.nlm.nih.gov/blast/db/"
		read -p "Which protein database to download? Select from the link above (DEFAULT: swissprot). Enter to use default. " ans2

                mkdir -p external_tools/psipred/BLAST+/blastdb
                cd external_tools/psipred/BLAST+/blastdb || exit 1

		if [ -n "$ans2" ]; then
			if ! wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/"$ans2".tar.gz; then
				echo "Download failed. Try downloading different database."
				exit 1
			else
				tar -xzf "$ans2".tar.gz
			fi
		else
			if ! wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/swissprot.tar.gz; then
				echo "Download failed. Try downloading different database."
				exit 1
			else
				tar -xzf swissprot.tar.gz
			fi
		fi

		# INSTALL tcsh
		if ! command -v tcsh >/dev/null 2>&1; then
			# install tcsh
			if command -v apt >/dev/null 2>&1; then
				echo "Searching for tcsh shell..."
				sudo apt update && sudo apt install -y tcsh
			else
				echo "Error installing tcsh..."
				exit 1
			fi
		else
			echo "tcsh is available."
		fi

		# ALLOW EXECUTIONS
		chmod -R +x .
	else
		echo "PSIPRED is not installed."
	fi

	# INSTALL PYTHON3 AND PYTHON3-PIP IF UNAVAILABLE
	if ! command -v python3 >/dev/null 2>&1 || ! dpkg -s python3-venv >/dev/null 2>&1; then
		if command -v apt-get >/dev/null 2>&1; then
			echo "Installing python3 and python3-venv..."
			sudo apt-get update && sudo apt-get install -y python3 python3-pip python3-venv
		else
			echo "Error installing python3. Please install directly on terminal."
			exit 1
		fi
	fi
	# INSTALL PY LIBRARIES
	python3 -m venv pypackages
	source pypackages/bin/activate
	pip install -r requirements.txt
fi

# RUN
source pypackages/bin/activate
python3 support_files/alignate.py
