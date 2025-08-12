#!/bin/bash

# Folder pypackages present?
if [ ! -d pypackages ]; then

	echo "First time opening file. Initial setup may take around 5-10 minutes..."
	read -p "Would you like to install local PSIPRED? (y/n) " ans

	# if user select y
	if [ "$ans" = "y" ] || [ "$ans" = "Y" ]; then

		#INSTALL tcsh shell
		if ! command -v tcsh >/dev/null 2>&1; then
			# install tcsh
			if command -v apt >/dev/null 2>&1; then
				echo "Searching for tcsh shell..."
                                sudo apt update
				sudo apt install tcsh
				echo "tcsh shell successfully installed."
			else
				echo "Error installing tcsh..."
				exit 1
			fi
		else
			echo "tcsh is available."
		fi


		#DOWNLOAD NCBI DB
		echo "NCBI Protein Databases: https://ftp.ncbi.nlm.nih.gov/blast/db/"
		read -p "Which protein database to download? Select from the link above (DEFAULT: swissprot. Enter to use): " ans2

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
		cd ../../../..

		# ALLOW EXECUTIONS
		chmod -R +x .
	else
		echo "PSIPRED is not installed."
	fi

        # INSTALL TRANALIGN
	echo "Checking python3 dependencies..."


	# INSTALL PYTHON3 AND PYTHON3-PIP IF UNAVAILABLE
	if ! command -v python3 >/dev/null 2>&1 || ! dpkg -s python3-venv >/dev/null 2>&1; then
		if command -v apt-get >/dev/null 2>&1; then
			echo "Installing python3..."
			sudo apt install python3

			if ! command -v pip3 >/dev/null 2>&1; then
				echo "Installing python3-pip..."
				sudo apt install python3-pip
			fi
			if ! dpkg -s python3-venv >/dev/null 2>&1; then
				echo "Installing python3-venv"
				sudo apt install python3-venv
			fi
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
echo "Opening Alignate.."
echo ""
echo "NOTE:"
echo "If you do not have local PSIPRED installed, PSIPRED web is available but, it may run for ~10 mins (Local <1 min)."
echo "To install local PSIPRED, see README.txt - (2) & (3)"


source pypackages/bin/activate
python3 support_files/alignate.py
