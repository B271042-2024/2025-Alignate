import sys, re, subprocess, tempfile, os, shutil, json, zipfile, tempfile, platform, uuid, requests, time, warnings
from PySide6.QtGui import QPixmap
from PySide6.QtWidgets import QToolTip, QGroupBox, QRadioButton, QProgressBar, QFileDialog, QMessageBox, QDialog, QTextEdit, QDialogButtonBox, QLayout, QScrollArea, QSizePolicy, QApplication, QMainWindow, QWidget, QCheckBox, QLabel, QLineEdit, QPushButton, QHBoxLayout, QVBoxLayout, QSlider
from PySide6.QtCore import QSize, Qt, QPoint, Signal, QThread
from Bio import SeqIO, AlignIO, BiopythonDeprecationWarning
from Bio.Align import MultipleSeqAlignment, AlignInfo
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrow
import matplotlib.colors as mcolors
from io import BytesIO
from collections import Counter

from drawingcanvas import DrawingCanvas
from ruler import ruler, ClickableLabel

warnings.filterwarnings("ignore", category=BiopythonDeprecationWarning)




class AlignmentWorker(QThread):
    finished = Signal(list)
    error = Signal(str)
    aligned_aa_signal = Signal(str)

    def __init__(self, uid, sequences, button_aln, parent=None):
        super().__init__(parent)
        self.sequences = sequences
        self.button_aln = button_aln
        self.output_file = None
        self.uid = uid

    def run(self):
        aligned_seq = []

        # Write FASTA file for alignment (Feed in sequences)
        base_path = os.path.dirname(os.path.abspath(__file__))
        output_folder = os.path.join(base_path, '..', 'output_files')
        os.makedirs(output_folder, exist_ok=True)

        # 2 create specific output folder
        session_folder = os.path.join(output_folder, self.uid)
        os.makedirs(session_folder, exist_ok=True)        

        # 3 create files
        fasta_file = os.path.join(output_folder, f"{self.uid}.fasta")
        aa_output_file = os.path.join(output_folder, f"{self.uid}_aa.fasta")
        self.aligned_aa_output_file = os.path.join(output_folder, f"{self.uid}_aa.aln")
        output_file = os.path.join(output_folder, f"{self.uid}.aln")
        
        with open(fasta_file, 'w') as fasta_out:
            for name, seq in self.sequences:
                if not name or not seq:
                    self.error.emit('Missing sequence name or sequence.')
                    #QMessageBox.warning(self, 'Error', 'Missing sequence name or sequence.')
                    return
                if not re.match(r'^[ACGTURYSWKMBDHVN\-]+$', seq.upper()):
                    self.error.emit(f'Invalid characters found in sequence: {name}. **For codon sequences only. Special characters are not allowed.')
                    #QMessageBox.warning(None, 'Error', f'Invalid characters found in sequence: {seq}. **For codon sequences only. Special characters are not allowed.')
                    return
                if '|' in name:
                    clean_id = name.split('|')[-1].split()[0]       # [-1] begin from the end
                    fasta_out.write(f'>{clean_id}\n{seq}\n')
                else:
                    fasta_out.write(f'>{name}\n{seq}\n')

# -------1 Translate codon sequences using biopython (output_output_files)
        translated_records = []
        for record in SeqIO.parse(fasta_file, "fasta"):
            aa_seq = record.seq.translate(to_stop=True)
            if '|' in record.id:
                clean_id = record.id.split('|')[-1].split()[0]                  # [-1] begin from the end
                translated_record = SeqRecord(aa_seq, id=clean_id, description="")
            else:
                translated_record = SeqRecord(aa_seq, id=record.id, description="")
            translated_records.append(translated_record)
        SeqIO.write(translated_records, aa_output_file, "fasta")

# -------2 Run Alignment: AA Sequences via mafft or clustalO
        try:
            # MAFFT 
            if self.button_aln.text() == 'MAFFT':
                print('#  Running MAFFT v7.526	(parameters: --anysymbol, --genafpair,--maxiterate 10000)...')
                mafft_path = "./external_tools/mafft_linux/mafft"
                with open(self.aligned_aa_output_file, 'w') as out:
                    subprocess.run([mafft_path, '--anysymbol', '--genafpair', '--maxiterate', '10000', aa_output_file], check=True, stdout=out, stderr=subprocess.DEVNULL)
            # ClustalO
            elif self.button_aln.text() == 'ClustalO':
                print('#  Running ClustalO v1.2.2 (parameters: default)...')
                clustalo_path = "./external_tools/clustalo_linux/clustalo"
                with open(self.aligned_aa_output_file, 'w') as out:
                    subprocess.run([clustalo_path, '-i', aa_output_file, '-o', self.aligned_aa_output_file, '--dealign', '--force'], check=True, stdout=out, stderr=subprocess.DEVNULL)
            # Neither
            else:
                self.error.emit('Error: Failed to select alignment method. Please contact Alignate.')
                #QMessageBox.warning(self, 'Error', 'Error: Failed to select alignment method. Please contact Alignate.')
                return
        except Exception as e:
            self.error.emit(f'{self.button_aln.text()} error: {str(e)}')
            #QMessageBox.critical(self, f'{self.button_aln.text()} error', str(e))
            return


# -------3 Run Alignment: Codon Sequences via tranalign
        tranalign_path = "./external_tools/tranalign/bin/tranalign"
        if not os.path.exists(self.aligned_aa_output_file) or os.path.getsize(self.aligned_aa_output_file) == 0:
            self.error.emit("Aligned amino acid output is missing or empty.")
            #QMessageBox.warning(self, "Error", "Aligned amino acid output is missing or empty.")
            return

        try:
            with open(output_file, 'w') as out:
                subprocess.run([tranalign_path, '-bsequence' , self.aligned_aa_output_file, '-asequence', fasta_file, '-outseq',  output_file], check=True, stdout=out, stderr=subprocess.DEVNULL)
        except Exception as e:
            self.error.emit(f'Tranalign alignment failed: {str(e)}')
            #QMessageBox.critical(self, "Tranalign alignment failed.", str(e))
            return
# ------4 Parse file and get aligned sequences
        for record in SeqIO.parse(output_file, 'fasta'):
            aligned_seq.append((record.id, str(record.seq)))                # extract name and seq for aligned sequences

        # emit signals
        self.aligned_aa_signal.emit(self.aligned_aa_output_file)
        self.finished.emit(aligned_seq)




class AlignmentDialog(QDialog):
    def __init__(self, uid, sequences, button_aln, callback_fn,
                 group=None, layout=None, output_file=None, return_only=False, seq_map=None,
                 parent=None):
        super().__init__(parent)
        self.setWindowTitle("Processing...")
        self.setModal(True)
        self.layout = QVBoxLayout(self)

        self.callback_fn = callback_fn
        self.group = group
        self.layout_param = layout
        self.output_file = output_file
        self.return_only = return_only
        self.seq_map = seq_map
        self.sequences = sequences
        self.uid = uid

        self.progress_bar = QProgressBar()
        self.progress_bar.setRange(0, 0)
        self.layout.addWidget(self.progress_bar)

        self.cancel_button = QPushButton("Cancel")
        self.cancel_button.clicked.connect(self.reject)
        self.layout.addWidget(self.cancel_button)

        self.worker = AlignmentWorker(self.uid, sequences, button_aln)
        self.worker.aligned_aa_signal.connect(self.receive_aa_aligned_file)
        self.worker.finished.connect(self.on_finished)
        self.worker.error.connect(self.show_error)
        self.worker.start()

    def receive_aa_aligned_file(self, aa_path: str):
        self.aa_aligned_path = aa_path

    def on_finished(self, aligned):
        self.accept()
        try:
            self.callback_fn(uid=self.uid,
                sequences=self.sequences,
                aligned_seq=aligned,
                group=self.group,
                layout=self.layout_param,
                button_aln=None,  # already used, optional here
                output_file=self.aa_aligned_path,
                return_only=self.return_only,
                seq_map=self.seq_map
            )
        except Exception as e:
            QMessageBox.critical(self, "Callback Error", str(e))

    def show_error(self, message):
        QMessageBox.critical(self, "Alignment Error", message)
        self.reject()





#_______________________________________________________________________________________________
#_______________________________________________________________________________________________
#_______________________________________________________________________________________________2 CLASS: Toolbar - codon
#_______________________________________________________________________________________________
#_______________________________________________________________________________________________





class codon(QWidget):
    def __init__(self, base_path, parent=None):
        super().__init__(parent)
        self.base_path = base_path

# --------------------------------------------Others
        # Initiation
        self.cancelled = False                                          # PROGRESS (DIALOG BOX)
        self.seq_rows = []                                              # SEQUENCE ROWS BY GROUP
        self.groups = []                                                # SEQUENCE DETAILS DICT
        self.widget_toggles = []                                        # FOR MENU2_HIDE TOGGLES
        self.is_alignall = False
        self.seq_map = None
        self.similarity_color = 'darkmagenta'

# --------------------------------------------Main
        # ---Layer 1
        # start window self
        layout_codon_l1 = QVBoxLayout()
        self.setLayout(layout_codon_l1)
        layout_codon_l1.setContentsMargins(0,0,0,0)

        # ---Layer 2
        widget_codon_l2 = QScrollArea()
        widget_codon_l2.setWidgetResizable(True)
        widget_codon_l2.setHorizontalScrollBarPolicy(Qt.ScrollBarAsNeeded)
        widget_codon_l2.setVerticalScrollBarPolicy(Qt.ScrollBarAsNeeded)
         # ---Add widgets to layout
        layout_codon_l1.addWidget(widget_codon_l2)

        # ---Layer 3
        widget_codon_l3 = QWidget()
        widget_codon_l3.setSizePolicy(QSizePolicy.MinimumExpanding, QSizePolicy.MinimumExpanding)
        self.layout_codon_l3 = QVBoxLayout()
        self.layout_codon_l3.setSpacing(0)
        self.layout_codon_l3.setContentsMargins(5,5,5,5)
        self.layout_codon_l3.setAlignment(Qt.AlignTop)
        widget_codon_l3.setLayout(self.layout_codon_l3)
         # ---Add widgets to parent widget
        widget_codon_l2.setWidget(widget_codon_l3)

# --------------------------------------------Sub-elements
        # ---Drawing Canvas (In Layer 3)
        self.canvas = DrawingCanvas()                                           # CONNECT TO FILE 2: drawingcanvas.py
        self.canvas.setToolTip('Left-click to draw, Right-click to erase.')
        self.canvas.setFixedHeight(40)
        self.canvas.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        self.layout_codon_l3.addWidget(self.canvas)

        # ---Ruler (In Layer 3)
        self.ruler = ruler(parent_context=self)                                        # CONNECT TO FILE 1: ruler.py
        self.layout_codon_l3.addWidget(self.ruler, alignment=Qt.AlignLeft)

        # ---Layer 4
        self.widget_codon_l4 = QWidget()
        self.layout_codon_l4 = QVBoxLayout()
        self.layout_codon_l4.setSpacing(0)
        self.layout_codon_l4.setContentsMargins(0, 0, 0, 0)
        self.widget_codon_l4.setLayout(self.layout_codon_l4)
        self.layout_codon_l3.addWidget(self.widget_codon_l4, alignment=Qt.AlignTop)

        # ---Drawing Canvas (In Layer 4)
#        self.canvas = DrawingCanvas()                                           # CONNECT TO FILE 2: drawingcanvas.py
#        self.canvas.setToolTip('Left-click to draw, Right-click to erase.')
#        self.canvas.setFixedHeight(40)
#        self.canvas.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
#        self.layout_codon_l4.addWidget(self.canvas)

        # ---Layer 4 Elements
        # Line 1
        # # 1 Button 1
        self.button1_addgroup = QPushButton('+Group')
        self.button1_addgroup.setFixedSize(QSize(60,28))
        self.button1_addgroup.setToolTip('Click to add more group(s).')
        self.button1_addgroup.setStyleSheet(
            """QPushButton {background-color: #00008B; color: white; font-weight: bold;}
            QPushButton:hover {background-color: #6495ED; color: white;}"""
        )
        # # 2 Button 2
        self.button2_alignall = QPushButton('Align all')
        self.button2_alignall.setFixedSize(QSize(60,28))
        self.button2_alignall.setToolTip(f'Click to align all sequences.\nPlease ensure that the reference group is correctly selected.')
        self.button2_alignall.setStyleSheet(
            """QPushButton {background-color: #00008B; color: white; font-weight: bold;}
            QPushButton:hover {background-color: #6495ED; color: white;}"""     
        )
        # # 3 Slider
        self.slidercon = QSlider(Qt.Horizontal)
        self.slidercon.setToolTip(f'Show %Conservation for each residue position: \nTick the right checkbox and move slider.')
        self.slidercon.setValue(100)
        self.slidercon.setMinimum(10)
        self.slidercon.setMaximum(100)
        self.slidercon.setSingleStep(10)
        self.slidercon.setTickPosition(QSlider.TicksBelow)
        self.slidercon.setFixedSize(120,20)
        self.slidercon.valueChanged.connect(self.slider_threshold)
        self.checkboxslider = QCheckBox()
        self.checkboxslider.setToolTip('Tick to activate slider.')
        self.checkboxslider.setChecked(False)
        self.checkboxslider.toggled.connect(self.handle_slider_mode_toggle)
        labelslider = QLabel('Tick to activate %Conservation slider')

        # # 4 Status
        self.status = QLabel('')

        # # Add widgets to parent widget
        self.widget_codon_buttons = QWidget()
        self.layout_codon_buttons = QHBoxLayout()
        self.layout_codon_buttons.setContentsMargins(0,2,0,2)
        self.widget_codon_buttons.setLayout(self.layout_codon_buttons)
        self.layout_codon_buttons.addWidget(self.button1_addgroup, alignment=Qt.AlignLeft)
        self.layout_codon_buttons.addWidget(self.button2_alignall, alignment=Qt.AlignLeft)
        self.layout_codon_buttons.addWidget(self.slidercon, alignment=Qt.AlignLeft)
        self.layout_codon_buttons.addWidget(self.checkboxslider, alignment=Qt.AlignLeft)
        self.layout_codon_buttons.addWidget(labelslider, alignment=Qt.AlignLeft)
        self.layout_codon_buttons.addWidget(self.status, alignment=Qt.AlignLeft)
        self.layout_codon_buttons.addStretch(1)
        self.layout_codon_l4.addWidget(self.widget_codon_buttons, alignment=Qt.AlignTop)

        # Line 2
        self.widget_codon_l4_2ndarystructure = QWidget()
        self.layout_codon_l4_2ndarystructure = QHBoxLayout()
        self.layout_codon_l4_2ndarystructure.setContentsMargins(5,0,0,0)
        self.layout_codon_l4_2ndarystructure.setSpacing(0)
        self.widget_codon_l4_2ndarystructure.setLayout(self.layout_codon_l4_2ndarystructure)
        self.widget_codon_l4_2ndarystructure.setObjectName("line_0")
        self.widget_codon_l4_2ndarystructure.setStyleSheet(
            """ #line_0 {
                padding: 2px;
            } """
        )    
        self.layout_codon_l4.addWidget(self.widget_codon_l4_2ndarystructure)

# --------------------------------------------Connect
        self.button1_addgroup.clicked.connect(self.button1_addgroup_clicked)
        self.button2_alignall.clicked.connect(self.button2_alignall_clicked)
        self.widget_toggles.append(self.widget_codon_buttons)

#_______________________________________________________________________________________________MAIN DEF___
#_______________________________________________________________________________________________1 DEF: QC System TCSH
#_______________________________________________________________________________________________

    def show_tcsh_warning(self):
        msg = QMessageBox(self)
        msg.setIcon(QMessageBox.Warning)
        msg.setWindowTitle('tcsh Not Found')
        msg.setText("Warning: Missing 'tcsh' shell. \nSome features like PSIPRED will not work. Please install using sudo apt install.")
        msg.setStandardButtons(QMessageBox.Ok)
        msg.show()


#_______________________________________________________________________________________________2 DEF: Set color for alignment
#_______________________________________________________________________________________________

    def set_similarity_color(self, color):
        self.similarity_color = color
        # Reapply coloring based on current alignment mode
        if getattr(self, 'is_alignall', False) and hasattr(self, 'seq_map'):
            self.color_code_seq(seq_map=self.seq_map, mode="all")
        else:
            for idx, group in enumerate(self.groups):
                self.color_code_seq(mode="group", group_idx=idx)

#_______________________________________________________________________________________________3 DEF: Search Sequence
#_______________________________________________________________________________________________

    def search_sequences(self, widget_search_seq):
        # 1 Define the target sequence
        search = widget_search_seq.strip().upper()

        if not search:                                                                              # if search is empty
            if self.is_searchseq:                                                                   # self.is_searchseq = True
               for group in self.groups:
                    for entry in group['widget_seq']:
                        for lbl in entry['seq_letters']:
                            bg_color = lbl.property("bg_color")
                            if bg_color is None:
                                bg_color = ""
                            lbl.setStyleSheet(f'background-color: {bg_color};')
            else:
                for group in self.groups:
                    for entry1 in group['widget_seq']:
                        for label in entry1['seq_letters']:
                            label.setStyleSheet("")  # clear search color
                self.is_searchseq = False
            return
        
        self.is_searchseq = True                                                                   # if not empty, self.is_searchseq = True

        for group in self.groups:
            for entry1 in group['widget_seq']:
                aligned_seq = ''.join(entry1['seq']).upper()
                label_row = entry1['seq_letters']

                # Reset previous highlights
                for label in label_row:
                    label.setStyleSheet("")

                # Remove gaps for logical match comparison, but keep index map
                non_gap_seq = ""
                pos_map = []  # index in aligned_seq → index in non_gap_seq
                for idx, char in enumerate(aligned_seq):
                    if char != "-":
                        non_gap_seq += char
                        pos_map.append(idx)

                # Search match in non-gap sequence
                for i in range(len(non_gap_seq) - len(search) + 1):
                    if non_gap_seq[i:i + len(search)] == search:
                        # Map match back to gapped sequence and highlight
                        for j in range(len(search)):
                            align_idx = pos_map[i + j]
                            label_row[align_idx].setStyleSheet("background: #FFA500;")


#_______________________________________________________________________________________________2-1 DEF: MENU2 VIEW
#_______________________________________________________________________________________________

    # 1 SHOW ALL
    def view_show_all(self):
        for widget in self.widget_toggles:
            widget.show()


    # 2 HIDE TOGGLES
    def view_hide_toggles(self):
        for widget in self.widget_toggles:
            widget.hide()


    # 3 ADJUST CONSENSUS MODE
    def adjust_consensusmode(self):
# --------------------------------------------Main
        dialog = QDialog(self)
        dialog.setStyleSheet('QDialog {background-color: white;}')
        dialog.move(500,500)
        dialog.setWindowTitle('Set Consensus Threshold')
        layout = QVBoxLayout()
        dialog.setLayout(layout)

# --------------------------------------------Sub-elements
        # ---1 label
        label = QLabel('Threshold (0.0 - 1.0):')
        layout.addWidget(label)

        # ---2 slider
        self.slider = QSlider(Qt.Horizontal)
        self.slider.setRange(0,100)
        self.slider.setValue(int(getattr(self, 'consensus_threshold', 1.0) * 100))
        layout.addWidget(self.slider)

        # ---3 threshold value
        value_label = QLabel(f"{self.slider.value()/100:.1f}")                   # :=start formatting, .1=1 decimal pts, f=fixed
        layout.addWidget(value_label)
        self.slider.valueChanged.connect(lambda v: value_label.setText(f"{v/100:.1f}"))

        # ---4 buttons
        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        layout.addWidget(buttons)
        buttons.accepted.connect(dialog.accept)                             # accept -> execution
        buttons.rejected.connect(dialog.reject)                             # cancel

# --------------------------------------------Connect
        if dialog.exec() == QDialog.Accepted:                               # execution
            threshold = self.slider.value() / 100.0
            self.consensus_threshold = threshold
            return threshold
        return None

#_______________________________________________________________________________________________2-2 DEF: MENU2 VIEW
#_______________________________________________________________________________________________

    def apply_new_consensus_threshold(self):
# . . .  CLEAR GUI . . . # . . . CLEAR DICT . . .
        for group in self.groups:
            layout = group['main_layout_seq']
            widget_to_remove = []
            target_names = ['consensus_row', 'conservation_block', 'custom_conservation_block', 'lbl_second', 'consensus_codon_row']
            for i in reversed(range(layout.count())):
                widget = layout.itemAt(i).widget()
                if widget and widget.objectName() in target_names:      # see in parent level
                    widget_to_remove.append(widget)
                else:                                                   
                    for name in target_names:
                        nested = widget.findChild(QWidget, name)        # else, see in children level
                        if nested:
                            widget_to_remove.append(nested)

            for widget in widget_to_remove:
                layout.removeWidget(widget)                             # remove from GUI
                widget.setParent(None)                                  # remove from hierarchy
                widget.deleteLater()                                    # schedule to delete later

# --------------------------------------------Connect
        # ---1 DEF: Adjust consensus mode
        threshold = self.adjust_consensusmode()

        # ---2 DEF: Update global consensus based on new threshold
        if threshold is not None:
            self.get_global_consensus(threshold=threshold, aa_aln_file=self.aligned_aa_output_file)

        # ---3 DEF: Update group consensus based on new threshold
            for group in self.groups:
                self.get_consensus_aln(group, threshold=threshold)

        # ---4 not DEF: Calculate % Base conservation
        # 1 get reference consensus
        ref_consensus = None
        for group in self.groups:
            if group['checkbox_setrefgroup'].isChecked():                                   
                ref_consensus = group.get('consensus_seq')
                break
        # 2 get group consensus and calculate %conservation
        for group in self.groups:
            layout = group['main_layout_seq']                                                   
            target_consensus = group.get('consensus_seq')
            if not target_consensus or target_consensus == ref_consensus:
                continue    # skip the for loop
            matches = sum(1 for a, b in zip(ref_consensus, target_consensus) if a == b)
            percent_conservation = (matches / len(ref_consensus)) * 100
            str_percent_conservation = f"{percent_conservation:.3g}%"        # 3sf
        # 3 update %conservation consensus against reference consensus
            for i in range(layout.count()):
                widget = layout.itemAt(i).widget()
                if widget and widget.objectName() == "consensus_row":
                    label = widget.findChild(QLabel, "percent_conservation")    # group['layout_seq'] > consensus_row > percent_conservation
                    if label:
                        label.setText(str_percent_conservation)

        # ---5 DEF: Update %conservation by region
        if hasattr(self, 'prediction_text'):
            self.compute_region_conservation(self.prediction_text)

#_______________________________________________________________________________________________3-1 DEF: Line 1 - Slider + Slider Checkbox
#_______________________________________________________________________________________________ (Called in DEF: init)

    def handle_slider_mode_toggle(self, checked):                                               # Turn the slider fxn ON/OFF                                           
        if checked:
            self.slider_threshold()
        else:
            if hasattr(self, 'seq_map') and self.seq_map:
                self.color_code_seq(seq_map=self.seq_map, mode="all")
            else:
                for idx, group in enumerate(self.groups):
                    self.color_code_seq(mode="group", group_idx=idx)


#_______________________________________________________________________________________________3-2 DEF: Slider
#_______________________________________________________________________________________________ (Called in DEF: run_alignment)

    def slider_threshold(self):
        if self.checkboxslider.isChecked():
            threshold = self.slidercon.value()
            lower = threshold - 5.5
            upper = threshold + 5.5

            all_seqs = []
            for group in self.groups:
                for entry in group['widget_seq']:
                    all_seqs.append(entry['seq'])

            if not all_seqs:
                return

            columns = list(zip(*all_seqs))
            conservation_scores = []
            for col in columns:
                most_common = Counter(col).most_common(1)[0][1]
                percent = round((most_common / len(col)) * 100)
                conservation_scores.append(percent)

            for group in self.groups:
                for entry in group['widget_seq']:
                    for i, lbl in enumerate(entry['seq_letters']):
                        if i >= len(conservation_scores):
                            continue
                        percent = conservation_scores[i]
                        bg_color = lbl.property("bg_color") or 'white'    # self.get_existing_bg_color(lbl)
                        if lower <= percent <= upper:
                            lbl.setStyleSheet(f'background-color: {bg_color}; color: black;')
                        else:
                            lbl.setStyleSheet('color: lightgray;')

        else:
            if getattr(self, 'is_alignall', False) and hasattr(self, 'seq_map'):
                self.color_code_seq(seq_map=self.seq_map, mode="all")
            else:
                for idx, group in enumerate(self.groups):
                    self.color_code_seq(mode="group", group_idx=idx)


#_______________________________________________________________________________________________4 DEF: Menu - Save Project
#_______________________________________________________________________________________________

    def save_project(self):
# --------------------------------------------Main
        # ---1 DIALOGBOX: SET FILE PATH (SAVE)
        file_path, _ = QFileDialog.getSaveFileName(self, 'Save Project', '', 'Alignate Project(*.alignate)')
        if not file_path:                                   # If no file, return
            return
        if not file_path.endswith('.alignate'):             # Automatically append extension if missing
            file_path += '.alignate'

        # ---2-1 MAIN STORAGE: STATE
        state = {
            'groups': [],
            'slider_value': self.slidercon.value(),
            'slider_checked': self.checkboxslider.isChecked()
        }
        for group in self.groups:
        # ---2-2 MAIN STORAGE: STATE < GROUP_DATA < [GROUP NAME, REFERENCE GROUP, SEQUENCES]
            group_data = {
                'name': group['lineedit_groupname'].text(),
                'is_reference': group['checkbox_setrefgroup'].isChecked(),
                'sequences': []
            }
        # ---2-3 MAIN STORAGE: STATE < GROUP_DATA < SEQUENCES NAME, REFERENCE GROUP, SEQUENCES: [SEQUENCE NAMES, SEQUENCES]
            for entry in group['widget_seq']:
                header = entry['seq_header'].text()
                seq = ''.join(label.text() for label in entry['seq_letters'])
                group_data['sequences'].append((header, seq))
            state['groups'].append(group_data)

        # ---3-1 OUTPUT FILE: JSON
        with tempfile.TemporaryDirectory() as temp_dir:         # Create temporary dir that is automatically deleted after with block ends
            json_path = os.path.join(temp_dir, 'state.json')    # Set path for file: state.json
            with open(json_path, 'w') as f:                         
                json.dump(state, f)                             # Open state.json to write and dump state dict as JSON
        # ---3-2 OUTPUT FILE: ZIP THE JSON FILE
            with zipfile.ZipFile(file_path, 'w') as zf:
                zf.write(json_path, 'state.json')               # Write state.json into Zip file, with the name: state.json

# --------------------------------------------Others
        # ---4 MESSAGE BOX: SAVED
        QMessageBox.information(self, 'Saved', f'Project saved to {file_path}')


#_______________________________________________________________________________________________5 DEF: Menu - Load Project
#_______________________________________________________________________________________________

    def load_project(self):
# --------------------------------------------Main
        # ---1 DIALOGBOX: OPEN FILE PATH
        file_path, _ = QFileDialog.getOpenFileName(self, 'Load Project', '', 'Alignate Project (*.alignate)')
        if not file_path:
            return
        
        # ---2-1 EXTRACT ZIP FILE
        with tempfile.TemporaryDirectory() as temp_dir:
            with zipfile.ZipFile(file_path, 'r') as zf:         
                zf.extractall(temp_dir)                                         # Extract Zip file into Temporary Directory: temp_dir
        # ---2-2 LOAD JSON FILE
            with open(os.path.join(temp_dir, 'state.json')) as f:
                state = json.load(f)                                            # Load json file

        # ---3 RESET WIDGET
            self.slidercon.setValue(state.get('slider_value', 100))                     # Reset slidercon
            self.checkboxslider.setChecked(state.get('slider_checked', False))          # Reset slider checked

# --------------------------------------------Connect
            # 1 Remove previous self.groups
            for group in self.groups[:]:                                                # Update Widget Group Data
                self.button2_removegroup_clicked(group['widget_group'])                 # ':' = create a copy of the self.groups

            for group_data in state["groups"]:
            # 2 Update self.groups with STATE (new)
                self.button1_addgroup_clicked()                                         # Both GUI and DICT
                group = self.groups[-1]                                                 # -1 = last item
                group['lineedit_groupname'].setText(group_data["name"])                 # Update: Group Name
                group['checkbox_setrefgroup'].setChecked(group_data["is_reference"])    # Update: Checkbox Set Reference Group
                for header, seq in group_data["sequences"]:
            # 3 Update GUI
                    self.add_sequences_toGUI(group, group["layout_seq"], header, seq)   # Add sequences to GUI

# --------------------------------------------Others
        # 4 MESSAGE BOX: LOADED
        QMessageBox.information(self, "Loaded", f"Project loaded from {file_path}")


#_______________________________________________________________________________________________OTHER DEF___
#_______________________________________________________________________________________________4-1 DEF: Line 1 - +Group
#_______________________________________________________________________________________________

    def button1_addgroup_clicked(self):
# . . .  CLEAR GUI . . .
        if hasattr(self, 'widget_global') and self.widget_global is not None:
            self.layout_codon_l3.removeWidget(self.widget_global)             # only remove from layout
            self.widget_global.setParent(None)                                  # detach from parent widget
            self.widget_global.deleteLater()                                    # schedule for deletion to free memory safely
            self.widget_global = None                                           # clear subject

        if hasattr(self, 'widget_global_codon') and self.widget_global_codon is not None:
            self.layout_codon_l3.removeWidget(self.widget_global_codon)
            self.widget_global_codon.setParent(None)
            self.widget_global_codon.deleteLater()
            self.widget_global_codon = None

# --------------------------------------------Main
        # MAIN WIDGET: GROUP
        self.widget_codon_l4_group_l1 = QWidget()
        self.widget_codon_l4_group_l1.setObjectName('group')         
        self.widget_codon_l4_group_l1.setStyleSheet(
            """ #group {
                border: 1px solid #A9A9A9;
                border-radius: 6px;
                padding: 2px;
            }"""
        )
        self.layout_codon_l4_group_l1 = QVBoxLayout()
        self.layout_codon_l4_group_l1.setContentsMargins(0,2,0,2)
        self.layout_codon_l4_group_l1.setSpacing(0)                                               
        self.widget_codon_l4_group_l1.setLayout(self.layout_codon_l4_group_l1)
        self.layout_codon_l4.addWidget(self.widget_codon_l4_group_l1)



# --------------------------------------------Sub-elements
        # ---1 LINE 1 (BUTTONS)
        # MAIN WIDGET: LINE 1
        widget_codon_l4_group_l1_line1 = QWidget()
        layout_codon_l4_group_l1_line1 = QHBoxLayout()
        widget_codon_l4_group_l1_line1.setLayout(layout_codon_l4_group_l1_line1)
        self.layout_codon_l4_group_l1.addWidget(widget_codon_l4_group_l1_line1, alignment=Qt.AlignTop)

        # ELEMENTS
        # 1 create elements
        totalgroups = len(self.groups)
        groupno = totalgroups + 1 if totalgroups > 0 else 1
        checkbox1_setrefgroup = QCheckBox()
        if len(self.groups) == 0:
            checkbox1_setrefgroup.setChecked(True)
        label2_setrefgroup = QLabel('Tick to set reference group (Default: Group 1)')

        lineedit_groupname = QLineEdit('Group'+str(groupno))
        lineedit_groupname.setFixedWidth(120)
        lineedit_groupname.setToolTip('Suggestion: Use scientific name!')
        button2_removegroup = QPushButton('-Group')
        button2_removegroup.setFixedWidth(60)
        button2_removegroup.setToolTip('Click to remove group.')
        button3_addseq = QPushButton('+')
        button3_addseq.setFixedWidth(30)
        button3_addseq.setToolTip('Click to add sequence(s) to this group.')
        button4_removeseq = QPushButton('-')
        button4_removeseq.setFixedWidth(30)
        button4_removeseq.setToolTip('Tick the row(s) below to mark for deletion and click to remove.')
        button5_align = QPushButton('Align')
        button5_align.setFixedWidth(60)
        button5_align.setToolTip('Click to align sequences in this group only.')
        button5_align.setStyleSheet("""
            QPushButton {background-color: #00008B; color: white; font-weight: bold;}
            QPushButton:hover {background-color: #6495ED; color: white;}""")

        # 2 add elements to layout
        layout_codon_l4_group_l1_line1.addWidget(lineedit_groupname, alignment=Qt.AlignLeft)
        layout_codon_l4_group_l1_line1.addWidget(button2_removegroup, alignment=Qt.AlignLeft)
        layout_codon_l4_group_l1_line1.addWidget(button3_addseq, alignment=Qt.AlignLeft)
        layout_codon_l4_group_l1_line1.addWidget(button4_removeseq, alignment=Qt.AlignLeft)
        layout_codon_l4_group_l1_line1.addWidget(button5_align, alignment=Qt.AlignLeft)
        layout_codon_l4_group_l1_line1.addWidget(checkbox1_setrefgroup, alignment=Qt.AlignLeft)
        layout_codon_l4_group_l1_line1.addWidget(label2_setrefgroup, alignment=Qt.AlignLeft)
        layout_codon_l4_group_l1_line1.addStretch(0)                      # to left-align all elements
        
        # ---2 LINE 2 (SEQUENCES)
        # MAIN WIDGET: SEQUENCES
                # MAIN WIDGET: SEQUENCES
        self.widget_codon_l4_group_l1_seq0 = QWidget()
        self.layout_codon_l4_group_l1_seq0 = QVBoxLayout()
        self.widget_codon_l4_group_l1_seq0.setFixedHeight(160)
        self.layout_codon_l4_group_l1_seq0.setSpacing(0)
        self.widget_codon_l4_group_l1_seq0.setLayout(self.layout_codon_l4_group_l1_seq0)
        self.layout_codon_l4_group_l1_seq0.setContentsMargins(5,0,0,0)
        self.layout_codon_l4_group_l1.addWidget(self.widget_codon_l4_group_l1_seq0)

        widget_codon_l4_group_l1_seq1 = QScrollArea()
        widget_codon_l4_group_l1_seq1.setWidgetResizable(True)
        widget_codon_l4_group_l1_seq1.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        widget_codon_l4_group_l1_seq1.setVerticalScrollBarPolicy(Qt.ScrollBarAsNeeded)
        self.layout_codon_l4_group_l1_seq0.addWidget(widget_codon_l4_group_l1_seq1)

        self.widget_codon_l4_group_l1_seq = QWidget()
        layout_codon_l4_group_l1_seq = QVBoxLayout()
        self.widget_codon_l4_group_l1_seq.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Minimum)
        layout_codon_l4_group_l1_seq.setSpacing(0)
        layout_codon_l4_group_l1_seq.setAlignment(Qt.AlignTop)
        self.widget_codon_l4_group_l1_seq.setLayout(layout_codon_l4_group_l1_seq)
        layout_codon_l4_group_l1_seq.setContentsMargins(5,0,0,0)
        widget_codon_l4_group_l1_seq1.setWidget(self.widget_codon_l4_group_l1_seq)   

        # INITIATE DICTIONARY
        group = {
            'widget_group': self.widget_codon_l4_group_l1,          # WIDGET GROUP
            'lineedit_groupname': lineedit_groupname,               # GROUP NAME
            'main_layout_seq': self.layout_codon_l4_group_l1_seq0,  # MAIN LAYOUT FOR LAYOUT SEQUENCE
            'layout_seq': layout_codon_l4_group_l1_seq,             # LAYOUT SEQUENCE
            'widget_seq': [],                                       # WIDGET SEQUENCE
            'checkbox_setrefgroup': checkbox1_setrefgroup           # CHECKBOX REFERENCE GROUP
        }

# --------------------------------------------Connect
        self.groups.append(group)
        lineedit_groupname.textChanged.connect(lambda text: group.update({'name': text}))
        self.widget_toggles.append(widget_codon_l4_group_l1_line1)
        button2_removegroup.clicked.connect(lambda _=None, w=self.widget_codon_l4_group_l1: self.button2_removegroup_clicked(w))
        button3_addseq.clicked.connect(lambda _=None, layout=layout_codon_l4_group_l1_seq: self.button3_addseq_clicked(layout))
        button4_removeseq.clicked.connect(lambda _=None, layout=self.layout_codon_l4_group_l1: self.button4_removeseq_clicked(layout))
        button5_align.clicked.connect(lambda _=None, layout=layout_codon_l4_group_l1_seq: self.button5_align_clicked(layout))
        checkbox1_setrefgroup.toggled.connect(lambda checked, this_box=checkbox1_setrefgroup: self.handle_reference_group_toggle(this_box))                                                 # only 1 group is allowed at a time


#_______________________________________________________________________________________________4-2 DEF: Line 1 - Exclusive Reference Checkbox
#_______________________________________________________________________________________________

    def handle_reference_group_toggle(self, selected_checkbox):
        if selected_checkbox.isChecked():
            for group in self.groups:
                box = group['checkbox_setrefgroup']
                if box != selected_checkbox:
                    box.setChecked(False)


#_______________________________________________________________________________________________5 DEF: Line 1 - Remove Group
#_______________________________________________________________________________________________

    def button2_removegroup_clicked(self, widget):
# . . .  CLEAR GUI . . .
        # 1 Remove global consensus
        if hasattr(self, 'widget_global') and self.widget_global is not None:
            self.layout_codon_l3.removeWidget(self.widget_global)     # remove from layout (optional)
            self.widget_global.setParent(None)                          # remove from parent GUI hierarchy
            self.widget_global.deleteLater()                            # schedule for safe deletion by Qt event loop
            self.widget_global = None                                   # no widget global

        if hasattr(self, 'widget_global_codon') and self.widget_global_codon is not None:
            self.layout_codon_l3.removeWidget(self.widget_global_codon)
            self.widget_global_codon.setParent(None)
            self.widget_global_codon.deleteLater()
            self.widget_global_codon = None

        # 2 Remove Widget Group
        widget.setParent(None)
        widget.deleteLater()

# . . . CLEAR DICT . . .
        # 3 If Widget Group is in dictionary groups, remove all related details
        for i, group in enumerate(self.groups):
            if group['widget_group'] == widget:
                del self.groups[i]
                break


#_______________________________________________________________________________________________6-1 DEF: Line 1 - Add Sequence
#_______________________________________________________________________________________________

    def button3_addseq_clicked(self, layout):
# . . .  CLEAR GUI . . . # . . . CLEAR DICT . . .
        for group in self.groups:
            if group['main_layout_seq'] == layout:
                widgets_to_remove = []
                target_names = ['consensus_row', 'conservation_block', 'custom_conservation_block', 'lbl_second', 'consensus_codon_row']
                for i in reversed(range(layout.count())):
                    widget = layout.itemAt(i).widget()
                    if widget and widget.objectName() in target_names:
                        widgets_to_remove.append(widget)
                    else:
                        for name in target_names:
                            nested = widget.findChild(QWidget, name)
                            if nested:
                                widgets_to_remove.append(nested)

                    for widget in widgets_to_remove:
                        layout.removeWidget(widget)
                        widget.setParent(None)
                        widget.deleteLater()

        if hasattr(self, 'widget_global') and self.widget_global is not None:
            self.layout_codon_l3.removeWidget(self.widget_global)
            self.widget_global.setParent(None)
            self.widget_global.deleteLater()
            self.widget_global = None

        if hasattr(self, 'widget_global_codon') and self.widget_global_codon is not None:
            self.layout_codon_l3.removeWidget(self.widget_global_codon)
            self.widget_global_codon.setParent(None)
            self.widget_global_codon.deleteLater()
            self.widget_global_codon = None

# --------------------------------------------Main
        self.widget_codon_14_group_l1_seq_dialoginput = QDialog(self)
        self.widget_codon_14_group_l1_seq_dialoginput.setStyleSheet('QDialog {background-color: white;}')
        self.widget_codon_14_group_l1_seq_dialoginput.move(500,500)
        self.layout_codon_l4_group_l1_seq_dialoginput = QVBoxLayout()
        self.widget_codon_14_group_l1_seq_dialoginput.setLayout(self.layout_codon_l4_group_l1_seq_dialoginput)

# --------------------------------------------Sub-elements
        # 1 label
        widget_codon_l4_group_l1_seq_dialoginput_label = QLabel()
        widget_codon_l4_group_l1_seq_dialoginput_label.setText('Please upload FASTA file containing codon sequences:')
        self.layout_codon_l4_group_l1_seq_dialoginput.addWidget(widget_codon_l4_group_l1_seq_dialoginput_label)

        # 2 buttons
        self.seqtext_button1file = QPushButton('Upload File')
        self.seqtext_button2cancel = QPushButton('Cancel')
        self.widget_codon_l4_group_l1_seq_dialoginput_seqbutton = QDialogButtonBox()
        self.widget_codon_l4_group_l1_seq_dialoginput_seqbutton.addButton(self.seqtext_button1file, QDialogButtonBox.ActionRole)
        self.widget_codon_l4_group_l1_seq_dialoginput_seqbutton.addButton(self.seqtext_button2cancel, QDialogButtonBox.RejectRole)
        self.layout_codon_l4_group_l1_seq_dialoginput.addWidget(self.widget_codon_l4_group_l1_seq_dialoginput_seqbutton, alignment=Qt.AlignCenter)

# --------------------------------------------Connect
        self.seqtext_button1file.clicked.connect(lambda _=None, layout=layout: self.seqtext_button1file_clicked(layout))
        self.seqtext_button2cancel.clicked.connect(self.widget_codon_14_group_l1_seq_dialoginput.reject)

# --------------------------------------------Others
        # Execution
        self.widget_codon_14_group_l1_seq_dialoginput.exec()




#_______________________________________________________________________________________________6-3 DEF: Line 1 - Add Codon & codon Files
#_______________________________________________________________________________________________

    def seqtext_button1file_clicked(self, layout):

# --------------------------------------------Initiation
        self.allinput_header = []
        self.allinput_seq = []
        self.dialog_seqlength = False

# --------------------------------------------Main
        # 1 File Upload Dialog Box
        file_path, _ = QFileDialog.getOpenFileName(self, 'Select File', '', 'FASTA files (*.fasta *.fa *.txt);;All Files (*)')
        if not file_path:
            return
        # Parse the uploaded file
        try:
            sequences = list(SeqIO.parse(file_path, 'fasta'))
            if not sequences:
                QMessageBox.warning(self, 'No sequences found', 'The uploaded file is empty or do not have sequences.')
                return
            for record in sequences:
                name = record.id or 'unnamed'
                seq = str(record.seq)

                if not re.match(r'^[ACGTURYSWKMBDHVN\-]+$', seq.upper()):
                    QMessageBox.warning(self, 'Invalid sequence', f'Invalid characters found in sequence: {name}. **For codon sequences only. Special characters are not allowed.')
                    return

                self.allinput_header.append(name)
                self.allinput_seq.append(seq)


    # --------------------------------------------QC self.allinput_seq/header
            # ---1 Get max length from the self.groups dict
            len_groupseq = []
            for group in self.groups:
                if group['widget_seq']:
                    for entry in group['widget_seq']:
                        seq = entry['seq']
                        len_groupseq.append(len(seq))                      # store length of all sequences
            max_len_groupseq = max(len_groupseq) if len_groupseq else 0
            min_len_groupseq = min(len_groupseq) if len_groupseq else 0

            # ---2 Get max length of the newly input sequences
            len_allinput_seq = []
            for seq in self.allinput_seq:
                if len(seq) > 900:
                    QMessageBox.warning(self, 'Sequence Length Exceeded', 'One or more of the uploaded sequences exceed the maximum allowed length of 900. Please remove them and re-upload.')
                    self.widget_codon_14_group_l1_seq_dialoginput.accept()
                    return
                len_allinput_seq.append(len(seq))
            max_len_allinput_seq = max(len_allinput_seq) if len_allinput_seq else 0

            # ---3 Compare 1 & 2 --> store the max out of 2
            max_len = max(max_len_groupseq, max_len_allinput_seq)

            # ---4 List down the list - possible to remove from allinput_seq/header using index
            list_index = []
            list_len = []
            list_header = []
            for len_index, length in enumerate(len_allinput_seq):
                if max_len - length > 10:
                    list_index.append(len_index)
                    list_len.append(length)
                    list_header.append(self.allinput_header[len_index])

            
            # ---5 Display with a checkbox for users to select and remove from the list
            # 1 main layout
            if list_index:
                widget_filter_seq_dialog = QDialog(self)
                widget_filter_seq_dialog.setStyleSheet('QDialog {background-color: white;}')
                widget_filter_seq_dialog.move(700,600)
                layout_filter_seq_dialog = QVBoxLayout()
                widget_filter_seq_dialog.setLayout(layout_filter_seq_dialog)
                widget_filter_seq_dialog.setWindowTitle('Filter Partial Sequences')

                if min_len_groupseq:
                    if max_len - min_len_groupseq > 10:
                        label_filter_min = QLabel(f'**Short sequence length detected in window (length={min_len_groupseq}). Please remove if required.')
                        layout_filter_seq_dialog.addWidget(label_filter_min)

                label_filter_max = QLabel(f'**Longest sequence: {max_len}')
                layout_filter_seq_dialog.addWidget(label_filter_max)

                label_filter_seq_dialog = QLabel('Please check sequences below and select to remove from the new list:')
                label_filter_seq_dialog.setStyleSheet('font-weight: bold;')
                layout_filter_seq_dialog.addWidget(label_filter_seq_dialog)

                # 2 main seq layout (lines of seqs)
                widget_filter_seq_scroll = QScrollArea()
                widget_filter_seq_scroll.setWidgetResizable(True)
                widget_filter_seq_scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAsNeeded)
                widget_filter_seq_scroll.setVerticalScrollBarPolicy(Qt.ScrollBarAsNeeded)
                layout_filter_seq_dialog.addWidget(widget_filter_seq_scroll)

                widget_filter_seq_dialog2 = QWidget()
                widget_filter_seq_dialog2.setStyleSheet('background-color: white;')
                layout_filter_seq_dialog2 = QVBoxLayout()
                widget_filter_seq_dialog2.setSizePolicy(QSizePolicy.MinimumExpanding, QSizePolicy.MinimumExpanding)
                widget_filter_seq_dialog2.setLayout(layout_filter_seq_dialog2)
                widget_filter_seq_scroll.setWidget(widget_filter_seq_dialog2)
    #            layout_filter_seq_dialog.addWidget(widget_filter_seq_dialog2)

                # 4 fill in each seq layout
                self.filter_checkboxes = []
                for i in range(len(list_index)):
                    checkbox_filter_seq_dialog3 = QCheckBox()  
                    label_filter = QLabel(f"{list_header[i]}    (length={list_len[i]})")

                    # 3 each seq layout (checkbox - header)
                    widget_filter_seq_dialog3 = QWidget()
                    layout_filter_seq_dialog3 = QHBoxLayout()
                    widget_filter_seq_dialog3.setLayout(layout_filter_seq_dialog3)
                    layout_filter_seq_dialog2.addWidget(widget_filter_seq_dialog3)

                    layout_filter_seq_dialog3.addWidget(checkbox_filter_seq_dialog3, alignment=Qt.AlignLeft)
                    layout_filter_seq_dialog3.addWidget(label_filter, alignment=Qt.AlignLeft)
                    layout_filter_seq_dialog3.addStretch(1)
                    self.filter_checkboxes.append((checkbox_filter_seq_dialog3, list_index[i]))

                # 5 buttons
                widget_buttons_filter = QWidget()
                layout_buttons_filter = QHBoxLayout()
                widget_buttons_filter.setLayout(layout_buttons_filter)
                layout_filter_seq_dialog.addWidget(widget_buttons_filter)
                button_filter_OK = QPushButton('OK')
                button_filter_cancel = QPushButton('Cancel')
                layout_buttons_filter.addWidget(button_filter_OK)
                layout_buttons_filter.addWidget(button_filter_cancel)
                layout_filter_seq_dialog.addStretch(1)

                # 6 if OK is clicked, remove all checked from all_inputheaader & all_inputseq
                button_filter_OK.clicked.connect(lambda: (self.filter_sequences(widget_filter_seq_dialog), setattr(self, 'dialog_seqlength', True)))
                button_filter_cancel.clicked.connect(widget_filter_seq_dialog.reject)

                widget_filter_seq_dialog.exec()
                if not self.dialog_seqlength:
                    return

    # --------------------------------------------Connect
            for seq_name, seq in zip(self.allinput_header, self.allinput_seq):
                for group in self.groups:
                    if group['layout_seq'] == layout:
                        self.add_sequences_toGUI(group, layout, seq_name, seq)
                        break

        except Exception as e:
            QMessageBox.critical(self, 'Error', f'Could not load file:\n{e}')

# --------------------------------------------Others
        self.widget_codon_14_group_l1_seq_dialoginput.accept()




#_______________________________________________________________________________________________6-4 DEF: Line 1 - Filter the added sequences before adding to GUI etc.
#_______________________________________________________________________________________________

    def filter_sequences(self, dialog):
        to_remove = [idx for checkbox, idx in self.filter_checkboxes if checkbox.isChecked()]
        for idx in sorted(to_remove, reverse=True):
            del self.allinput_header[idx]
            del self.allinput_seq[idx]

        dialog.accept()




#_______________________________________________________________________________________________7 DEF: Line 1 - Remove Sequence
#_______________________________________________________________________________________________

    def button4_removeseq_clicked(self, layout):
# . . .  CLEAR GUI . . .
        for group in self.groups:
            if group['main_layout_seq'] == layout:
                widgets_to_remove = []
                target_names = ['consensus_row', 'conservation_block', 'custom_conservation_block', 'lbl_second', 'consensus_codon_row']
                for i in reversed(range(layout.count())):
                    widget = layout.itemAt(i).widget()
                    if widget and widget.objectName() in target_names:
                        widgets_to_remove.append()
                    else:
                        for name in target_names:
                            nested = widget.findChild(QWidget, name)
                            if nested:
                                widgets_to_remove.append(nested)

                    for widget in widgets_to_remove:
                        layout.removeWidget(widget)
                        widget.setParent(None)
                        widget.deleteLater()

        if hasattr(self, 'widget_global') and self.widget_global is not None:
            self.layout_codon_l3.removeWidget(self.widget_global)
            self.widget_global.setParent(None)
            self.widget_global.deleteLater()
            self.widget_global = None

        if hasattr(self, 'widget_global_codon') and self.widget_global_codon is not None:
            self.layout_codon_l3.removeWidget(self.widget_global_codon)
            self.widget_global_codon.setParent(None)
            self.widget_global_codon.deleteLater()
            self.widget_global_codon = None


        seq_removed = []
        for row, checkbox in self.seq_rows:
            try:
                if checkbox.isChecked():
                    seq_removed.append((row, checkbox))
            except RuntimeError:
                continue    # in case that the widget has been deleted'

        for row, checkbox in seq_removed:
            row.setParent(None)
            row.deleteLater()

# . . . CLEAR DICT . . .
        rows_to_remove = set()
        for row, checkbox in seq_removed:
            try:
                rows_to_remove.add(row)
            except RuntimeError:
                continue

        # update group['widget_seq']
        for group in self.groups:
            if 'widget_seq' in group:
                group['widget_seq'] = [entry for entry in group['widget_seq'] if entry['widget_row'] not in rows_to_remove]

# --------------------------------------------Connect
        # update seq_rows to include only undeleted ones
        self.seq_rows = [(row, checkbox) for (row, checkbox) in self.seq_rows if (row, checkbox) not in seq_removed]    # used to mark for deletion


#_______________________________________________________________________________________________8 DEF: Line 1 - Align sequences
#__________________________________________________________________________________________ALIGN

    def button5_align_clicked(self, layout):
# --------------------------------------------Main
        widget_dialogbox_align = QDialog(self)
        widget_dialogbox_align.setStyleSheet('QDialog {background-color: white;}')
        widget_dialogbox_align.move(500,500)
        layout_dialogbox_align = QVBoxLayout()
        widget_dialogbox_align.setLayout(layout_dialogbox_align)

# --------------------------------------------Sub-elements
        # 1 Label
        widget_label_align = QLabel('ClustalO v1.2.2 [all platforms] (Settings: Default)\nMAFFT v7.526 [win & linux] v7.490 [mac] (Settings: --anysymbol, --genafpair, --maxiterate 10000')
        layout_dialogbox_align.addWidget(widget_label_align)

        # 2 Button
        widget_buttons_align = QWidget()
        layout_buttons_align = QHBoxLayout()
        widget_buttons_align.setLayout(layout_buttons_align)
        widget_button1_clustalo = QPushButton('ClustalO')
        widget_button2_mafft = QPushButton('MAFFT')
        widget_button3_cancel = QPushButton('Cancel')
        layout_buttons_align.addWidget(widget_button1_clustalo)
        layout_buttons_align.addWidget(widget_button2_mafft)
        layout_buttons_align.addWidget(widget_button3_cancel)
        layout_dialogbox_align.addWidget(widget_buttons_align)

        # 3 Radiobutton for PSIPRED (Run online or offline)
        widget_radiobtns_psipred = QGroupBox('PSIPRED Setting:')
        layout_radiobtns_psipred = QHBoxLayout()
        self.widget_psipred_offline = QRadioButton('Run Offline')
        self.widget_psipred_offline.setChecked(True)
        self.widget_psipred_online = QRadioButton('Run Online')
        layout_radiobtns_psipred.addWidget(self.widget_psipred_offline)
        layout_radiobtns_psipred.addWidget(self.widget_psipred_online)
        widget_radiobtns_psipred.setLayout(layout_radiobtns_psipred)
        layout_dialogbox_align.addWidget(widget_radiobtns_psipred)

        # 4 Visible only if run offline
        self.widget_radiobtns_psipred2 = QGroupBox('Run PSIPRED with/without BLAST+')
        layout_radiobtns_psipred2 = QHBoxLayout()
        self.widget_psipred_blast = QRadioButton('with BLAST+ (default)')
        self.widget_psipred_blast.setChecked(True)
        self.widget_psipred_single = QRadioButton('without BLAST+ (**less accurate)')
        layout_radiobtns_psipred2.addWidget(self.widget_psipred_blast)
        layout_radiobtns_psipred2.addWidget(self.widget_psipred_single)
        self.widget_radiobtns_psipred2.setLayout(layout_radiobtns_psipred2)
        layout_dialogbox_align.addWidget(self.widget_radiobtns_psipred2)

#        if self.widget_psipred_online.isChecked():
        self.qlineedit_email = QLineEdit()
        self.qlineedit_email.setPlaceholderText('Please enter your email')
        self.qlineedit_email.setVisible(False)                                      # defaule: not visible
        layout_dialogbox_align.addWidget(self.qlineedit_email)

# --------------------------------------------Connect
        # 1 if the group sequence layout is correct, extract group
        group = None
        for g in self.groups:
            if g['layout_seq'] == layout:
                group = g
                break

        widget_button3_cancel.clicked.connect(widget_dialogbox_align.reject)
        widget_button1_clustalo.clicked.connect(lambda _=None: (widget_dialogbox_align.accept(), self.run_alignment([(entry['seq_header'].text().strip(), ''.join(label.text() for label in entry['seq_letters']).strip()) for entry in group['widget_seq']], group, layout, widget_button1_clustalo)))
        widget_button2_mafft.clicked.connect(lambda _=None: (widget_dialogbox_align.accept(), self.run_alignment([(entry['seq_header'].text().strip(), ''.join(label.text() for label in entry['seq_letters']).strip()) for entry in group['widget_seq']], group, layout, widget_button2_mafft)))
        self.widget_psipred_offline.toggled.connect(self.toggle_psipred)
        self.widget_psipred_online.toggled.connect(self.toggle_psipred)

# --------------------------------------------Others
        widget_dialogbox_align.exec()


#_______________________________________________________________________________________________8-2 DEF: Line 1 - Align sequences
#__________________________________________________________________________________________ALIGN

    def toggle_psipred(self):
        if self.widget_psipred_online.isChecked():
            self.widget_radiobtns_psipred2.setVisible(False)
            self.qlineedit_email.setVisible(True)
        else:
            self.qlineedit_email.setVisible(False)
            self.widget_radiobtns_psipred2.setVisible(True)


#_______________________________________________________________________________________________9 DEF: 1 Align all sequences
#__________________________________________________________________________________________ALIGN
    def button2_alignall_clicked(self):
# --------------------------------------------Others
        # Initiation
        all_seq = []
        self.seq_map = []
        self.is_alignall = True

# --------------------------------------------Main
        widget_dialogbox_alignall = QDialog(self)
        widget_dialogbox_alignall.setStyleSheet('QDialog {background-color: white;}')
        widget_dialogbox_alignall.move(500,500)
        layout_dialogbox_alignall = QVBoxLayout()
        widget_dialogbox_alignall.setLayout(layout_dialogbox_alignall)

# --------------------------------------------Sub-elements
        # 1 Label
        widget_label_alignall = QLabel('ClustalO v1.2.2 [all platforms] (Settings: Default)\nMAFFT v7.526 [win & linux] v7.490 [mac] (Settings: --anysymbol, --genafpair, --maxiterate 10000')
        layout_dialogbox_alignall.addWidget(widget_label_alignall)

        # 2 Button
        widget_buttons_alignall = QDialogButtonBox()
        layout_buttons_alignall = QHBoxLayout()
        widget_buttons_alignall.setLayout(layout_buttons_alignall)
        widget_button1_clustalo = QPushButton('ClustalO')
        widget_button2_mafft = QPushButton('MAFFT')
        widget_button3_cancel = QPushButton('Cancel')
        widget_buttons_alignall.addButton(widget_button1_clustalo, QDialogButtonBox.ActionRole)
        widget_buttons_alignall.addButton(widget_button2_mafft, QDialogButtonBox.ActionRole)
        widget_buttons_alignall.addButton(widget_button3_cancel, QDialogButtonBox.RejectRole)
        layout_dialogbox_alignall.addWidget(widget_buttons_alignall)

        # 3 Radiobutton for PSIPRED (Run online or offline)
        widget_radiobtns_psipred = QGroupBox('PSIPRED Setting:')
        layout_radiobtns_psipred = QHBoxLayout()
        self.widget_psipred_offline = QRadioButton('Run Offline')
        self.widget_psipred_offline.setChecked(True)
        self.widget_psipred_online = QRadioButton('Run Online')
        layout_radiobtns_psipred.addWidget(self.widget_psipred_offline)
        layout_radiobtns_psipred.addWidget(self.widget_psipred_online)
        widget_radiobtns_psipred.setLayout(layout_radiobtns_psipred)
        layout_dialogbox_alignall.addWidget(widget_radiobtns_psipred)

        # 4 Visible only if run offline
        self.widget_radiobtns_psipred2 = QGroupBox('Run PSIPRED with/without BLAST+')
        layout_radiobtns_psipred2 = QHBoxLayout()
        self.widget_psipred_blast = QRadioButton('with BLAST+ (default)')
        self.widget_psipred_blast.setChecked(True)
        self.widget_psipred_single = QRadioButton('without BLAST+ (**less accurate)')
        layout_radiobtns_psipred2.addWidget(self.widget_psipred_blast)
        layout_radiobtns_psipred2.addWidget(self.widget_psipred_single)
        self.widget_radiobtns_psipred2.setLayout(layout_radiobtns_psipred2)
        layout_dialogbox_alignall.addWidget(self.widget_radiobtns_psipred2)

#        if self.widget_psipred_online.isChecked():
        self.qlineedit_email = QLineEdit()
        self.qlineedit_email.setPlaceholderText('Please enter your email')
        self.qlineedit_email.setVisible(False)                                      # defaule: not visible
        layout_dialogbox_alignall.addWidget(self.qlineedit_email)

# --------------------------------------------Action
        for group_idx, group in enumerate(self.groups):
            for entry_idx, entry in enumerate(group['widget_seq']):
                if not hasattr(entry['seq_header'], 'text'):
                    continue
                seq_header = entry.get('seq_header')
                if not seq_header or not seq_header.text():
                    continue
                name = entry['seq_header'].text().strip()                               # extract header
                seq = ''.join(label.text() for label in entry['seq_letters']).strip()   # extract sequence
                if name and seq:
                    all_seq.append((name, seq))                                         # extract all_seq (name, seq)
                    self.seq_map.append((group_idx, entry_idx))                         # extract seq_map (index)
                else:
                    QMessageBox.warning(self, 'Error', 'Missing sequence header or sequences.')
                    return
        if len(all_seq) < 2:
            QMessageBox.warning(self, 'Error', 'Need at least 2 sequences for alignment.')
            return

# --------------------------------------------Connect
        widget_button3_cancel.clicked.connect(widget_dialogbox_alignall.reject)                # cancel
        widget_button2_mafft.clicked.connect(lambda _=None: (widget_dialogbox_alignall.accept(), self.run_alignment(all_seq, None, None, widget_button2_mafft, return_only=False, seq_map=self.seq_map)))
        widget_button1_clustalo.clicked.connect(lambda _=None: (widget_dialogbox_alignall.accept(), self.run_alignment(all_seq, None, None, widget_button1_clustalo, return_only=False, seq_map=self.seq_map)))
        self.widget_psipred_offline.toggled.connect(self.toggle_psipred)
        self.widget_psipred_online.toggled.connect(self.toggle_psipred)

# --------------------------------------------Others
        widget_dialogbox_alignall.exec()



#_______________________________________________________________________________________________10-4 DEF: Run Alignment > etc.
#_______________________________________________________________________________________________

    def run_alignment(self, sequences, group=None, layout=None, button_aln=None, output_file=None, return_only=False, seq_map=None):
    
        print('')
        print('#1 Aligning amino acid sequence...')

        self.status.setText('Running . . .')
        self.status.setStyleSheet('color: red; font-weight: bold;')

    # --------------------------------------------Others
        # Initiation
        self.is_searchseq = False               # Retain the color after alignment in letter
        uid = uuid.uuid4().hex[:6]
        self.base_path = os.path.dirname(os.path.abspath(__file__))
        output_folder = os.path.join(self.base_path, '..', 'output_files')
        self.session_folder = os.path.join(output_folder, uid)

        # If no group is set as a reference
        any_checked = any(group['checkbox_setrefgroup'].isChecked() for group in self.groups)
        if not any_checked and self.groups:
            self.groups[0]['checkbox_setrefgroup'].setChecked(True)

        # QC 1: No. of sequences need to be > 2
        if len(sequences) < 2:
            QMessageBox.warning(self, 'Error', 'Need at least 2 sequences for alignment.')
            return

# --------------------------------------------Actions
# ------1 Run Alignment: MAFFT/ClustlO
        dlg = AlignmentDialog(uid, 
            sequences, button_aln, self.continue_run_alignment,
            group=group, layout=layout,
            output_file=output_file, return_only=return_only, seq_map=seq_map
        )
        dlg.exec()    



    def continue_run_alignment(self, uid, sequences, aligned_seq, group=None, layout=None, button_aln=None, output_file=None, return_only=False, seq_map=None):

        print('')
        print('#2 Generating consensus sequences and calculate %Conservation...')

        # ---Initiation
        self.uid = uid
        self.aligned_aa_output_file = output_file

        # ---1 Alignall
        if seq_map:                                                         # created in def button2_alignall_clicked(self)
            # -- 1 Clear GUI (unaligned seq/previously aligned seq) and DICT
            for group in self.groups:
                layout = group['layout_seq']
                for i in reversed(range(layout.count())):
                    widget = layout.itemAt(i).widget()
                    if widget:
                        widget.setParent(None)                              # REMOVE FROM GUI
                        widget.deleteLater()
                group['widget_seq'].clear()                                 # REMOVE FROM DICT: widget: sequence
                group['consensus_seq'] = None                               # REMOVE FROM DICT: consensus sequence 

            # -- 2 Connect - DEF: to add sequences to GUI
            for (group_idx, __), (aligned_name, aligned_seq) in zip(seq_map, aligned_seq):       
                group = self.groups[group_idx]                              # get each group
                layout = group['layout_seq']                                # get sequence layout
                if layout is None:
                    continue
                self.add_sequences_toGUI(group, layout, aligned_name, aligned_seq)

            print('#  Generate consensus sequences')

            # -- 3 Connect - DEF: Get and Display Consensus in each group
            for group in self.groups:
                consensus_str, codon_consensus_str = self.get_consensus_aln(group, threshold=None, button_aln=button_aln, sequences=sequences)

            print('#  Calculate %Conservation')

            # -- 4 not DEF: Calculate %Base conservation
            # 1 get reference consensus
            ref_consensus = None
            for group in self.groups:
                if group['checkbox_setrefgroup'].isChecked():                                   
                    ref_consensus = group.get('consensus_seq')
                    break
            # 2 get group consensus and calculate %conservation
            for group in self.groups:
                layout = group['main_layout_seq']                                                   
                target_consensus = group.get('consensus_seq')
                if not target_consensus or target_consensus == ref_consensus:
                    continue    # skip the for loop
                matches = sum(1 for a, b in zip(ref_consensus, target_consensus) if a == b)
                percent_conservation = (matches / len(ref_consensus)) * 100
                str_percent_conservation = f"{percent_conservation:.3g}%"        # 3sf
            # 3 update %conservation consensus against reference consensus
                for i in range(layout.count()):
                    widget = layout.itemAt(i).widget()
                    if widget and widget.objectName() == "consensus_row":
                        label = widget.findChild(QLabel, "percent_conservation")    # group['layout_seq'] > consensus_row > percent_conservation
                        if label:
                            label.setText(str_percent_conservation)

            # -- 5 Connect - DEF Color Code & Display on GUI (Aligned Sequences)
            self.color_code_seq(seq_map=self.seq_map, mode="all")


            # run global consensus to get both global consensus for codon and codon
            # -- 6 Connect - DEF Get & Display on GUI (Global Consensus)
            global_p_consensus_str = self.get_global_consensus(aa_aln_file=self.aligned_aa_output_file)
            if global_p_consensus_str:
                global_p_consensus_file = os.path.join(self.session_folder, f"{self.uid}_globalconsensus.fasta")
                with open(global_p_consensus_file, 'w') as f:
                    f.write(f">globalconsensus\n{global_p_consensus_str}\n")

                if not os.path.exists(global_p_consensus_file):
                    QMessageBox.critical(self, 'Error', f'Consensus FASTA file not found: {global_p_consensus_file}')
                    return
                elif os.path.getsize(global_p_consensus_file) == 0:
                    QMessageBox.critical(self, 'Error', f'Consensus FASTA file is empty: {global_p_consensus_file}')
                    return

                print('')
                print('#3 Generating secondary structure with PSIPRED...')

                if self.widget_psipred_offline.isChecked():
                    psipred_dir = os.path.join(self.base_path, '..', 'external_tools', 'psipred')
                    self.prediction_text = self.build_secondary_structure_offline(psipred_dir)
                else:
                    self.prediction_text = self.build_secondary_structure_online(global_p_consensus_file)

            # -- 8 Connect - DEF Display on GUI (PSIPRED Output)
                region_conservation = self.compute_region_conservation(self.prediction_text)
                self.draw_secondary_structure_to_gui(self.prediction_text, region_conservation)


        # ---2 Align
        else:
            if layout is not None:
            # -- 1 Clear GUI (unaligned seq/previously aligned seq) and DICT
                for i in reversed(range(layout.count())):
                    widget = layout.itemAt(i).widget()
                    if widget:
                        widget.setParent(None)                                          # removes from layout & GUI hierarchy but not delete
                        widget.deleteLater()                                            # REMOVE FROM GUI # ***** to amend ***** originally absent

            # -- 2 Connect - DEF: to add sequences to GUI
            if group is not None:                                                       # group is defined from def: button5_align_clicked
                group['widget_seq'].clear()                                             # REMOVE FROM DICT
                for name, seq in aligned_seq:
            # 1 Connect - DEF Add aligned seq on GUI
                    self.add_sequences_toGUI(group, layout, name, seq)                  

                print('#  Generate consensus sequences')

            # 2 Connect - DEF Add codon_consensus then, aa_consensus onto GUI
                consensus_str, codon_consensus_str = self.get_consensus_aln(group, threshold=None, button_aln=None)

            # -- 3 Connect - DEF Color Code & Display on GUI (Aligned Sequences)
                group_idx = self.groups.index(group)
                self.is_alignall = False
                self.color_code_seq(mode="group", group_idx=group_idx)

            # -- 4 Connect - DEF Run PSIPRED: Build Secondary Structure

                print('')
                print('#3 Generating secondary structure with PSIPRED...')

                if codon_consensus_str:
                    p_consensus_file = os.path.join(self.session_folder, f"{self.uid}_groupconsensus.fasta")
                    with open(p_consensus_file, 'w') as f:
                        f.write(f">consensus\n{codon_consensus_str}\n")

                    if not os.path.exists(p_consensus_file):
                        QMessageBox.critical(self, 'Error', f'Consensus FASTA file not found: {p_consensus_file}')
                        return
                    elif os.path.getsize(p_consensus_file) == 0:
                        QMessageBox.critical(self, 'Error', f'Consensus FASTA file is empty: {p_consensus_file}')
                        return

                    if self.widget_psipred_offline.isChecked():
                        psipred_dir = os.path.join(self.base_path, '..', 'external_tools', 'psipred')
                        self.prediction_text = self.build_secondary_structure_offline(psipred_dir)
                    else:
                        self.prediction_text = self.build_secondary_structure_online(p_consensus_file)

            # -- 5 Connect - DEF Display on GUI (PSIPRED Output)
                    self.draw_secondary_structure_to_gui(self.prediction_text)
                else:
                    QMessageBox.warning(self, 'Missing value', 'Consensus Sequences not found. Please contact Alignate.')
                    return

        print('')
        print('Finished.')

        self.status.setText('')

# --------------------------------------------Connect
        self.slider_threshold()

#_______________________________________________________________________________________________11 Add sequences to GUI (Unaligned & Aligned)
#_______________________________________________________________________________________________

    def add_sequences_toGUI(self, group, layout, seq_name='', seq=''):                 
# --------------------------------------------Main
        widget_seq = QWidget()
        layout_seq = QHBoxLayout()
        layout_seq.setContentsMargins(5,0,0,0)
        layout_seq.setSpacing(0)
        widget_seq.setLayout(layout_seq)
        layout.addWidget(widget_seq, alignment=Qt.AlignLeft)

# --------------------------------------------Sub-elements
        # 1 Checkbox
        seq_checkbox = QCheckBox()
        layout_seq.addWidget(seq_checkbox)
        layout_seq.addSpacing(5)

        # 2 Sequence name
        seq_header = QLineEdit(seq_name)
        seq_header.setFixedSize(120,20)
        layout_seq.addWidget(seq_header)
        layout_seq.addSpacing(5)
        layout_seq.setAlignment(Qt.AlignLeft)

        # 3 Sequence
        seq_letters = []
        for letter in seq:
            seq_letter = QLabel(letter)
            seq_letter.setFixedSize(15,20)
            seq_letter.setAlignment(Qt.AlignCenter)
            seq_letter.setStyleSheet('border: 1px solid #F8F8F8; padding: 1px;')
            layout_seq.addWidget(seq_letter)
            seq_letters.append(seq_letter)

# --------------------------------------------Connect
        group['widget_seq'].append({'seq_header': seq_header, 'seq': seq, 'seq_letters': seq_letters, 'widget_row': widget_seq})
        self.seq_rows.append((widget_seq, seq_checkbox))

# --------------------------------------------Others 
        # 4 Horizontal scrollbar: Dynamically resize parent container based on longest sequence
        avg_char_width = 8
        min_padding = 200
        seq_pixel_length = len(seq) * avg_char_width + min_padding
        current_width = self.widget_codon_l4.minimumWidth()
        if seq_pixel_length > current_width:
            self.widget_codon_l4.setMinimumWidth(seq_pixel_length)


#_______________________________________________________________________________________________12-1 Get & Display Consensus (by Group)
#_______________________________________________________________________________________________

    def get_consensus_aln(self, group, seq_map=None, threshold=None, button_aln=None, sequences=None):

        print('#  Generating group consensus') 

# . . .  CLEAR GUI . . .
        layout = group['main_layout_seq']
        if layout is None:
            return
        for i in reversed(range(layout.count())):
            widget = layout.itemAt(i).widget()
            if widget:
                name = widget.objectName()
                if name == "consensus_row":                             # Remove Group Consensus Row
                    layout.removeWidget(widget)
                    widget.setParent(None)
                    widget.deleteLater()

# --------------------------------------------Main
        widget_consensus = QWidget()
        widget_consensus.setObjectName('consensus_row')
        layout_consensus = QHBoxLayout()
        layout_consensus.setContentsMargins(5,0,0,0)
        layout_consensus.setSpacing(0)
        widget_consensus.setLayout(layout_consensus)
        layout.addWidget(widget_consensus, alignment=Qt.AlignLeft)       

# --------------------------------------------Sub-elements
        # 1 Checkbox (Only for spacing)
        invisible_checkbox = QCheckBox()
        invisible_checkbox.setEnabled(False)
        invisible_checkbox.setStyleSheet('background: transparent; border: none;')
        layout_consensus.addWidget(invisible_checkbox)

        # 2 Label
        label_consensus = QLabel('')
        label_consensus.setObjectName('percent_conservation')
        label_consensus.setFixedSize(120,20)
        layout_consensus.addWidget(label_consensus, alignment=Qt.AlignLeft)
        layout_consensus.addSpacing(5)  

# --------------------------------------------Action
        # 1 Get consensus from dict using SeqIO
        records = [SeqRecord(Seq(entry['seq']), id=entry['seq_header'].text()) for entry in group['widget_seq']]
        alignment = MultipleSeqAlignment(records)
        summary = AlignInfo.SummaryInfo(alignment)
        if threshold is None:
            threshold = getattr(self, 'consensus_threshold', 0.5)
        consensus = summary.dumb_consensus(threshold=threshold, ambiguous='N')
        consensus_str = str(consensus)

        for letter in consensus_str:
            lbl = QLabel(letter)
            lbl.setFixedSize(15,20)
            lbl.setAlignment(Qt.AlignCenter)
            lbl.setStyleSheet('color:gray;')
            layout_consensus.addWidget(lbl)


# --------------------------------------------Connect
        group['consensus_seq'] = consensus_str
        group_name = group['lineedit_groupname'].text()

# --------------------------------------------Output_files
        consensus_file = os.path.join(self.session_folder, f'{self.uid}_groupconsensus.fasta')
        with open(consensus_file, 'a') as c:
            c.write(f'>{group_name}\n{consensus_str}\n')



# --------------------------------------------
# --------------------------------------------Amino Acid Consensus: Translate codon consensus
        # Clean and chunk codon string into triplets
        codons = [consensus_str[i:i+3] for i in range(0, len(consensus_str), 3)]
        # Translate codons into amino acids
        codon_consensus_str = ""
        for codon in codons:
            if len(codon) < 3 or '-' in codon or 'N' in codon.upper():
                codon_consensus_str += 'X'  # Ambiguous or incomplete
            else:
                try:
                    aa = Seq(codon).translate()
                    codon_consensus_str += str(aa)
                except:
                    codon_consensus_str += 'X'


# --------------------------------------------Output_files
        codon_consensus_file = os.path.join(self.session_folder, f'{self.uid}_groupconsensus_codon.fasta')
        with open(codon_consensus_file, 'a') as c:
            c.write(f'>{group_name}\n{codon_consensus_str}\n')


                # --------------------------------------------Main_codon Consensus
        widget_consensus_codon = QWidget()
        widget_consensus_codon.setObjectName('consensus_codon_row')
        layout_consensus_codon = QHBoxLayout()
        layout_consensus_codon.setContentsMargins(5, 0, 0, 0)
        layout_consensus_codon.setSpacing(0)
        widget_consensus_codon.setLayout(layout_consensus_codon)
        layout.addWidget(widget_consensus_codon, alignment=Qt.AlignLeft)

        # Sub-elements
        invisible_checkbox_codon = QCheckBox()
        invisible_checkbox_codon.setEnabled(False)
        invisible_checkbox_codon.setStyleSheet('background: transparent; border: none;')
        layout_consensus_codon.addWidget(invisible_checkbox_codon)

        label_consensus_codon = QLabel('')
        label_consensus_codon.setFixedSize(120, 20)
        layout_consensus_codon.addWidget(label_consensus_codon, alignment=Qt.AlignLeft)
        layout_consensus_codon.addSpacing(5)

        for aa in codon_consensus_str.strip():
            lbl = QLabel(aa)
            lbl.setFixedSize(45, 20)
            lbl.setAlignment(Qt.AlignCenter)
            lbl.setStyleSheet('color:gray;')
            layout_consensus_codon.addWidget(lbl)
            group['consensus_seq_codon_codon'] = codon_consensus_str

# --------------------------------------------Others
        return consensus_str, codon_consensus_str



#_______________________________________________________________________________________________12-1 Get & Display Global Consensus
#_______________________________________________________________________________________________

    def get_global_consensus(self, threshold=None, aa_aln_file=None):

        print('#  Generating global consensus')

# . . .  CLEAR GUI . . .
        if hasattr(self, 'widget_global') and self.widget_global is not None:
            self.layout_codon_l3.removeWidget(self.widget_global)
            self.widget_global.setParent(None)
            self.widget_global.deleteLater()
            self.widget_global = None

        if hasattr(self, 'widget_global_codon') and self.widget_global_codon is not None:
            self.layout_codon_l3.removeWidget(self.widget_global_codon)
            self.widget_global_codon.setParent(None)
            self.widget_global_codon.deleteLater()
            self.widget_global_codon = None


# --------------------------------------------Main
        self.widget_global = QWidget()
        layout_global = QHBoxLayout()
        layout_global.setContentsMargins(5,0,0,0)
        layout_global.setSpacing(0)
        self.widget_global.setLayout(layout_global)
#        self.layout_codon_l3.addStretch(1)
        self.layout_codon_l3.addWidget(self.widget_global, alignment=Qt.AlignLeft)

# --------------------------------------------Sub-elements
        # 1 Checkbox (only for spacing)
        invisible_checkbox = QCheckBox()
        invisible_checkbox.setEnabled(False)
        invisible_checkbox.setStyleSheet('background: transparent; border: none;')
        layout_global.addWidget(invisible_checkbox)
        layout_global.addSpacing(5)

        # 2 Label
        invisible_label = QLabel('Consensus')
        invisible_label.setStyleSheet('font-weight: bold;')
        invisible_label.setFixedSize(120,20)
        layout_global.addWidget(invisible_label, alignment=Qt.AlignLeft)
        layout_global.addSpacing(5)

# --------------------------------------------Action
        all_records = []
        for group in self.groups:
            for entry in group['widget_seq']:
                name = entry['seq_header'].text().strip()
                seq = ''.join(label.text() for label in entry['seq_letters']).strip()
                all_records.append(SeqRecord(Seq(seq), id=name))
        if not all_records:
            return None
        alignment = MultipleSeqAlignment(all_records)
        summary = AlignInfo.SummaryInfo(alignment)
        if threshold is None:
            threshold = 0.5
        consensus = summary.dumb_consensus(threshold=threshold, ambiguous='N')
        consensus_str = str(consensus)

        for letter in str(consensus):
            lbl = QLabel(letter)
            lbl.setFixedSize(15,20)
            lbl.setAlignment(Qt.AlignCenter)
            lbl.setStyleSheet('color:gray;')
            layout_global.addWidget(lbl)



# --------------------------------------------Output_files
        global_consensus_file = os.path.join(self.session_folder, f'{self.uid}_globalconsensus.fasta')
        with open(global_consensus_file, 'a') as c:
            c.write(f'>global_consensus\n{consensus_str}\n')



# --------------------------------------------
# --------------------------------------------Amino Acid Consensus

        global_p_consensus_str = None  # Default to None for safe return
        if aa_aln_file and os.path.exists(aa_aln_file) and os.path.getsize(aa_aln_file) > 0:
            try:
                aa_alignment = AlignIO.read(aa_aln_file, 'fasta')
                summary = AlignInfo.SummaryInfo(aa_alignment)
                if threshold is None:
                    threshold = getattr(self, 'consensus_threshold', 0.5)
                codon_consensus = summary.dumb_consensus(threshold=threshold, ambiguous='X')
                global_p_consensus_str = str(codon_consensus)

                # --------------------------------------------Main_codon Consensus

                self.widget_global_codon = QWidget()  
                self.widget_global_codon.setObjectName('widget_global_codon')    
                layout_global_codon = QHBoxLayout()
                layout_global_codon.setContentsMargins(5,0,0,0)
                layout_global_codon.setSpacing(0)
                self.widget_global_codon.setLayout(layout_global_codon)
                self.layout_codon_l3.addWidget(self.widget_global_codon, alignment=Qt.AlignLeft)
                #layout_global_codon.setAlignment(Qt.AlignTop)

                # Sub-elements
                invisible_checkbox_codon = QCheckBox()
                invisible_checkbox_codon.setEnabled(False)
                invisible_checkbox_codon.setStyleSheet('background: transparent; border: none;')
                layout_global_codon.addWidget(invisible_checkbox_codon)
                layout_global_codon.addSpacing(5)

                label_consensus_codon = QLabel('')
                label_consensus_codon.setFixedSize(120, 20)
                layout_global_codon.addWidget(label_consensus_codon, alignment=Qt.AlignLeft)
                layout_global_codon.addSpacing(5)

                for aa in global_p_consensus_str:
                    if aa != "":
                        lbl = QLabel(aa)
                        lbl.setFixedSize(45, 20)
                        lbl.setAlignment(Qt.AlignCenter)
                        lbl.setStyleSheet('color:gray;')
                        layout_global_codon.addWidget(lbl)

        # --------------------------------------------Output_files
                global_consensus_codon_file = os.path.join(self.session_folder, f'{self.uid}_globalconsensus_codon.fasta')
                with open(global_consensus_codon_file, 'w') as c:
                    c.write(f'>global_consensus_codon\n{global_p_consensus_str}')

            except Exception as e:
                print(f"[Warning] Failed to compute AA consensus from {aa_aln_file}: {e}")

# --------------------------------------------Others
        return global_p_consensus_str



#_______________________________________________________________________________________________13 Color Code Aligned Sequences
#_______________________________________________________________________________________________

    def color_code_seq(self, seq_map=None, mode=None, group_idx=None):
# --------------------------------------------Define Inside Fxns
        def get_col_similarity(seqs):
            transposed = list(zip(*seqs))
            similarity = []
            for col in transposed:
                most_common = max(set(col), key=col.count)                              # get the most common base in a column
                similarity_score = col.count(most_common) / len(col)                    # most common base / total base in a column
                similarity.append(similarity_score)                                     # append the final value into similarity
            return similarity        

        def similarity_to_color(score):
            return mcolors.to_hex(mcolors.LinearSegmentedColormap.from_list('custom', ['white', self.similarity_color])(score))

# --------------------------------------------Main
        # ---1 Alignall
        if seq_map and mode == "all":
            all_seqs = []
            for group in self.groups:
                for entry in group['widget_seq']:
                    all_seqs.append(entry['seq'])                                       # Get all sequences

# --------------------------------------------Connect: inner DEF 1
            similarity = get_col_similarity(all_seqs)                                   # Calculate similarity
            for idx_col, score in enumerate(similarity):
                for (group_idx, entry_idx) in seq_map:
                    entry = self.groups[group_idx]['widget_seq'][entry_idx]             # grab specific element in widget_seq, for each group
                    if idx_col < len(entry['seq_letters']):
                        lbl = entry['seq_letters'][idx_col]
                        color = similarity_to_color(score)                             # grab letters from dict by using column index
# --------------------------------------------Connect: inner DEF 2
                        lbl.setProperty('bg_color', color)         # make a dynamic property global - can access later by object.property('bg_color'). setObjectName is static.
                        lbl.setStyleSheet(f'background-color: {color};')

        # ---2 Align
        elif mode == "group":
            selected_group = [self.groups[group_idx]] if group_idx is not None else self.groups
            for group in selected_group:
                group_seqs = [entry['seq'] for entry in group['widget_seq']]
                group_similarity = get_col_similarity(group_seqs)
                for idx_col, score in enumerate(group_similarity):
                    for entry in group['widget_seq']:
                        if idx_col < len(entry['seq_letters']):
                            lbl = entry['seq_letters'][idx_col]
                            color = similarity_to_color(score)
                            lbl.setProperty('bg_color', color)
                            lbl.setStyleSheet(f'background-color: {color};')


#_______________________________________________________________________________________________14 Custom Display % Conservation
#_______________________________________________________________________________________________

    def custom_display_perc_cons(self, pos1, pos2):
# --------------------------------------------Action
        # 1 Get Reference Consensus
        ref_consensus = None
        for group in self.groups:
            if group['checkbox_setrefgroup'].isChecked():
                ref_consensus = group.get('consensus_seq')
                break
        if not ref_consensus:
            print('No reference group selected.')

        # 2 Positions: MIN (Click 1) & MAX (Click 2)
        start = min(pos1, pos2)
        end = max(pos1, pos2)

        for group in self.groups:
        # 3 Get Target Consensus
            if group['checkbox_setrefgroup'].isChecked():
                continue
            target_consensus = group.get('consensus_seq')
            if not target_consensus:
                continue

        # 4 Remove existing Widget: Custom Display % Conservative
            layout = group['main_layout_seq']
            for i in reversed(range(layout.count())):
                widget = layout.itemAt(i).widget()
                if widget and widget.objectName() == 'custom_conservation_block':
                    layout.removeWidget(widget)
                    widget.deleteLater()

        # 5 Calculate % Similarity
            match_count = sum(
                1 for i in range(start, end + 1)
                if i < len(ref_consensus) and i < len(target_consensus) and ref_consensus[i] == target_consensus[i]
            )
            total = end - start + 1
            percent = (match_count / total * 100) if total else 0
            mid = (start + end) // 2                                                # // = division without decimals, / = normal division

# --------------------------------------------Main
            widget_result_main = QWidget()
            widget_result_main.setObjectName("custom_conservation_block")
            layout_result_main = QHBoxLayout()
            widget_result_main.setLayout(layout_result_main)
            layout_result_main.setContentsMargins(0, 0, 0, 0)
            layout_result_main.setSpacing(0)

# --------------------------------------------Sub-elements
            # ---1 2nd Layer
            widget_result = QWidget()
            layout_result = QHBoxLayout()
            layout_result.setContentsMargins(0, 0, 0, 0)
            layout_result.setSpacing(0)
            layout_result_main.addWidget(widget_result, alignment=Qt.AlignLeft)

            # ---2 Labels (Spacing)
            invisible_checkbox = QCheckBox()
            invisible_checkbox.setEnabled(False)
            invisible_checkbox.setStyleSheet('background: transparent; border: none;')
            layout_result.addWidget(invisible_checkbox)

            invisible_label = QLabel(f"(pos: {str(start)}-{str(end)}) {percent:.1f}%")
            invisible_label.setStyleSheet("font-size: 12px;")
            layout_result.addWidget(invisible_label)
            widget_result.setLayout(layout_result)
            layout.addWidget(widget_result_main)


#_______________________________________________________________________________________________15 PSIPRED: BUILD SECONDARY STRUCTURE
#________________________________________________________________________________________PSIPRED

    def build_secondary_structure_offline(self, psipred_dir):

        print('#  Offline: Running PSIPRED V4 with database: BLAST+')

# . . .  CLEAR GUI . . .
        # REMOVE EXISTING WIDGET SECONDARY STRUCTURE
        if hasattr(self, 'widget_horizontal') and self.widget_horizontal is not None:
            self.layout_codon_l4_2ndarystructure.removeWidget(self.widget_horizontal)     # remove from layout (optional)
            self.widget_horizontal.setParent(None)                                          # remove from parent GUI hierarchy
            self.widget_horizontal.deleteLater()                                            # schedule for safe deletion by Qt event loop
            self.widget_horizontal = None                                                   # no widget global

# --------------------------------------------Action
        # ---1 Set files

        alignedaa_file = self.aligned_aa_output_file
        base = os.path.splitext(os.path.basename(alignedaa_file))[0]                        # extract FASTA Filename
        os.makedirs(self.session_folder, exist_ok=True)                                         
        fasta_file = f"{self.session_folder}/{self.uid}_refaa.fasta"
        horiz_file = f"{self.session_folder}/{self.uid}_refaa.horiz"


# -------------------------- for aligned reference sequence ---------------------------

        refcodon = ""                               # ""=empty string, None=no value, nothing, missing
        refaa = ""
        # 1 Get ref codon seq
        for group in self.groups:
            if group['checkbox_setrefgroup'].isChecked():
                refcodon = group['widget_seq'][0]['seq']

        # 2 Translate to aa seq
        codons = [refcodon[i:i+3] for i in range(0, len(refcodon), 3)]
        for codon in codons:
            if len(codon) < 3 or '-' in codon or 'N' in codon.upper():
                refaa += 'X'  # Ambiguous or incomplete
            else:
                try:
                    aa = Seq(codon).translate()
                    refaa += str(aa)
                except:
                    refaa += 'X'

        # 3 Write to fasta_file
        with open(fasta_file, 'w') as f:
            f.write(f'>refaa\n{refaa}\n')

# -------------------------- for aligned reference sequence ---------------------------
        # ---2 Run PSIPRED
        # with PSI-BLAST
        if self.widget_psipred_blast.isChecked():
            try:
                runpsipred = os.path.join(psipred_dir, 'BLAST+', 'runpsipredplus')
                blastdb_path = os.path.join(psipred_dir, "BLAST+", "blastdb")
                env = os.environ.copy()
                env["BLASTDB"] = blastdb_path

                if shutil.which('tcsh') is None:
                    self.show_tcsh_warning()
                    return
                else:
                    subprocess.run([runpsipred, fasta_file], check=True, env=env, cwd=self.session_folder)
            except subprocess.CalledProcessError as e:
                QMessageBox.critical(self, 'PSIPRED error', f'PSIPRED with BLAST+ failed: {e}')
                raise e

        elif self.widget_psipred_single.isChecked():
            try:
                if shutil.which('tcsh') is None:
                    self.show_tcsh_warning()
                    return
                else:
                    runpsipred_single = os.path.join(psipred_dir, "runpsipred_single")
                    subprocess.run([runpsipred_single, fasta_file], check=True, cwd=self.session_folder)
            except subprocess.CalledProcessError as e:
                QMessageBox.critical(self, 'PSIPRED Error', f'PSIPRED Single failed: {e}')
                raise e    
         

        # ---3 Extract data from output: horiz file
        if not os.path.exists(horiz_file):
            raise FileNotFoundError(f"Expected output file {horiz_file} not found!")
        with open(horiz_file, "r") as f:
            self.prediction_text = f.read()
        return self.prediction_text




#_______________________________________________________________________________________________15 PSIPRED: BUILD SECONDARY STRUCTURE
#________________________________________________________________________________________PSIPRED

    def build_secondary_structure_online(self, fasta_file):

        print('#  Online: Running PSIPRED V4 online')

        with open(fasta_file, 'r') as fasta:
            lines = fasta.readlines()
        filtered_lines = [line for line in lines if not line.strip().startswith('>')]
        with open(fasta_file, 'w') as fasta:
            fasta.writelines(filtered_lines)

        # 1 Request from server
        user_email = self.qlineedit_email.text()
        if not user_email or '@' not in user_email:
            QMessageBox.warning(self, "Missing Email", "Please enter a valid email address to run PSIPRED online.")
            return

        fasta_file_name = os.path.basename(fasta_file)
        url = 'https://bioinf.cs.ucl.ac.uk/psipred/api/submission.json'

        with open(fasta_file, 'rb') as fasta:
            payload = {'input_data': (fasta_file_name, fasta)}      # (name for server to recognize as, open(inputfile, ..))
            data = {'job': 'psipred', 'submission_name': fasta_file_name, 'email': user_email}
#            print('PSIPRED: Sending Request')
            r = requests.post(url, data=data, files=payload)

            # Error trap
#            print("POST Code:", r.status_code)
#            print("POST response:", repr(r.text))
            if r.status_code != 200 and r.status_code != 201:
                raise Exception("PSIPRED submission failed: " + r.text)

            response_data = json.loads(r.text)

        # 2 Get from server
        while True:
            result_url = "https://bioinf.cs.ucl.ac.uk/psipred/api/submission/" + response_data['UUID']
            r = requests.get(result_url, headers={"Accept":"application/json"})

            if not r.text.strip():
                print("Empty response received from server.")
                print("Status code:", r.status_code)
                raise Exception("Server returned empty response while polling for results.")

            result_data = json.loads(r.text)
            print(result_data)
            if "Complete" in result_data["state"]:
                print(r.text)
                break
            else:
                time.sleep(30)

        # 3 Download horiz file from result
        horiz_path = None
        for result in result_data['submissions'][0].get("results", []):
            path = result.get("data_path", "")
            if path.endswith(".horiz"):
                horiz_path = path
                break
        if not horiz_path:
            raise ValueError("PSIPRED online: No .horiz file found in result data.")
            
        horiz_url = "https://bioinf.cs.ucl.ac.uk/psipred/api" + horiz_path
        r = requests.get(horiz_url)
        self.prediction_text = r.text

#---------------------------------------------Output files
        horiz_file = os.path.join(self.session_folder, f'{self.uid}_refseq.horiz')
        with open(horiz_file, 'w') as h:
            h.write(self.prediction_text)

        return self.prediction_text


#_______________________________________________________________________________________________16 PSIPRED: DISPLAY ON GUI
#________________________________________________________________________________________PSIPRED

    def draw_secondary_structure_to_gui(self, prediction_text, region_conservation=None):

        print('#  Drawing secondary structure from PSIPRED on GUI')

# --------------------------------------------Main
        self.widget_horizontal = QWidget()
        layout_horizontal = QHBoxLayout()
        layout_horizontal.setContentsMargins(0,0,0,0)
        layout_horizontal.addSpacing(0)
        self.widget_horizontal.setLayout(layout_horizontal)
        self.layout_codon_l4_2ndarystructure.addWidget(self.widget_horizontal, alignment=Qt.AlignLeft)

# --------------------------------------------Sub-elements
        # ---1 Labels (Spacing)
        invisible_checkbox = QCheckBox()
        invisible_checkbox.setEnabled(False)
        invisible_checkbox.setStyleSheet('background: transparent; border: none;')
        layout_horizontal.addWidget(invisible_checkbox)
   
        invisible_label = QLabel('')
        invisible_label.setFixedSize(117,20)
        layout_horizontal.addWidget(invisible_label, alignment=Qt.AlignLeft)

# --------------------------------------------Action
        # ---1 Parse horiz_rile: Prediction Text
        aa, pred = '', ''                                                   # Base, Secondary structure
        for line in prediction_text.splitlines():
            if 'AA' in line:
                parts = line.strip().split()
                if len(parts) > 1:
                    aa+= parts[1]
            elif 'Pred' in line:
                parts = line.strip().split()
                if len(parts) > 1:
                    pred += parts[1]


# -------------------------- for aligned reference sequence ---------------------------
        # 1 get aligned aa - find where it is
        refseq = ""
        for group in self.groups:
            if group['checkbox_setrefgroup'].isChecked():
                refseq = group['widget_seq'][0]['seq']

        # 2 if -, check aa b4 & after
        final_pred = []
        pred_idx = 0
        for i in range(0, len(refseq), 3):
            codon = refseq[i:i+3]
            if codon == '---':
                final_pred.append('C')  # fallback for gapped codon
            else:
                if pred_idx < len(pred):
                    final_pred.append(pred[pred_idx])
                    pred_idx += 1
                else:
                    final_pred.append('C')  # fallback if pred string runs out


 # -------------------------- for aligned reference sequence ---------------------------
        # ---2 Turn C, E, H into symbols: ---, Arrow, Box
        width_per_residue = 0.15 * 3                                            # Width per base: 15 pixels
        fig_width = len(final_pred) * width_per_residue                             # Figure width
        fig, ax = plt.subplots(figsize=(fig_width, 0.4), dpi=100)
        i=0
        troubleshoot = []
        while i < len(final_pred):
            ss = final_pred[i]                                                    # Secondary structure at pos i
            start = i                                                       # Start at pos i
            while i < len(final_pred) and final_pred[i] == ss:
                i += 1
            end = i

            if ss == 'H':
                rect = Rectangle((start, 0.1), end - start, 0.8, linewidth=1, edgecolor='red', facecolor='red', alpha=0.4)  # (x,y), width, height, border thickness, border color, fill color, semi-transparent      
                ax.add_patch(rect)   
                troubleshoot.append(f'H: {start}, {end}')                                                                                                       # add rect to plot
            elif ss == 'E':
                arrow = FancyArrow(start, 0.5, end - start - 0.2, 0, width=0.3, length_includes_head=True, head_width=0.5, head_length=0.3, color='blue') # x,y,x-length, y-change, ...
                ax.add_patch(arrow)
                troubleshoot.append(f'E: {start}, {end}')                
            else:
                ax.plot([start, end], [0.5, 0.5], color='gray', linewidth=1.2)
                troubleshoot.append(f'C: {start}, {end}')

        # ---3 Create Plot Figure
        ax.set_xlim(0, len(final_pred))
        ax.set_ylim(0, 1)
        ax.axis('off')
        plt.tight_layout(pad=0)
        buffer = BytesIO()                                                  # Save to buffer (FOR DISPLAY)
        fig.savefig(buffer, format='png', transparent=True, bbox_inches='tight', pad_inches=0)
        plt.close(fig)

        # ---4 Convert to QPixmap and Display on GUI
        pixmap = QPixmap()
        pixmap.loadFromData(buffer.getvalue())
        label = StructureLabel(pixmap, "\n".join(f"{reg['type']} ({reg['start']+1}-{reg['end']+1}):\t" + "\t".join(f"{g}: {v:.1f}%" for g, v in reg.get("group_scores", {}).items()) for reg in region_conservation) if region_conservation else "")
        layout_horizontal.addWidget(label, alignment=Qt.AlignLeft)



#-----------------------------------------------------Output files
        ss_conservation_file = os.path.join(self.session_folder, f'{self.uid}_ssconservation.txt')
        with open(ss_conservation_file, 'a') as ss:
            ss.write("\n".join(f"{reg['type']} ({reg['start']+1}-{reg['end']+1}):\t" + "\t".join(f"{g}: {v:.1f}%" for g, v in reg.get("group_scores", {}).items()) for reg in region_conservation) if region_conservation else "")




#_______________________________________________________________________________________________16 PSIPRED: COMPUTE REGION %CONSERVATIVE
#________________________________________________________________________________________PSIPRED

    def compute_region_conservation(self, prediction_text):

        print('#  Calculating %Conservation on secondary structure')

# --------------------------------------------Action
        # ---1 Get Reference Consensus Sequence
        ref_cons = None
        target_cons = None
        for group in self.groups:
            if group['checkbox_setrefgroup'].isChecked():
                ref_cons = group.get('consensus_seq')
                break
        if not ref_cons:
            return

# . . .  CLEAR GUI . . .
        for group in self.groups:
            if group['checkbox_setrefgroup'].isChecked():
                continue
            layout = group['main_layout_seq']
            for i in reversed(range(layout.count())):
                widget = layout.itemAt(i).widget()
                if widget and widget.objectName() == "conservation_block":
                    layout.removeWidget(widget)
                    widget.deleteLater()

        # ---3 Get Target Consensus Sequence
            target_cons = group.get('consensus_seq')
            if not target_cons:
                continue

            # ---4 Parse horiz_rile: Prediction Text
        aa, pred = '', ''
        for line in prediction_text.splitlines():
            if 'AA:' in line:
                parts = line.strip().split()
                if len(parts) > 1:
                    aa += parts[1]                                                  # extract the seq after 'AA:' (AA Seq)
            elif 'Pred:' in line:
                parts = line.strip().split()
                if len(parts) > 1:
                    pred += parts[1]  

# -------------------------- for aligned reference sequence ---------------------------
        # 1 get aligned aa - find where it is
        refseq = None
        for group in self.groups:
            if group['checkbox_setrefgroup'].isChecked():
                refseq = group['widget_seq'][0]['seq']

       # 2 if -, check aa b4 & after
        final_pred = []
        pred_idx = 0
        for i in range(0, len(refseq), 3):
            codon = refseq[i:i+3]
            if codon == '---':
                final_pred.append('C')  # fallback for gapped codon
            else:
                if pred_idx < len(pred):
                    final_pred.append(pred[pred_idx])
                    pred_idx += 1
                else:
                    final_pred.append('C')  # fallback if pred string runs out

# -------------------------- end: for aligned reference sequence ---------------------------
            # ---5 Collect data
        regions = []
        current_type = None
        start = None
        for i, ss in enumerate(final_pred):
            if ss in ['H', 'E']:
                if current_type != ss:                                              # If H/E but different from previous one
                    if current_type and start is not None:                          #
                        regions.append((current_type, start, i - 1))                # Save the previous region from start to i-1
                    current_type = ss                                               #
                    start = i
                # only add when it is H/E and different from the previous stored ss. If they are similar, keep looping (i++)
            else:                                                                   # If Coil, terminate the current region
                if current_type:                                                    #
                    regions.append((current_type, start, i - 1))                    #
                    current_type = None                                             # SS restarts
                    start = None                                                    # Start index restarts

        if current_type and start is not None:                                      # For the last base
            regions.append((current_type, start, len(final_pred) - 1))                    # e.g. ('H', 3, 10) = Helix from pos 3 to 10

# --------------------------------------------Action
        # ---7 Fill labels only at the midpoint of each region
        result = []
        for ss_type, start, end in regions:
            entry = {"type": ss_type, "start": start, "end": end, "group_scores": {}}
            for group in self.groups:
                group_name = group['lineedit_groupname'].text()
                if group['checkbox_setrefgroup'].isChecked():
                    continue
                target_cons = group.get('consensus_seq')
                if not target_cons:
                    continue

                codon_start = start * 3
                codon_end = (end  + 1) * 3 - 1
                match_count = sum(
                    1 for i in range(codon_start, codon_end + 1)
                    if i < len(target_cons) and i < len(ref_cons) and target_cons[i] == ref_cons[i]
                )
                total = codon_end - codon_start + 1
                percent = (match_count / total) * 100 if total else 0
                entry['group_scores'][group_name] = percent
            result.append(entry)

        return result



#_______________________________________________________________________________________________18 PSIPRED: For %Conservation per SS
#________________________________________________________________________________________PSIPRED

class StructureLabel(QLabel):
    def __init__(self, pixmap, tooltip_text, parent=None):
        super().__init__(parent)
        self.setPixmap(pixmap)
        self.tooltip_text = tooltip_text

    def enterEvent(self, event):
        QToolTip.showText(event.globalPosition().toPoint(), self.tooltip_text, self)

    def leaveEvent(self, event):
        QToolTip.hideText()

