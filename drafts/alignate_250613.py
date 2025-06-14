import sys, re, subprocess, tempfile, os, shutil, json, zipfile, tempfile, platform, requests, time, copy
from PySide6.QtGui import QIcon, QPixmap, QFont, QPainter, QPen, QMouseEvent, QAction, QKeySequence, QShortcut
from PySide6.QtWidgets import QGroupBox, QRadioButton, QStackedWidget, QProgressBar, QFileDialog, QMessageBox, QDialog, QTextEdit, QDialogButtonBox, QLayout, QScrollArea, QSizePolicy, QApplication, QMainWindow, QWidget, QCheckBox, QLabel, QLineEdit, QPushButton, QHBoxLayout, QVBoxLayout, QSlider
from PySide6.QtCore import QSize, Qt, QPoint, QTimer
from Bio import SeqIO, AlignIO
from Bio.Align import MultipleSeqAlignment, AlignInfo
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrow
import matplotlib.colors as mcolors
from io import BytesIO
from collections import Counter

from support_files.drawingcanvas import DrawingCanvas
from support_files.ruler import ruler, ClickableLabel
from support_files.codon import codon

app = QApplication(sys.argv)
app.setStyle('Fusion')
icon_logo = 'images/logo_ninja1'




################################################################################################
#FORMAT
#_______________________________________________________________________________________________
#_______________________________________________________________________________________________
#_______________________________________________________________________________________________1 CLASS: Main Window
#_______________________________________________________________________________________________
#_______________________________________________________________________________________________

#_______________________________________________________________________________________________1 DEF: ALIGNMENT
#_______________________________________________________________________________________________
# . . .  CLEAR GUI . . .
# . . . CLEAR DICT . . .
# . . .   BEGINS   . . .
# --------------------------------------------Main
# --------------------------------------------Sub-elements
# --------------------------------------------Connect
# --------------------------------------------Others
# ***** to amend *****
#
################################################################################################
#SUMMARY
# groups < widget group (header, seq, letters, widget row), group name, layout seq, widget seq, checkbox ref group
# group = {
# 'widget_group': # WIDGET GROUP: {'seq_header', 'seq', 'seq_letters', 'widget_row' (the whole row of each seq)}
# 'lineedit_groupname': # GROUP NAME
# 'layout_seq':  # LAYOUT SEQUENCE (The whole layout, not just 1 row)
# 'widget_seq': # WIDGET SEQUENCE (The whole widget, not just 1 row)
# 'checkbox1_setrefgroup': # CHECKBOX REFERENCE GROUP }
#
#
# 1 Delete temp files: tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix='.fasta'
#
################################################################################################
#DEFS
#    def show_tcsh_warning(self):
#    def view_show_all(self):
#    def view_hide_toggles(self):
#    def adjust_consensusmode(self):
#    def apply_new_consensus_threshold(self):
#    def handle_slider_mode_toggle(self, checked):       
#    def slider_threshold(self):
#    def save_project(self):
#    def load_project(self):
#    def button1_addgroup_clicked(self):
#    def handle_reference_group_toggle(self, selected_checkbox):
#    def button2_removegroup_clicked(self, widget):
#    def button3_addseq_clicked(self, layout):
#    def seqtext_button1text_clicked(self, layout):
#    def widget_seq_text_inputtext_dialogbox_button1_clicked(self, layout):
#    def seqtext_button2file_clicked(self, layout):
#    def button4_removeseq_clicked(self, layout):
#    def button5_align_clicked(self, layout):
#    def toggle_email_field(self):
#    def button2_alignall_clicked(self):
#    def get_mafft_path(self):
#    def get_clustalo_path(self):
#    def run_alignment(self, sequences, group=None, layout=None, button_aln=None, output_file=None, return_only=False, seq_map=None):
#    def add_sequences_toGUI(self, group, layout, seq_name='', seq=''):  
#    def get_consensus_aln(self, group, seq_map=None, threshold=None):
#    def get_global_consensus(self, threshold=None):
#    def color_code_seq(self, seq_map=None, mode=None, group_idx=None):
#    def custom_display_perc_cons(self, pos1, pos2):
#    def build_secondary_structure_offline(self, fasta_file, psipred_dir, base_path):
#    def build_secondary_structure_online(self, fasta_file):
#    def draw_secondary_structure_to_gui(self, prediction_text):
#    def compute_region_conservation(self, prediction_text):
#
################################################################################################




#_______________________________________________________________________________________________
#_______________________________________________________________________________________________
#_______________________________________________________________________________________________1 CLASS: main
#_______________________________________________________________________________________________
#_______________________________________________________________________________________________




class main(QMainWindow):
    def __init__(self):
        super().__init__()

# --------------------------------------------Main
        self.setWindowTitle('Alignate')
        self.setWindowIcon(QIcon(icon_logo))
        self.resize(1000, 500)                         # resize = initial size of widget (pixels)

# --------------------------------------------Sub-elements
        # ---Menu
        menu = self.menuBar()
        menu1 = menu.addMenu('File')
        menu1_load = menu1.addAction('Load Project')
        menu1_save = menu1.addAction('Save Project')
        menu1_saveas = menu1.addMenu('Save as')
        menu1_saveas_aln = menu1_saveas.addAction('.aln')           # ***** to amend 1 ***** connect           
        menu1_saveas_txt = menu1_saveas.addAction('.txt')           # ***** to amend 2 ***** connect
        menu1_saveas_png = menu1_saveas.addAction('.png')           # ***** to amend 3 ***** connect
        menu1_saveas_jpg = menu1_saveas.addAction('.jpg')           # ***** to amend 4 ***** connect
        menu1_saveas_pdf = menu1_saveas.addAction('.pdf')           # ***** to amend 5 ***** connect
        menu2 = menu.addMenu('View')
        menu2_all = menu2.addAction('All')
        menu2_hide = menu2.addAction('Hide toggles')
        menu2_consensus = menu2.addAction('Consensus mode (DEFAULT: 1)')
        menu3 = menu.addMenu('Help')

        # ---Toolbar
        self.stack = QStackedWidget()
        self.setCentralWidget(self.stack)
        self.window_about = about()
        self.window_protein = protein()
        base_path = os.path.dirname(os.path.abspath(__file__))
        self.window_codon = codon(base_path=base_path)                   # ***** to amend  ***** 1
        self.stack.addWidget(self.window_about)
        self.stack.addWidget(self.window_protein)
        self.stack.addWidget(self.window_codon)
        self.active_window = self.window_protein                        # set default: window_protein

        toolbar_spacer = QLabel()
        toolbar_spacer.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Preferred)

        widget_search_seq = QLineEdit()
        widget_search_seq.setFixedWidth(150)
        widget_search_seq.setPlaceholderText('Search sequence')

        toolbar = self.addToolBar('Main Toolbar')

# --------------------------------------------Connect   (triggered)
        # ---Menu
        menu1_load.triggered.connect(lambda: self.active_window.load_project())
        menu1_save.triggered.connect(lambda: self.active_window.save_project())
        menu2_all.triggered.connect(lambda: self.active_window.view_show_all())
        menu2_hide.triggered.connect(lambda: self.active_window.view_hide_toggles())
        menu2_consensus.triggered.connect(lambda: self.active_window.apply_new_consensus_threshold())

        self.action_about = QAction('About', self)
        self.action_about.setCheckable(True)
        self.action_protein = QAction('Protein', self)
        self.action_protein.setCheckable(True)
        self.action_codon = QAction('Codon', self)
        self.action_codon.setCheckable(True)

        toolbar.addAction(self.action_about)
        toolbar.addSeparator()
        toolbar.addAction(self.action_protein)
        toolbar.addSeparator()
        toolbar.addAction(self.action_codon)
        toolbar.addWidget(toolbar_spacer)
        toolbar.addWidget(widget_search_seq)

        self.action_about.triggered.connect(lambda: self.switch_page(self.window_about, self.action_about))
        self.action_protein.triggered.connect(lambda: self.switch_page(self.window_protein, self.action_protein))
        self.action_codon.triggered.connect(lambda: self.switch_page(self.window_codon, self.action_codon))


        def active_search_seq():                                                            # dynamically switch between protein & codon
            current_widget = self.stack.currentWidget()
            if hasattr(current_widget, 'search_sequences'):
                current_widget.search_sequences(widget_search_seq.text())
        widget_search_seq.textChanged.connect(active_search_seq)

        # Re-create Folder: tmp_files is absent
        tmp_folder = os.path.join(os.path.dirname(__file__), 'tmp_files')
        os.makedirs(tmp_folder, exist_ok=True)

    def switch_page(self, widget, active_action):
        self.stack.setCurrentWidget(widget)
        self.active_window = widget
        for action in [self.action_about, self.action_protein, self.action_codon]:
            action.setChecked(action == active_action)


    def closeEvent(self, event):                                                            # delete files in tmp_files when software is closed
        tmp_folder = os.path.join(os.path.dirname(__file__), 'tmp_files')
        if os.path.exists(tmp_folder):
            try:
                for filename in os.listdir(tmp_folder):
                    file_path = os.path.join(tmp_folder, filename)
                    if os.path.isfile(file_path):
                        os.remove(file_path)
                print(f'Folder cleared: {tmp_folder}')
            except Exception as e:
                print(f'Error when clearing {tmp_folder}: {e}')
        event.accept()




#_______________________________________________________________________________________________
#_______________________________________________________________________________________________
#_______________________________________________________________________________________________2 CLASS: Toolbar - about
#_______________________________________________________________________________________________
#_______________________________________________________________________________________________




class about(QWidget):
    def __init__(self):
        super().__init__()

# --------------------------------------------Main
        # ---Layer 1
        widget_about_l1 = QWidget()
        layout_about_l1 = QVBoxLayout()
        widget_about_l1.setLayout(layout_about_l1)
        layout_about_l1.setSizeConstraint(QLayout.SetFixedSize)

        # ---Layer 2: Widget for Texts only
        widget_about_l2 = QWidget()
        layout_about_l2 = QVBoxLayout()
        layout_about_l2.setSpacing(0)
        layout_about_l2.setSizeConstraint(QLayout.SetFixedSize)
        widget_about_l2.setLayout(layout_about_l2)  

# --------------------------------------------Sub-elements
        # ---Logo
        image_label1 = QLabel()
        image_label1.setPixmap(QPixmap(icon_logo).scaled(100, 100, Qt.KeepAspectRatio, Qt.SmoothTransformation))
        image_label1.setFixedSize(100,100)

        # ---Texts
        label_text1 = QLabel("Version: 1.0.0")
        label_text2 = QLabel("Source: githublink")
        label_text3 = QLabel("Developed by Adriana as part of her Masters Dissertation")
        label_text4 = QLabel("The University of Edinburgh")

        for label in [label_text1, label_text2, label_text3, label_text4]:
            label.setAlignment(Qt.AlignLeft | Qt.AlignTop)
            label.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Maximum)             # preferred = based on available space

        # ---Add widgets to layout
        layout_about_l1.addWidget(image_label1)
        layout_about_l1.addWidget(widget_about_l2)
        layout_about_l2.addWidget(label_text1)
        layout_about_l2.addWidget(label_text2)
        layout_about_l2.addWidget(label_text3)
        layout_about_l2.addWidget(label_text4)
        self.setLayout(layout_about_l1)




#_______________________________________________________________________________________________
#_______________________________________________________________________________________________
#_______________________________________________________________________________________________2 CLASS: Toolbar - protein
#_______________________________________________________________________________________________
#_______________________________________________________________________________________________




class protein(QWidget):
    def __init__(self):
        super().__init__()

# --------------------------------------------Others
        # Initiation
        self.cancelled = False                                          # PROGRESS (DIALOG BOX)
        self.seq_rows = []                                              # SEQUENCE ROWS BY GROUP
        self.groups = []                                                # SEQUENCE DETAILS DICT
        self.widget_toggles = []                                        # FOR MENU2_HIDE TOGGLES
        self.is_alignall = False
        self.is_searchseq = False
        self.widget_global = None

# --------------------------------------------Main
        # ---Layer 1
        # start window self
        layout_protein_l1 = QVBoxLayout()
        self.setLayout(layout_protein_l1)
        layout_protein_l1.setContentsMargins(0,0,0,0)

        # ---Layer 2
        widget_protein_l2 = QScrollArea()
        widget_protein_l2.setWidgetResizable(True)
        widget_protein_l2.setHorizontalScrollBarPolicy(Qt.ScrollBarAsNeeded)
        widget_protein_l2.setVerticalScrollBarPolicy(Qt.ScrollBarAsNeeded)
         # ---Add widgets to layout
        layout_protein_l1.addWidget(widget_protein_l2)

        # ---Layer 3
        widget_protein_l3 = QWidget()
        widget_protein_l3.setSizePolicy(QSizePolicy.MinimumExpanding, QSizePolicy.MinimumExpanding)
        self.layout_protein_l3 = QVBoxLayout()
        self.layout_protein_l3.setSpacing(0)
        self.layout_protein_l3.setContentsMargins(5,5,5,5)
        widget_protein_l3.setLayout(self.layout_protein_l3)
         # ---Add widgets to parent widget
        widget_protein_l2.setWidget(widget_protein_l3)

# --------------------------------------------Sub-elements
        # ---Ruler (In Layer 3)
        self.ruler = ruler(parent_context=self)                                        # CONNECT TO FILE 1: ruler.py
        self.layout_protein_l3.addWidget(self.ruler)

        # ---Layer 4
        self.widget_protein_l4 = QWidget()
        self.layout_protein_l4 = QVBoxLayout()
        self.layout_protein_l4.setSpacing(0)
        self.layout_protein_l4.setContentsMargins(0, 0, 0, 0)
        self.widget_protein_l4.setLayout(self.layout_protein_l4)
        self.layout_protein_l3.addWidget(self.widget_protein_l4, alignment=Qt.AlignTop)

        # ---Drawing Canvas (In Layer 4)
        self.canvas = DrawingCanvas()                                           # CONNECT TO FILE 2: drawingcanvas.py
        self.canvas.setFixedHeight(40)
        self.canvas.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        self.layout_protein_l4.addWidget(self.canvas)

        # ---Layer 4 Elements
        # Line 1
        # # 1 Button 1
        self.button1_addgroup = QPushButton('+Group')
        self.button1_addgroup.setFixedSize(QSize(60,28))
        self.button1_addgroup.setStyleSheet(
            """QPushButton {background-color: #00008B; color: white; font-weight: bold;}
            QPushButton:hover {background-color: #6495ED; color: white;}"""
        )
        # # 2 Button 2
        self.button2_alignall = QPushButton('Align all')
        self.button2_alignall.setFixedSize(QSize(60,28))
        self.button2_alignall.setStyleSheet(
            """QPushButton {background-color: #00008B; color: white; font-weight: bold;}
            QPushButton:hover {background-color: #6495ED; color: white;}"""     
        )
        # # 3 Slider
        self.slidercon = QSlider(Qt.Horizontal)
        self.slidercon.setValue(100)                                                            # Initial value: 10
        self.slidercon.setMinimum(10)
        self.slidercon.setMaximum(100)                                                          # Ticks: 1 - 10
        self.slidercon.setSingleStep(10)                                                        # Move every 10 units
        self.slidercon.setTickPosition(QSlider.TicksBelow)                                      # Tick below the slider
        self.slidercon.setFixedSize(120,20)
        self.slidercon.valueChanged.connect(self.slider_threshold)
        self.checkboxslider = QCheckBox()
        self.checkboxslider.setChecked(False)
        self.checkboxslider.toggled.connect(self.handle_slider_mode_toggle)
        # # Add widgets to parent widget
        self.widget_protein_buttons = QWidget()
        self.layout_protein_buttons = QHBoxLayout()
        self.layout_protein_buttons.setContentsMargins(0,2,0,2)
        self.widget_protein_buttons.setLayout(self.layout_protein_buttons)
        self.layout_protein_buttons.addWidget(self.button1_addgroup, alignment=Qt.AlignLeft)
        self.layout_protein_buttons.addWidget(self.button2_alignall, alignment=Qt.AlignLeft)
        self.layout_protein_buttons.addWidget(self.slidercon, alignment=Qt.AlignLeft)
        self.layout_protein_buttons.addWidget(self.checkboxslider, alignment=Qt.AlignLeft)
        self.layout_protein_buttons.addStretch(1)
        self.layout_protein_l4.addWidget(self.widget_protein_buttons, alignment=Qt.AlignTop)

        # Line 2
        self.widget_protein_l4_2ndarystructure = QWidget()
        self.layout_protein_l4_2ndarystructure = QHBoxLayout()
        self.layout_protein_l4_2ndarystructure.setContentsMargins(5,0,0,0)
        self.layout_protein_l4_2ndarystructure.setSpacing(0)
        self.widget_protein_l4_2ndarystructure.setLayout(self.layout_protein_l4_2ndarystructure)
        self.widget_protein_l4_2ndarystructure.setObjectName("line_0")
        self.widget_protein_l4_2ndarystructure.setStyleSheet(
            """ #line_0 {
                padding: 2px;
            } """
        )    
        self.layout_protein_l4.addWidget(self.widget_protein_l4_2ndarystructure)

# --------------------------------------------Connect
        self.button1_addgroup.clicked.connect(self.button1_addgroup_clicked)
        self.button2_alignall.clicked.connect(self.button2_alignall_clicked)
        self.widget_toggles.append(self.widget_protein_buttons)

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


#_______________________________________________________________________________________________2 DEF: Search Sequence
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


#_______________________________________________________________________________________________3-1 DEF: MENU2 VIEW
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
        dialog.setWindowTitle('Set Consensus Threshold')
        layout = QVBoxLayout()
        dialog.setLayout(layout)

# --------------------------------------------Sub-elements
        # ---1 label
        label = QLabel('Threshold (0.0 - 1.0) - for all consensus:')
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

#_______________________________________________________________________________________________3-2 DEF: MENU2 VIEW
#_______________________________________________________________________________________________

    def apply_new_consensus_threshold(self):
# . . .  CLEAR GUI . . . # . . . CLEAR DICT . . .
        for group in self.groups:
            layout = group['layout_seq']
            widget_to_remove = []
            target_names = ['consensus_row', 'conservation_block', 'custom_conservation_block', 'lbl_second']
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
            self.get_global_consensus(threshold=threshold)

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
            layout = group['layout_seq']                                                   
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


#_______________________________________________________________________________________________4-1 DEF: Line 1 - Slider + Slider Checkbox
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


#_______________________________________________________________________________________________4-2 DEF: Slider
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


#_______________________________________________________________________________________________5-1 DEF: Menu - Save Project
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


#_______________________________________________________________________________________________5-2 DEF: Menu - Load Project
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
            self.layout_protein_l3.removeWidget(self.widget_global)             # only remove from layout
            self.widget_global.setParent(None)                                  # detach from parent widget
            self.widget_global.deleteLater()                                    # schedule for deletion to free memory safely
            self.widget_global = None                                           # clear subject

# --------------------------------------------Main
        # MAIN WIDGET: GROUP
        self.widget_protein_l4_group_l1 = QWidget()
        self.widget_protein_l4_group_l1.setObjectName('group')         
        self.widget_protein_l4_group_l1.setStyleSheet(
            """ #group {
                border: 1px solid #A9A9A9;
                border-radius: 6px;
                padding: 2px;
            }"""
        )
        self.layout_protein_l4_group_l1 = QVBoxLayout()
        self.layout_protein_l4_group_l1.setContentsMargins(0,2,0,2)                                               
        self.widget_protein_l4_group_l1.setLayout(self.layout_protein_l4_group_l1)
        self.layout_protein_l4.addWidget(self.widget_protein_l4_group_l1)

# --------------------------------------------Sub-elements
        # ---1 LINE 1 (BUTTONS)
        # MAIN WIDGET: LINE 1
        widget_protein_l4_group_l1_line1 = QWidget()
        layout_protein_l4_group_l1_line1 = QHBoxLayout()
        widget_protein_l4_group_l1_line1.setLayout(layout_protein_l4_group_l1_line1)
        self.layout_protein_l4_group_l1.addWidget(widget_protein_l4_group_l1_line1, alignment=Qt.AlignTop)

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
        button2_removegroup = QPushButton('-Group')
        button2_removegroup.setFixedWidth(60)
        button3_addseq = QPushButton('+')
        button3_addseq.setFixedWidth(30)
        button4_removeseq = QPushButton('-')
        button4_removeseq.setFixedWidth(30)
        button5_align = QPushButton('Align')
        button5_align.setFixedWidth(60)
        button5_align.setStyleSheet("""
            QPushButton {background-color: #00008B; color: white; font-weight: bold;}
            QPushButton:hover {background-color: #6495ED; color: white;}""")

        # 2 add elements to layout
        layout_protein_l4_group_l1_line1.addWidget(lineedit_groupname, alignment=Qt.AlignLeft)
        layout_protein_l4_group_l1_line1.addWidget(button2_removegroup, alignment=Qt.AlignLeft)
        layout_protein_l4_group_l1_line1.addWidget(button3_addseq, alignment=Qt.AlignLeft)
        layout_protein_l4_group_l1_line1.addWidget(button4_removeseq, alignment=Qt.AlignLeft)
        layout_protein_l4_group_l1_line1.addWidget(button5_align, alignment=Qt.AlignLeft)
        layout_protein_l4_group_l1_line1.addWidget(checkbox1_setrefgroup, alignment=Qt.AlignLeft)
        layout_protein_l4_group_l1_line1.addWidget(label2_setrefgroup, alignment=Qt.AlignLeft)
        layout_protein_l4_group_l1_line1.addStretch(1)                      # to left-align all elements
        
        # ---2 LINE 2 (SEQUENCES)
        # MAIN WIDGET: SEQUENCES
        widget_protein_l4_group_l1_seq = QWidget()
        layout_protein_l4_group_l1_seq = QVBoxLayout()
        widget_protein_l4_group_l1_seq.setLayout(layout_protein_l4_group_l1_seq)
        layout_protein_l4_group_l1_seq.setSpacing(0)
        layout_protein_l4_group_l1_seq.setContentsMargins(5,0,0,0)
        self.layout_protein_l4_group_l1.addWidget(widget_protein_l4_group_l1_seq)

        # INITIATE DICTIONARY
        group = {
            'widget_group': self.widget_protein_l4_group_l1,        # WIDGET GROUP
            'lineedit_groupname': lineedit_groupname,               # GROUP NAME
            'layout_seq': layout_protein_l4_group_l1_seq,           # LAYOUT SEQUENCE
            'widget_seq': [],                                       # WIDGET SEQUENCE
            'checkbox_setrefgroup': checkbox1_setrefgroup,          # CHECKBOX REFERENCE GROUP
        }

# --------------------------------------------Connect
        self.groups.append(group)
        lineedit_groupname.textChanged.connect(lambda text: group.update({'name': text}))
        self.widget_toggles.append(widget_protein_l4_group_l1_line1)
        button2_removegroup.clicked.connect(lambda _=None, w=self.widget_protein_l4_group_l1: self.button2_removegroup_clicked(w))
        button3_addseq.clicked.connect(lambda _=None, layout=layout_protein_l4_group_l1_seq: self.button3_addseq_clicked(layout))
        button4_removeseq.clicked.connect(lambda _=None, layout=self.layout_protein_l4_group_l1: self.button4_removeseq_clicked(layout))
        button5_align.clicked.connect(lambda _=None, layout=layout_protein_l4_group_l1_seq: self.button5_align_clicked(layout))
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
            self.layout_protein_l3.removeWidget(self.widget_global)     # remove from layout (optional)
            self.widget_global.setParent(None)                          # remove from parent GUI hierarchy
            self.widget_global.deleteLater()                            # schedule for safe deletion by Qt event loop
            self.widget_global = None                                   # no widget global

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
            if group['layout_seq'] == layout:
                widgets_to_remove = []
                target_names = ['consensus_row', 'conservation_block', 'custom_conservation_block', 'lbl_second']
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
            self.layout_protein_l3.removeWidget(self.widget_global)
            self.widget_global.setParent(None)
            self.widget_global.deleteLater()
            self.widget_global = None

# --------------------------------------------Main
        self.widget_protein_l4_group_l1_seq_dialoginput = QDialog(self)
        self.layout_protein_l4_group_l1_seq_dialoginput = QVBoxLayout()
        self.widget_protein_l4_group_l1_seq_dialoginput.setLayout(self.layout_protein_l4_group_l1_seq_dialoginput)

# --------------------------------------------Sub-elements
        # 1 label
        widget_protein_l4_group_l1_seq_dialoginput_label = QLabel()
        widget_protein_l4_group_l1_seq_dialoginput_label.setText('Please specify method to add sequences:')
        self.layout_protein_l4_group_l1_seq_dialoginput.addWidget(widget_protein_l4_group_l1_seq_dialoginput_label)

        # 2 buttons
        self.seqtext_button1text = QPushButton('Input Text')
        self.seqtext_button2file = QPushButton('Upload File')
        self.seqtext_button3cancel = QPushButton('Cancel')
        self.widget_protein_l4_group_l1_seq_dialoginput_seqbutton = QDialogButtonBox()
        self.widget_protein_l4_group_l1_seq_dialoginput_seqbutton.addButton(self.seqtext_button1text, QDialogButtonBox.ActionRole)
        self.widget_protein_l4_group_l1_seq_dialoginput_seqbutton.addButton(self.seqtext_button2file, QDialogButtonBox.ActionRole)
        self.widget_protein_l4_group_l1_seq_dialoginput_seqbutton.addButton(self.seqtext_button3cancel, QDialogButtonBox.RejectRole)
        self.layout_protein_l4_group_l1_seq_dialoginput.addWidget(self.widget_protein_l4_group_l1_seq_dialoginput_seqbutton)

# --------------------------------------------Connect
        self.seqtext_button3cancel.clicked.connect(self.widget_protein_l4_group_l1_seq_dialoginput.reject)
        self.seqtext_button1text.clicked.connect(lambda _=None, layout=layout: self.seqtext_button1text_clicked(layout))
        self.seqtext_button2file.clicked.connect(lambda _=None, layout=layout: self.seqtext_button2file_clicked(layout))

# --------------------------------------------Others
        # Execution
        self.widget_protein_l4_group_l1_seq_dialoginput.move(500,500)
        self.widget_protein_l4_group_l1_seq_dialoginput.exec()


#_______________________________________________________________________________________________6-2 DEF: Line 1 - Add Sequence via text 1
#_______________________________________________________________________________________________

    def seqtext_button1text_clicked(self, layout):

        self.widget_protein_l4_group_l1_seq_dialoginput.hide()

# --------------------------------------------Main
        self.widget_seq_text_inputtext_dialogbox = QDialog(self)
        layout_seq_text_inputtext_dialogbox = QVBoxLayout()
        self.widget_seq_text_inputtext_dialogbox.setLayout(layout_seq_text_inputtext_dialogbox)

# --------------------------------------------Sub-elements
        # 1 label
        widget_seq_text_inputtext_dialogbox_label = QLabel()
        widget_seq_text_inputtext_dialogbox_label.setText('Sequences:')
        layout_seq_text_inputtext_dialogbox.addWidget(widget_seq_text_inputtext_dialogbox_label)

        # 2 input text box
        self.widget_seq_text_inputtext_dialogbox_textbox = QTextEdit()
        self.widget_seq_text_inputtext_dialogbox_textbox.setPlainText('>Sequence_name \nTHISISANEXAMPLEBASESEQUENCES')
        self.widget_seq_text_inputtext_dialogbox_textbox.setMinimumSize(500,400)
        layout_seq_text_inputtext_dialogbox.addWidget(self.widget_seq_text_inputtext_dialogbox_textbox)

        # 3 buttons
        widget_seq_text_inputtext_dialogbox_inputbuttons = QDialogButtonBox()
        widget_seq_text_inputtext_dialogbox_button1ok = QPushButton('OK')
        widget_seq_text_inputtext_dialogbox_button2cancel = QPushButton('Cancel')
        widget_seq_text_inputtext_dialogbox_inputbuttons.addButton(widget_seq_text_inputtext_dialogbox_button1ok, QDialogButtonBox.ActionRole)
        widget_seq_text_inputtext_dialogbox_inputbuttons.addButton(widget_seq_text_inputtext_dialogbox_button2cancel, QDialogButtonBox.RejectRole)
        layout_seq_text_inputtext_dialogbox.addWidget(widget_seq_text_inputtext_dialogbox_inputbuttons)

# --------------------------------------------Connect
        widget_seq_text_inputtext_dialogbox_button2cancel.clicked.connect(self.widget_seq_text_inputtext_dialogbox.reject)
        widget_seq_text_inputtext_dialogbox_button1ok.clicked.connect(lambda _=None, layout=layout: self.widget_seq_text_inputtext_dialogbox_button1_clicked(layout))

# --------------------------------------------Others
        # Execution
        self.widget_seq_text_inputtext_dialogbox.move(500,500)
        self.widget_seq_text_inputtext_dialogbox.exec()


#_______________________________________________________________________________________________6-3 DEF: Line 1 - Add Sequence via text 2
#_______________________________________________________________________________________________

    def widget_seq_text_inputtext_dialogbox_button1_clicked(self, layout):
# --------------------------------------------Main
        # 1 QC1: Messagebox Error if no input
        text = self.widget_seq_text_inputtext_dialogbox_textbox.toPlainText()
        if not text:
            QMessageBox.warning(self, 'Empty input', 'Please input sequences')
            return
        
        # 2 Remove spaces in text
        lines = [line.replace(' ','').strip() for line in text.splitlines() if line.strip()]    # Remove spaces

        # 3 Initiation
        self.allinput_header = []
        self.allinput_seq = []
        i = 0
        seq_count = 0

        # 4 Extract sequence header
        while i < len(lines):
            line = lines[i]
            if not line.startswith('>'):                            # if > is present
                QMessageBox.warning(self, 'Missing value', 'Missing header (>Sequence_name).')
                return
            seqname = line[1:].strip()                              # if > is present, but no name, name it as seq_#
            if not seqname:
                seqname = f"seq_{seq_count + 1}"
                seq_count += 1
            else:
                seq_count += 1
            self.allinput_header.append(seqname)                    # add to header list

        # 5 Extract sequence
            i += 1                                                  # add 1 line after header
            if i >= len(lines) or lines[i].startswith('>'):         # if > is present, cancel
                QMessageBox.warning(self, 'Missing value', f'Missing sequences: {seqname}')
                return
            sequence = []
            while i < len(lines) and not lines[i].startswith('>'):
                seq_line = lines[i].upper()                         # capitalize
                if not re.match(r"^[A-Z\>\-]+$", seq_line):          # QC: if special characters are present, cancel
                    QMessageBox.warning(self, 'Invalid sequence', f'Special characters are found in sequence: {sequence}')
                    return
                sequence.append(seq_line)                           # collect all lines before line with >
                i += 1
            full_sequence = ''.join(sequence)                       # form a full sequence
            self.allinput_seq.append(full_sequence)                 # add to sequence list
                
# --------------------------------------------Connect
        for seq_name, seq in zip(self.allinput_header, self.allinput_seq):
            for group in self.groups:
                if group['layout_seq'] == layout:
                    self.add_sequences_toGUI(group, layout, seq_name, seq)
                    break


#_______________________________________________________________________________________________6-3 DEF: Line 1 - Add Sequence via text 2
#_______________________________________________________________________________________________

    def seqtext_button2file_clicked(self, layout):

# --------------------------------------------Others
        self.widget_protein_l4_group_l1_seq_dialoginput.hide()

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
# --------------------------------------------Connect        
                for group in self.groups:
                    if group['layout_seq'] == layout:
                        self.add_sequences_toGUI(group, layout, name, seq)
        except Exception as e:
            QMessageBox.critical(self, 'Error', f'Could not load file:\n{e}')


#_______________________________________________________________________________________________7 DEF: Line 1 - Remove Sequence
#_______________________________________________________________________________________________

    def button4_removeseq_clicked(self, layout):
        
# . . .  CLEAR GUI . . .
        for group in self.groups:
            if group['layout_seq'] == layout:
                widgets_to_remove = []
                target_names = ['consensus_row', 'conservation_block', 'custom_conservation_block', 'lbl_second']
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
            self.layout_protein_l3.removeWidget(self.widget_global)
            self.widget_global.setParent(None)
            self.widget_global.deleteLater()
            self.widget_global = None

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


#_______________________________________________________________________________________________8-1 DEF: Line 1 - Align sequences
#__________________________________________________________________________________________ALIGN

    def button5_align_clicked(self, layout):
        
# --------------------------------------------Main
        widget_dialogbox_align = QDialog(self)
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
        self.widget_psipred_offline.toggled.connect(self.toggle_email_field)
        self.widget_psipred_online.toggled.connect(self.toggle_email_field)

# --------------------------------------------Others
        widget_dialogbox_align.move(500,500)
        widget_dialogbox_align.exec()


#_______________________________________________________________________________________________8-2 DEF: Line 1 - Align sequences
#__________________________________________________________________________________________ALIGN

    def toggle_email_field(self):
        if self.widget_psipred_online.isChecked():
            self.qlineedit_email.setVisible(True)
        else:
            self.qlineedit_email.setVisible(False)


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
        self.widget_psipred_offline.toggled.connect(self.toggle_email_field)
        self.widget_psipred_online.toggled.connect(self.toggle_email_field)

# --------------------------------------------Others
        widget_dialogbox_alignall.move(500,500)
        widget_dialogbox_alignall.exec()


#_______________________________________________________________________________________________10-1 DEF: Set path for MAFFT
#__________________________________________________________________________________________ALIGN

    def get_mafft_path(self):
        base_path = os.path.dirname(os.path.abspath(__file__))
        if sys.platform.startswith('win'):                                                  # win
            return os.path.join(base_path, 'external_tools', 'mafft_win64', 'mafft.bat')
        elif sys.platform.startswith('darwin'):                                             # macOS
            return os.path.join(base_path, 'external_tools', 'mafft_mac', 'mafft')
        elif sys.platform.startswith('linux'):                                              # linux
            return os.path.join(base_path, 'external_tools', 'mafft_linux', 'mafft')
        else:
            raise RuntimeError('Unsupported OS')


#_______________________________________________________________________________________________10-2 DEF: Set path for ClustalO
#__________________________________________________________________________________________ALIGN

    def get_clustalo_path(self):
        base_path = os.path.dirname(os.path.abspath(__file__))
        if sys.platform.startswith('win'):
            return os.path.join(base_path, 'external_tools', 'clustalo_win64', 'clustalo.exe')
        elif sys.platform.startswith('darwin'):                                 
            return os.path.join(base_path, 'external_tools', 'clustalo_mac', 'clustalo')
        elif sys.platform.startswith('linux'):
            return os.path.join(base_path, 'external_tools', 'clustalo_linux', 'clustalo')
        else:
            raise RuntimeError('Unsupported OS')


#_______________________________________________________________________________________________10-3 DEF: Run Alignment > etc.
#_______________________________________________________________________________________________

    def run_alignment(self, sequences, group=None, layout=None, button_aln=None, output_file=None, return_only=False, seq_map=None):

# --------------------------------------------Others
        # Initiation
        self.cancelled = False                  # To cancel when analysis is running
        aligned_seq = []
        self.is_searchseq = False               # Retain the color after alignment in letters

        # If no group is set as a reference
        any_checked = any(group['checkbox_setrefgroup'].isChecked() for group in self.groups)
        if not any_checked and self.groups:
            self.groups[0]['checkbox_setrefgroup'].setChecked(True)

        # QC 1: No. of sequences need to be > 2
        if len(sequences) < 2:
            QMessageBox.warning(self, 'Error', 'Need at least 2 sequences for alignment.')
            return
        
        # Write FASTA file for alignment (Feed in sequences)
        temp_fasta = tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix='.fasta')
        output_file = temp_fasta.name.replace('.fasta', '.aln')
        for name, seq in sequences:
            if not name or not seq:
                QMessageBox.warning(self, 'Error', 'Missing sequence name or sequence.')
                return
            if not re.match(r'^[A-Za-z\_\-]+$', seq):
                QMessageBox.warning(self, 'Error', 'Invalid characters found in sequence: {seq}. **Special characters are not allowed.')
                return
            temp_fasta.write(f'>{name}\n{seq}\n')
        temp_fasta.close()

# --------------------------------------------Main
        self.widget_progress = QDialog(self)
        parent_pos = self.geometry().center()
        self.widget_progress.move(parent_pos.x() - self.widget_progress.width() // 2, parent_pos.y() - self.widget_progress.height() // 2)
        self.widget_progress.setModal(True)
        self.widget_progress.setWindowTitle('Processing...')
        self.layout_progress = QVBoxLayout()
        self.widget_progress.setLayout(self.layout_progress)

# --------------------------------------------Sub-elements
        # 1 Progress Bar
        progress_bar = QProgressBar()
        progress_bar.setRange(0,0)
        self.layout_progress.addWidget(progress_bar)
        # 2 Button
        self.button_cancel = QPushButton('Cancel')
        self.button_cancel.setFixedWidth(60)
        self.layout_progress.addWidget(self.button_cancel)

# --------------------------------------------Connect
        # Button
        self.button_cancel.clicked.connect(lambda: [setattr(self, 'cancelled', True), self.widget_progress.reject()])
#           ------------------------------- Connect: Widget Progress 1 -------------------------------
        self.widget_progress.show()
        QApplication.processEvents()                  # make it real time
#           ------------------------------- Connect: Widget Progress 1 -------------------------------

# --------------------------------------------Actions
# ------1 Run Alignment: MAFFT/ClustlO
        try:
            # MAFFT 
            if button_aln.text() == 'MAFFT':
                mafft_path = self.get_mafft_path()
                with open(output_file, 'w') as out:
                    subprocess.run([mafft_path, '--anysymbol', '--genafpair', '--maxiterate', '10000', temp_fasta.name], check=True, stdout=out, cwd=os.path.dirname(mafft_path))
            # ClustalO
            elif button_aln.text() == 'ClustalO':
                clustalo_path = self.get_clustalo_path()
                with open(output_file, 'w') as out:
                    subprocess.run([clustalo_path, '-i', temp_fasta.name, '-o', output_file, '--dealign', '--force'], check=True, stdout=out)
            # Neither
            else:
                QMessageBox.warning(self, 'Error', 'Error: Failed to select alignment method. Please contact Alignate.')
                return
        except Exception as e:
            QMessageBox.critical(self, f'{button_aln.text()} error', str(e))
        finally:
            os.remove(temp_fasta.name)


# ------2 Parse file and get aligned sequences
        for record in SeqIO.parse(output_file, 'fasta'):
            aligned_seq.append((record.id, str(record.seq)))                # extract name and seq for aligned sequences

#           ------------------------------- Connect: Widget Progress 2 -------------------------------
        if self.cancelled:
            self.widget_progress.reject()
            return
#           ------------------------------- Connect: Widget Progress 2 -------------------------------

        if return_only:                                                     # ***** to amend *****  Possibly can delete                                             
            return aligned_seq
        

# ------3 Other actions (split by fxn: 1 Alignall & 2 Align)


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
#           ------------------------------- Connect: Widget Progress 3 -------------------------------
                if self.cancelled:
                    self.widget_progress.reject()
                    return
#           ------------------------------- Connect: Widget Progress 3 -------------------------------        
                group = self.groups[group_idx]                              # get each group
                layout = group['layout_seq']                                # get sequence layout
                if layout is None:
                    continue
                self.add_sequences_toGUI(group, layout, aligned_name, aligned_seq)

            # -- 3 Connect - DEF: Get and Display Consensus in each group
            for group in self.groups:
#           ------------------------------- Connect: Widget Progress 4 -------------------------------
                if self.cancelled:
                    self.widget_progress.reject()
                    return
#           ------------------------------- Connect: Widget Progress 4 -------------------------------
                self.get_consensus_aln(group, seq_map, threshold=None)

            # -- 4 not DEF: Calculate %Base conservation
            # 1 get reference consensus
            ref_consensus = None
            for group in self.groups:
                if group['checkbox_setrefgroup'].isChecked():                                   
                    ref_consensus = group.get('consensus_seq')
                    break
            # 2 get group consensus and calculate %conservation
            for group in self.groups:
                layout = group['layout_seq']                                                   
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

            # -- 6 Connect - DEF Get & Display on GUI (Global Consensus)
            global_consensus = self.get_global_consensus()
            if global_consensus:
#           ------------------------------- Connect: Widget Progress 5 -------------------------------
                if self.cancelled:
                    self.widget_progress.reject()
                    return
#           ------------------------------- Connect: Widget Progress 5 -------------------------------
                with tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix='.fasta') as f:
                    f.write(f">global_consensus\n{global_consensus.replace('-', '')}\n")
                    fasta_file = f.name

                with open(fasta_file, 'r') as debug_f:
                    print("FASTA contents:\n", debug_f.read())

            # -- 7 Connect - DEF Run PSIPRED: Build Secondary Structure
                base_path = os.path.dirname(os.path.abspath(__file__))
                if self.widget_psipred_offline.isChecked():
                    psipred_dir = os.path.join(base_path, 'external_tools', 'psipred')
                    self.prediction_text = self.build_secondary_structure_offline(fasta_file, psipred_dir, base_path)
                else:
                    self.prediction_text = self.build_secondary_structure_online(fasta_file)

            # -- 8 Connect - DEF Display on GUI (PSIPRED Output)
                self.draw_secondary_structure_to_gui(self.prediction_text)

            # -- 9 Connect - DEF Compute & Display (% Conservation based on PSIPRED Output)
                self.compute_region_conservation(self.prediction_text)


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
#           ------------------------------- Connect: Widget Progress 5 -------------------------------
                    if self.cancelled:
                        self.widget_progress.reject()
                        return
#           ------------------------------- Connect: Widget Progress 5 -------------------------------
                    self.add_sequences_toGUI(group, layout, name, seq)

            # -- 3 Connect - DEF: Get and Display Consensus in each group
                self.get_consensus_aln(group, threshold=None)

            # -- 4 Connect - DEF Color Code & Display on GUI (Aligned Sequences)
                group_idx = self.groups.index(group)
                self.is_alignall = False
                self.color_code_seq(mode="group", group_idx=group_idx)

            # -- 5 Connect - DEF Run PSIPRED: Build Secondary Structure
                if 'consensus_seq' in group:
                    seq = group['consensus_seq'].replace('-', '')
                    with tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix='.fasta') as f:
                        f.write(f">consensus\n{seq}\n")
                        fasta_file = f.name
                    base_path = os.path.dirname(os.path.abspath(__file__))

                    if self.widget_psipred_offline.isChecked():
                        psipred_dir = os.path.join(base_path, 'external_tools', 'psipred')
                        self.prediction_text = self.build_secondary_structure_offline(fasta_file, psipred_dir, base_path)
                    else:
                        self.prediction_text = self.build_secondary_structure_online(fasta_file)


            # -- 6 Connect - DEF Display on GUI (PSIPRED Output)
                    self.draw_secondary_structure_to_gui(self.prediction_text)
                else:
                    QMessageBox.warning(self, 'Missing value', 'Consensus Sequences not found. Please contact Alignate.')
                    return

# --------------------------------------------Others
        self.widget_progress.close()

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
#           ------------------------------- Connect: Widget Progress 5 -------------------------------
            if self.cancelled:
                self.widget_progress.reject()
                return
#           ------------------------------- Connect: Widget Progress 5 -------------------------------
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
        current_width = self.widget_protein_l4.minimumWidth()
        if seq_pixel_length > current_width:
            self.widget_protein_l4.setMinimumWidth(seq_pixel_length)


#_______________________________________________________________________________________________12-1 Get & Display Consensus (by Group)
#_______________________________________________________________________________________________

    def get_consensus_aln(self, group, seq_map=None, threshold=None):
        
# . . .  CLEAR GUI . . .
        layout = group['layout_seq']
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
            threshold = getattr(self, 'consensus_threshold', 1.0)
        consensus = summary.dumb_consensus(threshold=threshold, ambiguous='N')
        consensus_str = str(consensus)


        consensus_str_widget = []
        for letter in consensus_str:
            lbl = QLabel(letter)
            lbl.setFixedSize(15,20)
            lbl.setAlignment(Qt.AlignCenter)
            lbl.setStyleSheet('color:gray;')
            layout_consensus.addWidget(lbl)
            consensus_str_widget.append(lbl)

# --------------------------------------------Connect
        group['consensus_seq'] = consensus_str
        group['consensus_seq_widget'] = consensus_str_widget

# --------------------------------------------Others
        return consensus_str


#_______________________________________________________________________________________________12-1 Get & Display Global Consensus
#_______________________________________________________________________________________________

    def get_global_consensus(self, threshold=None):
        
# . . .  CLEAR GUI . . .
        if hasattr(self, 'widget_global') and self.widget_global is not None:
            self.layout_protein_l3.removeWidget(self.widget_global)
            self.widget_global.setParent(None)
            self.widget_global.deleteLater()
            self.widget_global = None

# --------------------------------------------Main
        self.widget_global = QWidget()
        layout_global = QHBoxLayout()
        layout_global.setContentsMargins(5,0,0,0)
        layout_global.setSpacing(0)
        self.widget_global.setLayout(layout_global)
        self.layout_protein_l3.addWidget(self.widget_global, alignment=Qt.AlignLeft)

# --------------------------------------------Sub-elements
        # 1 Checkbox (only for spacing)
        invisible_checkbox = QCheckBox()
        invisible_checkbox.setEnabled(False)
        invisible_checkbox.setStyleSheet('background: transparent; border: none;')
        layout_global.addWidget(invisible_checkbox)
        layout_global.addSpacing(5)

        # 2 Label
        invisible_label = QLabel('Consensus')
        invisible_label.setFixedSize(120,20)
        layout_global.addWidget(invisible_label, alignment=Qt.AlignLeft)
        layout_global.addSpacing(5)

# --------------------------------------------Action
        all_records = []
        for group in self.groups:
#           ------------------------------- Connect: Widget Progress 5 -------------------------------
            if self.cancelled:
                self.widget_progress.reject()
                return
#           ------------------------------- Connect: Widget Progress 5 -------------------------------
            for entry in group['widget_seq']:
                name = entry['seq_header'].text().strip()
                seq = ''.join(label.text() for label in entry['seq_letters']).strip()
                all_records.append(SeqRecord(Seq(seq), id=name))
        if not all_records:
            return None
        alignment = MultipleSeqAlignment(all_records)
        summary = AlignInfo.SummaryInfo(alignment)
        if threshold is None:
            threshold = 1.0
        consensus = summary.dumb_consensus(threshold=threshold, ambiguous='N')

        for letter in str(consensus):
            lbl = QLabel(letter)
            lbl.setFixedSize(15,20)
            lbl.setAlignment(Qt.AlignCenter)
            lbl.setStyleSheet('color:gray;')
            layout_global.addWidget(lbl)

# --------------------------------------------Others
        return str(consensus)


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
            return mcolors.to_hex(mcolors.LinearSegmentedColormap.from_list('custom', ['white', '#5b005b'])(score))

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
            layout = group['layout_seq']
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

            invisible_label = QLabel('Pos:'+ str(start) + ' - ' + str(end))
            invisible_label.setStyleSheet("color: gray; font-size: 9px;")
            layout_result.addWidget(invisible_label)

            # ---3 
            for i in range(len(ref_consensus)):
                lbl = QLabel('')
                lbl.setFixedSize(15, 20)
                lbl.setAlignment(Qt.AlignCenter)
                if i == mid:
                    lbl.setText(f"{percent:.0f}")
                    lbl.setStyleSheet("color: gray; font-size: 8px;")
                layout_result.addWidget(lbl)
            widget_result.setLayout(layout_result)
            layout.addWidget(widget_result_main)


#_______________________________________________________________________________________________15 PSIPRED: BUILD SECONDARY STRUCTURE
#________________________________________________________________________________________PSIPRED

    def build_secondary_structure_offline(self, fasta_file, psipred_dir, base_path):
# . . .  CLEAR GUI . . .
        # REMOVE EXISTING WIDGET SECONDARY STRUCTURE
        if hasattr(self, 'widget_horizontal') and self.widget_horizontal is not None:
            self.layout_protein_l4_2ndarystructure.removeWidget(self.widget_horizontal)     # remove from layout (optional)
            self.widget_horizontal.setParent(None)                                          # remove from parent GUI hierarchy
            self.widget_horizontal.deleteLater()                                            # schedule for safe deletion by Qt event loop
            self.widget_horizontal = None                                                   # no widget global

# --------------------------------------------Action
        # ---1 Set files
        base = os.path.splitext(os.path.basename(fasta_file))[0]                            # Extract FASTA Filename
        out_folder = os.path.join(base_path, 'tmp_files')                                   # Set output folder
        os.makedirs(out_folder, exist_ok=True)                                              # 
        horiz_file = f"{out_folder}/{base}.horiz"

        # ---2 Run PSIPRED
        # with PSI-BLAST
        try:
#           ------------------------------- Connect: Widget Progress 5 -------------------------------
            if self.cancelled:
                self.widget_progress.reject()
                return
#           ------------------------------- Connect: Widget Progress 5 -------------------------------
            base_path = os.path.dirname(os.path.abspath(__file__))
            shell_tcsh = os.path.join(base_path, 'external_tools', 'cygwin', 'bin', 'tcsh.exe')
            runpsipred = os.path.join(psipred_dir, 'BLAST+', 'runpsipredplus')
            blastdb_path = os.path.join(psipred_dir, "BLAST+", "blastdb")
            env = os.environ.copy()
            env["BLASTDB"] = blastdb_path
            system = platform.system()

            if system == 'Windows':
                subprocess.run([shell_tcsh, runpsipred, fasta_file], check=True, env=env, cwd=out_folder)
                print("Run runpsipredplus + PSI-BLAST")
            elif system == 'Linux' or system == 'Darwin':
                if shutil.which('tcsh') is None:
                    self.show_tcsh_warning()
                    return
                else:
                    subprocess.run([runpsipred, fasta_file], check=True, env=env, cwd=out_folder)

        # without PSI-BLAST
        except subprocess.CalledProcessError as e:
            print('Run runpsipred_single instead...')
            try:
#           ------------------------------- Connect: Widget Progress 5 -------------------------------
                if self.cancelled:
                    self.widget_progress.reject()
                    return
#           ------------------------------- Connect: Widget Progress 5 -------------------------------
                runpsipred_single = os.path.join(psipred_dir, "runpsipred_single")
                subprocess.run([runpsipred_single, fasta_file], check=True, cwd=out_folder)
                print("Run runpsipred_single")
            except subprocess.CalledProcessError as e:
                print("Both runpsipred and runpsipred_single failed.")
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
        print('online')

        with open(fasta_file, 'r') as fasta:
            lines = fasta.readlines()
        filtered_lines = [line for line in lines if not line.strip().startswith('>')]
        with open(fasta_file, 'w') as fasta:
            fasta.writelines(filtered_lines)
        with open(fasta_file, 'r') as fasta:
            print(fasta.read())

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
            print('PSIPRED: Sending Request')
            r = requests.post(url, data=data, files=payload)

            # Error trap
            print("POST Code:", r.status_code)
            print("POST response:", repr(r.text))
            if r.status_code != 200 and r.status_code != 201:
                raise Exception("PSIPRED submission failed: " + r.text)

            response_data = json.loads(r.text)

        # 2 Get from server
        while True:
            print('before requesting r')
            result_url = "https://bioinf.cs.ucl.ac.uk/psipred/api/submission/" + response_data['UUID']
            r = requests.get(result_url, headers={"Accept":"application/json"})
            print('after requesting r')

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
        for path in result_data.get("data_path", []):
            if path.endswith(".horiz"):
                horiz_path = path
                break
        if not horiz_path:
            raise ValueError("PSIPRED online: No .horiz file found in result data.")
            
        horiz_url = "https://bioinf.cs.ucl.ac.uk/psipred/api/submissions/" + horiz_path
        r = requests.get(horiz_url)
        self.prediction_text = r.text
        return self.prediction_text

#_______________________________________________________________________________________________16 PSIPRED: DISPLAY ON GUI
#________________________________________________________________________________________PSIPRED

    def draw_secondary_structure_to_gui(self, prediction_text):
        
# --------------------------------------------Main
        self.widget_horizontal = QWidget()
        layout_horizontal = QHBoxLayout()
        layout_horizontal.setContentsMargins(0,0,0,0)
        layout_horizontal.addSpacing(0)
        self.widget_horizontal.setLayout(layout_horizontal)
        self.layout_protein_l4_2ndarystructure.addWidget(self.widget_horizontal, alignment=Qt.AlignLeft)

# --------------------------------------------Sub-elements
        # ---1 Labels (Spacing)
        invisible_checkbox = QCheckBox()
        invisible_checkbox.setEnabled(False)
        invisible_checkbox.setStyleSheet('background: transparent; border: none;')
        layout_horizontal.addWidget(invisible_checkbox)
   
        invisible_label = QLabel('')
        invisible_label.setFixedSize(118,20)
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

        # ---2 Turn C, E, H into symbols: ---, Arrow, Box
        width_per_residue = 0.15                                            # Width per base: 15 pixels
        fig_width = len(aa) * width_per_residue                             # Figure width
        fig, ax = plt.subplots(figsize=(fig_width, 0.4), dpi=100)
        i=0
        while i < len(pred):
            ss = pred[i]                                                    # Secondary structure at pos i
            start = i                                                       # Start at pos i
            while i < len(pred) and pred[i] == ss:
                i += 1
            end = i

            if ss == 'H':
                rect = Rectangle((start, 0.1), end - start, 0.8, linewidth=1, edgecolor='red', facecolor='red', alpha=0.4)  # (x,y), width, height, border thickness, border color, fill color, semi-transparent      
                ax.add_patch(rect)                                                                                          # add rect to plot
            elif ss == 'E':
                arrow = FancyArrow(start, 0.5, end - start - 0.2, 0, width=0.3, length_includes_head=True, head_width=0.5, head_length=0.3, color='blue') # x,y,x-length, y-change, ...
                ax.add_patch(arrow)
            else:
                ax.plot([start, end], [0.5, 0.5], color='gray', linewidth=1.2)

        # ---3 Create Plot Figure
        ax.set_xlim(0, len(aa))
        ax.set_ylim(0, 1)
        ax.axis('off')
        plt.tight_layout(pad=0)
        buffer = BytesIO()                                                  # Save to buffer (FOR DISPLAY)
        fig.savefig(buffer, format='png', transparent=True, bbox_inches='tight', pad_inches=0)
        plt.close(fig)

        # ---4 Convert to QPixmap and Display on GUI
        pixmap = QPixmap()
        pixmap.loadFromData(buffer.getvalue())
        label = QLabel()
        label.setPixmap(pixmap)
        layout_horizontal.addWidget(label, alignment=Qt.AlignLeft)

#_______________________________________________________________________________________________16 PSIPRED: COMPUTE REGION %CONSERVATIVE
#________________________________________________________________________________________PSIPRED

    def compute_region_conservation(self, prediction_text):
        
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
            layout = group['layout_seq']
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

            # ---5 Collect data
            regions = []
            current_type = None
            start = None
            for i, ss in enumerate(pred):
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
                regions.append((current_type, start, len(pred) - 1))                    # e.g. ('H', 3, 10) = Helix from pos 3 to 10

    # --------------------------------------------Main
            widget_result = QWidget()
            widget_result.setObjectName("conservation_block")
            layout_result = QHBoxLayout()
            layout_result.setContentsMargins(5, 0, 0, 0)
            layout_result.setSpacing(0)

            invisible_checkbox = QCheckBox()
            invisible_checkbox.setEnabled(False)
            invisible_checkbox.setStyleSheet('background: transparent; border: none;')
            layout_result.addWidget(invisible_checkbox)

            lbl_second = QLabel('%Conservation')
            lbl_second.setObjectName("lbl_second")
            lbl_second.setFixedSize(120,20)
            layout_result.addWidget(lbl_second, alignment=Qt.AlignLeft)

    # --------------------------------------------Action
            # ---6 Pre-fill all columns with empty labels
            total_cols = len(ref_cons)
            labels = [QLabel('') for _ in range(total_cols)]
            for lbl in labels:
                lbl.setFixedSize(15, 20)
                lbl.setAlignment(Qt.AlignCenter)
                layout_result.addWidget(lbl)

            # ---7 Fill labels only at the midpoint of each region
            for ss_type, start, end in regions:
                match_count = sum(
                    1 for i in range(start, end + 1)
                    if i < len(target_cons) and i < len(ref_cons) and target_cons[i] == ref_cons[i]
                )
                total = end - start + 1
                percent = (match_count / total) * 100 if total else 0
                mid = (start + end) // 2
                if mid < len(labels):
                    labels[mid].setText(f"{int(percent)}")
                    labels[mid].setStyleSheet("color: gray; font-size: 8px; padding: 0px;")
                    labels[mid].setSizePolicy(QSizePolicy.Ignored, QSizePolicy.Fixed)
                    labels[mid].setAlignment(Qt.AlignCenter)

    # --------------------------------------------Main
            widget_result.setLayout(layout_result)
            main_widget_result = QWidget()
            main_layout_result = QHBoxLayout()
            main_layout_result.setContentsMargins(0,0,0,0)
            main_layout_result.setSpacing(0)
            main_widget_result.setLayout(main_layout_result)
            main_layout_result.addWidget(widget_result, alignment=Qt.AlignLeft)
            layout.addWidget(main_widget_result)

# Execute
window = main()
window.show()
app.exec()










# . . .  CLEAR GUI . . .
# . . . CLEAR DICT . . .
# . . .   BEGINS   . . .
# --------------------------------------------Main
# --------------------------------------------Sub-elements
# --------------------------------------------Connect
# --------------------------------------------Others
# ***** to amend *****


