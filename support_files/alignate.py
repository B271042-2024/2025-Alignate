import sys, re, subprocess, tempfile, os, shutil, json, zipfile, tempfile, platform, requests, time, uuid, warnings
from PySide6.QtGui import QIcon, QPixmap, QFont, QPainter, QPen, QMouseEvent, QAction, QKeySequence, QShortcut
from PySide6.QtWidgets import QGraphicsDropShadowEffect, QToolTip, QGroupBox, QRadioButton, QStackedWidget, QProgressBar, QFileDialog, QMessageBox, QDialog, QTextEdit, QDialogButtonBox, QLayout, QScrollArea, QSizePolicy, QApplication, QMainWindow, QWidget, QCheckBox, QLabel, QLineEdit, QPushButton, QHBoxLayout, QVBoxLayout, QSlider
from PySide6.QtCore import QSize, Qt, QPoint, QTimer, Signal, QThread
from Bio import SeqIO, AlignIO, BiopythonDeprecationWarning
from Bio.Align import MultipleSeqAlignment, AlignInfo
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrow
import matplotlib.colors as mcolors
from matplotlib.figure import Figure
from io import BytesIO
from collections import Counter

#from files
from drawingcanvas import DrawingCanvas
from ruler import ruler, ClickableLabel
from codon import codon

#begin system coding
warnings.filterwarnings("ignore", category=BiopythonDeprecationWarning)
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
        #self.setStyleSheet('background-color: #F8F8FF;')

# --------------------------------------------Sub-elements
        # ---Menu
        menu = self.menuBar()
        menu1 = menu.addMenu('File')
        menu1_load = menu1.addAction('Load Project')
        menu1_save = menu1.addAction('Save Project')
        menu1_others = menu1.addMenu('Others')
        menu1_others_saveconservationfile = menu1_others.addAction('Save Conservation (Only for %Conservation slider)')
        menu2 = menu.addMenu('View')
        menu2_color = menu2.addMenu('Color')
        menu2_color_orange = menu2_color.addAction('Orange')
        menu2_color_purple = menu2_color.addAction('Purple')
        menu2_color_green = menu2_color.addAction('Olive')
        menu2_color_blue = menu2_color.addAction('Red')

        menu2_toggles = menu2.addMenu('Toggles')
        menu2_all = menu2_toggles.addAction('Show')
        menu2_hide = menu2_toggles.addAction('Hide')
        menu2_consensus = menu2.addAction('Consensus mode (DEFAULT: 0.5)')
        menu3 = menu.addAction('Help')

        # ---Toolbar
        # window: about, protein, codon
        self.stack = QStackedWidget()
        self.setCentralWidget(self.stack)
        self.window_about = about()
        self.window_about.setStyleSheet('background-color: #F8F8FF')
        self.window_protein = protein()
        self.window_protein.setStyleSheet('background-color: #F8F8FF')
        base_path = os.path.dirname(os.path.abspath(__file__))
        self.window_codon = codon(base_path=base_path)                   
        self.window_codon.setStyleSheet('background-color: #F8F8FF')
        self.stack.addWidget(self.window_about)
        self.stack.addWidget(self.window_protein)
        self.stack.addWidget(self.window_codon)
        self.active_window = self.window_protein    
        toolbar_spacer = QLabel()
        toolbar_spacer.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Preferred)
        # search bar
        widget_search_seq = QLineEdit()
        widget_search_seq.setToolTip('!Only Sequence Match. Matches are in orange.')
        widget_search_seq.setFixedWidth(120)
        widget_search_seq.setPlaceholderText('Search sequence...')
        toolbar = self.addToolBar('Main Toolbar')

        # ---Style: MainWindow, Menu
        self.setStyleSheet("""
            QMainWindow {background-color: #F8F8FF;}
            QMenuBar::item {background-color: #eeeeee; padding: 4px 10px;}
            QMenuBar::item:selected {background-color: #87CEFA;}
            QMenu {background-color: #eeeeee; padding: 4px 1px; border: 1px gray;}
            QMenu::item {background-color: #eeeeee; padding: 4px 20px;}
            QMenu::item:selected {background-color: #87CEFA;}
        """)

# --------------------------------------------Connect   (triggered)
        # ---Menu 1
        menu1_load.triggered.connect(lambda: self.active_window.load_project())
        menu1_save.triggered.connect(lambda: self.active_window.save_project())
        menu2_all.triggered.connect(lambda: self.active_window.view_show_all())
        menu2_hide.triggered.connect(lambda: self.active_window.view_hide_toggles())
        menu2_consensus.triggered.connect(lambda: self.active_window.apply_new_consensus_threshold())
        menu3.triggered.connect(self.show_widget_help)

        # ---Toolbar
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

        # ---Menu 2: Color
        menu2_color_purple.triggered.connect(lambda: self.active_window.set_similarity_color("darkmagenta"))
        menu2_color_green.triggered.connect(lambda: self.active_window.set_similarity_color("olive"))
        menu2_color_orange.triggered.connect(lambda: self.active_window.set_similarity_color("darkorange"))
        menu2_color_blue.triggered.connect(lambda: self.active_window.set_similarity_color("DarkRed"))

        # ---Menu 3: Save residues at position with the pre-set conservation value
        menu1_others_saveconservationfile.triggered.connect(lambda: self.active_window.menu1_others_saveconservationfile_clicked())

        # ---Def: Search bar
        def active_search_seq():
            current_widget = self.stack.currentWidget()
            if hasattr(current_widget, 'search_sequences'):
                current_widget.search_sequences(widget_search_seq.text())
        widget_search_seq.textChanged.connect(active_search_seq)

        # ---Re-create folder if output_folder is missing
        output_folder = os.path.join(os.path.dirname(__file__), '..', 'output_files')
        os.makedirs(output_folder, exist_ok=True)


#_______________________________________________________________________________________________1 DEF: Switch windows_about, protein, codon
#_______________________________________________________________________________________________
    def switch_page(self, widget, active_action):
        self.stack.setCurrentWidget(widget)
        self.active_window = widget
        for action in [self.action_about, self.action_protein, self.action_codon]:
            action.setChecked(action == active_action)


#_______________________________________________________________________________________________2 DEF: Help Documentation
#_______________________________________________________________________________________________
    def show_widget_help(self):

        dialog_help = QDialog(self)
        dialog_help.setWindowFlags(Qt.Dialog | Qt.WindowTitleHint | Qt.WindowCloseButtonHint)
        dialog_help.setWindowTitle('Help')
        dialog_help.resize(1000, 300)
        dialog_help.move(50,400)
        layout_help = QVBoxLayout()
        dialog_help.setLayout(layout_help)

        dialog_help_l2 = QScrollArea()
        dialog_help_l2.setWidgetResizable(True)
        dialog_help_l2.setHorizontalScrollBarPolicy(Qt.ScrollBarAsNeeded)
        dialog_help_l2.setVerticalScrollBarPolicy(Qt.ScrollBarAsNeeded)
        layout_help.addWidget(dialog_help_l2)

        dialog_help_l3 = QWidget()
        dialog_help_l3.setStyleSheet('background-color: white;')
        layout_help2 = QVBoxLayout()
        dialog_help_l3.setLayout(layout_help2)
        dialog_help_l2.setWidget(dialog_help_l3)

        base_path = os.path.abspath(__file__)
        folder_help = os.path.join(os.path.dirname(base_path))
        file_help = os.path.join(folder_help, 'page_help.txt')
        with open(file_help, 'r') as hf:
            lines = hf.readlines()
            label_help = QLabel('\n'.join(lines))
            label_help.setTextFormat(Qt.TextFormat.RichText)                # Enable HTML format
            label_help.setWordWrap(True)
            layout_help2.addWidget(label_help)

        dialog_help.exec()




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

#_______________________________________________________________________________________________
#_______________________________________________________________________________________________
#_______________________________________________________________________________________________2-1 CLASS: Toolbar - protein
#_______________________________________________________________________________________________
#_______________________________________________________________________________________________
class AlignmentWorker(QThread):
    finished = Signal(list)
    error = Signal(str)

#_______________________________________________________________________________________________1 DEF: main
#_______________________________________________________________________________________________
    def __init__(self, uid, sequences, button_aln, parent=None):
        super().__init__(parent)
        self.sequences = sequences
        self.button_aln = button_aln
        self.uid = uid

#_______________________________________________________________________________________________2 DEF: connect from DEF: run_alignment
#_______________________________________________________________________________________________
    def run(self):
        aligned_seq = []
        try:
            # 1 go to output folder
            base_path = os.path.dirname(os.path.abspath(__file__))
            output_folder = os.path.join(base_path, '..', 'output_files')
            os.makedirs(output_folder, exist_ok=True)

            # 2 create specific output folder
            session_folder = os.path.join(output_folder, self.uid)
            os.makedirs(session_folder, exist_ok=True)

            # 3 create files
            fasta_file = os.path.join(session_folder, f'{self.uid}.fasta')
            aln_file = os.path.join(session_folder, f'{self.uid}.aln')
            with open(fasta_file, 'w') as f:
                for name, seq in self.sequences:
                    seq = seq.replace('-','')
                    if not name or not seq:
                        self.error.emit('Missing sequence name or sequence.')
                        return
                    if not re.match(r'^[A-Za-z]+$', seq):
                        self.error.emit(f'Invalid characters found in sequence: {name}. **Special characters are not allowed.')
                        return
                    if 'M' not in seq.upper():
                        self.error.emit(f'Sequence {name} does not contain Met. Ensure it is a valid protein sequence.')
                        return
                    f.write(f'>{name}\n{seq}\n')

            # 4 sequence alignment
            # 4-1 run alignment
            if self.button_aln.text() == 'MAFFT':
                print('#  Running MAFFT v7.526	(parameters: --anysymbol, --genafpair,--maxiterate 10000)...')
                mafft_path = "./external_tools/mafft_linux/mafft_7_525_withoutextensions/core/mafft"
                with open(aln_file, 'w') as out:
                    subprocess.run([mafft_path, '--anysymbol', '--genafpair', '--maxiterate', '10000', fasta_file], check=True, stdout=out, stderr=subprocess.DEVNULL)

            elif self.button_aln.text() == 'ClustalO':
                print('#  Running ClustalO v1.2.2 (parameters: default)...')
                clustalo_path = "./external_tools/clustalo_linux/clustalo"
                with open(aln_file, 'w') as out:
                    subprocess.run([clustalo_path, '-i', fasta_file, '-o', aln_file, '--dealign', '--force'], check=True, stdout=out, stderr=subprocess.DEVNULL)

            # 4-2 store aligned sequences from aln_file to variable aligned_seq
            for record in SeqIO.parse(aln_file, 'fasta'):
                aligned_seq.append((record.id, str(record.seq)))

            self.finished.emit(aligned_seq)

        except Exception as e:
            self.error.emit(str(e))




#_______________________________________________________________________________________________
#_______________________________________________________________________________________________
#_______________________________________________________________________________________________2-2 CLASS: Toolbar - protein
#_______________________________________________________________________________________________
#_______________________________________________________________________________________________
class AlignmentDialog(QDialog):

#_______________________________________________________________________________________________1 DEF: main
#_______________________________________________________________________________________________    
    def __init__(self, uid, sequences, button_aln, callback_fn, group=None, layout=None, output_file=None, return_only=False, seq_map=None, parent=None):
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
        self.uid = uid

        # ---1 Progress Bar
        self.progress_bar = QProgressBar()
        self.progress_bar.setRange(0, 0)
        self.layout.addWidget(self.progress_bar)
        self.cancel_button = QPushButton("Cancel")
        self.cancel_button.clicked.connect(self.reject)
        self.layout.addWidget(self.cancel_button)

        # ---2 Connect to Class: AlignmentWorker (Backend thread signaling)
        self.worker = AlignmentWorker(self.uid, sequences, button_aln)
        self.worker.finished.connect(self.on_finished)
        self.worker.error.connect(self.show_error)
        self.worker.start()

#_______________________________________________________________________________________________2 DEF: signal_Done
#_______________________________________________________________________________________________  
    def on_finished(self, aligned):
        self.accept()
        try:
            self.callback_fn(uid=self.uid, 
                aligned_seq=aligned,
                group=self.group,
                layout=self.layout_param,
                button_aln=None,  # already used, optional here
                output_file=self.output_file,
                return_only=self.return_only,
                seq_map=self.seq_map
            )
        except Exception as e:
            QMessageBox.critical(self, "Callback Error", str(e))

#_______________________________________________________________________________________________3 DEF: signal_Error
#_______________________________________________________________________________________________  
    def show_error(self, message):
        QMessageBox.critical(self, "Alignment Error", message)
        self.reject()




#_______________________________________________________________________________________________
#_______________________________________________________________________________________________
#_______________________________________________________________________________________________2-3 CLASS: Toolbar - protein
#_______________________________________________________________________________________________
#_______________________________________________________________________________________________
class protein(QWidget):

#_______________________________________________________________________________________________1 DEF: main
#_______________________________________________________________________________________________ 
    def __init__(self):
        super().__init__()

# --------------------------------------------Others
        # Initiation
        self.is_protein = True
        self.seq_rows = []                                              # SEQUENCE ROWS BY GROUP
        self.groups = []                                                # SEQUENCE DETAILS DICT
        self.widget_toggles = []                                        # FOR MENU2_HIDE TOGGLES
        self.is_alignall = False
        self.is_searchseq = False
        self.widget_global = None
        self.similarity_color = 'darkmagenta'
        self.align_blue = False
        self.prediction_text = ""
        self.saveconservation = False

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
        self.layout_protein_l3.setAlignment(Qt.AlignTop)
        widget_protein_l3.setLayout(self.layout_protein_l3)
         # ---Add widgets to parent widget
        widget_protein_l2.setWidget(widget_protein_l3)

# --------------------------------------------Sub-elements
        # ---Drawing Canvas (In Layer 3)
        self.canvas = DrawingCanvas()                                           # CONNECT TO FILE 2: drawingcanvas.py
        self.canvas.setToolTip('Left-click to draw, Right-click to erase.')
        self.canvas.setFixedHeight(40)
        self.canvas.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        self.layout_protein_l3.addWidget(self.canvas)

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
        self.slidercon.setToolTip(f'Show %Conservation for each residue position: \nTick checkbox and move slider.')
        self.slidercon.setValue(100)                                                            # Initial value: 10
        self.slidercon.setMinimum(10)
        self.slidercon.setMaximum(100)                                                          # Ticks: 1 - 10
        self.slidercon.setSingleStep(10)                                                        # Move every 10 units
        self.slidercon.setTickPosition(QSlider.TicksBelow)                                      # Tick below the slider
        self.slidercon.setFixedSize(120,20)
        self.slidercon.valueChanged.connect(self.slider_threshold)
        self.checkboxslider = QCheckBox()
        self.checkboxslider.setToolTip('Tick to activate slider.')
        self.checkboxslider.setChecked(False)
        self.checkboxslider.toggled.connect(self.handle_slider_mode_toggle)
        labelslider = QLabel('Tick to activate %Conservation slider.')

        # # 4 Status
        self.status = QLabel('')

        # # Add widgets to parent widget
        self.widget_protein_buttons = QWidget()
        self.layout_protein_buttons = QHBoxLayout()
        self.layout_protein_buttons.setContentsMargins(0,2,0,2)
        self.widget_protein_buttons.setLayout(self.layout_protein_buttons)
        self.layout_protein_buttons.addWidget(self.button1_addgroup, alignment=Qt.AlignLeft)
        self.layout_protein_buttons.addWidget(self.button2_alignall, alignment=Qt.AlignLeft)
        self.layout_protein_buttons.addWidget(self.slidercon, alignment=Qt.AlignLeft)
        self.layout_protein_buttons.addWidget(self.checkboxslider, alignment=Qt.AlignLeft)
        self.layout_protein_buttons.addWidget(labelslider, alignment=Qt.AlignLeft)
        self.layout_protein_buttons.addWidget(self.status, alignment=Qt.AlignLeft)
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


#_______________________________________________________________________________________________1 DEF: QC System TCSH
#_______________________________________________________________________________________________
    def show_tcsh_warning(self):
        msg = QMessageBox(self)
        msg.setIcon(QMessageBox.Warning)
        msg.setWindowTitle('tcsh Not Found')
        msg.setText("Warning: Missing 'tcsh' shell. \nSome features like local PSIPRED will not work. Please see README.txt to install or run PSIPRED online by selecting it when prompted after pressing the align button.")
        msg.setStandardButtons(QMessageBox.Ok)
        msg.show()

#_______________________________________________________________________________________________2 DEF: Set color for alignment
#_______________________________________________________________________________________________
    def set_similarity_color(self, color, group_idx=None):
        self.similarity_color = color
        # Reapply coloring based on current alignment mode
        if getattr(self, 'is_alignall', False) and hasattr(self, 'seq_map'):
            self.color_code_seq(seq_map=self.seq_map, mode="all")
        elif group_idx is not None:
            self.color_code_seq(mode="group", group_idx=group_idx)

#_______________________________________________________________________________________________3 DEF: Toolbar - Search Sequence
#_______________________________________________________________________________________________
    def search_sequences(self, widget_search_seq):
        # 1 Search from a search bar
        search = widget_search_seq.strip().upper()

        # 2 If search bar is empty
        if not search:
            # Restore the background color according to bg_color (align or alignall mode)                                                       
            if self.is_searchseq:
               for group in self.groups:
                    for entry in group['widget_seq']:
                        for lbl in entry['seq_letters']:
                            bg_color = lbl.property("bg_color")
                            if bg_color is None:
                                bg_color = ""
                            lbl.setStyleSheet(f'background-color: {bg_color};')
            # By default, self.is_searchseq is False (with no mode set, bg color is cleared)
            else:
                for group in self.groups:
                    for entry1 in group['widget_seq']:
                        for label in entry1['seq_letters']:
                            label.setStyleSheet("")  # clear search color
                self.is_searchseq = False
            return

        # 3 Set self.is_searchseq to True
        self.is_searchseq = True

        for group in self.groups:
        # 4 Clear background color
            for entry1 in group['widget_seq']:
                aligned_seq = ''.join(entry1['seq']).upper()
                label_row = entry1['seq_letters']
                for label in label_row:
                    label.setStyleSheet("")

        # 5 Skip gaps - still match even though separated by gaps
                non_gap_seq = ""
                pos_map = []
                # remove gaps (-)
                for idx, char in enumerate(aligned_seq):
                    if char != "-":
                        non_gap_seq += char
                        pos_map.append(idx)
                # search in sequences without gaps
                for i in range(len(non_gap_seq) - len(search) + 1):
                # map back the search
                    if non_gap_seq[i:i + len(search)] == search:
                        for j in range(len(search)):
                            align_idx = pos_map[i + j]
                            label_row[align_idx].setStyleSheet("background: #FFA500;")


#_______________________________________________________________________________________________4-1 DEF: MENU VIEW - Toggles
#_______________________________________________________________________________________________
    # 1 SHOW ALL
    def view_show_all(self):
        for widget in self.widget_toggles:
            widget.show()


    # 2 HIDE TOGGLES
    def view_hide_toggles(self):
        for widget in self.widget_toggles:
            widget.hide()


#_______________________________________________________________________________________________4-2 DEF: MENU VIEW - Consensus Mode
#_______________________________________________________________________________________________
    # 1 SET INTERFACE
    def adjust_consensusmode(self):
# --------------------------------------------Main
        dialog = QDialog(self)
        dialog.setStyleSheet('QDialog {background-color: white;}')
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
        self.slider.setValue(int(getattr(self, 'consensus_threshold', 0.5) * 100))
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


    # 2 APPLY THRESHOLD TO SEQUENCES
    def apply_new_consensus_threshold(self):
# . . .  CLEAR GUI . . . # . . . CLEAR DICT . . .
        for group in self.groups:
            layout = group['main_layout_seq']
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
            layout = group['main_layout_seq']                                                   
            target_consensus = group.get('consensus_seq')
            if not target_consensus or target_consensus == ref_consensus:
                continue    # skip the for loop
            matches = sum(
                1 for i in range(1, len(ref_consensus))  # start at 1
                if i < len(target_consensus) and ref_consensus[i] == target_consensus[i])
            percent_conservation = (matches / (len(ref_consensus))) * 100
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


#_______________________________________________________________________________________________5 DEF: Set Conservation Threshold Mode
#_______________________________________________________________________________________________
    # ACTIVATE/INACTIVATE CONSERVATION THRESHOLD
    def handle_slider_mode_toggle(self, checked):                                         
        if checked:
            self.slider_threshold()
        else:
            if hasattr(self, 'seq_map') and self.seq_map:
                self.color_code_seq(seq_map=self.seq_map, mode="all")
            else:
                for idx, group in enumerate(self.groups):
                    self.color_code_seq(mode="group", group_idx=idx)


    # SAVE CONSERVATION DATA INTO A FILE
    def menu1_others_saveconservationfile_clicked(self):
        self.saveconservation = True
        self.slider_threshold()


    # IF ACTIVATED, APPLY THRESHOLLD TO INTERFACE
    def slider_threshold(self):
        if self.checkboxslider.isChecked():
            threshold = self.slidercon.value()
            lower = threshold - 10
            upper = threshold + 10
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
                
            if self.saveconservation:
                selectconservation_file = os.path.join(self.session_folder, f'{self.uid}_conservation_{threshold}.txt')
                with open(selectconservation_file, 'w') as sc:
                    sc.write(f'# Conservation Threshold: {threshold} +/-10\n')
                    header = ['Position']
                    label_refs = []
                    for group in self.groups:
                        group_name = group['lineedit_groupname'].text()
                        for idx, entry in enumerate(group['widget_seq']):
                            seq_name = f'{group_name}_seq{idx+1}'
                            header.append(seq_name)
                            label_refs.append((seq_name, entry['seq_letters']))
                    sc.write('\t'.join(header) + '\n')

                    num_positions = len(conservation_scores)
                    for pos in range(num_positions):
                        percent = conservation_scores[pos]
                        if not (lower <= percent <= upper):
                            continue  # Skip positions if outside threshold

                        row = [str(pos + 1)]
                        keep_position = False
                        for seq_name, lbl_list in label_refs:
                            if pos < len(lbl_list):
                                lbl = lbl_list[pos]
                                residue = lbl.text()
                                row.append(residue)
                                keep_position = True
                                bg_color = lbl.property("bg_color") or 'white'
                                lbl.setStyleSheet(f'background-color: {bg_color}; color: black;')

                        if keep_position:
                            sc.write('\t'.join(row) + '\n')

                QMessageBox.information(self, 'Saved', f'Project saved to {self.session_folder}')

            else:
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

            # return to default if inactivated
            self.saveconservation = False




#_______________________________________________________________________________________________6-1 DEF: Menu - Save Project
#_______________________________________________________________________________________________
    def save_project(self):
# ------1 Set directory for saved file        
        file_path, _ = QFileDialog.getSaveFileName(self, 'Save Project', '', 'Alignate Project (*.alignate)')
        if not file_path:
            return
        if not file_path.endswith('.alignate'):
            file_path += '.alignate'

        # Set window protein
        window_tag = "protein"

# ------2 Save self.groups in group_data
        serializable_groups = []
        for group in self.groups:
            group_data = {'name': group['lineedit_groupname'].text(), 'is_reference': group['checkbox_setrefgroup'].isChecked(), 'sequences': [{'header': entry['seq_header'].text(), 'sequence': ''.join(lbl.text() for lbl in entry['seq_letters'])} for entry in group['widget_seq']]}
            serializable_groups.append(group_data)

# ------3 Save modes
        saved_mode = "group"
        if self.is_alignall:
            saved_mode = "all"
        else:
            saved_mode = "group"

# ------4 Save in states
        state = {
            'window': window_tag,
            'groups': serializable_groups,
            'prediction_text': self.prediction_text or "",
            'mode': saved_mode,
        }

# ------5 Save file
        with tempfile.TemporaryDirectory() as temp_dir:
            json_path = os.path.join(temp_dir, 'state.json')
            with open(json_path, 'w') as f:
                json.dump(state, f)
            with zipfile.ZipFile(file_path, 'w') as zf:
                zf.write(json_path, 'state.json')

# ------6 Saved file Message
        QMessageBox.information(self, 'Saved', f'Project saved to {file_path}')



#_______________________________________________________________________________________________6-2 DEF: Menu - Load Project
#_______________________________________________________________________________________________
    def load_project(self, file_path: str | None = None):

# ------1 Load file         
        file_path, _ = QFileDialog.getOpenFileName(self, 'Load Project', '', 'Alignate Project (*.alignate)')
        if not file_path:
            return

        # extract JSON
        with tempfile.TemporaryDirectory() as temp_dir:
            with zipfile.ZipFile(file_path, 'r') as zf:
                zf.extractall(temp_dir)
            with open(os.path.join(temp_dir, 'state.json')) as f:
                state = json.load(f)

        # window tag
        my_tag   = "protein"          #  in Codon.py use "codon"
        file_tag = state.get("window", my_tag)
        if file_tag != my_tag:
            main = self.parent()      # main window holds the other pages
            if file_tag == "protein" and hasattr(main, "window_protein"):
                main.window_protein.load_project(file_path)
                if hasattr(main, "switch_mode"):
                    main.switch_mode("protein")
            elif file_tag == "codon" and hasattr(main, "window_codon"):
                main.window_codon.load_project(file_path)
                if hasattr(main, "switch_mode"):
                    main.switch_mode("codon")
            else:
                QMessageBox.warning(self, "Unknown page", f"Open window “{file_tag}” and reload.")
            return

        # . . .  CLEAR GUI . . .
        # REMOVE EXISTING WIDGET SECONDARY STRUCTURE
        if hasattr(self, 'widget_horizontal') and self.widget_horizontal is not None:
            self.layout_protein_l4_2ndarystructure.removeWidget(self.widget_horizontal)     # remove from layout (optional)
            self.widget_horizontal.setParent(None)                                          # remove from parent GUI hierarchy
            self.widget_horizontal.deleteLater()                                            # schedule for safe deletion by Qt event loop
            self.widget_horizontal = None                                                   # no widget global

# ------2 Clear current groups
        for group in self.groups[:]:
            self.button2_removegroup_clicked(group['widget_group'])

# ------3 Restore objects
        saved_mode = state.get('mode')
        self.prediction_text = state.get('prediction_text', "")

        # self.groups
        for group_data in state['groups']:
            self.button1_addgroup_clicked()
            group = self.groups[-1]
            group['lineedit_groupname'].setText(group_data['name'])
            group['checkbox_setrefgroup'].setChecked(group_data['is_reference'])
            for seq in group_data['sequences']:

# ------4 Run: Display sequences on GUI
                self.add_sequences_toGUI(group, group['layout_seq'], seq['header'], seq['sequence'])

# ------5 Create output folder
        self.session_folder = os.path.join(os.path.dirname(__file__), '..', 'output_files', 'restored_session')
        os.makedirs(self.session_folder, exist_ok=True)
        self.uid = 'restored'

        # sequence mapping for Alignall
        seq_map = [(group_idx, entry_idx) for group_idx, group in enumerate(self.groups) for entry_idx, _ in enumerate(group['widget_seq'])]

        try:
            if self.prediction_text:
# ------6 Group mode - Run: Group consensus; Color sequences
                if saved_mode == "group":
                    self.is_alignall = False

                    # Group consensus
                    for group in self.groups:
                        self.get_consensus_aln(group=group)

                    # Color sequences
                    self.color_code_seq(mode="group")

            # . . .  CLEAR GUI . . .
                    # REMOVE EXISTING WIDGET SECONDARY STRUCTURE
                    if hasattr(self, 'widget_horizontal') and self.widget_horizontal is not None:
                        self.layout_protein_l4_2ndarystructure.removeWidget(self.widget_horizontal)     # remove from layout (optional)
                        self.widget_horizontal.setParent(None)                                          # remove from parent GUI hierarchy
                        self.widget_horizontal.deleteLater()                                            # schedule for safe deletion by Qt event loop
                        self.widget_horizontal = None                                                   # no widget global

                    print("This is group-based alignment file. Skipping secondary structure prediction and colors.")

# ------7 All mode - Color sequences; Group %Conservation, Global consensus, 2ndary structure & its %Conservation   
                else:
                    self.is_alignall = True

                    # Group consensus
                    for group in self.groups:
                        self.get_consensus_aln(group=group)

                    # Color sequences
                    self.color_code_seq(seq_map=seq_map, mode="all")

                    # Group %Conservation
                    ref_consensus = None
                    for group in self.groups:
                        if group['checkbox_setrefgroup'].isChecked():
                            ref_consensus = group.get('consensus_seq')
                            break

                    if ref_consensus:
                        for group in self.groups:
                            target_consensus = group.get('consensus_seq')
                            layout = group['main_layout_seq']
                            if not target_consensus or target_consensus == ref_consensus:
                                continue
                            # Calculate match %
                            matches = sum(
                                1 for i in range(1, len(ref_consensus))  # start at 1
                                if i < len(target_consensus) and ref_consensus[i] == target_consensus[i])
                            percent_conservation = (matches / (len(ref_consensus))) * 100
                            #matches = sum(1 for a, b in zip(ref_consensus, target_consensus) if a == b)
                            #percent_conservation = (matches / len(ref_consensus)) * 100
                            str_percent_conservation = f"{percent_conservation:.3g}%"        # 3sf

                            # Update GUI label
                            for i in range(layout.count()):
                                widget = layout.itemAt(i).widget()
                                if widget and widget.objectName() == "consensus_row":
                                    label = widget.findChild(QLabel, "percent_conservation")
                                    if label:
                                        label.setText(str_percent_conservation)

                    # Global consensus
                    self.get_global_consensus()

                    # 2ndary structure & its %Conservation
                    region_conservation = self.compute_region_conservation(self.prediction_text)
                    self.draw_secondary_structure_to_gui(
                        self.prediction_text,
                        region_conservation=region_conservation,
                        mode="all"
                    )

# ------9 Load file message
        except Exception as e:
            print(f"Error restoring secondary structure: {str(e)}")
        QMessageBox.information(self, "Loaded", f"Project loaded from {file_path}")


#_______________________________________________________________________________________________7-1 DEF: Line 1 - +Group
#_______________________________________________________________________________________________
    def button1_addgroup_clicked(self):
# --------------------------------------------Initiation
        self.align_blue = False
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
        self.widget_protein_l4_group_l1.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Preferred)
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
        self.widget_protein_l4_group_l1_seq0 = QWidget()
        self.layout_protein_l4_group_l1_seq0 = QVBoxLayout()
        self.widget_protein_l4_group_l1_seq0.setFixedHeight(160)
        self.layout_protein_l4_group_l1_seq0.setSpacing(0)
        self.widget_protein_l4_group_l1_seq0.setLayout(self.layout_protein_l4_group_l1_seq0)
        self.layout_protein_l4_group_l1_seq0.setContentsMargins(5,0,0,0)
        self.layout_protein_l4_group_l1.addWidget(self.widget_protein_l4_group_l1_seq0)

        widget_protein_l4_group_l1_seq1 = QScrollArea()
        widget_protein_l4_group_l1_seq1.setWidgetResizable(True)
        widget_protein_l4_group_l1_seq1.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        widget_protein_l4_group_l1_seq1.setVerticalScrollBarPolicy(Qt.ScrollBarAsNeeded)
        self.layout_protein_l4_group_l1_seq0.addWidget(widget_protein_l4_group_l1_seq1)

        self.widget_protein_l4_group_l1_seq = QWidget()
        layout_protein_l4_group_l1_seq = QVBoxLayout()
        self.widget_protein_l4_group_l1_seq.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Minimum)
        layout_protein_l4_group_l1_seq.setSpacing(0)
        layout_protein_l4_group_l1_seq.setAlignment(Qt.AlignTop)
        self.widget_protein_l4_group_l1_seq.setLayout(layout_protein_l4_group_l1_seq)
        layout_protein_l4_group_l1_seq.setContentsMargins(5,0,0,0)
        widget_protein_l4_group_l1_seq1.setWidget(self.widget_protein_l4_group_l1_seq)   

        # INITIATE DICTIONARY
        group = {
            'widget_group': self.widget_protein_l4_group_l1,        # WIDGET GROUP
            'lineedit_groupname': lineedit_groupname,               # GROUP NAME
            'layout_seq': layout_protein_l4_group_l1_seq,           # LAYOUT SEQUENCE
            'main_layout_seq': self.layout_protein_l4_group_l1_seq0,# MAIN LAYOUT FOR LAYOUT SEQUENCE
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


#_______________________________________________________________________________________________7-2 DEF: Line 1 - Exclusive Reference Checkbox
#_______________________________________________________________________________________________
    def handle_reference_group_toggle(self, selected_checkbox):
        if selected_checkbox.isChecked():
            for group in self.groups:
                box = group['checkbox_setrefgroup']
                if box != selected_checkbox:
                    box.setChecked(False)


#_______________________________________________________________________________________________7-3 DEF: Line 1 - -Group
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

        # 4 Reset color
        for group in self.groups:
            for entry in group['widget_seq']:
                for lbl in entry['seq_letters']:
                    lbl.setStyleSheet("")
                    lbl.setProperty("bg_color", None)


#_______________________________________________________________________________________________7-4 DEF: Line 1 - +Sequence
#_______________________________________________________________________________________________
    def button3_addseq_clicked(self, layout):
# --------------------------------------------Initiation
        self.align_blue = False

# . . .  CLEAR GUI . . . # . . . CLEAR DICT . . .
        for group in self.groups:
            if group['main_layout_seq'] == layout:
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
        self.widget_protein_l4_group_l1_seq_dialoginput.setStyleSheet('QDialog {background-color: white;}')
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


#_______________________________________________________________________________________________7-4-1 DEF: Line 1 - +Sequence via text
#_______________________________________________________________________________________________
    # ADD SEQUENCES TO A DIALOG BOX
    def seqtext_button1text_clicked(self, layout):
# --------------------------------------------Main
        self.widget_seq_text_inputtext_dialogbox = QDialog(self)
        self.widget_seq_text_inputtext_dialogbox.setStyleSheet('QDialog {background-color: white;}')
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


    # TRANSFER SEQUENCES FROM DIALOG BOX TO WINDOW
    def widget_seq_text_inputtext_dialogbox_button1_clicked(self, layout):
# --------------------------------------------Initiate
        self.dialog_seqlength = False

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
                    QMessageBox.warning(self, 'Invalid sequence', f'Special characters are found in sequence: {seq_line}')
                    return
                if 'M' not in seq_line.upper():
                    QMessageBox.warning(self, 'Invalid sequence', 'Sequences do not contain Met. Ensure it is a valid protein sequence.')
                sequence.append(seq_line)                           # collect all lines before line with >
                i += 1
            full_sequence = ''.join(sequence)                       # form a full sequence
            self.allinput_seq.append(full_sequence)                 # add to sequence list
                
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

        print('max_len_groupseq')
        print(len_groupseq)
        print(max_len_groupseq)
        print('---')

        # ---2 Get max length of the newly input sequences
        len_allinput_seq = []
        for seq in self.allinput_seq:
            if len(seq) > 900:
                QMessageBox.warning(self, 'Sequence Length Exceeded', 'One or more of the uploaded sequences exceed the maximum allowed length of 900. Please remove them and re-upload.')
                self.widget_seq_text_inputtext_dialogbox.accept()
                self.widget_protein_l4_group_l1_seq_dialoginput.accept()
                return
            len_allinput_seq.append(len(seq))
            print(len(seq))
        max_len_allinput_seq = max(len_allinput_seq) if len_allinput_seq else 0

        print('max_len_allinput_seq')
        print(len_allinput_seq)
        print(max_len_allinput_seq)
        print('---')

        # ---3 Compare 1 & 2 --> store the max out of 2
        max_len = max(max_len_groupseq, max_len_allinput_seq)
        print('max_len')
        print(max_len)

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
        if list_index:
            # 1 main layout
            widget_filter_seq_dialog = QDialog(self)
            widget_filter_seq_dialog.setStyleSheet('QDialog {background-color: white;}')
            widget_filter_seq_dialog.move(700,600)
            layout_filter_seq_dialog = QVBoxLayout()
            widget_filter_seq_dialog.setLayout(layout_filter_seq_dialog)
            widget_filter_seq_dialog.setWindowTitle('Filter: Potential partial sequences')

            if min_len_groupseq:
                if max_len - min_len_groupseq > 10:
                    label_filter_min = QLabel(f'**Short sequence length detected in GUI (length={min_len_groupseq}). Please remove if required.')
                    layout_filter_seq_dialog.addWidget(label_filter_min)

            label_filter_max = QLabel(f'**Longest sequence: {max_len}')
            layout_filter_seq_dialog.addWidget(label_filter_max)

            label_filter_seq_dialog = QLabel('Please check sequences below and select to remove from the new list:')
            label_filter_seq_dialog.setStyleSheet('font-weight: bold;')
            layout_filter_seq_dialog.addWidget(label_filter_seq_dialog)

            widget_filter_seq_scroll = QScrollArea()
            widget_filter_seq_scroll.setWidgetResizable(True)
            widget_filter_seq_scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAsNeeded)
            widget_filter_seq_scroll.setVerticalScrollBarPolicy(Qt.ScrollBarAsNeeded)
            layout_filter_seq_dialog.addWidget(widget_filter_seq_scroll)

            # 2 main seq layout (lines of seqs)
            widget_filter_seq_dialog2 = QWidget()
            widget_filter_seq_dialog2.setStyleSheet('background-color: white;')
            layout_filter_seq_dialog2 = QVBoxLayout()
            widget_filter_seq_dialog2.setSizePolicy(QSizePolicy.MinimumExpanding, QSizePolicy.MinimumExpanding)
            layout_filter_seq_dialog2.setSpacing(0)
            layout_filter_seq_dialog2.setContentsMargins(5,0,0,0)
            layout_filter_seq_dialog2.setAlignment(Qt.AlignTop)
            widget_filter_seq_dialog2.setLayout(layout_filter_seq_dialog2)
            widget_filter_seq_scroll.setWidget(widget_filter_seq_dialog2)
#            layout_filter_seq_dialog.addWidget(widget_filter_seq_dialog2)

            # 5 fill in each seq layout
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


            # 6 buttons
            widget_buttons_filter = QWidget()
            layout_buttons_filter = QHBoxLayout()
            widget_buttons_filter.setLayout(layout_buttons_filter)
            layout_filter_seq_dialog.addWidget(widget_buttons_filter)
            button_filter_OK = QPushButton('OK')
            button_filter_cancel = QPushButton('cancel')
            layout_buttons_filter.addWidget(button_filter_OK)
            layout_buttons_filter.addWidget(button_filter_cancel)
            layout_filter_seq_dialog.addStretch(1)

            # 7 if OK is clicked, remove all checked from all_inputheaader & all_inputseq
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
      
        self.widget_seq_text_inputtext_dialogbox.accept()
        self.widget_protein_l4_group_l1_seq_dialoginput.accept()


#_______________________________________________________________________________________________7-4-2 DEF: Line 1 - +Sequence from external File
#_______________________________________________________________________________________________
    def seqtext_button2file_clicked(self, layout):
        self.dialog_seqlength = False

# --------------------------------------------Main
        # 1 File Upload Dialog Box
        file_path, _ = QFileDialog.getOpenFileName(self, 'Select File', '', 'FASTA files (*.fasta *.fa *.txt);;All Files (*)')
        if not file_path:
            return
        # Parse the uploaded file
        try:
            self.allinput_header = []
            self.allinput_seq = []
            sequences = list(SeqIO.parse(file_path, 'fasta'))
            if not sequences:
                QMessageBox.warning(self, 'No sequences found', 'The uploaded file is empty or do not have sequences.')
                return
            for record in sequences:
                name = record.id or 'unnamed'
                seq = str(record.seq)

                if 'M' not in seq:
                    QMessageBox.warning(self, 'Invalid sequence', 'Met is not found in sequence. Please ensure only protein sequences are added.')
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
                    self.widget_protein_l4_group_l1_seq_dialoginput.accept()
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

                widget_filter_seq_scroll = QScrollArea()
                widget_filter_seq_scroll.setWidgetResizable(True)
                widget_filter_seq_scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAsNeeded)
                widget_filter_seq_scroll.setVerticalScrollBarPolicy(Qt.ScrollBarAsNeeded)
                layout_filter_seq_dialog.addWidget(widget_filter_seq_scroll)

                # 2 main seq layout (lines of seqs)
                widget_filter_seq_dialog2 = QWidget()
                widget_filter_seq_dialog2.setStyleSheet('background-color: white;')
                layout_filter_seq_dialog2 = QVBoxLayout()
                widget_filter_seq_dialog2.setSizePolicy(QSizePolicy.MinimumExpanding, QSizePolicy.MinimumExpanding)
                layout_filter_seq_dialog2.setSpacing(0)
                layout_filter_seq_dialog2.setContentsMargins(5,0,0,0)
                layout_filter_seq_dialog2.setAlignment(Qt.AlignTop)
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
        self.widget_protein_l4_group_l1_seq_dialoginput.accept()


#______________________________________________________________________________________________7-4-3 DEF: Line 1 - Filter seqs before adding to GUI
#_______________________________________________________________________________________________
    def filter_sequences(self, dialog):
        to_remove = [idx for checkbox, idx in self.filter_checkboxes if checkbox.isChecked()]
        for idx in sorted(to_remove, reverse=True):
            del self.allinput_header[idx]
            del self.allinput_seq[idx]
        dialog.accept()


#_______________________________________________________________________________________________7-4-4 DEF: Line 1 - -Sequences
#_______________________________________________________________________________________________
    def button4_removeseq_clicked(self, layout):
# . . .  CLEAR GUI . . .
        for group in self.groups:
            if group['main_layout_seq'] == layout:
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

        # 4 Reset color
        for group in self.groups:
            for entry in group['widget_seq']:
                for lbl in entry['seq_letters']:
                    lbl.setStyleSheet("")
                    lbl.setProperty("bg_color", None)

#_______________________________________________________________________________________________8-1 DEF: Line 1 - Align sequences
#__________________________________________________________________________________________
    def button5_align_clicked(self, layout):
# --------------------------------------------Initiation
        self.align_blue = True
# --------------------------------------------Main
        widget_dialogbox_align = QDialog(self)
        widget_dialogbox_align.setStyleSheet('QDialog {background-color: white;}')
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

        # 3 Radiobutton for PSIPRED (Run online or offline * with BLAST+ or single)
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
        widget_dialogbox_align.move(500,500)
        widget_dialogbox_align.exec()


#_______________________________________________________________________________________________8-2 DEF: 1 Align all sequences
#__________________________________________________________________________________________
    def button2_alignall_clicked(self):
# --------------------------------------------Others
        # Initiation
        all_seq = []
        self.seq_map = []
        self.is_alignall = True
        self.align_blue = False

# --------------------------------------------Main
        widget_dialogbox_alignall = QDialog(self)
        widget_dialogbox_alignall.setStyleSheet('QDialog {background-color: white;}')
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
        widget_dialogbox_alignall.move(500,500)
        widget_dialogbox_alignall.exec()


#_______________________________________________________________________________________________8-3 DEF: PSIPRED mode (select prior to alignment)
#__________________________________________________________________________________________
    def toggle_psipred(self):
        if self.widget_psipred_online.isChecked():
            self.widget_radiobtns_psipred2.setVisible(False)
            self.qlineedit_email.setVisible(True)
        else:
            self.qlineedit_email.setVisible(False)
            self.widget_radiobtns_psipred2.setVisible(True)


#_______________________________________________________________________________________________9 DEF: Run Alignment > etc.
#_______________________________________________________________________________________________
    # RUN ALIGNMENT
    def run_alignment(self, sequences, group=None, layout=None, button_aln=None, output_file=None, return_only=False, seq_map=None):
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
            self.status.setText('')
            return

        print('')
        print('#1 Aligning amino acid sequence...')

# --------------------------------------------Actions
# ------1 Run Alignment: MAFFT/ClustalO
        dlg = AlignmentDialog(uid, sequences, button_aln, self.continue_run_alignment, group=group, layout=layout, output_file=output_file, return_only=return_only, seq_map=seq_map)
        dlg.exec()


    # RUN ALIGNMENT DIALOG
    # CONTINUE RUN ALIGNMENT
    def continue_run_alignment(self, uid, aligned_seq, group=None, layout=None, button_aln=None, output_file=None, return_only=False, seq_map=None):
        print('')
        print('#2 Generating consensus sequences and calculate %Conservation...')

        # ---Initiation
        self.uid = uid

        # ---1 Mode: Alignall
        if seq_map:                                                         # created in def button2_alignall_clicked(self)
            # -- 1 check if reference group is empty
            for group in self.groups:
                if group['checkbox_setrefgroup'].isChecked():
                    if not group['widget_seq']:
                        QMessageBox.critical(self, 'Missing data', 'Void reference group. Discard all empty groups.')
                        self.status.setText('')
                        return

            # -- 2 Clear GUI (unaligned seq/previously aligned seq) and DICT
            for group in self.groups:
                layout = group['layout_seq']
                for i in reversed(range(layout.count())):
                    widget = layout.itemAt(i).widget()
                    if widget:
                        layout.removeWidget(widget)
                        widget.setParent(None)                              # REMOVE FROM GUI
                        widget.deleteLater()
                # Only clear these if we're keeping the group
                group['widget_seq'].clear()
                group['consensus_seq'] = None

            # -- 3 Connect - DEF: to add sequences to GUI
            for (group_idx, __), (aligned_name, aligned_seq) in zip(seq_map, aligned_seq):       
                group = self.groups[group_idx]                              # get each group
                layout = group['layout_seq']                                # get sequence layout
                if layout is None:
                    continue
                self.add_sequences_toGUI(group, layout, aligned_name, aligned_seq)

            # -- 4 Connect - DEF: Get and Display Consensus in each group
            print('#  Generate consensus sequences')
            for group in self.groups:
                self.get_consensus_aln(group, seq_map, threshold=None)

            # -- 5 not DEF: Calculate %Base conservation
            print('#  Calculate %Conservation')
            # 1 set the file for group conservation
            groupconservation_file = os.path.join(self.session_folder, f'{self.uid}_groupconservation.txt')

            # 2 get reference consensus
            ref_consensus = None
            for group in self.groups:
                if group['checkbox_setrefgroup'].isChecked():                                   
                    ref_consensus = group.get('consensus_seq')
                    break

            # 3 get group consensus and calculate %conservation
            for group in self.groups:
                group_name = group['lineedit_groupname'].text()
                layout = group['main_layout_seq']                                                   
                target_consensus = group.get('consensus_seq')
                if not target_consensus or target_consensus == ref_consensus:
                    continue    # skip the for loop

                matches = sum(
                    1 for i in range(1, len(ref_consensus))  # start at 1
                    if i < len(target_consensus) and ref_consensus[i] == target_consensus[i])
                percent_conservation = (matches / (len(ref_consensus))) * 100
                str_percent_conservation = f"{percent_conservation:.3g}%"        # 3sf

# --------------------------------------------Output_files
                with open(groupconservation_file, 'a') as c:
                    c.write(f'{group_name}\t{str_percent_conservation}\n')

            # 4 update %conservation consensus against reference consensus
                for i in range(layout.count()):
                    widget = layout.itemAt(i).widget()
                    if widget and widget.objectName() == "consensus_row":
                        label = widget.findChild(QLabel, "percent_conservation")    # group['layout_seq'] > consensus_row > percent_conservation
                        if label:
                            label.setText(str_percent_conservation)

            # -- 6 Connect - DEF Color Code & Display on GUI (Aligned Sequences)
            self.color_code_seq(seq_map=self.seq_map, mode="all")

            # -- 7 Connect - DEF Get & Display on GUI (Global Consensus)
            global_consensus = self.get_global_consensus()
            if global_consensus:
                print('')
                print('#3 Generating secondary structure with PSIPRED...')


            # -- 8 Connect - DEF Run PSIPRED: Build Secondary Structure
                if self.widget_psipred_offline.isChecked():
                    psipred_dir = os.path.join(self.base_path, '..', 'external_tools', 'psipred')
                    self.prediction_text = self.build_secondary_structure_offline(psipred_dir, mode="all")
                else:
                    self.prediction_text = self.build_secondary_structure_online(mode="all")

            # -- 9 Connect - DEF Compute & Display (% Conservation based on PSIPRED Output)
                region_conservation = self.compute_region_conservation(self.prediction_text, mode="all")

            # -- 10 Connect - DEF Display on GUI (PSIPRED Output)
                self.draw_secondary_structure_to_gui(self.prediction_text, region_conservation=region_conservation, mode="all")


        # ---2 Mode: Align
        else:
            if layout is not None:
            # -- 1 Clear GUI (unaligned seq/previously aligned seq) and DICT
                for i in reversed(range(layout.count())):
                    widget = layout.itemAt(i).widget()
                    if widget:
                        layout.removeWidget(widget)
                        widget.setParent(None)                                          # removes from layout & GUI hierarchy but not delete
                        widget.deleteLater()                                            # REMOVE FROM GUI # ***** to amend ***** originally absent

            # -- 2 Connect - DEF: to add sequences to GUI
            if group is not None:                                                       # group is defined from def: button5_align_clicked
                group['widget_seq'].clear()                                             # REMOVE FROM DICT
                for name, seq in aligned_seq:
                    self.add_sequences_toGUI(group, layout, name, seq)

            # -- 3 Connect - DEF: Get and Display Consensus in each group
                print('#  Generate consensus sequences')
                self.get_consensus_aln(group, threshold=None)

            # -- 4 Connect - DEF Color Code & Display on GUI (Aligned Sequences)
                group_idx = self.groups.index(group)
                self.is_alignall = False
                self.color_code_seq(mode="group", group_idx=group_idx)

            # -- 5 Connect - DEF Run PSIPRED: Build Secondary Structure
                if 'consensus_seq' in group:
                    print('')
                    print('#3 Generating secondary structure with PSIPRED...')
                    if self.widget_psipred_offline.isChecked():
                        psipred_dir = os.path.join(self.base_path, '..', 'external_tools', 'psipred')
                        self.prediction_text = self.build_secondary_structure_offline(psipred_dir, group, mode="group")
                    else:
                        self.prediction_text = self.build_secondary_structure_online(group, mode="group")
            # -- 6 Connect - DEF Display on GUI (PSIPRED Output)
                    self.draw_secondary_structure_to_gui(self.prediction_text, group=group, mode="group")

                else:
                    QMessageBox.warning(self, 'Missing value', 'Consensus Sequences not found. Please contact Alignate.')
                    self.status.setText('')
                    return

        print('')
        print('Finished.')
        self.status.setText('')

# --------------------------------------------Connect
        if self.checkboxslider.isChecked():
            self.slider_threshold()


#_______________________________________________________________________________________________10 Add sequences to GUI (Unaligned & Aligned)
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

        if self.align_blue == True:
            seq_header.setStyleSheet('border: 1px solid grey; border-radius: 2px; background-color: orange;')
        else:
            seq_header.setStyleSheet('border: 1px solid lightgray; border-radius: 1px; background-color: white;')

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
        current_width = self.widget_protein_l4.minimumWidth()
        if seq_pixel_length > current_width:
            self.widget_protein_l4.setMinimumWidth(seq_pixel_length)

#_______________________________________________________________________________________________11 Get & Display Consensus (by Group)
#_______________________________________________________________________________________________
    def get_consensus_aln(self, group, seq_map=None, threshold=None):
        print('#  Generating group consensus')

# . . .  CLEAR GUI . . .
        layout = group['main_layout_seq']
        if layout is None:
            return
            
        for i in reversed(range(layout.count())):
            widget = layout.itemAt(i).widget()
            if widget:
                name = widget.objectName()
                if name in ["consensus_row", "custom_conservation_block"]:
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
        layout_consensus.addSpacing(5)

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
        consensus = summary.dumb_consensus(threshold=threshold, ambiguous='X')
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
        group_name = group['lineedit_groupname'].text()

# --------------------------------------------Output_files
        consensus_file = os.path.join(self.session_folder, f'{self.uid}_groupconsensus.fasta')
        with open(consensus_file, 'a') as c:
            c.write(f'>{group_name}\n{consensus_str}\n')

# --------------------------------------------Others
        return consensus_str


#_______________________________________________________________________________________________12 Get & Display Global Consensus
#_______________________________________________________________________________________________
    def get_global_consensus(self, threshold=None):
        
        print('#  Generating global consensus')

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
        consensus = summary.dumb_consensus(threshold=threshold, ambiguous='X')
        consensus_str = str(consensus)

        for letter in str(consensus):
            lbl = QLabel(letter)
            lbl.setFixedSize(15,20)
            lbl.setAlignment(Qt.AlignCenter)
            lbl.setStyleSheet('color:gray;')
            layout_global.addWidget(lbl)

# --------------------------------------------Output_files
        global_consensus_file = os.path.join(self.session_folder, f'{self.uid}_globalconsensus.fasta')
        with open(global_consensus_file, 'w') as c:
            c.write(f'>global_consensus\n{consensus_str}')

# --------------------------------------------Others
        return str(consensus)


#_______________________________________________________________________________________________13 Color Code Aligned Sequences
#_______________________________________________________________________________________________
    def color_code_seq(self, seq_map=None, mode=None, group_idx=None):

# --------------------------------------------Define Inside Fxns
        # CALCULATE SIMILARITY %
        def get_col_similarity(seqs):
            transposed = list(zip(*seqs))
            similarity = []
            for col in transposed:
                most_common = max(set(col), key=col.count)                              # get the most common base in a column
                similarity_score = col.count(most_common) / len(col)                    # most common base / total base in a column
                similarity.append(similarity_score)                                     # append the final value into similarity
            return similarity        

        # COLOR BASED ON SIMILARITY %
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
        elif mode == "group" and group_idx is not None:
            # Color only the selected group
            group = self.groups[group_idx]
            group_seqs = [entry['seq'] for entry in group['widget_seq']]
            group_similarity = get_col_similarity(group_seqs)
            for idx_col, score in enumerate(group_similarity):
                for entry in group['widget_seq']:
                    if idx_col < len(entry['seq_letters']):
                        lbl = entry['seq_letters'][idx_col]
                        shade = similarity_to_color(score)
                        lbl.setProperty('bg_color', shade)
                        lbl.setStyleSheet(f'background-color: {shade};')


#_______________________________________________________________________________________________14 Custom Display % Conservation & Additional Residue Distribution Plot
#_______________________________________________________________________________________________
    # CALCULATE %CONSERVATION & RESIDUE DISTRIBUTION (BASED ON AA PROPERTIES)
    def custom_display_perc_cons(self, pos1, pos2, widget_checkbox):
        
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

            # Check start & end
            if start < 1:
                start = 1
            if end > len(target_consensus):
                end = len(target_consensus)

        # 4 Remove existing Widget: Custom Display % Conservative
            layout = group['main_layout_seq']
            for i in reversed(range(layout.count())):
                widget = layout.itemAt(i).widget()
                if widget and widget.objectName() == 'custom_conservation_block':
                    layout.removeWidget(widget)
                    widget.deleteLater()

        ## start & end positions are not index numbers
        # 5 Calculate %Similarity
            start_index = max(1, start)
            end_index = end
            match_count = sum(1 for i in range(start_index, end_index + 1) if i < len(ref_consensus) and i < len(target_consensus) and ref_consensus[i] == target_consensus[i])
            total = end - start + 1
            percent = (match_count / total * 100) if total else 0

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

# --------------------------------------------%Conservation based on amino acid properties
        if widget_checkbox.isChecked():
            aa_properties = [['P','G','C'],['A','V','L','I','M'],['F','Y','W'],['S','T','N','Q'],['E','D'],['R','H','K'],['X']]
            property_names = ["special","hydrophobic_aliphatic","hydrophobic_aromatic","polar","negative","positive","unknown"]

            # filter and name your groups
            if group in self.groups:
                group_names = [group['lineedit_groupname'].text() for group in self.groups]

            # initialize counts per group per property
            counts = [[0]*len(self.groups) for _ in aa_properties]

            # count
            for grp_idx, g in enumerate(self.groups):
                seq = g['consensus_seq']
                for i in range(start, end+1):
                    if i >= len(seq): 
                        continue
                    aa = seq[i]
                    for prop_idx, props in enumerate(aa_properties):
                        if aa in props:
                            counts[prop_idx][grp_idx] += 1
                            break

            # compute totals per group (sum over all properties)
            totals = [ sum(counts[prop][grp] for prop in range(len(aa_properties))) for grp in range(len(self.groups)) ]

            # build percent matrix
            percents = [[(counts[prop_idx][grp_idx] / totals[grp_idx] * 100) if totals[grp_idx] else 0.0 for grp_idx in range(len(self.groups))] for prop_idx in range(len(aa_properties))]

            # write file with header: aa_properties  group1  group2  ...
            aa_properties_perc_file = os.path.join(self.session_folder, f'{self.uid}_conservation_byaaproperties.txt')
            with open(aa_properties_perc_file, 'w') as afile:
                afile.write("# %Amino acid properties for each group\n")
                afile.write("aa_properties\t" + "\t".join(group_names) + "\n")
                for prop_idx, pname in enumerate(property_names):
                    row = [f"{pname}"] + [f"{percents[prop_idx][g]:.1f}" for g in range(len(self.groups))]
                    afile.write("\t".join(row) + "\n")

            png_path = self.save_aa_property_distribution_plot(percents, property_names, group_names)
            QMessageBox.information(self, 'Saved', f'%Conservation data saved:\nTable: {aa_properties_perc_file}\n\nBar Chart: {png_path}')


    # GENERATE DISTRIBUTION PLOT (BASED ON AA PROPERTIES)
    def save_aa_property_distribution_plot(self, percents, property_names, group_names):
        # Prepare species data from computed percents
        species_data = {group: [percents[prop_idx][grp_idx] for prop_idx in range(len(property_names))]
                        for grp_idx, group in enumerate(group_names)}

        species = list(species_data.keys())
        num_species = len(species)
        num_props = len(property_names)

        # Transpose data for plotting: props x species
        values = list(zip(*species_data.values()))

        custom_colors = [
            '#827717',  # special (purple)
            '#33691e',  # hydrophobic_aliphatic (blue)
            '#689138',  # hydrophobic_aromatic (green)
            '#c0ca33',  # polar (darker green)
            '#ffeb3b',  # negative (red)
            '#ffc107',  # positive (orange)
            '#999999'   # unknown (gray)
        ]

        # Generate Bar plot
        fig, ax = plt.subplots(figsize=(14, 6))
        width = 0.1
        x = list(range(num_species))
        for i, prop_vals in enumerate(values):
            ax.bar([pos + i * width for pos in x], prop_vals, width=width, color=custom_colors[i], label=property_names[i])

        ax.set_xticks([pos + (num_props / 2 - 0.5) * width for pos in x])
        ax.set_xticklabels(species, rotation=45, ha='right')
        ax.set_ylabel('% Amino Acid Property Composition')
        ax.set_title('Amino Acid Property Distribution')
        ax.legend(title='Amino Acid Properties', loc='center left', bbox_to_anchor=(1.02, 0.5), borderaxespad=0.)
        ax.grid(True, linestyle='--', alpha=0.5)
        plt.tight_layout()

        # Set plot as PNG file
        png_path = os.path.join(self.session_folder, f'{self.uid}_conservation_byaaproperties.png')
        plt.savefig(png_path, dpi=300)
        plt.close()

        return png_path


#_______________________________________________________________________________________________15 PSIPRED: BUILD SECONDARY STRUCTURE
#________________________________________________________________________________________
    # 1 RUN OFFLINE
    def build_secondary_structure_offline(self, psipred_dir, group=None, mode=None):

        print('#  Offline: Running PSIPRED V4 with database: BLAST+')

# . . .  CLEAR GUI . . .
        # REMOVE EXISTING WIDGET SECONDARY STRUCTURE
        if hasattr(self, 'widget_horizontal') and self.widget_horizontal is not None:
            self.layout_protein_l4_2ndarystructure.removeWidget(self.widget_horizontal)     # remove from layout (optional)
            self.widget_horizontal.setParent(None)                                          # remove from parent GUI hierarchy
            self.widget_horizontal.deleteLater()                                            # schedule for safe deletion by Qt event loop
            self.widget_horizontal = None                                                   # no widget global

# --------------------------------------------Action
        # ---1 Set files
        horiz_file = f"{self.session_folder}/{self.uid}_refseq.horiz"

# -------------------------- for aligned reference sequence ---------------------------
        if mode == "all":
            for group in self.groups:
                if group['checkbox_setrefgroup'].isChecked():   
                    refseq = group['widget_seq'][0]['seq']
                    fasta_file = os.path.join(self.session_folder, f'{self.uid}_refseq.fasta')
                    with open(fasta_file, 'w') as fasta:
                        fasta.write(f'>ref_seq\n{refseq}\n')    # write to fasta file
        elif mode == "group":
            refseq = group['widget_seq'][0]['seq']
            fasta_file = os.path.join(self.session_folder, f'{self.uid}_refseq.fasta')
            with open(fasta_file, 'w') as fasta:
                fasta.write(f'>ref_seq\n{refseq}\n')    # write to fasta file
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
                QMessageBox.critical(self, 'PSIPRED error', f'PSIPRED with BLAST+ failed. Please ensure PSIPRED is installed. If yes: {e}')
                self.status.setText('')
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
                QMessageBox.critical(self, 'PSIPRED Error', f'PSIPRED Single failed. Please ensure PSIPRED is installed. If yes: {e}')
                self.status.setText('')
                raise e                

        # ---3 Extract data from output: horiz file
        if not os.path.exists(horiz_file):
            #QMessageBox.critical(self, 'Error', f'Expected PSIPRED output file {horiz_file} not found!')
            self.status.setText('')
            raise RuntimeError(f"PSIPRED output file {horiz_file} not found!\n\nPlease ensure local PSIPRED is installed.")

        with open(horiz_file, "r") as f:
            self.prediction_text = f.read()

        return self.prediction_text



    # 2 RUN ONLINE
    def build_secondary_structure_online(self, group=None, mode=None):

        print('#  Online: Running PSIPRED V4 online')

# --------------------------------------------Action
        # ---1 Set files
        horiz_file = f"{self.session_folder}/{self.uid}_refseq.horiz"

# -------------------------- for aligned reference sequence ---------------------------
        if mode == "all":
            for group in self.groups:
                if group['checkbox_setrefgroup'].isChecked():   
                    refseq = group['widget_seq'][0]['seq']
                    fasta_file = os.path.join(self.session_folder, f'{self.uid}_refseq.fasta')
                    with open(fasta_file, 'w') as fasta:
                        fasta.write(f'>ref_seq\n{refseq}\n')    # write to fasta file
        else:
            refseq = group['widget_seq'][0]['seq']
            fasta_file = os.path.join(self.session_folder, f'{self.uid}_refseq.fasta')
            with open(fasta_file, 'w') as fasta:
                fasta.write(f'>ref_seq\n{refseq}\n')    # write to fasta file
# -------------------------- for aligned reference sequence ---------------------------

        with open(fasta_file, 'r') as fasta:
            lines = fasta.readlines()
        filtered_lines = [line for line in lines if not line.strip().startswith('>')]
        with open(fasta_file, 'w') as fasta:
            fasta.writelines(filtered_lines)

        # 1 Request from server
        user_email = self.qlineedit_email.text()
        if not user_email or '@' not in user_email:
            QMessageBox.warning(self, "Missing Email", "Please enter a valid email address to run PSIPRED online.")
            self.status.setText('')
            return

        fasta_file_name = os.path.basename(fasta_file)
        url = 'https://bioinf.cs.ucl.ac.uk/psipred/api/submission.json'

        with open(fasta_file, 'rb') as fasta:
            payload = {'input_data': (fasta_file_name, fasta)}      # (name for server to recognize as, open(inputfile, ..))
            data = {'job': 'psipred', 'submission_name': fasta_file_name, 'email': user_email}
            r = requests.post(url, data=data, files=payload)

            # Error trap
            if r.status_code != 200 and r.status_code != 201:
                print("Empty response received from server.")
                print("Status code:", r.status_code)
                raise Exception("PSIPRED submission failed: " + r.text)

            response_data = json.loads(r.text)

        # 2 Get from server
        while True:
            result_url = "https://bioinf.cs.ucl.ac.uk/psipred/api/submission/" + response_data['UUID']
            r = requests.get(result_url, headers={"Accept":"application/json"})

            if not r.text.strip():
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

#---------------------------------------------Others
        return self.prediction_text



    # 3 DRAW PSIPRED OUTPUT ONTO GUI
    def draw_secondary_structure_to_gui(self, prediction_text, region_conservation=None, group=None, mode=None):

        print('#  Drawing secondary structure from PSIPRED on GUI')

# . . .  CLEAR GUI . . .
        # REMOVE EXISTING WIDGET SECONDARY STRUCTURE
        if hasattr(self, 'widget_horizontal') and self.widget_horizontal is not None:
            self.layout_protein_l4_2ndarystructure.removeWidget(self.widget_horizontal)     # remove from layout (optional)
            self.widget_horizontal.setParent(None)                                          # remove from parent GUI hierarchy
            self.widget_horizontal.deleteLater()                                            # schedule for safe deletion by Qt event loop
            self.widget_horizontal = None                                                   # no widget global

# --------------------------------------------Main
        # --- 1 Widgets
        self.widget_horizontal = QWidget()
        layout_horizontal = QHBoxLayout()
        layout_horizontal.setContentsMargins(0, 0, 0, 0)
        layout_horizontal.addSpacing(0)
        self.widget_horizontal.setLayout(layout_horizontal)
        self.layout_protein_l4_2ndarystructure.addWidget(self.widget_horizontal, alignment=Qt.AlignLeft)

        invisible_checkbox = QCheckBox()
        invisible_checkbox.setEnabled(False)
        invisible_checkbox.setStyleSheet('background: transparent; border: none;')
        layout_horizontal.addWidget(invisible_checkbox)

        invisible_label = QLabel('')
        invisible_label.setFixedSize(117, 20)
        layout_horizontal.addWidget(invisible_label, alignment=Qt.AlignLeft)

        # --- 2 Get the final aligned secondary structure
        aa, pred = '', ''
        for line in prediction_text.splitlines():
            if 'AA' in line:
                parts = line.strip().split()
                if len(parts) > 1:
                    aa += parts[1]
            elif 'Pred' in line:
                parts = line.strip().split()
                if len(parts) > 1:
                    pred += parts[1]

        refseq = ''
        if mode == "group":
            if group.get('widget_seq') and len(group['widget_seq']) > 0:
                refseq = group['widget_seq'][0]['seq']
        elif mode == "all":
            for group in self.groups:
                if group['checkbox_setrefgroup'].isChecked():
                    refseq = group['widget_seq'][0]['seq']
                    break

        final_pred = []
        pred_idx = 0
        seq_len = len(refseq)

        for i, letter in enumerate(refseq):
            if letter == '-':
                # Look for valid neighbors
                left = final_pred[i - 1] if i > 0 else None
                right = pred[pred_idx] if i + 1 < seq_len and refseq[i + 1] != '-' else None

                # If both neighbors agree, use that structure
                if left and right and left == right:                # in the middle of similar states
                    final_pred.append(left)
                elif left:                                          # last residue
                    final_pred.append(left)
                elif right:                                         # first residue
                    final_pred.append(right)
                else:
                    final_pred.append('C')                          # fallback
            else:
                final_pred.append(pred[pred_idx])
                pred_idx += 1

        # --- 3 Draw the secondary structure
        fig_width = len(final_pred) * 0.15
        fig, ax = plt.subplots(figsize=(fig_width, 0.4), dpi=100)
        i = 0
        while i < len(final_pred):
            ss = final_pred[i]
            start = i
            while i < len(final_pred) and final_pred[i] == ss:
                i += 1
            end = i  # inclusive index

            region_length = end - start  # +1 for inclusive range

            # --- Draw the shapes ---
            if ss == 'H':
                rect = Rectangle((start, 0.1), region_length, 0.8, linewidth=1, edgecolor='orange', facecolor='orange')
                ax.add_patch(rect)
            elif ss == 'E':
                arrow = FancyArrow(start, 0.5, region_length - 0.2, 0, width=0.3, length_includes_head=True, head_width=0.5, head_length=0.3, color='blue')
                ax.add_patch(arrow)
            elif ss == 'C':
                ax.plot([start, end], [0.5, 0.5], color='gray', linewidth=1.2)

        ax.set_xlim(0, len(final_pred))
        ax.set_ylim(0, 1)
        ax.axis('off')
        fig.tight_layout(pad=0)
        buffer = BytesIO()
        fig.savefig(buffer, format='png', transparent=True)
        plt.close(fig)

        pixmap = QPixmap()
        pixmap.loadFromData(buffer.getvalue())
        header = "\t" + "\t".join(f"{g} (%)" for g in region_conservation[0].get("group_scores", {}).keys()) if region_conservation else ""
        body = "\n".join(f"{reg['type']} ({reg['start']+1}-{reg['end']+1}):\t" + "\t".join(f"{v:.1f}" for v in reg.get("group_scores", {}).values()) for reg in region_conservation) if region_conservation else ""
        if not self.align_blue:
            label = StructureLabel(pixmap, "Position for amino acid residues:\n" + header + "\n" + body)
        else:
            group_name = group['lineedit_groupname'].text()
            label = StructureLabel(pixmap, group_name)

        layout_horizontal.addWidget(label, alignment=Qt.AlignLeft)

#-----------------------------------------------------Output files
        ss_conservation_file = os.path.join(self.session_folder, f'{self.uid}_ssconservation.txt')
        with open(ss_conservation_file, 'w') as ss:
            ss.write("Position for amino acid residues:\n" + header + "\n" + body)


#_______________________________________________________________________________________________16 PSIPRED: Compute SS Region Conservation
#________________________________________________________________________________________
    def compute_region_conservation(self, prediction_text, mode=None):

        print('#  Calculating %Conservation on secondary structure')

        ref_cons = None
        for group in self.groups:
            if group['checkbox_setrefgroup'].isChecked():
                ref_cons = group.get('consensus_seq')
                break
        if not ref_cons:
            return []

        aa, pred = '', ''
        for line in prediction_text.splitlines():
            if 'AA:' in line:
                parts = line.strip().split()
                if len(parts) > 1:
                    aa += parts[1]
            elif 'Pred:' in line:
                parts = line.strip().split()
                if len(parts) > 1:
                    pred += parts[1]

        refseq = None
        for group in self.groups:
            if group['checkbox_setrefgroup'].isChecked():
                refseq = group['widget_seq'][0]['seq']

        final_pred = []
        pred_idx = 0
        for letter in refseq:
            if letter == '-':
                final_pred.append('C')
            else:
                final_pred.append(pred[pred_idx])
                pred_idx += 1

        regions = []
        current_type = None
        start = None
        for i, ss in enumerate(final_pred):
            if ss in ['H', 'E']:
                if current_type != ss:
                    if current_type and start is not None:
                        regions.append((current_type, start, i - 1))
                    current_type = ss
                    start = i
            else:
                if current_type:
                    regions.append((current_type, start, i - 1))
                    current_type = None
                    start = None
        if current_type and start is not None:
            regions.append((current_type, start, len(final_pred) - 1))

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
                match_count = sum(
                    1 for i in range(start, end + 1)
                    if i < len(target_cons) and i < len(ref_cons) and target_cons[i] == ref_cons[i]
                )
                total = end - start + 1
                percent = (match_count / total) * 100 if total else 0
                entry["group_scores"][group_name] = percent
            result.append(entry)

        return result



#_______________________________________________________________________________________________
#_______________________________________________________________________________________________
#_______________________________________________________________________________________________4 CLASS: PSIPRED: Display %Conservation on SS in interface
#_______________________________________________________________________________________________
#_______________________________________________________________________________________________
class StructureLabel(QLabel):
#_______________________________________________________________________________________________1 DEF: main to set tooltip
#________________________________________________________________________________________
    def __init__(self, pixmap, tooltip_text, parent=None):
        super().__init__(parent)
        self.setPixmap(pixmap)
        self.tooltip_text = tooltip_text

#_______________________________________________________________________________________________2 DEF: show tooltip
#________________________________________________________________________________________
    def enterEvent(self, event):
        QToolTip.showText(event.globalPosition().toPoint(), self.tooltip_text, self)

#_______________________________________________________________________________________________3 DEF: hide tooltip
#________________________________________________________________________________________
    def leaveEvent(self, event):
        QToolTip.hideText()



# Execute
window = main()
window.show()
app.exec()

