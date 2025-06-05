import sys, re, subprocess, tempfile, os, shutil, json, zipfile, tempfile
from PySide6.QtGui import QIcon, QPixmap, QFont, QPainter, QPen, QMouseEvent
from PySide6.QtWidgets import QToolButton, QStyleOptionSlider, QStackedWidget, QProgressBar, QFileDialog, QMessageBox, QDialog, QTextEdit, QDialogButtonBox, QLayout, QScrollArea, QSizePolicy, QApplication, QMainWindow, QWidget, QCheckBox, QLabel, QLineEdit, QPushButton, QHBoxLayout, QVBoxLayout, QSlider
from PySide6.QtCore import QSize, Qt, QPoint
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
from drawingcanvas import DrawingCanvas
from ruler import ruler, ClickableLabel

app = QApplication(sys.argv)
app.setStyle("Fusion")
icon_logo = 'images/logo_ninja1'




#_______________________________________________________________________________________________
#_______________________________________________________________________________________________1 Main Window
#_______________________________________________________________________________________________




class main(QMainWindow):
    def __init__(self):
        super().__init__()

        #1 setup basic interface
        self.setWindowTitle('Alignate')
        self.setWindowIcon(QIcon(icon_logo))
        self.resize(1000,500)   #setMinimumSize(800,500)

        #2 setup menubar and toolbar
        menu = self.menuBar()
        menu1 = menu.addMenu('File')
        menu1_load = menu1.addAction('Load Project')
        menu1_save = menu1.addAction('Save Project')
        menu1_saveas = menu1.addMenu('Save as')
        menu1_saveas_aln = menu1_saveas.addAction('.aln')
        menu1_saveas_txt = menu1_saveas.addAction('.txt')
        menu1_saveas_png = menu1_saveas.addAction('.png')
        menu1_saveas_jpg = menu1_saveas.addAction('.jpg')
        menu1_saveas_pdf = menu1_saveas.addAction('.pdf')
        menu2 = menu.addMenu('View')
        menu2_all = menu2.addAction('All')
        menu2_hide = menu2.addAction('Hide toggles')
        menu2_consensus = menu2.addAction('Consensus mode (DEFAULT: 1)')
        menu3 = menu.addMenu('Help')

        # Setup toolbar
        self.stack = QStackedWidget()
        self.setCentralWidget(self.stack)
        self.window_about = about()
        self.window_protein = protein()
        self.window_codon = QWidget()  # Placeholder
        self.stack.addWidget(self.window_about)    # index 0
        self.stack.addWidget(self.window_protein)  # index 1
        self.stack.addWidget(self.window_codon)    # index 2
        toolbar = self.addToolBar('Main Toolbar')
        toolbar.addAction('About').triggered.connect(lambda: self.stack.setCurrentWidget(self.window_about))
        toolbar.addAction('Protein').triggered.connect(lambda: self.stack.setCurrentWidget(self.window_protein))
        toolbar.addAction('Codon').triggered.connect(lambda: self.stack.setCurrentWidget(self.window_codon))

        menu1_load.triggered.connect(self.window_protein.load_project)
        menu1_save.triggered.connect(self.window_protein.save_project)
        menu2_all.triggered.connect(self.window_protein.view_show_all)
        menu2_hide.triggered.connect(self.window_protein.view_hide_toggles)
        menu2_consensus.triggered.connect(lambda: self.window_protein.apply_new_consensus_threshold())

#_______________________________________________________________________________________________
#_______________________________________________________________________________________________2 Toolbar: About
#_______________________________________________________________________________________________




class about(QWidget):
    def __init__(self):
        super().__init__()

        #1 setup elements
        ##1 texts
        text1 = "Version: 1.0.0"
        text2 = "Source: githublink"
        text3 = "Developed by Adriana as part of her Masters Dissertation"
        text4 = "The University of Edinburgh"
        label_text1 = QLabel(text1)

        label_text2 = QLabel(text2)
        label_text3 = QLabel(text3)
        label_text4 = QLabel(text4)

        for label in [label_text1, label_text2, label_text3, label_text4]:
            label.setAlignment(Qt.AlignLeft | Qt.AlignTop)
            label.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Maximum)

        ##2 logo
        image_label1 = QLabel()
        image_label1.setPixmap(QPixmap(icon_logo).scaled(100, 100, Qt.KeepAspectRatio, Qt.SmoothTransformation))
        image_label1.setFixedSize(100,100)

        # setup layout
        ##1 layer 1
        widget_about_l1 = QWidget()
        layout_about_l1 = QVBoxLayout()
        widget_about_l1.setLayout(layout_about_l1)
        layout_about_l1.setSizeConstraint(QLayout.SetFixedSize)

        ##2 layer 2
        widget_about_l2 = QWidget()
        layout_about_l2 = QVBoxLayout()
        layout_about_l2.setSpacing(0)
        layout_about_l2.setSizeConstraint(QLayout.SetFixedSize)
        widget_about_l2.setLayout(layout_about_l2)

        ##3 add widget layer 2 into layer 1 (layer 2: box for text)
        layout_about_l1.addWidget(image_label1)
        layout_about_l1.addWidget(widget_about_l2)
        layout_about_l2.addWidget(label_text1)
        layout_about_l2.addWidget(label_text2)
        layout_about_l2.addWidget(label_text3)
        layout_about_l2.addWidget(label_text4)
        self.setLayout(layout_about_l1)


#_______________________________________________________________________________________________
#_______________________________________________________________________________________________3 Toolbar: Protein
#_______________________________________________________________________________________________



class protein(QWidget):
    def __init__(self):
        super().__init__()      

        self.secondary_structure_widget = None
        self.cancelled = False
        self.widget_toggles = []

##______________________________________________________________________________________________1 Set main layout

        # 1 Set the main layout in protein widget window
        layout_protein_l1 = QVBoxLayout()               
        self.setLayout(layout_protein_l1)                  # layout_protein_layer1 (layout inside widget_l1)
        layout_protein_l1.setContentsMargins(0, 0, 0, 0)

        # 2 Set the scroll area to put in the widget protein layer 1
        widget_protein_l2 = QScrollArea()                               # widget_protein_layer2 (scroll)
        widget_protein_l2.setWidgetResizable(True)
        widget_protein_l2.setHorizontalScrollBarPolicy(Qt.ScrollBarAsNeeded)
        widget_protein_l2.setVerticalScrollBarPolicy(Qt.ScrollBarAsNeeded)
        layout_protein_l1.addWidget(widget_protein_l2)                  # add widget to layout layer1

        # 3 Set the widget 3 to put in the widget protein layer 2
        widget_protein_l3 = QWidget()                                   # widget_protein_layer3
        widget_protein_l3.setSizePolicy(QSizePolicy.MinimumExpanding, QSizePolicy.MinimumExpanding)
        self.layout_protein_l3 = QVBoxLayout()
        self.layout_protein_l3.setSpacing(0)
        self.layout_protein_l3.setContentsMargins(5, 5, 5, 5)
        widget_protein_l3.setLayout(self.layout_protein_l3)             # layout_protein_layer3
        widget_protein_l2.setWidget(widget_protein_l3)                  # add widget to widget layer2 (1 scroll can only has 1 widget)

        # 4 Add ruler
#        self.ruler = ruler()
        self.ruler = ruler(protein=self)
        self.layout_protein_l3.addWidget(self.ruler)

        # 5 Set the widget 4 to put in the widget protein layer 3
        self.widget_protein_l4 = QWidget()
        self.layout_protein_l4 = QVBoxLayout()
        self.layout_protein_l4.setSpacing(0)
        self.layout_protein_l4.setContentsMargins(0, 0, 0, 0)
        self.widget_protein_l4.setLayout(self.layout_protein_l4)
        self.layout_protein_l3.addWidget(self.widget_protein_l4, alignment=Qt.AlignTop)

##_______________________________________________________________________________________________2 Set drawing canvas


        self.canvas = DrawingCanvas()
        self.canvas.setFixedHeight(40)
        self.canvas.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        self.layout_protein_l4.addWidget(self.canvas)

##_______________________________________________________________________________________________3 


        # 5 Checkbox
        self.seq_rows = []      # store all checked rows in a tuple (used to remove sequences)

        #-- 6 Store objects (for alignments)
        self.groups = []

##______________________________________________________________________________________________2 Set main layout elements

        # 1 button1: +Group
        self.button1_addgroup = QPushButton('+Group')
        self.button1_addgroup.setFixedSize(QSize(60,28))
        self.button1_addgroup.setStyleSheet(
            """
            QPushButton {
                background-color: #00008B;
                color: white;
                font-weight: bold;                           
            }
            QPushButton:hover {
                background-color: #6495ED;
                color: white;
            }
            """
        )

        # 2 button2: Align all
        self.button2_alignall = QPushButton('Align all')
        self.button2_alignall.setFixedSize(QSize(60,28))
        self.button2_alignall.setStyleSheet(
            """
            QPushButton {
                background-color: #00008B;
                color: white;
                font-weight: bold;                           
            }
            QPushButton:hover {
                background-color: #6495ED;
                color: white;
            }
            """
        )

        # 3 Slider
        self.slidercon = QSlider(Qt.Horizontal)
        self.slidercon.setValue(100)                    # 100% = show all
        self.slidercon.setMinimum(10)
        self.slidercon.setMaximum(100)
        self.slidercon.setSingleStep(10)
        self.slidercon.setTickPosition(QSlider.TicksBelow)
        self.slidercon.setFixedSize(120,20)
        self.slidercon.valueChanged.connect(self.slider_threshold)
        self.checkboxslider = QCheckBox()
        self.checkboxslider.setChecked(False)
        self.checkboxslider.toggled.connect(self.handle_slider_mode_toggle)

        # 3 Setup layout for button1 & button2
        self.widget_protein_buttons = QWidget()
        self.layout_protein_buttons = QHBoxLayout()
        self.layout_protein_buttons.setContentsMargins(0,2,0,2)
        self.widget_protein_buttons.setLayout(self.layout_protein_buttons)
        self.layout_protein_buttons.addWidget(self.button1_addgroup, alignment=Qt.AlignLeft)
        self.layout_protein_buttons.addWidget(self.button2_alignall, alignment=Qt.AlignLeft)
        self.layout_protein_buttons.addWidget(self.slidercon, alignment=Qt.AlignLeft)
        self.layout_protein_buttons.addWidget(self.checkboxslider, alignment=Qt.AlignLeft)
        self.layout_protein_buttons.addStretch(1)

        # 4 Setup widget for secondary alignment
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

        # Check system requirement
        if shutil.which('tcsh') is None:
            self.show_tcsh_warning()

        # 5 Add elements to the layout
        self.layout_protein_l4.addWidget(self.widget_protein_buttons, alignment=Qt.AlignTop)
        self.widget_toggles.append(self.widget_protein_buttons)

        self.layout_protein_l4.addWidget(self.widget_protein_l4_2ndarystructure)

        # 6 Connect to DEF
        self.button1_addgroup.clicked.connect(self.button1_addgroup_clicked)
        self.button2_alignall.clicked.connect(self.button2_alignall_clicked)


##______________________________________________________________________________________________2 Check: tcsh (not bash) for PSIPRED

    def show_tcsh_warning(self):
        msg = QMessageBox(self)
        msg.setIcon(QMessageBox.Warning)
        msg.setWindowTitle("tcsh Not Found")
        msg.setText("Warning: 'tcsh' shell is not available.\nSome features like PSIPRED may not work.")
        msg.setStandardButtons(QMessageBox.Ok)
        msg.show()  # Non-blocking


##______________________________________________________________________________________________3 Slider checkbox

    def handle_slider_mode_toggle(self, checked):
        if checked:
            self.slider_threshold()
        else:
            self.color_code_seq(seq_map=self.seq_map)




    def view_show_all(self):
        for widget in self.widget_toggles:
            widget.show()

    def view_hide_toggles(self):
        for widget in self.widget_toggles:
            widget.hide()

    def todeladjust_consensusmode(self):
        dialog = QDialog(self)
        dialog.setWindowTitle("Set Consensus Threshold")
        layout = QVBoxLayout()
        dialog.setLayout(layout)
        label = QLabel("Threshold (0.0 – 1.0):")
        layout.addWidget(label)
        slider = QSlider(Qt.Horizontal)
        slider.setRange(0, 100)  # map 0–100 to 0.0–1.0
        slider.setValue(100)     # default 1.0
        layout.addWidget(slider)
        value_label = QLabel("1.00")
        layout.addWidget(value_label)

        def update_label(value):
            value_label.setText(f"{value / 100:.2f}")
        slider.valueChanged.connect(update_label)
        button_box = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        layout.addWidget(button_box)
        button_box.accepted.connect(dialog.accept)
        button_box.rejected.connect(dialog.reject)

        if dialog.exec() == QDialog.Accepted:
            return slider.value() / 100.0  # return float threshold
        else:
            return None  # or a default like 1.0


    def adjust_consensusmode(self):
        dialog = QDialog(self)
        dialog.setWindowTitle("Set Consensus Threshold")
        layout = QVBoxLayout(dialog)
        label = QLabel("Threshold (0.0 – 1.0):")
        layout.addWidget(label)
        slider = QSlider(Qt.Horizontal)
        slider.setRange(0, 100)
        slider.setValue(int(getattr(self, 'consensus_threshold', 1.0) * 100))
        layout.addWidget(slider)
        value_label = QLabel(f"{slider.value()/100:.2f}")
        layout.addWidget(value_label)
        slider.valueChanged.connect(lambda v: value_label.setText(f"{v/100:.2f}"))
        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        layout.addWidget(buttons)
        buttons.accepted.connect(dialog.accept)
        buttons.rejected.connect(dialog.reject)

        if dialog.exec() == QDialog.Accepted:
            threshold = slider.value() / 100.0
            self.consensus_threshold = threshold
            return threshold
        return None


    def apply_new_consensus_threshold(self):

        # 1. Remove consensus row widget (per group)
        for group in self.groups:
            layout = group['layout_seq']
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

        threshold = self.adjust_consensusmode()
        if threshold is not None:
            # global consensus
            self.get_global_consensus(threshold=threshold)
            # by group consensus
            for group in self.groups:
                self.get_consensus_aln(group, threshold=threshold)

            # calculate %base conservation
            ref_consensus = None
            for group in self.groups:
                if group['checkbox_setrefgroup'].isChecked():
                    ref_consensus = group.get('consensus_seq')
                    break
            for group in self.groups:
                layout = group['layout_seq']
                target_consensus = group.get('consensus_seq')
                if not target_consensus or target_consensus == ref_consensus:
                    continue    # Skip if reference group
                matches = sum(1 for a, b in zip(ref_consensus, target_consensus) if a == b)
                percent_conservation = (matches / len(ref_consensus)) * 100
                str_percent_conservation = f"{percent_conservation:.3g}%"
                # Find QLable inside consensus_row
                for i in range(layout.count()):
                    widget = layout.itemAt(i).widget()
                    if widget and widget.objectName() == "consensus_row":
                        label = widget.findChild(QLabel, "perc_conservation")
                        if label:
                            label.setText(str_percent_conservation)
            # calculate by region conservation
            if hasattr(self, 'prediction_text'):
                self.compute_region_conservation(self.prediction_text)










##______________________________________________________________________________________________3 Main Button: +Group

    def button1_addgroup_clicked(self):
        # 0 If global consensus is present, remove it from GUI
        if hasattr(self, 'widget_global') and self.widget_global is not None:
            self.layout_protein_l3.removeWidget(self.widget_global)     # remove from layout (optional)
            self.widget_global.setParent(None)                          # remove from parent GUI hierarchy
            self.widget_global.deleteLater()                            # schedule for safe deletion by Qt event loop
            self.widget_global = None                                   # no widget global

        # 1 Layout: Group box
        self.widget_protein_l4_group_l1 = QWidget()                     # widget_protein_l3_group_layer1
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
        self.widget_protein_l4_group_l1.setLayout(self.layout_protein_l4_group_l1)                   # layout_protein_l3_group_layer1 for widget layer1
        self.layout_protein_l4.addWidget(self.widget_protein_l4_group_l1)                            # add this widget to the main_widget_layout

        # 2 Layout: Group box - Line 1
        widget_protein_l4_group_l1_line1 = QWidget()                                            # widget for line 1
        layout_protein_l4_group_l1_line1 = QHBoxLayout()
        widget_protein_l4_group_l1_line1.setLayout(layout_protein_l4_group_l1_line1)
        self.layout_protein_l4_group_l1.addWidget(widget_protein_l4_group_l1_line1, alignment=Qt.AlignTop)

        # 3 Elements: Group box - Line 1
        totalgroups = len(self.groups)
        groupno = totalgroups + 1 if totalgroups > 0 else 1

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

        checkbox1_setrefgroup = QCheckBox()
        if len(self.groups) == 0:
            checkbox1_setrefgroup.setChecked(True)                                      # Default checked
        label2_setrefgroup = QLabel('Tick to set reference group (Default: Group 1)')
        layout_protein_l4_group_l1_line1.addWidget(lineedit_groupname, alignment=Qt.AlignLeft)
        layout_protein_l4_group_l1_line1.addWidget(button2_removegroup, alignment=Qt.AlignLeft)
        layout_protein_l4_group_l1_line1.addWidget(button3_addseq, alignment=Qt.AlignLeft)
        layout_protein_l4_group_l1_line1.addWidget(button4_removeseq, alignment=Qt.AlignLeft)
        layout_protein_l4_group_l1_line1.addWidget(button5_align, alignment=Qt.AlignLeft)
        layout_protein_l4_group_l1_line1.addWidget(checkbox1_setrefgroup, alignment=Qt.AlignLeft)
        layout_protein_l4_group_l1_line1.addWidget(label2_setrefgroup, alignment=Qt.AlignLeft)
        layout_protein_l4_group_l1_line1.addStretch(1)                      # to left-align all elements

        self.widget_toggles.append(widget_protein_l4_group_l1_line1)

        # 4 Layout: Group box - Line 2 (SEQ)
        widget_protein_l4_group_l1_seq = QWidget()
        layout_protein_l4_group_l1_seq = QVBoxLayout()
        widget_protein_l4_group_l1_seq.setLayout(layout_protein_l4_group_l1_seq)
        layout_protein_l4_group_l1_seq.setSpacing(0)
        layout_protein_l4_group_l1_seq.setContentsMargins(5,0,0,0)
        self.layout_protein_l4_group_l1.addWidget(widget_protein_l4_group_l1_seq)

        #-- 5 Store dictionary: Groups (1)
        group = {
            'widget_group': self.widget_protein_l4_group_l1,
            'lineedit_groupname': lineedit_groupname,
            'layout_seq': layout_protein_l4_group_l1_seq,
            'widget_seq': [],
            'checkbox_setrefgroup': checkbox1_setrefgroup
        }
        self.groups.append(group)

        # 7  Connect to DEF (lambda: to Specify Target)
        button2_removegroup.clicked.connect(lambda _=None, w=self.widget_protein_l4_group_l1: self.button2_removegroup_clicked(w))
        button3_addseq.clicked.connect(lambda _=None, layout=layout_protein_l4_group_l1_seq: self.button3_addseq_clicked(layout))
        button4_removeseq.clicked.connect(lambda _=None, layout=self.layout_protein_l4_group_l1: self.button4_removeseq_clicked(layout))
        button5_align.clicked.connect(lambda _=None, layout=layout_protein_l4_group_l1_seq, line1_layout=layout_protein_l4_group_l1_line1: self.button5_align_clicked(layout, line1_layout))
        checkbox1_setrefgroup.toggled.connect(lambda checked, this_box=checkbox1_setrefgroup: self.handle_reference_group_toggle(this_box))                                                 # only 1 group is allowed at a time


###_____________________________________________________________________________________________0 Checkbox: Exclusive

    def handle_reference_group_toggle(self, selected_checkbox):
        if selected_checkbox.isChecked():
            for group in self.groups:
                box = group['checkbox_setrefgroup']
                if box != selected_checkbox:
                    box.setChecked(False)

###_____________________________________________________________________________________________1 Button: Remove Group

    def button2_removegroup_clicked(self, widget):


        #--------------------------Cleaning Up Begins--------------------------#

        #### GUI: IF PRESENT, REMOVE CONSENSUS GLOBAL ###
        if hasattr(self, 'widget_global') and self.widget_global is not None:
            self.layout_protein_l3.removeWidget(self.widget_global)     # remove from layout (optional)
            self.widget_global.setParent(None)                          # remove from parent GUI hierarchy
            self.widget_global.deleteLater()                            # schedule for safe deletion by Qt event loop
            self.widget_global = None                                   # no widget global


        ### GROUPS ###
        for i, group in enumerate(self.groups):
            if group['widget_group'] == widget:
                del self.groups[i]
                break

        #--------------------------Cleaning Up Ends--------------------------#

        widget.setParent(None)
        widget.deleteLater()

###_____________________________________________________________________________________________2 Button: Add Sequences

    def button3_addseq_clicked(self, layout):


	#--------------------------Cleaning Up Begins--------------------------#
	#### GUI: IF PRESENT, REMOVE CONSENSUS ROW ###

        # 1. Remove consensus row widget (per group)
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


        #### GUI: IF PRESENT, REMOVE CONSENSUS GLOBAL ###
        if hasattr(self, 'widget_global') and self.widget_global is not None:
            self.layout_protein_l3.removeWidget(self.widget_global)     # remove from layout (optional)
            self.widget_global.setParent(None)                          # remove from parent GUI hierarchy
            self.widget_global.deleteLater()                            # schedule for safe deletion by Qt event loop
            self.widget_global = None                                   # no widget global


	### GROUPS ###
	# no action

	#--------------------------Cleaning Up Ends--------------------------#


        # EFFECTOR: LINE 2
        # 1 Dialog Box 1: Input Text or File: .txt/.fasta
        ## 1 Set layout
        self.widget_protein_l4_group_l1_seq_dialoginput = QDialog()
        self.layout_protein_l4_group_l1_seq_dialoginput = QVBoxLayout()
        self.widget_protein_l4_group_l1_seq_dialoginput.setLayout(self.layout_protein_l4_group_l1_seq_dialoginput)

        ## 2 Set elements
        ### 1 label & textbox
        widget_protein_l4_group_l1_seq_dialoginput_label = QLabel()
        widget_protein_l4_group_l1_seq_dialoginput_label.setText('Please specify method to insert sequence.')

        ### 2 buttons
        self.seqtext_button1text = QPushButton('Input Text')
        self.seqtext_button2file = QPushButton('Upload File')
        self.seqtext_button3cancel = QPushButton('Cancel')

        ## 3 Group buttons in a dialog button box
        self.widget_protein_l4_group_l1_seq_dialoginput_seqbutton = QDialogButtonBox()
        self.widget_protein_l4_group_l1_seq_dialoginput_seqbutton.addButton(self.seqtext_button1text, QDialogButtonBox.ActionRole)
        self.widget_protein_l4_group_l1_seq_dialoginput_seqbutton.addButton(self.seqtext_button2file, QDialogButtonBox.ActionRole)
        self.widget_protein_l4_group_l1_seq_dialoginput_seqbutton.addButton(self.seqtext_button3cancel, QDialogButtonBox.RejectRole)

        ## 4 Insert into layout
        self.layout_protein_l4_group_l1_seq_dialoginput.addWidget(widget_protein_l4_group_l1_seq_dialoginput_label)
        self.layout_protein_l4_group_l1_seq_dialoginput.addWidget(self.widget_protein_l4_group_l1_seq_dialoginput_seqbutton)

        ## 5 Connect to DEF
        self.seqtext_button3cancel.clicked.connect(self.widget_protein_l4_group_l1_seq_dialoginput.reject)     # cancel
        self.seqtext_button1text.clicked.connect(lambda _=None, layout= layout: self.seqtext_button1text_clicked(layout))
        self.seqtext_button2file.clicked.connect(lambda _=None, layout=layout: self.seqtext_button2file_clicked(layout))

        # 2 Dialog window: execute
        self.widget_protein_l4_group_l1_seq_dialoginput.exec()

###_____________________________________________________________________________________________2 Button: Add Sequences - Input Text

    def seqtext_button1text_clicked(self, layout):

        ###1 set main widget and layout
        self.widget_seq_text_inputtext_dialogbox = QDialog()
        layout_seq_text_inputtext_dialogbox = QVBoxLayout()
        self.widget_seq_text_inputtext_dialogbox.setLayout(layout_seq_text_inputtext_dialogbox)

        ###2 set elements: label, text input, buttons: OK and Cancel
        ####1 label
        widget_seq_text_inputtext_dialogbox_label = QLabel()
        widget_seq_text_inputtext_dialogbox_label.setText('Sequences:')
        
        ####2 input text box
        self.widget_seq_text_inputtext_dialogbox_textbox = QTextEdit()
        self.widget_seq_text_inputtext_dialogbox_textbox.setPlainText('>Sequence name\nATCGGKJKJ..')
        self.widget_seq_text_inputtext_dialogbox_textbox.setFixedSize(500,400)

        ####3 input buttons
        widget_seq_text_inputtext_dialogbox_inputbuttons = QDialogButtonBox()
        widget_seq_text_inputtext_dialogbox_button1ok = QPushButton('OK')
        widget_seq_text_inputtext_dialogbox_button2cancel = QPushButton('Cancel')
        widget_seq_text_inputtext_dialogbox_inputbuttons.addButton(widget_seq_text_inputtext_dialogbox_button1ok, QDialogButtonBox.ActionRole)
        widget_seq_text_inputtext_dialogbox_inputbuttons.addButton(widget_seq_text_inputtext_dialogbox_button2cancel, QDialogButtonBox.RejectRole)

        ###3 add elements to layout
        layout_seq_text_inputtext_dialogbox.addWidget(widget_seq_text_inputtext_dialogbox_label)
        layout_seq_text_inputtext_dialogbox.addWidget(self.widget_seq_text_inputtext_dialogbox_textbox)
        layout_seq_text_inputtext_dialogbox.addWidget(widget_seq_text_inputtext_dialogbox_inputbuttons)

        ###4 connect to def
        widget_seq_text_inputtext_dialogbox_button2cancel.clicked.connect(self.widget_seq_text_inputtext_dialogbox.reject)     # cancel
        widget_seq_text_inputtext_dialogbox_button1ok.clicked.connect(lambda _=None, layout=layout: self.widget_seq_text_inputtext_dialogbox_button1_clicked(layout))

        ###5 Dialog window: execute
        self.widget_seq_text_inputtext_dialogbox.exec()


    def widget_seq_text_inputtext_dialogbox_button1_clicked(self,layout):

        # 1 Cleanup 1: see if the input is empty
        text = self.widget_seq_text_inputtext_dialogbox_textbox.toPlainText()
        if not text:
            QMessageBox.warning(self, 'Empty input', 'Please input sequences.')
            return

        # 2 Extract sequences
        lines = [line.replace(' ', '').strip() for line in text.splitlines() if line.strip()]    # Remove all spaces, add all lines into a [] as is
        self.allinput_header = []
        self.allinput_seq = []
        i = 0                                               # begins at line 0
        seq_count = 0

        ## 1 Seq Header
        while i < len(lines):
            line = lines[i]                                 # line: header
            ### 1 Breakpoint 1: '>' (should have)
            if not line.startswith('>'):
                QMessageBox.warning(self, 'Missing value', 'Insert header (Format: >Name) for each sequence.')
                return
            ### 2 '>' has no name
            seqname = line[1:].strip()                      # var: letters after '>'
            if not seqname:                                 # if var is empty
                seqname = f"seq_{seq_count + 1}"       # name as seq_#
                seq_count += 1
            else:
                seq_count += 1
            ### 3 Collect seq header
            self.allinput_header.append(seqname)

        ## 2 Seq
            ### 1 Breakpoint 1: '>' (should not have)
            i += 1                                          # line: sequence
            if i >= len(lines) or lines[i].startswith('>'): # if no seq after
                QMessageBox.warning(self, 'Missing value', f'Missing sequences: {seqname}')
                return
            ### 2 Capitalize sequences
            sequence = []
            while i < len(lines) and not lines[i].startswith('>'):
                seq_line = lines[i].upper()
            ### 3 Breakpoint 2: No special character except (-/_) is allowed
                if not re.match(r"^[A-Z\_\-]+$", seq_line):
                    QMessageBox.warning(self, 'Invalid sequence', f'Special characters are found in sequence: {seqname}')
                    return
            ### 4 Join lines (no '>') to make a complete sequence
                sequence.append(seq_line)
                i += 1
            full_sequence = ''.join(sequence)
            ### 5 Collect seq
            self.allinput_seq.append(full_sequence)

        # 3 Hide windows
        self.widget_seq_text_inputtext_dialogbox.hide()
        self.widget_protein_l4_group_l1_seq_dialoginput.hide()

        # 4 Connect to DEF: add_sequences to GUI
        for seq_name, seq in zip(self.allinput_header, self.allinput_seq):
#            self.add_sequences_toGUI(layout, seq_name, seq)
            for group in self.groups:
                if group['layout_seq'] == layout:
                    self.add_sequences_toGUI(group, layout, seq_name, seq)
                    break

###_____________________________________________________________________________________________2 Button: Add Sequences - Input File

    def seqtext_button2file_clicked(self, layout):

        file_path, _ = QFileDialog.getOpenFileName(self, 'Seleft File', '', 'FASTA files (*.fasta *.fa *.txt);;All Files (*)')
        if not file_path:
            return
        try:
            sequences = list(SeqIO.parse(file_path, 'fasta'))                                   # SeqIO parse file as FASTA format -> convert into list
            if not sequences:                                                                   # Check 1: If there is no sequence
                QMessageBox.warning(self, 'No sequences found', 'The uploaded file is null.')
                return
            for record in sequences:
                name = record.id or 'unnamed'                                                   # Get the header. If null, name as unnamed
                seq = str(record.seq)
                for group in self.groups:
                    if group['layout_seq'] == layout:
                        self.add_sequences_toGUI(group, layout, name, seq)                      # Add sequences to GUI
        except Exception as e:
            QMessageBox.critical(self, 'Error', f'Could not load file:\n{e}')

        self.widget_protein_l4_group_l1_seq_dialoginput.hide()


###_____________________________________________________________________________________________3 Button: Remove Sequences

    def button4_removeseq_clicked(self, layout):
        #---------------------- GUI Cleanup Begins --------------------------#

        # 1. Remove consensus row widget (per group)
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

        # 2. Remove global consensus if exists
        if hasattr(self, 'widget_global') and self.widget_global is not None:
            self.layout_protein_l3.removeWidget(self.widget_global)
            self.widget_global.setParent(None)
            self.widget_global.deleteLater()
            self.widget_global = None

        #---------------------- GUI Cleanup Ends --------------------------#

        #---------------------- Sequence Row Deletion --------------------------#
        seq_removed = []

        # Identify sequences marked for deletion
        for row, checkbox in self.seq_rows:
            try:
                if checkbox.isChecked():
                    seq_removed.append((row, checkbox))
            except RuntimeError:
                continue  # widget might already be deleted

        # Actually delete rows from GUI
        for row, checkbox in seq_removed:
            row.setParent(None)
            row.deleteLater()

        # Update self.seq_rows list to keep only undeleted ones
        self.seq_rows = [
            (row, checkbox)
            for (row, checkbox) in self.seq_rows
            if (row, checkbox) not in seq_removed
        ]

        #---------------------- Group Data Cleanup --------------------------#

        # Remove deleted rows from group['widget_seq']
        rows_to_remove = set()
        for row, checkbox in seq_removed:
            try:
                rows_to_remove.add(row)
            except RuntimeError:
                continue

        for group in self.groups:
            if 'widget_seq' in group:
                group['widget_seq'] = [
                    entry for entry in group['widget_seq']
                    if entry['widget_row'] not in rows_to_remove
                ]






###_____________________________________________________________________________________________=Alignment=
###_____________________________________________________________________________________________4 Button: Align

    def button5_align_clicked(self, layout, layout_line1):

        # 1 Dialog box layout
        widget_dialogbox_align = QDialog()
        layout_dialogbox_align = QVBoxLayout()
        widget_dialogbox_align.setLayout(layout_dialogbox_align)

        # 2 Dialog box elements - Label, Buttons (method: mafft, clustalo)
        ## 1 Label
        widget_label_align = QLabel()
        widget_label_align.setText('ClustalO v1.2.2 [all platforms] (Settings: Default)\nMAFFT v7.526 [win & linux] v7.490 [mac] (Settings: --anysymbol, --genafpair, --maxiterate 10000')
        layout_dialogbox_align.addWidget(widget_label_align)

        ## 2 Buttons
        widget_buttons_align = QWidget() #QDialogButtonBox()
        layout_buttons_align = QHBoxLayout()
        widget_buttons_align.setLayout(layout_buttons_align)
        widget_button1_clustalo = QPushButton('ClustalO')
        widget_button2_mafft = QPushButton('MAFFT')
        widget_button3_cancel = QPushButton('Cancel')
        layout_buttons_align.addWidget(widget_button1_clustalo)
        layout_buttons_align.addWidget(widget_button2_mafft)
        layout_buttons_align.addWidget(widget_button3_cancel)
        layout_dialogbox_align.addWidget(widget_buttons_align)

        group = None
        for g in self.groups:                                                               # Only extract the correct group
            if g['layout_seq'] == layout:
                group = g
                break
        
        # 3 Connect to DEF
        widget_button3_cancel.clicked.connect(widget_dialogbox_align.reject)                # cancel
        widget_button2_mafft.clicked.connect(lambda _=None: (widget_dialogbox_align.accept(), self.run_alignment([(entry['seq_header'].text().strip(), ''.join(label.text() for label in entry['seq_letters']).strip()) for entry in group['widget_seq']], group, layout, widget_button2_mafft)))
        widget_button1_clustalo.clicked.connect(lambda _=None: (widget_dialogbox_align.accept(), self.run_alignment([(entry['seq_header'].text().strip(), ''.join(label.text() for label in entry['seq_letters']).strip()) for entry in group['widget_seq']], group, layout, widget_button1_clustalo)))

        widget_dialogbox_align.exec()

##_____________________________________________________________________________________________4 Main Button: Align all

    def button2_alignall_clicked(self):
        all_seq = []
        self.seq_map = []

        # 1 Collect all sequences and mapping info
        for group_idx, group in enumerate(self.groups):
            for entry_idx, entry in enumerate(group['widget_seq']):

                if not hasattr(entry['seq_header'], 'text'):
                    continue
                seq_header = entry.get('seq_header')

                if not seq_header or not seq_header.text():
                    continue  # skip deleted widgets safely
                name = entry['seq_header'].text().strip()
                seq = ''.join(label.text() for label in entry['seq_letters']).strip()
                if name and seq:
                    all_seq.append((name, seq))
                    self.seq_map.append((group_idx, entry_idx))
                else:
                    QMessageBox.warning(self, 'Error', 'Missing sequence header or sequences.')
                    return
                
        if len(all_seq) < 2:
            QMessageBox.warning(self, 'Error', 'You need at least 2 sequences to align.')
            return
        
        # 2 Dialog box
        ## 1 Layout
        widget_dialogbox_alignall = QDialog()
        layout_dialogbox_alignall = QVBoxLayout()
        widget_dialogbox_alignall.setLayout(layout_dialogbox_alignall)

        ## 2 Elements - Buttons (method: mafft, clustalo)
        ## 1 Label
        widget_label_alignall = QLabel()
        widget_label_alignall.setText('ClustalO v1.2.2 [all platforms] (Settings: Default)\nMAFFT v7.526 [win & linux] v7.490 [mac] (Settings: --anysymbol, --genafpair, --maxiterate 10000')
        layout_dialogbox_alignall.addWidget(widget_label_alignall)

        ## 2 Buttons
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

        # 3 Connect to DEF
        widget_button3_cancel.clicked.connect(widget_dialogbox_alignall.reject)                # cancel
        widget_button2_mafft.clicked.connect(lambda _=None: (widget_dialogbox_alignall.accept(), self.run_alignment(all_seq, None, None, widget_button2_mafft, return_only=False, seq_map=self.seq_map)))
        widget_button1_clustalo.clicked.connect(lambda _=None: (widget_dialogbox_alignall.accept(), self.run_alignment(all_seq, None, None, widget_button1_clustalo, return_only=False, seq_map=self.seq_map)))

        widget_dialogbox_alignall.exec()

# shared #______________________________________________________________________________________________________5 Additional fxns

## shared #_____________________________________________________________________________________________________1 Alignment

####____________________________________________________________________________________________________1 Match user's computer system - Alignment: MAFFT & ClustalO

    def get_mafft_path(self):
        base_path = os.path.dirname(os.path.abspath(__file__))
        if sys.platform.startswith('win'):                                      # windows
            return os.path.join(base_path, 'external_tools', 'mafft_win64', 'mafft.bat')
        elif sys.platform.startswith('darwin'):                                 # macOS
            return os.path.join(base_path, 'external_tools', 'mafft_mac', 'mafft')
        elif sys.platform.startswith('linux'):                                  # linux
            return os.path.join(base_path, 'external_tools', 'mafft_linux', 'mafft')
        else:
            raise RuntimeError('Unsupported OS')

    def get_clustalo_path(self):
        if sys.platform.startswith('win'):                                      # windows
            return os.path.join('external_tools', 'clustalo_win64', 'clustalo.exe')
        elif sys.platform.startswith('darwin'):                                 # macOS
            return os.path.join('external_tools', 'clustalo_mac', 'clustalo')
        elif sys.platform.startswith('linux'):                                  # linux
            return os.path.join('external_tools', 'clustalo_linux', 'clustalo')
        else:
            raise RuntimeError('Unsupported OS')

####____________________________________________________________________________________________________2 Run Alignment: MAFFT & ClustalO

    def run_alignment(self, sequences, group=None, layout=None, button_aln=None, output_file=None, return_only = False, seq_map=None):

        self.cancelled = False

        # 1 Check 1: Restriction on no. of sequences
        if len(sequences) < 2:
            QMessageBox.warning(self, 'Error', 'Need at least 2 sequences to align.')
            return

        # 2 Transfer sequences to file
        temp_fasta = tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix='.fasta')
        output_file = temp_fasta.name.replace('.fasta', '.aln')
        for name, seq in sequences:
            if not name or not seq:
                QMessageBox.warning(self, 'Error', 'Sequence name or sequence is missing.')
                return
            if not re.match(r'^[A-Za-z\-\_]+$', seq):
                QMessageBox.warning(self, 'Error', f'Invalid characters found sequence: {seq}')
                return
            temp_fasta.write(f'>{name}\n{seq}\n')
        temp_fasta.close()

        # 3 MAIN FUNCTION DialogBox: Progress
        ## 1 Layout
        self.widget_progress = QDialog(self)
        self.widget_progress.setModal(True)
        self.widget_progress.setWindowTitle('Processing...')                                  # block input until done
        self.layout_progress = QVBoxLayout()
        self.widget_progress.setLayout(self.layout_progress)
        ## 2 Progress Bar
        self.progress_bar = QProgressBar()
        self.progress_bar.setRange(0,0)
        self.layout_progress.addWidget(self.progress_bar)
        ## 3 Button: Cancel
        self.button_cancel = QPushButton('Cancel')
        self.button_cancel.setFixedWidth(60)
        self.button_cancel.clicked.connect(lambda: [setattr(self, 'cancelled', True), self.widget_progress.reject()])
        self.layout_progress.addWidget(self.button_cancel, alignment=Qt.AlignCenter)
        ## 4 Show
        self.widget_progress.show()
        QApplication.processEvents() 

        # 4 Run alignment: MAFFT or clustalo
        try:
            ## 1 MAFFT
            if button_aln.text() == 'MAFFT':
#                mafft_path = os.path.join('external_tools','mafft_win64', 'mafft.bat')
                mafft_path = self.get_mafft_path()
                with open(output_file, 'w') as out:
                    subprocess.run([mafft_path, '--anysymbol', '--genafpair', '--maxiterate', '10000', temp_fasta.name], check=True, stdout=out, cwd=os.path.dirname(mafft_path))
            ## 2 ClustalO
            elif button_aln.text() == 'ClustalO':
                clustalo_path = self.get_clustalo_path()
                with open(output_file, 'w') as out:
                    subprocess.run([clustalo_path, '-i', temp_fasta.name, '-o', output_file, '--dealign', '--force'], check=True, stdout=out)            
            ## 3 Neither
            else:
                QMessageBox.warning(self, 'Error', 'Error: Selecting alignment method.')
                return
        except Exception as e:
            QMessageBox.critical(self, f'{button_aln.text()} error', str(e))
        finally:
            os.remove(temp_fasta.name)

        # 4 parse the file and get the sequences
        aligned_seq = []
        for record in SeqIO.parse(output_file, 'fasta'):
            aligned_seq.append((record.id, str(record.seq))) 

        # --- Now check for cancel before starting alignment ---
        if self.cancelled:
            self.widget_progress.reject()
            return

        # 6 Clear seqs from GUI, then add sequences and consensus
        if return_only:
            return aligned_seq

        ## 1 for Alignall (all groups)
        if seq_map:
            for group in self.groups:
                layout = group['layout_seq']
                for i in reversed(range(layout.count())):
                    widget = layout.itemAt(i).widget()
                    if widget:
                        widget.setParent(None)
                        widget.deleteLater()
                group['widget_seq'].clear()
                group['consensus_seq'] = None


            # 5 Re-add aligned sequences to correct groups
            for (group_idx, _), (aligned_name, aligned_seq) in zip(seq_map, aligned_seq):

                # --- Now check for cancel before starting alignment ---
                if self.cancelled:
                    self.widget_progress.reject()
                    return

                group = self.groups[group_idx]
                layout = group['layout_seq']
                if layout is None:
                    continue
                self.add_sequences_toGUI(group, layout, aligned_name, aligned_seq)              # this is for each seq

            # get consensus for each group
            for group in self.groups:                                                           # after fill up seq, do consensus next
                # --- Now check for cancel before starting alignment ---
                if self.cancelled:
                    self.widget_progress.reject()
                    return
                self.get_consensus_aln(group, seq_map, threshold=None)



            # calculate %base conservation
            ref_consensus = None
            for group in self.groups:
                if group['checkbox_setrefgroup'].isChecked():
                    ref_consensus = group.get('consensus_seq')
                    break

            for group in self.groups:
                layout = group['layout_seq']
                target_consensus = group.get('consensus_seq')

                if not target_consensus or target_consensus == ref_consensus:
                    continue	# Skip if reference group

                matches = sum(1 for a, b in zip(ref_consensus, target_consensus) if a == b)
                percent_conservation = (matches / len(ref_consensus)) * 100
                str_percent_conservation = f"{percent_conservation:.3g}%"

                # Find QLable inside consensus_row
                for i in range(layout.count()):
                    widget = layout.itemAt(i).widget()
                    if widget and widget.objectName() == "consensus_row":
                        label = widget.findChild(QLabel, "perc_conservation")
                        if label:
                            label.setText(str_percent_conservation)

            # color code sequences
            self.color_code_seq(seq_map=self.seq_map)

            # get global consensus
            global_consensus = self.get_global_consensus()
            if global_consensus:
                # --- Now check for cancel before starting alignment ---
                if self.cancelled:
                    self.widget_progress.reject()
                    return
                with tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix='.fasta') as f:
                    f.write(f">global_consensus\n{global_consensus.replace('-', '')}\n")
                    fasta_file = f.name

                # run PSIPRED
                base_path = os.path.dirname(os.path.abspath(__file__))
                psipred_dir = os.path.join(base_path, 'external_tools', 'psipred')
                self.prediction_text = self.build_secondary_structure(fasta_file, psipred_dir, base_path)
                self.draw_secondary_structure_to_gui(self.prediction_text)
                self.compute_region_conservation(self.prediction_text)

        ## 2 for Align (by group)
        else:

            if layout is not None:
                for i in reversed(range(layout.count())):
                    # --- Now check for cancel before starting alignment ---
                    if self.cancelled:
                        self.widget_progress.reject()
                        return

                    widget = layout.itemAt(i).widget()
                    if widget:
                        widget.setParent(None)
            if group is not None:
                group['widget_seq'].clear()
                for name, seq in aligned_seq:
                    # --- Now check for cancel before starting alignment ---
                    if self.cancelled:
                        self.widget_progress.reject()
                        return

                    #group['widget_seq'].append({'seq_header': QLineEdit(name), 'seq': seq, 'seq_letters': [], 'widget_row': None})
                    self.add_sequences_toGUI(group, layout, name, seq)
                self.get_consensus_aln(group, threshold=None)


                # color code sequences
                self.color_code_seq()

                # Use consensus sequence for current group
                if 'consensus_seq' in group:
                    # --- Now check for cancel before starting alignment ---
                    if self.cancelled:
                        self.widget_progress.reject()
                        return
                    seq = group['consensus_seq'].replace('-', '')                                   # remove gaps for PSIPRED
                    with tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix='.fasta') as f:
                        f.write(f">consensus\n{seq}\n")
                        fasta_file = f.name
                    # run PSIPRED
                    base_path = os.path.dirname(os.path.abspath(__file__))
                    psipred_dir = os.path.join(base_path, 'external_tools', 'psipred')
                    self.prediction_text = self.build_secondary_structure(fasta_file, psipred_dir, base_path)
                    self.draw_secondary_structure_to_gui(self.prediction_text)
                else:
                    QMessageBox.warning(self, "Missing consensus", "Consensus sequence not found.")
                    return

        self.widget_progress.close()
        self.slider_threshold()


####____________________________________________________________________________________________________3 Get consensus sequences


    # 1 Global consensus (Alignall)
    def get_global_consensus(self, threshold=None):

        #### GUI: IF PRESENT, REMOVE CONSENSUS GLOBAL ###
        if hasattr(self, 'widget_global') and self.widget_global is not None:
            self.layout_protein_l3.removeWidget(self.widget_global)     # remove from layout (optional)
            self.widget_global.setParent(None)                          # remove from parent GUI hierarchy
            self.widget_global.deleteLater()                            # schedule for safe deletion by Qt event loop
            self.widget_global = None                                   # no widget global

        # run
        all_records = []
        for group in self.groups:

            ### to cancel half-way
            if self.cancelled:
                self.widget_progress.reject()
                return

            for entry in group['widget_seq']:
                name = entry['seq_header'].text().strip()
                seq = ''.join(label.text() for label in entry['seq_letters']).strip()
                all_records.append(SeqRecord(Seq(seq), id=name))
        if not all_records:
            return None
        alignment = MultipleSeqAlignment(all_records)
        summary = AlignInfo.SummaryInfo(alignment)
#        consensus = summary.dumb_consensus(threshold=1.0, ambiguous='X')  # strict: even 1 difference → 'N'
        if threshold is None:
            threshold = 1.0
        consensus = summary.dumb_consensus(threshold=threshold, ambiguous='X')

        # Layout: Global consensus
        ## 1 Set main layout
        self.widget_global = QWidget()
        layout_global = QHBoxLayout()
        layout_global.setContentsMargins(5,0,0,0)
        layout_global.setSpacing(0)
        self.widget_global.setLayout(layout_global)

        ## 2 Invisible checkbox
        invisible_checkbox = QCheckBox()
        invisible_checkbox.setEnabled(False)
        invisible_checkbox.setStyleSheet('background: transparent; border: none;')
        layout_global.addWidget(invisible_checkbox)
        layout_global.addSpacing(5)

        ## 3 Invisible label
        invisible_label = QLabel('Consensus')
        invisible_label.setFixedSize(120,20)
        layout_global.addWidget(invisible_label, alignment=Qt.AlignLeft)
        layout_global.addSpacing(5)

        ## 4 Set layout for sequences (by letter)
        for letter in str(consensus):
            lbl = QLabel(letter)
            lbl.setFixedSize(15,20)
            lbl.setAlignment(Qt.AlignCenter)
            lbl.setStyleSheet('color: gray;')
            layout_global.addWidget(lbl)

        ## 5 Add to main layout
        self.layout_protein_l3.addWidget(self.widget_global, alignment=Qt.AlignLeft)
        return str(consensus)                                           # for secondary structure


####____________________________________________________________________________________________________5 Consensus per group


    # 2 Consensus per group
    def get_consensus_aln(self, group, seq_map=None, threshold=None):

        # 0 If global consensus is present, remove it from GUI
#        if hasattr(self, 'widget_global') and self.widget_global is not None:
#            self.layout_protein_l3.removeWidget(self.widget_global)     # remove from layout (optional)
#            self.widget_global.setParent(None)                          # remove from parent GUI hierarchy
#            self.widget_global.deleteLater()                            # schedule for safe deletion by Qt event loop
#            self.widget_global = None                                   # no widget global

        # 1 Get sequence layout
        layout = group['layout_seq']
        if layout is None:                                              # skip layouts that are not ready yet
            return

        # 2 Clear previous consensus row
        for i in reversed(range(layout.count())):
            widget = layout.itemAt(i).widget()                          # i = last row
            if widget:
                name = widget.objectName()
                if name == "consensus_row":
                    layout.removeWidget(widget)
                    widget.deleteLater()

        # 3 Get aligned seq from group
        records = [SeqRecord(Seq(entry['seq']), id=entry['seq_header'].text()) for entry in group['widget_seq']]
        alignment = MultipleSeqAlignment(records)
        summary = AlignInfo.SummaryInfo(alignment)
#        consensus = summary.dumb_consensus(threshold=0.5, ambiguous='X')
        if threshold is None:
#            threshold = 1.0
            threshold = getattr(self, 'consensus_threshold', 1.0)
        consensus = summary.dumb_consensus(threshold=threshold, ambiguous='X')

        consensus_str = str(consensus)

        # 4 Store consensus seq for PSIRED (without gaps)
        group['consensus_seq'] = consensus_str

        # 5 Set Layout: Consensus
        widget_consensus = QWidget()
        widget_consensus.setObjectName('consensus_row')
        layout_consensus = QHBoxLayout()
        layout_consensus.setContentsMargins(5,0,0,0)
        layout_consensus.setSpacing(0)
        widget_consensus.setLayout(layout_consensus)
        layout.addWidget(widget_consensus, alignment=Qt.AlignLeft)

        # 6 Checkbox (invisible: to add spaces to align with the sequences)
        invisible_checkbox = QCheckBox()
        invisible_checkbox.setEnabled(False)
        invisible_checkbox.setStyleSheet('background: transparent; border: none;')
        layout_consensus.addWidget(invisible_checkbox)
        #layout_consensus.addSpacing(5)

        # 7 Label: Consensus
        label_consensus = QLabel('')
        label_consensus.setObjectName('perc_conservation')
        label_consensus.setFixedSize(120,20)
        layout_consensus.addWidget(label_consensus, alignment=Qt.AlignLeft)
        layout_consensus.addSpacing(5)

#        if seq_map:
#            self.consensus_calculate_conservation(label_consensus)

        # 8 Add concensus onto GUI
        for letter in consensus_str:
            lbl = QLabel(letter)
            lbl.setFixedSize(15,20)
            lbl.setAlignment(Qt.AlignCenter)
            lbl.setStyleSheet('color: gray;')
            layout_consensus.addWidget(lbl)

        return consensus_str

####____________________________________________________________________________________________________4 Add sequences from file to GUI

    def add_sequences_toGUI(self, group, layout, seq_name='', seq=''):
        # 1 set layout for sequence (horizontal)
        widget_seq = QWidget()
        layout_seq = QHBoxLayout()
        layout_seq.setContentsMargins(5,0,0,0)
        layout_seq.setSpacing(0)
        widget_seq.setLayout(layout_seq)
        layout.addWidget(widget_seq, alignment=Qt.AlignLeft)

        # 2 +Checkbox to sequence layout
        seq_checkbox = QCheckBox()
        layout_seq.addWidget(seq_checkbox)
        layout_seq.addSpacing(5)
        self.seq_rows.append((widget_seq, seq_checkbox))        # save the rows and checkbox

        # 3 +Seq header to sequence layout
        seq_header = QLineEdit(seq_name)
        seq_header.setFixedSize(120,20)
        layout_seq.addWidget(seq_header)
        layout_seq.addSpacing(5)
        layout_seq.setAlignment(Qt.AlignLeft)

        seq_letters = []
        for letter in seq:
            ### MAIN FUNCTION: to cancel half-way
            if self.cancelled:
                self.widget_progress.reject()
                return

            ## to add letters onto GUI
            seq_letter = QLabel(letter)
            seq_letter.setFixedSize(15,20)
            seq_letter.setAlignment(Qt.AlignCenter)
            seq_letter.setStyleSheet('border: 1px solid #F8F8F8; padding: 1px;')
            layout_seq.addWidget(seq_letter)
            seq_letters.append(seq_letter)

        #-- 5 Create dictionary to store objects related to sequences
        group['widget_seq'].append({'seq_header': seq_header, 'seq': seq, 'seq_letters': seq_letters, 'widget_row': widget_seq})

        # 6 Horizontal scrollbar: Dynamically resize parent container width based on longest sequence
        avg_char_width = 8                                                  # Approximate pixel width per character (monospace)
        min_padding = 200                                                   # Extra space for name + spacing
        seq_pixel_length = len(seq) * avg_char_width + min_padding          # Box length applied
        current_width = self.widget_protein_l4.minimumWidth()               # Get the width of parent layout
        if seq_pixel_length > current_width:                                # Apply it to self.widget_protein_l4 if wider than current
            self.widget_protein_l4.setMinimumWidth(seq_pixel_length)        # Set the length




#_______________________________________________________________________________________________
#_______________________________________________________________________________________________4 Secondary Structure
#_______________________________________________________________________________________________




    def build_secondary_structure(self, fasta_file, psipred_dir, base_path):

        # 0 If global consensus is present, remove it from GUI
        if hasattr(self, 'widget_horizontal') and self.widget_horizontal is not None:
            self.layout_protein_l4_2ndarystructure.removeWidget(self.widget_horizontal)     # remove from layout (optional)
            self.widget_horizontal.setParent(None)                          # remove from parent GUI hierarchy
            self.widget_horizontal.deleteLater()                            # schedule for safe deletion by Qt event loop
            self.widget_horizontal = None                                   # no widget global

        base = os.path.splitext(os.path.basename(fasta_file))[0]                                        # get the file name
        out_folder = os.path.join(base_path, 'tmp_files')
        os.makedirs(out_folder, exist_ok=True)
        horiz_file = f"{out_folder}/{base}.horiz"

        # 1 Run PSIPRED with PSI-BLAST
        try:
            ### MAIN FUNCTION: to cancel half-way
            if self.cancelled:
                self.widget_progress.reject()
                return

            ## run
            runpsipred = os.path.join(psipred_dir, "BLAST+", "runpsipredplus")
            blastdb_path = os.path.join(psipred_dir, "BLAST+", "blastdb")
            env = os.environ.copy()
            env["BLASTDB"] = blastdb_path
            subprocess.run([runpsipred, fasta_file], check=True, env=env, cwd=out_folder)
            print("Run runpsipredplus + PSI-BLAST")

        # 2 Run PSIPRED without BLAST
        except subprocess.CalledProcessError as e:
            print("runpsipredplus failed, trying runpsipred_single instead...")
            try:
                ### MAIN FUNCTION: to cancel half-way
                if self.cancelled:
                    self.widget_progress.reject()
                    return
                ## run
                runpsipred_single = os.path.join(psipred_dir, "runpsipred_single")
                subprocess.run([runpsipred_single, fasta_file], check=True, cwd=out_folder)
                print("Run runpsipred_single")
            except subprocess.CalledProcessError as e:
                print("Both runpsipred and runpsipred_single failed.")
                raise e

        # Make sure .horiz file exists
        if not os.path.exists(horiz_file):
            raise FileNotFoundError(f"Expected output file {horiz_file} not found!")
        with open(horiz_file, "r") as f:
            self.prediction_text = f.read()
        return self.prediction_text


##
    def draw_secondary_structure_to_gui(self, prediction_text):
        # Step 1: Parse structure and sequence
        aa, pred = '', ''
        for line in prediction_text.splitlines():
            if 'AA:' in line:
                parts = line.strip().split()
                if len(parts) > 1:
                    aa += parts[1]                                                  # extract the seq after 'AA:' (AA Seq)
            elif 'Pred:' in line:
                parts = line.strip().split()
                if len(parts) > 1:
                    pred += parts[1]                                                # extract the seq after 'Pred:' (Conformation)

        # Step 2: Translate prediction to symbolic characters
        width_per_residue = 0.15
        fig_width = len(aa) * width_per_residue
        fig, ax = plt.subplots(figsize=(fig_width, 0.4), dpi=100)

        i=0
        while i < len(pred):
            ss = pred[i]
            start = i
            while i < len(pred) and pred[i] == ss:
                i += 1
            end = i

            if ss == 'H':
                rect = Rectangle((start, 0.1), end - start, 0.8, linewidth=1, edgecolor='red', facecolor='red', alpha=0.4)
                ax.add_patch(rect)
            elif ss == 'E':
                arrow = FancyArrow(start, 0.5, end - start - 0.2, 0, width=0.3, length_includes_head=True, head_width=0.5, head_length=0.3, color='blue')
                ax.add_patch(arrow)
            else:
                ax.plot([start, end], [0.5, 0.5], color='gray', linewidth=1.2)

        # Step 3: Create figure
#        fig, ax = plt.subplots(figsize=(len(aa) * 0.15, 0.4))  # 0.15 in ~ 15 pixels
        ax.set_xlim(0, len(aa))
        ax.set_ylim(0, 1)
        ax.axis('off')
        plt.tight_layout(pad=0)

        # Step 5: Save to buffer (tight bounding box)
        buffer = BytesIO()
        fig.savefig(buffer, format='png', transparent=True, bbox_inches='tight', pad_inches=0)
        plt.close(fig)

        # Step 4: Convert to QPixmap and display in QLabel
        pixmap = QPixmap()
        pixmap.loadFromData(buffer.getvalue())
        label = QLabel()
        label.setPixmap(pixmap)

        # 5 Build layout aligned with seuqence rows
        self.widget_horizontal = QWidget()
        layout_horizontal = QHBoxLayout()
        layout_horizontal.setContentsMargins(0,0,0,0)
        layout_horizontal.addSpacing(0)
        self.widget_horizontal.setLayout(layout_horizontal)

        ## 6 Checkbox (invisible: to add spaces to align with the sequences)
        invisible_checkbox = QCheckBox()
        invisible_checkbox.setEnabled(False)
        invisible_checkbox.setStyleSheet('background: transparent; border: none;')
        layout_horizontal.addWidget(invisible_checkbox)
        #layout_horizontal.addSpacing(5)

        ## 2 Label: Consensus
        invisible_label = QLabel('')
        invisible_label.setFixedSize(118,20)
        layout_horizontal.addWidget(invisible_label, alignment=Qt.AlignLeft)
        #layout_horizontal.addSpacing(5)

        ## 3 Label: Secondary structure
        layout_horizontal.addWidget(label, alignment=Qt.AlignLeft)

        ## 4 Add Widget to main widget
        self.layout_protein_l4_2ndarystructure.addWidget(self.widget_horizontal, alignment=Qt.AlignLeft)

        # 00 MAIN FUNCTION Progress
        label_progress_sec = QLabel('Secondary structure is generated.')
        self.layout_progress.addWidget(label_progress_sec)
        QApplication.processEvents()


    def compute_region_conservation(self, prediction_text):

        # 1: Get ref_consensus & target_consensus
        ref_cons = None
        target_cons = None
        for group in self.groups:
            if group['checkbox_setrefgroup'].isChecked():
                ref_cons = group.get('consensus_seq')
                break
        if not ref_cons:
            return

        # Step 2: Parse structure and sequence
        aa, pred = '', ''
        for line in prediction_text.splitlines():
            if 'AA:' in line:
                parts = line.strip().split()
                if len(parts) > 1:
                    aa += parts[1]                                                  # extract the seq after 'AA:' (AA Seq)
            elif 'Pred:' in line:
                parts = line.strip().split()
                if len(parts) > 1:
                    pred += parts[1]                                                # extract the seq after 'Pred:' (Conformation)

        regions = []
        current_type = None
        start = None
        for i, ss in enumerate(pred):
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
            regions.append((current_type, start, len(pred) - 1))

        for group in self.groups:
            if group['checkbox_setrefgroup'].isChecked():
                continue
            layout = group['layout_seq']
            for i in reversed(range(layout.count())):
                widget = layout.itemAt(i).widget()
                if widget and widget.objectName() == "conservation_block":
                    layout.removeWidget(widget)
                    widget.deleteLater()

            target_cons = group.get('consensus_seq')
            if not target_cons:
                continue

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

            # Pre-fill all columns with empty labels
            total_cols = len(ref_cons)
            labels = [QLabel('') for _ in range(total_cols)]
            for lbl in labels:
                lbl.setFixedSize(15, 20)
                lbl.setAlignment(Qt.AlignCenter)
                layout_result.addWidget(lbl)

            # Step 3: Fill labels only at the midpoint of each region
            for ss_type, start, end in regions:
                target_cons = group.get('consensus_seq')
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

            widget_result.setLayout(layout_result)
            main_widget_result = QWidget()
            main_layout_result = QHBoxLayout()
            main_layout_result.setContentsMargins(0,0,0,0)
            main_layout_result.setSpacing(0)
            main_widget_result.setLayout(main_layout_result)
            main_layout_result.addWidget(widget_result, alignment=Qt.AlignLeft)
            layout.addWidget(main_widget_result)




#_______________________________________________________________________________________________
#_______________________________________________________________________________________________5 Color code alignment
#_______________________________________________________________________________________________




    def color_code_seq(self, seq_map=None):
        def get_col_similarity(seqs):
            transposed = list(zip(*seqs))
            similarity = []
            for col in transposed:
                most_common = max(set(col), key=col.count)
                similarity_score = col.count(most_common) / len(col)
                similarity.append(similarity_score)
            return similarity

        def similarity_to_color(score):
            return mcolors.to_hex(mcolors.LinearSegmentedColormap.from_list('custom', ['white', '#5b005b'])(score))

        # 1 Alignall fxn
        if seq_map:
            all_seqs = []
            for group in self.groups:
                for entry in group['widget_seq']:
                    all_seqs.append(entry['seq'])
    
            similarity = get_col_similarity(all_seqs)

            for idx_col, score in enumerate(similarity):
                for (group_idx, entry_idx) in seq_map:
                    entry = self.groups[group_idx]['widget_seq'][entry_idx]
                    if idx_col < len(entry['seq_letters']):
                        lbl = entry['seq_letters'][idx_col]
                        lbl.setProperty("bg_color", similarity_to_color(score))
                        lbl.setStyleSheet(f'background-color: {similarity_to_color(score)};')

#        else:
#            for group in self.groups:
#                if 'widget_seq' not in group or not group['widget_seq']:
#                    continue
#                seqs = [entry['seq'] for entry in group['widget_seq']]
#                similarity = get_col_similarity(seqs)
#                for entry in group['widget_seq']:
#                    seq_len = len(entry['seq_letters'])
#                    for idx in range(seq_len):
#                        if idx < len(similarity):
#                            lbl = entry['seq_letters'][idx]
#                            lbl.setStyleSheet(f'background-color: {similarity_to_color(similarity[idx])};')




#_______________________________________________________________________________________________
#_______________________________________________________________________________________________5 Slider
#_______________________________________________________________________________________________




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
                        bg_color = lbl.property("bg_color") or 'white'    #self.get_existing_bg_color(lbl)
                        if lower <= percent <= upper:
                            lbl.setStyleSheet(f'background-color: {bg_color}; color: black;')
                        else:
                            lbl.setStyleSheet('color: lightgray;')

        else:
            self.color_code_seq(seq_map=self.seq_map)




#_______________________________________________________________________________________________
#_______________________________________________________________________________________________6 Custom display % conservation
#_______________________________________________________________________________________________





    def custom_display_perc_cons(self, pos1, pos2):
       # 1. Get reference consensus
        ref_consensus = None
        for group in self.groups:

            if group['checkbox_setrefgroup'].isChecked():
                ref_consensus = group.get('consensus_seq')
                break

        if not ref_consensus:
            print("No reference group selected.")
            return

        # Ensure positions are sorted
        start = min(pos1, pos2)
        end = max(pos1, pos2)

        # 2. For each group, compare consensus in the selected region
        for group in self.groups:
            if group['checkbox_setrefgroup'].isChecked():
                continue  # skip reference group

            target_consensus = group.get('consensus_seq')
            if not target_consensus:
                continue

            # Clean up existing blocks
            layout = group['layout_seq']
            for i in reversed(range(layout.count())):
                widget = layout.itemAt(i).widget()
                if widget and widget.objectName() == "custom_conservation_block":
                    layout.removeWidget(widget)
                    widget.deleteLater()

            # 3. Calculate percent similarity in selected region
            match_count = sum(
                1 for i in range(start, end + 1)
                if i < len(ref_consensus) and i < len(target_consensus) and ref_consensus[i] == target_consensus[i]
            )
            total = end - start + 1
            percent = (match_count / total * 100) if total else 0

            # 4. Create a row of labels with % displayed at center
            mid = (start + end) // 2

            widget_result_main = QWidget()
            widget_result_main.setObjectName("custom_conservation_block")
            layout_result_main = QHBoxLayout()
            widget_result_main.setLayout(layout_result_main)
            layout_result_main.setContentsMargins(0, 0, 0, 0)
            layout_result_main.setSpacing(0)

            widget_result = QWidget()
            layout_result = QHBoxLayout()
            layout_result.setContentsMargins(0, 0, 0, 0)
            layout_result.setSpacing(0)
            layout_result_main.addWidget(widget_result, alignment=Qt.AlignLeft)

            # Add spacing label (to align left indent)
            invisible_checkbox = QCheckBox()
            invisible_checkbox.setEnabled(False)
            invisible_checkbox.setStyleSheet('background: transparent; border: none;')
            layout_result.addWidget(invisible_checkbox)

            invisible_label = QLabel('Pos:'+ str(start) + ' - ' + str(end))
            invisible_label.setStyleSheet("color: gray; font-size: 9px;")
            layout_result.addWidget(invisible_label)

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






#_______________________________________________________________________________________________
#_______________________________________________________________________________________________6 Save Project
#_______________________________________________________________________________________________



    def save_project(self):
        file_path, _ = QFileDialog.getSaveFileName(self, "Save Project", "", "Alignate Project (*.alignate)")
        if not file_path:
            return

        # Automatically append extension if missing
        if not file_path.endswith(".alignate"):
            file_path += ".alignate"

        state = {
            "groups": [],
            "slider_value": self.slidercon.value(),
            "slider_checked": self.checkboxslider.isChecked()
        }

        for group in self.groups:
            group_data = {
                "name": group['lineedit_groupname'].text(),
                "is_reference": group['checkbox_setrefgroup'].isChecked(),
                "sequences": []
            }
            for entry in group['widget_seq']:
                header = entry['seq_header'].text()
                seq = ''.join(label.text() for label in entry['seq_letters'])
                group_data['sequences'].append((header, seq))
            state["groups"].append(group_data)

        with tempfile.TemporaryDirectory() as temp_dir:
            json_path = os.path.join(temp_dir, "state.json")
            with open(json_path, "w") as f:
                json.dump(state, f)

#            canvas_path = os.path.join(temp_dir, "canvas.png")
#            self.canvas.image.save(canvas_path)

            with zipfile.ZipFile(file_path, "w") as zf:
                zf.write(json_path, "state.json")
#                zf.write(canvas_path, "canvas.png")

        QMessageBox.information(self, "Saved", f"Project saved to {file_path}")



#_______________________________________________________________________________________________
#_______________________________________________________________________________________________6 Load Project
#_______________________________________________________________________________________________




    def load_project(self):
        file_path, _ = QFileDialog.getOpenFileName(self, "Load Project", "", "Alignate Project (*.alignate)")
        if not file_path:
            return

        with tempfile.TemporaryDirectory() as temp_dir:
            with zipfile.ZipFile(file_path, "r") as zf:
                zf.extractall(temp_dir)

            with open(os.path.join(temp_dir, "state.json")) as f:
                state = json.load(f)

            self.slidercon.setValue(state.get("slider_value", 100))
            self.checkboxslider.setChecked(state.get("slider_checked", False))

            # Clear current state
            for group in self.groups[:]:
                self.button2_removegroup_clicked(group['widget_group'])

            for group_data in state["groups"]:
                self.button1_addgroup_clicked()
                group = self.groups[-1]
                group['lineedit_groupname'].setText(group_data["name"])
                group['checkbox_setrefgroup'].setChecked(group_data["is_reference"])
                for header, seq in group_data["sequences"]:
                    self.add_sequences_toGUI(group, group["layout_seq"], header, seq)

            # Load canvas
#            self.canvas.image.load(os.path.join(temp_dir, "canvas.png"))
#            self.canvas.update()

        QMessageBox.information(self, "Loaded", f"Project loaded from {file_path}")





#_______________________________________________________________________________________________
#_______________________________________________________________________________________________6 General
#_______________________________________________________________________________________________





#_______________________________________________________________________________________________
#_______________________________________________________________________________________________6 Action: File
#_______________________________________________________________________________________________




#_______________________________________________________________________________________________
#_______________________________________________________________________________________________7 Action: View
#_______________________________________________________________________________________________




#_______________________________________________________________________________________________
#_______________________________________________________________________________________________8 Action: Help
#_______________________________________________________________________________________________




#_______________________________________________________________________________________________
#_______________________________________________________________________________________________9 Main Execution
#_______________________________________________________________________________________________




window = main()
window.show()
app.exec()


# To give thought
## 1 To limit the seq length to <1k
