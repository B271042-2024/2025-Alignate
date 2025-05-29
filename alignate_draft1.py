import sys, re, subprocess, tempfile, os
from PySide6.QtGui import QIcon, QPixmap, QFont
from PySide6.QtWidgets import QStackedWidget, QProgressBar, QFileDialog, QMessageBox, QDialog, QTextEdit, QDialogButtonBox, QLayout, QScrollArea, QSizePolicy, QApplication, QMainWindow, QWidget, QCheckBox, QLabel, QLineEdit, QPushButton, QGroupBox, QToolBar, QHBoxLayout, QVBoxLayout, QScrollBar, QSlider
from PySide6.QtCore import QSize, Qt
from Bio import SeqIO
app = QApplication(sys.argv)
app.setStyle("Fusion")
icon_logo = 'logo_ninja1'


#_______________________________________________________________________________________________1 Main Window

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
        menu2.addAction('All')
        menu2.addAction('Selection')
        menu2.addAction('Unselected')
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


    #1 def Toolbar_action: window about
    def trigger_toolbar1_about(self):
        self.setCentralWidget(self.window_about)


    #2 def Toolbar_action: window protein
    def trigger_toolbar2_protein(self):
        self.setCentralWidget(self.window_protein)


    



#_______________________________________________________________________________________________2 Toolbar: About

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



#_______________________________________________________________________________________________3 Toolbar: Protein

class protein(QWidget):
    def __init__(self):
        super().__init__()      

#_______1 SET PROTEIN MAIN LAYOUT___1 Set main layout

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

        # 4 Set the widget 4 to put in the widget protein layer 3
        self.widget_protein_l4 = QWidget()
        self.layout_protein_l4 = QVBoxLayout()
        self.layout_protein_l4.setSpacing(0)
        self.layout_protein_l4.setContentsMargins(0, 0, 0, 0)
        self.widget_protein_l4.setLayout(self.layout_protein_l4)
        self.layout_protein_l3.addWidget(self.widget_protein_l4, alignment=Qt.AlignTop)

        # 5 Checkbox
        self.seq_rows = []      # store all checked rows in a tuple (used to remove sequences)

        #-- 6 Store objects (for alignments)
        self.groups = []



#_______1 SET PROTEIN MAIN LAYOUT___2 Set main elements

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

        # 3 Setup layout for button1 & button2
        self.widget_protein_buttons = QWidget()
        self.layout_protein_buttons = QHBoxLayout()
        self.layout_protein_buttons.setContentsMargins(0,2,0,2)
        self.widget_protein_buttons.setLayout(self.layout_protein_buttons)
        self.layout_protein_buttons.addWidget(self.button1_addgroup, alignment=Qt.AlignLeft)
        self.layout_protein_buttons.addWidget(self.button2_alignall, alignment=Qt.AlignLeft)
        self.layout_protein_buttons.addStretch(1)

        # 4 Setup widget for secondary alignment
        self.widget_protein_l4_2ndarystructure = QWidget()
        self.layout_protein_l4_2ndarystructure = QHBoxLayout()
        self.layout_protein_l4_2ndarystructure.setContentsMargins(0,2,0,2)
        self.widget_protein_l4_2ndarystructure.setLayout(self.layout_protein_l4_2ndarystructure)
        self.widget_protein_l4_2ndarystructure.setObjectName("line_0")
        self.widget_protein_l4_2ndarystructure.setStyleSheet(
            """ #line_0 {
                border: 1px solid #A9A9A9;
                border-radius: 6px;
                padding: 2px;
            } """
        )

        # 5 Add elements to the layout
        self.layout_protein_l4.addWidget(self.widget_protein_buttons, alignment=Qt.AlignTop)
        self.layout_protein_l4.addWidget(self.widget_protein_l4_2ndarystructure)

#_______1 SET PROTEIN MAIN LAYOUT___3 Connect to DEF

        self.button1_addgroup.clicked.connect(self.button1_addgroup_clicked)
        self.button2_alignall.clicked.connect(self.button2_alignall_clicked)


#_______2 SET PROTEIN MAIN DEF___DEF 1: Add New Group
    def button1_addgroup_clicked(self):

#_______2 SET PROTEIN MAIN DEF___DEF 1: Add New Group___1 Set main layout
        # 1 Set main layout
        self.widget_protein_l4_group_l1 = QWidget()                                                  # widget_protein_l3_group_layer1
        self.widget_protein_l4_group_l1.setObjectName("group")
        self.widget_protein_l4_group_l1.setStyleSheet(
            """ #group {
                border: 1px solid #A9A9A9;
                border-radius: 6px;
                padding: 2px;
            } """
        )
        self.layout_protein_l4_group_l1 = QVBoxLayout()
        self.layout_protein_l4_group_l1.setContentsMargins(0,2,0,2)                                               
        self.widget_protein_l4_group_l1.setLayout(self.layout_protein_l4_group_l1)                   # layout_protein_l3_group_layer1 for widget layer1
        self.layout_protein_l4.addWidget(self.widget_protein_l4_group_l1)                            # add this widget to the main_widget_layout


#_______2 SET PROTEIN MAIN DEF___DEF 1: Add New Group___2 Set LINE 1
        # 1 Set layout (for LINE 1 - buttons)
        widget_protein_l4_group_l1_line1 = QWidget()                                            # widget for line 1
        layout_protein_l4_group_l1_line1 = QHBoxLayout()
        widget_protein_l4_group_l1_line1.setLayout(layout_protein_l4_group_l1_line1)
        self.layout_protein_l4_group_l1.addWidget(widget_protein_l4_group_l1_line1, alignment=Qt.AlignTop)

        # 1 Set layout (for LINE 2 - sequences)
        widget_protein_l4_group_l1_seq = QWidget()
        layout_protein_l4_group_l1_seq = QVBoxLayout()
        widget_protein_l4_group_l1_seq.setLayout(layout_protein_l4_group_l1_seq)
        layout_protein_l4_group_l1_seq.setSpacing(0)
        layout_protein_l4_group_l1_seq.setContentsMargins(5,0,0,0)
        self.layout_protein_l4_group_l1.addWidget(widget_protein_l4_group_l1_seq)

        # 2 Set elements
        lineedit_groupname = QLineEdit('edit group name')
        lineedit_groupname.setFixedWidth(120)
        button2_removegroup = QPushButton('-Group')
        button2_removegroup.setFixedWidth(60)
        button3_addseq = QPushButton('+')
        button3_addseq.setFixedWidth(30)
        button4_removeseq = QPushButton('-')
        button4_removeseq.setFixedWidth(30)
        button5_align = QPushButton('Align')
        button5_align.setFixedWidth(60)

        button5_align.setStyleSheet(
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


        # 3 Add elements to the layout
        layout_protein_l4_group_l1_line1.addWidget(lineedit_groupname)
        layout_protein_l4_group_l1_line1.addWidget(button2_removegroup)
        layout_protein_l4_group_l1_line1.addWidget(button3_addseq)
        layout_protein_l4_group_l1_line1.addWidget(button4_removeseq)
        layout_protein_l4_group_l1_line1.addWidget(button5_align)
        layout_protein_l4_group_l1_line1.addStretch(1)                      # to left-align all elements

        #-- 4 Create dictionary to store the layout (for alignment)
        group = {
            'widget_group': self.widget_protein_l4_group_l1,
            'lineedit_groupname': lineedit_groupname,
            'layout_seq': layout_protein_l4_group_l1_seq,
            'widget_seq': []
        }
        self.groups.append(group)


        # 5  Connect to DEF (lambda: to Specify Target)
        button2_removegroup.clicked.connect(lambda _, w=self.widget_protein_l4_group_l1: self.button2_removegroup_clicked(w))
        button3_addseq.clicked.connect(lambda _, layout=layout_protein_l4_group_l1_seq: self.button3_addseq_clicked(layout))
        button4_removeseq.clicked.connect(lambda _, layout=self.layout_protein_l4_group_l1: self.button4_removeseq_clicked(layout))
#        button5_align.clicked.connect(lambda _, layout=layout_protein_l4_group_l1_seq: self.button5_align_clicked(layout))

        layout = layout_protein_l4_group_l1_seq
        layout_line1 = layout_protein_l4_group_l1_line1
        button5_align.clicked.connect(lambda _, layout=layout, line1_layout=layout_line1: self.button5_align_clicked(layout, line1_layout))


#_______2 SET PROTEIN MAIN DEF___DEF 1: Add New Group___2 Set LINE 1__DEF 1: Remove Group
    def button2_removegroup_clicked(self, widget):
        widget.setParent(None)
        widget.deleteLater()
        

#_______2 SET PROTEIN MAIN DEF___DEF 1: Add New Group___2 Set LINE 1__DEF 2: Add Sequences
#_______2 SET PROTEIN MAIN DEF___DEF 1: Add New Group___2 Set LINE 1__DEF 2: Add Sequences___1 Method: Input text/Upload file
    def button3_addseq_clicked(self, layout):
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
        self.seqtext_button1text.clicked.connect(lambda _, layout= layout: self.seqtext_button1text_clicked(layout))
        self.seqtext_button2file.clicked.connect(lambda _, layout=layout: self.seqtext_button2file_clicked(layout))

        # 2 Dialog window: execute
        self.widget_protein_l4_group_l1_seq_dialoginput.exec()

#_______2 SET PROTEIN MAIN DEF___DEF 1: Add New Group___2 Set LINE 1__DEF 2: Add Sequences___1 Method: Input text___Dialog box
    ##2 input sequences from a dialog box
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
        widget_seq_text_inputtext_dialogbox_button1ok.clicked.connect(lambda _, layout=layout: self.widget_seq_text_inputtext_dialogbox_button1_clicked(layout))

        ###5 Dialog window: execute
        self.widget_seq_text_inputtext_dialogbox.exec()

#_______2 SET PROTEIN MAIN DEF___DEF 1: Add New Group___2 Set LINE 1__DEF 2: Add Sequences___1 Method: Input text___Extract Sequences
    def widget_seq_text_inputtext_dialogbox_button1_clicked(self,layout):       # layout: widget for sequences

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

        # 4 Connect to def: add_sequences to GUI
        for seq_name, seq in zip(self.allinput_header, self.allinput_seq):
#            self.add_sequences_toGUI(layout, seq_name, seq)
            for group in self.groups:
                if group['layout_seq'] == layout:
                    self.add_sequences_toGUI(group, layout, seq_name, seq)
                    break


#_______2 SET PROTEIN MAIN DEF___DEF 1: Add New Group___2 Set LINE 1__DEF 2: Add Sequences___2 Method: Upload file
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
#                self.add_sequences_toGUI(layout, name, seq)                                     # Transfer to function GUI
                for group in self.groups:
                    if group['layout_seq'] == layout:
                        self.add_sequences_toGUI(group, layout, name, seq)


        except Exception as e:
            QMessageBox.critical(self, 'Error', f'Could not load file:\n{e}')

        # 3 Hide windows
        self.widget_protein_l4_group_l1_seq_dialoginput.hide()


#_______2 SET PROTEIN MAIN DEF___DEF 1: Add New Group___2 Set LINE 2: Add Seq to GUI___DEF 1: Add Sequences to GUI
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



#_______2 SET PROTEIN MAIN DEF___DEF 1: Add New Group___2 Set LINE 1__DEF 3: Remove Sequences
    def button4_removeseq_clicked(self, layout):
        seq_removed = []
        for row, checkbox in self.seq_rows:                     # keep the rows in a new list
            if checkbox.isChecked():
                seq_removed.append((row, checkbox))

        for row, checkbox in seq_removed:                       # delete them
            row.setParent(None)
            row.deleteLater()
 
        self.seq_rows = [item for item in self.seq_rows if item not in seq_removed]



#_______2 SET PROTEIN MAIN DEF___DEF 1: Add New Group___1 Set LINE 1__DEF 3: Align
    def button5_align_clicked(self, layout, layout_line1):

#___________________# 1 Dialog Box
        # 1 Set layout
        widget_dialogbox_align = QDialog()
        layout_dialogbox_align = QVBoxLayout()
        widget_dialogbox_align.setLayout(layout_dialogbox_align)

        # 2 Set elements: 1 buttons (method: mafft, clustalo)
        widget_buttons_align = QDialogButtonBox()
        layout_buttons_align = QHBoxLayout()
        widget_buttons_align.setLayout(layout_buttons_align)
        widget_button1_clustalo = QPushButton('ClustalO')
        widget_button2_mafft = QPushButton('MAFFT')
        widget_button3_cancel = QPushButton('Cancel')
        widget_buttons_align.addButton(widget_button1_clustalo, QDialogButtonBox.ActionRole)
        widget_buttons_align.addButton(widget_button2_mafft, QDialogButtonBox.ActionRole)
        widget_buttons_align.addButton(widget_button3_cancel, QDialogButtonBox.RejectRole)
        layout_dialogbox_align.addWidget(widget_buttons_align)

#___________________# 2 Extract the correct group using the layout
        group = None
        for g in self.groups:                                                               # Only extract the correct group
            if g['layout_seq'] == layout:
                group = g
                break

#___________________# 3 Connect to DEF
        widget_button3_cancel.clicked.connect(widget_dialogbox_align.reject)                # cancel
        widget_button2_mafft.clicked.connect(lambda: (
            widget_dialogbox_align.accept(),
            #self.run_mafft([(entry['seq_header'].text().strip(), entry['seq'].toPlainText().strip()) for entry in group['widget_seq']], group, layout, layout_line1)))
            self.run_mafft([(entry['seq_header'].text().strip(), ''.join(label.text() for label in entry['seq_letters']).strip()) for entry in group['widget_seq']], group, layout)))


#___________________# 4 Add Dialog Box to GUI
        widget_dialogbox_align.exec()


#___________________# 5 Parse .aln file
    def parse_aln(self, output_file):
        aligned = []
        for record in SeqIO.parse(output_file, 'fasta'):
            aligned.append((record.id, str(record.seq)))
        return aligned


#___________________# 6 Run MAAFT and add to GUI (by group)
    def run_mafft(self, sequences, group=None, layout=None, output_file=None, return_only = False):

        print(sequences)

        # 1 Check 1: Restriction on no. of sequences
        if len(sequences) < 2:
            QMessageBox.warning(self, 'Error', 'You need to have at least 2 sequences to align.')
            return
        
        # 2 Set temporary output file
        temp_fasta = tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix='.fasta')
        output_file = temp_fasta.name.replace('.fasta', '.aln')

        # 3 Write sequences to fasta
        for name, seq in sequences:
        ## 1 Check 2: if missing sequences & if sequence has special characters
            if not name or not seq:
                QMessageBox.warning(self, 'Error', 'Sequence name or sequence is missing.')
                return
            if not re.match(r'^[A-Za-z\-\_]+$', seq):
                QMessageBox.warning(self, 'Error', f'Invalid characters found sequence: {seq}')
                return
        ## 2 write to file    
            temp_fasta.write(f'>{name}\n{seq}\n')
        temp_fasta.close()

        # 4 Run MAFFT
        mafft_path = os.path.join('external_tools','mafft', 'mafft.bat')
        try:
            with open(output_file, 'w') as out:
                subprocess.run([mafft_path, '--anysymbol', '--genafpair', '--maxiterate', '10000', temp_fasta.name], check=True, stdout=out)
        except Exception as e:
            QMessageBox.critical(self, 'MAFFT error', str(e))
        finally:
            os.remove(temp_fasta.name)

        # 5 parse the file and get the sequences
        aligned_seq = self.parse_aln(output_file)

        if return_only:
            return aligned_seq

        # 6 clear seq in GUI
        for i in reversed(range(layout.count())):
            widget = layout.itemAt(i).widget()
            if widget:
                widget.setParent(None)
        group['widget_seq'].clear()

        # 7 add aligned seq
        for name, seq in aligned_seq:
            self.add_sequences_toGUI(group, layout, name, seq)



#_______2 SET PROTEIN MAIN DEF___DEF 2: Align all sequences
    def button2_alignall_clicked(self):
        all_seq = []
        seq_map = []

        # 1 Collect all sequences and mapping info
        for group_idx, group in enumerate(self.groups):
            for entry_idx, entry in enumerate(group['widget_seq']):
                name = entry['seq_header'].text().strip()
                seq = ''.join(label.text() for label in entry['seq_letters']).strip()

                if name and seq:
                    all_seq.append((name, seq))
                    seq_map.append((group_idx, entry_idx))
                else:
                    QMessageBox.warning(self, 'Error', 'Missing sequence header or sequences.')
                    return
                
        if len(all_seq) < 2:
            QMessageBox.warning(self, 'Error', 'You need at least 2 sequences to align.')
            return
        
        # 2 Run MAFFT as one batch
        try:
            aligned = self.run_mafft(all_seq, return_only=True)
        except Exception as e:
            QMessageBox.critical(self, 'MAFFT error', str(e))
            return

        # 3 Clear all groups' GUI and widget_seq
        for group in self.groups:
            layout = group['layout_seq']
            for i in reversed(range(layout.count())):
                widget = layout.itemAt(i).widget()
                if widget:
                    widget.setParent(None)
            group['widget_seq'].clear()

        # 4 Re-add aligned sequences to correct groups
        for (group_idx, _), (aligned_name, aligned_seq) in zip(seq_map, aligned):
            group = self.groups[group_idx]
            layout = group['layout_seq']

            if layout is None:
                continue

            self.add_sequences_toGUI(group, layout, aligned_name, aligned_seq)







#7 Action: File

#8 Action: View

#9 Action: Help





    #3 def Toolbar_action: window codon
        

# Execute
window = main()
window.show()
app.exec()





# To give thought
## 1 To limit the seq length to <1k


