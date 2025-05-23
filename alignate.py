import sys
from PySide6.QtGui import QIcon, QPixmap
from PySide6.QtWidgets import QApplication, QMainWindow, QWidget, QScrollArea, QGridLayout, QLabel, QLineEdit, QPushButton, QGroupBox, QToolBar, QHBoxLayout, QVBoxLayout, QScrollBar, QSlider
from PySide6.QtCore import QSize, Qt

app = QApplication(sys.argv)
app.setStyle("Fusion")

icon_logo = 'alignate_logo'

# Sections & Toolbar
#1 Section: Alignment box

#2 Section: Secondary structure box

#3 Section: %Conservation

#4 Toolbar: About
class about(QWidget):
    def __init__(self):
        super().__init__()

        #1 setup text
        text1 = "Alignate v1.0 was created on 22/05/2025 by Adrian.\nnothing is impossible..."
        label_text1 = QLabel(text1)
        label_text1.setAlignment(Qt.AlignLeft | Qt.AlignTop)

        # setup layout
        v_layout = QVBoxLayout()
        v_layout.addWidget(label_text1)
        self.setLayout(v_layout)

#5 Toolbar: Protein & Codon

# Actions
#4 Action: Add new sequences (multiple selections)

#5 Action: Align 1 group    (protein & codon)

#6 Action: Align all groups (protein & codon)

#7 Action: File

#8 Action: View

#9 Action: Help

# Class: main
class main(QMainWindow):
    def __init__(self):
        super().__init__()

        #1 setup basic interface
        self.setWindowTitle('Alignate')
        self.setWindowIcon(QIcon(icon_logo))
        self.setMinimumSize(800,500)

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

        toolbar = QToolBar()
        toolbar1_about = toolbar.addAction('About')
        toolbar1_about.triggered.connect(self.trigger_toolbar1_about)
        toolbar2_protein = toolbar.addAction('Protein')
        toolbar2_protein.triggered.connect(self.trigger_toolbar2_protein)
        toolbar3_codon = toolbar.addAction('Codon')
#        toolbar3_codon.triggered.connect(self.trigger_toolbar3_codon)
        self.addToolBar(toolbar)


    #1 def Toolbar_action: window about
    def trigger_toolbar1_about(self):
        self.window_about = about()
        self.setCentralWidget(self.window_about)

    #2 def Toolbar_action: window protein
    def trigger_toolbar2_protein(self):

        #1 Set protein base widget and layout for sub-widgets
        ##1 Set the mainwidget to put in the mainwindow
        widget_protein_l1 = QWidget()                                   # widget_protein_layer1
        self.setCentralWidget(widget_protein_l1)   
        layout_protein_l1 = QVBoxLayout()               
        widget_protein_l1.setLayout(layout_protein_l1)                  # layout_protein_layer1 (layout inside widget_l1)
        layout_protein_l1.setContentsMargins(0, 0, 0, 0)

        ##2 Set the scroll area to put in the widget protein layer 1
        widget_protein_l2 = QScrollArea()                               # widget_protein_layer2 (scroll)
        widget_protein_l2.setWidgetResizable(True)                 
        layout_protein_l1.addWidget(widget_protein_l2)                  # add widget to layout layer1

        ##3 Set the widget 3 to put in the widget protein layer 2
        widget_protein_l3 = QWidget()                                   # widget_protein_layer3
        self.layout_protein_l3 = QVBoxLayout()
        self.layout_protein_l3.setSpacing(0)
        self.layout_protein_l3.setContentsMargins(5, 5, 5, 5)
        widget_protein_l3.setLayout(self.layout_protein_l3)             # layout_protein_layer3
        widget_protein_l2.setWidget(widget_protein_l3)                  # add widget to widget layer2 (1 scroll can only has 1 widget)

        ##4 Set the widget 4 to put in the widget protein layer 3
        self.widget_protein_l4 = QWidget()
        self.layout_protein_l4 = QVBoxLayout()
        self.layout_protein_l4.setSpacing(0)
        self.layout_protein_l4.setContentsMargins(0, 0, 0, 0)
        self.widget_protein_l4.setLayout(self.layout_protein_l4)
        self.layout_protein_l3.addWidget(self.widget_protein_l4, alignment=Qt.AlignTop)

        #2 Setup elements
        ##1 button1
        self.button1_addgroup = QPushButton('+Group')
        self.button1_addgroup.setFixedSize(QSize(80,28))
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

        ##2 button2
        self.button2_alignall = QPushButton('Align all')
        self.button2_alignall.setFixedSize(QSize(80,28))
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

        ##3 Setup layout for buttons
        self.widget_protein_buttons = QWidget()
        self.layout_protein_buttons = QHBoxLayout()
        self.layout_protein_buttons.setContentsMargins(0,2,0,2)
        self.widget_protein_buttons.setLayout(self.layout_protein_buttons)
        self.layout_protein_buttons.addWidget(self.button1_addgroup, alignment=Qt.AlignLeft)
        self.layout_protein_buttons.addWidget(self.button2_alignall, alignment=Qt.AlignLeft)
        self.layout_protein_buttons.addStretch(1)

        ##4 Setup widget for secondary alignment
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
        
        #3 Add elements to the layout
        self.layout_protein_l4.addWidget(self.widget_protein_buttons, alignment=Qt.AlignTop)
        self.layout_protein_l4.addWidget(self.widget_protein_l4_2ndarystructure)

        #4 Connect elements to def
        self.button1_addgroup.clicked.connect(self.button1_addgroup_clicked)
        self.button2_alignall.clicked.connect(self.button2_alignall_clicked)

    #3 def [continue #2 def] to add new group
    def button1_addgroup_clicked(self):

        #1 Setup main alignment box --> add to main layout after button1
        widget_protein_l4_group_l1 = QWidget()                                                  # widget_protein_l3_group_layer1
        layout_protein_l4_group_l1 = QVBoxLayout()
        layout_protein_l4_group_l1.setContentsMargins(0,2,0,2)                                               
        widget_protein_l4_group_l1.setLayout(layout_protein_l4_group_l1)                        # layout_protein_l3_group_layer1 for widget layer1
        self.layout_protein_l4.addWidget(widget_protein_l4_group_l1)                            # add this widget to the main_widget_layout

        #2 Setup LINE 1
        ##1 Set widget and layout
        widget_protein_l4_group_l1_line1 = QWidget()                                            # widget for line 1
        widget_protein_l4_group_l1_line1.setObjectName("line_1")
        widget_protein_l4_group_l1_line1.setStyleSheet(
            """ #line_1 {
                border: 1px solid #A9A9A9;
                border-radius: 6px;
                padding: 2px;
            } """
        )
        layout_protein_l4_group_l1_line1 = QHBoxLayout()
        widget_protein_l4_group_l1_line1.setLayout(layout_protein_l4_group_l1_line1)
        layout_protein_l4_group_l1.addWidget(widget_protein_l4_group_l1_line1, alignment=Qt.AlignTop)

        ##2 Setup elements
        lineedit_groupname = QLineEdit('edit group name')
        lineedit_groupname.setFixedWidth(120)
        button2_removegroup = QPushButton('-Group')
        button2_removegroup.setFixedWidth(100)
        button3_addseq = QPushButton('+')
        button3_addseq.setFixedWidth(50)
        button4_removeseq = QPushButton('-')
        button4_removeseq.setFixedWidth(50)
        button5_align = QPushButton('Align')
        button5_align.setFixedWidth(100)

        ##3 Add elements to the layout
        layout_protein_l4_group_l1_line1.addWidget(lineedit_groupname)
        layout_protein_l4_group_l1_line1.addWidget(button2_removegroup)
        layout_protein_l4_group_l1_line1.addWidget(button3_addseq)
        layout_protein_l4_group_l1_line1.addWidget(button4_removeseq)
        layout_protein_l4_group_l1_line1.addWidget(button5_align)
        layout_protein_l4_group_l1_line1.addStretch(1)                      # to left-align all elements





    #3 def Toolbar_action: window codon
        

# Execute
window = main()
window.show()
app.exec()