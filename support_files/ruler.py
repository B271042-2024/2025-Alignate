from PySide6.QtWidgets import QWidget, QHBoxLayout, QVBoxLayout, QCheckBox, QLabel, QSizePolicy, QDialog, QPushButton
from PySide6.QtCore import Qt, Signal
from PySide6.QtGui import QMouseEvent, QPainter, QPen


##---------------------------------------------------------clickable labels


class ClickableLabel(QLabel):
    clicked = Signal()

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setCursor(Qt.PointingHandCursor)
        self.marked = False                                            # 1 red mark


    def mousePressEvent(self, event: QMouseEvent):
        if event.button() == Qt.LeftButton:
            self.clicked.emit()

    def paintEvent(self, event):                                       # 2 red mark (the whole def)
        super().paintEvent(event)

        if self.marked:
            painter = QPainter(self)
            painter.setRenderHint(QPainter.Antialiasing)
            pen = QPen(Qt.red, 2)
            painter.setPen(pen)
            radius = 6
            center = self.rect().center()
            painter.drawEllipse(center, radius, radius)


##---------------------------------------------------------ruler


class ruler(QWidget):
    def __init__(self, parent_context, parent=None):
#    def __init__(self, context, parent=None):
        super().__init__(parent)

#        self.context = context
        self.context = parent_context

        # call fxn 1: build the ruler
        self.build_ruler()

        # call Class 2, fxn : set clickable
        self.clickable = ClickableLabel()
        self.layout_ruler.addWidget(self.clickable)
        #self.clickable.clicked.connect(self.handle_clickable)
        self.clickable.clicked.connect(lambda: self.handle_clickable(-1))
        self.click_count = 0
        self.positions = []

    def build_ruler(self):
        # 1 set main widget & layout
        self.setFixedHeight(20)
        self.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        self.setObjectName('ruler')
        layout_ruler_main = QHBoxLayout()
        layout_ruler_main.setContentsMargins(0,0,0,0)
        layout_ruler_main.setSpacing(0)
        self.setLayout(layout_ruler_main)

        # 2 add checkbox
        invisible_checkbox = QCheckBox()
        invisible_checkbox.setEnabled(False)
        invisible_checkbox.setStyleSheet('border: none; background: transparent')
        layout_ruler_main.addWidget(invisible_checkbox, alignment=Qt.AlignLeft)

        # 3 add label
        invisible_label = QLabel('')
        invisible_label.setFixedSize(115,20)
        layout_ruler_main.addWidget(invisible_label, alignment=Qt.AlignLeft)

        # 4 add ruler
        widget_ruler = QWidget()
        self.layout_ruler = QHBoxLayout()
        self.layout_ruler.setContentsMargins(0,0,0,0)
        self.layout_ruler.setSpacing(0)
        widget_ruler.setLayout(self.layout_ruler)
        layout_ruler_main.addWidget(widget_ruler, alignment=Qt.AlignLeft)

        self.ruler_labels = []                                                  # 3 red mark
        approx_width = 15000        # max 1000 bases
        num_cols = int(approx_width / 15)
        for i in range(num_cols):
            label = ClickableLabel(str(i) if i % 10 == 0 else '.')
            label.clicked.connect(lambda i=i: self.handle_clickable(i))
            label.setFixedSize(15,20)
            label.setAlignment(Qt.AlignCenter)
            label.setStyleSheet('color: gray; font-size: 8px')
            self.layout_ruler.addWidget(label)
            self.ruler_labels.append(label)                                     # 4 red mark


    def handle_clickable(self, position):

        label = self.ruler_labels[position]                                  # 5 red mark
        label.marked = True                                                  # 6 red mark
        label.update()                                                       # 7 red mark

        self.click_count += 1
        self.positions.append(position)

        if self.click_count == 2:
            # positions
            pos1, pos2 = self.positions

            # create a dialog
            widget_dialog = QDialog(self)
            layout_dialog = QVBoxLayout()
            widget_dialog.setLayout(layout_dialog)

            widget_label1 = QLabel("WARNING! Only for aligned sequences.")
            widget_label2 = QLabel("---")
            widget_label3 = QLabel(f"Show '%'conservation for selected position: {pos1}, {pos2}?")
            widget_button = QWidget()
            layout_button = QHBoxLayout()
            widget_button.setLayout(layout_button)
            button_yes = QPushButton('Yes')
            button_cancel = QPushButton('Cancel')
            layout_button.addWidget(button_yes)
            layout_button.addWidget(button_cancel)

            layout_dialog.addWidget(widget_label1)
            layout_dialog.addWidget(widget_label2)
            layout_dialog.addWidget(widget_label3)
            layout_dialog.addWidget(widget_button)

            def reset_state():
                for idx in [pos1, pos2]:
                    self.ruler_labels[idx].marked = False
                    self.ruler_labels[idx].update()
                self.click_count = 0
                self.positions.clear()

            print(f'test: {position}')

            button_cancel.clicked.connect(widget_dialog.reject)
            button_cancel.clicked.connect(reset_state)

            button_yes.clicked.connect(lambda: (widget_dialog.accept(), self.context.custom_display_perc_cons(pos1, pos2)))
            button_yes.clicked.connect(reset_state)

            widget_dialog.exec()

            for label in self.ruler_labels:                                      # 8 red mark
                label.marked = False
                label.update()



    def tedelhandle_clickable(self, position):

        label = self.ruler_labels[position]                                  # 5 red mark
        label.marked = True                                                  # 6 red mark
        label.update()                                                       # 7 red mark

        self.click_count += 1
        self.positions.append(position)

        if self.click_count == 2:
            # positions
            pos1, pos2 = self.positions
            self.click_count = 0
            self.positions.clear()
            # create a dialog
            widget_dialog = QDialog(self)
            layout_dialog = QVBoxLayout()
            widget_dialog.setLayout(layout_dialog)

            widget_label = QLabel(f"Show '%'conservation for selected position: {pos1}, {pos2}?")
            widget_button = QWidget()
            layout_button = QHBoxLayout()
            widget_button.setLayout(layout_button)
            button_yes = QPushButton('Yes')
            button_cancel = QPushButton('Cancel')
            layout_button.addWidget(button_yes)
            layout_button.addWidget(button_cancel)

            layout_dialog.addWidget(widget_label)
            layout_dialog.addWidget(widget_button)

            count = 0       # clears count
            print(f'test: {position}')

            button_cancel.clicked.connect(widget_dialog.reject)
            button_yes.clicked.connect(lambda: (widget_dialog.accept(), self.context.custom_display_perc_cons(pos1, pos2)))

            label = self.ruler_labels[position]                                  # 5 red mark
            label.marked = True                                                  # 6 red mark
            label.update()                                                       # 7 red mark

            widget_dialog.exec()

            for label in self.ruler_labels:                                      # 8 red mark
                label.marked = False
                label.update()

