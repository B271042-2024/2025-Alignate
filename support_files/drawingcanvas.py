from PySide6.QtWidgets import QWidget, QSizePolicy
from PySide6.QtGui import QPainter, QPen, QImage
from PySide6.QtCore import Qt, QPoint, QSize


#_______________________________________________________________________________________________
#_______________________________________________________________________________________________
#_______________________________________________________________________________________________1 CLASS: Drawing pad
#_______________________________________________________________________________________________
#_______________________________________________________________________________________________
class DrawingCanvas(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setAttribute(Qt.WA_StaticContents)
        self.drawing = False
        self.last_point = QPoint()
        self.setFixedHeight(40)  # Only fix height
        self.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        self.image = None

    def resizeEvent(self, event):
        new_width = self.width()
        new_height = self.height()

        if self.image is None:
            self.image = QImage(new_width, new_height, QImage.Format_RGB32)
            self.image.fill(Qt.white)
        else:
            # Resize and preserve existing content
            new_image = QImage(new_width, new_height, QImage.Format_RGB32)
            new_image.fill(Qt.white)
            painter = QPainter(new_image)
            painter.drawImage(0, 0, self.image)
            self.image = new_image

        super().resizeEvent(event)

    def paintEvent(self, event):
        painter = QPainter(self)
        painter.drawImage(0, 0, self.image)

    def mousePressEvent(self, event):
        if event.button() in (Qt.LeftButton, Qt.RightButton):
            self.drawing = True
            self.last_point = event.position().toPoint()

    def mouseMoveEvent(self, event):
        if self.drawing:
            current_point = event.position().toPoint()
            painter = QPainter(self.image)
            if event.buttons() & Qt.LeftButton:
                color = Qt.black  # Draw normally
                pen = QPen(Qt.black, 2, Qt.SolidLine)
            elif event.buttons() & Qt.RightButton:
                pen = QPen(Qt.white, 10, Qt.SolidLine)
            else:
                return  # Ignore other buttons
            
            painter.setPen(pen)
            painter.drawLine(self.last_point, current_point)
            self.last_point = current_point
            self.update()

    def mouseReleaseEvent(self, event):
        if event.button() == Qt.LeftButton:
            self.drawing = False
