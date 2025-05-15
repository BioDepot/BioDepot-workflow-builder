"""
Preview item model.
"""

import logging

from AnyQt.QtWidgets import QApplication, QStyleOption
from AnyQt.QtGui import (
    QStandardItemModel,
    QStandardItem,
    QIcon,
    QIconEngine,
    QPainter,
    QPixmap,
)
from AnyQt.QtSvg import QSvgRenderer

# pylint: disable=unused-import
from AnyQt.QtCore import Qt, QTimer, QRectF, QRect, QSize


from . import scanner

log = logging.getLogger(__name__)

# Preview Data Roles
####################

# Name of the item, (same as `Qt.DisplayRole`)
NameRole = Qt.DisplayRole

# Items description (items long description)
DescriptionRole = Qt.UserRole + 1

# Items url/path (where previewed resource is located).
PathRole = Qt.UserRole + 2

# Items preview SVG contents string
ThumbnailSVGRole = Qt.UserRole + 3


UNKNOWN_SVG = """<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<svg width="161.8mm" height="100.0mm"
 xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink"
 version="1.2" baseProfile="tiny">
</svg>
"""


class PreviewModel(QStandardItemModel):
    """A model for preview items.
    """

    def __init__(self, parent=None, items=None):
        QStandardItemModel.__init__(self, parent)

        if items is not None:
            self.insertColumn(0, items)

        self.__timer = QTimer(self)

    def delayedScanUpdate(self, delay=10):
        """Run a delayed preview item scan update.
        """

        def iter_update(items):
            for item in items:
                try:
                    scanner.scan_update(item)
                except Exception:
                    log.error(
                        "An unexpected error occurred while " "scanning %r.",
                        str(item.text()),
                        exc_info=True,
                    )
                    item.setEnabled(False)
                yield

        items = [self.item(i) for i in range(self.rowCount())]

        iter_scan = iter_update(items)

        def process_one():
            try:
                next(iter_scan)
            except StopIteration:
                self.__timer.timeout.disconnect(process_one)
                self.__timer.stop()

        self.__timer.timeout.connect(process_one)
        self.__timer.start(delay)


class PreviewItem(QStandardItem):
    """A preview item.
    """

    def __init__(
        self, name=None, description=None, thumbnail=None, icon=None, path=None
    ):
        QStandardItem.__init__(self)

        self.__name = ""

        if name is None:
            name = "Untitled"

        self.setName(name)

        if description is None:
            description = "No description."
        self.setDescription(description)

        if thumbnail is None:
            thumbnail = UNKNOWN_SVG
        self.setThumbnail(thumbnail)

        if icon is not None:
            self.setIcon(icon)

        if path is not None:
            self.setPath(path)

    def name(self):
        """Return the name (title) of the item (same as `text()`.
        """
        return self.__name

    def setName(self, value):
        """Set the item name. `value` if not empty will be used as
        the items DisplayRole otherwise an 'untitled' placeholder will
        be used.

        """
        self.__name = value

        if not value:
            self.setText("untitled")
        else:
            self.setText(value)

    def description(self):
        """Return the detailed description for the item.

        This is stored as `DescriptionRole`, if no data is set then
        return the string for `WhatsThisRole`.

        """
        desc = self.data(DescriptionRole)

        if desc is not None:
            return desc

        whatsthis = self.data(Qt.WhatsThisRole)
        return whatsthis

    def setDescription(self, description):
        self.setData(description, DescriptionRole)
        self.setWhatsThis(description)

    def thumbnail(self):
        """Return the thumbnail SVG string for the preview item.

        This is stored as `ThumbnailSVGRole`
        """
        thumb = self.data(ThumbnailSVGRole)
        if thumb is not None:
            return thumb

    def setThumbnail(self, thumbnail):
        """Set the thumbnail SVG contents as a string.

        When set it also overrides the icon role.

        """
        self.setData(thumbnail, ThumbnailSVGRole)
        engine = SvgIconEngine(thumbnail.encode("utf-8"))
        self.setIcon(QIcon(engine))

    def path(self):
        """Return the path item data.
        """
        return self.data(PathRole)

    def setPath(self, path):
        """Set the path data of the item.

        .. note:: This also sets the Qt.StatusTipRole

        """
        self.setData(path, PathRole)
        self.setStatusTip(path)
        self.setToolTip(path)


class SvgIconEngine(QIconEngine):
    def __init__(self, contents):
        # type: (bytes) -> None
        super().__init__()
        self.__contents = contents
        self.__generator = QSvgRenderer(contents)

    def paint(self, painter, rect, mode, state):
        # type: (QPainter, QRect, QIcon.Mode, QIcon.State) -> None
        if self.__generator.isValid():
            size = rect.size()
            dpr = 1.0
            try:
                dpr = painter.device().devicePixelRatioF()
            except AttributeError:
                pass
            if dpr != 1.0:
                size = size * dpr
            painter.drawPixmap(rect, self.pixmap(size, mode, state))

    def pixmap(self, size, mode, state):
        # type: (QSize, QIcon.Mode, QIcon.State) -> QPixmap
        if not self.__generator.isValid():
            return QPixmap()

        dsize = self.__generator.defaultSize()  # type: QSize
        if not dsize.isNull():
            dsize.scale(size, Qt.KeepAspectRatio)
            size = dsize

        pm = QPixmap(size)
        pm.fill(Qt.transparent)
        painter = QPainter(pm)
        try:
            self.__generator.render(painter, QRectF(0, 0, size.width(), size.height()))
        finally:
            painter.end()
        style = QApplication.style()
        if style is not None:
            opt = QStyleOption()
            opt.palette = QApplication.palette()
            pm = style.generatedIconPixmap(mode, pm, opt)
        return pm
