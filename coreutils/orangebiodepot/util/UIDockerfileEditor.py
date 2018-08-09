from PyQt5.QtGui import (
    QBrush, QFont, QSyntaxHighlighter, QTextCharFormat
)
from PyQt5.QtCore import Qt, QRegExp

keyword_list = ['ADD', 'COPY', 'RUN', 'FROM', 'CMD', 'ENTRYPOINT',
                'VOLUME', 'LABEL', 'MAINTAINER', 'EXPOSE', 'ENV',
                'USER', 'WORKDIR', 'ARG', 'STOPSIGNAL', 'HEALTHCHECK', 'SHELL']

class DockerSyntaxHighlighter(QSyntaxHighlighter):
    def __init__(self, parent=None):

        self.keywordFormat = self._text_format(Qt.blue)
        self.stringFormat = self._text_format(Qt.darkRed)
        self.commentFormat = self._text_format(Qt.darkGreen)
        self.decoratorFormat = self._text_format(Qt.darkGray)

        self.keywords = list(keyword_list)

        self.rules = [(QRegExp(r"\b%s\b" % kwd), self.keywordFormat)
                      for kwd in self.keywords] + \
                     [(QRegExp(r"'.*'"), self.stringFormat),
                      (QRegExp(r'".*"'), self.stringFormat),
                      (QRegExp(r"#.*"), self.commentFormat),
                      (QRegExp(r"@[A-Za-z_]+[A-Za-z0-9_]+"),
                       self.decoratorFormat)]

        self.multilineStart = QRegExp(r"(''')|" + r'(""")')
        self.multilineEnd = QRegExp(r"(''')|" + r'(""")')

        super().__init__(parent)

    def highlightBlock(self, text):
        for pattern, format in self.rules:
            exp = QRegExp(pattern)
            index = exp.indexIn(text)
            while index >= 0:
                length = exp.matchedLength()
                if exp.captureCount() > 0:
                    self.setFormat(exp.pos(1), len(str(exp.cap(1))), format)
                else:
                    self.setFormat(exp.pos(0), len(str(exp.cap(0))), format)
                index = exp.indexIn(text, index + length)

        # Multi line strings
        start = self.multilineStart
        end = self.multilineEnd

        self.setCurrentBlockState(0)
        startIndex, skip = 0, 0
        if self.previousBlockState() != 1:
            startIndex, skip = start.indexIn(text), 3
        while startIndex >= 0:
            endIndex = end.indexIn(text, startIndex + skip)
            if endIndex == -1:
                self.setCurrentBlockState(1)
                commentLen = len(text) - startIndex
            else:
                commentLen = endIndex - startIndex + 3
            self.setFormat(startIndex, commentLen, self.stringFormat)
            startIndex, skip = (start.indexIn(text,startIndex + commentLen + 3), 3)

    def _text_format(self, foreground=Qt.black, weight=QFont.Normal):
        fmt = QTextCharFormat()
        fmt.setForeground(QBrush(foreground))
        fmt.setFontWeight(weight)
        return fmt