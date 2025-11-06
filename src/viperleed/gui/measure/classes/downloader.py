"""Module downloader of viperleed.gui.measure.classes.

Contains NetworkHelper, a downloader that uses QtNetwork or requests as
fallback to download data from a given url.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    'Florian Dörr (@FlorianDoerr)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2025-11-04'
__license__ = 'GPLv3+'

import importlib

from PyQt5 import QtCore as qtc
from PyQt5 import QtNetwork as qtn


class _RequestsWorker(qtc.QObject):
    """Downloder using requests module."""

    # Emitted if download finished successfully.
    finished = qtc.pyqtSignal(bytes)
    # Emitted if download failed.
    failed = qtc.pyqtSignal()

    def __init__(self, url):
        """Initialise downloader."""
        super().__init__()
        self._url = url

    @qtc.pyqtSlot()
    def run(self):
        """Download from self._url."""
        try:
            requests = importlib.import_module('requests')
        except ImportError:
            self.failed.emit()
            return

        try:
            response = requests.get(self._url, timeout=2)
        except requests.exceptions.RequestException:
            self.failed.emit()
            return

        try:
            response.raise_for_status()
        except requests.exceptions.RequestException:
            self.failed.emit()
            return

        self.finished.emit(response.content)


class NetworkHelper(qtc.QObject):
    """Abstraction for QtNetwork with requests fallback."""

    # Emitted if the download was successful. Carries downloaded data.
    download_finished = qtc.pyqtSignal(bytes)

    # Emitted if the download failed.
    download_failed = qtc.pyqtSignal()

    # Emitted if the requests module must be installed. This is usually
    # the case if Qt5 QtNetwork is incompatible with the OS.
    install_requests = qtc.pyqtSignal()

    def __init__(self, parent=None):
        """Initialise the NetworkHelper."""
        super().__init__(parent=parent)
        # Notice that we do not give self as a parent for self._network.
        # That's because moving a QNetworkAccessManager automatically
        # to a new thread together with its parent (at is normally the
        # case with parent-child relationships) seems to be broken in
        # Qt5. We move self._network explicitly in self.moveToThread.
        # Assigning a parent here would break the overridden
        # self.moveToThread, as children with a parent cannot be moved
        # 'independently'.
        self._network = qtn.QNetworkAccessManager()
        self._network.setTransferTimeout(timeout=2000)
        self._thread = None
        self._worker = None

    def download(self, url):
        """Start a download from url.

        Choose QtNetwork or requests dynamically depending on system
        availability. QtNetwork is preferred. Missing QtNetwork support
        is detected automatically.

        Parameters
        ----------
        url : str
            Link to data to download.

        Returns
        -------
        None.
        """
        if not qtn.QSslSocket.supportsSsl():
            # Cannot use QtNetwork, use requests module instead.
            self._download_with_requests(url)
            return

        request = qtn.QNetworkRequest(qtc.QUrl(url))
        # Since Qt does not come with SSL support the RedirectAttribute
        # must be set in order to get the Arduino CLI file via http.
        request.setAttribute(qtn.QNetworkRequest.FollowRedirectsAttribute,
                             True)
        self._network.finished.connect(self._on_qt_finished)
        self._network.get(request)

    def moveToThread(self, thread):     # pylint: disable=invalid-name
        """Move self and self._network to thread."""
        super().moveToThread(thread)
        self._network.moveToThread(thread)

    @qtc.pyqtSlot(bytes)
    @qtc.pyqtSlot()
    def _cleanup_thread(self, *_):
        """Delete thread and worker objects."""
        if self._thread:
            self._thread.quit()
            self._thread.wait()
            self._thread.deleteLater()

        if self._worker:
            self._worker.deleteLater()

        self._worker = None
        self._thread = None

    def _download_with_requests(self, url):
        """Download data with requests module.

        Parameters
        ----------
        url : str
            Link to data to download.

        Emits
        -----
        install_requests
            If requests module is not installed.
        """
        try:
            importlib.import_module('requests')
        except ImportError:
            self.install_requests.emit()
            return
        self._thread = qtc.QThread(self)
        self._worker = _RequestsWorker(url)
        self._worker.moveToThread(self._thread)

        # Wire signals
        self._worker.finished.connect(self._cleanup_thread)
        self._worker.failed.connect(self._cleanup_thread)
        self._worker.finished.connect(self.download_finished.emit)
        self._worker.failed.connect(self.download_failed.emit)
        self._thread.started.connect(self._worker.run)

        self._thread.start()

    @qtc.pyqtSlot(qtn.QNetworkReply)
    def _on_qt_finished(self, reply):
        """Check reply and emit.

        Parameters
        ----------
        reply : qtn.QNetworkReply
            Contains downloaded data or error message.

        Emits
        -----
        download_failed
            If the download failed.
        download_finished(data)
            If the download succeeded. Carries data as bytes.
        """
        self._network.finished.disconnect(self._on_qt_finished)
        if reply.error():
            self.download_failed.emit()
        else:
            self.download_finished.emit(bytes(reply.readAll()))
        reply.deleteLater()
