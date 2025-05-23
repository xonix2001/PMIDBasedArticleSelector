# thread_worker.py

from PySide6.QtCore import QObject, Signal, Slot
from extend_test import PubMedDataProcessorExtended

class ExtendedWorker(QObject):
    finished = Signal()
    error = Signal(str)
    progress_signal = Signal(int, int)
    eta_signal = Signal(float)

    def __init__(self, email, api_key, filepath):
        super().__init__()
        self.email = email
        self.api_key = api_key
        self.filepath = filepath

    @Slot()
    def run(self):
        try:
            proc = PubMedDataProcessorExtended(
                email=self.email,
                api_key=self.api_key,
                progress_signal=self.progress_signal,
                eta_signal=self.eta_signal
            )
            proc.run(self.filepath)
            self.finished.emit()
        except Exception as e:
            self.error.emit(str(e))
