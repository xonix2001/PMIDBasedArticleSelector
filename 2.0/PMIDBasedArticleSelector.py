import sys
from PySide6.QtWidgets import QApplication, QMainWindow, QFileDialog, QMessageBox
from PySide6.QtUiTools import QUiLoader
from PySide6.QtCore import QFile, QIODevice, QObject, Signal
from PySide6.QtGui import QTextCursor, QIcon
import pandas as pd
from base import PubMedDataProcessorBase
from sort import PmidFilter
from extend import PubMedDataProcessorExtended

class EmittingStream(QObject):
    textWritten = Signal(str)  # Signal to emit the captured text

    def write(self, text):
        self.textWritten.emit(str(text))

    def flush(self):
        pass  # Placeholder to avoid errors when calling flush()

class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.load_ui()
        # Connecting UI components
        self.ui.pushButton.clicked.connect(self.open_file_dialog)
        self.ui.pushButton_2.clicked.connect(self.execute_core_logic)
        self.ui.pushButton_4.clicked.connect(self.select_file1)
        self.ui.pushButton_6.clicked.connect(self.run_filter)
        self.ui.pushButton_7.clicked.connect(self.select_file2)
        self.ui.pushButton_8.clicked.connect(self.run_processing)
        self.ui.lineEdit_3.textChanged.connect(self.text_changed)
        self.ui.comboBox.activated.connect(self.onComboBoxChanged)
        self.selected_file_path = None
        self.file_path1 = None
        self.file_path2 = None
        self.filter = None
        self.emit_stream = EmittingStream()
        self.emit_stream.textWritten.connect(self.appendTextToTextEdit)
        sys.stdout = self.emit_stream  # Redirect sys.stdout to custom stream
        self.show()

    def appendTextToTextEdit(self, text):
        cleaned_text = text.strip()
        self.ui.plainTextEdit.moveCursor(QTextCursor.End)
        self.ui.plainTextEdit.appendPlainText(cleaned_text)

    def load_ui(self):
        loader = QUiLoader()
        file = QFile("get_shortinfo.ui")
        if not file.open(QIODevice.ReadOnly):
            QMessageBox.critical(self, "Error", "Cannot load UI file.")
            sys.exit(-1)
        self.ui = loader.load(file, self)
        file.close()

    # 第一段逻辑代码
    def open_file_dialog(self):
        file_name, _ = QFileDialog.getOpenFileName(self, "Open File", "", "Text Files (*.txt)")
        if file_name:
            self.selected_file_path = file_name
            self.ui.plainTextEdit.appendPlainText(f"待处理文件: {file_name}"+ '\n')

    def execute_core_logic(self):
        if self.selected_file_path:
            email = ""
            api_key = ""
            processor = PubMedDataProcessorBase(email, api_key)
            processor.run(input_file_path=self.selected_file_path, 
                        output_file_path=None)
            QMessageBox.information(self, "完成", "数据处理完成。")
        else:
            QMessageBox.warning(self, "警告", "请先选择一个文件。")

    # 第二段逻辑代码
    def select_file1(self):
        file_path1, _ = QFileDialog.getOpenFileName(self, "Open File", "", "Excel Files (*.xlsx)")
        if file_path1:
            self.ui.lineEdit.setText(file_path1)
            self.file_path1 = file_path1
            bad_journal_list = ['CANCERS', 
                    'DIAGNOSTICS', 
                    'ENVIRONMENTAL SCIENCE AND POLLUTION RESEARCH', 
                    'FUEL', 
                    'JOURNAL OF CLINICAL MEDICINE', 
                    'JOURNAL OF PERSONALIZED MEDICINE', 
                    'RADIOLOGIA MEDICA',
                    'BIOENGINEERED', 
                    'CONNECTION SCIENCE', 
                    'MULTIMEDIA TOOLS AND APPLICATIONS', 
                    'PSYCHIATRIA DANUBINA', 
                    'JOURNAL OF BIOBASED MATERIALS AND BIOENERGY', 
                    'JOURNAL OF BIOMATERIALS AND TISSUE ENGINEERING', 
                    'JOURNAL OF BIOMEDICAL NANOTECHNOLOGY', 
                    'JOURNAL OF NANOELECTRONICS AND OPTOELECTRONICS', 
                    'JOURNAL OF SENSORS', 
                    'MATERIALS EXPRESS', 
                    'SCIENCE OF ADVANCED MATERIALS', 
                    'ALTERNATIVE THERAPIES IN HEALTH AND MEDICINE', 
                    'CMES-COMPUTER MODELING IN ENGINEERING & SCIENCES', 
                    'EXPERIMENTAL AND THERAPEUTIC MEDICINE', 
                    'FRONTIERS IN ENERGY RESEARCH', 
                    'MATHEMATICAL BIOSCIENCES AND ENGINEERING', 
                    'TROPICAL JOURNAL OF PHARMACEUTICAL RESEARCH']
            self.filter = PmidFilter(file_path1, bad_journal_list)

    def run_filter(self):
        if not self.filter:
            QMessageBox.warning(self, "Warning", "Please load a file first.")
            return
        # 读取用户在GUI中的输入
        citation_filter = self.ui.spinBox.value()
        journal_filter = self.ui.lineEdit_3.text()
        interested_journals = [journal_filter.strip()] if journal_filter else []

        # 映射checkBox选择到output_choice
        if self.ui.checkBox.isChecked() and self.ui.checkBox_2.isChecked():
            output_choice = '3'
        elif self.ui.checkBox.isChecked():
            output_choice = '2'
        elif self.ui.checkBox_2.isChecked():
            output_choice = '1'
        else:
            QMessageBox.warning(self, "Warning", "Please select an output option.")
            return

        # 运行PubMedFilter逻辑
        if output_choice == '1':
            df_interested_journals = self.filter.filter_by_journal(interested_journals)
            self.filter.output_result(df_interested_journals)
        elif output_choice == '2':
            df_final = self.filter.filter_by_citation(citation_filter)
            self.filter.output_result(df_final)
        elif output_choice == '3':
            df_interested_journals = self.filter.filter_by_journal(interested_journals)
            df_final_journals_citations = self.filter.filter_by_citation(citation_filter, df_interested_journals)
            self.filter.output_result(df_final_journals_citations)

    def onComboBoxChanged(self,index):
        self.combobox_item=self.ui.comboBox.itemText(index)
        self.ui.lineEdit_3.setText(self.combobox_item)

    def text_changed(self):
        self.df=pd.read_excel(self.file_path1, engine='openpyxl')
        self.journal_list=self.df['Journal'].tolist()
        self.edit_text=self.ui.lineEdit_3.text()
        matched_element=[item for item in self.journal_list if self.edit_text.lower() in item.lower()]
        self.ui.comboBox.setItems(matched_element)

    # 第三段逻辑代码
    def select_file2(self):
        file_name, _ = QFileDialog.getOpenFileName(self, "Open File", "", "Text Files (*.txt)")
        if file_name:
            self.file_path2 = file_name
            self.ui.lineEdit_2.setText(file_name)

    def run_processing(self):
        if self.file_path2:
            email = ""
            api_key = ""
            processor = PubMedDataProcessorExtended(email, api_key)
            try:
                processor.run(input_file_path=self.file_path2)
                QMessageBox.information(self, "Completed", "Processing completed successfully.")
            except Exception as e:
                QMessageBox.critical(self, "Error", f"An error occurred: {e}")
        else:
            QMessageBox.warning(self, "Warning", "Please select a file first.")

    def closeEvent(self, event):
        sys.stdout = sys.__stdout__  # Restore sys.stdout upon closing the app
        super().closeEvent(event)

if __name__ == "__main__":
    app = QApplication(sys.argv)
    app.setWindowIcon(QIcon('logo.png'))
    window = MainWindow()
    window.show()
    sys.exit(app.exec())
