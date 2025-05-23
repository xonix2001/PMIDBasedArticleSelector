# main.py

import sys
import json
import os
import time
from pathlib import Path
from datetime import datetime

from PySide6.QtWidgets import (
    QApplication, QWidget, QFileDialog, QMessageBox
)
from PySide6.QtUiTools import QUiLoader
from PySide6.QtCore import QFile, QIODevice, QObject, Signal, QThread
from PySide6.QtGui import QTextCursor, QIcon

from base import PubMedDataProcessorBase
from sort import PmidFilter
from extend_test import PubMedDataProcessorExtended
from translator import translator
from thread_worker import ExtendedWorker
# 配置文件路径
CONFIG_PATH = Path.home() / ".pmid_selector_config.json"

class EmittingStream(QObject):
    """把 stdout 重定向到 GUI 文本区"""
    textWritten = Signal(str)
    def write(self, text):      self.textWritten.emit(text)
    def flush(self):            pass

class MainWindow(QWidget):
    # 供 Extended 模块发射
    progress_signal = Signal(int, int)
    eta_signal      = Signal(float)

    def __init__(self):
        super().__init__()
        self._load_ui()
        self._setup_logging()
        self._connect_signals()
        self._load_config()
        #self.show()

        # 存储三个模块的文件路径
        self.file1 = None
        self.file2 = None
        self.file3 = None

    def _load_ui(self):
        loader = QUiLoader()
        ui_file = QFile("test.ui")
        if not ui_file.open(QIODevice.ReadOnly):
            QMessageBox.critical(self, "错误", "无法加载 UI！")
            sys.exit(-1)
        self.ui = loader.load(ui_file, self)
        ui_file.close()
        self.setWindowTitle(self.ui.windowTitle())

    def _setup_logging(self):
        self.stream = EmittingStream()
        self.stream.textWritten.connect(self._append_log)
        sys.stdout = self.stream

    def _connect_signals(self):
        # 全局配置
        self.ui.pushButton_save_config.clicked.connect(self._save_config)

        # -------- 模块一：基础信息 --------
        self.ui.pushButton.clicked.connect(self._browse_file1)      # 选择文件1
        self.ui.pushButton_2.clicked.connect(self._run_basic_info)  # 获取1

        # -------- 模块二：筛选高被引 --------
        self.ui.pushButton_4.clicked.connect(self._browse_file2)    # 选择文件2
        self.ui.pushButton_6.clicked.connect(self._run_filter)      # 获取2

        # -------- 模块三：摘要 & 链接 --------
        self.ui.pushButton_7.clicked.connect(self._browse_file3)    # 选择文件3
        self.ui.pushButton_8.clicked.connect(self._run_extended)    # 获取3

        # Extended 进度 & ETA
        self.progress_signal.connect(self._update_progress)
        self.eta_signal.connect(self._update_eta)

    def _load_config(self):
        if CONFIG_PATH.exists():
            try:
                cfg = json.loads(CONFIG_PATH.read_text("utf-8"))
                self.ui.lineEdit_cfg_email.setText(cfg.get("email",""))
                self.ui.lineEdit_cfg_pubmed_key.setText(cfg.get("pubmed_api_key",""))
                self.ui.lineEdit_cfg_trans_base.setText(cfg.get("translator_api_base",""))
                self.ui.lineEdit_cfg_trans_key.setText(cfg.get("translator_api_key",""))
                model = cfg.get("translator_model","gpt-4o-mini")
                idx = self.ui.comboBox_cfg_trans_model.findText(model)
                if idx>=0:
                    self.ui.comboBox_cfg_trans_model.setCurrentIndex(idx)
            except Exception as e:
                QMessageBox.warning(self, "配置加载失败", str(e))

    def _save_config(self):
        cfg = {
            "email": self.ui.lineEdit_cfg_email.text().strip(),
            "pubmed_api_key": self.ui.lineEdit_cfg_pubmed_key.text().strip(),
            "translator_api_base": self.ui.lineEdit_cfg_trans_base.text().strip(),
            "translator_api_key": self.ui.lineEdit_cfg_trans_key.text().strip(),
            "translator_model": self.ui.comboBox_cfg_trans_model.currentText(),
        }
        try:
            CONFIG_PATH.write_text(json.dumps(cfg, ensure_ascii=False, indent=2), "utf-8")
            QMessageBox.information(self, "配置已保存", str(CONFIG_PATH))
        except Exception as e:
            QMessageBox.critical(self, "保存失败", str(e))

    # ------ 模块一：基础信息 ------
    def _browse_file1(self):
        fn, _ = QFileDialog.getOpenFileName(self, "选择文件1", "", "Text Files (*.txt)")
        if fn:
            self.file1 = fn
            print(f"[Module1] 选定文件1: {fn}")

    def _run_basic_info(self):
        if not self.file1:
            QMessageBox.warning(self, "提示", "请先选择文件1")
            return
        email = self.ui.lineEdit_cfg_email.text().strip()
        api_key = self.ui.lineEdit_cfg_pubmed_key.text().strip()
        processor = PubMedDataProcessorBase(email, api_key)
        print(f"[{datetime.now():%H:%M:%S}] 开始模块1：基础信息爬取")
        try:
            processor.run(input_file_path=self.file1, output_file_path=None)
            QMessageBox.information(self, "完成", "模块1 执行完成")
        except Exception as e:
            QMessageBox.critical(self, "模块1 出错", str(e))

    # ------ 模块二：筛选高被引 ------
    def _browse_file2(self):
        fn, _ = QFileDialog.getOpenFileName(self, "选择文件2", "", "Excel Files (*.xlsx)")
        if fn:
            self.file2 = fn
            self.ui.lineEdit.setText(fn)

    def _run_filter(self):
        if not self.file2:
            QMessageBox.warning(self, "提示", "请先选择文件2")
            return
        # 构建 PmidFilter
        bad_list = [j.strip() for j in self.ui.lineEdit_3.text().split(";") if j.strip()]
        filter_obj = PmidFilter(self.file2, bad_list)
        # 读取用户选择
        cit_per = self.ui.spinBox.value()
        journal = self.ui.lineEdit_3.text().strip()
        # 判断选项
        if self.ui.checkBox.isChecked() and self.ui.checkBox_2.isChecked():
            mode = "both"
        elif self.ui.checkBox.isChecked():
            mode = "citation"
        elif self.ui.checkBox_2.isChecked():
            mode = "journal"
        else:
            QMessageBox.warning(self, "提示", "请至少选择一种筛选方式")
            return
        print(f"[{datetime.now():%H:%M:%S}] 开始模块2：筛选 (mode={mode})")
        try:
            if mode=="citation":
                df = filter_obj.filter_by_citation(cit_per)
            elif mode=="journal":
                df = filter_obj.filter_by_journal([journal])
            else:
                df0 = filter_obj.filter_by_journal([journal])
                df = filter_obj.filter_by_citation(cit_per, df0)
            filter_obj.output_result(df)
            QMessageBox.information(self, "完成", "模块2 执行完成")
        except Exception as e:
            QMessageBox.critical(self, "模块2 出错", str(e))

    # ------ 模块三：摘要 & 链接 (Extended) ------
    def _browse_file3(self):
        fn, _ = QFileDialog.getOpenFileName(self, "选择文件3", "", "Text Files (*.txt)")
        if fn:
            self.file3 = fn
            self.ui.lineEdit_2.setText(fn)

    def _run_extended(self):
        if not self.file3:
            QMessageBox.warning(self, "提示", "请先选择文件3")
            return

        self.ui.plainTextEdit.clear()
        self.ui.label_progress.setText("[0/0]")
        self.ui.label_eta.setText("ETA: --:--")

        email = self.ui.lineEdit_cfg_email.text().strip()
        api_key = self.ui.lineEdit_cfg_pubmed_key.text().strip()

        # 更新翻译参数
        translator.api_base = self.ui.lineEdit_cfg_trans_base.text().strip()
        translator.api_key  = self.ui.lineEdit_cfg_trans_key.text().strip()
        translator.model    = self.ui.comboBox_cfg_trans_model.currentText()

        print(f"[{datetime.now():%H:%M:%S}] 开始模块3：摘要 & 链接爬取")

        # 启动 Worker 线程
        self.thread = QThread()
        self.worker = ExtendedWorker(email, api_key, self.file3)
        self.worker.moveToThread(self.thread)

        # 连接信号槽
        self.thread.started.connect(self.worker.run)
        self.worker.progress_signal.connect(self._update_progress)
        self.worker.eta_signal.connect(self._update_eta)
        self.worker.finished.connect(self.thread.quit)
        self.worker.finished.connect(lambda: QMessageBox.information(self, "完成", "模块3 执行完成"))
        self.worker.error.connect(lambda msg: QMessageBox.critical(self, "模块3 出错", msg))
        self.worker.error.connect(self.thread.quit)
        self.thread.finished.connect(self.thread.deleteLater)

        self.thread.start()
    # ---- 日志 & 进度回调 ----
    def _append_log(self, txt):
        self.ui.plainTextEdit.moveCursor(QTextCursor.End)
        self.ui.plainTextEdit.appendPlainText(txt.rstrip())

    def _update_progress(self, done, total):
        self.ui.label_progress.setText(f"[{done}/{total}]")

    def _update_eta(self, sec):
        m, s = divmod(int(sec), 60)
        self.ui.label_eta.setText(f"ETA: {m:02d}:{s:02d}")

    def closeEvent(self, e):
        sys.stdout = sys.__stdout__
        super().closeEvent(e)


if __name__ == "__main__":
    app = QApplication(sys.argv)
    app.setWindowIcon(QIcon('logo.png'))
    w = MainWindow()
    w.show()
    sys.exit(app.exec())
