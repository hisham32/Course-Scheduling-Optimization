import sys
import base64
import json
import requests
from PySide6.QtCore import Qt, QThread, Signal
from PySide6.QtWidgets import (
    QApplication, QWidget, QVBoxLayout, QHBoxLayout, QLabel, QPushButton,
    QFileDialog, QComboBox, QTableWidget, QTableWidgetItem, QMessageBox,
    QCheckBox, QLineEdit, QHeaderView, QTextEdit
)

API_URL_DEFAULT = "http://localhost:8000"


class LogStreamThread(QThread):
    line = Signal(str)          # emits each log line
    finishedEvt = Signal(dict)  # emits final {"type":"done", ...} payload

    def __init__(self, url: str, payload: dict):
        super().__init__()
        self.url = url
        self.payload = payload

    def run(self):
        try:
            with requests.post(
                self.url,
                json=self.payload,
                stream=True,
                headers={"Accept": "text/event-stream"}
            ) as r:
                r.raise_for_status()
                # small chunk_size to avoid client-side buffering
                for raw in r.iter_lines(chunk_size=1, decode_unicode=True):
                    if raw is None:
                        continue
                    if raw.startswith(":"):  # heartbeat/comment
                        continue
                    if not raw.startswith("data:"):
                        continue
                    data = raw[len("data:"):].strip()
                    if not data:
                        continue
                    try:
                        obj = json.loads(data)
                    except Exception:
                        continue
                    typ = obj.get("type")
                    if typ == "log":
                        self.line.emit(obj.get("text", ""))
                    elif typ == "done":
                        self.finishedEvt.emit(obj)
                        break
        except Exception as e:
            self.finishedEvt.emit({"type": "done", "status": "error", "message": str(e)})


class SchedulerUI(QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("CEE Scheduler")
        self.resize(980, 720)

        self.offerings_path = ""
        self.input_data_path = ""
        self.current_schedule = []
        self.stream_thread: LogStreamThread | None = None

        root = QVBoxLayout(self)

        # === API Row: URL + Health + Optimize ===
        api_row = QHBoxLayout()
        api_row.addWidget(QLabel("API URL:"))
        self.api_url_input = QLineEdit(API_URL_DEFAULT)
        self.api_url_input.setPlaceholderText("http://localhost:8000")
        api_row.addWidget(self.api_url_input)

        self.health_btn = QPushButton("Health")
        self.health_btn.clicked.connect(self.check_health)
        api_row.addWidget(self.health_btn)

        self.optimize_btn = QPushButton("Optimize")
        self.optimize_btn.clicked.connect(self.run_opt)  # streaming
        api_row.addWidget(self.optimize_btn)

        root.addLayout(api_row)

        # === Files Row ===
        file_row = QHBoxLayout()
        self.offerings_btn = QPushButton("Load Hx_Data_Offerings.xlsx")
        self.offerings_btn.clicked.connect(self.pick_offerings)
        file_row.addWidget(self.offerings_btn)

        self.input_btn = QPushButton("Load input_data.xlsx")
        self.input_btn.clicked.connect(self.pick_input)
        file_row.addWidget(self.input_btn)

        root.addLayout(file_row)

        # === Options Row ===
        opt_row = QHBoxLayout()
        self.semester = QComboBox()
        self.semester.addItems(["Fall", "Spring"])
        opt_row.addWidget(QLabel("Semester:"))
        opt_row.addWidget(self.semester)

        self.return_excel_chk = QCheckBox("Return Excel in response")
        opt_row.addWidget(self.return_excel_chk)
        opt_row.addStretch(1)
        root.addLayout(opt_row)

        # === Table ===
        self.table = QTableWidget(0, 4)
        self.table.setHorizontalHeaderLabels(["Course", "Day", "Start", "End"])
        header = self.table.horizontalHeader()
        header.setSectionResizeMode(QHeaderView.Stretch)
        root.addWidget(self.table)

        # === Optimization Log ===
        log_label = QLabel("Optimization Log")
        root.addWidget(log_label)

        self.log_text = QTextEdit()
        self.log_text.setReadOnly(True)
        self.log_text.setMinimumHeight(200)
        self.log_text.setLineWrapMode(QTextEdit.NoWrap)
        self.log_text.setStyleSheet("font-family: Consolas, 'Courier New', monospace; font-size: 12px;")
        root.addWidget(self.log_text)

        # === Save Excel ===
        self.save_btn = QPushButton("Save table to Excel...")
        self.save_btn.clicked.connect(self.save_excel)
        self.save_btn.setEnabled(False)
        root.addWidget(self.save_btn)

    # ---------- UI Actions ----------

    def pick_offerings(self):
        path, _ = QFileDialog.getOpenFileName(self, "Select Hx_Data_Offerings.xlsx", "", "Excel (*.xlsx)")
        if path:
            self.offerings_path = path
            self.offerings_btn.setText(path)

    def pick_input(self):
        path, _ = QFileDialog.getOpenFileName(self, "Select input_data.xlsx", "", "Excel (*.xlsx)")
        if path:
            self.input_data_path = path
            self.input_btn.setText(path)

    def api_base(self) -> str:
        return (self.api_url_input.text() or API_URL_DEFAULT).strip().rstrip("/")

    def check_health(self):
        url = self.api_base() + "/health"
        try:
            r = requests.get(url, timeout=10)
            r.raise_for_status()
            data = r.json()
            QMessageBox.information(self, "Health", f"{url}\n\n{data}")
        except Exception as e:
            QMessageBox.critical(self, "Health check failed", f"{url}\n\n{e}")

    def run_opt(self):
        if not self.offerings_path or not self.input_data_path:
            QMessageBox.warning(self, "Missing files", "Pick both Hx_Data_Offerings.xlsx and input_data.xlsx")
            return

        def to_container_path(host_path):
            return "/data/" + host_path.split("/")[-1]

        payload = {
            "plan_semester": self.semester.currentText(),
            "offerings_path": to_container_path(self.offerings_path),
            "input_data_path": to_container_path(self.input_data_path),
            # Provide if you want the API to also write a file; else rely on return_excel_bytes
            "output_excel_path": None,
            "return_excel_bytes": self.return_excel_chk.isChecked(),
        }

        # reset UI
        self.log_text.clear()
        self.table.setRowCount(0)
        self.save_btn.setEnabled(False)
        self.optimize_btn.setEnabled(False)
        self.optimize_btn.setText("Optimizing...")

        url = self.api_base() + "/optimize/stream"
        self.stream_thread = LogStreamThread(url, payload)
        self.stream_thread.line.connect(self._append_log_line)
        self.stream_thread.finishedEvt.connect(self._stream_done)
        self.stream_thread.finished.connect(self._stream_cleanup)
        self.stream_thread.start()

    def _append_log_line(self, text: str):
        if not text:
            return
        self.log_text.append(text)
        cursor = self.log_text.textCursor()
        cursor.movePosition(cursor.End)
        self.log_text.setTextCursor(cursor)

    def _stream_done(self, result: dict):
        if result.get("status") != "success":
            msg = result.get("message", "Unknown error")
            QMessageBox.warning(self, "Optimization error", msg)
            return

        rows = result.get("rows")
        if isinstance(rows, list) and rows:
            self.current_schedule = rows
            self.populate_table(rows)
            self.save_btn.setEnabled(True)

        excel_b64 = result.get("excel_base64")
        if excel_b64:
            save_path, _ = QFileDialog.getSaveFileName(self, "Save returned Excel", "Optimal_Solution.xlsx", "Excel (*.xlsx)")
            if save_path:
                try:
                    with open(save_path, "wb") as f:
                        f.write(base64.b64decode(excel_b64))
                    QMessageBox.information(self, "Saved", f"Excel written to:\n{save_path}")
                except Exception as e:
                    QMessageBox.critical(self, "Save failed", str(e))

        QMessageBox.information(self, "Done", "Optimization finished. See the log and table above.")

    def _stream_cleanup(self):
        self.optimize_btn.setEnabled(True)
        self.optimize_btn.setText("Optimize")

    def populate_table(self, rows):
        self.table.setRowCount(len(rows))
        for i, row in enumerate(rows):
            self.table.setItem(i, 0, QTableWidgetItem(row.get("Course", "")))
            self.table.setItem(i, 1, QTableWidgetItem(row.get("Day", "")))
            self.table.setItem(i, 2, QTableWidgetItem(str(row.get("Start", ""))))
            self.table.setItem(i, 3, QTableWidgetItem(str(row.get("End", ""))))

    def save_excel(self):
        if not self.current_schedule:
            return
        import pandas as pd
        save_path, _ = QFileDialog.getSaveFileName(self, "Save table to Excel", "Optimal_Solution.xlsx", "Excel (*.xlsx)")
        if not save_path:
            return
        try:
            df = pd.DataFrame(self.current_schedule)
            with pd.ExcelWriter(save_path, engine="xlsxwriter") as writer:
                df.to_excel(writer, sheet_name="Optimal_Solution", index=False)
            QMessageBox.information(self, "Saved", f"Wrote {len(df)} rows to:\n{save_path}")
        except Exception as e:
            QMessageBox.critical(self, "Save failed", str(e))


def main():
    app = QApplication(sys.argv)
    ui = SchedulerUI()
    ui.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
