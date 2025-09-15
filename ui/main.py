import sys
import base64
import requests
from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QApplication, QWidget, QVBoxLayout, QHBoxLayout, QLabel, QPushButton,
    QFileDialog, QComboBox, QTableWidget, QTableWidgetItem, QMessageBox,
    QCheckBox, QLineEdit, QHeaderView
)

API_URL_DEFAULT = "http://localhost:8000"

class SchedulerUI(QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("CEE Scheduler")
        self.resize(980, 640)

        self.offerings_path = ""
        self.input_data_path = ""
        self.current_schedule = []

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
        self.optimize_btn.clicked.connect(self.run_opt)
        api_row.addWidget(self.optimize_btn)

        root.addLayout(api_row)

        # === Files Row: pick Excel files ===
        file_row = QHBoxLayout()
        self.offerings_btn = QPushButton("Pick Hx_Data_Offerings.xlsx")
        self.offerings_btn.clicked.connect(self.pick_offerings)
        file_row.addWidget(self.offerings_btn)

        self.input_btn = QPushButton("Pick input_data.xlsx")
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
        header.setSectionResizeMode(QHeaderView.Stretch)  # equal width columns
        root.addWidget(self.table)

        # === Save Excel button ===
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

        # Map host path to container path (assumes ./data mounted to /data in Docker)
        def to_container_path(host_path):
            return "/data/" + host_path.split("/")[-1]

        payload = {
            "plan_semester": self.semester.currentText(),
            "offerings_path": to_container_path(self.offerings_path),
            "input_data_path": to_container_path(self.input_data_path),
            "output_excel_path": "/data/{}_Schedule.xlsx".format(self.semester.currentText()),
            "return_excel_bytes": self.return_excel_chk.isChecked(),
        }

        url = self.api_base() + "/optimize"
        try:
            r = requests.post(url, json=payload, timeout=600)
            r.raise_for_status()
            data = r.json()
        except Exception as e:
            QMessageBox.critical(self, "Request failed", f"{url}\n\n{e}")
            return

        if data.get("status") != "success":
            QMessageBox.warning(self, "Optimization error", data.get("message", "Unknown error"))
            return

        schedule = data.get("schedule", [])
        self.current_schedule = schedule
        self.populate_table(schedule)
        self.save_btn.setEnabled(True)

        if data.get("excel_base64"):
            save_path, _ = QFileDialog.getSaveFileName(self, "Save returned Excel", "Optimal_Solution.xlsx", "Excel (*.xlsx)")
            if save_path:
                with open(save_path, "wb") as f:
                    f.write(base64.b64decode(data["excel_base64"]))
                QMessageBox.information(self, "Saved", f"Excel written to:\n{save_path}")

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
