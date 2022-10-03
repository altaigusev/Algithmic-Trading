#include "mainwindow.h"

#include "ui_mainwindow.h"

MainWindow::MainWindow(QWidget* parent)
    : QMainWindow(parent), ui(new Ui::MainWindow) {
  ui->setupUi(this);
  graph_count = 1;
  graphCountAprox = 1;
  y_min_ = 0;
  x_min_ = 0;
  y_max_ = 10;
  x_max_ = 100;
  ui->interpolationView->setInteraction(QCP::iRangeZoom, true);
  ui->interpolationView->setInteraction(QCP::iRangeDrag, true);
}

MainWindow::~MainWindow() { delete ui; }

void MainWindow::on_addGraph_clicked() {
  if (ui->tableWidget->rowCount() < 5) {
    if (filename.isEmpty() == false) {
      try {
        GraphLine new_line = AddGraphDataToInterpolationTable();
        CalculateAndPlot(new_line);
      } catch (std::exception& e) {
        QString error = e.what();
        QMessageBox msg;
        msg.setText("Ошибка! " + error);
        msg.exec();
      }
    }
  } else {
    QMessageBox msg;
    msg.setText("Количество графиков не может быть больше 5!");
    msg.exec();
  }
}

void MainWindow::on_switchInterpolation_clicked() {
  ui->stackedWidget->setCurrentIndex(1);
}

MainWindow::GraphLine MainWindow::AddGraphDataToInterpolationTable() {
  GraphLine new_line = GetLineParametersDialog(filename);
  AddStructToTable(new_line);
  return new_line;
}

void MainWindow::AddStructToTable(GraphLine& line) {
  line.linecolor = GetNewRndColor();
  ui->tableWidget->setColumnWidth(0, 100);
  ui->tableWidget->setColumnWidth(2, 100);
  ui->tableWidget->insertRow(ui->tableWidget->rowCount());
  ui->tableWidget->setItem(ui->tableWidget->rowCount() - 1, 0,
                           new QTableWidgetItem(line.filename));
  ui->tableWidget->setItem(ui->tableWidget->rowCount() - 1, 1,
                           new QTableWidgetItem);
  ui->tableWidget->item(ui->tableWidget->rowCount() - 1, 1)
      ->setBackground(line.linecolor);
  if (line.type == InterType::kCubic) {
    ui->tableWidget->setItem(ui->tableWidget->rowCount() - 1, 2,
                             new QTableWidgetItem(line.type_string));
  } else if (line.type == InterType::kNewton) {
    ui->tableWidget->setItem(
        ui->tableWidget->rowCount() - 1, 2,
        new QTableWidgetItem(line.type_string + " со степенью " +
                             QString::number(line.polynome_pow)));
  }
}

MainWindow::GraphLine MainWindow::GetLineParametersDialog(const QString file) {
  GraphLine new_line;
  new_line.filename = file.split('/').last();
  QDialog interpol;
  QHBoxLayout lt(&interpol);
  QLabel lbl("Тип интерполяции", &interpol);
  QComboBox bx(&interpol);
  QPushButton btn("Добавить график", &interpol);
  QLabel lbl2("Степень полинома", &interpol);
  QSpinBox pow_newton;
  pow_newton.setMinimum(3);
  pow_newton.setMaximum(9);
  connect(&btn, SIGNAL(clicked()), &interpol, SLOT(close()));
  connect(&bx, &QComboBox::currentIndexChanged, &lbl2,
          [&lbl2](int param) { lbl2.setHidden(!bool(param)); });
  connect(&bx, &QComboBox::currentIndexChanged, &pow_newton,
          [&pow_newton](int param) { pow_newton.setHidden(!bool(param)); });
  lt.insertWidget(0, &lbl);
  lt.insertWidget(1, &bx);
  lt.insertWidget(2, &lbl2);
  lt.insertWidget(3, &pow_newton);
  lt.insertWidget(4, &btn);
  bx.addItem("Кубическая");
  bx.addItem("Полиномиальная");
  interpol.setModal(true);
  interpol.exec();

  new_line.type_string = bx.currentText();
  if (bx.currentIndex() == 0) {
    new_line.type = InterType::kCubic;
  } else {
    new_line.type = InterType::kNewton;
  }
  new_line.polynome_pow = pow_newton.value();
  return new_line;
}

void MainWindow::AddGraphDataToApproximationTable() {}

void MainWindow::CalculateAndPlot(const GraphLine& line) {
  QDateTime new_date = QDateTime::fromString(start_date, "yyyy-MM-dd");
  QVector<double> result;
  QVector<double> result_x;
  std::vector<double> x_std;
  for (auto xx : X_interpol) {
    x_std.push_back(xx);
  }

  std::vector<double> y_std;
  for (auto yy : Y_interpol) {
    y_std.push_back(yy);
  }
  double iter_step = x_max_ / ui->spinBoxPointsInter->value();
  if (line.type == InterType::kCubic) {
    for (double i = 0; i < x_max_; i += iter_step) {
      result_x.push_back(
          QCPAxisTickerDateTime::dateTimeToKey(new_date.addSecs(i * 86400)));
      result.push_back(interpol.GetValueCubicSpline(x_std, y_std, i));
    }
    result_x.push_back(QCPAxisTickerDateTime::dateTimeToKey(
        new_date.addSecs(*(--x_std.end()) * 86400)));
    result.push_back(
        interpol.GetValueCubicSpline(x_std, y_std, *(--x_std.end())));
  } else if (line.type == InterType::kNewton) {
    int divider = x_max_ / line.polynome_pow;
    std::vector<double> x_new;
    std::vector<double> y_new;
    for (size_t k = 0; k < x_std.size(); ++k) {
      if (k % divider == 0) {
        x_new.push_back(x_std[k]);
        y_new.push_back(y_std[k]);
      }
    }
    x_new.push_back(*(--x_std.end()));
    y_new.push_back(*(--y_std.end()));
    for (double i = 0.0; i <= x_max_; i += iter_step) {
      result_x.push_back(
          QCPAxisTickerDateTime::dateTimeToKey(new_date.addSecs(i * 86400)));
      double new_val = interpol.GetValueNewtonPolynome(x_new, y_new, i);
      if (new_val > y_max_ && new_val < 1000) y_max_ = new_val;
      result.push_back(new_val);
    }
    result_x.push_back(QCPAxisTickerDateTime::dateTimeToKey(
        new_date.addSecs(*(--x_std.end()) * 86400)));
    result.push_back(
        interpol.GetValueNewtonPolynome(x_new, y_new, *(--x_std.end())));
  }
  ui->interpolationView->addGraph();
  ui->interpolationView->graph(graph_count)->setData(result_x, result);
  ui->interpolationView->graph(graph_count)->setPen(QPen(line.linecolor, 2));
  ++graph_count;
  ui->interpolationView->yAxis->setRange(y_min_, y_max_);
  ui->interpolationView->rescaleAxes();
  ui->interpolationView->replot();
}

QColor MainWindow::GetNewRndColor() {
  return QColor::fromRgb(QRandomGenerator::global()->generate());
  ;
}

void MainWindow::ClearCommon() {
  ui->interpolationView->clearGraphs();
  ui->interpolationView->replot();
  for (int i = ui->tableWidget->rowCount(); i >= 0; --i) {
    ui->tableWidget->removeRow(i);
  }
  graph_count = 1;
  x_max_ = 100;
  y_max_ = 10;
}

void MainWindow::on_clearInterpolationButton_clicked() {
  ClearCommon();
  filename.clear();
}

void MainWindow::on_deleteSelectedButton_clicked() {
  QModelIndex currentIndex = ui->tableWidget->currentIndex();
  if (currentIndex.row() != -1) {
    int line_number = currentIndex.row();
    ui->tableWidget->removeRow(line_number);
    if (ui->interpolationView->removeGraph(line_number + 1) == true) {
      graph_count--;
    }
    ui->interpolationView->replot();
  }
}

void MainWindow::on_fileOpenButton_clicked() {
  X_interpol.clear();
  Y_interpol.clear();

  QFileDialog dlg;
  dlg.setWindowTitle("Open csv");
  QString searchpath;
  if (filename.isEmpty() == false) {
    searchpath = filename;
  } else {
    searchpath =
        QStandardPaths::writableLocation(QStandardPaths::DesktopLocation);
  }
  filename.clear();
  filename = dlg.getOpenFileName(this, "open file", searchpath, "*.csv");
  if (filename.isEmpty() == false && filename != searchpath) {
    ClearCommon();
    try {
      auto vec = s21::Loader::getDataNorm(filename.toStdString());
      start_date = QString::fromStdString(s21::Loader::getFirstDate(
          filename.toStdString()));  // добавить начальную дату
      QDate new_date = QDate::fromString(start_date, "yyyy-MM-dd");
      ui->dateTimeEdit->setMinimumDate(new_date);
      QMap<double, double> XY;
      double firstVal = vec.at(0).second;
      for (auto val : vec) {
        XY.insert(val.second, val.first);
      }
      for (auto it = XY.begin(); it != XY.end(); ++it) {
        if ((it.key() - firstVal) > x_max_) x_max_ = it.key() - firstVal;
        if (it.value() > y_max_) y_max_ = it.value();
        X_interpol.push_back(it.key() - firstVal);
        Y_interpol.push_back(it.value());
      }
      ui->spinBoxPointsInter->setMinimum(int(XY.size()));
      QVector<double> tmpX, tmpY;
      for (auto x_val : X_interpol) {
        tmpX.push_back(
            QCPAxisTickerDateTime::dateTimeToKey(new_date.addDays(x_val)));
      }
      for (auto y_val : Y_interpol) {
        tmpY.push_back(y_val);
      }
      ui->interpolationView->xAxis->setLabel("Дата");
      ui->interpolationView->yAxis->setLabel("Цена");
      ui->interpolationView->addGraph();
      ui->interpolationView->graph(0)->setData(tmpX, tmpY);
      ui->interpolationView->graph(0)->setScatterStyle(
          QCPScatterStyle::ssCircle);
      ui->interpolationView->graph(0)->setPen(QPen(Qt::black, 1));
      QSharedPointer<QCPAxisTickerDateTime> dateTimeTicker(
          new QCPAxisTickerDateTime);
      ui->interpolationView->xAxis->setTicker(dateTimeTicker);
      dateTimeTicker->setDateTimeFormat("d. MMM\nyyyy");
      ui->interpolationView->xAxis->setRange(
          QCPAxisTickerDateTime::dateTimeToKey(new_date),
          QCPAxisTickerDateTime::dateTimeToKey(new_date.addDays(x_max_)));
      ui->interpolationView->yAxis->setRange(y_min_, y_max_);
      ui->interpolationView->replot();
    } catch (std::exception& e) {
      QString error = QString::fromStdString(e.what());
      QMessageBox msg;
      msg.setText("ERROR: " + error);
      msg.exec();
    }
  }
}

void MainWindow::on_calculateInterpolatedValue_clicked() {
  if (X_interpol.empty() == false) {
    double date_time =
        double(ui->dateTimeEdit->dateTime().toSecsSinceEpoch() -
               QDateTime(QDate::fromString("2021-03-22", "yyyy-MM-dd"), QTime())
                   .toSecsSinceEpoch()) /
        double(86400);
    double cubic_value =
        interpol.GetValueCubicSpline(X_interpol, Y_interpol, date_time);
    QString cubic_price = QString::number(cubic_value);
    ui->cubicLabel->setText(cubic_price);

    int polynome_pow = ui->polynomialSpinBox->value();
    int divider = x_max_ / polynome_pow;
    std::vector<double> x_new;
    std::vector<double> y_new;
    for (size_t k = 0; k < X_interpol.size(); ++k) {
      if (k % divider == 0) {
        x_new.push_back(X_interpol[k]);
        y_new.push_back(Y_interpol[k]);
      }
    }
    x_new.push_back(*(--X_interpol.end()));
    y_new.push_back(*(--Y_interpol.end()));
    double newton_value =
        interpol.GetValueNewtonPolynome(x_new, y_new, date_time);
    QString newton_price = QString::number(newton_value);
    ui->polynomialLabel->setText(newton_price);
  }
}

//*******************************************APPROXIMATION*******************************************

void MainWindow::on_switchApprox_clicked() {
  ui->stackedWidget->setCurrentIndex(0);
}

void MainWindow::on_addGraphApprox_clicked() {
  if (ui->tableWidgetApprox->rowCount() < 5) {
    if (filenameApprox.isEmpty() == false) {
      try {
        GraphLine new_line = AddGraphDataToApproxTable();
        CalculateAndPlotApprox(new_line);
      } catch (std::exception& e) {
        QString error = QString::fromStdString(e.what());
        QMessageBox msg;
        msg.setText("ERROR: " + error);
        msg.exec();
      }
    }
  } else {
    QMessageBox msg;
    msg.setText("Количество графиков не может быть больше 5!");
    msg.exec();
  }
}

MainWindow::GraphLine MainWindow::AddGraphDataToApproxTable() {
  GraphLine new_line = GetLineParametersDialogApprox(filenameApprox);
  new_line.linecolor = GetNewRndColor();
  AddStructToTableApprox(new_line);
  return new_line;
}

void MainWindow::AddStructToTableApprox(GraphLine line) {
  ui->tableWidgetApprox->setColumnWidth(0, 100);
  ui->tableWidgetApprox->setColumnWidth(2, 100);
  ui->tableWidgetApprox->insertRow(ui->tableWidgetApprox->rowCount());
  ui->tableWidgetApprox->setItem(ui->tableWidgetApprox->rowCount() - 1, 0,
                                 new QTableWidgetItem(line.filenameApprox));
  ui->tableWidgetApprox->setItem(ui->tableWidgetApprox->rowCount() - 1, 1,
                                 new QTableWidgetItem);
  ui->tableWidgetApprox->item(ui->tableWidgetApprox->rowCount() - 1, 1)
      ->setBackground(line.linecolor);

  if (line.type == InterType::approx) {
    QString w;
    if (loadWeights == true) {
      w = " С весами ";
    } else {
      w = " Без весов ";
    }
    ui->tableWidgetApprox->setItem(
        ui->tableWidgetApprox->rowCount() - 1, 2,
        new QTableWidgetItem("полином со степенью " +
                             QString::number(line.polynome_pow) + w));
  }
}
MainWindow::GraphLine MainWindow::GetLineParametersDialogApprox(
    const QString file) {
  GraphLine new_line;
  new_line.type = InterType::approx;
  new_line.filenameApprox = file.split('/').last();
  QDialog approx;
  QHBoxLayout lt(&approx);
  QLabel lbl2("Степень полинома", &approx);
  QPushButton btn("Добавить график", &approx);
  QSpinBox powP;
  powP.setMinimum(1);
  powP.setMaximum(6);
  connect(&btn, SIGNAL(clicked()), &approx, SLOT(close()));
  lt.insertWidget(0, &lbl2);
  lt.insertWidget(1, &powP);
  lt.insertWidget(2, &btn);
  approx.setModal(true);
  approx.exec();
  new_line.polynome_pow = powP.value();
  new_line.xPlus = ui->spinBox_2->value();
  return new_line;
}

void MainWindow::CalculateAndPlotApprox(const GraphLine& line) {
  std::vector<std::pair<double, double>> data =
      s21::Loader::getDataNorm(filenameApprox.toStdString());
  start_date = QString::fromStdString(
      s21::Loader::getFirstDate(filenameApprox.toStdString()));
  QDate new_date = QDate::fromString(start_date, "yyyy-MM-dd");

  QVector<double> X, newX, xPlus, yPlus, newY;

  double step = data[data.size() - 1].second / countPointsApprox;
  for (double i = 1.0; i < data[data.size() - 1].second; i += step) {
    newX.push_back(QCPAxisTickerDateTime::dateTimeToKey(new_date) + i * 86400);
  }
  y_max_app_ = 0.0;
  x_max_app_ = 0.0;
  int powP = line.polynome_pow;
  for (auto val : data) {
    X.push_back(val.second);
  }
  std::vector<double> tmp, weights;
  std::vector<std::pair<double, double>> data2;
  if (filenameApproxWeights.isEmpty() == false && loadWeights == true) {
    weights = s21::Loader::loadWeights(filenameApproxWeights.toStdString());
    tmp = s21::Approximation::calcAllNewYWithW(weights, powP, data,
                                               countPointsApprox);
    data2 = s21::Approximation::getNewData(weights, powP, data);
  } else {
    tmp = s21::Approximation::calcALLNewY(powP, data, countPointsApprox);
  }
  if (line.xPlus > 0) {
    for (int i = 1; i < line.xPlus; ++i) {
      double t = X[X.size() - 1] + i;
      xPlus.push_back(QCPAxisTickerDateTime::dateTimeToKey(new_date) +
                      t * 86400);
      if (filenameApproxWeights.isEmpty() == false && loadWeights == true) {
        yPlus.push_back(s21::Approximation::getResult(powP, data2, t));
      } else {
        yPlus.push_back(s21::Approximation::getResult(powP, data, t));
      }
    }
  }
  for (auto x : tmp) {
    newY.push_back(x);
  }
  ui->approxView->addGraph();
  ui->approxView->graph(graphCountAprox)->setData(newX, newY);
  ui->approxView->graph(graphCountAprox)->addData(xPlus, yPlus);
  ui->approxView->graph(graphCountAprox)->setPen(QPen(line.linecolor, 3));
  ++graphCountAprox;
  ui->approxView->replot();
}

void MainWindow::on_loadDataset_clicked() {
  on_clearApproximationButton_clicked();
  QFileDialog dlg;
  dlg.setWindowTitle("Open csv");
  QString searchpath;
  if (filenameApprox.isEmpty() == false) {
    searchpath = filenameApprox;
  } else {
    searchpath =
        QStandardPaths::writableLocation(QStandardPaths::DesktopLocation);
  }
  try {
    filenameApprox =
        dlg.getOpenFileName(this, "open file", searchpath, "*.csv");
    if (filenameApprox.isEmpty() == false) {
      std::vector<std::pair<double, double>> data =
          s21::Loader::getDataNorm(filenameApprox.toStdString());
      start_date = QString::fromStdString(
          s21::Loader::getFirstDate(filenameApprox.toStdString()));
      QDate new_date = QDate::fromString(start_date, "yyyy-MM-dd");

      countPointsApprox = data.size();
      ui->point_box->setMinimum(countPointsApprox);
      QVector<double> x, y;
      for (auto val : data) {
        x.push_back(QCPAxisTickerDateTime::dateTimeToKey(new_date) +
                    val.second * 86400);
        y.push_back(val.first);
        y_max_app_ = val.first >= y_max_app_ ? val.first : y_max_app_;
        x_max_app_ = (double)val.second >= x_max_app_ ? val.second : x_max_app_;
      }
      ui->approxView->addGraph();
      ui->approxView->graph(0)->setData(x, y);
      ui->approxView->graph(0)->setScatterStyle(QCPScatterStyle::ssCircle);
      ui->approxView->graph(0)->setPen(QPen(Qt::black, 1));
      ui->approxView->xAxis->setLabel("Дата");
      ui->approxView->yAxis->setLabel("Цена");

      QSharedPointer<QCPAxisTickerDateTime> dateTimeTicker(
          new QCPAxisTickerDateTime);
      ui->approxView->xAxis->setTicker(dateTimeTicker);
      dateTimeTicker->setDateTimeFormat("d. MMM\nyyyy");
      ui->approxView->xAxis->setRange(
          QCPAxisTickerDateTime::dateTimeToKey(new_date),
          QCPAxisTickerDateTime::dateTimeToKey(
              new_date.addDays(x_max_app_ * 1.2)));

      //      ui->approxView->xAxis->setRange(0, x_max_app_ * 1.2);
      ui->approxView->yAxis->setRange(0, y_max_app_ * 1.5);
      ui->approxView->setInteraction(QCP::iRangeZoom, true);
      ui->approxView->setInteraction(QCP::iRangeDrag, true);
      ui->approxView->replot();
    }
  } catch (std::exception& e) {
    QString error = QString::fromStdString(e.what());
    QMessageBox msg;
    msg.setText("ERROR: " + error);
    msg.exec();
  }
}

void MainWindow::on_clearApproximationButton_clicked() {
  ui->approxView->clearGraphs();
  ui->approxView->replot();
  for (int i = ui->tableWidgetApprox->rowCount(); i >= 0; --i) {
    ui->tableWidgetApprox->removeRow(i);
  }
  graphCountAprox = 1;
  loadWeights = false;
  filename.clear();
  filenameApprox.clear();
  filenameApproxWeights.clear();
}

void MainWindow::on_clearRowApprox_clicked() {
  QModelIndex currentIndex = ui->tableWidgetApprox->currentIndex();
  if (currentIndex.row() != -1) {
    int line_number = currentIndex.row();
    ui->tableWidgetApprox->removeRow(line_number);
    if (ui->approxView->removeGraph(line_number + 1) == true) {
      graphCountAprox--;
    }
    ui->approxView->replot();
  }
}

void MainWindow::on_calculationOnePoint_clicked() {
  if (filenameApprox.isEmpty() == false) {
    try {
      std::vector<std::pair<double, double>> data =
          s21::Loader::getDataNorm(filenameApprox.toStdString());
      int pow = ui->spinBox_3->value();
      double date_time =
          double(QDateTime::fromString(ui->dateTimeEdit_2->dateTime().toString(
                                           "yyyy-MM-dd-hh-mm"),
                                       "yyyy-MM-dd-hh-mm")
                     .toSecsSinceEpoch() -
                 QDateTime::fromString(
                     QString::fromStdString(s21::Loader::getFirstDate(
                         filenameApprox.toStdString())),
                     "yyyy-MM-dd")
                     .toSecsSinceEpoch()) /
          double(86400);
      if (date_time < 0) {
        QString error = "ERROR: The given date is less than the given one";
        QMessageBox msg;
        msg.setText(error);
        msg.exec();
      }

      double rez = s21::Approximation::getResult(pow, data, date_time);
      ui->label_14->setText(QString::number(rez));
    } catch (std::exception& e) {
      QString error = QString::fromStdString(e.what());
      QMessageBox msg;
      msg.setText("ERROR: " + error);
      msg.exec();
    }
  } else {
    QMessageBox msg;
    msg.setText("Датасет не загружен!");
    msg.exec();
  }
}

void MainWindow::on_loadWeights_clicked() {
  try {
    QFileDialog dlg;
    dlg.setWindowTitle("Open csv");
    QString searchpath;
    if (filenameApproxWeights.isEmpty() == false) {
      searchpath = filenameApproxWeights;
    } else {
      searchpath =
          QStandardPaths::writableLocation(QStandardPaths::DesktopLocation);
    }
    filenameApproxWeights =
        dlg.getOpenFileName(this, "open file", searchpath, "*.csv");
    if (filenameApproxWeights.isEmpty() == false) loadWeights = true;

  } catch (std::exception& e) {
    QString error = QString::fromStdString(e.what());
    QMessageBox msg;
    msg.setText("ERROR: " + error);
    msg.exec();
  }
}

void MainWindow::on_clearWeights_clicked() {
  loadWeights = false;
  filenameApproxWeights.clear();
}

void MainWindow::on_spinBox_2_valueChanged() {
  while (graphCountAprox > 1) {
    ui->tableWidgetApprox->removeRow(ui->tableWidgetApprox->rowCount() - 1);
    ui->approxView->removeGraph(ui->tableWidgetApprox->rowCount() + 1);
    --graphCountAprox;
  }
  ui->approxView->replot();
}

void MainWindow::on_point_box_valueChanged() {
  countPointsApprox = ui->point_box->value();
}
