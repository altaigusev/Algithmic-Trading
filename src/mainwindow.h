#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QComboBox>
#include <QDialog>
#include <QDialogButtonBox>
#include <QHBoxLayout>
#include <QMainWindow>
#include <QMessageBox>
#include <QPushButton>
#include <QRandomGenerator>
#include <QStandardPaths>
#include <QTableWidgetItem>

#include "QtWidgets/QFileDialog"
#include "approximation.h"
#include "interpolation.h"
#include "loader.h"
#include "qcustomplot.h"

QT_BEGIN_NAMESPACE
namespace Ui {
class MainWindow;
}
QT_END_NAMESPACE

class MainWindow : public QMainWindow {
  Q_OBJECT

 public:
  MainWindow(QWidget *parent = nullptr);
  ~MainWindow();

 private slots:
  void on_addGraph_clicked();

  void on_switchInterpolation_clicked();
  void on_clearInterpolationButton_clicked();
  void on_deleteSelectedButton_clicked();
  void on_fileOpenButton_clicked();
  void on_calculateInterpolatedValue_clicked();

  void on_switchApprox_clicked();
  void on_addGraphApprox_clicked();
  void on_clearApproximationButton_clicked();
  void on_clearRowApprox_clicked();
  void on_loadDataset_clicked();
  void on_calculationOnePoint_clicked();
  void on_loadWeights_clicked();
  void on_clearWeights_clicked();
  void on_spinBox_2_valueChanged();
  void on_point_box_valueChanged();

 private:
  enum class InterType { kCubic, kNewton, approx };
  struct GraphLine {
   public:
    QString filename;
    QString filenameApprox;
    QColor linecolor;
    int polynome_pow;
    int xPlus;
    InterType type;
    QString type_string;
  };

  Ui::MainWindow *ui;
  QString filename;
  QString filenameApprox;
  QString filenameApproxWeights;
  bool loadWeights = false;
  int graph_count;
  int graphCountAprox;
  int countPointsApprox;
  double y_max_, y_min_, x_max_, x_min_, y_max_app_, x_max_app_;
  GraphLine AddGraphDataToInterpolationTable();
  GraphLine AddGraphDataToApproxTable();
  void AddStructToTable(GraphLine &line);
  void AddStructToTableApprox(GraphLine line);
  GraphLine GetLineParametersDialog(const QString file);
  GraphLine GetLineParametersDialogApprox(const QString file);
  void AddGraphDataToApproximationTable();
  void CalculateAndPlot(const GraphLine &line);
  void CalculateAndPlotApprox(const GraphLine &line);
  std::vector<double> X_interpol, Y_interpol;
  Interpolation interpol;
  int num_of_points_inter_;
  QColor GetNewRndColor();
  QSet<int> all_colors_;
  QString start_date;
  void ClearCommon();
};
#endif  // MAINWINDOW_H
