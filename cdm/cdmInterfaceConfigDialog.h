#ifndef INTERFACECONFIGDIALOG_H
#define INTERFACECONFIGDIALOG_H

#include <QtSql>
#include <QtGui>
#if QT_VERSION < 0x050000
#else
#include <QtWidgets>
#endif
#include "QsLog.h"
#include "cdmOptions.h"

class InterfaceConfigDialog : public QDialog
{
  
  Q_OBJECT
    
    public:
  
  InterfaceConfigDialog( QWidget *_parentwidget = 0, Qt::WindowFlags fl = Qt::SubWindow,
		    Options *_options = 0);
  virtual ~InterfaceConfigDialog();

 public slots:
 void ok();
 void onMaterialInterfaceChanged(int index);

 signals:

  void handleSelectionChanged(int index);
 
 protected:
  void closeEvent(QCloseEvent *event);
  void showWarning(QString message);

 private:

  void         updateOptions();

  Options      *options;
  QVector<QLineEdit*>   epsilonr;
  QVector<QLineEdit*>   epsiloni;
  QVector<QLineEdit*>   position;
  QVector<QComboBox*>   materiallayer;
  QPushButton  *okButton;
};



#endif
