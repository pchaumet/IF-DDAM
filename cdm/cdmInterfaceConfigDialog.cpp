#include "cdmInterfaceConfigDialog.h"

InterfaceConfigDialog::InterfaceConfigDialog(QWidget *parent,Qt::WindowFlags fl,
				   Options *_options)
: QDialog(parent,fl)
{
   options = _options;
   okButton = new QPushButton("Ok");
   connect(okButton, SIGNAL(clicked()), this, SLOT(ok()));
   QSignalMapper *signalMapper = new QSignalMapper(this);
   
   QFormLayout *layout = new QFormLayout(this);
  
   for (int i = 0 ; i < options->getInterfaceNumber()+1; i++) {
     epsilonr.push_back(new QLineEdit(QString::number(real(options->getEpsilonlayer().at(i)))));
     epsiloni.push_back(new QLineEdit(QString::number(imag(options->getEpsilonlayer().at(i)))));

      QComboBox *materiallayerBox = new QComboBox(this);
      connect( materiallayerBox, SIGNAL(currentIndexChanged(int)), signalMapper, SLOT(map()));
      signalMapper->setMapping(materiallayerBox, i);
      materiallayerBox->addItems(options->materialList);
      materiallayerBox->setCurrentIndex(0);
      materiallayer.push_back(materiallayerBox);


     
     QBoxLayout   *epsilonlayout = new QBoxLayout(QBoxLayout::LeftToRight);
     epsilonlayout->addWidget(epsilonr.at(i));
     epsilonlayout->addWidget(epsiloni.at(i));
     epsilonlayout->addWidget(materiallayer.at(i));
     if (i ==0) {
       layout->addRow("Epsilon substrat:",epsilonlayout);
     }
     else if (i == options->getInterfaceNumber()) {
       layout->addRow("Epsilon superstrat:",epsilonlayout);
     }
     else {
       layout->addRow("Epsilon:",epsilonlayout);
     }

     if ( i < options->getInterfaceNumber() ) {
       position.push_back(new QLineEdit(QString::number(options->getPositioninterface().at(i))));
       QString title = "Interface " + QString::number(i+1) + " ---- Position (nm):";
       layout->addRow(title,position.at(i));
     }
   }
   layout->addRow("",okButton);
   this->setWindowTitle(tr("Multilayer properties"));


   connect(signalMapper, SIGNAL(mapped(int)), this, SIGNAL(handleSelectionChanged(int)));
   connect(this, SIGNAL(handleSelectionChanged(int)), this, SLOT(onMaterialInterfaceChanged(int)));
     for (int i = 0 ; i < options->getInterfaceNumber()+1; i++) {
        QComboBox *materiallayerBox = materiallayer.at(i);
        QLOG_DEBUG() << "InterfaceEpsilonConfigDialog::InterfaceEpsilonConfigDialog " << i+1 
		    << " material layer: " <<  options->getMaterialInterface().at(i);
        materiallayerBox->setCurrentIndex(materiallayerBox->findText(options->getMaterialInterface().at(i)));
     }
   
}



InterfaceConfigDialog::~InterfaceConfigDialog()
{
  QLOG_DEBUG ( ) << "Deleting InterfaceConfigDialog";
}

void 
InterfaceConfigDialog::onMaterialInterfaceChanged(int index) {
 QLOG_DEBUG() << "InterfaceEpsilonConfigDialog::onMaterialInterfaceChanged";
 QComboBox *materiallayerBox = materiallayer.at(index);
 if ( materiallayerBox->currentText() != "xx" )  {
   QLineEdit *epsilonrEdit = epsilonr.at(index);
   QLineEdit *epsiloniEdit = epsiloni.at(index);
   epsilonrEdit->setEnabled( false );
   epsiloniEdit->setEnabled( false );
 }
 else {
   QLineEdit *epsilonrEdit = epsilonr.at(index);
   QLineEdit *epsiloniEdit = epsiloni.at(index);
   epsilonrEdit->setEnabled( true );
   epsiloniEdit->setEnabled( true );
 }
}



void
InterfaceConfigDialog::ok()
{
  this->updateOptions();
  this->close();
}
void
InterfaceConfigDialog::updateOptions()
{
   QVector<dcmplx> epsilonlayer;
   QVector<double> positioninterface;
   QVector<QString> _materiallayer;
   
   for (int i = 0 ; i < MAX_LAYER_NUMBER; i++) {
     if (i < options->getInterfaceNumber() + 1 ) {
       epsilonlayer.push_back(dcmplx(epsilonr.at(i)->text().toDouble(),epsiloni.at(i)->text().toDouble()));
       _materiallayer.push_back(materiallayer.at(i)->currentText());
     }
     else {
       epsilonlayer.push_back(dcmplx(1,0));
     }
     if ( i < options->getInterfaceNumber() )
       positioninterface.push_back(position.at(i)->text().toDouble());
     else if ( i == options->getInterfaceNumber() )
       positioninterface.push_back(position.at(i-1)->text().toDouble());
     else 
       positioninterface.push_back(0);
   }
   options->setEpsilonlayer(epsilonlayer);
   options->setPositioninterface(positioninterface);  
   options->setMaterialInterface(_materiallayer); 
}
void 
InterfaceConfigDialog::closeEvent(QCloseEvent* event)
{
  event->accept();  
  QLOG_DEBUG ( ) << "Closing InterfaceConfigDialog";
  this->close();
}
void 
InterfaceConfigDialog::showWarning(QString message) {
  QMessageBox::warning(this, "Error:", message);
}
