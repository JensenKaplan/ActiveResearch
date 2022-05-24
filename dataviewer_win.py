#!/usr/bin/python3
from PyQt5 import QtGui, QtCore, QtWidgets
import os, sys, io
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from misc_lib import *

import matplotlib.cm as cm
from matplotlib import container

#import scipy.signal as sig

# import matplotlib
# matplotlib.rcParams['mathtext.fontset'] = 'stixsans'
# matplotlib.rcParams['mathtext.it'] = 'Arial:italic'
# matplotlib.rcParams['mathtext.rm'] = 'Arial'
# matplotlib.rcParams['mathtext.bf'] = 'Arial:bold'
# matplotlib.rcParams['mathtext.default'] = 'it'
# matplotlib.rcParams['font.family'] = 'Arial'

class PPMSDataViewer(QtWidgets.QWidget):

  def __init__(self):
    super(PPMSDataViewer, self).__init__()
    self.initUI()
    self.dfs = []
    self.processed = False
    self.plotted = False

  def initUI(self):
    self.setGeometry(600, 300, 750, 800)
    self.center()
    self.setWindowTitle("PPMS Data Viewr v0.2.0")

    #Grid Layout
    grid = QtWidgets.QGridLayout()
    self.setLayout(grid)
    
    #Canvas and Toolbar
    self.figure = plt.figure(figsize=(10, 8))
    self.figure.set_tight_layout(True)
    self.canvas = FigureCanvas(self.figure)
    self.canvas.get_default_filename = self.getFileName
    self.toolbar = NavigationToolbar(self.canvas, self)
    grid.addWidget(self.canvas, 1, 0, 1, 5)
    grid.addWidget(self.toolbar, 0, 0, 1, 5)
    self.canvas.draw()

    #Add group boxes
    groupbox1 = QtWidgets.QGroupBox("", self)
    groupbox1.resize(groupbox1.sizeHint())
    grid.addWidget(groupbox1, 2, 0, 1, 2)

    self.title_input = QtWidgets.QLineEdit(self)
    self.title_input.resize(self.title_input.sizeHint())
    self.compound_input = QtWidgets.QLineEdit(self)
    self.compound_input.resize(self.compound_input.sizeHint())
    self.checkbox_grid = QtWidgets.QCheckBox('Grid On', self)
    self.checkbox_grid.resize(self.checkbox_grid.sizeHint())
    self.checkbox_xlog = QtWidgets.QCheckBox('XLogScale', self)
    self.checkbox_xlog.resize(self.checkbox_xlog.sizeHint())
    self.checkbox_ylog = QtWidgets.QCheckBox('YLogScale', self)
    self.checkbox_ylog.resize(self.checkbox_ylog.sizeHint())
    self.checkbox_legendoff = QtWidgets.QCheckBox('Legend Off', self)
    self.checkbox_legendoff.resize(self.checkbox_legendoff.sizeHint())
    #Combo box to select data mode
    self.combobox_mode = QtWidgets.QComboBox(self)
    self.combobox_mode.resize(self.combobox_mode.sizeHint())
    self.combobox_mode.addItem("HC Cp vs T")
    self.combobox_mode.addItem("HC Cp/T vs T")
    self.combobox_mode.addItem("VSM: M vs H")
    self.combobox_mode.addItem("VSM: M vs T")
    self.combobox_mode.addItem("VSM: Chi vs T")
    self.combobox_mode.addItem("VSM: 1/Chi vs T")

    button_plot = QtWidgets.QPushButton('Plot', self)
    button_plot.resize(button_plot.sizeHint())
    button_plot.clicked.connect(self.plot)

    formlayout_mode = QtWidgets.QFormLayout()
    formlayout_mode.addRow("Data Mode: ", self.combobox_mode)
    formlayout_mode.addRow("Compound: ", self.compound_input)
    formlayout_mode.addRow("Title:", self.title_input)
    formlayout_mode.addWidget(self.checkbox_grid)
    formlayout_mode.addWidget(self.checkbox_xlog)
    formlayout_mode.addWidget(self.checkbox_ylog)
    formlayout_mode.addWidget(self.checkbox_legendoff)
    formlayout_mode.addWidget(button_plot)
    groupbox1.setLayout(formlayout_mode)

    groupbox2 = QtWidgets.QGroupBox("", self)
    groupbox2.resize(groupbox2.sizeHint())
    grid.addWidget(groupbox2, 2, 2, 1, 2)

    self.totmass_input = QtWidgets.QLineEdit(self)
    self.totmass_input.resize(self.totmass_input.sizeHint())
    self.masserr_input = QtWidgets.QLineEdit(self)
    self.masserr_input.resize(self.masserr_input.sizeHint())
    self.molmass_input = QtWidgets.QLineEdit(self)
    self.molmass_input.resize(self.totmass_input.sizeHint())
    self.silver_input = QtWidgets.QLineEdit(self)
    self.silver_input.resize(self.silver_input.sizeHint())
    self.field_input = QtWidgets.QLineEdit(self)
    self.field_input.resize(self.field_input.sizeHint())
    #self.checkbox_silver = QtWidgets.QCheckBox('Silver', self)
    #self.checkbox_silver.resize(self.checkbox_silver.sizeHint())
    #self.filetable = QtWidgets.QTableWidget(self)
    #self.filetable.setColumnCount(1)
    #self.filetable.resize(self.filetable.sizeHint())
    #self.filetable.resizeRowsToContents()
    #self.filetable.resizeColumnsToContents()

    formlayout_pars = QtWidgets.QFormLayout()
    formlayout_pars.addRow("Total Mass (mg): ", self.totmass_input)
    formlayout_pars.addRow("Molar Mass (g/mol):", self.molmass_input)
    formlayout_pars.addRow("Total Mass Error:", self.masserr_input)
    #formlayout_pars.addWidget(self.checkbox_silver)
    formlayout_pars.addRow("Silver Ratio (HC only):", self.silver_input)
    formlayout_pars.addRow("Field (Oe) (Chi only):", self.field_input)
    #formlayout_pars.addRow(self.filetable)
    #formlayout_pars.setAlignment(QtCore.Qt.AlignLeft)
    groupbox2.setLayout(formlayout_pars)

    groupbox3 = QtWidgets.QGroupBox("", self)
    groupbox3.resize(groupbox3.sizeHint())
    grid.addWidget(groupbox3, 2, 4, 1, 1)

    #Load CSV Button
    button_importfiles = QtWidgets.QPushButton('Load Data Files', self)
    button_importfiles.resize(button_importfiles.sizeHint())
    button_importfiles.clicked.connect(self.getCSV)

    #Add CSV Button
    button_addfiles = QtWidgets.QPushButton('Add Data Files', self)
    button_addfiles.resize(button_addfiles.sizeHint())
    button_addfiles.clicked.connect(self.addCSV)

    #Plot Button
    button_clearplot = QtWidgets.QPushButton('Clear Plot', self)
    button_clearplot.resize(button_clearplot.sizeHint())
    button_clearplot.clicked.connect(self.clearPlot)

    #Update Parameters Button
    button_updatepars = QtWidgets.QPushButton('Update Parameters', self)
    button_updatepars.resize(button_updatepars.sizeHint())
    button_updatepars.clicked.connect(self.updateParameters)

    #Export Plotted Data
    button_exportdata = QtWidgets.QPushButton('ExportPlotData', self)
    button_exportdata.resize(button_updatepars.sizeHint())
    button_exportdata.clicked.connect(self.exportPlotData)

    #Clear Files Button
    button_clearfiles = QtWidgets.QPushButton('Clear Files', self)
    button_clearfiles.resize(button_clearfiles.sizeHint())
    button_clearfiles.clicked.connect(self.clearFiles)

    #Exit Program Button
    button_exit = QtWidgets.QPushButton('Exit', self)
    button_exit.resize(button_exit.sizeHint())
    button_exit.clicked.connect(self.quit)

    vbox3 = QtWidgets.QVBoxLayout()
    #vbox3.addStretch(1)
    vbox3.addWidget(button_importfiles)
    vbox3.addWidget(button_addfiles)
    vbox3.addWidget(button_clearfiles)
    vbox3.addWidget(button_exportdata)
    vbox3.addWidget(button_clearplot)
    vbox3.addWidget(button_updatepars)
    vbox3.addWidget(button_exit)
    groupbox3.setAlignment(QtCore.Qt.AlignTop)
    vbox3.setAlignment(QtCore.Qt.AlignTop)
    groupbox3.setLayout(vbox3)

    #set default path as home
    self.filepath = "C:\Data"
    self.savepath = "C:\Data\Exports"

    self.show()

  def center(self):
    qr = self.frameGeometry()
    cp = QtWidgets.QDesktopWidget().availableGeometry().center()
    qr.moveCenter(cp)
    self.move(qr.topLeft())

  def getCSV(self):
    filelist = QtWidgets.QFileDialog.getOpenFileNames(self,
                                              'Select one or more files',
                                              self.filepath,
                                              '(*.dat *.DAT)')
    if filelist:
      files = [str(f) for f in filelist[0]]
      try:
        self.filepath = QtCore.QFileInfo(filelist[0][0]).path()
      except:
        self.filepath = "C:\Data"
      self.readCSV(files)

  def addCSV(self):
    filelist = QtWidgets.QFileDialog.getOpenFileNames(self,
                                              'Select one or more files',
                                              self.filepath,
                                              '(*.dat *.DAT)')
    if filelist:
      self.filepath = QtCore.QFileInfo(filelist[0][0]).path()
      tmp_files = [str(f) for f in filelist[0]]
      files = self.files + tmp_files
    self.processed = False
    self.readCSV(files)

  def readCSV(self, files):
    self.files = []
    self.dfs = []
    for n, file in enumerate(files):
      print(file)
      tmp = io.open(file, 'r', encoding = 'latin1')
      nlines = 50
      while nlines > 0:
        tmp_str = tmp.readline()
        if tmp_str.find("Heat") != -1:
          num_rows_skip = 13
          self.combobox_mode.setCurrentIndex(0)
          tmp.close()
          break
        if tmp_str.find("VSM") != -1:
          num_rows_skip = 30
          self.combobox_mode.setCurrentIndex(2)
          tmp.close()
          break
        nlines -= 1
      if nlines > 0:
        self.dfs.append(pd.read_csv(file, skiprows = num_rows_skip,
                        encoding = 'latin1', sep = 'delimiter',
                        delimiter = ',', error_bad_lines = False))
        self.files.append(file)
    self.labels = [s.split('.')[0].split('/')[-1] for s in self.files]
    mode = str(self.combobox_mode.currentText())
    if mode[0:2] == "HC" and self.processed == False:
      self.processed = True
      self.processHC()
    if mode[0:2] == "VS" and self.processed == False:
      self.processed = True
      self.processVSM()

  def updateParameters(self):
    self.processed = False
    self.readCSV(self.files)

  def processHC(self):
    try:
      total_mass = float(self.totmass_input.text())/1e3
    except:
      total_mass = 1
    try:
      mol_mass = float(self.molmass_input.text())
    except:
      mol_mass = 1e6
    try:
      mass_err = float(self.masserr_input.text())
    except:
      mass_err = 0.02
    try:
      silver_ratio = float(self.silver_input.text())
    except:
      silver_ratio = 0.
    if silver_ratio > 0.99:
      silver_ratio = 0.
    if silver_ratio > 0.:
      sample_mass = total_mass*(1-silver_ratio)
      silver_mass = total_mass*silver_ratio
      for n in range(len(self.dfs)):
        self.dfs[n] = self.dfs[n][pd.notnull(self.dfs[n][hc_cols[9]]) \
                        #& (self.dfs[n][hc_cols[18]] > 80) \
                        & pd.notnull(self.dfs[n][hc_cols[7]])]
        self.dfs[n] = self.dfs[n].sort_values([hc_cols[7]])
        self.dfs[n][hc_cols[9]] = self.dfs[n][hc_cols[9]]/1e6
        silver_cp = cp_ag(self.dfs[n][hc_cols[7]], silver_mass)
        self.dfs[n][hc_cols[9]] = (self.dfs[n][hc_cols[9]]-silver_cp)/sample_mass*mol_mass
        self.dfs[n][hc_cols[10]] = self.dfs[n][hc_cols[10]]/1e6/sample_mass*mol_mass
        self.dfs[n][hc_cols[10]] = self.dfs[n][hc_cols[10]] + self.dfs[n][hc_cols[9]]*mass_err
    else:
      for n in range(len(self.dfs)):
        self.dfs[n] = self.dfs[n][pd.notnull(self.dfs[n][hc_cols[9]]) \
                        #& (self.dfs[n][hc_cols[18]] > 80) \
                        & pd.notnull(self.dfs[n][hc_cols[7]])]
        self.dfs[n] = self.dfs[n].sort_values([hc_cols[7]])
        self.dfs[n][hc_cols[9]] = self.dfs[n][hc_cols[9]]/1e6/total_mass*mol_mass
        self.dfs[n][hc_cols[10]] = self.dfs[n][hc_cols[10]]/1e6/total_mass*mol_mass
        self.dfs[n][hc_cols[10]] = self.dfs[n][hc_cols[10]] + self.dfs[n][hc_cols[9]]*mass_err

  def processVSM(self):
    try:
      sample_mass = float(self.totmass_input.text())/1e3
    except:
      sample_mass = 1
    try:
      mol_mass = float(self.molmass_input.text())
    except:
      mol_mass = 1e6
    for n in range(len(self.dfs)):
      self.dfs[n] = self.dfs[n][pd.notnull(self.dfs[n][vsm_cols[4]])]
      self.dfs[n][vsm_cols[4]] = self.dfs[n][vsm_cols[4]]*mol_mass/sample_mass
      self.dfs[n][vsm_cols[5]] = self.dfs[n][vsm_cols[5]]*mol_mass/sample_mass

  def plot(self):
    mode = str(self.combobox_mode.currentText())
    if mode == "HC Cp vs T":
      self.plotCp()
    if mode == "HC Cp/T vs T":
      self.plotCpoT()
    if mode == "VSM: Chi vs T":
      self.plotChi()
    if mode == "VSM: 1/Chi vs T":
      self.plotInverseChi()
    if mode == "VSM: M vs T":
      self.plotMvT()
    if mode == "VSM: M vs H":
      self.plotMvH()
    self.setTitle()
    self.setPlotOptions()
    self.canvas.draw()
    self.plotted = True

  def plotCp(self):
    colors = cm.rainbow(np.linspace(0, 1, len(self.dfs)*10))
    colors = colors[0:-1:10]
    plt.cla()
    for n, df in enumerate(self.dfs):
      plt.errorbar(df[hc_cols[7]].values, df[hc_cols[9]].values,
          df[hc_cols[10]].values, linestyle = '-',
          marker = 'o', color = colors[n],
          markeredgecolor = colors[n],
          markerfacecolor = colors[n],
          capsize = 0, errorevery = 1,
          label = self.labels[n])
    plt.xlabel("Temperature (K)", fontsize = 14)
    plt.ylabel("$C_p$ (J/K/mol)", fontsize = 14)

  def plotCpoT(self):
    colors = cm.rainbow(np.linspace(0, 1, len(self.dfs)*10))
    colors = colors[0:-1:10]
    plt.cla()
    for n, df in enumerate(self.dfs):
      plt.errorbar(df[hc_cols[7]].values,
          df[hc_cols[9]].values/df[hc_cols[7]].values,
          df[hc_cols[10]].values/df[hc_cols[7]].values, linestyle = '-',
          marker = 'o', color = colors[n],
          markeredgecolor = colors[n],
          markerfacecolor = colors[n],
          capsize = 0, errorevery = 1,
          label = self.labels[n])
    plt.xlabel("Temperature (K)", fontsize = 14)
    plt.ylabel("$C_p/T$ (J/K$^2$/mol)", fontsize = 14)

  def plotChi(self):
    colors = cm.rainbow(np.linspace(0, 1, len(self.dfs)*10))
    colors = colors[0:-1:10]
    plt.cla()
    for n, df in enumerate(self.dfs):
      try:
        field = float(self.field_input.text())
      except:
        field = np.mean(df[vsm_cols[3]])
      df.sort_values([vsm_cols[2]])
      plt.errorbar(df[vsm_cols[2]].values, df[vsm_cols[4]].values/field,
        df[vsm_cols[5]].values/field, linestyle = '-',
        marker = 'o', color = colors[n], markeredgecolor = colors[n],
        markerfacecolor = colors[n],
        capsize = 0, errorevery = 1,
        label = self.labels[n])
    plt.xlabel("Temperature (K)", fontsize = 14)
    plt.ylabel("$\chi$ (emu.mol$_\mathrm{fu}^{-1}$)", fontsize = 14)

  def plotInverseChi(self):
    colors = cm.rainbow(np.linspace(0, 1, len(self.dfs)*10))
    colors = colors[0:-1:10]
    plt.cla()
    for n, df in enumerate(self.dfs):
      try:
        field = float(self.field_input.text())
      except:
        field = np.mean(df[vsm_cols[3]])
      df.sort_values([vsm_cols[2]])
      plt.errorbar(df[vsm_cols[2]].values, field/df[vsm_cols[4]].values,
        field/df[vsm_cols[4]].values/df[vsm_cols[4]].values*df[vsm_cols[5]].values,
        linestyle = '-',
        marker = 'o', color = colors[n], markeredgecolor = colors[n],
        markerfacecolor = colors[n],
        capsize = 0, errorevery = 1,
        label = self.labels[n])
    plt.xlabel("Temperature (K)", fontsize = 14)
    plt.ylabel("$\chi^{-1}$ (mol$_{\mathrm{fu}}$.emu$^{-1}$)", fontsize = 14)

  def plotMvH(self):
    colors = cm.rainbow(np.linspace(0, 1, len(self.dfs)*10))
    colors = colors[0:-1:10]
    plt.cla()
    mu_b = 5584.93
    for n, df in enumerate(self.dfs):
      df.sort_values([vsm_cols[3]])
      plt.errorbar(df[vsm_cols[3]].values/1e4, df[vsm_cols[4]].values/mu_b,
        df[vsm_cols[5]].values/mu_b, linestyle = '-',
        marker = 'o', color = colors[n], markeredgecolor = colors[n],
        markerfacecolor = colors[n],
        capsize = 0, errorevery = 1,
        label = self.labels[n])
    plt.xlabel("Magnetic Field (T)", fontsize = 14)
    plt.ylabel("Moment ($\mu_{\mathrm{B}}$.mol$_{\mathrm{fu}}^{-1}$)", fontsize = 14)

  def plotMvT(self):
    colors = cm.rainbow(np.linspace(0, 1, len(self.dfs)*10))
    colors = colors[0:-1:10]
    plt.cla()
    mu_b = 5584.93
    for n, df in enumerate(self.dfs):
      df.sort_values([vsm_cols[2]])
      plt.errorbar(df[vsm_cols[2]].values, df[vsm_cols[4]].values/mu_b,
        df[vsm_cols[5]].values/mu_b, linestyle = '-',
        marker = 'o', color = colors[n], markeredgecolor = colors[n],
        markerfacecolor = colors[n],
        capsize = 0, errorevery = 1,
        label = self.labels[n])
    plt.xlabel("Temperature (K)", fontsize = 14)
    plt.ylabel("Moment ($\mu_{\mathrm{B}}$.mol$_\mathrm{fu}^{-1}$)", fontsize = 14)

  def setTitle(self):
    compound = convertCompoundString(str(self.compound_input.text()))
    plt.title(compound + " " + str(self.title_input.text()), fontsize = 14)
    ax = self.figure.gca()
    if self.checkbox_legendoff.isChecked():
      return
    handles, labels = ax.get_legend_handles_labels()
    #print(labels)
    handles = [h[0] if isinstance(h, container.ErrorbarContainer) else h for h in handles]
    ax.legend(handles, labels, numpoints = 1, frameon = False,
              loc = 'best', fontsize = 10)

  def exportPlotData(self):
    filepath = QtWidgets.QFileDialog.getExistingDirectory(self, 
                           'Export plotted data at',
                           self.savepath)
    print(filepath)
    self.savepath = filepath
    if self.plotted and (len(self.dfs) > 0):
      mode = str(self.combobox_mode.currentText())
      if mode == "HC Cp vs T" or mode == "HC Cp/T vs T":
        self.exportHC()
      if mode == "VSM: Chi vs T" or mode == "VSM: 1/Chi vs T":
        self.exportChi()
      if mode == "VSM: M vs T":
        self.exportMvT()
      if mode == "VSM: M vs H":
        self.exportMvH()
      else:
        return

  def exportHC(self):
    for n, df in enumerate(self.dfs):
      f = self.savepath + '/' + self.labels[n] + '.dat'
      data = np.array([df[hc_cols[7]].values, df[hc_cols[9]].values, df[hc_cols[10]].values]).T
      np.savetxt(f, data, fmt='%.6e', delimiter = ' ', header = 'Temperature (Kelvin), Specific Heat (J/K/mol), Specific Heat Error)')

  def exportChi(self):
    for n, df in enumerate(self.dfs):
      f = self.savepath + '/' + self.labels[n] + '.dat'
      try:
        field = float(self.field_input.text())
      except:
        field = np.mean(df[vsm_cols[3]])
      data = np.array([df[vsm_cols[2]].values, df[vsm_cols[4]].values/field, df[vsm_cols[5]].values/field]).T
      np.savetxt(f, data, fmt='%.6e', delimiter = ' ', header = 'Temperature (Kelvin), Chi (emu/mol), Chi Error')

  def exportMvT(self):
    for n, df in enumerate(self.dfs):
      mu_b = 5584.93
      f = self.savepath + '/' + self.labels[n] + '.dat'
      data = np.array([df[vsm_cols[2]].values, df[vsm_cols[4]].values/mu_b, df[vsm_cols[5]].values/mu_b]).T
      np.savetxt(f, data, fmt='%.6e', delimiter = ' ', header = 'Temperature (Kelvin), Moment (muB/mol), Moment Error')

  def exportMvH(self):
    for n, df in enumerate(self.dfs):
      mu_b = 5584.93
      f = self.savepath + '/' + self.labels[n] + '.dat'
      data = np.array([df[vsm_cols[3]].values/1e4, df[vsm_cols[4]].values/mu_b, df[vsm_cols[5]].values/mu_b]).T
      np.savetxt(f, data, fmt='%.6e', delimiter = ' ', header = 'Magnetic Field (Tesla), Moment (muB/mol), Moment Error')

  def getFileName(self):
    compound = str(self.compound_input.text())
    if compound == "":
      compound = "Secret"
    return compound + "_PPMS_" + getTimeString()

  def setPlotOptions(self):
    if self.checkbox_grid.isChecked():
      plt.grid(True)
    if self.checkbox_xlog.isChecked():
      plt.xscale('log')
    if self.checkbox_ylog.isChecked():
      plt.yscale('log')

  def clearPlot(self):
    plt.cla()
    self.canvas.draw()
    self.plotted = False

  def clearFiles(self):
    self.processed = False
    self.files = []
    self.dfs = []

  def quit(self):
    sys.exit();

def main():
  app = QtWidgets.QApplication(sys.argv)
  w = PPMSDataViewer()
  app.exec_()

if __name__ == '__main__':
  main()
