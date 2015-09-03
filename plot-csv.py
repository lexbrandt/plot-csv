#! /usr/bin/env python
#! -*- coding:utf-8 -*-

import sys
import csv
import re
import os
import math
import random
import argparse
import re


import prettyplotlib as ppl
import numpy as np

# This is "import matplotlib.pyplot as plt" from the prettyplotlib library
import matplotlib.pyplot as plt

# This is "import matplotlib as mpl" from the prettyplotlib library
import matplotlib as mpl

import brewer2mpl


# import pylab as P


######################################################################
## input parser
######################################################################

def setParser():
    parser = argparse.ArgumentParser()
    parser.add_argument("inputFiles", nargs='+',
                        type=argparse.FileType('rb'),
                        help="Input files")
    parser.add_argument("--solsPath",
                        help="Path for instances' solutions",
                        type=str, default='')
    parser.add_argument("--solsExt",
                        help="Extension for instances' solutions",
                        type=str, default='.grasp')
    parser.add_argument("--add_time", nargs='+', type=argparse.FileType('rb'),
                        help="File(s) with additional time to consider")
    parser.add_argument("-b", "--block",
                        help="Block for interaction",
                        action="store_true")
    parser.add_argument("-v", "--verbose",
                        help="Increase output verbosity",
                        action='count', default=0)

    return parser

######################################################################
## process
######################################################################

def read_table(iFile):
    reader = csv.reader(iFile, delimiter=";", skipinitialspace=True)
    header = reader.next()

    table = []
    for i, row in enumerate(reader):
        table.append(row)

    return header, table
    

def per_type(table, header):
    table.sort()

    type_table = {}
    for row in table:
        # instance identification
        tmp = re.search("(^.*/)([a-z0-9_]+)_([0-9]+)_([0-9]+).in", row[0])
        instancePath, instanceType, instanceSize, instanceNumber = tmp.groups()
        instanceSize = int(instanceSize)
        instanceNumber = int(instanceNumber)

        config_index = header.index("comment")
        config = row[config_index]
        if config not in type_table:
            type_table[config] = {}
        if instanceType not in type_table[config]:
            type_table[config][instanceType] = {}
        if instanceSize not in type_table[config][instanceType]:
            type_table[config][instanceType][instanceSize] = {}

        type_header = ["instanceNumber"] + header[1:]
        for i,col in enumerate(row):
            if type_header[i] not in type_table[config][instanceType][instanceSize].keys():
                type_table[config][instanceType][instanceSize][type_header[i]] = []

            if type_header[i] == "instanceNumber":
              type_table[config][instanceType][instanceSize][type_header[i]].append(instanceNumber)
            else:
              type_table[config][instanceType][instanceSize][type_header[i]].append(col)

    return type_header, type_table


def process(iFile):
    header, table = read_table(iFile)

    type_header, type_table = per_type(table, header)

    return header, table, type_header, type_table

######################################################################
## plots
######################################################################

meancolor = brewer2mpl.get_map('Set2', 'qualitative', 3).mpl_colors[0]
meanprops = {'marker':'o', 'markerfacecolor':meancolor, 'markeredgecolor':'w'}

def createTableStatsReduction(type_header, type_table, redConfig, comConfig, instanceType,
                              fileName = None):

  if redConfig not in type_table:
    raise Exception("Config \""+redConfig+"\" not found!")
  if comConfig not in type_table:
    raise Exception("Config \""+comConfig+"\" not found!")
  if instanceType not in type_table[redConfig]:
    raise Exception("Instance type \""+instanceType+"\" not found!")

  with open(fileName, "w") as oFile:

    sizes = sorted(type_table[comConfig][instanceType].iterkeys())
    for size in sizes:
      if size not in type_table[redConfig][instanceType]:
        continue
      
      nConstrs = []
      pConstrs = []

      comNConstrs = map(int, type_table[comConfig][instanceType][size]['nConstr'])
      comInstNumbers = type_table[comConfig][instanceType][size]['instanceNumber']

      for comNConstrs_inst, comInstNumber in zip(comNConstrs, comInstNumbers):

        redInstIndex = type_table[redConfig][instanceType][size]['instanceNumber'].index(comInstNumber)
        redNConstrs_inst = int(type_table[redConfig][instanceType][size]['nConstr'][redInstIndex])

        nConstrs.append(redNConstrs_inst)
        pConstrs.append(float(redNConstrs_inst)/float(comNConstrs_inst) * 100.0)

      nConstrs = np.array(nConstrs)
      pConstrs = np.array(pConstrs)


      totalVars = np.array( map(float, type_table[redConfig][instanceType][size]['nVarsTotal']) )
      totalArcs = totalVars -1
  
      # fixedArcs   = np.array( map(int, type_table[redConfig][instanceType][size]['preproc.fixedArcs']) )
      # blockedArcs = np.array( map(int, type_table[redConfig][instanceType][size]['preproc.blockedArcs']) )
      # freeArcs = totalArcs - fixedArcs - blockedArcs

      freeVars = np.array( map(int, type_table[redConfig][instanceType][size]['nVarsFree']) )
      freeArcs = freeVars -1

      pFreeArcs = (freeArcs/totalArcs) * 100.0
  
      row = "%d" % (size)
      places = 2
      zerof = "%.*f" % (places, 0)
      zeroi = "%*d"  % (places, 0)

      row += (" & %.*f" % (places, nConstrs.mean())).replace(zerof, zeroi)
      row += (" & %.*f" % (places, pConstrs.mean())).replace(zerof, zeroi)

      # row +=  " & %d"   % (totalVars.mean())
      # row += (" & %.*f" % (places, fixedArcs.mean())).replace(zerof, zeroi)
      row += (" & %.*f" % (places, freeArcs.mean())).replace(zerof, zeroi)
      row += (" & %.*f" % (places, pFreeArcs.mean())).replace(zerof, zeroi)
      # row += (" & %.*f" % (places, blockedArcs.mean())).replace(zerof, zeroi)

      row +=  " \\\\\n"
      
      oFile.write(row.replace('_', '\_').encode("utf8"))
      oFile.write("\\hline\n".encode("utf8"))


def getNEdges(solPath):
  with open(solPath, 'r') as s:
    n = int(s.readline())
    return n


def boxplotFixedEdgesSols(type_header, type_table, config, instanceType,
                          solsPath, solsExt, figName = None):
  
  if config not in type_table:
    raise Exception("Config \""+config+"\" not found!")
  if instanceType not in type_table[config]:
    raise Exception("Instance type \""+instanceType+"\" not found!")
  
  fig, ax = plt.subplots(1)

  data = []
  for size in sorted(type_table[config][instanceType].iterkeys()):
    fixedEdges = np.array(map(int, type_table[config][instanceType][size]["preproc.fixedEdges"]))
    instanceNumbers = np.array(type_table[config][instanceType][size]["instanceNumber"])

    fixedEdgesProp = []

    for fEdges, instanceNumber in zip(fixedEdges, instanceNumbers):
      instanceName = "%s_%03d_%02d%s" % (instanceType, size, instanceNumber, solsExt)
      solEdges = getNEdges(os.path.join(solsPath, instanceName))

      fEdgesProp = float(fEdges)/float(solEdges) * 100.0
      fixedEdgesProp.append(fEdgesProp)

    data.append(fixedEdgesProp)

  ret = ppl.boxplot(ax, data, xticklabels=map(str, sorted(type_table[config][instanceType].iterkeys())),
                    widths=0.3, showmeans=True, meanprops=meanprops)

  ax.set_xlabel(u'# Pontos')
  ax.set_ylabel(u'% Arestas Fixadas')

  if figName != None:
    fig.savefig(figName, bbox_inches='tight')


def scatterFreeByMissingEdges(type_header, type_table, config, instanceType,
                              solsPath, solsExt, figName = None):
  if config not in type_table:
    raise Exception("Config \""+config+"\" not found!")
  if instanceType not in type_table[config]:
    raise Exception("Instance type \""+instanceType+"\" not found!")
  
  fig, ax = plt.subplots(1)

  for size in sorted(type_table[config][instanceType].iterkeys()): # instanceType
    fixedEdges = np.array(map(int, type_table[config][instanceType][size]["preproc.fixedEdges"]))
    blockedEdges = np.array(map(int, type_table[config][instanceType][size]["preproc.blockedEdges"]))
    instanceNumbers = np.array(type_table[config][instanceType][size]["instanceNumber"])

    if size == 80:
      solsPath += '_80-90'
    missingEdges = []
    freeEdges = []
    for fEdges, bEdges, instanceNumber in zip(fixedEdges, blockedEdges, instanceNumbers):
      instanceName = "%s_%03d_%02d%s" % (instanceType, size, instanceNumber, solsExt)
      solEdges = getNEdges(os.path.join(solsPath, instanceName))

      mEdges = solEdges - fEdges
      missingEdges.append(mEdges)

      frEdges = ((size-1)*size)/2 - (fEdges + bEdges)
      freeEdges.append(frEdges)

    ppl.scatter(ax, missingEdges, freeEdges, label=str(size))
  
  ppl.legend(ax, loc="lower right")
  ax.set_xlabel(u'Arestas faltantes')
  ax.set_ylabel(u'Arestas livres')
  # ax.set_aspect('equal')

  ax.set_xlim((0, ax.get_xlim()[1]))
  ax.set_ylim((0, ax.get_ylim()[1]))

  # ax.set_title('prettyplotlib `scatter` example\nshowing default color cycle and scatter params')
  if figName != None:
    fig.savefig(figName, bbox_inches='tight')
  

def boxplotFixedPaths(type_header, type_table, config, instanceType, figName = None):
  
  if config not in type_table:
    raise Exception("Config \""+config+"\" not found!")
  if instanceType not in type_table[config]:
    raise Exception("Instance type \""+instanceType+"\" not found!")
  
  fig, ax = plt.subplots(1)

  data = []
  for size in sorted(type_table[config][instanceType].iterkeys()):
    fixedpaths = np.array(map(int, type_table[config][instanceType][size]["preproc.fixedPaths"]))
    nPaths = (size*(size-1))/2.0
    fixedPathsProp = fixedpaths/nPaths * 100.0

    data.append(fixedPathsProp)

  ret = ppl.boxplot(ax, data, xticklabels=map(str, sorted(type_table[config][instanceType].iterkeys())),
                    widths=0.3, showmeans=True, meanprops=meanprops)

  ax.set_xlabel(u'# Pontos')
  ax.set_ylabel(u'% Caminhos Fixados')

  if figName != None:
    fig.savefig(figName, bbox_inches='tight')


def stackedBarsBlockedArcsByRoutine(type_header, type_table, config, instanceType,
                                    figName = None, routine_map = None):

  if config not in type_table:
    raise Exception("Config \""+config+"\" not found!")
  if instanceType not in type_table[config]:
    raise Exception("Instance type \""+instanceType+"\" not found!")

  fig, ax = plt.subplots(1)

  routines = []
  rgexs = r'[^.]+\.([^.]+).blockedArcs'
  routineR = re.compile(rgexs)
  for h in type_header:
    routine = routineR.search(h)
    if routine != None:
      routines.append(routine.group(1))

  data = []
  sizes = sorted(type_table[config][instanceType].iterkeys())
  for size in sizes:
  
    totalVars = np.array( map(int, type_table[config][instanceType][size]['nVarsTotal']) )
    totalArcs = totalVars -1
    for index, routine in enumerate(routines):
      
      if index >= len(data):
        data.append([])

      columnName = 'preproc.'+routine+'.blockedArcs'
      blockedArcs = np.array( map(int, type_table[config][instanceType][size][columnName]) )

      propBlockedArcs = np.array( map(lambda r, t: float(r)/float(t), blockedArcs, totalArcs) )
      avgPropBlockedArcs = propBlockedArcs.mean() * 100.0

      data[index].append(avgPropBlockedArcs)

  # Get "Set2" colors from ColorBrewer (all colorbrewer scales: http://bl.ocks.org/mbostock/5577023)
  colors = brewer2mpl.get_map('Set2', 'qualitative', len(routines)).mpl_colors

  clean_data     = []
  clean_routines = []
  clean_colors   = []
  for index, routine in enumerate(routines):
    if data[index].count(0) != len(data[index]):
      clean_data.append(data[index])
      clean_routines.append(routine)
      clean_colors.append(colors[index])

  data = np.array(clean_data)
  routines = clean_routines
  colors = clean_colors

  if routine_map != None:
    # labels = np.array( [routine_map[r] for r in routines] )
    labels = np.array( map(routine_map.get, routines) )
  else:
    labels = np.array(routines)

  bottom = np.vstack((np.zeros((data.shape[1],), dtype=data.dtype),
                      np.cumsum(data, axis=0)[:-1]))

  width = 0.8
  ind = [x-width/2 for x in range(1, len(data[0])+1)]

  for dat, lab, bot, col in reversed(zip(data, labels, bottom, colors)):
    ppl.bar(ax, ind, dat, width, grid='y',
             bottom=bot, label=lab, color=col)

  ax.set_xlabel(u'# Pontos')
  ax.set_xticks(range(1, len(sizes)+1))
  ax.set_xticklabels(map(str, sizes))
  ax.set_xlim(0.5, len(sizes)+0.5)
  ax.set_ylabel(u'% Arcos bloqueados')

  ppl.legend(ax, loc="lower right")

  if figName != None and type(figName) == str:
    fig.savefig(figName, bbox_inches='tight')


def getNEdgesArray(type_table, type_header, config, instanceType, size, solsPath, solsExt):
  instanceNumbers = np.array(type_table[config][instanceType][size]["instanceNumber"])
  if size == 80:
    solsPath += '_80-90'
  solsNEdges = []
  for instanceNumber in instanceNumbers:
    instanceName = "%s_%03d_%02d%s" % (instanceType, size, instanceNumber, solsExt)
    solEdges = getNEdges(os.path.join(solsPath, instanceName))
    solsNEdges.append(solEdges)

  return np.array(solsNEdges)

def stackedBarsFixedEdgesByRoutine(type_header, type_table, config, instanceType,
                                   solsPath = None, solsExt = None, useTotalEdges = False,
                                   figName = None, routine_map = None):

  if config not in type_table:
    raise Exception("Config \""+config+"\" not found!")
  if instanceType not in type_table[config]:
    raise Exception("Instance type \""+instanceType+"\" not found!")

  fig, ax = plt.subplots(1)

  routines = []
  rgexs = r'[^.]+\.([^.]+).fixedEdges'
  routineR = re.compile(rgexs)
  for h in type_header:
    routine = routineR.search(h)
    if routine != None:
      routines.append(routine.group(1))

  data = []
  sizes = sorted(type_table[config][instanceType].iterkeys())
  for size in sizes:
  
    if useTotalEdges:
      totalVars = np.array( map(int, type_table[config][instanceType][size]['nVarsTotal']) )
      totalEdges = [size*(size-1)/2] * len(totalVars)
      refNEdges = totalEdges
    else:
      solsNEdges = getNEdgesArray(type_table, type_header, config, instanceType, size, solsPath, solsExt)
      refNEdges = solsNEdges

    for index, routine in enumerate(routines):
      
      if index >= len(data):
        data.append([])

      columnName = 'preproc.'+routine+'.fixedEdges'
      fixedEdges = np.array( map(int, type_table[config][instanceType][size][columnName]) )

      propFixedEdges = np.array( map(lambda r, t: float(r)/float(t), fixedEdges, refNEdges) )
      avgPropFixedEdges = propFixedEdges.mean() * 100.0

      data[index].append(avgPropFixedEdges)

  # Get "Set2" colors from ColorBrewer (all colorbrewer scales: http://bl.ocks.org/mbostock/5577023)
  colors = brewer2mpl.get_map('Set2', 'qualitative', max(3, len(routines))).mpl_colors

  clean_data     = []
  clean_routines = []
  clean_colors   = []
  for index, routine in enumerate(routines):
    if data[index].count(0) != len(data[index]):
      clean_data.append(data[index])
      clean_routines.append(routine)
      clean_colors.append(colors[index])

  data = np.array(clean_data)
  routines = clean_routines
  colors = clean_colors

  if routine_map != None:
    # labels = np.array( [routine_map[r] for r in routines] )
    labels = np.array( map(routine_map.get, routines) )
  else:
    labels = np.array(routines)

  bottom = np.vstack((np.zeros((data.shape[1],), dtype=data.dtype),
                      np.cumsum(data, axis=0)[:-1]))

  width = 0.8
  ind = [x-width/2 for x in range(1, len(data[0])+1)]

  for dat, lab, bot, col in reversed(zip(data, labels, bottom, colors)):
    ppl.bar(ax, ind, dat, width, grid='y',
             bottom=bot, label=lab, color=col)

  ax.set_xlabel(u'# Pontos')
  ax.set_xticks(range(1, len(sizes)+1))
  ax.set_xticklabels(map(str, sizes))
  ax.set_xlim(0.5, len(sizes)+0.5)
  ax.set_ylabel(u'% Arestas fixadas')

  ppl.legend(ax, loc="upper right")

  if figName != None and type(figName) == str:
    fig.savefig(figName, bbox_inches='tight')





def stackedBarsBlockedEdgesByRoutine(type_header, type_table, config, instanceType,
                                     solsPath = None, solsExt = None, useTotalEdges = False,
                                     figName = None, routine_map = None):

  if config not in type_table:
    raise Exception("Config \""+config+"\" not found!")
  if instanceType not in type_table[config]:
    raise Exception("Instance type \""+instanceType+"\" not found!")

  fig, ax = plt.subplots(1)

  routines = []
  rgexs = r'[^.]+\.([^.]+).blockedEdges'
  routineR = re.compile(rgexs)
  for h in type_header:
    routine = routineR.search(h)
    if routine != None:
      routines.append(routine.group(1))

  data = []
  sizes = sorted(type_table[config][instanceType].iterkeys())
  for size in sizes:

    if useTotalEdges:
      totalVars = np.array( map(int, type_table[config][instanceType][size]['nVarsTotal']) )
      totalEdges = [size*(size-1)/2] * len(totalVars)
      refNEdges = totalEdges
    else:
      solsNEdges = getNEdgesArray(type_table, type_header, config, instanceType, size,
                                  solsPath, solsExt)
      totalEdges = [size*(size-1)/2] * len(solsNEdges)
      refNEdges = totalEdges - solsNEdges

    for index, routine in enumerate(routines):
      
      if index >= len(data):
        data.append([])

      columnName = 'preproc.'+routine+'.blockedEdges'
      blockedEdges = np.array( map(int, type_table[config][instanceType][size][columnName]) )

      propBlockedEdges = np.array( map(lambda r, t: float(r)/float(t), blockedEdges, refNEdges) )
      avgPropBlockedEdges = propBlockedEdges.mean() * 100.0

      if figName == 'bEdges-exclReg-routines-mdtp.pdf':
        print size, ": bEdges:", routine, ":", avgPropBlockedEdges

      data[index].append(avgPropBlockedEdges)


  # Get "Set2" colors from ColorBrewer (all colorbrewer scales: http://bl.ocks.org/mbostock/5577023)
  colors = brewer2mpl.get_map('Set2', 'qualitative', max(len(routines), 3)).mpl_colors

  clean_data     = []
  clean_routines = []
  clean_colors   = []
  for index, routine in enumerate(routines):
    if data[index].count(0) != len(data[index]):
      clean_data.append(data[index])
      clean_routines.append(routine)
      clean_colors.append(colors[index])

  data = np.array(clean_data)
  routines = clean_routines
  colors = clean_colors

  if routine_map != None:
    # labels = np.array( [routine_map[r] for r in routines] )
    labels = np.array( map(routine_map.get, routines) )
  else:
    labels = np.array(routines)

  bottom = np.vstack((np.zeros((data.shape[1],), dtype=data.dtype),
                      np.cumsum(data, axis=0)[:-1]))

  width = 0.8
  ind = [x-width/2 for x in range(1, len(data[0])+1)]

  for dat, lab, bot, col in reversed(zip(data, labels, bottom, colors)):
    ppl.bar(ax, ind, dat, width, grid='y',
             bottom=bot, label=lab, color=col)

  ax.set_xlabel(u'# Pontos')
  ax.set_xticks(range(1, len(sizes)+1))
  ax.set_xticklabels(map(str, sizes))
  ax.set_xlim(0.5, len(sizes)+0.5)
  ax.set_ylabel(u'% Arestas bloqueadas')

  ppl.legend(ax, loc="lower right")

  if figName != None and type(figName) == str:
    fig.savefig(figName, bbox_inches='tight')


def createResultsTableConfig(type_header, type_table, fileName, config, instanceType,
                             configs_map = None, addt_table = None, addt_configs = [],
                             addt_config = None):
  data = {}

  if config not in type_table:
    raise Exception("Config \""+config+"\" not found!")
  if instanceType not in type_table[config]:
    print type_table[config].keys()
    raise Exception("Instance type \""+instanceType+"\" not found!")

  for size in type_table[config][instanceType]:

    if size not in data.iterkeys():
      data[size] = []

    # Optimal solutions
    nOptimalSols = 0
    timeOptimalSolsSum = 0.0

    # Feasible solutions
    nFeasibleSols = 0
    gapSum = 0.0

    # General
    nInstances   = 0.0
    nNodesSum = 0.0

    for i, index in enumerate(type_table[config][instanceType][size]["instanceNumber"]):
      
      status = type_table[config][instanceType][size]["status"][i]

      if status == "Optimal":
        nOptimalSols += 1
        time = float(type_table[config][instanceType][size]["totalTime"][i])
        if config in addt_configs:
          assert addt_table != None, "Missing necessary argument: addt_table"
          assert addt_config != None, "Missing necessary argument: addt_config"
          addt_i = addt_table[addt_config][instanceType][size]["instanceNumber"].index(index)
          addt = float(addt_table[addt_config][instanceType][size]["time"][addt_i])
          time += addt/1000.0 # ms -> s
        timeOptimalSolsSum += time

        # nNodes = int(type_table[config][instanceType][size]["nNodes"][i])
        # nNodesSum += nNodes

      elif status == "Feasible":
        nFeasibleSols += 1
        gap = float(type_table[config][instanceType][size]["Gap"][i])
        gapSum += gap

      nInstances += 1

    avgTimeOptimalSols = 0.0
    if nOptimalSols != 0:
      avgTimeOptimalSols = timeOptimalSolsSum / nOptimalSols

    avgGap = 0.0
    if nFeasibleSols != 0:
      avgGap = gapSum / nFeasibleSols

    avgNNodes = 0.0
    if nOptimalSols != 0:
      avgNNodes = nNodesSum / nOptimalSols

    data[size].append(nOptimalSols)
    data[size].append(avgTimeOptimalSols)
    data[size].append(avgNNodes)

    data[size].append(nFeasibleSols)
    data[size].append(avgGap)


  with open(fileName, "w") as oFile:
  
    for size in sorted(data.iterkeys()):
      if configs_map != None:
        config_label = configs_map[config]
      else:
        config_label = config

      row = "%d" % (size)
      row += " & %d" % (data[size][0])

      if data[size][0] != 0:
        row += " & %.1f" % (data[size][1])
        # row += " & %.1f" % (data[size][2])
      else:
        row += " & --"
        # row += " & --"
        assert data[size][1] == 0, "Non zero run time " + str(data[size][1])
        assert data[size][2] == 0, "Non zero nNodes " + str(data[size][2])

      if data[size][0] != nInstances:
        row += " & %d" % (data[size][3])
        row += " & %.1f" % (data[size][4]*100.0)
      else:
        row += ' & --'
        row += ' & --'
        assert data[size][3] == 0, "Non feasible sols"
        assert data[size][4] == 0.0, "Non zero gap"

      row +=  " \\\\\n"
      
      oFile.write(row.replace('_', '\_').encode("utf8"))
      oFile.write("\\hline\n".encode("utf8"))

def boxplotEllipseVsExclReg(type_header, type_table, config, instanceType, figName = None):
  
  if config not in type_table:
    raise Exception("Config \""+config+"\" not found!")
  if instanceType not in type_table[config]:
    raise Exception("Instance type \""+instanceType+"\" not found!")
  
  fig, ax = plt.subplots(1)

  data = []
  for size in sorted(type_table[config][instanceType].iterkeys()):
    byEllipse = np.array(map(float, type_table[config][instanceType][size]["excludedByEllipse"]))
    exclugReg = np.array(map(float, type_table[config][instanceType][size]["exclusionReg"]))

    blockedEdgesPropInc = ((byEllipse - exclugReg) / exclugReg) * 100.0

    if size == 70:
      print "box:", byEllipse.mean(), exclugReg.mean()

    data.append(blockedEdgesPropInc)

  ret = ppl.boxplot(ax, data, xticklabels=map(str, sorted(type_table[config][instanceType].iterkeys())),
                    widths=0.3, showmeans=True, meanprops=meanprops)

  ax.set_xlabel(u'# Pontos')
  ax.set_ylabel(u'Incremento % Arestas Bloqueadas')

  if figName != None:
    fig.savefig(figName, bbox_inches='tight')


def main():
  parser = setParser()
  args = parser.parse_args()
  
  block = args.block
  
  type_header = None
  type_table = {}

  for f in args.inputFiles:
    f_header, f_table, f_type_header, f_type_table = process(f)
    if type_header == None:
      type_header = f_type_header
    elif type_header != f_type_header:
      raise Exception("File {} has different header from the previous one(s)"
                      .format(f.name))
    type_table.update(f_type_table)
    f.close()
    
    addt_header = None
    addt_table = {}

  if args.add_time != None:
    for f in args.add_time:
      f_header, f_table, f_addt_header, f_addt_table = process(f)
      if addt_header == None:
        addt_header = f_addt_header
      elif addt_header != f_addt_header:
        raise Exception("File {} has different header from the previous one(s)"
                        .format(f.name))
      addt_table.update(f_addt_table)
      f.close()


  configs_map_alter = {
    'gps' : r'\eqref{eq:strongCouplingXX}',
    'gps_pathCoupl' : r'\eqref{eq:pathCoupl}',
    'gps_pathLB' : r'\eqref{eq:pathLB}',
    }
  
  routine_map = {
    'cleanExtremes' : u'Limpa extremos',
    'outOfEllipse'  : u'Eliminação por Elipse',
    'bigEdges'      : u'Eliminação de Arcos Longos',
    }

  createTableStatsReduction(type_header, type_table, 'gps_dryRun', 'g_s_dryRun', 'uniform',
                            fileName = 'stats-preproc-lc-triang-data.tex')


  boxplotFixedEdgesSols(type_header, type_table, 'gps_dryRun', 'uniform',
                        args.solsPath, args.solsExt, figName = 'preproc-stats-fe-triang.pdf')

  scatterFreeByMissingEdges(type_header, type_table, 'gps_dryRun', 'uniform',
                            args.solsPath, args.solsExt, figName = 'free-edges-mdtp-stats.pdf')

  boxplotFixedPaths(type_header, type_table, 'gps_dryRun', 'uniform',
                    figName = 'preproc-stats-fp-triang.pdf')

  routine_map = {
    # 'cleanExtremes'   : u'Limpa extremos',
    'outOfEllipse'      : u'Eliminação baseada em Elipse',
    # 'bigEdges'        : u'Eliminação de Arcos Longos',
    'excludedByEllipse' : u'Eliminação por Elipse',
    'fixUniquePath'     : u'Único caminho possível',
    'fixNotCrossed'     : u'Aresta livre de Cruzamentos'
    }


  stackedBarsBlockedArcsByRoutine(type_header, type_table,
                                  'gps_dryRun', 'uniform',
                                  figName = 'bArcs-routines-mdtp.pdf',
                                  routine_map = routine_map
                                  )

  stackedBarsFixedEdgesByRoutine(type_header, type_table,
                                 'gps_dryRun', 'uniform',
                                 useTotalEdges = True,
                                 figName = 'fEdges-routines-mdtp.pdf',
                                 routine_map = routine_map
                                 )

  stackedBarsFixedEdgesByRoutine(type_header, type_table,
                                 'gps_dryRun', 'uniform',
                                 args.solsPath, args.solsExt,
                                 figName = 'fEdges-sol-routines-mdtp.pdf',
                                 routine_map = routine_map
                                 )

  stackedBarsBlockedEdgesByRoutine(type_header, type_table,
                                   'gps_dryRun', 'uniform',
                                   useTotalEdges = True,
                                   figName = 'bEdges-routines-mdtp.pdf',
                                   routine_map = routine_map
                                   )

  stackedBarsBlockedEdgesByRoutine(type_header, type_table,
                                   'gps_dryRun', 'uniform',
                                   solsPath = args.solsPath, solsExt = args.solsExt,
                                   useTotalEdges = False,
                                   figName = 'bEdges-sol-routines-mdtp.pdf',
                                   routine_map = routine_map
                                   )

  ## exclReg

  scatterFreeByMissingEdges(type_header, type_table, 'gps_exclReg', 'uniform',
                            args.solsPath, args.solsExt, figName = 'free-edges-exclReg-mdtp-stats.pdf')



  stackedBarsBlockedArcsByRoutine(type_header, type_table,
                                  'gps_exclReg', 'uniform',
                                  figName = 'bArcs-exclReg-routines-mdtp.pdf',
                                  routine_map = routine_map
                                  )

  stackedBarsBlockedEdgesByRoutine(type_header, type_table,
                                  'gps_exclReg', 'uniform',
                                   solsPath = args.solsPath, solsExt = args.solsExt,
                                   figName = 'bEdges-exclReg-routines-mdtp.pdf',
                                   routine_map = routine_map
                                   )

  stackedBarsFixedEdgesByRoutine(type_header, type_table,
                                 'gps_exclReg', 'uniform',
                                 args.solsPath, args.solsExt,
                                 figName = 'fEdges-exclReg-routines-mdtp.pdf',
                                 routine_map = routine_map
                                 )

  createResultsTableConfig(type_header, type_table, 'stats-GPS-exclReg-triang-data.tex',
                           'gps_exclReg', 'uniform',
                           configs_map = None,
                           addt_table=addt_table,
                           addt_configs=['gps_exclReg'],
                           addt_config = 'grasp')

  boxplotEllipseVsExclReg(type_header, type_table, 'gps_exclReg', 'uniform',
                          figName = 'ellipse-vs-exclReg-be-triang.pdf')
  
if __name__ == "__main__":
  main()
