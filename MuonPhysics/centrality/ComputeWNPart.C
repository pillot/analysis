/*
 *  ComputeWNPart.C
 *  aliroot_dev
 *
 *  Created by philippe pillot on 11/12/11.
 *  Copyright 2011 SUBATECH. All rights reserved.
 *
 */

void ComputeWNPart(Bool_t weight = kTRUE)
{
  /// - launch the centrality analysis to get the weighted Npart value for each set of parameters
  /// - launch the ComputeErrors macro to get the error and print the results for each centrality bin
  
  gROOT->ProcessLine(Form(".x AnalyzeCentralityB.C(%d,%d,%d,%d,%d)", weight, 64, 4, 662, 546));
  gROOT->ProcessLine(Form(".x AnalyzeCentralityB.C(%d,%d,%d,%d,%d)", weight, 69, 4, 662, 546));
  gROOT->ProcessLine(Form(".x AnalyzeCentralityB.C(%d,%d,%d,%d,%d)", weight, 59, 4, 662, 546));
  gROOT->ProcessLine(Form(".x AnalyzeCentralityB.C(%d,%d,%d,%d,%d)", weight, 64, 0, 662, 546));
  gROOT->ProcessLine(Form(".x AnalyzeCentralityB.C(%d,%d,%d,%d,%d)", weight, 64, 8, 662, 546));
  gROOT->ProcessLine(Form(".x AnalyzeCentralityB.C(%d,%d,%d,%d,%d)", weight, 64, 4, 668, 546));
  gROOT->ProcessLine(Form(".x AnalyzeCentralityB.C(%d,%d,%d,%d,%d)", weight, 64, 4, 656, 546));
  gROOT->ProcessLine(Form(".x AnalyzeCentralityB.C(%d,%d,%d,%d,%d)", weight, 64, 4, 662, 556));
  gROOT->ProcessLine(Form(".x AnalyzeCentralityB.C(%d,%d,%d,%d,%d)", weight, 64, 4, 662, 536));
  gROOT->ProcessLine(Form(".x AnalyzeCentralityB.C(%d,%d,%d,%d,%d)", weight, 64, 4, 668, 556));
  gROOT->ProcessLine(Form(".x AnalyzeCentralityB.C(%d,%d,%d,%d,%d)", weight, 64, 4, 668, 536));
  gROOT->ProcessLine(Form(".x AnalyzeCentralityB.C(%d,%d,%d,%d,%d)", weight, 64, 4, 656, 556));
  gROOT->ProcessLine(Form(".x AnalyzeCentralityB.C(%d,%d,%d,%d,%d)", weight, 64, 4, 656, 536));
  
  gROOT->ProcessLine(Form(".x ComputeErrors.C(%d)",weight));
  
}

