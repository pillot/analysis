#ifndef ALIANALYSISTASKTEST_H
#define ALIANALYSISTASKTEST_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/// \ingroup muondep
/// \class AliAnalysisTaskTest
/// \brief basic task to start
//Author: Philippe Pillot - SUBATECH Nantes

#include "AliAnalysisTaskSE.h"

class TH1F;

class AliAnalysisTaskTest : public AliAnalysisTaskSE {
public:
  
  AliAnalysisTaskTest();
  AliAnalysisTaskTest(const char *name);
  virtual ~AliAnalysisTaskTest();
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *);
  virtual void   Terminate(Option_t *);
  
private:
  
  /// Not implemented
  AliAnalysisTaskTest(const AliAnalysisTaskTest& rhs);
  /// Not implemented
  AliAnalysisTaskTest& operator = (const AliAnalysisTaskTest& rhs);
  
private:
  
  TH1F *fhPt; //!< output histogram
  
  ClassDef(AliAnalysisTaskTest, 1);
};

#endif

