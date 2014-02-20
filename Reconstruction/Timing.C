/*
 *  Timing.C
 *  aliroot_dev
 *
 *  Created by philippe pillot on 21/10/13.
 *  Copyright 2013 Subatech. All rights reserved.
 *
 */


// -------------------------------------------------------------------
void Timing(TString log)
{
  
  AliCounterCollection c;
  c.AddRubric("timer", "Process/RunLocalEventReconstruction/FillTreeR/RunMuonTracking");
  c.AddRubric("time", "Real/CPU/slices");
  c.Init(kTRUE);
  
  TObjArray *timers = c.GetKeyWords("timer").Tokenize(",");
  TIter next(timers);
  TObjString *timer = 0x0;
  TString line;
  while ((timer = static_cast<TObjString*>(next()))) {
    
    gSystem->Exec(Form("grep -i \" %s  R:\" %s > __%s__", timer->GetName(), log.Data(), timer->GetName()));
    
    ifstream inFile(Form("__%s__", timer->GetName()));
    if (inFile.is_open()) while (! inFile.eof() ) {
      
      line.ReadLine(inFile,kTRUE);
      if (line.IsNull()) continue;
      
      Double_t realTime = GetRealTime(line);
      if (realTime < 0) continue;
      
      Double_t CPUTime = GetCPUTime(line);
      if (CPUTime < 0) continue;
      
      Int_t slices = GetSlices(line);
      if (slices < 0) continue;
      
      c.Count(Form("timer:%s/time:Real", timer->GetName()), realTime);
      c.Count(Form("timer:%s/time:CPU", timer->GetName()), CPUTime);
      c.Count(Form("timer:%s/time:slices", timer->GetName()), slices);
      
    }
    inFile.close();
    
    gSystem->Exec(Form("rm -f __%s__", timer->GetName()));
    
  }
  
  c.Print("timer/time");
  
}

// -------------------------------------------------------------------
Double_t GetRealTime(TString line)
{
  line.Remove(0, line.Index("R:")+2);
  line.Remove(line.Index("s"), line.Length());
  if (!line.IsFloat()) return -1.;
  return line.Atof();
}

// -------------------------------------------------------------------
Double_t GetCPUTime(TString line)
{
  line.Remove(0, line.Index("C:")+2);
  line.Remove(line.Index("s"), line.Length());
  if (!line.IsFloat()) return -1.;
  return line.Atof();
}

// -------------------------------------------------------------------
Int_t GetSlices(TString line)
{
  line.Remove(0, line.Index("(")+1);
  line.Remove(line.Index(" s"), line.Length());
  if (!line.IsDec()) return -1;
  return line.Atoi();
}

