/*
 *  CountCandidates.C
 *  aliroot_dev
 *
 *  Created by philippe pillot on 17/10/13.
 *  Copyright 2013 Subatech. All rights reserved.
 *
 */


// -------------------------------------------------------------------
void CountCandidates(TString log)
{
  
  AliCounterCollection c;
  c.AddRubric("candidates", "all/good");
  c.AddRubric("level", "first/more/st3/st2/st1");
  TString cutKeys = "any";
  for (Int_t j = 1; j <= 10; j++) cutKeys += Form("/>%d", 10000*j);
  c.AddRubric("cut", cutKeys.Data());
  c.Init();
  
  gSystem->Exec(Form("grep \"Number of candidates before cleaning\" %s > __ALL__", log.Data()));
  gSystem->Exec(Form("grep \"Number of good candidates\" %s > __GOOD__", log.Data()));
  
  TString line;
  TString file[2] = {"__ALL__", "__GOOD__"};
  TString candidates[2] = {"candidates:all", "candidates:good"};
  for (Int_t i = 0; i < 2; i++) {
    
    ifstream inFile(file[i].Data());
    if (inFile.is_open()) while (! inFile.eof() ) {
      
      line.ReadLine(inFile,kTRUE);
      if (line.IsNull()) continue;
      
      TString level(GetLevel(line));
      if (level == "level:") continue;
      
      Int_t n = GetNCounts(line);
      if (n <= 0) continue;
      
      c.Count(Form("%s/%s/cut:any", candidates[i].Data(), level.Data()), n);
      for (Int_t j = 1; j <= 10; j++) if (n > 10000*j)
	c.Count(Form("%s/%s/cut:>%d", candidates[i].Data(), level.Data(), 10000*j));
      
    }
    inFile.close();
    
  }
  
  gSystem->Exec("rm -f __ALL__");
  gSystem->Exec("rm -f __GOOD__");
  
  c.Print("level/candidates/cut");
  
}

// -------------------------------------------------------------------
TString GetLevel(TString line)
{
  TString level = "level:";
  if (line.Contains("MakeTrackCandidates")) level += "first";
  else if (line.Contains("MakeMoreTrackCandidates")) level += "more";
  else if (line.Contains("In stations(1..) 3")) level += "st3";
  else if (line.Contains("In stations(1..) 2")) level += "st2";
  else if (line.Contains("In stations(1..) 1")) level += "st1";
  return level;
}

// -------------------------------------------------------------------
Int_t GetNCounts(TString line)
{
  line.Remove(0, line.Index("=")+2);
  line.Remove(line.Index(" "), line.Length());
  if (!line.IsDec()) return -1;
  return line.Atoi();
}

