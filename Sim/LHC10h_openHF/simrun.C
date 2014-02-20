// #define VERBOSEARGS
// simrun.C
{
// set job and simulation variables as :
// root.exe -b -q simrun.C  --run <x> --chunk <y> --event <n> --process <kPythia6/kPhojet/kPythia6ATLAS_Flat/kPythia6D6T> --field <kNoField/k5kG> --energy <900/2760/10000>

  int nrun = 0;
  int nchunk = 0;
  int nevent = 0;
  int seed = 0;

  char sseed[1024];
  char srun[1024];
  char schunk[1024];
  char sprocess[1024];
  char sfield[1024];
  char senergy[1024];
  char spolar[1024];

  sprintf(srun,"");
  sprintf(schunk,"");
  sprintf(sprocess,"");
  sprintf(sfield,"");
  sprintf(senergy,"");
  sprintf(spolar,"");

  for (int i=0; i< gApplication->Argc();i++){
#ifdef VERBOSEARGS
    printf("Arg  %d:  %s\n",i,gApplication->Argv(i));
#endif
    if (!(strcmp(gApplication->Argv(i),"--run")))
      nrun = atoi(gApplication->Argv(i+1));
    sprintf(srun,"%d",nrun);

    if (!(strcmp(gApplication->Argv(i),"--chunk")))
      nchunk = atoi(gApplication->Argv(i+1));
    sprintf(schunk,"%d",nchunk);

    if (!(strcmp(gApplication->Argv(i),"--event")))
      nevent = atoi(gApplication->Argv(i+1));

    if (!(strcmp(gApplication->Argv(i),"--process")))
      sprintf(sprocess, gApplication->Argv(i+1));

    if (!(strcmp(gApplication->Argv(i),"--field")))
      sprintf(sfield,gApplication->Argv(i+1));

    if (!(strcmp(gApplication->Argv(i),"--energy")))
      sprintf(senergy,gApplication->Argv(i+1));

    if (!(strcmp(gApplication->Argv(i),"--polar")))
      sprintf(spolar,gApplication->Argv(i+1));
  }

  seed = nrun * 10000 + nchunk;
  sprintf(sseed,"%d",seed);

  if (seed==0) {
    fprintf(stderr,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
    fprintf(stderr,"!!!!  WARNING! Seeding variable for MC is 0          !!!!\n");
    fprintf(stderr,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
  } else {
    fprintf(stdout,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
    fprintf(stdout,"!!!  MC Seed is %d \n",seed);
    fprintf(stdout,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
  }
  
// set the seed environment variable
  gSystem->Setenv("CONFIG_SEED",sseed);
  gSystem->Setenv("CONFIG_RUN_TYPE",sprocess); // kPythia6 or kPhojet
  gSystem->Setenv("CONFIG_FIELD",sfield);      // kNoField or k5kG
  gSystem->Setenv("CONFIG_ENERGY",senergy);    // 900 or 10000 (GeV)
  gSystem->Setenv("DC_RUN",srun); // used in AliSimulation.cxx
  gSystem->Setenv("DC_EVENT",schunk); // Not used in Config.C (used for setting seed)
  gSystem->Setenv("CONFIG_POLAR",spolar); // JPsi polarization
  
// Needed to produce simulated RAW data
  gSystem->Setenv("ALIMDC_RAWDB1","./mdc1");
  gSystem->Setenv("ALIMDC_RAWDB2","./mdc2");
  gSystem->Setenv("ALIMDC_TAGDB","./mdc1/tag");
  gSystem->Setenv("ALIMDC_RUNDB","./mdc1/meta");
  cout<< "SIMRUN:: Run " << gSystem->Getenv("DC_RUN")
          << " Chunk " << gSystem->Getenv("DC_EVENT")
	  << " Generator " << gSystem->Getenv("CONFIG_RUN_TYPE")
	  << " Field " << gSystem->Getenv("CONFIG_FIELD")
	  << " Energy " << gSystem->Getenv("CONFIG_ENERGY")
	  << endl;

  cout<<">>>>> SIMULATION <<<<<"<<endl;
  gSystem->Exec(Form("aliroot -b -q sim.C\\(%d\\) > sim.log 2>&1",nevent));
  gSystem->Exec("mv syswatch.log simwatch.log");
  cout<<">>>>> RECONSTRUCTION <<<<<"<<endl;
  gSystem->Exec("aliroot -b -q rec.C > rec.log 2>&1");
  gSystem->Exec("mv syswatch.log recwatch.log");
  cout<<">>>>> TAG <<<<<"<<endl;
  gSystem->Exec("aliroot -b -q tag.C > tag.log 2>&1");
  cout<<">>>>> CHECK ESD <<<<<"<<endl;
  gSystem->Exec("aliroot -b -q CheckESD.C > check.log 2>&1");
  cout<<">>>>> MAKE AOD <<<<<"<<endl;
  gSystem->Exec("aliroot -b -q AODtrain.C > aod.log 2>&1");
  cout<<">>>>> CHECK AOD <<<<<"<<endl;
  gSystem->Exec("aliroot -b -q CheckAOD.C >> check.log 2>&1");
}
