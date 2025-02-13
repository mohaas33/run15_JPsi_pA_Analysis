//void uDstSkim(Int_t nEvents = 1000, Int_t nFiles = 10, TString InputFileList = "./fileLists/test.list", TString OutputDir = "./outputFiles", TString JobIdName = "mudst_EEMC_matching_test" )
void uDstSkim(Int_t nEvents, Int_t nFiles, TString InputFileList, TString OutputDir, TString JobIdName ) 
{
 
  // Load libraries
  gROOT   -> Macro("loadMuDst.C");
  gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
  loadSharedLibraries();

  gSystem -> Load("uDstSkim.so") ;

  // for BEMC a la $STAR/StRoot/macros/mudst/exampleEmc.C
  // Load St_db_Maker and co
  gSystem->Load("StDbLib.so");
  gSystem->Load("StDbBroker.so");
  gSystem->Load("St_db_Maker");
  gSystem->Load("StMaker");
  // Load Emc libraries
  gSystem->Load("StDaqLib");
  gSystem->Load("StEmcRawMaker");
  gSystem->Load("StEmcADCtoEMaker");
  gSystem->Load("StPreEclMaker");
  gSystem->Load("StEpcMaker");
  // Load TOF stuff (extracted from StPeCMAker doMuEvent_picoFromEmbedding.C)
  gSystem->Load("StBTofUtil");
  gSystem->Load("StBTofCalibMaker");
  // Load TOF stuff (was commented out in Leszek's pp2ppAnalysis.C?)
  gSystem->Load("StBTofHitMaker");
  gSystem->Load("StBTofMatchMaker");
  // gSystem->Load("StBTofUtil");
  gSystem->Load("StBTofGeometry");
  
  //EEMC libraries
  gSystem->Load("StEEmcUtil");
  gSystem->Load("$STAR/StRoot/StEEmcPool");
  gSystem->Load("StEEmcDbMaker");
  gSystem->Load("StEEmcA2EMaker");
  gSystem->Load("StEEmcClusterMaker");
  gSystem->Load("EemcGeomSimple");  
  
  gSystem->Load("StSpinDbMaker");   

  // List of member links in the chain
  StChain*                    chain  =  new StChain ;
  StMuDstMaker*          muDstMaker  =  new StMuDstMaker(0,0,"",InputFileList,"MuDst",nFiles) ;
  // Turn off everything but Primary tracks in order to speed up the analysis and eliminate IO
  //muDstMaker -> SetStatus("*",1) ;                // Turn on all branches
  muDstMaker -> SetStatus("*",0) ;                // Turn off all branches
  //muDstMaker -> SetStatus("MuEventAll",1) ;          // Turn on
  muDstMaker -> SetStatus("MuEvent",1) ;          // Turn on the Event data (esp. Event number)
  muDstMaker -> SetStatus("PrimaryTracks",1) ;    // Turn on the primary track data
  muDstMaker -> SetStatus("PrimaryVertices",1) ;    // Turn on the primary vertex data
  //muDstMaker -> SetStatus("GlobalTracks",1) ;    // Turn on the global track data
  muDstMaker -> SetStatus("EmcAll",1) ;    // Turn on all standard Emc branches
  muDstMaker -> SetStatus("TofAll" ,1) ; //  Turn on all standard Tof branches
  muDstMaker -> SetStatus("BTofAll" ,1) ; //  Turn on all standard BTof branches
  muDstMaker -> SetStatus("TofHit",1) ;    // Turn on some TOF stuff
  muDstMaker -> SetStatus("TofData",1) ;    // Turn on some TOF stuff
  muDstMaker -> SetStatus("TofRawData",1) ;    // Turn on some TOF stuff
  muDstMaker -> SetStatus("BTofHit",1) ;    // Turn on some TOF stuff
  muDstMaker -> SetStatus("BTofRawHit",1) ;    // Turn on some TOF stuff
  muDstMaker -> SetStatus("BTofHeader",1) ;    // Turn on some TOF stuff
  muDstMaker -> SetStatus("pp2pp",1); // RP stuff
  //muDstMaker -> SetStatus("McEvent",1); // MC ?
  muDstMaker -> SetStatus("MCAll",1); // This gets MC vertices & Tracks

  // for BEMC a la $STAR/StRoot/macros/mudst/exampleEmc.C
  // Need St_db_Maker for Emc calibration
  St_db_Maker    *db1  = new St_db_Maker("db","$HOME/StarDb","MySQL:StarDb","$STAR/StarDb"); // from exampleEmc.C
  StEEmcDbMaker  *myMk = new StEEmcDbMaker("eemcDb"); //for the EEMC pedestals
  

  // Database for Spin data
  StSpinDbMaker *spDb=new StSpinDbMaker("spinDb");

  // Maker to apply calibration -- for the BEMC
  StEmcADCtoEMaker *adc_to_e=new StEmcADCtoEMaker(); adc_to_e->setPrint(kFALSE);
  
  // Makers for clusterfinding
  StPreEclMaker *pre_ecl=new StPreEclMaker(); pre_ecl->setPrint(kFALSE);
  StEpcMaker *epc=new StEpcMaker(); epc->setPrint(kFALSE);


  //-----> Code
  uDstSkimMaker* AnalysisCode  =  new uDstSkimMaker(muDstMaker) ;

 
  // Miscellaneous things we need before starting the chain
  TString Name = JobIdName ;
  Name.Append(".histograms.root") ;
  AnalysisCode -> SetOutputFileName(Name) ;       // Name the output file for histograms
  if ( nEvents == 0 )  nEvents = muDstMaker->chain()->GetEntries();      // Take all events in nFiles if nEvents = 0


  
  // Loop over the links in the chain
  chain -> Init() ;


  // Reset clusterfinding params
  cout << "RESETTING CLUSTERING PARAMS" << endl;
  Int_t sizeMax = 4; Float_t energySeed = 0.7; Float_t energyAdd  = 0.001; Float_t minTotE  = 0.1; // contructor defaults
  energySeed = 0.3; // lower Seed energy
  pre_ecl->SetClusterConditions("bemc", sizeMax, energySeed, energyAdd, minTotE, kFALSE);

  //chain -> EventLoop(1,nEvents) ;
  for (Int_t ievt = 1; ievt <= nEvents; ievt++) {
    chain->Clear(); // DON'T FORGET THIS!!!
    Int_t stat = chain->Make();
    if (ievt%10000==0) {cout << "ievt: " << ievt << " stat: " << stat << endl;}
  }
  chain->Finish(); 

  // Cleanup
  delete chain ;

}

