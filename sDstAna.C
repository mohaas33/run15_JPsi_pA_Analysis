
void sDstAna( TString InputFileList = "./input_file_list_local.list", TString OutputDir = "./", TString JobIdName = "./interactive_output")//Int_t nEvents, Int_t nFiles, TString InputFileList, TString OutputDir, TString JobIdName ) 
{
 
  // Load libraries
  gROOT   -> Macro("loadMuDst.C");
  gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
  loadSharedLibraries();

  gSystem -> Load("sDstAna.so") ;

	std::cout << "Input FileList: " << InputFileList << std::endl;
	std::cout << "Output dir: " << OutputDir << std::endl;
	std::cout << "JobIdName: " << JobIdName << std::endl;

	TString fileType_ROOT = ".root";
  TString OutputFileName = JobIdName + fileType_ROOT;	
  cout << "Output file: " << OutputFileName << endl;

	std::cout << "OutputFileName: " << OutputFileName << std::endl;


  //-----> Code
  sDstAnaMaker* AnalysisCode  =  new sDstAnaMaker() ;
  //TString Name = "/gpfs01/star/pwg/eshulga/Files/skim_from_MuDST/B15A5902FE4781FFE5DF2FC5E4298B08_12066.histograms.root";
  AnalysisCode->SetInputFileName(InputFileList);
  AnalysisCode->SetOutputFileName(OutputFileName);
  AnalysisCode->Make();
  /*
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
  */
  // Cleanup
  //delete chain ;

}

