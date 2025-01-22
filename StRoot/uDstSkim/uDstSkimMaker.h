// Based on the MuDST tools written by Frank Laue.
// Based on the DST Tutorial by Dan Magestro on the STAR Computing/Tutorials page.
// Updated 9/4/2006 by Jim Thomas to include the latest DST format, and Scheduler techniques.

#ifndef uDstSkimMaker_def
#define uDstSkimMaker_def

#include "StMaker.h"
#include "TString.h"
// my additions:
#include "TTree.h"
#include "StMuDSTMaker/COMMON/StMuEvent.h"

class StMuDstMaker ;
class TFile        ;
class TH1F         ;
// my additions:
class TTree        ;

#define MaxNumberOfTH1F      10

class uDstSkimMaker : public StMaker
{
  
 private:

  StMuDstMaker* mMuDstMaker ;                      //  Make MuDst pointer available to member functions

  TH1F*         histogram[MaxNumberOfTH1F] ;       //  1D Histograms
  TFile*        root_output ;                 //  Histograms outputfile pointer

  UInt_t        mEventsProcessed ;                 //  Number of Events read and processed
  TString       mRootOutputFileName ;         //  Name of the histogram output file 

  // start my stuff header private:

  // the tree:
  TTree* T;

  // end my stuff header private


 protected:


 public:

  uDstSkimMaker(StMuDstMaker* maker) ;       //  Constructor
  virtual          ~uDstSkimMaker( ) ;       //  Destructor

  Int_t Init    ( ) ;                              //  Initiliaze the analysis tools ... done once
  Int_t Make    ( ) ;                              //  The main analysis that is done on each event
  Int_t Finish  ( ) ;                              //  Finish the analysis, close files, and clean up.
  Int_t TPL     ( ) ;                              //  Call T->Fill() w/o X_br_fill, for template maker

  void SetOutputFileName(TString name) {mRootOutputFileName = name;} // Make name available to member functions

  void event_br_init();
  void event_br_fill();
  void MC_br_init();
  void MC_br_fill();
  void emccl_br_init();
  void emccl_br_fill();
  void eemc_clusters_init();
  void eemc_clusters_fill();
  void fdtrk_br_init();
  void fdtrk_br_fill();
  void fdvtx_br_init();
  void fdvtx_br_fill();
  void fdpair_br_init();
  void fdpair_br_fill();
  void rptrk_br_init();
  void rptrk_br_fill();
  
  ClassDef(uDstSkimMaker,1)                  //  Macro for CINT compatability
    
};

#endif

