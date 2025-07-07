// Based on the MuDST tools written by Frank Laue.
// Based on the DST Tutorial by Dan Magestro on the STAR Computing/Tutorials page.
// Updated 9/4/2006 by Jim Thomas to include the latest DST format, and Scheduler techniques.

#include "uDstSkimMaker.h"

#include <iostream>

#include "StMuDSTMaker/COMMON/StMuDstMaker.h"
#include "StMuDSTMaker/COMMON/StMuTrack.h"
#include "StMuDSTMaker/COMMON/StMuEvent.h"
// my additions:
#include "StEventTypes.h" // THIS HAS >150 #include's
// MC info
#include "StMuDSTMaker/COMMON/StMuMcVertex.h"
#include "StMuDSTMaker/COMMON/StMuMcTrack.h"
#include <StMuDSTMaker/COMMON/StMuPrimaryVertex.h>
// EMC: // no more needed
// track->BEMC extrap.:
#include "StEmcUtil/geometry/StEmcGeom.h"
#include "StEmcUtil/projection/StEmcPosition.h"

//StSpinDbMaker
#include "StSpinPool/StSpinDbMaker/StSpinDbMaker.h"

//#include "StEvent/StL0Trigger.h"


//EEMC stuff
#include "StEEmcUtil/EEmcGeom/EEmcGeomSimple.h"
#include "StEEmcUtil/database/StEEmcDb.h"   
#include "StEEmcUtil/database/EEmcDbItem.h"
#include "StEEmcPool/StEEmcA2EMaker/StEEmcA2EMaker.h"

// TOF hit:
#include "StMuDSTMaker/COMMON/StMuBTofHit.h"
// RP stuff
#include "StMuDSTMaker/COMMON/StMuRpsCollection.h"

#include "event_br.h"
#include "MC_br.h"
#include "emccl_br.h"
#include "eemc_clusters.h"
#include "fdtrk_br.h"
#include "fdvtx_br.h"
#include "fdpair_br.h"
#include "rptrk_br.h"

#include "TH1.h"
#include "TFile.h"
#include "TObjArray.h"
#include "TClonesArray.h"
#include "TList.h"
// my additions:
#include "TTree.h"

#define NumberOfTH1F      0                     // Number of Histograms

ClassImp(uDstSkimMaker)                   // Macro for CINT compatibility


uDstSkimMaker::uDstSkimMaker( StMuDstMaker* maker ) : StMaker("uDstSkimMaker")

{ // Initialize and/or zero all public/private data members here.

  for ( Int_t i = 0 ; i < NumberOfTH1F ; i++ )  // Zero the histogram pointers
    {
      histogram[i] = NULL ;
    }

  mMuDstMaker      = maker ;                    // Pass MuDst pointer to DstAnlysisMaker Class member functions
  root_output = NULL  ;                    // Zero the Pointer to histogram output file
  mEventsProcessed = 0     ;                    // Zero the Number of Events processed by the maker 
  mRootOutputFileName = "" ;               // Histogram Output File Name will be set inside the "analysis".C macro

}


uDstSkimMaker::~uDstSkimMaker() 

{ // Destroy and/or zero out all public/private data members here.

}


Int_t uDstSkimMaker::Init( )

{ // Do once at the start of the analysis

  // Create Histogram output file

  root_output = new TFile( mRootOutputFileName, "recreate" ) ;  // Name was set in "analysis".C macro

  // Create Histograms

  // the tree:
  T = new TTree("T","uDst");

  // start my stuff init:

  event_br_init();
  MC_br_init();
  emccl_br_init();
  eemc_clusters_init();
  fdtrk_br_init();
  fdvtx_br_init();
  fdpair_br_init();
  rptrk_br_init();

  // end my stuff init

  return kStOK ;

}


Int_t uDstSkimMaker::Make( )
 
{ // Do each event

  // start my stuff event:

  event_br_fill();


  //cout << mEventsProcessed <<" "<< runN <<" "<< eventN << endl;

  mEventsProcessed++ ;

  // for DATA select only these trigs; ignore for MC
  if ( (runN >= 11000000) && // assume MC for Run < 2010
       !(
  	 trig_Zerobias || trig_Zdcmon || trig_Bbcmon ||
  	 trig_UPCmain || trig_UPCtopo || trig_UPChighG || trig_UPCjpsiB ||
  	 trig_RP2E || trig_2E || trig_RPUPC ||
  	 trig_UPCJPsi || trig_UPCJPsizdc || trig_UPCinc || trig_JPsiHTTP
	  ) ) {return kStOK ;}

  // MC vertex & track info
  MC_br_fill();

  // BEMC clusters
  emccl_br_fill();

  // EEMC cluster
  eemc_clusters_fill();

  // track->FD extrap.
  fdtrk_br_fill();

  // vertices info
  fdvtx_br_fill();

  // pair info
  fdpair_br_fill();

  // RP info
  rptrk_br_fill();

  T->Fill();

  // end my stuff event
  return kStOK ;
  
}


Int_t uDstSkimMaker::Finish( )

{ // Do once at the end of the analysis

  // start my stuff finish

  // end my stuff finish

  // Write histograms to disk, output miscellaneous other information

  root_output -> Write() ;   // Write all histograms to disk 

  cout << "Total Events Processed in DstMaker " << mEventsProcessed << endl ;

  return kStOK ;

}


Int_t uDstSkimMaker::TPL( )

{ // for making template

  T->Fill(); // one entry

  return kStOK ;

}

