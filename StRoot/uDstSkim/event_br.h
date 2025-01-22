// event branch header

#ifndef event_br_h
#define event_br_h

// tree general event vars:
Int_t runN;
Int_t eventN;

Int_t filln;
Int_t bid;
Int_t bid7;

Double_t bfield;

// tree trigger bits:
// generic monitors
Int_t trig_Zerobias;
Int_t trig_Zdcmon;
Int_t trig_Bbcmon;
// Run10-14 vintage
Int_t trig_UPCmain;
Int_t trig_UPCtopo;
Int_t trig_UPChighG;
Int_t trig_UPCjpsiB;
// Run15 new
Int_t trig_RP2E;
Int_t trig_2E;
Int_t trig_RPUPC;
// Run16 new
Int_t trig_UPCJPsi;
Int_t trig_UPCJPsizdc;
Int_t trig_UPCinc;
// Run17 new
Int_t trig_JPsiHTTP;

// last DSM words; DSMs: TOF&RP, BBC&ZDC, EMC
Int_t lastDSM_TOFRP;
Int_t lastDSM_BBCZDC;
Int_t lastDSM_EMC;

// ZDC info:
Int_t zdce;
Int_t zdcw;
Int_t zdceadc[3];
Int_t zdcwadc[3];
Int_t zdcetdc;
Int_t zdcwtdc;
Int_t zdctdcdiff;

// BBC, TOFtrig info:
static Int_t bbce, bbcw;
static Int_t ntoftrig;

// some rates
Double_t zdcerate;
Double_t zdcwrate;
Double_t zdccrate;

// some global event quantities:
Int_t nglobtrk;
Int_t nprimtrk;
Int_t nvtx;


// event branch code

void uDstSkimMaker::event_br_init()
{

  // tree general event vars:
  T->Branch("runN",&runN,"runN/I");
  T->Branch("eventN",&eventN,"eventN/I");

  T->Branch("filln",&filln,"filln/I");
  T->Branch("bid",&bid,"bid/I");
  T->Branch("bid7",&bid7,"bid7/I");

  T->Branch("bfield",&bfield,"bfield/D");

  // tree trigger bits:
  // generic monitors
  T->Branch("trig_Zerobias",&trig_Zerobias,"trig_Zerobias/I");
  T->Branch("trig_Zdcmon",&trig_Zdcmon,"trig_Zdcmon/I");
  T->Branch("trig_Bbcmon",&trig_Bbcmon,"trig_Bbcmon/I");
  // Run10-14 vintage
  T->Branch("trig_UPCmain",&trig_UPCmain,"trig_UPCmain/I");
  T->Branch("trig_UPCtopo",&trig_UPCtopo,"trig_UPCtopo/I");
  T->Branch("trig_UPChighG",&trig_UPChighG,"trig_UPChighG/I");
  T->Branch("trig_UPCjpsiB",&trig_UPCjpsiB,"trig_UPCjpsiB/I");
  // Run15 new
  T->Branch("trig_RP2E",&trig_RP2E,"trig_RP2E/I");
  T->Branch("trig_2E",&trig_2E,"trig_2E/I");
  T->Branch("trig_RPUPC",&trig_RPUPC,"trig_RPUPC/I");
  // Run16 new
  T->Branch("trig_UPCJPsi",&trig_UPCJPsi,"trig_UPCJPsi/I");
  T->Branch("trig_UPCJPsizdc",&trig_UPCJPsizdc,"trig_UPCJPsizdc/I");
  T->Branch("trig_UPCinc",&trig_UPCinc,"trig_UPCinc/I");
  // Run17 new
  T->Branch("trig_JPsiHTTP",&trig_JPsiHTTP,"trig_JPsiHTTP/I");

  // last DSM words; DSMs: TOF&RP, BBC&ZDC, EMC
  T->Branch("lastDSM_TOFRP",&lastDSM_TOFRP,"lastDSM_TOFRP/I");
  T->Branch("lastDSM_BBCZDC",&lastDSM_BBCZDC,"lastDSM_BBCZDC/I");
  T->Branch("lastDSM_EMC",&lastDSM_EMC,"lastDSM_EMC/I");

  // ZDC info:
  T->Branch("zdce",&zdce,"zdce/I");
  T->Branch("zdcw",&zdcw,"zdcw/I");
  T->Branch("zdceadc",&zdceadc,"zdceadc[3]/I");
  T->Branch("zdcwadc",&zdcwadc,"zdcwadc[3]/I");
  T->Branch("zdcetdc",&zdcetdc,"zdcetdc/I");
  T->Branch("zdcwtdc",&zdcwtdc,"zdcwtdc/I");
  T->Branch("zdctdcdiff",&zdctdcdiff,"zdctdcdiff/I");

  // BBC, TOFtrig info:
  T->Branch("bbce",&bbce,"bbce/I");
  T->Branch("bbcw",&bbcw,"bbcw/I");
  T->Branch("ntoftrig",&ntoftrig,"ntoftrig/I");

  // some rates
  T->Branch("zdcerate",&zdcerate,"zdcerate/D");
  T->Branch("zdcwrate",&zdcwrate,"zdcwrate/D");
  T->Branch("zdccrate",&zdccrate,"zdccrate/D");

  // some global event quantities:
  T->Branch("nglobtrk",&nglobtrk,"nglobtrk/I");
  T->Branch("nprimtrk",&nprimtrk,"nprimtrk/I");
  T->Branch("nvtx",&nvtx,"nvtx/I");

}

void uDstSkimMaker::event_br_fill()
{

  // defaults:
  runN = -9999; eventN = -9999;
  filln = -9999; bid = -9999; bid7 = -9999;
  bfield = -9999;

  trig_Zerobias = -9999; trig_Zdcmon = -9999; trig_Bbcmon = -9999;
  trig_UPCmain = -9999; trig_UPCtopo = -9999; trig_UPChighG = -9999; trig_UPCjpsiB = -9999;
  trig_RP2E = -9999; trig_2E = -9999; trig_RPUPC = -9999;
  trig_UPCJPsi = -9999; trig_UPCJPsizdc = -9999; trig_UPCinc = -9999;
  trig_JPsiHTTP = -9999;

  lastDSM_TOFRP = -9999; lastDSM_BBCZDC = -9999; lastDSM_EMC = -9999;
  zdce = -9999.; zdcw = -9999.;
  zdceadc[0] =  -9999.; zdceadc[1] =  -9999.; zdceadc[2] =  -9999.;
  zdcwadc[0] =  -9999.; zdcwadc[1] =  -9999.; zdcwadc[2] =  -9999.;
  zdcetdc =  -9999.;  zdcwtdc =  -9999.; zdctdcdiff =  -9999.;
  bbce =  -9999.; bbcw =  -9999.;
  ntoftrig =  -9999.;
  zdcerate = -9999.; zdcwrate = -9999.; zdccrate = -9999.;
  nglobtrk = -9999; nprimtrk = -9999; nvtx = -9999;

  StMuEvent* muEvent = mMuDstMaker->muDst()->event();
  if (!muEvent) {
    cout << "WARNING: event_br_fill no muEvent" << endl;
    return;
  }

  // tree general event vars:
  runN = muEvent->runNumber();
  eventN = muEvent->eventNumber();

  double fy = muEvent->runInfo().beamFillNumber(yellow);
  double fb = muEvent->runInfo().beamFillNumber(blue);
  if (fy != fb) {cout << "WARNING: event_br_fill no muEvent" << endl;}
  filln = (int)fy;

  // these seem to be cumulative bunch #, adding through revolutions
  // bcn1 seems to be higher order bit(s) of unsigned int bcn0
  // int bcn0 = muEvent->eventInfo().bunchCrossingNumber(0);
  // int bcn1 = muEvent->eventInfo().bunchCrossingNumber(1);
  // cout << "bcn0: " << bcn0  << " bcn1: " << bcn1 << endl;

  bid = muEvent->l0Trigger().bunchCrossingId();
  bid7 = muEvent->l0Trigger().bunchCrossingId7bit(runN);

  bfield = muEvent->eventSummary().magneticField();

  // Trigger bits:
  StTriggerId ttid = muEvent->triggerIdCollection().nominal();
  // each trig. Tid & run ranges from Jamie's 'lum_pertriggerid' files; checked OK with lum_perrun files
  trig_Zerobias = (
		     ((runN>=16128030 && runN<=16159024) && // Run15 pAu RP_Zerobias
		      (ttid.isTrigger(500712))
		      ) ||
		     ((runN>=16160033 && runN<=16169094) && // Run15 pAu RP_Zerobias
		      (ttid.isTrigger(510712))
		      ) ||
		     ((runN>=17038047 && runN<=17179012) && // Run16 AuAu200 Zerobias
		      (ttid.isTrigger(9300))
		      ) ||
		     ((runN>=18000000 && runN<=19000000) && // Run17 RP_Zerobias
		      (ttid.isTrigger(570704))
		      )
		   );
  trig_Zdcmon = 0; // placeholder for now
  trig_Bbcmon = 0; // placeholder for now
  trig_UPCmain = (
		  ((runN<11039046) && // Run10 AuAu pre-PHYSICS, exact run#s?
		 (ttid.isTrigger(1))
		   ) ||
		  ((runN>=11039046 && runN<=11077018) && // Run10 AuAu PHYSICS
		 (ttid.isTrigger(260750))
		   ) ||
		  ((runN<12146003) && // Run11 AuAu pre-PHYSICS, exact run#s?
		 (ttid.isTrigger(4))
		   ) ||
		  ((runN>=12146003 && runN<=12171017) && // Run11 AuAu PHYSICS
		 (ttid.isTrigger(350007) || ttid.isTrigger(350017))
		   ) ||
		  ((runN>=13139060 && runN<=13177009) && // Run12 CuAu PHYSICS
		 (ttid.isTrigger(410601) || ttid.isTrigger(410605))
		   ) ||
		  ((runN>=15078073 && runN<=15167014) && // Run14 AuAu PHYSICS
		 (ttid.isTrigger(450701) || ttid.isTrigger(450711))
		   ) ||
		  ((runN>=17041033 && runN<=17043045) && // Run16 AuAu PHYSICS
		   (ttid.isTrigger(520701))
		   )
		  );
  trig_UPCtopo = 0; // placeholder for now
  trig_UPChighG = (
		   ((runN>=17070048 && runN<=17179012) && // Run16 AuAu 200
		    (ttid.isTrigger(520741) || ttid.isTrigger(520751) ||
		     ttid.isTrigger(570701))
		    )
		   );
  trig_UPCjpsiB = (
		   ((runN>=15084052 && runN<=15167014) && // Run14 AuAu PHYSICS
		    (ttid.isTrigger(450705) || ttid.isTrigger(450725))
		    )
		   );
  trig_RP2E = (
	       ((runN>=16058079 && runN<=16117019) && // Run15 pp
		 (ttid.isTrigger(470710) || ttid.isTrigger(470730) ||
		  ttid.isTrigger(480710) || ttid.isTrigger(480730) || 
		  ttid.isTrigger(490710))
		) ||
	       ((runN>=16138044 && runN<=16148022) && // Run15 pAu
		 (ttid.isTrigger(500710))
		)
	       );
  trig_2E = (
	     ((runN>=16149001 && runN<=16159024) && // Run15 pAu
	      (ttid.isTrigger(500730) || ttid.isTrigger(500750))
	      ) ||
	     ((runN>=16160033 && runN<=16169094) && // Run15 pAl
	      (ttid.isTrigger(510710))
	      )
	     );
  trig_RPUPC = (
		((runN>=16128030 && runN<=16159024) && // Run15 pAu
		 (ttid.isTrigger(500000) || ttid.isTrigger(500020))
		 ) ||
		((runN>=16160033 && runN<=16169094) && // Run15 pAl
		 (ttid.isTrigger(510714))
		 ) ||
		((runN>=18000000 && runN<=19000000) && // Run17 pp
		 (ttid.isTrigger(570702) || ttid.isTrigger(570712))
		 )
		);
  trig_UPCJPsi = (
		  ((runN>=17041033 && runN<=17179012) && // Run16 AuAu 200
		   (ttid.isTrigger(520703) || ttid.isTrigger(520713) ||
		    ttid.isTrigger(520723) || ttid.isTrigger(520733) ||
		    ttid.isTrigger(570703))
		   )
		  );
  trig_UPCJPsizdc = (
		     ((runN>=17049051 && runN<=17179012) && // Run16 AuAu 200
		      (ttid.isTrigger(520722) || ttid.isTrigger(520732) ||
		       ttid.isTrigger(570702))
		      )
		     );
  trig_UPCinc = (
		 ((runN>=17043046 && runN<=17070019) && // Run16 AuAu 200
		  (ttid.isTrigger(520711) || ttid.isTrigger(520721) ||
		   ttid.isTrigger(520731))
		  )
		 );
  trig_JPsiHTTP = (
		   ((runN>=18000000 && runN<=19000000) && // Run17 pp
		    (ttid.isTrigger(570209) || ttid.isTrigger(570219) || ttid.isTrigger(570229))
		    )
		   );

  // some trigger info
  // trigger data object
  //const StTriggerData* trigData = muEvent->triggerData(); // not clear why not this???
  const StTriggerData * trigData;
  trigData =  const_cast< StTriggerData *> (muEvent->triggerData());
  if ((!trigData) && (runN >= 11000000)) {cout << "WARNING: event_br_fill no trigData" << endl;}
  if (trigData) {

    // last DSM words from DSMs: TOF&RP, BBC&ZDC, EMC;
    // DSM->channel assignment seems stable since ~Run9; CHECK AGAIN
    lastDSM_TOFRP = trigData->lastDSM(0);
    lastDSM_BBCZDC = trigData->lastDSM(1);
    lastDSM_EMC = trigData->lastDSM(3);

    // ZDC
    zdce = trigData->zdcUnAttenuated(east);
    zdcw = trigData->zdcUnAttenuated(west);
    for (int ipmt=0; ipmt<3; ipmt++) {
      zdceadc[ipmt] = trigData->zdcADC(east,ipmt+1); // guessed from StEvent/StZdcTriggerDetector.cxx
      zdcwadc[ipmt] = trigData->zdcADC(west,ipmt+1);
    }
    zdcetdc = trigData->zdcTDC(east);
    zdcwtdc = trigData->zdcTDC(west);
    zdctdcdiff = trigData->zdcTimeDifference();

    // BBC
    // truncated sum over BBC small tile chans. 1-16
    bbce = 0;  bbcw = 0;
    // BBC TAC,ADC cuts from Run10; CHECK!
    //int bbcTACmin = 100; int bbcTACmax = 2300; int bbcADCmin = 20;
    // BBC TAC,ADC cuts: TACmin,ADCmin from C. Perkins email 19.02.16, TACmax from trig. web page
    int bbcTACmin = 100; int bbcTACmax = 2400; int bbcADCmin = 20;
    for (int i=0; i<16; i++) { // trunc. sum over small tiles
      //cout << "BBCe " << i <<" "<< trigData->bbcTDC(east,i+1) <<" "<< trigData->bbcADC(east,i+1) << endl;
      //cout << "BBCw " << i <<" "<< trigData->bbcTDC(west,i+1) <<" "<< trigData->bbcADC(west,i+1) << endl;
      if ( trigData->bbcTDC(east,i+1)>bbcTACmin &&
           trigData->bbcTDC(east,i+1)<bbcTACmax &&
           trigData->bbcADC(east,i+1)>bbcADCmin ) {bbce += trigData->bbcADC(east,i+1);}
      if ( trigData->bbcTDC(west,i+1)>bbcTACmin &&
           trigData->bbcTDC(west,i+1)<bbcTACmax &&
           trigData->bbcADC(west,i+1)>bbcADCmin ) {bbcw += trigData->bbcADC(west,i+1);}
    }

    // TOF trig. info:
    ntoftrig =trigData->tofMultiplicity();

  }

  // some rates
  zdcerate = muEvent->runInfo().zdcEastRate();
  zdcwrate = muEvent->runInfo().zdcWestRate();
  zdccrate = muEvent->runInfo().zdcCoincidenceRate();

  // some global event quantities:
  nglobtrk = muEvent->eventSummary().numberOfGoodTracks();
  nprimtrk = muEvent->eventSummary().numberOfGoodPrimaryTracks();
  nvtx = mMuDstMaker->muDst()->numberOfPrimaryVertices();

  //cout << "BBCsums " << bbce <<" "<< bbcw << endl;

}

#endif

// spin states for this event, for tree
//Int_t sblu, syel;

/*
// spin pattern tables & info
const Int_t maxnfilltab = 320; // Run15 had 301 fills in polar data
Int_t nfillstab;
Int_t fillstab[maxnfilltab];
Int_t sblustab[maxnfilltab][120];
Int_t syelstab[maxnfilltab][120];
Int_t stab_lastfill = -9999;
Int_t stab_fillindex = -9999;
Int_t bid7_bluoffset = 0;
Int_t bid7_yeloffset = 80;

>>> TREE BRANCHES:

Te->Branch("sblu",&sblu,"sblu/I");
Te->Branch("syel",&syel,"syel/I");

>>> FILL TABLES BEFORE EVENT PROCESSING:

// read in spin pattern tables
int spinyr = 15; // non-spin years will have 2015 patterns
if (tag == 17) {spinyr = 17;}
sprintf(cvar,"/star/u/wschmidk/polana/spin_patterns/dat/spin_patterns_%i.dat",spinyr);
ifstream ifile(cvar);
nfillstab = 0;
while (1) {
  ifile >> fillstab[nfillstab];
  for (int ib=0; ib<120; ib++) {ifile >> sblustab[nfillstab][ib];}
  for (int ib=0; ib<120; ib++) {ifile >> syelstab[nfillstab][ib];}
  if (ifile.eof()) {break;}
  nfillstab++;
 }
ifile.close();
cout << "nfillstab: " << nfillstab << endl;

// check tables
// ofstream* ofile = new ofstream("spin_patterns_15.dat");
// for (int ifill = 0; ifill<nfillstab; ifill++) {
//   *ofile << fillstab[ifill] << endl;
//   for (int ib=0; ib<120; ib++) {*ofile <<" "<< sblustab[ifill][ib];} *ofile << endl;
//   for (int ib=0; ib<120; ib++) {*ofile <<" "<< syelstab[ifill][ib];} *ofile << endl;
// }
// ofile->close();

>>> EACH EVENT GET SPIN STATES:

  // get the spins
  if (filln != stab_lastfill) { // find new fill index
    stab_fillindex = -9999;
    for (int ifill=0; ifill<nfillstab; ifill++) {
      if (fillstab[ifill] == filln) {
        stab_lastfill = filln;
        stab_fillindex = ifill;
        break;
      }
    }
  }
  if (stab_fillindex<0) {
    if (yr==15) {cout << "WARNING: filln not in spin pattern table" << endl;}
    sblu = -9999;
    syel = -9999;
  }
  else {
    sblu = sblustab[stab_fillindex][(bid7+bid7_bluoffset)%120];
    syel = syelstab[stab_fillindex][(bid7+bid7_yeloffset)%120];
  }
  */
