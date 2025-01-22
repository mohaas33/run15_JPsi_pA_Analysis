// Roman Pot track branch header

#ifndef rptrk_br_h
#define rptrk_br_h


// RP tree vars
const Int_t max_rptrk = 100;
Int_t ntot_rptrk; // total # RP tracks
Int_t n_rptrk; // # RP tracks limted <= max_rptrk (array size)
Int_t br_rptrk[max_rptrk]; // branch: EU=0, ED=1, WU=2, WD=3
Int_t typ_rptrk[max_rptrk]; // my track type: 1=1st station, 2=2nd, 3=both
Int_t np_rptrk[max_rptrk]; // total # points

Int_t npst1_rptrk[max_rptrk]; // # points station 1
Int_t qualst1_rptrk[max_rptrk]; // quality station 1: 0=rpsNormal, 1=rpsGolden, 2=rpsNotSet
Double_t xst1_rptrk[max_rptrk]; // X station 1                     ^= single cluster all 4 planes 
Double_t yst1_rptrk[max_rptrk]; // Y station 1                     ^= single cluster all 4 planes 
Double_t t1st1_rptrk[max_rptrk]; // t PMT1 station 1 
Double_t t2st1_rptrk[max_rptrk]; // t PMT2 station 1

Int_t npst2_rptrk[max_rptrk];
Int_t qualst2_rptrk[max_rptrk];
Double_t xst2_rptrk[max_rptrk];
Double_t yst2_rptrk[max_rptrk]; 
Double_t t1st2_rptrk[max_rptrk];
Double_t t2st2_rptrk[max_rptrk];

Double_t pt_rptrk[max_rptrk]; // Pt
Double_t pz_rptrk[max_rptrk]; // Pz
Double_t phi_rptrk[max_rptrk]; // phi

Double_t t_rptrk[max_rptrk]; // t
Double_t xi_rptrk[max_rptrk]; // xi = 1-xL


// Roman Pot track branch code

void uDstSkimMaker::rptrk_br_init()
{

  // rptrk vars:
  T->Branch("ntot_rptrk",&ntot_rptrk,"ntot_rptrk/I");
  T->Branch("n_rptrk",&n_rptrk,"n_rptrk/I");
  T->Branch("br_rptrk",&br_rptrk,"br_rptrk[n_rptrk]/I");
  T->Branch("typ_rptrk",&typ_rptrk,"typ_rptrk[n_rptrk]/I");
  T->Branch("np_rptrk",&np_rptrk,"np_rptrk[n_rptrk]/I");

  T->Branch("npst1_rptrk",&npst1_rptrk,"npst1_rptrk[n_rptrk]/I");
  T->Branch("qualst1_rptrk",&qualst1_rptrk,"qualst1_rptrk[n_rptrk]/I");
  T->Branch("xst1_rptrk",&xst1_rptrk,"xst1_rptrk[n_rptrk]/D");
  T->Branch("yst1_rptrk",&yst1_rptrk,"yst1_rptrk[n_rptrk]/D");
  T->Branch("t1st1_rptrk",&t1st1_rptrk,"t1st1_rptrk[n_rptrk]/D");
  T->Branch("t2st1_rptrk",&t2st1_rptrk,"t2st1_rptrk[n_rptrk]/D");

  T->Branch("npst2_rptrk",&npst2_rptrk,"npst2_rptrk[n_rptrk]/I");
  T->Branch("qualst2_rptrk",&qualst2_rptrk,"qualst2_rptrk[n_rptrk]/I");
  T->Branch("xst2_rptrk",&xst2_rptrk,"xst2_rptrk[n_rptrk]/D");
  T->Branch("yst2_rptrk",&yst2_rptrk,"yst2_rptrk[n_rptrk]/D");
  T->Branch("t1st2_rptrk",&t1st2_rptrk,"t1st2_rptrk[n_rptrk]/D");
  T->Branch("t2st2_rptrk",&t2st2_rptrk,"t2st2_rptrk[n_rptrk]/D");

  T->Branch("pt_rptrk",&pt_rptrk,"pt_rptrk[n_rptrk]/D");
  T->Branch("pz_rptrk",&pz_rptrk,"pz_rptrk[n_rptrk]/D");
  T->Branch("phi_rptrk",&phi_rptrk,"phi_rptrk[n_rptrk]/D");

  T->Branch("t_rptrk",&t_rptrk,"t_rptrk[n_rptrk]/D");
  T->Branch("xi_rptrk",&xi_rptrk,"xi_rptrk[n_rptrk]/D");

  // max. array for template
  n_rptrk = max_rptrk;

}

void uDstSkimMaker::rptrk_br_fill()
{

  // defaults:
  ntot_rptrk = 0;
  n_rptrk = 0;

  StMuRpsCollection* pp2pp = mMuDstMaker->muDst()->RpsCollection();
  if (!pp2pp) {
    Int_t yr = (Int_t)(mMuDstMaker->muDst()->event()->runNumber()/1000000) - 1;
    if (yr==15 || yr==17) {cout << "WARNING: rptrk_br_fill no StMuRpsCollection" << endl;}
    return;
  }

  ntot_rptrk = pp2pp->numberOfTracks();
  if (ntot_rptrk <= max_rptrk) {
    n_rptrk = ntot_rptrk;
  }
  else {
    cout <<" WARNING: rptrk_br_fill ntot_rptrk limit exceeded" << endl;
    n_rptrk = max_rptrk;
  }

  Double_t beamMomentum = 100.; // SHOULD GET FROM StRunInfo OR SUCH 

  for (int i=0; i<n_rptrk; i++) { // start loop over tracks

    StMuRpsTrack* ptrack = pp2pp->track(i);
    if (!ptrack) {
      cout << "WARNING: rptrk_br_fill missing StMuRpsTrack" << endl;

      pt_rptrk[i] = -9999.;
      pz_rptrk[i] = -9999.;
      phi_rptrk[i] = -9999.;
      t_rptrk[i] = -9999.;
      xi_rptrk[i] = -9999.;

    }
    else { // start track OK

      br_rptrk[i] = ptrack->branch(); // branch: (0,1,2,3) = (EU,ED,WU,WD)

      // my track type: 1 = first station, 2 = second station, 3 = both stations
      //                    rpsLocal           rpsLocal            rpsGlobal
      typ_rptrk[i] = 0;
      for (int ist=0; ist<=1; ist++) {
	if (ptrack->trackPoint(ist) != 0) {typ_rptrk[i] += ist+1;}
      }

      np_rptrk[i] = ptrack->planesUsed(); // total # planes used (0-8)

      //const StMuRpsTrackPoint* point1 = ptrack->trackPoint(0); // couldn't use this correctly???
      if (!ptrack->trackPoint(0)) {
	npst1_rptrk[i] = 0;
	qualst1_rptrk[i] = 2;
	xst1_rptrk[i] = -9999.;
	yst1_rptrk[i] = -9999.;
	t1st1_rptrk[i] = -9999.;
	t2st1_rptrk[i] = -9999.;
      }
      else {
	npst1_rptrk[i] = ptrack->trackPoint(0)->planesUsed();
	qualst1_rptrk[i] = ptrack->trackPoint(0)->quality();
	xst1_rptrk[i] = ptrack->trackPoint(0)->x();
	yst1_rptrk[i] = ptrack->trackPoint(0)->y();
	t1st1_rptrk[i] = ptrack->trackPoint(0)->time(0);
	t2st1_rptrk[i] = ptrack->trackPoint(0)->time(1);
      }
      if (!ptrack->trackPoint(1)) {
	npst2_rptrk[i] = 0;
	qualst2_rptrk[i] = 2;
	xst2_rptrk[i] = -9999.;
	yst2_rptrk[i] = -9999.;
	t1st2_rptrk[i] = -9999.;
	t2st2_rptrk[i] = -9999.;
      }
      else {
	npst2_rptrk[i] = ptrack->trackPoint(1)->planesUsed();
	qualst2_rptrk[i] = ptrack->trackPoint(1)->quality();
	xst2_rptrk[i] = ptrack->trackPoint(1)->x();
	yst2_rptrk[i] = ptrack->trackPoint(1)->y();
	t1st2_rptrk[i] = ptrack->trackPoint(1)->time(0);
	t2st2_rptrk[i] = ptrack->trackPoint(1)->time(1);
      }

      pt_rptrk[i] = ptrack->pt();
      pz_rptrk[i] = ptrack->pVec().z();
      phi_rptrk[i] = ptrack->phi();
      t_rptrk[i] = ptrack->t(beamMomentum);
      xi_rptrk[i] = ptrack->xi(beamMomentum);

    } // end track OK

  } // end loop over tracks

}

#endif
