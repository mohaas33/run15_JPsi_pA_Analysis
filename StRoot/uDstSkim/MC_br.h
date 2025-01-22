// MC vertex & track branch header

#ifndef MC_br_h
#define MC_br_h

// MC vertex tree vars:
const Int_t max_MCvtx = 1000;
Int_t n_MCvtx;
Double_t x_MCvtx[max_MCvtx];
Double_t y_MCvtx[max_MCvtx];
Double_t z_MCvtx[max_MCvtx];

// MC track tree vars:
const Int_t max_MCtrk = 1000;
Int_t n_MCtrk;
Int_t iMCvtx_MCtrk[max_MCtrk]; // index in MCvtx list
Int_t pid_MCtrk[max_MCtrk];
Int_t q_MCtrk[max_MCtrk];
Double_t px_MCtrk[max_MCtrk];
Double_t py_MCtrk[max_MCtrk];
Double_t pz_MCtrk[max_MCtrk];

// MC vertex & track branch code

void uDstSkimMaker::MC_br_init()
{

  // treeMC vertex vars:
  T->Branch("n_MCvtx",&n_MCvtx,"n_MCvtx/I");
  T->Branch("x_MCvtx",&x_MCvtx,"x_MCvtx[n_MCvtx]/D");
  T->Branch("y_MCvtx",&y_MCvtx,"y_MCvtx[n_MCvtx]/D");
  T->Branch("z_MCvtx",&z_MCvtx,"z_MCvtx[n_MCvtx]/D");

  // treeMC track vars:
  T->Branch("n_MCtrk",&n_MCtrk,"n_MCtrk/I");
  T->Branch("iMCvtx_MCtrk",&iMCvtx_MCtrk,"iMCvtx_MCtrk[n_MCtrk]/I");
  T->Branch("pid_MCtrk",&pid_MCtrk,"pid_MCtrk[n_MCtrk]/I");
  T->Branch("q_MCtrk",&q_MCtrk,"q_MCtrk[n_MCtrk]/I");
  T->Branch("px_MCtrk",&px_MCtrk,"px_MCtrk[n_MCtrk]/D");
  T->Branch("py_MCtrk",&py_MCtrk,"py_MCtrk[n_MCtrk]/D");
  T->Branch("pz_MCtrk",&pz_MCtrk,"pz_MCtrk[n_MCtrk]/D");

  // max. array for template
  n_MCvtx = max_MCvtx;
  n_MCtrk = max_MCtrk;

}

void uDstSkimMaker::MC_br_fill()
{

  // defaults:
  n_MCvtx = 0;
  n_MCtrk = 0;

  TClonesArray *MuMcVertices = mMuDstMaker->muDst()->mcArray(0);
  TClonesArray *MuMcTracks = mMuDstMaker->muDst()->mcArray(1);

  n_MCvtx = MuMcVertices->GetEntriesFast();
  n_MCtrk = MuMcTracks->GetEntriesFast();

  if (n_MCvtx > max_MCvtx) {
    //cout << " WARNING: MC_br_fill n_MCvtx limit exceeded" << endl;
    n_MCvtx = max_MCvtx;
  }
  if (n_MCtrk > max_MCtrk) {
    //cout << " WARNING: MC_br_fill n_MCtrk limit exceeded" << endl;
    n_MCtrk = max_MCtrk;
  }

  for (int ivtx=0; ivtx<n_MCvtx; ivtx++) { // start loop over MC vertices

    StMuMcVertex *mcVertex = (StMuMcVertex*) MuMcVertices->UncheckedAt(ivtx);

    x_MCvtx[ivtx] = mcVertex->XyzV().x();
    y_MCvtx[ivtx] = mcVertex->XyzV().y();
    z_MCvtx[ivtx] = mcVertex->XyzV().z();

  } // end loop over MC vertices

  for (int itrk=0; itrk<n_MCtrk; itrk++) { // start loop over MC tracks

    StMuMcTrack *mcTrack = (StMuMcTrack *) MuMcTracks->UncheckedAt(itrk);

    // index in MCvtx array = 1 (primary), 2, 3, ...; avoid overflows
    if (mcTrack->IdVx() <= max_MCvtx) {iMCvtx_MCtrk[itrk] = mcTrack->IdVx() - 1;}
    else {iMCvtx_MCtrk[itrk] = -9999;}
    pid_MCtrk[itrk] = mcTrack->GePid();
    q_MCtrk[itrk] = mcTrack->Charge();
    px_MCtrk[itrk] = mcTrack->Pxyz().x();
    py_MCtrk[itrk] = mcTrack->Pxyz().y();
    pz_MCtrk[itrk] = mcTrack->Pxyz().z();

  } // end loop over MC tracks

}

#endif
