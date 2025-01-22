// FD pair branch header

#ifndef fdpair_br_h
#define fdpair_br_h

// pair tree vars:
const Int_t max_fdpair = 1000;
Int_t n_fdpair;
Int_t ifdtrk1_fdpair[max_fdpair]; // fdtrk list index = 0,1,...
Int_t ifdtrk2_fdpair[max_fdpair];

Int_t q_fdpair[max_fdpair];
Double_t pt_fdpair[max_fdpair];
Double_t pz_fdpair[max_fdpair];
Double_t p_fdpair[max_fdpair];
Double_t phi_fdpair[max_fdpair];
Double_t acolin_fdpair[max_fdpair]; // track p-vectors acolinearity

Double_t mee_fdpair[max_fdpair];
Double_t mmumu_fdpair[max_fdpair];
Double_t mpipi_fdpair[max_fdpair];
Double_t mkk_fdpair[max_fdpair];
Double_t mpp_fdpair[max_fdpair];
Double_t rapee_fdpair[max_fdpair];
Double_t rapmumu_fdpair[max_fdpair];
Double_t rappipi_fdpair[max_fdpair];
Double_t rapkk_fdpair[max_fdpair];
Double_t rappp_fdpair[max_fdpair];

// pair branch code

void uDstSkimMaker::fdpair_br_init()
{

  // tree pair vars:
  T->Branch("n_fdpair",&n_fdpair,"n_fdpair/I");
  T->Branch("ifdtrk1_fdpair",&ifdtrk1_fdpair,"ifdtrk1_fdpair[n_fdpair]/I");
  T->Branch("ifdtrk2_fdpair",&ifdtrk2_fdpair,"ifdtrk2_fdpair[n_fdpair]/I");

  T->Branch("q_fdpair",&q_fdpair,"q_fdpair[n_fdpair]/I");
  T->Branch("pt_fdpair",&pt_fdpair,"pt_fdpair[n_fdpair]/D");
  T->Branch("pz_fdpair",&pz_fdpair,"pz_fdpair[n_fdpair]/D");
  T->Branch("p_fdpair",&p_fdpair,"p_fdpair[n_fdpair]/D");
  T->Branch("phi_fdpair",&phi_fdpair,"phi_fdpair[n_fdpair]/D");
  T->Branch("acolin_fdpair",&acolin_fdpair,"acolin_fdpair[n_fdpair]/D");

  T->Branch("mee_fdpair",&mee_fdpair,"mee_fdpair[n_fdpair]/D");
  T->Branch("mmumu_fdpair",&mmumu_fdpair,"mmumu_fdpair[n_fdpair]/D");
  T->Branch("mpipi_fdpair",&mpipi_fdpair,"mpipi_fdpair[n_fdpair]/D");
  T->Branch("mkk_fdpair",&mkk_fdpair,"mkk_fdpair[n_fdpair]/D");
  T->Branch("mpp_fdpair",&mpp_fdpair,"mpp_fdpair[n_fdpair]/D");
  T->Branch("rapee_fdpair",&rapee_fdpair,"rapee_fdpair[n_fdpair]/D");
  T->Branch("rapmumu_fdpair",&rapmumu_fdpair,"rapmumu_fdpair[n_fdpair]/D");
  T->Branch("rappipi_fdpair",&rappipi_fdpair,"rappipi_fdpair[n_fdpair]/D");
  T->Branch("rapkk_fdpair",&rapkk_fdpair,"rapkk_fdpair[n_fdpair]/D");
  T->Branch("rappp_fdpair",&rappp_fdpair,"rappp_fdpair[n_fdpair]/D");

  // max. array for template
  n_fdpair = max_fdpair;

}

void uDstSkimMaker::fdpair_br_fill()
{

  // defaults:
  n_fdpair = 0;

  // loop over pairs of FD matched tracks
  for (int itrk1=0; itrk1<n_fdtrk-1; itrk1++) { // start loop overt trk1
    for (int itrk2=itrk1+1; itrk2<n_fdtrk; itrk2++) { // start loop overt trk2

      if (
	  ivtx_fdtrk[itrk1] == ivtx_fdtrk[itrk2] // same vertex
	  //&& iemccl_fdtrk[itrk1] != iemccl_fdtrk[itrk2] // not same emc cluster
	  )
	{ // start new pair

	  ifdtrk1_fdpair[n_fdpair] = itrk1;
	  ifdtrk2_fdpair[n_fdpair] = itrk2;

	  q_fdpair[n_fdpair] = q_fdtrk[itrk1] + q_fdtrk[itrk2];

	  // get the tracks (can do w/o TObjArray?)
	  StMuDst::setVertexIndex(ivtx_fdtrk[itrk1]); //<-BOOBOO
	  TObjArray* tracks = mMuDstMaker->muDst()->primaryTracks();
	  StMuTrack* trk1 = (StMuTrack*)tracks->At(idx_fdtrk[itrk1]);
	  StMuTrack* trk2 = (StMuTrack*)tracks->At(idx_fdtrk[itrk2]);
	  if (!tracks || !trk1 || !trk2)
	    {cout << "WARNING: fdpair_br_fill missing tracks pointer(s)" << endl;}
	  else { // start track pointers OK

	    // 3-momenta
	    StThreeVectorF p31 = trk1->p();
	    StThreeVectorF p32 = trk2->p();
	    StThreeVectorF p3pair = p31 + p32;

	    pt_fdpair[n_fdpair] = p3pair.perp();
	    pz_fdpair[n_fdpair] = p3pair.z();
	    p_fdpair[n_fdpair] = p3pair.mag();
	    phi_fdpair[n_fdpair] = p3pair.phi();

	    acolin_fdpair[n_fdpair] = TMath::Pi()-p31.angle(p32); // track p-vectors acolinearity

	    // where're the SCL constants?
	    //cout << "electron mass: " << electron_mass_c2*GeV << endl;

	    // CHECK UNIT CONVERSIONS, e.g. M_e*MeV ???
	    //const Double_t M_e = 140.; // checked p4.m() OK

	    // ee hypothesis
	    const Double_t M_e = 0.510998928; // from PDG 2014/2015 update
	    StLorentzVectorF p4e1(p31,p31.massHypothesis(M_e*MeV));
	    StLorentzVectorF p4e2(p32,p32.massHypothesis(M_e*MeV));
	    //cout << " Mtrk1: " << p4e1.m() << " Mtrk2: " << p4e2.m() << endl;
	    StLorentzVectorF p4ee = p4e1 + p4e2;
	    mee_fdpair[n_fdpair] = p4ee.m();
	    rapee_fdpair[n_fdpair] = p4ee.rapidity();

	    // mumu hypothesis
	    const Double_t M_mu = 105.6583715; // from PDG 2014/2015 update
	    StLorentzVectorF p4mu1(p31,p31.massHypothesis(M_mu*MeV));
	    StLorentzVectorF p4mu2(p32,p32.massHypothesis(M_mu*MeV));
	    StLorentzVectorF p4mumu = p4mu1 + p4mu2;
	    mmumu_fdpair[n_fdpair] = p4mumu.m();
	    rapmumu_fdpair[n_fdpair] = p4mumu.rapidity();

	    // pipi hypothesis
	    const Double_t M_pi = 139.57018; // from PDG 2014/2015 update
	    StLorentzVectorF p4pi1(p31,p31.massHypothesis(M_pi*MeV));
	    StLorentzVectorF p4pi2(p32,p32.massHypothesis(M_pi*MeV));
	    StLorentzVectorF p4pipi = p4pi1 + p4pi2;
	    mpipi_fdpair[n_fdpair] = p4pipi.m();
	    rappipi_fdpair[n_fdpair] = p4pipi.rapidity();

	    // kk hypothesis
	    const Double_t M_k = 493.677; // from PDG 2014/2015 update
	    StLorentzVectorF p4k1(p31,p31.massHypothesis(M_k*MeV));
	    StLorentzVectorF p4k2(p32,p32.massHypothesis(M_k*MeV));
	    StLorentzVectorF p4kk = p4k1 + p4k2;
	    mkk_fdpair[n_fdpair] = p4kk.m();
	    rapkk_fdpair[n_fdpair] = p4kk.rapidity();

	    // pp hypothesis
	    const Double_t M_p = 938.272046; // from PDG 2014/2015 update
	    StLorentzVectorF p4p1(p31,p31.massHypothesis(M_p*MeV));
	    StLorentzVectorF p4p2(p32,p32.massHypothesis(M_p*MeV));
	    StLorentzVectorF p4pp = p4p1 + p4p2;
	    mpp_fdpair[n_fdpair] = p4pp.m();
	    rappp_fdpair[n_fdpair] = p4pp.rapidity();

	    n_fdpair++;
	    if (n_fdpair >= max_fdpair) {
	      cout <<" WARNING: fdpair_br_fill n_fdpair limit reached" << endl;
	      return;
	    }

	  } // end track pointers OK

	} // end new pair

    } // end loop overt trk2
  } // end loop overt trk1

}

#endif
