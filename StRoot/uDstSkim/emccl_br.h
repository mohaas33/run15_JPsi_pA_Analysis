// EMC cluster branch header

#ifndef emccl_br_h
#define emccl_br_h
  
// table of cells in clusters, and energies
Int_t emccl_cells[120][2][20];
Double_t emccl_ecells[120][2][20];

// BEMC 'clusters' tree vars
const Int_t max_emccl = 4800;
Int_t n_emccl;
Double_t e_emccl[max_emccl];
Double_t eta_emccl[max_emccl];
Double_t phi_emccl[max_emccl];
Double_t sigeta_emccl[max_emccl];
Double_t sigphi_emccl[max_emccl];
Int_t nhits_emccl[max_emccl];
Int_t idxmax_emccl[max_emccl]; // index (1-4800) of max. E cell
Double_t emax_emccl[max_emccl]; // E of max. E cell
Int_t adcmax_emccl[max_emccl]; // ADC of max. ADC cell
Int_t etabinmax_emccl[max_emccl]; // max. E cell eta bin 1-40 from (-1,+1)
Int_t phibinmax_emccl[max_emccl]; // max. E cell phi bin 1-120 from (-Pi,+Pi)
Int_t phiwdgmax_emccl[max_emccl]; // max. E cell phi wedge 1-6 from (-Pi,+Pi)
//Int_t n_emcclOK;
//Int_t emcclOK[max_emccl];

int mod2phi(int mod); // prototype for map BEMC modules to phi bins 1-60
int smod2phi(int smod); // prototype for map BEMC sub modules to phi bins 1-120

// EMC cluster branch code

void uDstSkimMaker::emccl_br_init()
{

  // emccl vars:
  T->Branch("n_emccl",&n_emccl,"n_emccl/I");
  T->Branch("e_emccl",e_emccl,"e_emccl[n_emccl]/D");
  T->Branch("eta_emccl",eta_emccl,"eta_emccl[n_emccl]/D");
  T->Branch("phi_emccl",phi_emccl,"phi_emccl[n_emccl]/D");
  T->Branch("sigeta_emccl",sigeta_emccl,"sigeta_emccl[n_emccl]/D");
  T->Branch("sigphi_emccl",sigphi_emccl,"sigphi_emccl[n_emccl]/D");
  T->Branch("nhits_emccl",&nhits_emccl,"nhits_emccl[n_emccl]/I");
  T->Branch("idxmax_emccl",&idxmax_emccl,"idxmax_emccl[n_emccl]/I");
  T->Branch("emax_emccl",&emax_emccl,"emax_emccl[n_emccl]/D");
  T->Branch("adcmax_emccl",&adcmax_emccl,"adcmax_emccl[n_emccl]/I");
  T->Branch("etabinmax_emccl",&etabinmax_emccl,"etabinmax_emccl[n_emccl]/I");
  T->Branch("phibinmax_emccl",&phibinmax_emccl,"phibinmax_emccl[n_emccl]/I");
  T->Branch("phiwdgmax_emccl",&phiwdgmax_emccl,"phiwdgmax_emccl[n_emccl]/I");
  //T->Branch("n_emcclOK",&n_emcclOK,"n_emcclOK/I");
  //T->Branch("_emcclOK",_emcclOK,"_emcclOK[n_emccl]/I");

  // max. array for template
  n_emccl = max_emccl;

}

void uDstSkimMaker::emccl_br_fill()
{

  // defaults:
  n_emccl = 0;

  // clear cells used table
  for (int i=0; i<120; i++) { for(int j=0; j<2; j++) { for (int k=0; k<20; k++) {
	emccl_cells[i][j][k] = -9999;
	emccl_ecells[i][j][k] = 0.;
      }}}

  // for BEMC a la $STAR/StRoot/macros/mudst/exampleEmc.C
  // get BEMC clusters (several layers deep), get cluster info & max. E cell

  StEmcCollection *emcCollection = mMuDstMaker->muDst()->emcCollection(); // EMC

  if (!emcCollection) {cout << "WARNING: emccl_br_fill no emc collection" << endl;}
  else { // start emc collection OK

    StEmcDetector *barrel = emcCollection->detector(kBarrelEmcTowerId); // BEMC
    if(!barrel){/*cout << "WARNING: emccl_br_fill no barrel" << endl;*/}
    else { // start barrel OK

      if (!barrel->cluster()) {/*cout << "WARNING: emccl_br_fill no barrel cluster collection" << endl;*/}
      else { // start barrel cluster collection OK

	StSPtrVecEmcCluster& clusters = barrel->cluster()->clusters(); // BEMC clusters
	n_emccl = clusters.size();
	if (n_emccl > max_emccl) {
	  cout << "WARNING: emccl_br_fill n_emccl > max_emccl" << endl;
	  n_emccl = max_emccl;
	}
	for (int icl=0; icl<n_emccl; icl++) { // start loop over BEMC clusters

	  StEmcCluster* emccl = clusters[icl]; // BEMC cluster
	  if (!emccl) {
	    cout << "WARNING: emccl_br_fill no emccl" << endl;
	    e_emccl[icl] = -9999.; // flag for missing emccl
	  }
	  else { // start emccl OK

	    e_emccl[icl] = emccl->energy();
	    eta_emccl[icl] = emccl->eta();
	    phi_emccl[icl] = emccl->phi();
	    sigeta_emccl[icl] = emccl->sigmaEta();
	    sigphi_emccl[icl] = emccl->sigmaPhi();

	    StPtrVecEmcRawHit& hits = emccl->hit(); // BEMC cluster hits
	    nhits_emccl[icl] = emccl->nHits();
	    Int_t idxmax = -9999; Double_t emax=-9999.; Int_t adcmax = -9999;
	    Int_t modmax = -9999; Int_t submax = -9999; Int_t etamax = -9999; 
	    for (int ih=0; ih<nhits_emccl[icl]; ih++) { // start loop over cluster cells
	      StEmcRawHit* h = hits[ih]; // BEMC cluster hit
	      if (!h) {cout << "WARNING: emccl_br_fill no h (cluster hits)" << endl;}
	      else {
		if (emccl_cells[h->module()-1][h->sub()-1][h->eta()-1] >= 0) { // check reused cell?
		  cout << "WARNING: emccl_br_fill BEMC HIT USED > 1 CLUSTER" << endl;
		}
		emccl_cells[h->module()-1][h->sub()-1][h->eta()-1] = icl;
		emccl_ecells[h->module()-1][h->sub()-1][h->eta()-1] = h->energy();
		if (h->energy() > emax) { // get max. E cell in cluster
		  idxmax = 20*(2*(h->module()-1) + (h->sub()-1)) + h->eta();
		  emax = h->energy();
		  modmax = h->module(); submax = h->sub(); etamax = h->eta();
		}
		if ((Int_t)h->adc() > adcmax) { // get max. cell ADC in cluster
		  adcmax = (Int_t)h->adc();
		}
	      }
	    } // end loop over cluster cells
	    idxmax_emccl[icl] = idxmax;
	    emax_emccl[icl] = emax;
	    adcmax_emccl[icl] = adcmax;
	    if (modmax >= 61) {etabinmax_emccl[icl] = 21 - etamax;} // eta<0
	    else {etabinmax_emccl[icl] = 20 + etamax;} // eta>0
	    Int_t smodmax = 2*(modmax-1) + submax; // sub-mod. # 1-240
	    phibinmax_emccl[icl] = smod2phi(smodmax); // complicated mapping
	    phiwdgmax_emccl[icl] = (int)((phibinmax_emccl[icl]-1)/20) + 1; // phi wedge 1-6

	  } // end emccl OK
	} // end loop over BEMC clusters
      } // end barrel cluster collection OK
    } // end barrel OK
  } // end emc collection OK

}

int mod2phi(int mod) // map BEMC modules to phi bins 1-60
{
  int phibin = 0; // 0 if bad module #
  if (mod>= 1 && mod <= 43) {phibin =  44-mod;}
  if (mod>=44 && mod <= 60) {phibin = 104-mod;}
  if (mod>=61 && mod <= 72) {phibin =  mod-12;}
  if (mod>=73 && mod <=120) {phibin =  mod-72;}
  return phibin;
}

int smod2phi(int smod) // map BEMC sub modules to phi bins 1-120
{
  int phibin = 0; // 0 if bad sub module #
  if (smod>=  1 && smod <= 85) {phibin =  86-smod;}
  if (smod>= 86 && smod <=120) {phibin = 206-smod;}
  if (smod>=121 && smod <=145) {phibin =  smod-25;}
  if (smod>=146 && smod <=240) {phibin =  smod-145;}
  return phibin;
}

#endif
