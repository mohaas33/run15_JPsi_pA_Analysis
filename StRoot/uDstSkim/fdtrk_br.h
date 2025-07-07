// Fast Detecotr (FD) matched track branch header
// (matched == track points to a cell in BEMC cluster and/or matched to a TOF hit)

#ifndef fdtrk_br_h
#define fdtrk_br_h

// projection & geometry objects:
StEmcPosition* project;
StEmcGeom* mGeom;
EEmcGeomSimple *endcapGeom;

// Fast Detector matched track tree vars
const Int_t max_fdtrk = 4800;
Int_t n_fdtrk;
Int_t ifdvtx_fdtrk[max_fdtrk]; // prim. vtx. index = 0,1,... in fdvtx list; SET IN VTX CODE
// only for local use, not saved in tree:
Int_t ivtx_fdtrk[max_fdtrk]; // prim. vtx. index = 0,1,... in MuDst list
Int_t idx_fdtrk[max_fdtrk]; // prim. trk index in MuDst list for this vertex = 0,1,...

// track quantities
Int_t mtd_match_fdtrk[max_fdtrk];
Int_t q_fdtrk[max_fdtrk];
Double_t pt_fdtrk[max_fdtrk];
Double_t pz_fdtrk[max_fdtrk];
Double_t p_fdtrk[max_fdtrk];
Double_t eta_fdtrk[max_fdtrk];
Double_t phi_fdtrk[max_fdtrk];
Int_t nh_fdtrk[max_fdtrk];
Double_t dcar_fdtrk[max_fdtrk];
Double_t dcaz_fdtrk[max_fdtrk];

Double_t vtxz_fdtrk[max_fdtrk];
Int_t    vtxi_fdtrk[max_fdtrk];

// track dE/dx info
Double_t dedx_fdtrk[max_fdtrk];
Int_t nhdedx_fdtrk[max_fdtrk];
Double_t sigel_fdtrk[max_fdtrk];
Double_t sigpi_fdtrk[max_fdtrk];
Double_t sigk_fdtrk[max_fdtrk];
Double_t sigp_fdtrk[max_fdtrk];

// track->BEMC info
Int_t emcext_fdtrk[max_fdtrk]; // 0/1 NO/YES track extrapolates to BEMC
Double_t etaext_fdtrk[max_fdtrk]; // track extrap. eta @ BEMC
Double_t phiext_fdtrk[max_fdtrk]; // track extrap. phi @ BEMC
Int_t iemccl_fdtrk[max_fdtrk]; // <0 no BEMC match; >=0 matched emc cluster index

// track->EEMC info
Int_t eemc_ext_fdtrk[max_fdtrk]; // 0/1 NO/YES track extrapolates to EEMC
Double_t eemc_etaext_fdtrk[max_fdtrk]; // track extrap. eta @ EEMC
Double_t eemc_phiext_fdtrk[max_fdtrk]; // track extrap. phi @ EEMC
Int_t ieemc_cl_fdtrk[max_fdtrk]; // <0 no EEMC match; >=0 matched eemc tower index
Double_t eemc_adc_fdtrk[max_fdtrk]; //Associates EEMC ADC value to fdtrack
Double_t eemc_energy_fdtrk[max_fdtrk]; //Associates EEMC ADC->Energy value to fdtrack
Double_t eemc_energy_prs1_fdtrk[max_fdtrk];
Double_t eemc_energy_prs2_fdtrk[max_fdtrk];
Double_t eemc_energy_tow_fdtrk[max_fdtrk];
Double_t eemc_energy_post_fdtrk[max_fdtrk];
Double_t eemc_adc_prs1_fdtrk[max_fdtrk];
Double_t eemc_adc_prs2_fdtrk[max_fdtrk];
Double_t eemc_adc_tow_fdtrk[max_fdtrk];
Double_t eemc_adc_post_fdtrk[max_fdtrk];
Int_t eemc_num_towers_in_track[max_fdtrk];

// track TOF info
Int_t tofhit_fdtrk[max_fdtrk];
Double_t tofpatlen_fdtrk[max_fdtrk];
Double_t toftof_fdtrk[max_fdtrk];
Double_t tofbeta_fdtrk[max_fdtrk];
Double_t toflet_fdtrk[max_fdtrk];
Double_t toftet_fdtrk[max_fdtrk];

//EEMC track projection function

Bool_t trackOnEndcap(StThreeVectorD* position, StThreeVectorD* momentum, const StMuTrack* const track, double magField, double endcapDist, int &sector, int &sub, int &etabin );
Bool_t projectEEMCTrack(StThreeVectorD* atFinal, StThreeVectorD* momentumAtFinal, const StMuTrack* const track, double magField, double radius, int option);


// Fast Detector matched track branch code

void uDstSkimMaker::fdtrk_br_init()
{

  // instantiate projection & geometry:
  project = new StEmcPosition();
  mGeom = StEmcGeom::instance("bemc");
  endcapGeom = new EEmcGeomSimple();

  // fdtrk vars:
  T->Branch("n_fdtrk",&n_fdtrk,"n_fdtrk/I");
  T->Branch("ifdvtx_fdtrk",&ifdvtx_fdtrk,"ifdvtx_fdtrk[n_fdtrk]/I");
  //T->Branch("ivtx_fdtrk",&ivtx_fdtrk,"ivtx_fdtrk[n_fdtrk]/I"); // only for local use
  //T->Branch("idx_fdtrk",&idx_fdtrk,"idx_fdtrk[n_fdtrk]/I"); // only for local use

  // track quantities
  T->Branch("mtd_match_fdtrk",&mtd_match_fdtrk,"mtd_match_fdtrk[n_fdtrk]/I");
  T->Branch("q_fdtrk",&q_fdtrk,"q_fdtrk[n_fdtrk]/I");
  T->Branch("pt_fdtrk",&pt_fdtrk,"pt_fdtrk[n_fdtrk]/D");
  T->Branch("pz_fdtrk",&pz_fdtrk,"pz_fdtrk[n_fdtrk]/D");
  T->Branch("p_fdtrk",&p_fdtrk,"p_fdtrk[n_fdtrk]/D");
  T->Branch("eta_fdtrk",&eta_fdtrk,"eta_fdtrk[n_fdtrk]/D");
  T->Branch("phi_fdtrk",&phi_fdtrk,"phi_fdtrk[n_fdtrk]/D");
  T->Branch("nh_fdtrk",&nh_fdtrk,"nh_fdtrk[n_fdtrk]/I");
  T->Branch("dcar_fdtrk",&dcar_fdtrk,"dcar_fdtrk[n_fdtrk]/D");
  T->Branch("dcaz_fdtrk",&dcaz_fdtrk,"dcaz_fdtrk[n_fdtrk]/D");
  T->Branch("vtxz_fdtrk",&vtxz_fdtrk,"vtxz_fdtrk/D");
  T->Branch("vtxi_fdtrk",&vtxi_fdtrk,"vtxi_fdtrk/I");
  
  // track dE/dx info
  T->Branch("dedx_fdtrk",&dedx_fdtrk,"dedx_fdtrk[n_fdtrk]/D");
  T->Branch("nhdedx_fdtrk",&nhdedx_fdtrk,"nhdedx_fdtrk[n_fdtrk]/I");
  T->Branch("sigel_fdtrk",&sigel_fdtrk,"sigel_fdtrk[n_fdtrk]/D");
  T->Branch("sigpi_fdtrk",&sigpi_fdtrk,"sigpi_fdtrk[n_fdtrk]/D");
  T->Branch("sigk_fdtrk",&sigk_fdtrk,"sigk_fdtrk[n_fdtrk]/D");
  T->Branch("sigp_fdtrk",&sigp_fdtrk,"sigp_fdtrk[n_fdtrk]/D");

  // track->BEMC info
  T->Branch("emcext_fdtrk",&emcext_fdtrk,"emcext_fdtrk[n_fdtrk]/I");
  T->Branch("etaext_fdtrk",&etaext_fdtrk,"etaext_fdtrk[n_fdtrk]/D");
  T->Branch("phiext_fdtrk",&phiext_fdtrk,"phiext_fdtrk[n_fdtrk]/D");
  T->Branch("iemccl_fdtrk",&iemccl_fdtrk,"iemccl_fdtrk[n_fdtrk]/I");
  
  // track->EEMC info
  T->Branch("eemc_ext_fdtrk",&eemc_ext_fdtrk,"eemc_ext_fdtrk[n_fdtrk]/I");
  T->Branch("eemc_etaext_fdtrk",&eemc_etaext_fdtrk,"eemc_etaext_fdtrk[n_fdtrk]/D");
  T->Branch("eemc_phiext_fdtrk",&eemc_phiext_fdtrk,"eemc_phiext_fdtrk[n_fdtrk]/D");
  T->Branch("ieemc_cl_fdtrk",&ieemc_cl_fdtrk,"ieemc_cl_fdtrk[n_fdtrk]/I");
  T->Branch("eemc_adc_fdtrk",&eemc_adc_fdtrk,"eemc_adc_fdtrk[n_fdtrk]/D");
  T->Branch("eemc_energy_fdtrk",&eemc_energy_fdtrk,"eemc_energy_fdtrk[n_fdtrk]/D");
  T->Branch("eemc_energy_prs1_fdtrk",&eemc_energy_prs1_fdtrk, "eemc_energy_prs1_fdtrk[n_fdtrk]/D");
  T->Branch("eemc_energy_prs2_fdtrk",&eemc_energy_prs2_fdtrk, "eemc_energy_prs2_fdtrk[n_fdtrk]/D");
  T->Branch("eemc_energy_tow_fdtrk",&eemc_energy_tow_fdtrk, "eemc_energy_tow_fdtrk[n_fdtrk]/D");
  T->Branch("eemc_energy_post_fdtrk",&eemc_energy_post_fdtrk, "eemc_energy_post_fdtrk[n_fdtrk]/D");
  T->Branch("eemc_adc_prs1_fdtrk",&eemc_adc_prs1_fdtrk, "eemc_adc_prs1_fdtrk[n_fdtrk]/D");
  T->Branch("eemc_adc_prs2_fdtrk",&eemc_adc_prs2_fdtrk, "eemc_adc_prs2_fdtrk[n_fdtrk]/D");
  T->Branch("eemc_adc_tow_fdtrk",&eemc_adc_tow_fdtrk, "eemc_adc_tow_fdtrk[n_fdtrk]/D");
  T->Branch("eemc_adc_post_fdtrk",&eemc_adc_post_fdtrk, "eemc_adc_post_fdtrk[n_fdtrk]/D");
  T->Branch("eemc_num_towers_in_track", &eemc_num_towers_in_track, "eemc_num_towers_in_track[n_fdtrk]/I");

  // track TOF info
  T->Branch("tofhit_fdtrk",&tofhit_fdtrk,"tofhit_fdtrk[n_fdtrk]/I");
  T->Branch("tofpatlen_fdtrk",&tofpatlen_fdtrk,"tofpatlen_fdtrk[n_fdtrk]/D");
  T->Branch("toftof_fdtrk",&toftof_fdtrk,"toftof_fdtrk[n_fdtrk]/D");
  T->Branch("tofbeta_fdtrk",&tofbeta_fdtrk,"tofbeta_fdtrk[n_fdtrk]/D");
  T->Branch("toflet_fdtrk",&toflet_fdtrk,"toflet_fdtrk[n_fdtrk]/D");
  T->Branch("toftet_fdtrk",&toftet_fdtrk,"toftet_fdtrk[n_fdtrk]/D");

  // max. array for template
  n_fdtrk = max_fdtrk;

}

void uDstSkimMaker::fdtrk_br_fill()
{

  // defaults:
  n_fdtrk = 0;

  //cout << "\n" << endl;

  StMuEvent* muEvent = mMuDstMaker->muDst()->event();
  if (!muEvent) {
    cout << "WARNING: fdtrk_br_fill no muEvent" << endl;
    return;
  }

  // assume MC for Run < 2010
  int isMC = runN < 11000000;

  // loop over vertices,tracks, look for fast detector matches

  Int_t Nvert = mMuDstMaker->muDst()->numberOfPrimaryVertices();
  //cout << "Number of vertices = " << Nvert << endl;
  for (int ivert = 0; ivert<Nvert; ivert++) { // start loop over vertices

    StMuDst::setVertexIndex(ivert);

    // MC matching: require Zvtx close to generated
    const double MCZvtxcut = 1.; // cut 1 cm
    int MC_vtxmat = 0;
	double Zvtxrec = mMuDstMaker->muDst()->event()->primaryVertexPosition().z();
    if (isMC) {
       // this Zvtx
      for (int iMCvtx=0; iMCvtx<n_MCvtx; iMCvtx++) { // loop over gen. vertices
	double Zvtxgen = z_MCvtx[iMCvtx];
	if (abs(Zvtxrec-Zvtxgen) < MCZvtxcut) {MC_vtxmat = 1;}
      }
    }
    if (isMC &&!MC_vtxmat) {continue;} // MC: skip vertex if not matched

    TObjArray* tracks = mMuDstMaker->muDst()->primaryTracks() ; // Create a TObject array containing the primary tracks
    if (!tracks) {cout << "WARNING: fdtrk_br_fill missing tracks" << endl;}
    else { // start tracks OK

      for (int itrk=0; itrk<= tracks->GetLast(); itrk++) { // start loop over tracks

		//track = (StMuTrack*)tracks[itrk]; // error: invalid cast from type 'TObjArray' to type 'StMuTrack*'
		StMuTrack* track = (StMuTrack*)tracks->At(itrk);
		if (!track) {cout << "WARNING: fdtrk_br_fill missing track" << endl;}
		else { // start track OK

	  		// get B-field (not done once, can change run to run)
	  		Double_t bFld = muEvent->eventSummary().magneticField()/10.; // CHECK x10???

	  		// start track->BEMC cluster matching:
	  		Int_t BEMC_ext = 0;
	  		Int_t BEMC_match = 0;
	  		Int_t MTD_match = 0;
	  
	  		StThreeVectorD BEMCposition, BEMCmomentum;
	  		Int_t BEMCmod,BEMCetabin,BEMCsub;
	  		Bool_t ok = project->trackOnEmc(&BEMCposition,&BEMCmomentum,track,bFld);
	  		if (!ok) {/*cout <<" WARNING: fdtrk_br_fill projection failed" << endl;*/}
	  		else { // projection OK
	    		BEMC_ext = !(mGeom->getBin(BEMCposition.phi(),BEMCposition.pseudoRapidity(),BEMCmod,BEMCetabin,BEMCsub));
	    		//cout << "Barrel track: " << " phi = " << BEMCposition.phi() << "  eta = " << BEMCposition.pseudoRapidity() << endl;
	    		if (BEMC_ext && emccl_cells[BEMCmod-1][BEMCsub-1][BEMCetabin-1] >= 0) {BEMC_match = 1;}
	  		}
	  
			if(track->index2MtdHit() >=0){MTD_match=1;}

	  		//cout << "Begin EEMC to track matching..." << endl;
	  		// start track->EEMC cluster matching:
	  		Int_t EEMC_ext = 0;
	  		Int_t EEMC_match = 0;
	  		Int_t distToEEMCFromOriginZ = 279.542; //270.0; //2.7 meters
	  		double clusterEnergy = 0.0;
	  		//double prs1Energy = 0.0;
	  		//double prs2Energy = 0.0;
	  		//double postEnergy = 0.0;
			double maxTowEnergy = 0.0;
			Int_t nTowInCluster = 0;	
 
			//if(track->eta() < 0.95) { continue; }
 
	  		StThreeVectorD EEMCposition, EEMCmomentum;
	  		Int_t EEMCsec = 0, EEMCsub = 0, EEMCetabin = 0;
			Int_t EEMCsec_max = 0, EEMCsub_max = 0, EEMCetabin_max = 0;
	  		Bool_t eemc_ok = trackOnEndcap(&EEMCposition, &EEMCmomentum, track, bFld, distToEEMCFromOriginZ, EEMCsec, EEMCsub, EEMCetabin);
	  		//eemc_ok = eemc_cluster_adc[EEMCsec][EEMCsub][EEMCetabin] > 10;
	  		//if (!eemc_ok) {cout <<" WARNING: eemc projection failed" << endl;}
	  		if (eemc_ok){ // projection OK
		  		EEMC_ext = eemc_ok;
		  		//cout << "EEMC track projected to tower - sector = " << EEMCsec << " subsector = " << EEMCsub << "  etabin = " << EEMCetabin << endl;
	      
		  		//fine adjacent towers bins
		  		//NB: This is a very in-effecient cluster-finder

		  		//int nTowInCluster = 0;

		  		int etaStart = 0; int etaEnd = 0;

		  		int subSecZero[3] = {4, 0, 1};
		  		int subSecFour[3] = {3, 4, 0};
		  		int secNormalZero[3] = {-1, 0, 0};
		  		int secNormalFour[3] = {0, 0, 1};
		  		int secZeroZero[3] = {11, 0, 0};
		  		int secFourFour[3] = {0, 0, -11};

          		//clusters not near a sector boundary
		  		if(EEMCsub > 0 && EEMCsub < 4){

			  		if(EEMCetabin == 0){ etaStart = 0; etaEnd = 2;}
              		if(EEMCetabin == 11){ etaStart = -1; etaEnd = 1;}
			  		if(EEMCetabin > 0 && EEMCetabin < 11){ etaStart = -1; etaEnd = 2;}	    

		      		for(int subRange = -1; subRange < 2; subRange++){
                  		for(int etaRange = etaStart; etaRange < etaEnd; etaRange++){
					  		//cout << "checking - " << EEMCsec << " subsector = " << EEMCsub+subRange << "  etabin = " << EEMCetabin+etaRange << endl;
				      		if (eemc_ok && eemc_tower_energy[EEMCsec][EEMCsub+subRange][EEMCetabin+etaRange] > 0.0) {
                          		//EEMC_match = 1;
                          		//EEMCsec = EEMCsec;
                          		//EEMCsub = EEMCsub+subRange;
                          		//EEMCetabin = EEMCetabin+etaRange;
						  		//EEMC_match = 1;
						  		nTowInCluster++;	
								if(eemc_tower_energy[EEMCsec][EEMCsub+subRange][EEMCetabin+etaRange] > maxTowEnergy){
									
									maxTowEnergy   = eemc_tower_energy[EEMCsec][EEMCsub+subRange][EEMCetabin+etaRange];
									EEMCsec_max    = EEMCsec;
									EEMCsub_max    = EEMCsub+subRange;
									EEMCetabin_max = EEMCetabin+etaRange;

								}
								
								clusterEnergy += eemc_tower_energy[EEMCsec][EEMCsub+subRange][EEMCetabin+etaRange];
						  		//prs1Energy    = eemc_prs1_energy[EEMCsec][EEMCsub+subRange][EEMCetabin+etaRange];
      					  		//prs2Energy 	  = eemc_prs2_energy[EEMCsec][EEMCsub+subRange][EEMCetabin+etaRange];
      					  		//postEnergy    = eemc_post_energy[EEMCsec][EEMCsub+subRange][EEMCetabin+etaRange];
 						  		
								//cout << "Track >>MATCHED<< to EEMC tower - sector = " << EEMCsec << " subsector = " << EEMCsub+subRange << "  etabin = " 
						  	   	//	 << EEMCetabin+etaRange << " energy = " << eemc_tow_energy[EEMCsec][EEMCsub+subRange][EEMCetabin+etaRange] << endl;
								//cout << "Preshower1 Energy = " << eemc_prs1_energy[EEMCsec][EEMCsub+subRange][EEMCetabin+etaRange] << " GeV " << endl;
								//cout << "Preshower2 Energy = " << eemc_prs2_energy[EEMCsec][EEMCsub+subRange][EEMCetabin+etaRange] << " GeV " << endl;
								//cout << "Postshower Energy = " << eemc_post_energy[EEMCsec][EEMCsub+subRange][EEMCetabin+etaRange] << " GeV " << endl;		      
					  		}
				 		}
			  		}
		 		}			
		
				//if(nTowInCluster > 0 && EEMC_match != 1) { 
                //    cout << "-----Cluster Energy = " << clusterEnergy << " GeV with  "<< nTowInCluster << " towers in cluster." <<  endl;
                //    cout << "Max tower energy = " << maxTowEnergy << endl;
                //}
		  		if(nTowInCluster > 0 && clusterEnergy > 0.0 && EEMC_match != 1){ EEMC_match = 1; }          

          		//cluster near 0 edge of sector
          		int subOffset = 0; int secOffset = 0;
	      		if(!EEMC_match && EEMCsub == 0){

              		if(EEMCetabin == 0){ etaStart = 0; etaEnd = 2;}
              		if(EEMCetabin == 11){ etaStart = -1; etaEnd = 1;}
              		if(EEMCetabin > 0 && EEMCetabin < 11){ etaStart = -1; etaEnd = 2;}

			  		for(int i = 0; i < 3; i++){
                  		subOffset = subSecZero[i];
			      		if(EEMCsec > 0) {secOffset = secNormalZero[i];}
                  		if(EEMCsec == 0){secOffset = secZeroZero[i];}
				  		for(int etaRange = etaStart; etaRange < etaEnd; etaRange++){
					  		//cout << "checking - " << EEMCsec+secOffset << " subsector = " << subOffset << "  etabin = " << EEMCetabin+etaRange << endl;
                      		if (eemc_ok && eemc_tower_energy[EEMCsec+secOffset][subOffset][EEMCetabin+etaRange] > 0.0) {
                          		//EEMC_match = 1;
                          		//EEMCsec = EEMCsec+secOffset;
                          		//EEMCsub = subOffset;
                          		//EEMCetabin = EEMCetabin+etaRange;
						  		//EEMC_match = 1;
						  		nTowInCluster++;
								if(eemc_tower_energy[EEMCsec+secOffset][subOffset][EEMCetabin+etaRange] > maxTowEnergy){
                                    
                                    maxTowEnergy   = eemc_tower_energy[EEMCsec+secOffset][subOffset][EEMCetabin+etaRange];
                                    EEMCsec_max    = EEMCsec+secOffset;
                                    EEMCsub_max    = subOffset;
                                    EEMCetabin_max = EEMCetabin+etaRange;

                                }
                          		clusterEnergy += eemc_tower_energy[EEMCsec+secOffset][subOffset][EEMCetabin+etaRange];
						  		//prs1Energy    = eemc_prs1_energy[EEMCsec+secOffset][subOffset][EEMCetabin+etaRange];
                          		//prs2Energy    = eemc_prs2_energy[EEMCsec+secOffset][subOffset][EEMCetabin+etaRange];
                          		//postEnergy    = eemc_post_energy[EEMCsec+secOffset][subOffset][EEMCetabin+etaRange];
                          		

								//cout << "Track >>MATCHED<< to EEMC tower - sector = " << EEMCsec+secOffset << " subsector = " << subOffset << "  etabin = "
                               	//	<< EEMCetabin+etaRange << " energy = " <<  eemc_tow_energy[EEMCsec+secOffset][subOffset][EEMCetabin+etaRange] << endl; 
								//cout << "Preshower1 Energy = " << 			  eemc_prs1_energy[EEMCsec+secOffset][subOffset][EEMCetabin+etaRange] << " GeV " << endl;    
                                //cout << "Preshower2 Energy = " << 			  eemc_prs2_energy[EEMCsec+secOffset][subOffset][EEMCetabin+etaRange] << " GeV " << endl;
                                //cout << "Postshower Energy = " << 			  eemc_post_energy[EEMCsec+secOffset][subOffset][EEMCetabin+etaRange] << " GeV " << endl;

                      		}
                  		}
              		}
		  		}
			
				//if(nTowInCluster > 0 && EEMC_match != 1) { 
                //    cout << "-----Cluster Energy = " << clusterEnergy << " GeV with  "<< nTowInCluster << " towers in cluster." <<  endl;
                //    cout << "Max tower energy = " << maxTowEnergy << endl;
                //}
          		if(nTowInCluster > 0 && clusterEnergy > 0.0 && EEMC_match != 1){ EEMC_match = 1; }

	 	  		//cluster near 4 edge of sector
		  		if(!EEMC_match && EEMCsub == 4){

              		if(EEMCetabin == 0){ etaStart = 0; etaEnd = 2;}
              		if(EEMCetabin == 11){ etaStart = -1; etaEnd = 1;}
              		if(EEMCetabin > 0 && EEMCetabin < 11){ etaStart = -1; etaEnd = 2;}

              		for(int i = 0; i < 3; i++){
                  		subOffset = subSecFour[i];
                  		if(EEMCsec < 11) {secOffset = secNormalFour[i];}
                  		if(EEMCsec == 11){secOffset = secFourFour[i];}
                  		for(int etaRange = etaStart; etaRange < etaEnd; etaRange++){
                  	  		//cout << "checking - " << EEMCsec+secOffset << " subsector = " << subOffset << "  etabin = " << EEMCetabin+etaRange << endl;
					  		if (eemc_ok && eemc_tower_energy[EEMCsec+secOffset][subOffset][EEMCetabin+etaRange] > 0.0) {
                          		//EEMC_match = 1;
                          		//EEMCsec = EEMCsec+secOffset;
                          		//EEMCsub = subOffset;
                          		//EEMCetabin = EEMCetabin+etaRange;
                          		//EEMC_match = 1;
                          		nTowInCluster++;
								if(eemc_tower_energy[EEMCsec+secOffset][subOffset][EEMCetabin+etaRange] > maxTowEnergy){

                                    maxTowEnergy   = eemc_tower_energy[EEMCsec+secOffset][subOffset][EEMCetabin+etaRange];
                                    EEMCsec_max    = EEMCsec+secOffset;
                                    EEMCsub_max    = subOffset;
                                    EEMCetabin_max = EEMCetabin+etaRange;

                                }
                          		clusterEnergy  += eemc_tower_energy[EEMCsec+secOffset][subOffset][EEMCetabin+etaRange];
						  		//prs1Energy     = eemc_prs1_energy[EEMCsec+secOffset][subOffset][EEMCetabin+etaRange];
                          		//prs2Energy     = eemc_prs2_energy[EEMCsec+secOffset][subOffset][EEMCetabin+etaRange];
                          		//postEnergy     = eemc_post_energy[EEMCsec+secOffset][subOffset][EEMCetabin+etaRange];
                          		
								//cout << "Track >>MATCHED<< to EEMC tower - sector = " << EEMCsec+secOffset << " subsector = " << subOffset << "  etabin = "
                               	//	<< EEMCetabin+etaRange << " energy = " << eemc_tow_energy[EEMCsec+secOffset][subOffset][EEMCetabin+etaRange] << endl;
								//cout << "Preshower1 Energy = " 			   << eemc_prs1_energy[EEMCsec+secOffset][subOffset][EEMCetabin+etaRange] << " GeV " << endl;    
                                //cout << "Preshower2 Energy = " 			   << eemc_prs2_energy[EEMCsec+secOffset][subOffset][EEMCetabin+etaRange] << " GeV " << endl;
                                //cout << "Postshower Energy = " 			   << eemc_post_energy[EEMCsec+secOffset][subOffset][EEMCetabin+etaRange] << " GeV " << endl;
                      		}
                  		}
              		}
		  		}
	
		  		//if(nTowInCluster > 0 && EEMC_match != 1) { 
				//	cout << "-----Cluster Energy = " << clusterEnergy << " GeV with  "<< nTowInCluster << " towers in cluster." <<  endl;
				//	cout << "Max tower energy = " << maxTowEnergy << endl;
				//}
		  		if(nTowInCluster > 0 && clusterEnergy > 0.0 && EEMC_match != 1 ){ EEMC_match = 1; }
			
		  		//if(EEMC_match){

				//cout << "Final matched track eta = " << EEMCposition.pseudoRapidity() << endl;
            	//cout << "Final matched track phi = " << EEMCposition.phi() << endl;

		  		//}

	  	}//end EEMC projection OKAY

	  // start track->TOF hit matching:
	  Int_t TOF_match = 0;
	  const StMuBTofPidTraits btof = track->btofPidTraits();
	  if (btof.matchFlag() != 0) { // track->TOF hit match
	    TOF_match = 1;
	  }

	  // start track->MC gen. matching: track 3-vector close to generated and same sign
	  //const double MCtrkcut = 5.; // no cut
	  //const double MCtrkcut = 0.5; // cut 0.5 rad
	  //int MC_trkmat = 0;
	  if (isMC) {
	    //int Qrec = track->charge();
	    TVector3 P3rec; P3rec.SetPtEtaPhi(track->pt(),track->eta(),track->phi());
	    for (int iMCtrk=0; iMCtrk<n_MCtrk; iMCtrk++) { // loop over gen. tracks
	      TVector3 P3gen; P3gen.SetXYZ(px_MCtrk[iMCtrk],py_MCtrk[iMCtrk],pz_MCtrk[iMCtrk]);
	      //if (Qrec==q_MCtrk[iMCtrk] && P3rec.Angle(P3gen)<MCtrkcut) {MC_trkmat = 1;}
	    }
	  }

	  // define fast detector match:
	  //Int_t FD_match = BEMC_match;
	  //Int_t FD_match = TOF_match;
	  //Int_t FD_match = BEMC_match && TOF_match;
	  //Int_t FD_match = BEMC_match || TOF_match;
	  //Int_t FD_match = BEMC_match;
	  //Int_t FD_match = EEMC_match;
	  Int_t FD_match = BEMC_match || EEMC_match || TOF_match || MTD_match;
	  //Int_t FD_match = BEMC_match || EEMC_match;

	  // FD or MC match
	  //if (FD_match || MC_trkmat) { // start fast detector or MC match
	  if (FD_match) { // start fast detector or MC match

	    // ifdvtx_fdtrk SET IN VTX CODE
	    ivtx_fdtrk[n_fdtrk] = ivert; // only for local use
	    idx_fdtrk[n_fdtrk] = itrk; // only for local use

	    // track quantities
	    q_fdtrk[n_fdtrk] = track->charge();
	    pt_fdtrk[n_fdtrk] = track->pt();
	    pz_fdtrk[n_fdtrk] = track->p().z();
	    p_fdtrk[n_fdtrk] = track->p().mag();
	    eta_fdtrk[n_fdtrk] = track->eta();
	    phi_fdtrk[n_fdtrk] = track->phi();
	    nh_fdtrk[n_fdtrk] = track->nHits();
	    dcar_fdtrk[n_fdtrk] = track->dcaGlobal().perp();
	    dcaz_fdtrk[n_fdtrk] = track->dcaGlobal().z();
		vtxz_fdtrk[n_fdtrk] = Zvtxrec;
		vtxi_fdtrk[n_fdtrk] = ivert;
		
	    // track dE/dx info
	    dedx_fdtrk[n_fdtrk] = track->dEdx();
	    nhdedx_fdtrk[n_fdtrk] = track->nHitsDedx();
	    sigel_fdtrk[n_fdtrk] = track->nSigmaElectron();
	    sigpi_fdtrk[n_fdtrk] = track->nSigmaPion();
	    sigk_fdtrk[n_fdtrk] = track->nSigmaKaon();
	    sigp_fdtrk[n_fdtrk] = track->nSigmaProton();

	    // track->BEMC info
	    //emcext_fdtrk[n_fdtrk] = BEMC_ext;
	    //etaext_fdtrk[n_fdtrk] = BEMCposition.pseudoRapidity();
	    //phiext_fdtrk[n_fdtrk] = BEMCposition.phi();
	    if (BEMC_match) {
	        emcext_fdtrk[n_fdtrk] = BEMC_ext;
        	etaext_fdtrk[n_fdtrk] = BEMCposition.pseudoRapidity();
       		phiext_fdtrk[n_fdtrk] = BEMCposition.phi();
			iemccl_fdtrk[n_fdtrk] = emccl_cells[BEMCmod-1][BEMCsub-1][BEMCetabin-1];
	    }
	    else {
	      iemccl_fdtrk[n_fdtrk] = -9999;
	    }
		mtd_match_fdtrk[n_fdtrk] = MTD_match;	
	    // track->EEMC info
	    //eemc_ext_fdtrk[n_fdtrk] = EEMC_ext;
	    //eemc_etaext_fdtrk[n_fdtrk] = EEMCposition.pseudoRapidity();
	    //eemc_phiext_fdtrk[n_fdtrk] = EEMCposition.phi();
	    if (EEMC_match) {

			//maxTowEnergy
            //EEMCsec_max   
            //EEMCsub_max
            //EEMCetabin_max

	      	eemc_ext_fdtrk[n_fdtrk]    = EEMC_ext;
        	eemc_etaext_fdtrk[n_fdtrk] = EEMCposition.pseudoRapidity();
        	eemc_phiext_fdtrk[n_fdtrk] = EEMCposition.phi();
			ieemc_cl_fdtrk[n_fdtrk]    = eemc_tower_cells[EEMCsec_max][EEMCsub_max][EEMCetabin_max];
			eemc_adc_fdtrk[n_fdtrk]    = eemc_tower_adc[EEMCsec_max][EEMCsub_max][EEMCetabin_max];
			eemc_energy_fdtrk[n_fdtrk] = clusterEnergy; //emc_cluster_energy[EEMCsec][EEMCsub][EEMCetabin];
			
			//high tower energies for core tower and pre/post-shower layers!
			eemc_energy_prs1_fdtrk[n_fdtrk] = eemc_prs1_energy[EEMCsec_max][EEMCsub_max][EEMCetabin_max];
            eemc_energy_prs2_fdtrk[n_fdtrk] = eemc_prs2_energy[EEMCsec_max][EEMCsub_max][EEMCetabin_max];
			eemc_energy_tow_fdtrk[n_fdtrk]  = maxTowEnergy; //eemc_tow_energy[EEMCsec_max][EEMCsub_max][EEMCetabin_max];
            eemc_energy_post_fdtrk[n_fdtrk] = eemc_post_energy[EEMCsec_max][EEMCsub_max][EEMCetabin_max];

			//adc information 
			eemc_adc_prs1_fdtrk[n_fdtrk] = eemc_prs1_adc[EEMCsec_max][EEMCsub_max][EEMCetabin_max]; //raw ADC
            eemc_adc_prs2_fdtrk[n_fdtrk] = eemc_prs2_adc[EEMCsec_max][EEMCsub_max][EEMCetabin_max]; //raw ADC
            eemc_adc_tow_fdtrk[n_fdtrk]  = eemc_tower_adc[EEMCsec_max][EEMCsub_max][EEMCetabin_max]; //raw ADC
            eemc_adc_post_fdtrk[n_fdtrk] = eemc_post_adc[EEMCsec_max][EEMCsub_max][EEMCetabin_max]; //raw ADC
		
			//some other important information to help with cuts
			eemc_num_towers_in_track[n_fdtrk] = nTowInCluster;	

			//cout << "matched EEMC track saved" << endl;
	    }
	    else {
	      ieemc_cl_fdtrk[n_fdtrk] = -9999;
	    }
		

	    // track TOF info
	    if (TOF_match) {
	      tofhit_fdtrk[n_fdtrk] = 1;
	      tofpatlen_fdtrk[n_fdtrk] = btof.pathLength();
	      toftof_fdtrk[n_fdtrk] = btof.timeOfFlight();
	      tofbeta_fdtrk[n_fdtrk] = btof.beta();
	      const StMuBTofHit * btofHit = track->tofHit();
	      //const StMuBTofHit btofHit = track->tofHit(); DON'T WORK
	      if (btofHit) {
		toflet_fdtrk[n_fdtrk] = btofHit->leadingEdgeTime();
		toftet_fdtrk[n_fdtrk] = btofHit->trailingEdgeTime();
	      }
	      else {
		toflet_fdtrk[n_fdtrk] = -99;
		toftet_fdtrk[n_fdtrk] = -99;
	      }
	    }
	    else {
	      tofhit_fdtrk[n_fdtrk] = 0;
	      tofpatlen_fdtrk[n_fdtrk] = -9999;
	      toflet_fdtrk[n_fdtrk] = -9999;
	      toftet_fdtrk[n_fdtrk] = -9999;
	    }
	    //cout << "BTOF: " << tofhit_fdtrk[n_fdtrk]
	    //   <<" "<< btof.position().x() << endl; ALWAYS 0 IF NO TOF HIT

	    n_fdtrk++;
	    if (n_fdtrk >= max_fdtrk) {
	      cout <<" WARNING: fdtrk_br_fill nfdtrk limit reached" << endl;
	      return;
	    }

	  } // end track->fast detector match
	} // end track OK
      } // end loop over tracks
    } // end tracks OK
  } // end loop over vertices
	
	//cout << "number of fd_tracks = " << n_fdtrk << endl;

}

Bool_t projectEEMCTrack(StThreeVectorD* atFinal, StThreeVectorD* momentumAtFinal, 
								const StMuTrack* const track, double magField, double radius, int option)
{
    StThreeVectorD Zero(0,0,0);
    *atFinal=Zero;
    *momentumAtFinal=Zero;
	
    /* this was for StTrack
		const StThreeVectorF& origin = track->geometry()->origin();
	const StThreeVectorF& momentum = track->geometry()->momentum();
	double charge = track->geometry()->charge();
	StPhysicalHelixD helix(momentum, origin, magField*tesla, charge);
    */
    StPhysicalHelixD helix = track->outerHelix();
    const StThreeVectorF momentum = track->momentum();
	
	const StThreeVectorF rVect(0,0,radius);
	const StThreeVectorF nVect(0,0,1);
	
    double pathLength = helix.pathLength(rVect, nVect);
    double charge = track->charge();
	
    double s;//,s1,s2; 
    //s=0;
    s = pathLength;//.first;
    //s2 = pathLength.second;
	
	//cout << "Helix path length = " << s << endl;
	
    Bool_t goProj;
    goProj = kFALSE;
	
    if (finite(s) <= 0 ) { return kFALSE;} // Track couldn't be projected!
	else goProj = kTRUE;
	
    /*if (option == 1)  // Selects positive path lenght to project track forwards along its helix relative to
				  // first point of track. The smaller solution is taken when both are positive
	{
	    if (s1 >= 0 && s2 >= 0) {s = s1; goProj = kTRUE; }
	    if (s1 >= 0 && s2 < 0) { s = s1; goProj = kTRUE; }
	    if (s1 < 0 && s2 >= 0) { s = s2; goProj = kTRUE; }
	}
	
    if (option == -1) // Selects negative path lenght to project track backwards along its helix relative to
					  // first point of track. The smaller absolute solution is taken when both are negative 
	{
	    if (s1 <= 0 && s2 <= 0) { s = s2; goProj = kTRUE; }
	    if (s1 <= 0 && s2 > 0) { s = s1; goProj = kTRUE; }
	    if (s1 > 0 && s2 <= 0) { s = s2; goProj = kTRUE; }
	}
	*/
    if (goProj) 
	{
	    *atFinal = helix.at( s );
	    *momentumAtFinal = helix.momentumAt( s, magField*tesla );
		//*momentumAtFinal = helix.momentum(magField*tesla);
	    if (charge == 0) *momentumAtFinal = momentum;
	}
    return goProj;
}

Bool_t trackOnEndcap(StThreeVectorD* position, StThreeVectorD* momentum, const StMuTrack* const track, double magField, double endcapDist, int &sector, int &sub, int &etabin ){
	
	StPhysicalHelixD helix = track->outerHelix();
	const StThreeVectorD& origin = helix.origin();
	//StEmcPosition * proj;
	
	//cout << "begin track projection to EEMC..." << endl;
	
	StThreeVectorF tpcMom = track->momentum();

	//cout << " Eta 'track' = " << track->eta() << "  Eta 'momentum' = " << tpcMom.pseudoRapidity() << endl;

	if(track->eta() < 0.7){ return kFALSE; }

    float xO = origin.x();
	float yO = origin.y();
	float distToOrigin = sqrt( pow(xO, 2) + pow(yO, 2) );    
	//cout << "dist to origin = " << distToOrigin << endl;
	if ( distToOrigin < endcapDist ){
		
		Bool_t projTrackOk = projectEEMCTrack( position, momentum, track, magField, endcapDist, 1) && position->perp() > 0.0;
		//cout << "track projected?  " << projTrackOk << " transerse component = " << position->perp() <<  endl;
		if ( projTrackOk && position->pseudoRapidity() > 0.8 && position->pseudoRapidity() < 2.0 ){//0.7 3.2
			int sectortmp = 0, subtmp = 0, etabintmp = 0;
			Float_t dphi = 0, deta = 0;
			
			float phi = position->phi();
			float eta = position->pseudoRapidity();
			
			//if(eta < 1.0){ cout << "invalid eta value" << endl; return kFALSE;}
			
			Double_t z   = endcapDist;
			Double_t rho = z/sinh(eta); 
			//cout << "z = " << z << endl;	
			//cout << "rho = " << rho << endl;	
			//cout << "eta = " << eta << endl;
			//cout << "extracted projected vector parameters...   phi = " << phi << ", eta = " << eta <<  endl;
			
			TVector3 emcDir(rho*cos(phi), rho*sin(phi), z);
			
			//cout << "emcDir.Pt() = " << emcDir.Pt() << endl;	
			//cout << "emc Vector - x = " << emcDir.X() << "  y = " << emcDir.y() << "  z = " << emcDir.Z() << "  eta = " << eta << endl;	
			
			if(tpcMom.pseudoRapidity() < 1.0){ return kFALSE; }
	
			if ( endcapGeom->getTower(emcDir, sectortmp, subtmp, etabintmp, dphi, deta) ){ 
			
				sector = sectortmp;
				sub = subtmp;
				etabin = etabintmp;
			
				//cout << "Tower found!!!  sector = " << sector << " sub = " << sub << " etabin = " << etabin << endl;
				//cout << "deltaPhi = " << dphi <<  "  deltaEta = " << deta << endl;
				//cout << "TPC track eta =  " << tpcMom.pseudoRapidity() << endl;				

				return kTRUE;
			}      
		}
	} 
	
	    return kFALSE;
}


#endif
