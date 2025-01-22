// Endcap EMC cluster branch header
//
//
// Author: Alex Jentsch
//
//


#ifndef eemc_clusters_h
#define eemc_clusters_h
  
// table of cells in clusters, and energies
Int_t eemc_tower_cells[12][5][12];
Double_t eemc_prs1_energy[12][5][12];
Double_t eemc_prs2_energy[12][5][12];
Double_t eemc_tower_energy[12][5][12];
Double_t eemc_post_energy[12][5][12];
Double_t eemc_prs1_adc[12][5][12];
Double_t eemc_prs2_adc[12][5][12];
Double_t eemc_tower_adc[12][5][12];
Double_t eemc_post_adc[12][5][12];

// EEMC 'clusters' tree vars
const Int_t max_eemc_towers = 720;
const Double_t towEThreshold = 0.1;
//Int_t n_eemc_clusters;


StEEmcDb        *mDbE; // access to EEMC database 


// EMC cluster branch code

void uDstSkimMaker::eemc_clusters_init()
{

    // eemc cluster vars:
	T->Branch("eemc_tower_cells",eemc_tower_cells,"eemc_tower_cells[12][5][12]/I");
	T->Branch("eemc_prs1_energy",eemc_prs1_energy,"eemc_prs1_energy[12][5][12]/D");
	T->Branch("eemc_prs2_energy",eemc_prs2_energy,"eemc_prs2_energy[12][5][12]/D");
	T->Branch("eemc_tower_energy",eemc_tower_energy,"eemc_tower_energy[12][5][12]/D");
	T->Branch("eemc_post_energy",eemc_post_energy,"eemc_post_energy[12][5][12]/D");	
	T->Branch("eemc_prs1_adc",eemc_prs1_adc,"eemc_prs1_adc[12][5][12]/D");
    T->Branch("eemc_prs2_adc",eemc_prs2_adc,"eemc_prs2_adc[12][5][12]/D");
    T->Branch("eemc_tower_adc",eemc_tower_adc,"eemc_tower_adc[12][5][12]/D");
    T->Branch("eemc_post_adc",eemc_post_adc,"eemc_post_adc[12][5][12]/D");

	
	mDbE = (StEEmcDb*)GetDataSet("StEEmcDb"); 
	assert(mDbE);

}

void uDstSkimMaker::eemc_clusters_fill()
{

	// defaults:

	// clear cells used table
	for (int i=0; i<12; i++) { 
		for(int j=0; j<5; j++) { 
			for(int k=0; k<12; k++) {
			
				eemc_tower_cells[i][j][k] = -9999;
				eemc_prs1_energy[i][j][k] = 0.;
				eemc_prs2_energy[i][j][k] = 0.;
				eemc_tower_energy[i][j][k] = 0.;
				eemc_post_energy[i][j][k] = 0.;
				eemc_prs1_adc[i][j][k] = 0.;
                eemc_prs2_adc[i][j][k] = 0.;
                eemc_tower_adc[i][j][k] = 0.;
                eemc_post_adc[i][j][k] = 0.;
			}
		}
	}

	// for BEMC a la $STAR/StRoot/macros/mudst/exampleEmc.C
	// get BEMC clusters (several layers deep), get cluster info & max. E cell

  	//cout << "...begin EEMC cluster analysis..." << endl;
  	//StEmcCollection *emcCollectionEndcap = mMuDstMaker->muDst()->emcCollection(); // EMC
  
  	//EEMC doesn't seem to have emcCollection - try StMuEmcCollection
  	StMuEmcCollection *muEmcCollection = mMuDstMaker->muDst()->muEmcCollection();
	//EEmcGeomSimple *endcapGeomMap = new EEmcGeomSimple();
	//mDbE = (StEEmcDb*)GetDataSet("StEEmcDb"); //database object to obtain gains

  	if (!muEmcCollection) {cout << "WARNING: eemc_cluster_fill no emc collection" << endl;}
  	else { // start emc collection OK

		//This first loop is for the preshower information
        //keeping them separate is a bit inefficient (two loops), but wiser to keep information cleanly separated

		for(int prsHit = 0; prsHit < muEmcCollection->getNEndcapPrsHits(); prsHit++){

			int prsSec, prsSub, prsEtabin, pre;
			double energy = 0.0;			

			StMuEmcHit* preShowerMuDstHit = muEmcCollection->getEndcapPrsHit(prsHit, prsSec, prsSub, prsEtabin, pre);
			
			double rawadc = preShowerMuDstHit->getAdc();
			if(rawadc < 0){continue;}

			static const Char_t detectors[] = { 'T', 'P', 'Q', 'R' };
			const EEmcDbItem *x=mDbE->getTile(prsSec,'A'+prsSub-1,prsEtabin, detectors[pre]);
            assert(x); // it should never happened for muDst
            if(x->fail ) continue; // drop not working channels

            double adc = rawadc - (x->ped);
			
			if(adc<3.0*x->sigPed) { adc = 0.0; }
			//if(adc < 0){ continue; }


			if(x->gain <=0) {energy = 0.0;} //drop channels with no gain
			else energy = (adc/(x->gain));

			if(energy < 0.0 || adc < 0) { continue;}

			//eemc_tower_adc[prsSec-1][prsSub-1][prsEtabin-1] += rawadc;    // may need to deprecate this
			if(pre == 1){
				eemc_prs1_energy[prsSec-1][prsSub-1][prsEtabin-1] = energy;
				eemc_prs1_adc[prsSec-1][prsSub-1][prsEtabin-1] = adc;
				//cout << "preshower 1 rawADC = " << rawadc << " and pedestal is = " << x->ped << endl;
			}
            if(pre == 2){
				eemc_prs2_energy[prsSec-1][prsSub-1][prsEtabin-1] = energy;
				eemc_prs2_adc[prsSec-1][prsSub-1][prsEtabin-1] = adc;
				//cout << "preshower 2 rawADC = " << rawadc << " and pedestal is = " << x->ped << endl;
			}
            if(pre == 3){
				eemc_post_energy[prsSec-1][prsSub-1][prsEtabin-1] = energy;
				eemc_post_adc[prsSec-1][prsSub-1][prsEtabin-1] = adc;
				//cout << "posthower rawADC = " << rawadc << " and pedestal is = " << x->ped << endl;
			}

		}//end loop voer EEMC pre- and post-shower

		//loop to get normal, total EEMC tower ADC information, and calibrate with gains
		for(int towId = 0; towId < muEmcCollection->getNEndcapTowerADC(); ++towId) {
	  		
			int rawadc, sec, sub, etabin;
			//int prsSec, prsSub, prsEtabin, pre;
			muEmcCollection->getEndcapTowerADC(towId, rawadc, sec, sub, etabin);
            
			double energy = 0.0;
		    const EEmcDbItem *x=mDbE->getTile(sec,'A'+sub-1,etabin,'T');
		    assert(x); // it should never happened for muDst
		    if(x->fail ) continue; // drop not working channels
	
			double adc = rawadc - (x->ped);
	
			if(adc<3.0*x->sigPed) { adc = 0.0; } //reject if adc < 3Sigma above pedestal
				
			if(x->gain <=0) {energy = 0.0;} //drop channels with no gain
			else energy = adc/(x->gain);

			//if(adc < 0 || energy < 0){ adc = 0.0; energy = 0.0;} //set negative (unphysical) values to zero -- hopefully avoid actual errors since E < 0.2 dropped anyhow
	
			eemc_tower_cells[sec-1][sub-1][etabin-1] = towId; //just a tower Id from 0 to 719
			eemc_tower_adc[sec-1][sub-1][etabin-1]  = adc;	//pedestal-subtracted ADC values, stored for debugging purposes	
			eemc_tower_energy[sec-1][sub-1][etabin-1] = energy;  //calibrated tower energy - ACTUAL, USED ENERGY

	
		} // end loop over EEMC towers
	} // end emc collection OK
} // end main function



#endif
////////////////////////
//reference ONLY
////////////////////////

			//adc_eemc_cluster[towId] = rawadc; // deprecated
			//e_eemc_cluster[towId] = energy; //deprecated
	   	   	//eta_eemc_cluster[towId] = endcapGeomMap->getEtaMean(etabin-1); // deprecated
	    	//phi_eemc_cluster[towId] = endcapGeomMap->getPhiMean(sec-1, sub-1); // deprecated
				
	    	//sigeta_eemc_cluster[towId] = 2.0*endcapGeomMap->getEtaHalfWidth(etabin-1); // deprecated
	    	//sigphi_eemc_cluster[towId] = 2.0*endcapGeomMap->getPhiHalfWidth(etabin-1); // deprecated

/*		
for (int i=0; i< emc->getNEndcapTowerADC(); i++) {
	int sec,eta,sub,rawAdc; //muDst  ranges:sec:1-12, sub:1-5, eta:1-12
	emc->getEndcapTowerADC(i,rawAdc,sec,sub,eta);

	const EEmcDbItem *x=mDbE->getTile(sec,'A'+sub-1,eta,'T');
	assert(x); // it should never happened for muDst
	if(x->fail ) continue; // drop not working channels
	int isec=x->sec-1;
	int isub=x->sub-'A';
	int ieta=x->eta-1;

	assert(isec>=0 && isec<mxEtowSec); // check input is ok
	assert(isub>=0 && isub<mxEtowSub);
	assert(ieta>=0 && ieta<mxEtowEta);

	float adc=rawAdc-x->ped; // ped subtracted ADC
	if(adc<par_kSigPed*x->sigPed) continue;

	wEve->etow.adc[isec*mxEtowSub+isub][ieta]=adc;

	if(x->gain<=0) continue;// drop channels w/o gains
	float ene=adc/x->gain;

	//method for shifting energy scale 
	ene*=par_etowScale;//(default is par_etowScale=1)
    wEve->etow.ene[isec*mxEtowSub+isub][ieta]=ene;
    wEve->etow.stat[isec*mxEtowSub+isub][ieta]=0;

    if(maxADC<adc) { maxIdName=x->name; maxADC=adc; maxSec=isec; maxSub=isub; maxEta=ieta;}
    adcSum+=adc;
}

wEve->etow.maxAdc=maxADC;
wEve->etow.maxSec=maxSec; wEve->etow.maxSub=maxSub; wEve->etow.maxEta=maxEta;
hE[31]->Fill(maxADC);
hE[32]->Fill(adcSum);

if(maxADC<par_maxADC)  return -2 ;  // not enough energy

return 0;
*/
