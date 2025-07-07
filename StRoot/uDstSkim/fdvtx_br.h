// FD vertex branch header

#ifndef fdvtx_br_h
#define fdvtx_br_h

// vertex tree vars:
const Int_t max_fdvtx = 1000;
Int_t n_fdvtx;
Int_t ntrk_fdvtx[max_fdvtx];
Int_t nfdtrk_fdvtx[max_fdvtx];
Double_t x_fdvtx[max_fdvtx];
Double_t y_fdvtx[max_fdvtx];
Double_t z_fdvtx[max_fdvtx];
Double_t rank_fdvtx[max_fdvtx];

// vertex branch code

void uDstSkimMaker::fdvtx_br_init()
{

  // tree vertex vars:
  T->Branch("n_fdvtx",&n_fdvtx,"n_fdvtx/I");
  T->Branch("ntrk_fdvtx",&ntrk_fdvtx,"ntrk_fdvtx[n_fdvtx]/I");
  T->Branch("nfdtrk_fdvtx",&nfdtrk_fdvtx,"nfdtrk_fdvtx[n_fdvtx]/I");
  T->Branch("x_fdvtx",&x_fdvtx,"x_fdvtx[n_fdvtx]/D");
  T->Branch("y_fdvtx",&y_fdvtx,"y_fdvtx[n_fdvtx]/D");
  T->Branch("z_fdvtx",&z_fdvtx,"z_fdvtx[n_fdvtx]/D");
  T->Branch("rank_fdvtx",&rank_fdvtx,"rank_fdvtx[n_fdvtx]/D");

  // max. array for template
  n_fdvtx = max_fdvtx;

}

void uDstSkimMaker::fdvtx_br_fill()
{

  // defaults:
  n_fdvtx = 0;

  Int_t nvtx = mMuDstMaker->muDst()->numberOfPrimaryVertices();
  for (int ivtx=0; ivtx<nvtx; ivtx++) { // start loop over all vertices

    // count FD tracks this vertex, and reset fdtrk vertex pointer to fdvtx list
    nfdtrk_fdvtx[n_fdvtx] = 0;
    for (int itrk=0; itrk<n_fdtrk; itrk++) {
      if (ivtx_fdtrk[itrk] == ivtx) {
	      ifdvtx_fdtrk[itrk] = n_fdvtx; // will be index in fdvtx list
	      nfdtrk_fdvtx[n_fdvtx]++;
      }
    }

    // keep vertices with FD track(s)
    //if (nfdtrk_fdvtx[n_fdvtx]>0) { 
    // keep all vertices
    if (nfdtrk_fdvtx[n_fdvtx]>=0) {
      //StMuPrimaryVertex* V= mMuDstMaker->muDst()->primaryVertex(ivtx);
      //assert(V);

      //StMuDst::setVertexIndex(ivtx);



      // checked: consistent with length prim. tracks array each vtx (fdtrk_br.h):
      ntrk_fdvtx[n_fdvtx] = mMuDstMaker->muDst()->primaryTracks()->GetEntries();
      x_fdvtx[n_fdvtx] = mMuDstMaker->muDst()->event()->primaryVertexPosition().x();
      y_fdvtx[n_fdvtx] = mMuDstMaker->muDst()->event()->primaryVertexPosition().y();
      z_fdvtx[n_fdvtx] = mMuDstMaker->muDst()->event()->primaryVertexPosition().z();

      //rank_fdvtx[n_fdvtx] = mMuDstMaker->muDst()->primaryVertex(ivtx)->ranking();
    	StMuPrimaryVertex* V= mMuDstMaker->muDst()->primaryVertex(ivtx);
	    assert(V);
	    mMuDstMaker->muDst()->setVertexIndex(ivtx);
	    rank_fdvtx[n_fdvtx]=V->ranking();

      n_fdvtx++;
      if (n_fdvtx >= max_fdvtx) {
	      cout <<" WARNING: fdvtx_br_fill n_fdvtx limit reached" << endl;
	      return;
      }

    }

  } // end loop over all vertices

}

#endif
