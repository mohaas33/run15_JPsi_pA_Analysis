#define sDstAnaMaker_cxx
#include "sDstAnaMaker.h"

#include "TFile.h"
#include "TH1F.h"
#include <TH2F.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TObjArray.h"


Int_t sDstAnaMaker::Init( )
{
   std::cout<< "sDstAnaMaker::Init" <<std::endl;
   //InitTree(tree);
   // Variable binning
   Double_t m_bins[] = {
    1.00000000e-01, 1.03608624e-01, 1.07347470e-01, 1.11221236e-01,
    1.15234793e-01, 1.19393183e-01, 1.23701634e-01, 1.28165561e-01,
    1.32790575e-01, 1.37582487e-01, 1.42547322e-01, 1.47691319e-01,
    1.53020943e-01, 1.58542894e-01, 1.64264111e-01, 1.70191785e-01,
    1.76333367e-01, 1.82696575e-01, 1.89289408e-01, 1.96120151e-01,
    2.03197390e-01, 2.10530020e-01, 2.18127257e-01, 2.25998649e-01,
    2.34154091e-01, 2.42603832e-01, 2.51358492e-01, 2.60429075e-01,
    2.69826981e-01, 2.79564022e-01, 2.89652437e-01, 3.00104905e-01,
    3.10934562e-01, 3.22155022e-01, 3.33780385e-01, 3.45825265e-01,
    3.58304798e-01, 3.71234672e-01, 3.84631135e-01, 3.98511027e-01,
    4.12891792e-01, 4.27791504e-01, 4.43228891e-01, 4.59223356e-01,
    4.75795000e-01, 4.92964653e-01, 5.10753894e-01, 5.29185082e-01,
    5.48281382e-01, 5.68066796e-01, 5.88566191e-01, 6.09805332e-01,
    6.31810914e-01, 6.54610595e-01, 6.78233030e-01, 7.02707911e-01,
    7.28065998e-01, 7.54339162e-01, 7.81560427e-01, 8.09764004e-01,
    8.38985343e-01, 8.69261170e-01, 9.00629538e-01, 9.33129872e-01,
    9.66803021e-01, 1.00169131e+00, 1.03783858e+00, 1.07529027e+00,
    1.11409346e+00, 1.15429690e+00, 1.19595114e+00, 1.23910852e+00,
    1.28382329e+00, 1.33015164e+00, 1.37815181e+00, 1.42788413e+00,
    1.47941110e+00, 1.53279749e+00, 1.58811039e+00, 1.64541932e+00,
    1.70479632e+00, 1.76631601e+00, 1.83005571e+00, 1.89609554e+00,
    1.96451850e+00, 2.03541059e+00, 2.10886091e+00, 2.18496177e+00,
    2.26380882e+00, 2.34550117e+00, 2.43014149e+00, 2.51783616e+00,
    2.60869540e+00, 2.70283341e+00, 2.80036851e+00, 2.90142328e+00,
    3.00612474e+00, 3.11460448e+00, 3.22699885e+00, 3.34344911e+00,
    3.46410162e+00, 3.58910802e+00, 3.71862543e+00, 3.85281665e+00,
    3.99185032e+00, 4.13590119e+00, 4.28515031e+00, 4.43978528e+00,
    4.60000044e+00, 4.76599716e+00, 4.93798408e+00, 5.11617736e+00,
    5.30080097e+00, 5.49208695e+00, 5.69027572e+00, 5.89561638e+00,
    6.10836701e+00, 6.32879501e+00, 6.55717743e+00, 6.79380131e+00,
    7.03896406e+00, 7.29297381e+00, 7.55614982e+00, 7.82882286e+00,
    8.11133564e+00, 8.40404325e+00, 8.70731358e+00, 9.02152779e+00,
    9.34708081e+00, 9.68438182e+00, 1.00338548e+01, 1.03959388e+01,
    1.07710892e+01, 1.11597773e+01, 1.15624917e+01, 1.19797386e+01,
    1.24120423e+01, 1.28599463e+01, 1.33240134e+01, 1.38048269e+01,
    1.43029912e+01, 1.48191324e+01, 1.53538992e+01, 1.59079637e+01,
    1.64820223e+01, 1.70767965e+01, 1.76930339e+01, 1.83315090e+01,
    1.89930242e+01, 1.96784111e+01, 2.03885309e+01, 2.11242764e+01,
    2.18865721e+01, 2.26763762e+01, 2.34946814e+01, 2.43425161e+01,
    2.52209460e+01, 2.61310751e+01, 2.70740474e+01, 2.80510480e+01,
    2.90633048e+01, 3.01120902e+01, 3.11987224e+01, 3.23245670e+01,
    3.34910391e+01, 3.46996048e+01, 3.59517830e+01, 3.72491477e+01,
    3.85933294e+01, 3.99860176e+01, 4.14289627e+01, 4.29239782e+01,
    4.44729432e+01, 4.60778045e+01, 4.77405792e+01, 4.94633573e+01,
    5.12483039e+01, 5.30976625e+01, 5.50137575e+01, 5.69989972e+01,
    5.90558767e+01, 6.11869813e+01, 6.33949894e+01, 6.56826763e+01,
    6.80529171e+01, 7.05086911e+01, 7.30530847e+01, 7.56892959e+01,
    7.84206380e+01, 8.12505440e+01, 8.41825707e+01, 8.72204032e+01,
    9.03678597e+01, 9.36288960e+01, 9.70076108e+01, 1.00508251e+02,
    1.04135216e+02, 1.07893064e+02, 1.11786519e+02, 1.15820474e+02,
    1.20000000e+02
    };
   Int_t   m_binN   = sizeof(m_bins)/sizeof(Double_t) - 1; // # of Bins


   h_mee_OS  = new TH1F( "h_mee_OS", "Pair OS mass for e", m_binN, m_bins ) ; 
   h_mee_SS  = new TH1F( "h_mee_SS", "Pair SS mass for e", m_binN, m_bins ) ; 
   //h_mee_OS  = new TH1F( "h_mee_OS", "Pair OS mass for e", 1000, 0,100) ; 
   //h_mee_SS  = new TH1F( "h_mee_SS", "Pair SS mass for e", 1000, 0,100) ; 


   return kStOK ;

}
Int_t sDstAnaMaker::Make( )
{
   std::cout<< "sDstAnaMaker::Make" <<std::endl;

   // if parameter tree is not specified (or zero), connect the file
   // used to generate this class and read the Tree.
   if (!mRootInputFileName) {
      std::cout<< "ERROR: mRootInputFileName Not defined !!!" << std::endl;
      return 0;
   }
	std::cout << "Input FileList: " << mRootInputFileName << std::endl;
   ifstream fileListStream;
   fileListStream.open(mRootInputFileName);
   if(!fileListStream) { 
      cout << "NO_LIST_FILE " << mRootInputFileName << endl; 
      return 0;
   }
  string fileName =  "";

   //TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(mRootInputFileName);
   //if (!f || !f->IsOpen()) {
   //   f = new TFile(mRootInputFileName);
   //}
   //TTree *tree;
   //f->GetObject("T",tree);

   TChain *tree_chain = new TChain("T");
	while(getline(fileListStream, fileName)){        
    TString tmp = fileName;
    std::cout << "Input file: " << fileName << std::endl;
    tree_chain->Add(fileName.c_str());
   }   

   InitTree(tree_chain);   

   if (fChain == 0) return 0;

   Long64_t nentries = fChain->GetEntries();

   std::cout<< "-2" <<std::endl;
   std::cout<< "nentries = " << nentries <<std::endl;
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      //std::cout<< "-1" <<std::endl;
      //std::cout<< "jentry = " << jentry <<std::endl;

      //Long64_t ientry = LoadTree(jentry);
      //if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      //tree->GetEvent(jentry);
      // if (Cut(ientry) < 0) continue;
      
      if(n_fdpair[0]>0){
         //std::cout<< "runN = " << runN << " n_fdpair = "<< n_fdpair[0] <<std::endl;
         //std::cout<< "0" <<std::endl;
         for(int p=0; p<n_fdpair[0]; p++){
            //std::cout<< "1 p = " << p <<std::endl;
            double mee = mee_fdpair[p];
            if(q_fdpair[p]==0){
               //std::cout<< "2" <<std::endl;
               //std::cout<< "q = "<< q_fdpair[p] <<std::endl;
               //std::cout<< "Mee = "<<mee_fdpair[p] <<std::endl;
               
               h_mee_SS->Fill(mee);
            }else{
               //std::cout<< "3" <<std::endl;
               h_mee_OS->Fill(mee);
            }
         }
      }
      //else{
      //   continue;
      //}
   }

   return kStOK ;
 
}

Int_t sDstAnaMaker::Finish( )

{ // Do once at the end of the analysis
   std::cout<< "sDstAnaMaker::Finish" <<std::endl;
   printf("writting out\n");
   std::cout << "OUTPUT FILE NAME:" << mRootOutputFileName << std::endl;
   histogram_output = new TFile(mRootOutputFileName,"recreate");
  //TFile *outf = new TFile("./hist_collection.root","recreate");
  // Write histograms to disk, output miscellaneous other information
   h_mee_SS->Write();
   h_mee_OS->Write();
   histogram_output -> Write();
   histogram_output -> Close();
  //outf -> Close();// Write all histograms to disk 
  //cout << "Total Events Processed in DstMaker " << mEventsProcessed << endl ;

  return kStOK ;

}


Int_t sDstAnaMaker::TPL( )

{ // for making template
   std::cout<< "sDstAnaMaker::TPL" <<std::endl;

  //T->Fill(); // one entry

  return kStOK ;

}


