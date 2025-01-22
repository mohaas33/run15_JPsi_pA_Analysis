//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Jan 14 14:52:34 2025 by ROOT version 5.34/38
// from TTree T/uDst
// found on file: /gpfs01/star/pwg/eshulga/Files/skim_from_MuDST/9D615F903AC82A224460400731353AE6_21878.histograms.root
//////////////////////////////////////////////////////////

#ifndef sDstAnaMaker_def
#define sDstAnaMaker_def

#include "StMaker.h"
#include <TROOT.h>
#include <TChain.h>
//#include <TFile.h>

class TFile        ;
class TH1F         ;
class TH2F         ;

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class sDstAnaMaker : public StMaker
{
private:
   TString       mRootOutputFileName ;         //  Name of the histogram output file 
   TString       mRootInputFileName ;//  Name of the histogram output file 
   
   TFile*        histogram_output ;                  //  Histograms outputfile pointer
   
   TH1F* h_mee_OS;
   TH1F* h_mee_SS;

public :

   void SetOutputFileName(TString name) {
      mRootOutputFileName = name;
      std::cout << "Setting mRootOutputFileName to: " << mRootOutputFileName << std::endl;
   } // Make name available to member functions
   void SetInputFileName(TString name) {mRootInputFileName = name;} // Make name available to member functions


   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           runN[1];
   Int_t           eventN[1];
   Int_t           filln[1];
   Int_t           bid[1];
   Int_t           bid7[1];
   Double_t        bfield;
   Int_t           trig_Zerobias;
   Int_t           trig_Zdcmon;
   Int_t           trig_Bbcmon;
   Int_t           trig_UPCmain;
   Int_t           trig_UPCtopo;
   Int_t           trig_UPChighG;
   Int_t           trig_UPCjpsiB;
   Int_t           trig_RP2E;
   Int_t           trig_2E;
   Int_t           trig_RPUPC;
   Int_t           trig_UPCJPsi;
   Int_t           trig_UPCJPsizdc;
   Int_t           trig_UPCinc;
   Int_t           trig_JPsiHTTP;
   Int_t           lastDSM_TOFRP;
   Int_t           lastDSM_BBCZDC;
   Int_t           lastDSM_EMC;
   Int_t           zdce;
   Int_t           zdcw;
   Int_t           zdceadc[3];
   Int_t           zdcwadc[3];
   Int_t           zdcetdc;
   Int_t           zdcwtdc;
   Int_t           zdctdcdiff;
   Int_t           bbce;
   Int_t           bbcw;
   Int_t           ntoftrig;
   Double_t        zdcerate;
   Double_t        zdcwrate;
   Double_t        zdccrate;
   Int_t           nglobtrk;
   Int_t           nprimtrk;
   Int_t           nvtx;
   Int_t           n_MCvtx;
   Double_t        x_MCvtx[1];   //[n_MCvtx]
   Double_t        y_MCvtx[1];   //[n_MCvtx]
   Double_t        z_MCvtx[1];   //[n_MCvtx]
   Int_t           n_MCtrk;
   Int_t           iMCvtx_MCtrk[1];   //[n_MCtrk]
   Int_t           pid_MCtrk[1];   //[n_MCtrk]
   Int_t           q_MCtrk[1];   //[n_MCtrk]
   Double_t        px_MCtrk[1];   //[n_MCtrk]
   Double_t        py_MCtrk[1];   //[n_MCtrk]
   Double_t        pz_MCtrk[1];   //[n_MCtrk]
   Int_t           n_emccl;
   Double_t        e_emccl[36];   //[n_emccl]
   Double_t        eta_emccl[36];   //[n_emccl]
   Double_t        phi_emccl[36];   //[n_emccl]
   Double_t        sigeta_emccl[36];   //[n_emccl]
   Double_t        sigphi_emccl[36];   //[n_emccl]
   Int_t           nhits_emccl[36];   //[n_emccl]
   Int_t           idxmax_emccl[36];   //[n_emccl]
   Double_t        emax_emccl[36];   //[n_emccl]
   Int_t           adcmax_emccl[36];   //[n_emccl]
   Int_t           etabinmax_emccl[36];   //[n_emccl]
   Int_t           phibinmax_emccl[36];   //[n_emccl]
   Int_t           phiwdgmax_emccl[36];   //[n_emccl]
   Int_t           eemc_tower_cells[12][5][12];
   Double_t        eemc_prs1_energy[12][5][12];
   Double_t        eemc_prs2_energy[12][5][12];
   Double_t        eemc_tower_energy[12][5][12];
   Double_t        eemc_post_energy[12][5][12];
   Double_t        eemc_prs1_adc[12][5][12];
   Double_t        eemc_prs2_adc[12][5][12];
   Double_t        eemc_tower_adc[12][5][12];
   Double_t        eemc_post_adc[12][5][12];
   Int_t           n_fdtrk;
   Int_t           ifdvtx_fdtrk[50];   //[n_fdtrk]
   Int_t           q_fdtrk[50];   //[n_fdtrk]
   Double_t        pt_fdtrk[50];   //[n_fdtrk]
   Double_t        pz_fdtrk[50];   //[n_fdtrk]
   Double_t        p_fdtrk[50];   //[n_fdtrk]
   Double_t        eta_fdtrk[50];   //[n_fdtrk]
   Double_t        phi_fdtrk[50];   //[n_fdtrk]
   Int_t           nh_fdtrk[50];   //[n_fdtrk]
   Double_t        dcar_fdtrk[50];   //[n_fdtrk]
   Double_t        dcaz_fdtrk[50];   //[n_fdtrk]
   Double_t        dedx_fdtrk[50];   //[n_fdtrk]
   Int_t           nhdedx_fdtrk[50];   //[n_fdtrk]
   Double_t        sigel_fdtrk[50];   //[n_fdtrk]
   Double_t        sigpi_fdtrk[50];   //[n_fdtrk]
   Double_t        sigk_fdtrk[50];   //[n_fdtrk]
   Double_t        sigp_fdtrk[50];   //[n_fdtrk]
   Int_t           emcext_fdtrk[50];   //[n_fdtrk]
   Double_t        etaext_fdtrk[50];   //[n_fdtrk]
   Double_t        phiext_fdtrk[50];   //[n_fdtrk]
   Int_t           iemccl_fdtrk[50];   //[n_fdtrk]
   Int_t           eemc_ext_fdtrk[50];   //[n_fdtrk]
   Double_t        eemc_etaext_fdtrk[50];   //[n_fdtrk]
   Double_t        eemc_phiext_fdtrk[50];   //[n_fdtrk]
   Int_t           ieemc_cl_fdtrk[50];   //[n_fdtrk]
   Double_t        eemc_adc_fdtrk[50];   //[n_fdtrk]
   Double_t        eemc_energy_fdtrk[50];   //[n_fdtrk]
   Double_t        eemc_energy_prs1_fdtrk[50];   //[n_fdtrk]
   Double_t        eemc_energy_prs2_fdtrk[50];   //[n_fdtrk]
   Double_t        eemc_energy_tow_fdtrk[50];   //[n_fdtrk]
   Double_t        eemc_energy_post_fdtrk[50];   //[n_fdtrk]
   Double_t        eemc_adc_prs1_fdtrk[50];   //[n_fdtrk]
   Double_t        eemc_adc_prs2_fdtrk[50];   //[n_fdtrk]
   Double_t        eemc_adc_tow_fdtrk[50];   //[n_fdtrk]
   Double_t        eemc_adc_post_fdtrk[50];   //[n_fdtrk]
   Int_t           eemc_num_towers_in_track[50];   //[n_fdtrk]
   Int_t           tofhit_fdtrk[50];   //[n_fdtrk]
   Double_t        tofpatlen_fdtrk[50];   //[n_fdtrk]
   Double_t        toftof_fdtrk[50];   //[n_fdtrk]
   Double_t        tofbeta_fdtrk[50];   //[n_fdtrk]
   Double_t        toflet_fdtrk[50];   //[n_fdtrk]
   Double_t        toftet_fdtrk[50];   //[n_fdtrk]
   Int_t           n_fdvtx;
   Int_t           ntrk_fdvtx[8];   //[n_fdvtx]
   Int_t           nfdtrk_fdvtx[8];   //[n_fdvtx]
   Double_t        x_fdvtx[8];   //[n_fdvtx]
   Double_t        y_fdvtx[8];   //[n_fdvtx]
   Double_t        z_fdvtx[8];   //[n_fdvtx]
   Int_t           n_fdpair[1];
   Int_t           ifdtrk1_fdpair[354];   //[n_fdpair]
   Int_t           ifdtrk2_fdpair[354];   //[n_fdpair]
   Int_t           q_fdpair[354];   //[n_fdpair]
   Double_t        pt_fdpair[354];   //[n_fdpair]
   Double_t        pz_fdpair[354];   //[n_fdpair]
   Double_t        p_fdpair[354];   //[n_fdpair]
   Double_t        phi_fdpair[354];   //[n_fdpair]
   Double_t        acolin_fdpair[354];   //[n_fdpair]
   Double_t        mee_fdpair[354];   //[n_fdpair]
   Double_t        mmumu_fdpair[354];   //[n_fdpair]
   Double_t        mpipi_fdpair[354];   //[n_fdpair]
   Double_t        mkk_fdpair[354];   //[n_fdpair]
   Double_t        mpp_fdpair[354];   //[n_fdpair]
   Double_t        rapee_fdpair[354];   //[n_fdpair]
   Double_t        rapmumu_fdpair[354];   //[n_fdpair]
   Double_t        rappipi_fdpair[354];   //[n_fdpair]
   Double_t        rapkk_fdpair[354];   //[n_fdpair]
   Double_t        rappp_fdpair[354];   //[n_fdpair]
   Int_t           ntot_rptrk;
   Int_t           n_rptrk;
   Int_t           br_rptrk[100];   //[n_rptrk]
   Int_t           typ_rptrk[100];   //[n_rptrk]
   Int_t           np_rptrk[100];   //[n_rptrk]
   Int_t           npst1_rptrk[100];   //[n_rptrk]
   Int_t           qualst1_rptrk[100];   //[n_rptrk]
   Double_t        xst1_rptrk[100];   //[n_rptrk]
   Double_t        yst1_rptrk[100];   //[n_rptrk]
   Double_t        t1st1_rptrk[100];   //[n_rptrk]
   Double_t        t2st1_rptrk[100];   //[n_rptrk]
   Int_t           npst2_rptrk[100];   //[n_rptrk]
   Int_t           qualst2_rptrk[100];   //[n_rptrk]
   Double_t        xst2_rptrk[100];   //[n_rptrk]
   Double_t        yst2_rptrk[100];   //[n_rptrk]
   Double_t        t1st2_rptrk[100];   //[n_rptrk]
   Double_t        t2st2_rptrk[100];   //[n_rptrk]
   Double_t        pt_rptrk[100];   //[n_rptrk]
   Double_t        pz_rptrk[100];   //[n_rptrk]
   Double_t        phi_rptrk[100];   //[n_rptrk]
   Double_t        t_rptrk[100];   //[n_rptrk]
   Double_t        xi_rptrk[100];   //[n_rptrk]

   // List of branches
   TBranch        *b_runN;   //!
   TBranch        *b_eventN;   //!
   TBranch        *b_filln;   //!
   TBranch        *b_bid;   //!
   TBranch        *b_bid7;   //!
   TBranch        *b_bfield;   //!
   TBranch        *b_trig_Zerobias;   //!
   TBranch        *b_trig_Zdcmon;   //!
   TBranch        *b_trig_Bbcmon;   //!
   TBranch        *b_trig_UPCmain;   //!
   TBranch        *b_trig_UPCtopo;   //!
   TBranch        *b_trig_UPChighG;   //!
   TBranch        *b_trig_UPCjpsiB;   //!
   TBranch        *b_trig_RP2E;   //!
   TBranch        *b_trig_2E;   //!
   TBranch        *b_trig_RPUPC;   //!
   TBranch        *b_trig_UPCJPsi;   //!
   TBranch        *b_trig_UPCJPsizdc;   //!
   TBranch        *b_trig_UPCinc;   //!
   TBranch        *b_trig_JPsiHTTP;   //!
   TBranch        *b_lastDSM_TOFRP;   //!
   TBranch        *b_lastDSM_BBCZDC;   //!
   TBranch        *b_lastDSM_EMC;   //!
   TBranch        *b_zdce;   //!
   TBranch        *b_zdcw;   //!
   TBranch        *b_zdceadc;   //!
   TBranch        *b_zdcwadc;   //!
   TBranch        *b_zdcetdc;   //!
   TBranch        *b_zdcwtdc;   //!
   TBranch        *b_zdctdcdiff;   //!
   TBranch        *b_bbce;   //!
   TBranch        *b_bbcw;   //!
   TBranch        *b_ntoftrig;   //!
   TBranch        *b_zdcerate;   //!
   TBranch        *b_zdcwrate;   //!
   TBranch        *b_zdccrate;   //!
   TBranch        *b_nglobtrk;   //!
   TBranch        *b_nprimtrk;   //!
   TBranch        *b_nvtx;   //!
   TBranch        *b_n_MCvtx;   //!
   TBranch        *b_x_MCvtx;   //!
   TBranch        *b_y_MCvtx;   //!
   TBranch        *b_z_MCvtx;   //!
   TBranch        *b_n_MCtrk;   //!
   TBranch        *b_iMCvtx_MCtrk;   //!
   TBranch        *b_pid_MCtrk;   //!
   TBranch        *b_q_MCtrk;   //!
   TBranch        *b_px_MCtrk;   //!
   TBranch        *b_py_MCtrk;   //!
   TBranch        *b_pz_MCtrk;   //!
   TBranch        *b_n_emccl;   //!
   TBranch        *b_e_emccl;   //!
   TBranch        *b_eta_emccl;   //!
   TBranch        *b_phi_emccl;   //!
   TBranch        *b_sigeta_emccl;   //!
   TBranch        *b_sigphi_emccl;   //!
   TBranch        *b_nhits_emccl;   //!
   TBranch        *b_idxmax_emccl;   //!
   TBranch        *b_emax_emccl;   //!
   TBranch        *b_adcmax_emccl;   //!
   TBranch        *b_etabinmax_emccl;   //!
   TBranch        *b_phibinmax_emccl;   //!
   TBranch        *b_phiwdgmax_emccl;   //!
   TBranch        *b_eemc_tower_cells;   //!
   TBranch        *b_eemc_prs1_energy;   //!
   TBranch        *b_eemc_prs2_energy;   //!
   TBranch        *b_eemc_tower_energy;   //!
   TBranch        *b_eemc_post_energy;   //!
   TBranch        *b_eemc_prs1_adc;   //!
   TBranch        *b_eemc_prs2_adc;   //!
   TBranch        *b_eemc_tower_adc;   //!
   TBranch        *b_eemc_post_adc;   //!
   TBranch        *b_n_fdtrk;   //!
   TBranch        *b_ifdvtx_fdtrk;   //!
   TBranch        *b_q_fdtrk;   //!
   TBranch        *b_pt_fdtrk;   //!
   TBranch        *b_pz_fdtrk;   //!
   TBranch        *b_p_fdtrk;   //!
   TBranch        *b_eta_fdtrk;   //!
   TBranch        *b_phi_fdtrk;   //!
   TBranch        *b_nh_fdtrk;   //!
   TBranch        *b_dcar_fdtrk;   //!
   TBranch        *b_dcaz_fdtrk;   //!
   TBranch        *b_dedx_fdtrk;   //!
   TBranch        *b_nhdedx_fdtrk;   //!
   TBranch        *b_sigel_fdtrk;   //!
   TBranch        *b_sigpi_fdtrk;   //!
   TBranch        *b_sigk_fdtrk;   //!
   TBranch        *b_sigp_fdtrk;   //!
   TBranch        *b_emcext_fdtrk;   //!
   TBranch        *b_etaext_fdtrk;   //!
   TBranch        *b_phiext_fdtrk;   //!
   TBranch        *b_iemccl_fdtrk;   //!
   TBranch        *b_eemc_ext_fdtrk;   //!
   TBranch        *b_eemc_etaext_fdtrk;   //!
   TBranch        *b_eemc_phiext_fdtrk;   //!
   TBranch        *b_ieemc_cl_fdtrk;   //!
   TBranch        *b_eemc_adc_fdtrk;   //!
   TBranch        *b_eemc_energy_fdtrk;   //!
   TBranch        *b_eemc_energy_prs1_fdtrk;   //!
   TBranch        *b_eemc_energy_prs2_fdtrk;   //!
   TBranch        *b_eemc_energy_tow_fdtrk;   //!
   TBranch        *b_eemc_energy_post_fdtrk;   //!
   TBranch        *b_eemc_adc_prs1_fdtrk;   //!
   TBranch        *b_eemc_adc_prs2_fdtrk;   //!
   TBranch        *b_eemc_adc_tow_fdtrk;   //!
   TBranch        *b_eemc_adc_post_fdtrk;   //!
   TBranch        *b_eemc_num_towers_in_track;   //!
   TBranch        *b_tofhit_fdtrk;   //!
   TBranch        *b_tofpatlen_fdtrk;   //!
   TBranch        *b_toftof_fdtrk;   //!
   TBranch        *b_tofbeta_fdtrk;   //!
   TBranch        *b_toflet_fdtrk;   //!
   TBranch        *b_toftet_fdtrk;   //!
   TBranch        *b_n_fdvtx;   //!
   TBranch        *b_ntrk_fdvtx;   //!
   TBranch        *b_nfdtrk_fdvtx;   //!
   TBranch        *b_x_fdvtx;   //!
   TBranch        *b_y_fdvtx;   //!
   TBranch        *b_z_fdvtx;   //!
   TBranch        *b_n_fdpair;   //!
   TBranch        *b_ifdtrk1_fdpair;   //!
   TBranch        *b_ifdtrk2_fdpair;   //!
   TBranch        *b_q_fdpair;   //!
   TBranch        *b_pt_fdpair;   //!
   TBranch        *b_pz_fdpair;   //!
   TBranch        *b_p_fdpair;   //!
   TBranch        *b_phi_fdpair;   //!
   TBranch        *b_acolin_fdpair;   //!
   TBranch        *b_mee_fdpair;   //!
   TBranch        *b_mmumu_fdpair;   //!
   TBranch        *b_mpipi_fdpair;   //!
   TBranch        *b_mkk_fdpair;   //!
   TBranch        *b_mpp_fdpair;   //!
   TBranch        *b_rapee_fdpair;   //!
   TBranch        *b_rapmumu_fdpair;   //!
   TBranch        *b_rappipi_fdpair;   //!
   TBranch        *b_rapkk_fdpair;   //!
   TBranch        *b_rappp_fdpair;   //!
   TBranch        *b_ntot_rptrk;   //!
   TBranch        *b_n_rptrk;   //!
   TBranch        *b_br_rptrk;   //!
   TBranch        *b_typ_rptrk;   //!
   TBranch        *b_np_rptrk;   //!
   TBranch        *b_npst1_rptrk;   //!
   TBranch        *b_qualst1_rptrk;   //!
   TBranch        *b_xst1_rptrk;   //!
   TBranch        *b_yst1_rptrk;   //!
   TBranch        *b_t1st1_rptrk;   //!
   TBranch        *b_t2st1_rptrk;   //!
   TBranch        *b_npst2_rptrk;   //!
   TBranch        *b_qualst2_rptrk;   //!
   TBranch        *b_xst2_rptrk;   //!
   TBranch        *b_yst2_rptrk;   //!
   TBranch        *b_t1st2_rptrk;   //!
   TBranch        *b_t2st2_rptrk;   //!
   TBranch        *b_pt_rptrk;   //!
   TBranch        *b_pz_rptrk;   //!
   TBranch        *b_phi_rptrk;   //!
   TBranch        *b_t_rptrk;   //!
   TBranch        *b_xi_rptrk;   //!

   sDstAnaMaker();
   virtual ~sDstAnaMaker();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     InitTree(TChain *tree);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   Int_t Init    ( ) ;                              //  Initiliaze the analysis tools ... done once
   Int_t Make    ( ) ;                              //  The main analysis that is done on each event
   Int_t Finish  ( ) ;                              //  Finish the analysis, close files, and clean up.
   Int_t TPL     ( ) ;                              //  Call T->Fill() w/o X_br_fill, for template maker


   ClassDef(sDstAnaMaker,1)                  //  Macro for CINT compatability

};

#endif

#ifdef sDstAnaMaker_cxx
sDstAnaMaker::sDstAnaMaker() : StMaker("uDstSkimMaker")
{
   h_mee_OS = NULL ;
   h_mee_SS = NULL ;
   histogram_output = NULL  ;
   
   mRootOutputFileName = "" ;               // Histogram Output File Name will be set inside the "analysis".C macro
   mRootInputFileName = "" ;               // Histogram Input File Name will be set inside the "analysis".C macro

   Init();
}

sDstAnaMaker::~sDstAnaMaker()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t sDstAnaMaker::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t sDstAnaMaker::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void sDstAnaMaker::InitTree(TChain *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("runN", &runN, &b_runN);
   fChain->SetBranchAddress("eventN", &eventN, &b_eventN);
   fChain->SetBranchAddress("filln", &filln, &b_filln);
   fChain->SetBranchAddress("bid", &bid, &b_bid);
   fChain->SetBranchAddress("bid7", &bid7, &b_bid7);
   fChain->SetBranchAddress("bfield", &bfield, &b_bfield);
   fChain->SetBranchAddress("trig_Zerobias", &trig_Zerobias, &b_trig_Zerobias);
   fChain->SetBranchAddress("trig_Zdcmon", &trig_Zdcmon, &b_trig_Zdcmon);
   fChain->SetBranchAddress("trig_Bbcmon", &trig_Bbcmon, &b_trig_Bbcmon);
   fChain->SetBranchAddress("trig_UPCmain", &trig_UPCmain, &b_trig_UPCmain);
   fChain->SetBranchAddress("trig_UPCtopo", &trig_UPCtopo, &b_trig_UPCtopo);
   fChain->SetBranchAddress("trig_UPChighG", &trig_UPChighG, &b_trig_UPChighG);
   fChain->SetBranchAddress("trig_UPCjpsiB", &trig_UPCjpsiB, &b_trig_UPCjpsiB);
   fChain->SetBranchAddress("trig_RP2E", &trig_RP2E, &b_trig_RP2E);
   fChain->SetBranchAddress("trig_2E", &trig_2E, &b_trig_2E);
   fChain->SetBranchAddress("trig_RPUPC", &trig_RPUPC, &b_trig_RPUPC);
   fChain->SetBranchAddress("trig_UPCJPsi", &trig_UPCJPsi, &b_trig_UPCJPsi);
   fChain->SetBranchAddress("trig_UPCJPsizdc", &trig_UPCJPsizdc, &b_trig_UPCJPsizdc);
   fChain->SetBranchAddress("trig_UPCinc", &trig_UPCinc, &b_trig_UPCinc);
   fChain->SetBranchAddress("trig_JPsiHTTP", &trig_JPsiHTTP, &b_trig_JPsiHTTP);
   fChain->SetBranchAddress("lastDSM_TOFRP", &lastDSM_TOFRP, &b_lastDSM_TOFRP);
   fChain->SetBranchAddress("lastDSM_BBCZDC", &lastDSM_BBCZDC, &b_lastDSM_BBCZDC);
   fChain->SetBranchAddress("lastDSM_EMC", &lastDSM_EMC, &b_lastDSM_EMC);
   fChain->SetBranchAddress("zdce", &zdce, &b_zdce);
   fChain->SetBranchAddress("zdcw", &zdcw, &b_zdcw);
   fChain->SetBranchAddress("zdceadc", zdceadc, &b_zdceadc);
   fChain->SetBranchAddress("zdcwadc", zdcwadc, &b_zdcwadc);
   fChain->SetBranchAddress("zdcetdc", &zdcetdc, &b_zdcetdc);
   fChain->SetBranchAddress("zdcwtdc", &zdcwtdc, &b_zdcwtdc);
   fChain->SetBranchAddress("zdctdcdiff", &zdctdcdiff, &b_zdctdcdiff);
   fChain->SetBranchAddress("bbce", &bbce, &b_bbce);
   fChain->SetBranchAddress("bbcw", &bbcw, &b_bbcw);
   fChain->SetBranchAddress("ntoftrig", &ntoftrig, &b_ntoftrig);
   fChain->SetBranchAddress("zdcerate", &zdcerate, &b_zdcerate);
   fChain->SetBranchAddress("zdcwrate", &zdcwrate, &b_zdcwrate);
   fChain->SetBranchAddress("zdccrate", &zdccrate, &b_zdccrate);
   fChain->SetBranchAddress("nglobtrk", &nglobtrk, &b_nglobtrk);
   fChain->SetBranchAddress("nprimtrk", &nprimtrk, &b_nprimtrk);
   fChain->SetBranchAddress("nvtx", &nvtx, &b_nvtx);
   fChain->SetBranchAddress("n_MCvtx", &n_MCvtx, &b_n_MCvtx);
   fChain->SetBranchAddress("x_MCvtx", &x_MCvtx, &b_x_MCvtx);
   fChain->SetBranchAddress("y_MCvtx", &y_MCvtx, &b_y_MCvtx);
   fChain->SetBranchAddress("z_MCvtx", &z_MCvtx, &b_z_MCvtx);
   fChain->SetBranchAddress("n_MCtrk", &n_MCtrk, &b_n_MCtrk);
   fChain->SetBranchAddress("iMCvtx_MCtrk", &iMCvtx_MCtrk, &b_iMCvtx_MCtrk);
   fChain->SetBranchAddress("pid_MCtrk", &pid_MCtrk, &b_pid_MCtrk);
   fChain->SetBranchAddress("q_MCtrk", &q_MCtrk, &b_q_MCtrk);
   fChain->SetBranchAddress("px_MCtrk", &px_MCtrk, &b_px_MCtrk);
   fChain->SetBranchAddress("py_MCtrk", &py_MCtrk, &b_py_MCtrk);
   fChain->SetBranchAddress("pz_MCtrk", &pz_MCtrk, &b_pz_MCtrk);
   fChain->SetBranchAddress("n_emccl", &n_emccl, &b_n_emccl);
   fChain->SetBranchAddress("e_emccl", e_emccl, &b_e_emccl);
   fChain->SetBranchAddress("eta_emccl", eta_emccl, &b_eta_emccl);
   fChain->SetBranchAddress("phi_emccl", phi_emccl, &b_phi_emccl);
   fChain->SetBranchAddress("sigeta_emccl", sigeta_emccl, &b_sigeta_emccl);
   fChain->SetBranchAddress("sigphi_emccl", sigphi_emccl, &b_sigphi_emccl);
   fChain->SetBranchAddress("nhits_emccl", nhits_emccl, &b_nhits_emccl);
   fChain->SetBranchAddress("idxmax_emccl", idxmax_emccl, &b_idxmax_emccl);
   fChain->SetBranchAddress("emax_emccl", emax_emccl, &b_emax_emccl);
   fChain->SetBranchAddress("adcmax_emccl", adcmax_emccl, &b_adcmax_emccl);
   fChain->SetBranchAddress("etabinmax_emccl", etabinmax_emccl, &b_etabinmax_emccl);
   fChain->SetBranchAddress("phibinmax_emccl", phibinmax_emccl, &b_phibinmax_emccl);
   fChain->SetBranchAddress("phiwdgmax_emccl", phiwdgmax_emccl, &b_phiwdgmax_emccl);
   fChain->SetBranchAddress("eemc_tower_cells", eemc_tower_cells, &b_eemc_tower_cells);
   fChain->SetBranchAddress("eemc_prs1_energy", eemc_prs1_energy, &b_eemc_prs1_energy);
   fChain->SetBranchAddress("eemc_prs2_energy", eemc_prs2_energy, &b_eemc_prs2_energy);
   fChain->SetBranchAddress("eemc_tower_energy", eemc_tower_energy, &b_eemc_tower_energy);
   fChain->SetBranchAddress("eemc_post_energy", eemc_post_energy, &b_eemc_post_energy);
   fChain->SetBranchAddress("eemc_prs1_adc", eemc_prs1_adc, &b_eemc_prs1_adc);
   fChain->SetBranchAddress("eemc_prs2_adc", eemc_prs2_adc, &b_eemc_prs2_adc);
   fChain->SetBranchAddress("eemc_tower_adc", eemc_tower_adc, &b_eemc_tower_adc);
   fChain->SetBranchAddress("eemc_post_adc", eemc_post_adc, &b_eemc_post_adc);
   fChain->SetBranchAddress("n_fdtrk", &n_fdtrk, &b_n_fdtrk);
   fChain->SetBranchAddress("ifdvtx_fdtrk", ifdvtx_fdtrk, &b_ifdvtx_fdtrk);
   fChain->SetBranchAddress("q_fdtrk", q_fdtrk, &b_q_fdtrk);
   fChain->SetBranchAddress("pt_fdtrk", pt_fdtrk, &b_pt_fdtrk);
   fChain->SetBranchAddress("pz_fdtrk", pz_fdtrk, &b_pz_fdtrk);
   fChain->SetBranchAddress("p_fdtrk", p_fdtrk, &b_p_fdtrk);
   fChain->SetBranchAddress("eta_fdtrk", eta_fdtrk, &b_eta_fdtrk);
   fChain->SetBranchAddress("phi_fdtrk", phi_fdtrk, &b_phi_fdtrk);
   fChain->SetBranchAddress("nh_fdtrk", nh_fdtrk, &b_nh_fdtrk);
   fChain->SetBranchAddress("dcar_fdtrk", dcar_fdtrk, &b_dcar_fdtrk);
   fChain->SetBranchAddress("dcaz_fdtrk", dcaz_fdtrk, &b_dcaz_fdtrk);
   fChain->SetBranchAddress("dedx_fdtrk", dedx_fdtrk, &b_dedx_fdtrk);
   fChain->SetBranchAddress("nhdedx_fdtrk", nhdedx_fdtrk, &b_nhdedx_fdtrk);
   fChain->SetBranchAddress("sigel_fdtrk", sigel_fdtrk, &b_sigel_fdtrk);
   fChain->SetBranchAddress("sigpi_fdtrk", sigpi_fdtrk, &b_sigpi_fdtrk);
   fChain->SetBranchAddress("sigk_fdtrk", sigk_fdtrk, &b_sigk_fdtrk);
   fChain->SetBranchAddress("sigp_fdtrk", sigp_fdtrk, &b_sigp_fdtrk);
   fChain->SetBranchAddress("emcext_fdtrk", emcext_fdtrk, &b_emcext_fdtrk);
   fChain->SetBranchAddress("etaext_fdtrk", etaext_fdtrk, &b_etaext_fdtrk);
   fChain->SetBranchAddress("phiext_fdtrk", phiext_fdtrk, &b_phiext_fdtrk);
   fChain->SetBranchAddress("iemccl_fdtrk", iemccl_fdtrk, &b_iemccl_fdtrk);
   fChain->SetBranchAddress("eemc_ext_fdtrk", eemc_ext_fdtrk, &b_eemc_ext_fdtrk);
   fChain->SetBranchAddress("eemc_etaext_fdtrk", eemc_etaext_fdtrk, &b_eemc_etaext_fdtrk);
   fChain->SetBranchAddress("eemc_phiext_fdtrk", eemc_phiext_fdtrk, &b_eemc_phiext_fdtrk);
   fChain->SetBranchAddress("ieemc_cl_fdtrk", ieemc_cl_fdtrk, &b_ieemc_cl_fdtrk);
   fChain->SetBranchAddress("eemc_adc_fdtrk", eemc_adc_fdtrk, &b_eemc_adc_fdtrk);
   fChain->SetBranchAddress("eemc_energy_fdtrk", eemc_energy_fdtrk, &b_eemc_energy_fdtrk);
   fChain->SetBranchAddress("eemc_energy_prs1_fdtrk", eemc_energy_prs1_fdtrk, &b_eemc_energy_prs1_fdtrk);
   fChain->SetBranchAddress("eemc_energy_prs2_fdtrk", eemc_energy_prs2_fdtrk, &b_eemc_energy_prs2_fdtrk);
   fChain->SetBranchAddress("eemc_energy_tow_fdtrk", eemc_energy_tow_fdtrk, &b_eemc_energy_tow_fdtrk);
   fChain->SetBranchAddress("eemc_energy_post_fdtrk", eemc_energy_post_fdtrk, &b_eemc_energy_post_fdtrk);
   fChain->SetBranchAddress("eemc_adc_prs1_fdtrk", eemc_adc_prs1_fdtrk, &b_eemc_adc_prs1_fdtrk);
   fChain->SetBranchAddress("eemc_adc_prs2_fdtrk", eemc_adc_prs2_fdtrk, &b_eemc_adc_prs2_fdtrk);
   fChain->SetBranchAddress("eemc_adc_tow_fdtrk", eemc_adc_tow_fdtrk, &b_eemc_adc_tow_fdtrk);
   fChain->SetBranchAddress("eemc_adc_post_fdtrk", eemc_adc_post_fdtrk, &b_eemc_adc_post_fdtrk);
   fChain->SetBranchAddress("eemc_num_towers_in_track", eemc_num_towers_in_track, &b_eemc_num_towers_in_track);
   fChain->SetBranchAddress("tofhit_fdtrk", tofhit_fdtrk, &b_tofhit_fdtrk);
   fChain->SetBranchAddress("tofpatlen_fdtrk", tofpatlen_fdtrk, &b_tofpatlen_fdtrk);
   fChain->SetBranchAddress("toftof_fdtrk", toftof_fdtrk, &b_toftof_fdtrk);
   fChain->SetBranchAddress("tofbeta_fdtrk", tofbeta_fdtrk, &b_tofbeta_fdtrk);
   fChain->SetBranchAddress("toflet_fdtrk", toflet_fdtrk, &b_toflet_fdtrk);
   fChain->SetBranchAddress("toftet_fdtrk", toftet_fdtrk, &b_toftet_fdtrk);
   fChain->SetBranchAddress("n_fdvtx", &n_fdvtx, &b_n_fdvtx);
   fChain->SetBranchAddress("ntrk_fdvtx", ntrk_fdvtx, &b_ntrk_fdvtx);
   fChain->SetBranchAddress("nfdtrk_fdvtx", nfdtrk_fdvtx, &b_nfdtrk_fdvtx);
   fChain->SetBranchAddress("x_fdvtx", x_fdvtx, &b_x_fdvtx);
   fChain->SetBranchAddress("y_fdvtx", y_fdvtx, &b_y_fdvtx);
   fChain->SetBranchAddress("z_fdvtx", z_fdvtx, &b_z_fdvtx);
   fChain->SetBranchAddress("n_fdpair", &n_fdpair, &b_n_fdpair);
   fChain->SetBranchAddress("ifdtrk1_fdpair", ifdtrk1_fdpair, &b_ifdtrk1_fdpair);
   fChain->SetBranchAddress("ifdtrk2_fdpair", ifdtrk2_fdpair, &b_ifdtrk2_fdpair);
   fChain->SetBranchAddress("q_fdpair", q_fdpair, &b_q_fdpair);
   fChain->SetBranchAddress("pt_fdpair", pt_fdpair, &b_pt_fdpair);
   fChain->SetBranchAddress("pz_fdpair", pz_fdpair, &b_pz_fdpair);
   fChain->SetBranchAddress("p_fdpair", p_fdpair, &b_p_fdpair);
   fChain->SetBranchAddress("phi_fdpair", phi_fdpair, &b_phi_fdpair);
   fChain->SetBranchAddress("acolin_fdpair", acolin_fdpair, &b_acolin_fdpair);
   fChain->SetBranchAddress("mee_fdpair", mee_fdpair, &b_mee_fdpair);
   fChain->SetBranchAddress("mmumu_fdpair", mmumu_fdpair, &b_mmumu_fdpair);
   fChain->SetBranchAddress("mpipi_fdpair", mpipi_fdpair, &b_mpipi_fdpair);
   fChain->SetBranchAddress("mkk_fdpair", mkk_fdpair, &b_mkk_fdpair);
   fChain->SetBranchAddress("mpp_fdpair", mpp_fdpair, &b_mpp_fdpair);
   fChain->SetBranchAddress("rapee_fdpair", rapee_fdpair, &b_rapee_fdpair);
   fChain->SetBranchAddress("rapmumu_fdpair", rapmumu_fdpair, &b_rapmumu_fdpair);
   fChain->SetBranchAddress("rappipi_fdpair", rappipi_fdpair, &b_rappipi_fdpair);
   fChain->SetBranchAddress("rapkk_fdpair", rapkk_fdpair, &b_rapkk_fdpair);
   fChain->SetBranchAddress("rappp_fdpair", rappp_fdpair, &b_rappp_fdpair);
   fChain->SetBranchAddress("ntot_rptrk", &ntot_rptrk, &b_ntot_rptrk);
   fChain->SetBranchAddress("n_rptrk", &n_rptrk, &b_n_rptrk);
   fChain->SetBranchAddress("br_rptrk", br_rptrk, &b_br_rptrk);
   fChain->SetBranchAddress("typ_rptrk", typ_rptrk, &b_typ_rptrk);
   fChain->SetBranchAddress("np_rptrk", np_rptrk, &b_np_rptrk);
   fChain->SetBranchAddress("npst1_rptrk", npst1_rptrk, &b_npst1_rptrk);
   fChain->SetBranchAddress("qualst1_rptrk", qualst1_rptrk, &b_qualst1_rptrk);
   fChain->SetBranchAddress("xst1_rptrk", xst1_rptrk, &b_xst1_rptrk);
   fChain->SetBranchAddress("yst1_rptrk", yst1_rptrk, &b_yst1_rptrk);
   fChain->SetBranchAddress("t1st1_rptrk", t1st1_rptrk, &b_t1st1_rptrk);
   fChain->SetBranchAddress("t2st1_rptrk", t2st1_rptrk, &b_t2st1_rptrk);
   fChain->SetBranchAddress("npst2_rptrk", npst2_rptrk, &b_npst2_rptrk);
   fChain->SetBranchAddress("qualst2_rptrk", qualst2_rptrk, &b_qualst2_rptrk);
   fChain->SetBranchAddress("xst2_rptrk", xst2_rptrk, &b_xst2_rptrk);
   fChain->SetBranchAddress("yst2_rptrk", yst2_rptrk, &b_yst2_rptrk);
   fChain->SetBranchAddress("t1st2_rptrk", t1st2_rptrk, &b_t1st2_rptrk);
   fChain->SetBranchAddress("t2st2_rptrk", t2st2_rptrk, &b_t2st2_rptrk);
   fChain->SetBranchAddress("pt_rptrk", pt_rptrk, &b_pt_rptrk);
   fChain->SetBranchAddress("pz_rptrk", pz_rptrk, &b_pz_rptrk);
   fChain->SetBranchAddress("phi_rptrk", phi_rptrk, &b_phi_rptrk);
   fChain->SetBranchAddress("t_rptrk", t_rptrk, &b_t_rptrk);
   fChain->SetBranchAddress("xi_rptrk", xi_rptrk, &b_xi_rptrk);
   Notify();
}

Bool_t sDstAnaMaker::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void sDstAnaMaker::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t sDstAnaMaker::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef sDstAnaMaker_cxx
