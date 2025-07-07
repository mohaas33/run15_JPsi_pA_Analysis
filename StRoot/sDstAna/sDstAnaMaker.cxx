#define sDstAnaMaker_cxx
#include "sDstAnaMaker.h"

#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TObjArray.h"

int bemc_phi_wedge(double phi_ext);
int eemc_phi_wedge(double phi_ext);
int both_phi_wedge(double phi_ext);
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

   char name[100];

   //h_muu_OS  = new TH1F( "h_muu_OS", "Pair OS mass for #mu", m_binN, m_bins ) ; 
   //h_muu_SS  = new TH1F( "h_muu_SS", "Pair SS mass for #mu", m_binN, m_bins ) ; 
   h_muu_OS  = new TH1F( "h_muu_OS", "Pair OS mass for #mu", 1000, 0, 10) ; 
   h_muu_SS  = new TH1F( "h_muu_SS", "Pair SS mass for #mu", 1000, 0, 10) ; 
   h_mpiP_OS  = new TH1F( "h_mpiP_OS", "Pair OS mass for #pip", 1000, 0, 10) ; 
   h_mpiP_SS  = new TH1F( "h_mpiP_SS", "Pair SS mass for #pip", 1000, 0, 10) ; 

   h_mkk_OS    = new TH1F( "h_mkk_OS", "Pair OS mass for K", 1000, 0, 10) ; 
   h_mkk_SS    = new TH1F( "h_mkk_SS", "Pair SS mass for K", 1000, 0, 10) ; 
   h_mpp_OS    = new TH1F( "h_mpp_OS", "Pair OS mass for p", 1000, 0, 10) ; 
   h_mpp_SS    = new TH1F( "h_mpp_SS", "Pair SS mass for p", 1000, 0, 10) ; 
   h_mpipi_OS  = new TH1F( "h_mpipi_OS", "Pair OS mass for #pi", 1000, 0, 10) ; 
   h_mpipi_SS  = new TH1F( "h_mpipi_SS", "Pair SS mass for #pi", 1000, 0, 10) ; 

   h_mee_OS  = new TH1F( "h_mee_OS", "Pair OS mass for e", m_binN, m_bins ) ; 
   h_mee_SS  = new TH1F( "h_mee_SS", "Pair SS mass for e", m_binN, m_bins ) ; 
   h_y_OS  = new TH1F( "h_y_OS", "Pair OS y for e", 200, -2, 2 ) ; 
   h_y_SS  = new TH1F( "h_y_SS", "Pair SS y for e", 200, -2, 2 ) ; 

   // Main parameters 
	// - number of TPC hits
   // - number of fit hits De/dx
   // - BEMC Clusters
   // - various Chi2
   // - tower energy   

   h_pair_bemc_wedge_b2b = new TH1F("h_pair_bemc_wedge_b2b","h_pair_bemc_wedge_b2b", 100, -4.5, 95.5);
   for(int i=0;i<4;i++){
      sprintf(name,"h_spin_mee_trkCuts_%d",i);
      h_spin_mee_trkCuts[i]       = new TH2F(name,name, 20, -3.5, 16.5, 40, 2,6);
      sprintf(name,"h_spin_mee_SS_trkCuts_%d",i);
      h_spin_mee_SS_trkCuts[i]       = new TH2F(name,name, 20, -3.5, 16.5, 40, 2,6);
   }
   h_evt = new TH1F("h_evt","h_evt",20,-2.5,17.5);
   h_nvtx = new TH1F("h_nvtx","h_nvtx",20,-0.5,19.5);
   h_nvtx_rank = new TH1F("h_nvtx_rank","h_nvtx_rank",20,-0.5,19.5);
   h_evt_z_rank = new TH1F("h_evt_z_rank","h_evt_z_rank",800,-200,200);
   h_evt_z_high_rank = new TH1F("h_evt_z_high_rank","h_evt_z_high_rank",800,-200,200);
   h_evt_z = new TH1F("h_evt_z","h_evt_z",800,-200,200);
   h_evt_x = new TH1F("h_evt_x","h_evt_x",800,-200,200);
   h_evt_y = new TH1F("h_evt_y","h_evt_y",800,-200,200);
   h_pair_evt_z = new TH1F("h_pair_evt_z","h_pair_evt_z",800,-200,200);

   h_evt_ntoftrig = new TH1F("h_evt_ntoftrig","h_evt_ntoftrig",200,-0.5,199.5);
   h_zb_evt_ntoftrig = new TH1F("h_zb_evt_ntoftrig","h_zb_evt_ntoftrig",200,-0.5,199.5);


   h_vtx1_vtx2 = new TH2F("h_vtx1_vtx2","h_vtx1_vtx2",15,-1.5,13.5,15,-1.5,13.5);

   h_evt_spin_z   = new TH2F("h_evt_spin_z","h_evt_spin_z", 20, -3.5, 16.5, 800, -200,200);
   h_evt_spin     = new TH2F("h_evt_spin","h_evt_spin", 20, -3.5, 16.5, 200, -0.5,199.5);
   h_evt_spin_mw  = new TH2F("h_evt_spin_mw","h_evt_spin_mw", 20, -3.5, 16.5, 200, -0.5,199.5);

   h_spin_mee_cuts       = new TH2F("h_spin_mee_cuts","h_spin_mee_cuts", 20, -3.5, 16.5, 40, 2,6);
   h_spin_muu_cuts       = new TH2F("h_spin_muu_cuts","h_spin_muu_cuts", 20, -3.5, 16.5, 40, 2,6);
   h_spin_mee       = new TH2F("h_spin_mee","h_spin_mee", 20, -3.5, 16.5, 40, 2,6);
   h_spin_mee_R       = new TH2F("h_spin_mee_R","h_spin_mee_R", 20, -3.5, 16.5, 40, 2,6);
   h_spin_mee_L       = new TH2F("h_spin_mee_L","h_spin_mee_L", 20, -3.5, 16.5, 40, 2,6);
   h_spin_mee_SS       = new TH2F("h_spin_mee_SS","h_spin_mee_SS", 20, -3.5, 16.5, 40, 2,6);
   h_spin_mee_SS_R       = new TH2F("h_spin_mee_SS_R","h_spin_mee_SS_R", 20, -3.5, 16.5, 40, 2,6);
   h_spin_mee_SS_L       = new TH2F("h_spin_mee_SS_L","h_spin_mee_SS_L", 20, -3.5, 16.5, 40, 2,6);

   h_dca_spin_mee       = new TH2F("h_dca_spin_mee",     "h_dca_spin_mee", 20, -3.5, 16.5, 40, 2,6);
   h_dca_spin_mee_R     = new TH2F("h_dca_spin_mee_R",   "h_dca_spin_mee_R", 20, -3.5, 16.5, 40, 2,6);
   h_dca_spin_mee_L     = new TH2F("h_dca_spin_mee_L",   "h_dca_spin_mee_L", 20, -3.5, 16.5, 40, 2,6);
   h_dca_spin_mee_SS    = new TH2F("h_dca_spin_mee_SS",  "h_dca_spin_mee_SS", 20, -3.5, 16.5, 40, 2,6);
   h_dca_spin_mee_SS_R  = new TH2F("h_dca_spin_mee_SS_R","h_dca_spin_mee_SS_R", 20, -3.5, 16.5, 40, 2,6);
   h_dca_spin_mee_SS_L  = new TH2F("h_dca_spin_mee_SS_L","h_dca_spin_mee_SS_L", 20, -3.5, 16.5, 40, 2,6);

   h_RP2E_dca_spin_mee       = new TH2F("h_RP2E_dca_spin_mee",     "h_RP2E_dca_spin_mee", 20, -3.5, 16.5, 40, 2,6);
   h_RP2E_dca_spin_mee_R     = new TH2F("h_RP2E_dca_spin_mee_R",   "h_RP2E_dca_spin_mee_R", 20, -3.5, 16.5, 40, 2,6);
   h_RP2E_dca_spin_mee_L     = new TH2F("h_RP2E_dca_spin_mee_L",   "h_RP2E_dca_spin_mee_L", 20, -3.5, 16.5, 40, 2,6);
   h_RP2E_dca_spin_mee_SS    = new TH2F("h_RP2E_dca_spin_mee_SS",  "h_RP2E_dca_spin_mee_SS", 20, -3.5, 16.5, 40, 2,6);
   h_RP2E_dca_spin_mee_SS_R  = new TH2F("h_RP2E_dca_spin_mee_SS_R","h_RP2E_dca_spin_mee_SS_R", 20, -3.5, 16.5, 40, 2,6);
   h_RP2E_dca_spin_mee_SS_L  = new TH2F("h_RP2E_dca_spin_mee_SS_L","h_RP2E_dca_spin_mee_SS_L", 20, -3.5, 16.5, 40, 2,6);

   h_spin_dcaz       = new TH2F("h_spin_dcaz","h_spin_dcaz", 20, -3.5, 16.5, 40, -2,2);
   h_spin_dcar       = new TH2F("h_spin_dcar","h_spin_dcar", 20, -3.5, 16.5, 40, 0,4);
   h_spin_dcaz_SS    = new TH2F("h_spin_dcaz_SS","h_spin_dcaz_SS", 20, -3.5, 16.5, 40, -2,2);
   h_spin_dcar_SS    = new TH2F("h_spin_dcar_SS","h_spin_dcar_SS", 20, -3.5, 16.5, 40, 0,4);
   h_mass_window_spin_dcaz       = new TH2F("h_mass_window_spin_dcaz","h_mass_window_spin_dcaz", 20, -3.5, 16.5, 40, -2,2);
   h_mass_window_spin_dcar       = new TH2F("h_mass_window_spin_dcar","h_mass_window_spin_dcar", 20, -3.5, 16.5, 40, 0,4);
   h_mass_window_spin_dcaz_SS    = new TH2F("h_mass_window_spin_dcaz_SS","h_mass_window_spin_dcaz_SS", 20, -3.5, 16.5, 40, -2,2);
   h_mass_window_spin_dcar_SS    = new TH2F("h_mass_window_spin_dcar_SS","h_mass_window_spin_dcar_SS", 20, -3.5, 16.5, 40, 0,4);
   h_spin_phi       = new TH2F("h_spin_phi","h_spin_phi", 20, -3.5, 16.5, 160, -4,4);
   h_spin_phi_SS       = new TH2F("h_spin_phi_SS","h_spin_phi_SS", 20, -3.5, 16.5, 160, -4,4);
   h_spin_y       = new TH2F("h_spin_y","h_spin_y", 20, -3.5, 16.5, 60, -3,3);
   h_spin_y_R       = new TH2F("h_spin_y_R","h_spin_y_R", 20, -3.5, 16.5, 60, -3,3);
   h_spin_y_L       = new TH2F("h_spin_y_L","h_spin_y_L", 20, -3.5, 16.5, 60, -3,3);
   h_spin_y_SS       = new TH2F("h_spin_y_SS","h_spin_y_SS", 20, -3.5, 16.5, 60, -3,3);
   h_spin_y_SS_R       = new TH2F("h_spin_y_SS_R","h_spin_y_SS_R", 20, -3.5, 16.5, 60, -3,3);
   h_spin_y_SS_L       = new TH2F("h_spin_y_SS_L","h_spin_y_SS_L", 20, -3.5, 16.5, 60, -3,3);
   h_spin_pT       = new TH2F("h_spin_pT",     "h_spin_pT",      20, -3.5, 16.5, 30, 0,3);
   h_spin_pT_R     = new TH2F("h_spin_pT_R",   "h_spin_pT_R",    20, -3.5, 16.5, 30, 0,3);
   h_spin_pT_L     = new TH2F("h_spin_pT_L",   "h_spin_pT_L",    20, -3.5, 16.5, 30, 0,3);
   h_spin_pT_SS    = new TH2F("h_spin_pT_SS",  "h_spin_pT_SS",   20, -3.5, 16.5, 30, 0,3);
   h_spin_pT_SS_R  = new TH2F("h_spin_pT_SS_R","h_spin_pT_SS_R", 20, -3.5, 16.5, 30, 0,3);
   h_spin_pT_SS_L  = new TH2F("h_spin_pT_SS_L","h_spin_pT_SS_L", 20, -3.5, 16.5, 30, 0,3);

   h_dca_spin_pT       = new TH2F("h_dca_spin_pT",     "h_dca_spin_pT",      20, -3.5, 16.5, 30, 0,3);
   h_dca_spin_pT_R     = new TH2F("h_dca_spin_pT_R",   "h_dca_spin_pT_R",    20, -3.5, 16.5, 30, 0,3);
   h_dca_spin_pT_L     = new TH2F("h_dca_spin_pT_L",   "h_dca_spin_pT_L",    20, -3.5, 16.5, 30, 0,3);
   h_dca_spin_pT_SS    = new TH2F("h_dca_spin_pT_SS",  "h_dca_spin_pT_SS",   20, -3.5, 16.5, 30, 0,3);
   h_dca_spin_pT_SS_R  = new TH2F("h_dca_spin_pT_SS_R","h_dca_spin_pT_SS_R", 20, -3.5, 16.5, 30, 0,3);
   h_dca_spin_pT_SS_L  = new TH2F("h_dca_spin_pT_SS_L","h_dca_spin_pT_SS_L", 20, -3.5, 16.5, 30, 0,3);

   h_dca_spin_y         = new TH2F("h_dca_spin_y",     "h_dca_spin_y",      20, -3.5, 16.5, 60, -3,3);
   h_dca_spin_y_L       = new TH2F("h_dca_spin_y_L",   "h_dca_spin_y_L",    20, -3.5, 16.5, 60, -3,3);
   h_dca_spin_y_R       = new TH2F("h_dca_spin_y_R",   "h_dca_spin_y_R",    20, -3.5, 16.5, 60, -3,3);
   h_dca_spin_y_SS      = new TH2F("h_dca_spin_y_SS",  "h_dca_spin_y_SS",   20, -3.5, 16.5, 60, -3,3);
   h_dca_spin_y_SS_L    = new TH2F("h_dca_spin_y_SS_L","h_dca_spin_y_SS_L", 20, -3.5, 16.5, 60, -3,3);
   h_dca_spin_y_SS_R    = new TH2F("h_dca_spin_y_SS_R","h_dca_spin_y_SS_R", 20, -3.5, 16.5, 60, -3,3);

   h_dca_spin_sy_y         = new TH2F("h_dca_spin_sy_y",     "h_dca_spin_sy_y",      20, -3.5, 16.5, 60, -3,3);
   h_dca_spin_sy_y_L       = new TH2F("h_dca_spin_sy_y_L",   "h_dca_spin_sy_y_L",    20, -3.5, 16.5, 60, -3,3);
   h_dca_spin_sy_y_R       = new TH2F("h_dca_spin_sy_y_R",   "h_dca_spin_sy_y_R",    20, -3.5, 16.5, 60, -3,3);
   h_dca_spin_sy_y_SS      = new TH2F("h_dca_spin_sy_y_SS",  "h_dca_spin_sy_y_SS",   20, -3.5, 16.5, 60, -3,3);
   h_dca_spin_sy_y_SS_L    = new TH2F("h_dca_spin_sy_y_SS_L","h_dca_spin_sy_y_SS_L", 20, -3.5, 16.5, 60, -3,3);
   h_dca_spin_sy_y_SS_R    = new TH2F("h_dca_spin_sy_y_SS_R","h_dca_spin_sy_y_SS_R", 20, -3.5, 16.5, 60, -3,3);

   h_RP2E_dca_spin_y         = new TH2F("h_RP2E_dca_spin_y",     "h_RP2E_dca_spin_y",      20, -3.5, 16.5, 60, -3,3);
   h_RP2E_dca_spin_y_L       = new TH2F("h_RP2E_dca_spin_y_L",   "h_RP2E_dca_spin_y_L",    20, -3.5, 16.5, 60, -3,3);
   h_RP2E_dca_spin_y_R       = new TH2F("h_RP2E_dca_spin_y_R",   "h_RP2E_dca_spin_y_R",    20, -3.5, 16.5, 60, -3,3);
   h_RP2E_dca_spin_y_SS      = new TH2F("h_RP2E_dca_spin_y_SS",  "h_RP2E_dca_spin_y_SS",   20, -3.5, 16.5, 60, -3,3);
   h_RP2E_dca_spin_y_SS_L    = new TH2F("h_RP2E_dca_spin_y_SS_L","h_RP2E_dca_spin_y_SS_L", 20, -3.5, 16.5, 60, -3,3);
   h_RP2E_dca_spin_y_SS_R    = new TH2F("h_RP2E_dca_spin_y_SS_R","h_RP2E_dca_spin_y_SS_R", 20, -3.5, 16.5, 60, -3,3);

   h_dca_spin_phi       = new TH2F("h_dca_spin_phi",   "h_dca_spin_phi",    20, -3.5, 16.5, 160, -4,4);
   h_dca_spin_phi_SS    = new TH2F("h_dca_spin_phi_SS","h_dca_spin_phi_SS", 20, -3.5, 16.5, 160, -4,4);
   h_dca_spin_cosphi       = new TH2F("h_dca_spin_cosphi",   "h_dca_spin_cosphi",    20, -3.5, 16.5, 44, -1.1,1.1);
   h_dca_spin_cosphi_SS    = new TH2F("h_dca_spin_cosphi_SS","h_dca_spin_cosphi_SS", 20, -3.5, 16.5, 44, -1.1,1.1);
   h_spin_cosphi       = new TH2F("h_spin_cosphi",   "h_spin_cosphi",    20, -3.5, 16.5, 44, -1.1,1.1);
   h_spin_cosphi_SS    = new TH2F("h_spin_cosphi_SS","h_spin_cosphi_SS", 20, -3.5, 16.5, 44, -1.1,1.1);

   h_dca_cosphi_sinphi_u       = new TH2F("h_dca_cosphi_sinphi_u",   "h_dca_cosphi_sinphi_u", 44, -1.1,1.1, 44, -1.1,1.1);
   h_dca_cosphi_sinphi_d       = new TH2F("h_dca_cosphi_sinphi_d",   "h_dca_cosphi_sinphi_d", 44, -1.1,1.1, 44, -1.1,1.1);


   h_spin_w       = new TH2F("h_spin_w","h_spin_w", 20, -3.5, 16.5, 60, 0,180);
   h_spin_w_R       = new TH2F("h_spin_w_R","h_spin_w_R", 20, -3.5, 16.5, 60, 0,180);
   h_spin_w_L       = new TH2F("h_spin_w_L","h_spin_w_L", 20, -3.5, 16.5, 60, 0,180);
   h_spin_w_SS       = new TH2F("h_spin_w_SS","h_spin_w_SS", 20, -3.5, 16.5, 60, 0,180);
   h_spin_w_SS_R       = new TH2F("h_spin_w_SS_R","h_spin_w_SS_R", 20, -3.5, 16.5, 60, 0,180);
   h_spin_w_SS_L       = new TH2F("h_spin_w_SS_L","h_spin_w_SS_L", 20, -3.5, 16.5, 60, 0,180);


   h_zdceadc_sum       = new TH1F("h_zdceadc_sum","h_zdceadc_sum", 1000, -0.5, 999.5);
   h_zdcwadc_sum       = new TH1F("h_zdcwadc_sum","h_zdcwadc_sum", 1000, -0.5, 999.5);
   h_zdceadc_sum_mee       = new TH2F("h_zdceadc_sum_mee","h_zdceadc_sum_mee", 1000, -0.5, 999.5, 40, 2,6);
   h_zdcwadc_sum_mee       = new TH2F("h_zdcwadc_sum_mee","h_zdcwadc_sum_mee", 1000, -0.5, 999.5, 40, 2,6);
   h_zdceadc_sum_mee_SS       = new TH2F("h_zdceadc_sum_mee_SS","h_zdceadc_sum_mee_SS", 1000, -0.5, 999.5, 40, 2,6);
   h_zdcwadc_sum_mee_SS       = new TH2F("h_zdcwadc_sum_mee_SS","h_zdcwadc_sum_mee_SS", 1000, -0.5, 999.5, 40, 2,6);

   h_n_fdtrk       = new TH1F("h_n_fdtrk","h_n_fdtrk",             200, -4.5, 195.5);
   h_nh_fdtrk       = new TH1F("h_nh_fdtrk","h_nh_fdtrk",             200, -4.5, 195.5);
   h_nhdedx_fdtrk   = new TH1F("h_nhdedx_fdtrk","h_nhdedx_fdtrk",     200, -4.5, 195.5);
   h_emax_emccl     = new TH1F("h_emax_emccl","h_emax_emccl",         1000, -4.5, 95.5);
   h_pair_nh_fdtrk       = new TH1F("h_pair_nh_fdtrk"    ,"h_pair_nh_fdtrk",         200, -4.5, 195.5);
   h_pair_nhdedx_fdtrk   = new TH1F("h_pair_nhdedx_fdtrk","h_pair_nhdedx_fdtrk",     200, -4.5, 195.5);
   h_pair_emax_emccl     = new TH1F("h_pair_emax_emccl"  ,"h_pair_emax_emccl",      1000, -4.5, 95.5);

   h_pair_pt_emax_emccl          = new TH2F("h_pair_pt_emax_emccl"  ,"h_pair_pt_emax_emccl",  30, 0,3,    1000, -4.5, 95.5);
   h_pair_pt_emax_emccl_mass     = new TH2F("h_pair_pt_emax_emccl_mass"  ,"h_pair_pt_emax_emccl_mass",  30, 0,3,    1000, -4.5, 95.5);

   h_eta_vs_phi_fdtrack             = new TH2F("h_eta_vs_phi_fdtrack",             "h_eta_vs_phi_fdtrack"              , 200 , -2 , 2, 600, -3, 3);
   h_eta_vs_phi_fdtrack_pairs_bemc  = new TH2F("h_eta_vs_phi_fdtrack_pairs_bemc",  "h_eta_vs_phi_fdtrack_pairs_bemc"   , 200 , -2 , 2, 600, -3, 3);
   h_phi_track_1_vs_phi_track_2_bemc= new TH2F("h_phi_track_1_vs_phi_track_2_bemc","h_phi_track_1_vs_phi_track_2_bemc" , 600, -3, 3, 600, -3, 3);

   h_chiSquare_pipi_ee = new TH2F("h_chiSquare_pipi_ee","h_chiSquare_pipi_ee",1500,-10,140,1500,-10,140);
   h_chiSquare_kk_ee   = new TH2F("h_chiSquare_kk_ee","h_chiSquare_kk_ee",1500,-10,140,1500,-10,140);
   h_chiSquare_piK_ee  = new TH2F("h_chiSquare_piK_ee","h_chiSquare_piK_ee",1500,-10,140,1500,-10,140);
   h_chiSquare_kPi_ee  = new TH2F("h_chiSquare_kPi_ee","h_chiSquare_kPi_ee",1500,-10,140,1500,-10,140);
   h_chiSquare_kP_ee   = new TH2F("h_chiSquare_kP_ee","h_chiSquare_kP_ee",1500,-10,140,1500,-10,140);
   h_chiSquare_pK_ee   = new TH2F("h_chiSquare_pK_ee","h_chiSquare_pK_ee",1500,-10,140,1500,-10,140);
   h_chiSquare_pPi_ee  = new TH2F("h_chiSquare_pPi_ee","h_chiSquare_pPi_ee",1500,-10,140,1500,-10,140);
   h_chiSquare_piP_ee  = new TH2F("h_chiSquare_piP_ee","h_chiSquare_piP_ee",1500,-10,140,1500,-10,140);

   h_chiSquare_pipi_ee_SS = new TH2F("h_chiSquare_pipi_ee_SS","h_chiSquare_pipi_ee_SS",1500,-10,140,1500,-10,140);
   h_chiSquare_kk_ee_SS   = new TH2F("h_chiSquare_kk_ee_SS","h_chiSquare_kk_ee_SS",1500,-10,140,1500,-10,140);
   h_chiSquare_piK_ee_SS  = new TH2F("h_chiSquare_piK_ee_SS","h_chiSquare_piK_ee_SS",1500,-10,140,1500,-10,140);
   h_chiSquare_kPi_ee_SS  = new TH2F("h_chiSquare_kPi_ee_SS","h_chiSquare_kPi_ee_SS",1500,-10,140,1500,-10,140);
   h_chiSquare_kP_ee_SS   = new TH2F("h_chiSquare_kP_ee_SS","h_chiSquare_kP_ee_SS",1500,-10,140,1500,-10,140);
   h_chiSquare_pK_ee_SS   = new TH2F("h_chiSquare_pK_ee_SS","h_chiSquare_pK_ee_SS",1500,-10,140,1500,-10,140);
   h_chiSquare_pPi_ee_SS  = new TH2F("h_chiSquare_pPi_ee_SS","h_chiSquare_pPi_ee_SS",1500,-10,140,1500,-10,140);
   h_chiSquare_piP_ee_SS  = new TH2F("h_chiSquare_piP_ee_SS","h_chiSquare_piP_ee_SS",1500,-10,140,1500,-10,140);

   h_chiSquare_sel_pipi_ee = new TH2F("h_chiSquare_sel_pipi_ee","h_chiSquare_sel_pipi_ee",1500,-10,140,1500,-10,140);
   h_chiSquare_sel_kk_ee   = new TH2F("h_chiSquare_sel_kk_ee",  "h_chiSquare_sel_kk_ee",1500,-10,140,1500,-10,140);
   h_chiSquare_sel_piK_ee  = new TH2F("h_chiSquare_sel_piK_ee", "h_chiSquare_sel_piK_ee",1500,-10,140,1500,-10,140);
   h_chiSquare_sel_kPi_ee  = new TH2F("h_chiSquare_sel_kPi_ee", "h_chiSquare_sel_kPi_ee",1500,-10,140,1500,-10,140);
   h_chiSquare_sel_kP_ee   = new TH2F("h_chiSquare_sel_kP_ee",  "h_chiSquare_sel_kP_ee",1500,-10,140,1500,-10,140);
   h_chiSquare_sel_pK_ee   = new TH2F("h_chiSquare_sel_pK_ee",  "h_chiSquare_sel_pK_ee",1500,-10,140,1500,-10,140);
   h_chiSquare_sel_pPi_ee  = new TH2F("h_chiSquare_sel_pPi_ee", "h_chiSquare_sel_pPi_ee",1500,-10,140,1500,-10,140);
   h_chiSquare_sel_piP_ee  = new TH2F("h_chiSquare_sel_piP_ee", "h_chiSquare_sel_piP_ee",1500,-10,140,1500,-10,140);

   h_chiSquare_sel_pipi_ee_SS = new TH2F("h_chiSquare_sel_pipi_ee_SS","h_chiSquare_sel_pipi_ee_SS",1500,-10,140,1500,-10,140);
   h_chiSquare_sel_kk_ee_SS   = new TH2F("h_chiSquare_sel_kk_ee_SS",  "h_chiSquare_sel_kk_ee_SS",1500,-10,140,1500,-10,140);
   h_chiSquare_sel_piK_ee_SS  = new TH2F("h_chiSquare_sel_piK_ee_SS", "h_chiSquare_sel_piK_ee_SS",1500,-10,140,1500,-10,140);
   h_chiSquare_sel_kPi_ee_SS  = new TH2F("h_chiSquare_sel_kPi_ee_SS", "h_chiSquare_sel_kPi_ee_SS",1500,-10,140,1500,-10,140);
   h_chiSquare_sel_kP_ee_SS   = new TH2F("h_chiSquare_sel_kP_ee_SS",  "h_chiSquare_sel_kP_ee_SS",1500,-10,140,1500,-10,140);
   h_chiSquare_sel_pK_ee_SS   = new TH2F("h_chiSquare_sel_pK_ee_SS",  "h_chiSquare_sel_pK_ee_SS",1500,-10,140,1500,-10,140);
   h_chiSquare_sel_pPi_ee_SS  = new TH2F("h_chiSquare_sel_pPi_ee_SS", "h_chiSquare_sel_pPi_ee_SS",1500,-10,140,1500,-10,140);
   h_chiSquare_sel_piP_ee_SS  = new TH2F("h_chiSquare_sel_piP_ee_SS", "h_chiSquare_sel_piP_ee_SS",1500,-10,140,1500,-10,140);

   //h_mee_OS  = new TH1F( "h_mee_OS", "Pair OS mass for e", 1000, 0,100) ; 
   //h_mee_SS  = new TH1F( "h_mee_SS", "Pair SS mass for e", 1000, 0,100) ; 


   return kStOK ;

}
Int_t sDstAnaMaker::Make( )
{
   std::cout<< "sDstAnaMaker::Make" <<std::endl;


	//cut variables
	
	int nHits_fdtrk_cut = 15; //15 is default for BEMC-only optimized extraction

	int nHitsFit_bemc = 11; //11 is the default for BEMC-only optimized extraction
	//int nHitsFit_eemc = 8;  //6 is the default

	//double energyCut_bemc = 0.5; //0.5 is the default for BEMC-only optimized extraction
	
	
	//double adcCut_eemc = 5;     // 10 is the default
	//double adcPOST_eemc = 20;     // 10 is the default
	//double energyCut_eemc = 0.3; //default value I used in online code

	//int num_EEMC_towers_in_cluster_cut = 1; //cut to decide if a EEMC track is really a potential lepton candidate


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

   TChain *tree_chain = new TChain("T");
	while(getline(fileListStream, fileName)){        
    TString tmp = fileName;
    std::cout << "Input file: " << fileName << std::endl;
    tree_chain->Add(fileName.c_str());
   }   

   InitTree(tree_chain);   

   if (fChain == 0) return 0;

   Long64_t nentries = fChain->GetEntries();

   std::cout<< "nentries = " << nentries <<std::endl;
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {

      //Long64_t ientry = LoadTree(jentry);
      //if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      //tree->GetEvent(jentry);
      // if (Cut(ientry) < 0) continue;
      double zdceadc_sum = zdceadc[0]+zdceadc[1]+zdceadc[2];
      double zdcwadc_sum = zdcwadc[0]+zdcwadc[1]+zdcwadc[2];
      h_zdceadc_sum       ->Fill(zdceadc_sum);
      h_zdcwadc_sum       ->Fill(zdcwadc_sum);



      if((trig_JPsiHTTP[0]||trig_RP2E[0]||trig_2E[0]) ){
         h_evt->Fill(0);//J/psi trigger
         h_evt_ntoftrig->Fill(ntoftrig[0]);
      }
      if(trig_Zerobias[0]){
         h_evt->Fill(-1);//Zero bias trigger
         h_zb_evt_ntoftrig->Fill(ntoftrig[0]);
      }
      if (ntoftrig[0]>4) continue;
      int nvtx_rank=0;
      double running_rank = -100;         
      int high_rank_vx = -100;  
      if((trig_JPsiHTTP[0]||trig_RP2E[0]||trig_2E[0]) ){
         h_evt->Fill(1);//tof
       
         //if(trig_JPsiHTTP[0]||trig_RP2E[0]||trig_2E[0]){
            //h_evt->Fill(0);
            //std::cout<<" n_fdvtx[0] "<<n_fdvtx[0]<<std::endl;
            h_nvtx->Fill(n_fdvtx[0]);

            for(int vx=0;vx<n_fdvtx[0];vx++){
               //std::cout<<" z_fdvtx[vx] "<<z_fdvtx[vx]<<std::endl;
               h_evt_spin_z->Fill(run_sb[bid7[0]],z_fdvtx[vx]);
               if(rank_fdvtx[vx]>0){
                  h_evt_z_rank->Fill(z_fdvtx[vx]);
                  if(rank_fdvtx[vx]>running_rank){
                     running_rank = rank_fdvtx[vx];
                     high_rank_vx = vx;
                  }
                  nvtx_rank++;
               }
               h_evt_z->Fill(z_fdvtx[vx]);
               h_evt_x->Fill(x_fdvtx[vx]);
               h_evt_y->Fill(y_fdvtx[vx]);
               //std::cout<<" "<<vx<<": "<< z_fdvtx[vx]<<std::endl;
            }
            h_nvtx_rank->Fill(nvtx_rank);
            if(running_rank>0)h_evt_z_high_rank->Fill(z_fdvtx[high_rank_vx]);
   
         //}
         for(int si=0;si<120;si++){
            h_evt_spin->Fill(run_sb[si],si);
         }
      }
      //if((trig_JPsiHTTP[0]||trig_RP2E[0]||trig_2E[0]) )

      for(int trk=0;trk< n_fdtrk[0];trk++){
         if(trig_JPsiHTTP[0]||trig_RP2E[0]||trig_2E[0]){
            h_nh_fdtrk       ->Fill(   nh_fdtrk[trk]       );
            h_nhdedx_fdtrk   ->Fill(   nhdedx_fdtrk[trk]   );
            h_emax_emccl     ->Fill(   emax_emccl[trk]     );
            h_eta_vs_phi_fdtrack -> Fill(eta_fdtrk[trk], phi_fdtrk[trk]);
         }
      }

      
      if(n_fdpair[0]>0){

         //if((trig_JPsiHTTP[0]||trig_RP2E[0]||trig_2E[0]) )
         for(int p=0; p<n_fdpair[0]; p++){
            int pair_trk_idx_1 = ifdtrk1_fdpair[p];  //fd_trk index for pair_particle 1
				int pair_trk_idx_2 = ifdtrk2_fdpair[p];  //fd_trk index for pair_particle 2
            
            double mee = mee_fdpair[p];
            double pt_pair = pt_fdpair[p];
            double muu = mmumu_fdpair[p];
            double y = rapee_fdpair[p];
            double phi = phi_fdpair[p];
            double Ep = 100;
            double Mjpsi = 3.0969;
            double w = sqrt(2*Ep*Mjpsi*exp(y));

				//bool eemcONLY = false;
				bool bemcONLY = false;
				//bool bothONLY = false;
			
					

		      int bemc_wedge_1 = 0;
		      int bemc_wedge_2 = 0;
		      int bemc_wedge_b2b = 0;

		      //int eemc_wedge_1 = 0;
            //int eemc_wedge_2 = 0;
            //int eemc_wedge_b2b = 0;

		      //int both_wedge_1 = 0;
		      //int both_wedge_2 = 0;
		      //int both_wedge_b2b = 0;

				pair_trk_idx_1 = ifdtrk1_fdpair[p];  //fd_trk index for pair_particle 1
				pair_trk_idx_2 = ifdtrk2_fdpair[p];  //fd_trk index for pair_particle 2

				int bemc_ht_idx_trk_1 = iemccl_fdtrk[pair_trk_idx_1]; //high tower BEMC for fd_trk 1 -> if == -9999, there was no matched cluster
     			int bemc_ht_idx_trk_2 = iemccl_fdtrk[pair_trk_idx_2]; //
	
				int eemc_ht_idx_trk_1 = ieemc_cl_fdtrk[pair_trk_idx_1]; //matched EEMC tower -> if == -9999, no matched tower
            int eemc_ht_idx_trk_2 = ieemc_cl_fdtrk[pair_trk_idx_2]; //
				

            
                  		//double prsSum_eemc = 0.0;
      		//double prsSum_trk1_eemc = 0.0;
      		//double prsSum_trk2_eemc = 0.0;

				//prsSum_trk1_eemc = eemc_energy_prs1_fdtrk[pair_trk_idx_1]+eemc_energy_prs2_fdtrk[pair_trk_idx_1];
				//prsSum_trk2_eemc = eemc_energy_prs1_fdtrk[pair_trk_idx_2]+eemc_energy_prs2_fdtrk[pair_trk_idx_2];

	         //identical particles
	         double chiSquare_pipi = 0;
	         double chiSquare_pp   = 0;
	         double chiSquare_ee   = 0;
	         double chiSquare_kk   = 0;	

	         //non-identical particles
	         double chiSquare_piK  = 0;
	         double chiSquare_kPi  = 0;
	         double chiSquare_kP   = 0;
	         double chiSquare_pK   = 0;
	         double chiSquare_pPi  = 0;
	         double chiSquare_piP  = 0;

	         //double chiSquare_ePi  = 0;
	         //double chiSquare_piE  = 0;
	         //double chiSquare_eK  = 0;
            //double chiSquare_kE  = 0;
	         //double chiSquare_eP  = 0;
            //double chiSquare_pE  = 0;	

				chiSquare_pipi = sigpi_fdtrk[pair_trk_idx_1]*sigpi_fdtrk[pair_trk_idx_1] + sigpi_fdtrk[pair_trk_idx_2]*sigpi_fdtrk[pair_trk_idx_2];
				chiSquare_ee   = sigel_fdtrk[pair_trk_idx_1]*sigel_fdtrk[pair_trk_idx_1] + sigel_fdtrk[pair_trk_idx_2]*sigel_fdtrk[pair_trk_idx_2];
				chiSquare_kk   = sigk_fdtrk[pair_trk_idx_1]*sigk_fdtrk[pair_trk_idx_1] + sigk_fdtrk[pair_trk_idx_2]*sigk_fdtrk[pair_trk_idx_2];
				chiSquare_pp   = sigp_fdtrk[pair_trk_idx_1]*sigp_fdtrk[pair_trk_idx_1] + sigp_fdtrk[pair_trk_idx_2]*sigp_fdtrk[pair_trk_idx_2];	

				chiSquare_piK  = sigpi_fdtrk[pair_trk_idx_1]*sigpi_fdtrk[pair_trk_idx_1] + sigk_fdtrk[pair_trk_idx_2]*sigk_fdtrk[pair_trk_idx_2];
    			chiSquare_kPi  = sigk_fdtrk[pair_trk_idx_1]*sigk_fdtrk[pair_trk_idx_1] + sigpi_fdtrk[pair_trk_idx_2]*sigpi_fdtrk[pair_trk_idx_2];
    			chiSquare_kP   = sigk_fdtrk[pair_trk_idx_1]*sigk_fdtrk[pair_trk_idx_1] + sigp_fdtrk[pair_trk_idx_2]*sigp_fdtrk[pair_trk_idx_2];
    			chiSquare_pK   = sigp_fdtrk[pair_trk_idx_1]*sigp_fdtrk[pair_trk_idx_1] + sigk_fdtrk[pair_trk_idx_2]*sigk_fdtrk[pair_trk_idx_2];
    			chiSquare_pPi  = sigp_fdtrk[pair_trk_idx_1]*sigp_fdtrk[pair_trk_idx_1] + sigpi_fdtrk[pair_trk_idx_2]*sigpi_fdtrk[pair_trk_idx_2];
    			chiSquare_piP  = sigpi_fdtrk[pair_trk_idx_1]*sigpi_fdtrk[pair_trk_idx_1] + sigp_fdtrk[pair_trk_idx_2]*sigp_fdtrk[pair_trk_idx_2];

    			//chiSquare_ePi  = sigel_fdtrk[pair_trk_idx_1]*sigel_fdtrk[pair_trk_idx_1] + sigpi_fdtrk[pair_trk_idx_2]*sigpi_fdtrk[pair_trk_idx_2];
    			//chiSquare_piE  = sigpi_fdtrk[pair_trk_idx_1]*sigpi_fdtrk[pair_trk_idx_1] + sigel_fdtrk[pair_trk_idx_2]*sigel_fdtrk[pair_trk_idx_2];
				
    			//chiSquare_eK  = sigel_fdtrk[pair_trk_idx_1]*sigel_fdtrk[pair_trk_idx_1] + sigk_fdtrk[pair_trk_idx_2]*sigk_fdtrk[pair_trk_idx_2];
    			//chiSquare_kE  = sigk_fdtrk[pair_trk_idx_1]*sigk_fdtrk[pair_trk_idx_1] + sigel_fdtrk[pair_trk_idx_2]*sigel_fdtrk[pair_trk_idx_2];
				
    			//chiSquare_eP  = sigel_fdtrk[pair_trk_idx_1]*sigel_fdtrk[pair_trk_idx_1] + sigp_fdtrk[pair_trk_idx_2]*sigp_fdtrk[pair_trk_idx_2];
    			//chiSquare_pE  = sigp_fdtrk[pair_trk_idx_1]*sigp_fdtrk[pair_trk_idx_1] + sigel_fdtrk[pair_trk_idx_2]*sigel_fdtrk[pair_trk_idx_2];

            // 3/5 for blue up/down
            // run_sb[bid7[0]]

            //if((trig_JPsiHTTP[0]||trig_RP2E[0]||trig_2E[0]) ){
         	   if(n_fdtrk[0] <20 ){ // number of BEMC clusters
				      if((nh_fdtrk[pair_trk_idx_1] >= nHits_fdtrk_cut && nh_fdtrk[pair_trk_idx_2] >= nHits_fdtrk_cut)){ //hits per track
				         if((nhdedx_fdtrk[pair_trk_idx_1] >= nHitsFit_bemc && nhdedx_fdtrk[pair_trk_idx_2] >= nHitsFit_bemc)){ //hits used for de/dx calculation  
                        if(dcaz_fdtrk[pair_trk_idx_1]<1 && dcaz_fdtrk[pair_trk_idx_2]<1 && dcar_fdtrk[pair_trk_idx_1]<1.5 && dcar_fdtrk[pair_trk_idx_2]<1.5){
                           if(q_fdpair[p]==0){
                              h_spin_mee_trkCuts[0] ->Fill(run_sb[bid7[0]], mee);
            			   	   if(n_fdtrk[0] <15 ){h_spin_mee_trkCuts[1] ->Fill(run_sb[bid7[0]], mee); }
            			   	   if(n_fdtrk[0] <10 ){h_spin_mee_trkCuts[2] ->Fill(run_sb[bid7[0]], mee); }
            			   	   if(n_fdtrk[0] <5 ){ h_spin_mee_trkCuts[3] ->Fill(run_sb[bid7[0]], mee); }
                           }else{
                              h_spin_mee_SS_trkCuts[0] ->Fill(run_sb[bid7[0]], mee);
            			   	   if(n_fdtrk[0] <15 ){h_spin_mee_SS_trkCuts[1] ->Fill(run_sb[bid7[0]], mee); }
            			   	   if(n_fdtrk[0] <10 ){h_spin_mee_SS_trkCuts[2] ->Fill(run_sb[bid7[0]], mee); }
            			   	   if(n_fdtrk[0] <5 ){ h_spin_mee_SS_trkCuts[3] ->Fill(run_sb[bid7[0]], mee); }
                           }
                        }
                     }
                  }
               }
            //}
            //if(trig_RP2MU[0] ){//&& z_fdvtx[high_rank_vx]<=100){
				   if((nh_fdtrk[pair_trk_idx_1] >= nHits_fdtrk_cut && nh_fdtrk[pair_trk_idx_2] >= nHits_fdtrk_cut)){ //hits per track
				      //if((nhdedx_fdtrk[pair_trk_idx_1] >= nHitsFit_bemc && nhdedx_fdtrk[pair_trk_idx_2] >= nHitsFit_bemc)){ //hits used for de/dx calculation  
                  if (std::fabs(vtxz_fdtrk[pair_trk_idx_1])<100 && std::fabs(vtxz_fdtrk[pair_trk_idx_2])<100){
                     //if(chiSquare_kk<10 && chiSquare_ee>10 && chiSquare_pipi>10 && dcaz_fdtrk[pair_trk_idx_1]<1 && dcaz_fdtrk[pair_trk_idx_2]<1 && dcar_fdtrk[pair_trk_idx_1]<1.5 && dcar_fdtrk[pair_trk_idx_2]<1.5){
                     //if(chiSquare_ee>10 && dcaz_fdtrk[pair_trk_idx_1]<1 && dcaz_fdtrk[pair_trk_idx_2]<1 && dcar_fdtrk[pair_trk_idx_1]<1.5 && dcar_fdtrk[pair_trk_idx_2]<1.5){
                        h_spin_muu_cuts             ->Fill(run_sb[bid7[0]], muu); 
                        if(q_fdpair[p]==0){
                           h_muu_OS    ->Fill(muu);
                           if(chiSquare_kk<9 
                              && chiSquare_pipi  > 10.0  
                              && chiSquare_piK > 10.0 
										&& chiSquare_kPi > 10.0 
                              && chiSquare_kP  > 10.0
                              && chiSquare_pK  > 10.0  
										&& chiSquare_pPi > 10.0 
                              && chiSquare_piP > 10.0 
                              && trig_RPUPC[0])   h_mkk_OS    ->Fill(mkk_fdpair[p]);
                           if(chiSquare_pp<10 && trig_RPUPC[0])   h_mpp_OS    ->Fill(mpp_fdpair[p]);
                           if(chiSquare_pipi<9
                              && chiSquare_kk  > 10.0  
                              && chiSquare_piK > 10.0 
										&& chiSquare_kPi > 10.0 
                              && chiSquare_kP  > 10.0
                              && chiSquare_pK  > 10.0  
										&& chiSquare_pPi > 10.0 
                              && chiSquare_piP > 10.0 
                              && trig_RPUPC[0]) h_mpipi_OS  ->Fill(mpipi_fdpair[p]);
                           if(q_fdtrk[pair_trk_idx_1]<0){
                              if (chiSquare_piP<10)h_mpiP_OS->Fill(mpiP_fdpair[p]);
                           }else{
                              if(chiSquare_pPi<10)h_mpiP_OS->Fill(mPpi_fdpair[p]);
                           }
                        }else{
                           h_muu_SS    ->Fill(muu);
                           if(chiSquare_kk<9 
                              && chiSquare_pipi  > 10.0  
                              && chiSquare_piK > 10.0 
										&& chiSquare_kPi > 10.0 
                              && chiSquare_kP  > 10.0
                              && chiSquare_pK  > 10.0  
										&& chiSquare_pPi > 10.0 
                              && chiSquare_piP > 10.0 
                              && trig_RPUPC[0])   h_mkk_SS    ->Fill(mkk_fdpair[p]);
                           if(chiSquare_pp<10 && trig_RPUPC[0])   h_mpp_SS    ->Fill(mpp_fdpair[p]);
                           if(chiSquare_pipi<9
                              && chiSquare_kk  > 10.0  
                              && chiSquare_piK > 10.0 
										&& chiSquare_kPi > 10.0 
                              && chiSquare_kP  > 10.0
                              && chiSquare_pK  > 10.0  
										&& chiSquare_pPi > 10.0 
                              && chiSquare_piP > 10.0 
                              && trig_RPUPC[0]) h_mpipi_SS  ->Fill(mpipi_fdpair[p]);
                        
                           if(q_fdtrk[pair_trk_idx_1]<0){
                              if (chiSquare_piP<10)h_mpiP_SS->Fill(mpiP_fdpair[p]);
                           }else{
                              if(chiSquare_pPi<10)h_mpiP_SS->Fill(mPpi_fdpair[p]);
                           }
                        }
                     //}
                  }
                  //}
               }    
            //}  
            if((trig_JPsiHTTP[0]||trig_RP2E[0]||trig_2E[0]) ){
               h_evt->Fill(2);//|nfd pair exists
               h_pair_evt_z->Fill(vtxz_fdtrk[pair_trk_idx_1]);
               h_pair_evt_z->Fill(vtxz_fdtrk[pair_trk_idx_2]);
               
            }
            
            //if (vtxi_fdtrk[pair_trk_idx_1]!=vtxi_fdtrk[pair_trk_idx_2]) continue; 
            if (std::fabs(vtxz_fdtrk[pair_trk_idx_1])>100 || std::fabs(vtxz_fdtrk[pair_trk_idx_2])>100) continue;
            if((trig_JPsiHTTP[0]||trig_RP2E[0]||trig_2E[0]) )h_evt->Fill(3);//z vertex
            //if (std::fabs(z_fdvtx[high_rank_vx])>100) continue;

            if(bemc_ht_idx_trk_1 >= 0 && bemc_ht_idx_trk_2 >= 0 && eemc_ht_idx_trk_1 >= 0 && eemc_ht_idx_trk_2 >= 0 ){continue;}
            if(bemc_ht_idx_trk_1 >= 0 && bemc_ht_idx_trk_2 >= 0 && eemc_ht_idx_trk_1 >= 0){continue;}
            if(bemc_ht_idx_trk_1 >= 0 && bemc_ht_idx_trk_2 >= 0 && eemc_ht_idx_trk_2 >= 0){continue;}
            
            //if(n_fdvtx[0]>2)continue;

				bemc_wedge_1 = bemc_phi_wedge(phiext_fdtrk[pair_trk_idx_1]);
				bemc_wedge_2 = bemc_phi_wedge(phiext_fdtrk[pair_trk_idx_2]);

				bemc_wedge_b2b = TMath::Abs(bemc_wedge_2 - bemc_wedge_1);

				//eemc_wedge_1 = eemc_phi_wedge(eemc_phiext_fdtrk[pair_trk_idx_1]);
				//eemc_wedge_2 = eemc_phi_wedge(eemc_phiext_fdtrk[pair_trk_idx_2]);

				//eemc_wedge_b2b = TMath::Abs(eemc_wedge_2 - eemc_wedge_1);
				
				//if(bemc_ht_idx_trk_1 >= 0 && eemc_ht_idx_trk_2 >= 0){
				//	both_wedge_1 = both_phi_wedge(phiext_fdtrk[pair_trk_idx_1]);
            //   both_wedge_2 = both_phi_wedge(eemc_phiext_fdtrk[pair_trk_idx_2]);
				//}				
				//if(eemc_ht_idx_trk_1 >= 0 && bemc_ht_idx_trk_2 >= 0){
            //   both_wedge_1 = both_phi_wedge(eemc_phiext_fdtrk[pair_trk_idx_1]);
            //   both_wedge_2 = both_phi_wedge(phiext_fdtrk[pair_trk_idx_2]);
            //}

            //both_wedge_b2b = TMath::Abs(both_wedge_2 - both_wedge_1);

				


				

				if(bemc_ht_idx_trk_1 >= 0 && bemc_ht_idx_trk_2 >= 0 && eemc_ht_idx_trk_1 < 0 && eemc_ht_idx_trk_2 < 0 ){ bemcONLY = true;}
				//if(bemc_ht_idx_trk_1 < 0 && bemc_ht_idx_trk_2 < 0 && eemc_ht_idx_trk_1 >= 0 && eemc_ht_idx_trk_2 >= 0 ){ eemcONLY = true;}
				//if((bemc_ht_idx_trk_1 >= 0 && eemc_ht_idx_trk_2 >= 0) ||
				//   (eemc_ht_idx_trk_1 >= 0 && bemc_ht_idx_trk_2 >= 0)){bothONLY = true;}
            
            if((trig_JPsiHTTP[0]||trig_RP2E[0]||trig_2E[0]) && bemcONLY)h_evt->Fill(4);//trk is there and it is matched to BEMC
	            
				//--------------------------------------------------------
				//-----------------track quality cuts---------------------
				//--------------------------------------------------------
				
				
				//--------------------------------------------------------------------------------------------------
				//BEMC cuts application - THESE HAVE BEEN OPTIMIZED AS OF JAN 24th, 2023 - DO NOT CHANGE LIGHTLY!!!
				//--------------------------------------------------------------------------------------------------
            h_pair_bemc_wedge_b2b ->Fill(   bemc_wedge_b2b );
            h_pair_nh_fdtrk       ->Fill(   nh_fdtrk[pair_trk_idx_1]       );
            h_pair_nhdedx_fdtrk   ->Fill(   nhdedx_fdtrk[pair_trk_idx_1]   );
            h_pair_emax_emccl     ->Fill(   emax_emccl[bemc_ht_idx_trk_1]     );
            h_pair_nh_fdtrk       ->Fill(   nh_fdtrk[pair_trk_idx_2]       );
            h_pair_nhdedx_fdtrk   ->Fill(   nhdedx_fdtrk[pair_trk_idx_2]   );
            h_pair_emax_emccl     ->Fill(   emax_emccl[bemc_ht_idx_trk_2]     );

             
            if((trig_JPsiHTTP[0]||trig_RP2E[0]||trig_2E[0]) ){
               if( bemc_wedge_b2b == 3){
                  h_pair_pt_emax_emccl     ->Fill(pt_fdtrk[pair_trk_idx_1],emax_emccl[bemc_ht_idx_trk_1]);     
                  h_pair_pt_emax_emccl     ->Fill(pt_fdtrk[pair_trk_idx_2],emax_emccl[bemc_ht_idx_trk_2]);     


                  h_eta_vs_phi_fdtrack_pairs_bemc -> Fill(eta_fdtrk[pair_trk_idx_1], phi_fdtrk[pair_trk_idx_1]);
                  h_eta_vs_phi_fdtrack_pairs_bemc -> Fill(eta_fdtrk[pair_trk_idx_2], phi_fdtrk[pair_trk_idx_2]);
                  h_phi_track_1_vs_phi_track_2_bemc -> Fill(phi_fdtrk[pair_trk_idx_1], phi_fdtrk[pair_trk_idx_2]);
                  
				      if((nh_fdtrk[pair_trk_idx_1] >= nHits_fdtrk_cut && nh_fdtrk[pair_trk_idx_2] >= nHits_fdtrk_cut)){ //hits per track
				         if((nhdedx_fdtrk[pair_trk_idx_1] >= nHitsFit_bemc && nhdedx_fdtrk[pair_trk_idx_2] >= nHitsFit_bemc)){ //hits used for de/dx calculation  
                        
                        h_spin_muu_cuts             ->Fill(run_sb[bid7[0]], muu);                        
                        if(mee>2){



				               if(chiSquare_ee <= 10.0){

				                  if((chiSquare_pipi > 10.0 && chiSquare_kk > 10.0 && chiSquare_piK > 10.0 && 
				                  		chiSquare_kPi > 10.0  && chiSquare_kP  > 10.0 && chiSquare_pK > 10.0  &&
				                  		chiSquare_pPi > 10.0  && chiSquare_piP > 10.0) ) { 
                                       if(dcaz_fdtrk[pair_trk_idx_1]<1 && dcaz_fdtrk[pair_trk_idx_2]<1 && dcar_fdtrk[pair_trk_idx_1]<1.5 && dcar_fdtrk[pair_trk_idx_2]<1.5){
                                          h_spin_mee_cuts             ->Fill(run_sb[bid7[0]], mee);
                                       }
                              }
                           }
                        }
                     }
                  }
               }

            }

				if(bemcONLY && n_emccl[0] > 3){continue;} // number of BEMC clusters
				if(bemcONLY && bemc_wedge_b2b != 3) {continue;} //BEMC back-to-back cut
            if((trig_JPsiHTTP[0]||trig_RP2E[0]||trig_2E[0])  && bemcONLY)h_evt->Fill(5);//wedge            
				if(bemcONLY && (nh_fdtrk[pair_trk_idx_1] < nHits_fdtrk_cut || nh_fdtrk[pair_trk_idx_2] < nHits_fdtrk_cut)){continue;} //hits per track
				if(bemcONLY && (nhdedx_fdtrk[pair_trk_idx_1] < nHitsFit_bemc || nhdedx_fdtrk[pair_trk_idx_2] < nHitsFit_bemc)){continue;}   //hits used for de/dx calculation  
            if((trig_JPsiHTTP[0]||trig_RP2E[0]||trig_2E[0])  && bemcONLY)h_evt->Fill(6);//hits cut
            if((trig_JPsiHTTP[0]||trig_RP2E[0]||trig_2E[0])  && bemcONLY && mee>2)h_evt->Fill(7);//mee>2
            if(bemcONLY && (emax_emccl[bemc_ht_idx_trk_1] < 0.5 || emax_emccl[bemc_ht_idx_trk_2] < 0.5)){ continue;}
            if((trig_JPsiHTTP[0]||trig_RP2E[0]||trig_2E[0])  && bemcONLY  && mee>2)h_evt->Fill(8);//emax emccl ee > 0.5
            if((trig_JPsiHTTP[0]||trig_RP2E[0]||trig_2E[0]) && bemcONLY){
               if(q_fdpair[p]==0){
                  h_chiSquare_pipi_ee -> Fill(chiSquare_pipi, chiSquare_ee);
                  h_chiSquare_kk_ee   -> Fill(chiSquare_kk, chiSquare_ee);
                  h_chiSquare_piK_ee  -> Fill(chiSquare_piK, chiSquare_ee);
                  h_chiSquare_kPi_ee  -> Fill(chiSquare_kPi, chiSquare_ee);
                  h_chiSquare_kP_ee   -> Fill(chiSquare_kP, chiSquare_ee);
                  h_chiSquare_pK_ee   -> Fill(chiSquare_pK, chiSquare_ee);
                  h_chiSquare_pPi_ee  -> Fill(chiSquare_pPi, chiSquare_ee);
                  h_chiSquare_piP_ee  -> Fill(chiSquare_piP, chiSquare_ee);
               }else{
                  h_chiSquare_pipi_ee_SS -> Fill(chiSquare_pipi, chiSquare_ee);
                  h_chiSquare_kk_ee_SS   -> Fill(chiSquare_kk, chiSquare_ee);
                  h_chiSquare_piK_ee_SS  -> Fill(chiSquare_piK, chiSquare_ee);
                  h_chiSquare_kPi_ee_SS  -> Fill(chiSquare_kPi, chiSquare_ee);
                  h_chiSquare_kP_ee_SS   -> Fill(chiSquare_kP, chiSquare_ee);
                  h_chiSquare_pK_ee_SS   -> Fill(chiSquare_pK, chiSquare_ee);
                  h_chiSquare_pPi_ee_SS  -> Fill(chiSquare_pPi, chiSquare_ee);
                  h_chiSquare_piP_ee_SS  -> Fill(chiSquare_piP, chiSquare_ee);
               }
            }
            
            if(bemcONLY && chiSquare_ee > 10.0){continue;} //|| (chiSquare_ee < 10.0 && chiSquare_pipi < 10.0)) ) {continue;}
            if((trig_JPsiHTTP[0]||trig_RP2E[0]||trig_2E[0])  && bemcONLY  && mee>2)h_evt->Fill(9);//Chi2 ee
				
            if(bemcONLY && chiSquare_ee < 10.0 && (chiSquare_pipi < 10.0 || chiSquare_kk < 10.0 || chiSquare_piK < 10.0 || 
														   chiSquare_kPi < 10.0  || chiSquare_kP < 10.0 || chiSquare_pK < 10.0  ||
														   chiSquare_pPi < 10.0  || chiSquare_piP < 10.0) ) { continue; }
            
            if((trig_JPsiHTTP[0]||trig_RP2E[0]||trig_2E[0])  && bemcONLY  && mee>2)h_evt->Fill(10);//Chi2 other


				//--------------------------------------------------------------------------------------------------
				//BEMC + EEMC cuts application
				//--------------------------------------------------------------------------------------------------
            /*	
				//-------GOOD CUTS-------
				if(bothONLY && bemc_ht_idx_trk_1 >= 0 && eemc_ht_idx_trk_2 >= 0 && (nhdedx_fdtrk[pair_trk_idx_1] < nHitsFit_bemc || nhdedx_fdtrk[pair_trk_idx_2] < nHitsFit_eemc)){continue;}
				if(bothONLY && eemc_ht_idx_trk_1 >= 0 && bemc_ht_idx_trk_2 >= 0 && (nhdedx_fdtrk[pair_trk_idx_1] < nHitsFit_eemc || nhdedx_fdtrk[pair_trk_idx_2] < nHitsFit_bemc)){continue;}
				//---------------	-------

				//cout << "---both TPC nSigma cuts---" << endl;
				//cout << "TPC nSigma_e_1 = " << sigel_fdtrk[pair_trk_idx_1] << endl;
				//cout << "TPC nSigma_e_2 = " << sigel_fdtrk[pair_trk_idx_2] << endl;
			
				//topological momentum cut

				//double mom1 = p_fdtrk[pair_trk_idx_1];
				//double mom2 = p_fdtrk[pair_trk_idx_2];
	
				//if(bothONLY && TMath::Abs(mom1 - mom2) < 2.5){continue;}

				//nSigma cuts for pions and electrons in BEMC
				
				
				if(bothONLY && 
                ( (TMath::Abs(sigel_fdtrk[pair_trk_idx_1]) > 3.15 && bemc_ht_idx_trk_1 >= 0) || 
                  (TMath::Abs(sigel_fdtrk[pair_trk_idx_2]) > 3.15 && bemc_ht_idx_trk_2 >= 0) )) {continue;}
				
				if(bothONLY &&  
                ( (TMath::Abs(sigpi_fdtrk[pair_trk_idx_1]) < 3.15 && bemc_ht_idx_trk_1 >= 0) ||  
                  (TMath::Abs(sigpi_fdtrk[pair_trk_idx_2]) < 3.15 && bemc_ht_idx_trk_2 >= 0) )) {continue;}

				
				if(bothONLY &&  
                ( (TMath::Abs(sigk_fdtrk[pair_trk_idx_1]) < 3.15 && bemc_ht_idx_trk_1 >= 0) ||  
                  (TMath::Abs(sigk_fdtrk[pair_trk_idx_2]) < 3.15 && bemc_ht_idx_trk_2 >= 0) )) {continue;}	
				
				
				//if(bothONLY &&
                //( (TMath::Abs(sigp_fdtrk[pair_trk_idx_1]) < 4.0 && bemc_ht_idx_trk_1 >= 0) ||                                      
                //  (TMath::Abs(sigp_fdtrk[pair_trk_idx_2]) < 4.0 && bemc_ht_idx_trk_2 >= 0) )) {continue;}

				//if(bothONLY && chiSquare_ee < 10.0 && (chiSquare_pipi < 10.0 || chiSquare_kk < 10.0 || chiSquare_piK < 10.0 || chiSquare_kPi < 10.0 || 
				//									   chiSquare_kP < 10.0   || chiSquare_pK < 10.0 || chiSquare_pPi < 10.0 || chiSquare_piP < 10.0) ) { continue; }

				//energy cuts
				
				if(bothONLY && bemc_ht_idx_trk_1 >= 0 && emax_emccl[bemc_ht_idx_trk_1] < 0.5){continue;}
				if(bothONLY && bemc_ht_idx_trk_2 >= 0 && emax_emccl[bemc_ht_idx_trk_2] < 0.5){continue;}
		
				if(bothONLY && both_wedge_b2b != 3) {continue;}

				if(bothONLY && n_emccl[0] > 2){continue;}
				//if(bothONLY && numEEMCClusters > 2){ continue;} //2 is the nominal for current results
				
				//eemc_adc_fdtrk

				//if(bothONLY && bemc_ht_idx_trk_2 >= 0 && eemc_adc_fdtrk[pair_trk_idx_1] < adcCut_eemc){ continue; }
				//if(bothONLY && bemc_ht_idx_trk_1 >= 0 && eemc_adc_fdtrk[pair_trk_idx_2] < adcCut_eemc){ continue; }
				
				if(bothONLY && bemc_ht_idx_trk_2 >= 0 && eemc_num_towers_in_track[pair_trk_idx_1] > num_EEMC_towers_in_cluster_cut){ continue; }
				if(bothONLY && bemc_ht_idx_trk_1 >= 0 && eemc_num_towers_in_track[pair_trk_idx_2] > num_EEMC_towers_in_cluster_cut){ continue; }
				//if(bothONLY && tofMult > 2){ continue; }
				//-------------------------------------------------------------------------------------------------
				//EEMC cuts application
				//-------------------------------------------------------------------------------------------------

				//if(eemcONLY){ continue; }	
				//if(eemcONLY){ cout << "EEMC triggered event..." << endl;}	
				//if(eemcONLY && (eemc_adc_fdtrk[pair_trk_idx_1] < adcCut_eemc || eemc_adc_fdtrk[pair_trk_idx_2] < adcCut_eemc)){ continue; }

				//if(eemcONLY && (eemc_adc_post_fdtrk[pair_trk_idx_1] > 10 || eemc_adc_post_fdtrk[pair_trk_idx_2] > 10)){ continue; }
	
				//if(eemcONLY && (nhdedx_fdtrk[pair_trk_idx_1] < nHitsFit_eemc || nhdedx_fdtrk[pair_trk_idx_2] < nHitsFit_eemc)){continue;}
				//if(eemcONLY && eemc_wedge_b2b != 3) {continue;}				
				//if(eemcONLY && numEEMCClusters > 3) {continue;}
				//if(eemcONLY && chiSquare_ee > 10.0){ continue; }	
			
				//if(eemcONLY && ( (p_fdtrk[pair_trk_idx_1] > 1.0 && eemc_energy_fdtrk[pair_trk_idx_1] < 0.5) || 
                //                 (p_fdtrk[pair_trk_idx_2] > 1.0 && eemc_energy_fdtrk[pair_trk_idx_2] < 0.5) )){continue;}
				
				//if(eemcONLY && (eemc_energy_tow_fdtrk[pair_trk_idx_1]/p_fdtrk[pair_trk_idx_1] < 0.35 || eemc_energy_tow_fdtrk[pair_trk_idx_2]/p_fdtrk[pair_trk_idx_2] < 0.35)){ continue; }
				//if(eemcONLY && (eemc_energy_tow_fdtrk[pair_trk_idx_1] < 0.6 || eemc_energy_tow_fdtrk[pair_trk_idx_2] < 0.6)){ continue; }
				//if(eemcONLY && (prsSum_trk1_eemc == 0.0 || prsSum_trk2_eemc == 0.0)){continue;}
				//if(eemcONLY && eemc_energy_prs1_fdtrk[pair_trk_idx_1] < 0.001 && eemc_energy_prs2_fdtrk[pair_trk_idx_1] < 0.001){ continue; }
				//if(eemcONLY && eemc_energy_prs1_fdtrk[pair_trk_idx_2] < 0.001 && eemc_energy_prs2_fdtrk[pair_trk_idx_2] < 0.001){ continue; }
				//if(eemcONLY && (eemc_energy_prs1_fdtrk[pair_trk_idx_1] < 0.0001 || eemc_energy_prs1_fdtrk[pair_trk_idx_2] < 0.0001)){ continue; }
                //if(eemcONLY && (eemc_energy_prs2_fdtrk[pair_trk_idx_1] < 0.0001 || eemc_energy_prs2_fdtrk[pair_trk_idx_2] < 0.0001)){ continue; }
				//if(eemcONLY && (eemc_energy_post_fdtrk[pair_trk_idx_1] > 0.0 || eemc_energy_post_fdtrk[pair_trk_idx_2] > 0.0)){ continue; }
				//if(eemcONLY && totalEnergyEEMC < 2.0){continue;}
            


            if(eemcONLY && TMath::Abs(pt_fdtrk[pair_trk_idx_1] - pt_fdtrk[pair_trk_idx_2]) > 1.0) { continue; } 

            */

            if((trig_JPsiHTTP[0]||trig_RP2E[0]||trig_2E[0]) && bemcONLY){
               h_vtx1_vtx2->Fill(vtxi_fdtrk[pair_trk_idx_1],vtxi_fdtrk[pair_trk_idx_2]);
               if(q_fdpair[p]==0){
                  h_mee_OS->Fill(mee);
                  h_y_OS ->Fill(y);

                  h_spin_dcaz              ->Fill(run_sb[bid7[0]],dcaz_fdtrk[pair_trk_idx_1]);   
                  h_spin_dcaz              ->Fill(run_sb[bid7[0]],dcaz_fdtrk[pair_trk_idx_2]);   
                  h_spin_dcar              ->Fill(run_sb[bid7[0]],dcar_fdtrk[pair_trk_idx_1]);   
                  h_spin_dcar              ->Fill(run_sb[bid7[0]],dcar_fdtrk[pair_trk_idx_2]);   

                  h_zdceadc_sum_mee       ->Fill(zdceadc_sum, mee);
                  h_zdcwadc_sum_mee       ->Fill(zdcwadc_sum, mee);
                  h_spin_mee              ->Fill(run_sb[bid7[0]], mee);
                  h_n_fdtrk               ->Fill(   n_fdtrk[0]       );
                  if(mee > 2.8 && mee < 3.2){
                     h_pair_pt_emax_emccl_mass->Fill(pt_fdtrk[pair_trk_idx_1],emax_emccl[bemc_ht_idx_trk_1]);
                     h_pair_pt_emax_emccl_mass->Fill(pt_fdtrk[pair_trk_idx_2],emax_emccl[bemc_ht_idx_trk_2]);
                     
                     h_spin_y              ->Fill(run_sb[bid7[0]], y);
                     h_spin_pT              ->Fill(run_sb[bid7[0]], pt_pair);
                     h_spin_phi            ->Fill(run_sb[bid7[0]], phi);
                     h_spin_w              ->Fill(run_sb[bid7[0]], w);

                     
                     h_mass_window_spin_dcaz              ->Fill(run_sb[bid7[0]],dcaz_fdtrk[pair_trk_idx_1]);   
                     h_mass_window_spin_dcaz              ->Fill(run_sb[bid7[0]],dcaz_fdtrk[pair_trk_idx_2]);   
                     h_mass_window_spin_dcar              ->Fill(run_sb[bid7[0]],dcar_fdtrk[pair_trk_idx_1]);   
                     h_mass_window_spin_dcar              ->Fill(run_sb[bid7[0]],dcar_fdtrk[pair_trk_idx_2]);   

                     h_chiSquare_sel_pipi_ee -> Fill(chiSquare_pipi, chiSquare_ee);
                     h_chiSquare_sel_kk_ee   -> Fill(chiSquare_kk, chiSquare_ee);
                     h_chiSquare_sel_piK_ee  -> Fill(chiSquare_piK, chiSquare_ee);
                     h_chiSquare_sel_kPi_ee  -> Fill(chiSquare_kPi, chiSquare_ee);
                     h_chiSquare_sel_kP_ee   -> Fill(chiSquare_kP, chiSquare_ee);
                     h_chiSquare_sel_pK_ee   -> Fill(chiSquare_pK, chiSquare_ee);
                     h_chiSquare_sel_pPi_ee  -> Fill(chiSquare_pPi, chiSquare_ee);
                     h_chiSquare_sel_piP_ee  -> Fill(chiSquare_piP, chiSquare_ee);



                  }


                  if(dcaz_fdtrk[pair_trk_idx_1]<1 && dcaz_fdtrk[pair_trk_idx_2]<1 && dcar_fdtrk[pair_trk_idx_1]<1.5 && dcar_fdtrk[pair_trk_idx_2]<1.5){
                     if((trig_JPsiHTTP[0]||trig_RP2E[0]||trig_2E[0])  && mee>2)h_evt->Fill(11);//DCA

                     if(mee > 2.8 && mee < 3.2){
                        h_dca_spin_y              ->Fill(run_sb[bid7[0]], y);
                        h_dca_spin_pT              ->Fill(run_sb[bid7[0]], pt_pair);
                        h_dca_spin_sy_y           ->Fill(run_sy[bid7[0]], y);
                        if(trig_RP2E[0]) h_RP2E_dca_spin_y         ->Fill(run_sb[bid7[0]], y);
                        h_dca_spin_phi            ->Fill(run_sb[bid7[0]], phi);   
                        h_dca_spin_cosphi         ->Fill(run_sb[bid7[0]], cos(phi));   
                        if(run_sb[bid7[0]]==3){
                           h_dca_cosphi_sinphi_u         ->Fill(cos(phi), sin(phi));   
                        }
                        if(run_sb[bid7[0]]==5){
                           h_dca_cosphi_sinphi_d         ->Fill(cos(phi), sin(phi));   
                        }                        
                        for(int si=0;si<120;si++){
                           h_evt_spin_mw->Fill(run_sb[si],si);
                        }
                        if(cos(phi)>=0){
                           h_dca_spin_y_L              ->Fill(run_sb[bid7[0]], y); 
                           h_dca_spin_pT_L              ->Fill(run_sb[bid7[0]], pt_pair);
                           h_dca_spin_sy_y_L              ->Fill(run_sy[bid7[0]], y); 
                           if(trig_RP2E[0]) h_RP2E_dca_spin_y_L              ->Fill(run_sb[bid7[0]], y);
                        }else{
                           h_dca_spin_y_R              ->Fill(run_sb[bid7[0]], y); 
                           h_dca_spin_pT_R              ->Fill(run_sb[bid7[0]], pt_pair);
                           h_dca_spin_sy_y_R              ->Fill(run_sy[bid7[0]], y); 
                           if(trig_RP2E[0]) h_RP2E_dca_spin_y_R              ->Fill(run_sb[bid7[0]], y);
                        }
                     }
                     h_dca_spin_mee              ->Fill(run_sb[bid7[0]], mee);
                     if(trig_RP2E[0]) h_RP2E_dca_spin_mee              ->Fill(run_sb[bid7[0]], mee);
                     if(cos(phi)>=0){
                        h_dca_spin_mee_L              ->Fill(run_sb[bid7[0]], mee); 
                        if(trig_RP2E[0]) h_RP2E_dca_spin_mee_L              ->Fill(run_sb[bid7[0]], mee); 
                     }else{
                        h_dca_spin_mee_R              ->Fill(run_sb[bid7[0]], mee); 
                        if(trig_RP2E[0]) h_RP2E_dca_spin_mee_R              ->Fill(run_sb[bid7[0]], mee); 
                     }

                  }

                  if(cos(phi)>=0){
                     h_spin_mee_L              ->Fill(run_sb[bid7[0]], mee);
                     if(mee > 2.8 && mee < 3.2){
                        h_spin_y_L              ->Fill(run_sb[bid7[0]], y); 
                        h_spin_pT_L              ->Fill(run_sb[bid7[0]], pt_pair);
                        h_spin_w_L              ->Fill(run_sb[bid7[0]], w); 
                     }
                  }else{
                     h_spin_mee_R              ->Fill(run_sb[bid7[0]], mee);
                     if(mee > 2.8 && mee < 3.2){
                        h_spin_y_R              ->Fill(run_sb[bid7[0]], y); 
                        h_spin_pT_R              ->Fill(run_sb[bid7[0]], pt_pair);
                        h_spin_w_R              ->Fill(run_sb[bid7[0]], w); 
                     }
                  }
               }else{
                  h_mee_SS->Fill(mee);
                  h_y_SS ->Fill(y);
                  h_zdceadc_sum_mee_SS       ->Fill(zdceadc_sum, mee);
                  h_zdcwadc_sum_mee_SS       ->Fill(zdcwadc_sum, mee);
                  h_spin_mee_SS              ->Fill(run_sb[bid7[0]], mee);
                  h_spin_dcaz_SS              ->Fill(run_sb[bid7[0]],dcaz_fdtrk[pair_trk_idx_1]);   
                  h_spin_dcaz_SS              ->Fill(run_sb[bid7[0]],dcaz_fdtrk[pair_trk_idx_2]);                     
                  h_spin_dcar_SS              ->Fill(run_sb[bid7[0]],dcar_fdtrk[pair_trk_idx_1]);   
                  h_spin_dcar_SS              ->Fill(run_sb[bid7[0]],dcar_fdtrk[pair_trk_idx_2]);                     
                  if(mee > 2.8 && mee < 3.2){
                     h_spin_phi_SS              ->Fill(run_sb[bid7[0]], phi);
                     h_spin_y_SS                ->Fill(run_sb[bid7[0]], y);
                     h_spin_pT_SS              ->Fill(run_sb[bid7[0]], pt_pair);
                     h_spin_w_SS                ->Fill(run_sb[bid7[0]], w);
                     h_spin_cosphi_SS         ->Fill(run_sb[bid7[0]], cos(phi));   

                     h_mass_window_spin_dcaz_SS              ->Fill(run_sb[bid7[0]],dcaz_fdtrk[pair_trk_idx_1]);   
                     h_mass_window_spin_dcaz_SS              ->Fill(run_sb[bid7[0]],dcaz_fdtrk[pair_trk_idx_2]);                     
                     h_mass_window_spin_dcar_SS              ->Fill(run_sb[bid7[0]],dcar_fdtrk[pair_trk_idx_1]);   
                     h_mass_window_spin_dcar_SS              ->Fill(run_sb[bid7[0]],dcar_fdtrk[pair_trk_idx_2]);                     
    
                     h_chiSquare_sel_pipi_ee_SS -> Fill(chiSquare_pipi, chiSquare_ee);
                     h_chiSquare_sel_kk_ee_SS   -> Fill(chiSquare_kk, chiSquare_ee);
                     h_chiSquare_sel_piK_ee_SS  -> Fill(chiSquare_piK, chiSquare_ee);
                     h_chiSquare_sel_kPi_ee_SS  -> Fill(chiSquare_kPi, chiSquare_ee);
                     h_chiSquare_sel_kP_ee_SS   -> Fill(chiSquare_kP, chiSquare_ee);
                     h_chiSquare_sel_pK_ee_SS   -> Fill(chiSquare_pK, chiSquare_ee);
                     h_chiSquare_sel_pPi_ee_SS  -> Fill(chiSquare_pPi, chiSquare_ee);
                     h_chiSquare_sel_piP_ee_SS  -> Fill(chiSquare_piP, chiSquare_ee);
                  }
                  if(cos(phi)>=0){
                     h_spin_mee_SS_L              ->Fill(run_sb[bid7[0]], mee);
                     if(mee > 2.8 && mee < 3.2){
                        h_spin_pT_SS_L              ->Fill(run_sb[bid7[0]], pt_pair);
                        h_spin_y_SS_L              ->Fill(run_sb[bid7[0]], y); 
                        h_spin_w_SS_L              ->Fill(run_sb[bid7[0]], w); 
                     }
                  }else{
                     h_spin_mee_SS_R              ->Fill(run_sb[bid7[0]], mee);
                     if(mee > 2.8 && mee < 3.2){
                        h_spin_y_SS_R              ->Fill(run_sb[bid7[0]], y); 
                        h_spin_pT_SS_R              ->Fill(run_sb[bid7[0]], pt_pair);
                        h_spin_w_SS_R              ->Fill(run_sb[bid7[0]], w); 
                     }
                  }
                  if(dcaz_fdtrk[pair_trk_idx_1]<1 && dcaz_fdtrk[pair_trk_idx_2]<1 && dcar_fdtrk[pair_trk_idx_1]<1.5 && dcar_fdtrk[pair_trk_idx_2]<1.5){
                     if(mee > 2.8 && mee < 3.2){
                        h_dca_spin_y_SS              ->Fill(run_sb[bid7[0]], y);
                        h_dca_spin_pT_SS              ->Fill(run_sb[bid7[0]], pt_pair);
                        h_dca_spin_sy_y_SS              ->Fill(run_sy[bid7[0]], y);
                        if(trig_RP2E[0]) h_RP2E_dca_spin_y_SS         ->Fill(run_sb[bid7[0]], y);
                        h_dca_spin_phi_SS            ->Fill(run_sb[bid7[0]], phi);   
                        h_dca_spin_cosphi_SS         ->Fill(run_sb[bid7[0]], cos(phi));   

                        if(cos(phi)>=0){
                           h_dca_spin_y_SS_L              ->Fill(run_sb[bid7[0]], y); 
                           h_dca_spin_pT_SS_L              ->Fill(run_sb[bid7[0]], pt_pair);
                           h_dca_spin_sy_y_SS_L              ->Fill(run_sy[bid7[0]], y); 
                           if(trig_RP2E[0]) h_RP2E_dca_spin_y_SS_L              ->Fill(run_sb[bid7[0]], y);
                        }else{
                           h_dca_spin_y_SS_R              ->Fill(run_sb[bid7[0]], y); 
                           h_dca_spin_pT_SS_R              ->Fill(run_sb[bid7[0]], pt_pair);
                           h_dca_spin_sy_y_SS_R              ->Fill(run_sy[bid7[0]], y); 
                           if(trig_RP2E[0]) h_RP2E_dca_spin_y_SS_R              ->Fill(run_sb[bid7[0]], y);
                        }
                     }
                     h_dca_spin_mee_SS              ->Fill(run_sb[bid7[0]], mee);
                     if(trig_RP2E[0]) h_RP2E_dca_spin_mee_SS              ->Fill(run_sb[bid7[0]], mee);
                     if(cos(phi)>=0){
                        h_dca_spin_mee_SS_L              ->Fill(run_sb[bid7[0]], mee); 
                        if(trig_RP2E[0]) h_RP2E_dca_spin_mee_SS_L              ->Fill(run_sb[bid7[0]], mee); 
                     }else{
                        h_dca_spin_mee_SS_R              ->Fill(run_sb[bid7[0]], mee); 
                        if(trig_RP2E[0]) h_RP2E_dca_spin_mee_SS_R              ->Fill(run_sb[bid7[0]], mee); 
                     }


                  }                  
                  
               }
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

   h_evt_spin    ->Write();
   h_evt_spin_mw  ->Write();
   
   h_evt_spin_z ->Write();
   h_vtx1_vtx2 ->Write();
   h_evt   ->Write();
   h_nvtx  ->Write(); 
   h_nvtx_rank  ->Write(); 
   h_evt_z ->Write();
   h_evt_z_rank ->Write();
   h_evt_z_high_rank ->Write();
   h_pair_evt_z ->Write();
   h_evt_x ->Write();
   h_evt_y ->Write();
   h_evt_ntoftrig ->Write();
   h_zb_evt_ntoftrig ->Write();

   h_muu_SS -> Write();
   h_muu_OS -> Write();
   h_mpiP_OS -> Write();
   h_mpiP_SS -> Write();
   h_mkk_OS   -> Write();
   h_mkk_SS   -> Write();
   h_mpp_OS   -> Write();
   h_mpp_SS   -> Write();
   h_mpipi_OS -> Write();
   h_mpipi_SS -> Write();
   h_mee_SS -> Write();
   h_mee_OS -> Write();
   h_zdceadc_sum -> Write();
   h_zdcwadc_sum -> Write();

   h_spin_mee -> Write();
   h_spin_mee_cuts -> Write();
   h_spin_muu_cuts -> Write();
   h_spin_mee_R -> Write();
   h_spin_mee_L -> Write();
   h_spin_mee_SS -> Write();
   h_spin_mee_SS_R -> Write();
   h_spin_mee_SS_L -> Write();
   for(int i=0;i<4;i++){
      h_spin_mee_trkCuts[i]    -> Write();
      h_spin_mee_SS_trkCuts[i] -> Write();
   }

   h_spin_dcaz    -> Write();
   h_spin_dcar    -> Write();
   h_spin_dcaz_SS -> Write();
   h_spin_dcar_SS -> Write();
   h_mass_window_spin_dcaz    -> Write();
   h_mass_window_spin_dcar    -> Write();
   h_mass_window_spin_dcaz_SS -> Write();
   h_mass_window_spin_dcar_SS -> Write();
   
   h_dca_spin_mee       -> Write();
   h_dca_spin_mee_R     -> Write();
   h_dca_spin_mee_L     -> Write();
   h_dca_spin_mee_SS    -> Write();
   h_dca_spin_mee_SS_R  -> Write();
   h_dca_spin_mee_SS_L  -> Write();

   h_RP2E_dca_spin_mee       -> Write();
   h_RP2E_dca_spin_mee_R     -> Write();
   h_RP2E_dca_spin_mee_L     -> Write();
   h_RP2E_dca_spin_mee_SS    -> Write();
   h_RP2E_dca_spin_mee_SS_R  -> Write();
   h_RP2E_dca_spin_mee_SS_L  -> Write();

   h_dca_spin_phi -> Write();
   h_dca_spin_phi_SS -> Write();

   h_dca_spin_y -> Write();
   h_dca_spin_y_L -> Write();
   h_dca_spin_y_R -> Write();
   h_dca_spin_y_SS -> Write();
   h_dca_spin_y_SS_L -> Write();
   h_dca_spin_y_SS_R -> Write();

   h_dca_spin_pT -> Write();
   h_dca_spin_pT_L -> Write();
   h_dca_spin_pT_R -> Write();
   h_dca_spin_pT_SS -> Write();
   h_dca_spin_pT_SS_L -> Write();
   h_dca_spin_pT_SS_R -> Write();

   h_RP2E_dca_spin_y -> Write();
   h_RP2E_dca_spin_y_L -> Write();
   h_RP2E_dca_spin_y_R -> Write();
   h_RP2E_dca_spin_y_SS -> Write();
   h_RP2E_dca_spin_y_SS_L -> Write();
   h_RP2E_dca_spin_y_SS_R -> Write();

   h_dca_spin_sy_y -> Write();
   h_dca_spin_sy_y_L -> Write();
   h_dca_spin_sy_y_R -> Write();
   h_dca_spin_sy_y_SS -> Write();
   h_dca_spin_sy_y_SS_L -> Write();
   h_dca_spin_sy_y_SS_R -> Write();

   h_spin_cosphi -> Write();
   h_spin_cosphi_SS -> Write();
   h_dca_spin_cosphi -> Write();
   h_dca_spin_cosphi_SS -> Write();

   h_dca_cosphi_sinphi_u -> Write();
   h_dca_cosphi_sinphi_d -> Write();

   h_spin_phi -> Write();
   h_spin_phi_SS -> Write();

   h_pair_pt_emax_emccl -> Write();
   h_pair_pt_emax_emccl_mass -> Write();

   h_spin_y -> Write();
   h_spin_y_R -> Write();
   h_spin_y_L -> Write();
   h_spin_y_SS -> Write();
   h_spin_y_SS_R -> Write();
   h_spin_y_SS_L -> Write();
   h_spin_pT -> Write();
   h_spin_pT_R -> Write();
   h_spin_pT_L -> Write();
   h_spin_pT_SS -> Write();
   h_spin_pT_SS_R -> Write();
   h_spin_pT_SS_L -> Write();

   h_spin_w -> Write();
   h_spin_w_R -> Write();
   h_spin_w_L -> Write();
   h_spin_w_SS -> Write();
   h_spin_w_SS_R -> Write();
   h_spin_w_SS_L -> Write();

   h_zdceadc_sum_mee -> Write();
   h_zdcwadc_sum_mee -> Write();
   h_zdceadc_sum_mee_SS -> Write();
   h_zdcwadc_sum_mee_SS -> Write();

   h_y_OS -> Write();
   h_y_SS -> Write();
   h_pair_bemc_wedge_b2b -> Write();
   h_pair_nh_fdtrk       -> Write();
   h_pair_nhdedx_fdtrk   -> Write();
   h_pair_emax_emccl     -> Write();
   h_nh_fdtrk       -> Write();
   h_n_fdtrk       -> Write();
   h_nhdedx_fdtrk   -> Write();
   h_emax_emccl     -> Write();

   h_eta_vs_phi_fdtrack ->Write();
   h_eta_vs_phi_fdtrack_pairs_bemc ->Write();
   h_phi_track_1_vs_phi_track_2_bemc->Write();

   h_chiSquare_pipi_ee -> Write();
   h_chiSquare_kk_ee   -> Write();
   h_chiSquare_piK_ee  -> Write();
   h_chiSquare_kPi_ee  -> Write();
   h_chiSquare_kP_ee   -> Write();
   h_chiSquare_pK_ee   -> Write();
   h_chiSquare_pPi_ee  -> Write();
   h_chiSquare_piP_ee  -> Write();

   h_chiSquare_pipi_ee_SS -> Write();
   h_chiSquare_kk_ee_SS   -> Write();
   h_chiSquare_piK_ee_SS  -> Write();
   h_chiSquare_kPi_ee_SS  -> Write();
   h_chiSquare_kP_ee_SS   -> Write();
   h_chiSquare_pK_ee_SS   -> Write();
   h_chiSquare_pPi_ee_SS  -> Write();
   h_chiSquare_piP_ee_SS  -> Write();

   h_chiSquare_sel_pipi_ee -> Write();
   h_chiSquare_sel_kk_ee   -> Write();
   h_chiSquare_sel_piK_ee  -> Write();
   h_chiSquare_sel_kPi_ee  -> Write();
   h_chiSquare_sel_kP_ee   -> Write();
   h_chiSquare_sel_pK_ee   -> Write();
   h_chiSquare_sel_pPi_ee  -> Write();
   h_chiSquare_sel_piP_ee  -> Write();

   h_chiSquare_sel_pipi_ee_SS -> Write();
   h_chiSquare_sel_kk_ee_SS   -> Write();
   h_chiSquare_sel_piK_ee_SS  -> Write();
   h_chiSquare_sel_kPi_ee_SS  -> Write();
   h_chiSquare_sel_kP_ee_SS   -> Write();
   h_chiSquare_sel_pK_ee_SS   -> Write();
   h_chiSquare_sel_pPi_ee_SS  -> Write();
   h_chiSquare_sel_piP_ee_SS  -> Write();

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


int bemc_phi_wedge(double phi_ext){
	
	double pi = TMath::Pi();

	if(-pi <= phi_ext && phi_ext < -(2./3.)*pi) {return 1;}
	else if(-(2./3.)*pi <= phi_ext && phi_ext < -(1./3.)*pi) {return 2;}
	else if(-(1./3.)*pi <= phi_ext && phi_ext < 0.0) {return 3;}
	else if(0.0 <= phi_ext && phi_ext < pi/3.) {return 4;}
	else if(pi/3. <= phi_ext && phi_ext < (2./3.)*pi) {return 5;}
	else if((2./3.)*pi <= phi_ext && phi_ext < pi) {return 6;}

	return 0;

}



int eemc_phi_wedge(double phi_ext){

    double pi = TMath::Pi();

    if(-pi <= phi_ext && phi_ext < -(2./3.)*pi) {return 1;}
    else if(-(2./3.)*pi <= phi_ext && phi_ext < -(1./3.)*pi) {return 2;}
    else if(-(1./3.)*pi <= phi_ext && phi_ext < 0.0) {return 3;}
    else if(0.0 <= phi_ext && phi_ext < pi/3.) {return 4;}
    else if(pi/3. <= phi_ext && phi_ext < (2./3.)*pi) {return 5;}
    else if((2./3.)*pi <= phi_ext && phi_ext < pi) {return 6;}

    return 0;

}
     

int both_phi_wedge(double phi_ext){

    double pi = TMath::Pi();

    if(-pi <= phi_ext && phi_ext < -(2./3.)*pi) {return 1;}
    else if(-(2./3.)*pi <= phi_ext && phi_ext < -(1./3.)*pi) {return 2;}
    else if(-(1./3.)*pi <= phi_ext && phi_ext < 0.0) {return 3;}
    else if(0.0 <= phi_ext && phi_ext < pi/3.) {return 4;}
    else if(pi/3. <= phi_ext && phi_ext < (2./3.)*pi) {return 5;}
    else if((2./3.)*pi <= phi_ext && phi_ext < pi) {return 6;}

    return 0;

}

