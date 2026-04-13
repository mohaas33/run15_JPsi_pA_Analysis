

//#include "AtlasStyle.h"
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TSystem.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <vector>

#include "TFractionFitter.h"

const int p_bin = 3; // spin +1
const int n_bin = 5; // spin -1
const float mee_min = 2;
const float mee_max = 4.5;
float y_min = -1.2;
float y_max = 1.2;
const double Mjpsi = 3.0969;

double P_avg = 0.0;
double P_avg_err = 0.0;
double sqrt_cosPhiLR = 0.0;
double correction = 0.0;
double y_Mean = 0.154;



void draw_Mark(TString Mark, float x1, float y1, float x2, float y2, float transp=0){
    TLegend *leg_mark = new TLegend(x1, y1, x2, y2);
    leg_mark->SetTextSize(0.048);
    leg_mark-> SetTextFont(42);
    if(transp>0){
        leg_mark->SetFillColorAlpha(kWhite, transp);
     }else{
        leg_mark->SetFillColor(0);
        leg_mark->SetFillStyle(0);
    }
    leg_mark->SetBorderSize(0);

    leg_mark->AddEntry((TObject*)0, Mark, ""); // Independent text
    leg_mark->Draw();

}

double crystalball(const double *x, const double *p) {
  // if ((!x) || (!p)) return 0.; // just a precaution
  // [Constant] * ROOT::Math::crystalball_function(x, [Alpha], [N], [Sigma], [Mean])
  return (p[0] * crystalball_function(x[0], p[3], p[4], p[2], p[1]));
}

// Define the Crystal Ball function
Double_t crystalBall(Double_t *x, Double_t *par) {
    Double_t alpha = par[0];  // Tail parameter
    Double_t n     = par[1];  // Tail power
    Double_t mean  = par[2];  // Mean
    Double_t sigma = par[3];  // Width
    Double_t norm  = par[4];  // Normalization

    Double_t t = (x[0] - mean) / sigma;
    if (alpha < 0) t = -t;

    Double_t absAlpha = fabs(alpha);
    Double_t A = pow(n / absAlpha, n) * exp(-0.5 * absAlpha * absAlpha);
    Double_t B = n / absAlpha - absAlpha;

    //if (t > -absAlpha)
    //    return norm * exp(-0.5 * t * t);  // Gaussian core
    //else
    //    return norm * A * pow(B - t, -n);  // Power-law tail
    // Normalization factor for unit area
    Double_t sqrtPiOver2 = sqrt(TMath::Pi()/2);
    Double_t C = (n/absAlpha) * (1.0/(n-1)) * exp(-0.5*absAlpha*absAlpha);
    Double_t D = sqrtPiOver2 * (1 + TMath::Erf(absAlpha / sqrt(2.0)));

    Double_t N = 1.0 / (sigma * (C + D)); // total area = 1

    // Shape
    Double_t shape = (t > -absAlpha) ?
                     exp(-0.5 * t * t) :
                     A * pow(B - t, -n);

    return norm * N * shape;
}
// Normalized Crystal Ball: par[0]=alpha, par[1]=n, par[2]=mean, par[3]=sigma, par[4]=yield (area)
Double_t crystalBallNorm(Double_t *x, Double_t *par) {
    Double_t alpha = par[0];
    Double_t n     = par[1];
    Double_t mean  = par[2];
    Double_t sigma = par[3];
    Double_t yield = par[4];

    Double_t t = (x[0] - mean) / sigma;
    // allow negative alpha convention if you wish:
    if (alpha < 0) t = -t;

    Double_t absA = fabs(alpha);
    Double_t A = TMath::Power(n/absA, n) * TMath::Exp(-0.5 * absA * absA);
    Double_t B = n/absA - absA;

    // Integral over t (dimensionless) of shape:
    Double_t sqrtPiOver2 = TMath::Sqrt(TMath::Pi()/2.0);
    Double_t D = sqrtPiOver2 * (1.0 + TMath::Erf(absA / TMath::Sqrt(2.0)));
    Double_t C = (n/absA) * TMath::Exp(-0.5 * absA * absA) / (n - 1.0);

    Double_t I_t = D + C;           // integral over t of the shape
    Double_t N = 1.0 / (sigma * I_t); // factor so that Integral_x [ N * shape(t) ] = 1.0

    Double_t shape;
    if (t > -absA) {
        shape = TMath::Exp(-0.5 * t * t);
    } else {
        shape = A * TMath::Power(B - t, -n);
    }

    // Return yield * normalized shape (so integrated area = yield)
    return yield * N * shape;
}

// Full fit function: CB (signal) + exponential (background)
Double_t totalFit(Double_t *x, Double_t *par) {
    // CB params: par[0] - par[4]
    Double_t cb = crystalBallNorm(x, par);

    // Exponential background: par[5], par[6]
    Double_t bkg = par[5] * exp(par[6] * x[0]);

    return cb + bkg;
}
Double_t totalFitPol(Double_t *x, Double_t *par) {
    // CB params: par[0] - par[4]
    Double_t cb = crystalBallNorm(x, par);

    // Polinomial background: par[5], par[6]
    Double_t bkg  = par[5] + x[0] / par[6] + x[0]*x[0] / (par[7]*par[7]);

    return cb + bkg;
}


Double_t GaussExp(Double_t *x, Double_t *par) {
    // Parameters
    // par[0] = Gaussian amplitude
    // par[1] = Gaussian mean
    // par[2] = Gaussian sigma
    // par[3] = Exponential amplitude
    // par[4] = Exponential decay constant

    Double_t gauss = par[0] * TMath::Gaus(x[0], par[1], par[2], kTRUE);
    Double_t expo  = par[3] * TMath::Exp(-x[0] / par[4]);
    return gauss + expo;
}
Double_t GaussPol(Double_t *x, Double_t *par) {
    // Parameters
    // par[0] = Gaussian amplitude
    // par[1] = Gaussian mean
    // par[2] = Gaussian sigma
    // par[3] = Exponential amplitude
    // par[4] = Exponential decay constant

    Double_t gauss = par[0] * TMath::Gaus(x[0], par[1], par[2], kTRUE);
    Double_t pol  = par[3] + x[0] / par[4] + x[0]*x[0] / (par[5]*par[5]);
    return gauss + pol;
}
void fit_Jpsi_with_templates(TString canvasName,TH1F *hData, TH1F *hMC1, TH1F *hMC2, TH1F *hMC3, TString titleX, TString titleY, TString Mark="empty", TString Mark2="empty" ) {
    // Load histograms (could also be created in the same script)
    //TFile *f = TFile::Open("histos.root");
    //TH1D *hData = (TH1D*)f->Get("hData");     // data histogram
    //TH1D *hMC1  = (TH1D*)f->Get("hProc1");    // MC: process 1
    //TH1D *hMC2  = (TH1D*)f->Get("hProc2");    // MC: process 2
    //TH1D *hMC3  = (TH1D*)f->Get("hProc3");    // MC: process 3

    // Make sure binning matches exactly
    //if (!hData || !hMC1 || !hMC2 || !hMC3) {
    //    Error("fit_Jpsi_with_templates", "Missing histograms!");
    //    return;
    //}
    TCanvas *c_MC = new TCanvas("c_MC", "MC", 800, 600);

    hMC1  ->SetLineColor(2);
    hMC2  ->SetLineColor(4);
    hMC3  ->SetLineColor(6);

    hMC1->Draw("HIST");
    hMC2->Draw("HISTsame");
    hMC3->Draw("HISTsame");

    c_MC -> Print("./Plots/"+canvasName+"_templates.png");

    float n_hMC1 = hMC1->Integral();
    float n_hMC2 = hMC2->Integral();
    float n_hMC3 = hMC3->Integral();
    float n_total = n_hMC1 + n_hMC2 + n_hMC3;

    hMC1->Scale(1.0/n_hMC1);
    hMC2->Scale(1.0/n_hMC2);
    hMC3->Scale(1.0/n_hMC3);
    // Create template list
    TObjArray *mcTemplates = new TObjArray(3);
    mcTemplates->Add(hMC1);
    mcTemplates->Add(hMC2);
    mcTemplates->Add(hMC3);

    // Create fraction fitter
    TFractionFitter *fit = new TFractionFitter(hData, mcTemplates);
    fit->Constrain(0, 0.0, 1); // fraction between 0 and 1
    fit->Constrain(1, 0.0, 0.2);
    fit->Constrain(2, 0.0, 0.1);

    // Perform fit
    Int_t status = fit->Fit();
    if (status != 0) {
        Error("fit_Jpsi_with_templates", "Fit failed, status = %d", status);
        return;
    }

    // Get fractions and errors
    Double_t frac[3], err[3];
    for (int i = 0; i < 3; ++i) {
        fit->GetResult(i, frac[i], err[i]);
        printf("Process %d: fraction = %.4f ± %.4f\n", i+1, frac[i], err[i]);
        printf("Yield: %.2f ± %.2f\n", frac[i] * hData->Integral(), err[i] * hData->Integral());
        //
    }

    // Clone MC histos so we don't overwrite originals
    TH1D *hMC1_fit = (TH1D*)hMC1->Clone("hMC1_fit");
    TH1D *hMC2_fit = (TH1D*)hMC2->Clone("hMC2_fit");
    TH1D *hMC3_fit = (TH1D*)hMC3->Clone("hMC3_fit");

    // Scale each by fraction × total data yield
    double totalData = hData->Integral();
    hMC1_fit->Scale(frac[0] * totalData / hMC1_fit->Integral());
    hMC2_fit->Scale(frac[1] * totalData / hMC2_fit->Integral());
    hMC3_fit->Scale(frac[2] * totalData / hMC3_fit->Integral());

    double Yield1 = hMC1_fit->Integral(hMC1_fit->FindBin(2.8),hMC1_fit->FindBin(3.2));
    double Yield2 = hMC2_fit->Integral(hMC2_fit->FindBin(2.8),hMC2_fit->FindBin(3.2));
    double Yield3 = hMC3_fit->Integral(hMC3_fit->FindBin(2.8),hMC3_fit->FindBin(3.2));

    printf("Yield 1 [2.8, 3.2]: %.2f; Fraction: %.2f \n", Yield1, Yield1/(Yield1+Yield2+Yield3));
    printf("Yield 2 [2.8, 3.2]: %.2f; Fraction: %.2f \n", Yield2, Yield2/(Yield1+Yield2+Yield3));
    printf("Yield 3 [2.8, 3.2]: %.2f; Fraction: %.2f \n", Yield3, Yield3/(Yield1+Yield2+Yield3));


    // Get total fit histogram
    TH1 *hFit = (TH1*)fit->GetPlot();
    hFit->SetLineColor(kRed);
    hFit->SetLineWidth(2);

    // Styling

    hData->GetXaxis()->SetTitle(titleX);
    hData->GetYaxis()->SetTitle(titleY);

    hFit->SetLineColor(1);
    hMC1_fit  ->SetLineColor(2);
    hMC2_fit  ->SetLineColor(4);
    hMC3_fit  ->SetLineColor(6);

    //hMC1_fit->SetLineColor(kBlue);
    hMC1_fit->SetLineWidth(2);

    //hMC2_fit->SetLineColor(kGreen+2);
    hMC2_fit->SetLineWidth(2);

    //hMC3_fit->SetLineColor(kOrange+1);
    hMC3_fit->SetLineWidth(2);


    // Draw
    TCanvas *c_fit = new TCanvas("c", "Fit result", 800, 600);
    hData->SetMarkerStyle(20);
    hData->Draw("E");           // Data points
    hFit->Draw("HIST SAME");    // Total fit
    hMC1_fit->Draw("HIST SAME");
    hMC2_fit->Draw("HIST SAME");
    hMC3_fit->Draw("HIST SAME");

    // Legend
    TLegend *leg = new TLegend(0.6, 0.55, 0.88, 0.85);
    leg->AddEntry(hData, "Data", "ep");
    leg->AddEntry(hFit, "Total Fit", "l");
    leg->AddEntry(hMC1_fit, "#gammap#uparrow#rightarrowJ/#psip","l");
    leg->AddEntry(hMC2_fit, "#gammaAu#rightarrowJ/#psiAu","l");
    leg->AddEntry(hMC3_fit, "#gamma#gamma#rightarrowe^{#plus}e^{#minus}","l");
    leg->SetTextSize(0.048);
    leg-> SetTextFont(42);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->Draw();

    if(Mark!="empty"){
        draw_Mark(Mark, 0.45, 0.85, 0.6, 0.9);
    }
    if(Mark2!="empty"){
        draw_Mark(Mark2, 0.45, 0.8, 0.6, 0.83);
    }
    char name[100];
    sprintf(name,"Total Fit = %d ", hFit->Integral(hFit->FindBin(2.8),hFit->FindBin(3.2)));
    draw_Mark(name, 0.45, 0.35, 0.6, 0.4);
    sprintf(name,"#gammap#uparrow#rightarrowJ/#psip = %d ", hMC1_fit->Integral(hMC1_fit->FindBin(2.8),hMC1_fit->FindBin(3.2)));
    draw_Mark(name, 0.45, 0.3, 0.6, 0.35);
    sprintf(name,"#gammaAu#rightarrowJ/#psiAu = %d ", hMC2_fit->Integral(hMC2_fit->FindBin(2.8),hMC2_fit->FindBin(3.2)));
    draw_Mark(name, 0.45, 0.25, 0.6, 0.30);
    sprintf(name,"#gamma#gamma#rightarrowe^{#plus}e^{#minus} = %d ", hMC3_fit->Integral(hMC3_fit->FindBin(2.8),hMC3_fit->FindBin(3.2)));
    draw_Mark(name, 0.45, 0.2, 0.6, 0.25);

    c_fit -> Print("./Plots/"+canvasName+"_fitting_peak.png");
    c_fit -> Print("./Plots/"+canvasName+"_fitting_peak.pdf");

}
void fit_Jpsi_with_templates_pp(TString canvasName,TH1F *hData, TH1F *hMC1, TH1F *hMC3, TString titleX, TString titleY, TString Mark="empty", TString Mark2="empty" ) {
    // Load histograms (could also be created in the same script)
    //TFile *f = TFile::Open("histos.root");
    //TH1D *hData = (TH1D*)f->Get("hData");     // data histogram
    //TH1D *hMC1  = (TH1D*)f->Get("hProc1");    // MC: process 1
    //TH1D *hMC2  = (TH1D*)f->Get("hProc2");    // MC: process 2
    //TH1D *hMC3  = (TH1D*)f->Get("hProc3");    // MC: process 3

    // Make sure binning matches exactly
    //if (!hData || !hMC1 || !hMC2 || !hMC3) {
    //    Error("fit_Jpsi_with_templates", "Missing histograms!");
    //    return;
    //}

    // Create template list
    TObjArray *mcTemplates = new TObjArray(3);
    mcTemplates->Add(hMC1);
    //mcTemplates->Add(hMC2);
    mcTemplates->Add(hMC3);

    // Create fraction fitter
    TFractionFitter *fit = new TFractionFitter(hData, mcTemplates);
    fit->Constrain(0, 0.0, 0.90); // fraction between 0 and 1
    fit->Constrain(1, 0.0, 0.1);
    //fit->Constrain(2, 0.0, 0.20);

    // Perform fit
    Int_t status = fit->Fit();
    if (status != 0) {
        Error("fit_Jpsi_with_templates", "Fit failed, status = %d", status);
        return;
    }

    // Get fractions and errors
    Double_t frac[2], err[2];
    for (int i = 0; i < 2; ++i) {
        fit->GetResult(i, frac[i], err[i]);
        printf("Process %d: fraction = %.4f ± %.4f\n", i+1, frac[i], err[i]);
        printf("Yield: %.2f ± %.2f\n", frac[i] * hData->Integral(), err[i] * hData->Integral());
        //
    }

    // Clone MC histos so we don't overwrite originals
    TH1D *hMC1_fit = (TH1D*)hMC1->Clone("hMC1_fit");
    //TH1D *hMC2_fit = (TH1D*)hMC2->Clone("hMC2_fit");
    TH1D *hMC3_fit = (TH1D*)hMC3->Clone("hMC3_fit");

    // Scale each by fraction × total data yield
    double totalData = hData->Integral();
    hMC1_fit->Scale(frac[0] * totalData / hMC1_fit->Integral());
    //hMC2_fit->Scale(frac[1] * totalData / hMC2_fit->Integral());
    hMC3_fit->Scale(frac[1] * totalData / hMC3_fit->Integral());

    double Yield1 = hMC1_fit->Integral(hMC1_fit->FindBin(2.8),hMC1_fit->FindBin(3.2));
    //double Yield2 = hMC2_fit->Integral(hMC2_fit->FindBin(2.8),hMC2_fit->FindBin(3.2));
    double Yield3 = hMC3_fit->Integral(hMC3_fit->FindBin(2.8),hMC3_fit->FindBin(3.2));

    printf("Yield 1 [2.8, 3.2]: %.2f; Fraction: %.2f \n", Yield1, Yield1/(Yield1+Yield3));
    //printf("Yield 2 [2.8, 3.2]: %.2f; Fraction: %.2f \n", Yield2, Yield2/(Yield1+Yield2+Yield3));
    printf("Yield 3 [2.8, 3.2]: %.2f; Fraction: %.2f \n", Yield3, Yield3/(Yield1+Yield3));


    // Get total fit histogram
    TH1 *hFit = (TH1*)fit->GetPlot();
    hFit->SetLineColor(kRed);
    hFit->SetLineWidth(2);

    // Styling

    hData->GetXaxis()->SetTitle(titleX);
    hData->GetYaxis()->SetTitle(titleY);

    hFit->SetLineColor(1);
    hMC1_fit  ->SetLineColor(2);
    //hMC2_fit  ->SetLineColor(4);
    hMC3_fit  ->SetLineColor(6);

    //hMC1_fit->SetLineColor(kBlue);
    hMC1_fit->SetLineWidth(2);

    //hMC2_fit->SetLineColor(kGreen+2);
    //hMC2_fit->SetLineWidth(2);

    //hMC3_fit->SetLineColor(kOrange+1);
    hMC3_fit->SetLineWidth(2);


    // Draw
    TCanvas *c_fit = new TCanvas("c", "Fit result", 800, 600);
    hData->SetMarkerStyle(20);
    hData->Draw("E");           // Data points
    hFit->Draw("HIST SAME");    // Total fit
    hMC1_fit->Draw("HIST SAME");
    //hMC2_fit->Draw("HIST SAME");
    hMC3_fit->Draw("HIST SAME");

    // Legend
    TLegend *leg = new TLegend(0.6, 0.55, 0.88, 0.85);
    leg->AddEntry(hData, "Data", "ep");
    leg->AddEntry(hFit, "Total Fit", "l");
    leg->AddEntry(hMC1_fit, "#gammap#uparrow#rightarrowJ/#psip","l");
    //leg->AddEntry(hMC2_fit, "#gammaAu#rightarrowJ/#psiAu","l");
    leg->AddEntry(hMC3_fit, "#gamma#gamma#rightarrowe^{#plus}e^{#minus}","l");
    leg->SetTextSize(0.048);
    leg-> SetTextFont(42);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->Draw();

    if(Mark!="empty"){
        draw_Mark(Mark, 0.45, 0.85, 0.6, 0.9);
    }
    if(Mark2!="empty"){
        draw_Mark(Mark2, 0.45, 0.8, 0.6, 0.83);
    }
    char name[100];
    sprintf(name,"Total Fit = %d ", hFit->Integral(hFit->FindBin(2.8),hFit->FindBin(3.2)));
    draw_Mark(name, 0.45, 0.35, 0.6, 0.4);
    sprintf(name,"#gammap#uparrow#rightarrowJ/#psip = %d ", hMC1_fit->Integral(hMC1_fit->FindBin(2.8),hMC1_fit->FindBin(3.2)));
    //draw_Mark(name, 0.45, 0.3, 0.6, 0.35);
    //sprintf(name,"#gammaAu#rightarrowJ/#psiAu = %d ", hMC2_fit->Integral(hMC2_fit->FindBin(2.8),hMC2_fit->FindBin(3.2)));
    draw_Mark(name, 0.45, 0.25, 0.6, 0.30);
    sprintf(name,"#gamma#gamma#rightarrowe^{#plus}e^{#minus} = %d ", hMC3_fit->Integral(hMC3_fit->FindBin(2.8),hMC3_fit->FindBin(3.2)));
    draw_Mark(name, 0.45, 0.2, 0.6, 0.25);

    c_fit -> Print("./Plots/"+canvasName+"_fitting_peak.png");
    c_fit -> Print("./Plots/"+canvasName+"_fitting_peak.pdf");

}

void fit_background_with_gap(TH2F *h, TString canvasName, TString titleX, TString titleY, TString Mark="empty", TString Mark2="empty", float minx=0, float maxx=-10000, float miny=0, float maxy=-10000, int Nrebin=-1){
    h->GetXaxis()->SetRangeUser(2,6);
    TH1F *hist = (TH1F*)h->ProjectionY()->Clone("h_"+canvasName+"_p");

    TCanvas *c = new TCanvas("c_"+canvasName,"c_"+canvasName,800,600);

    // Simulate histogram with exponential background + gaussian peak
    //TH1F *hist = new TH1F("hist", "Histogram with Exponential + Gaussian", 40, 2, 6);

    // Fill exponential background
    //TF1 *expTrue = new TF1("expTrue", "expo", 2, 6);
    //expTrue->SetParameters(3, -0.1); // amplitude, exponent
    //hist->FillRandom("expTrue", 1000);

    // Add gaussian peak
    //TF1 *gaus = new TF1("gaus", "gaus", 2.5, 4.5);
    //gaus->SetParameters(200, 3, 0.1); // amplitude, mean, sigma
    //hist->FillRandom("gaus", 400);
    if(Nrebin>0){
        hist   ->Rebin(Nrebin);
    }
    if(maxx!=-10000){
        hist->GetXaxis()->SetRangeUser(minx, maxx);
    }
    hist->GetYaxis()->SetRangeUser(hist->GetMinimum(), 1.3*hist->GetMaximum());
    if(maxy!=-10000){
        hist->GetYaxis()->SetRangeUser(miny, maxy);
    }
    hist->GetXaxis()->SetTitle(titleX);
    hist->GetYaxis()->SetTitle(titleY);

    // Draw the histogram
    hist->Draw("E1X0");

    // --- Fit the background excluding the peak region ---
    // Define exponential function
    TF1 *expFit = new TF1("expFit", "expo", 1, 10);

    // Fit only sidebands: from 0–4.5 and 5.5–10
    hist->Fit("expFit", "REM0", "", 3.4, 6); // right side (use "+" to append to previous fit)
    hist->Fit("expFit", "+REM0", "", 1.9, 2.5);   // left side

    // Draw the fitted exponential
    expFit->SetLineColor(kRed);
    expFit->SetLineStyle(7);
    expFit->Draw("same");


    // Define full fit function
    TF1 *fitFunc = new TF1("fitFunc", totalFit, 2.0, 6.0, 7);
    fitFunc->SetParNames("alpha", "n", "mean", "sigma", "CBnorm", "BKGnorm", "BKGslope");

    // Initial parameter guesses
    fitFunc->SetParameters(1e-2, 3.0, 3.1, 1.1e-3, 400, 5, -0.8);

    // Fit
    hist->Fit(fitFunc, "ER", "", 2.1, 5.0);
    //fitFunc->SetParameters(fitFunc->GetParameters())
    //hist->Fit(fitFunc, "ERM", "", 2.1, 5.0);
    //fitFunc->SetParameters(fitFunc->GetParameters())
    //hist->Fit(fitFunc, "ERM", "", 2.1, 5.0);

    // Background estimate in window (using background part only)
    TF1 *bkgOnly = new TF1("bkgOnly", "([0]*exp([1]*x))", 2.8, 3.2);
    bkgOnly->SetParameters(fitFunc->GetParameter(5), fitFunc->GetParameter(6));
    double bkgCounts = bkgOnly->Integral(2.8, 3.2);
    

    // Draw the fitted fitFunc function
    fitFunc->SetLineColor(kBlue);
    fitFunc->Draw("same");

    // Optional: draw Gaussian separately
    //gaus->SetLineColor(kBlue);
    //gaus->Draw("same");

    // Add legend
    TLegend *leg = new TLegend(0.55, 0.65, 0.85, 0.85);
    leg->AddEntry(hist, "Data", "p");
    leg->AddEntry(expFit, "Exp Fit", "l");
    leg->AddEntry(fitFunc, "Crystal Ball + Exp Fit", "l");
    //leg->AddEntry(gaus, "Injected Gaussian Peak", "l");
    leg->SetTextSize(0.048);
    leg-> SetTextFont(42);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->Draw();
    char name[100];
    sprintf(name,"N = %d - %d(%f)", hist->Integral(hist->FindBin(2.8),hist->FindBin(3.2)),expFit->Integral(2.8,3.2), bkgCounts);
    draw_Mark(name, 0.55, 0.55, 0.8, 0.6);
    if(Mark!="empty"){
        draw_Mark(Mark, 0.2, 0.85, 0.3, 0.9);
    }
    if(Mark2!="empty"){
        draw_Mark(Mark2, 0.2, 0.8, 0.3, 0.83);
    }


    c -> Print("./Plots/"+canvasName+".png");

}

void draw_1d(TH1F *h, TString canvasName, TString titleX, TString titleY, TString Mark="empty", TString Mark2="empty", float minx=0, float maxx=-10000, float miny=0, float maxy=-10000, bool fLogY=false){
    TCanvas *c = new TCanvas("c_"+canvasName,"c_"+canvasName,800,600);
    if(fLogY){
        c->SetLogy();
    }
    h -> GetXaxis() -> SetTitle(titleX);
    h -> GetYaxis() -> SetTitle(titleY);

    if(maxx!=-10000){
        h->GetXaxis()->SetRangeUser(minx, maxx);
    }
    if(maxy!=-10000){
        h->GetYaxis()->SetRangeUser(miny, maxy);
    }

    h -> Draw("E1X0");

    if(Mark!="empty"){
        draw_Mark(Mark, 0.2, 0.85, 0.3, 0.9);
    }
    if(Mark2!="empty"){
        draw_Mark(Mark2, 0.2, 0.8, 0.3, 0.83);
    }
    c -> Print("./Plots/"+canvasName+".png");
    c -> Print("./Plots/"+canvasName+".pdf");
}

void draw_1d_z(TH1F *h, TString canvasName, TString titleX, TString titleY, TString Mark="empty", TString Mark2="empty", float minx=0, float maxx=-10000, float miny=0, float maxy=-10000, bool fLogY=false){
    TCanvas *c = new TCanvas("c_"+canvasName,"c_"+canvasName,800,600);
    if(fLogY){
        c->SetLogy();
    }
    h -> GetXaxis() -> SetTitle(titleX);
    h -> GetYaxis() -> SetTitle(titleY);

    if(maxx!=-10000){
        h->GetXaxis()->SetRangeUser(minx, maxx);
    }
    if(maxy!=-10000){
        h->GetYaxis()->SetRangeUser(miny, maxy);
    }
    TH1F *h_copy = (TH1F*)h->Clone("h_zero");
    h -> Draw("E1X0");
    h_copy -> GetXaxis()->SetRangeUser(-100, 100);
    h_copy -> SetFillColorAlpha(kBlue, 0.5);
    h_copy -> SetLineColor(kBlue);
    h_copy -> SetLineWidth(2);
    h_copy -> SetLineStyle(7);
    h_copy -> Draw("HISTsame");

    if(Mark!="empty"){
        draw_Mark(Mark, 0.2, 0.85, 0.3, 0.9);
    }
    if(Mark2!="empty"){
        draw_Mark(Mark2, 0.2, 0.8, 0.3, 0.83);
    }
    c -> Print("./Plots/"+canvasName+".png");
    c -> Print("./Plots/"+canvasName+".pdf");
}

void MC_background(TH2F *h, TH2F *h_SS, TString canvasName, TString titleX, TString titleY, TString Mark="empty", TString Mark2="empty", float minx=0, float maxx=-10000, float miny=0, float maxy=-10000, int Nrebin=-1){
    h->GetXaxis()->SetRangeUser(2,6);
    TH1F *hist = (TH1F*)h->ProjectionY()->Clone("h_"+canvasName+"_p");
    TH1F *hist_SS = (TH1F*)h_SS->ProjectionY()->Clone("h_"+canvasName+"_SS_p");
    TFile *inputFile_QED = new TFile("./JPsi_pAu_QED_p_100_Au_100.root");
    TFile *inputFile_p = new TFile("./JPsi_pAu_pSrc_100_AuTgt_100.root");
    TFile *inputFile_A = new TFile("./JPsi_pAu_AuSrc_100_pTgt_100.root");

    //TFile *inputFile = new TFile("/star/u/wschmidk/starana/UPC/pA_2015/MCplay/root/fit_mpthists.root");
    TH2D*     h2MCsum = (TH2D*)inputFile_A  ->Get("h_dca_spin_mee")->Clone("h2MCsum");
    TH2D*     h2mpt_1 = (TH2D*)inputFile_A  ->Get("h_dca_spin_mee")->Clone("h2mpt_1");
    TH2D*     h2mpt_2 = (TH2D*)inputFile_p  ->Get("h_dca_spin_mee")->Clone("h2mpt_2");
    TH2D*     h2mpt_3 = (TH2D*)inputFile_QED->Get("h_dca_spin_mee")->Clone("h2mpt_3");

    h2MCsum ->Sumw2();
    h2MCsum ->Add(h2mpt_2);
    h2MCsum ->Add(h2mpt_3);

    TH1F *h_m_sum  = (TH1F*)h2MCsum -> ProjectionY()->Clone("h_m_sum");
    TH1F *h_m_1    = (TH1F*)h2mpt_1 -> ProjectionY()->Clone("h_m_1");
    TH1F *h_m_2    = (TH1F*)h2mpt_2 -> ProjectionY()->Clone("h_m_2");
    TH1F *h_m_3    = (TH1F*)h2mpt_3 -> ProjectionY()->Clone("h_m_3");


    fit_Jpsi_with_templates(canvasName, hist, h_m_1 , h_m_2, h_m_3, titleX, titleY, Mark, Mark2);

    TCanvas *c = new TCanvas("c_"+canvasName,"c_"+canvasName,800,600);

    // Simulate histogram with exponential background + gaussian peak
    //TH1F *hist = new TH1F("hist", "Histogram with Exponential + Gaussian", 40, 2, 6);

    // Fill exponential background
    //TF1 *expTrue = new TF1("expTrue", "expo", 2, 6);
    //expTrue->SetParameters(3, -0.1); // amplitude, exponent
    //hist->FillRandom("expTrue", 1000);

    // Add gaussian peak
    //TF1 *gaus = new TF1("gaus", "gaus", 2.5, 4.5);
    //gaus->SetParameters(200, 3, 0.1); // amplitude, mean, sigma
    //hist->FillRandom("gaus", 400);
    if(Nrebin>0){
        hist   ->Rebin(Nrebin);
    }
    if(maxx!=-10000){
        hist->GetXaxis()->SetRangeUser(minx, maxx);
    }
    hist->GetYaxis()->SetRangeUser(hist->GetMinimum(), 1.3*hist->GetMaximum());
    if(maxy!=-10000){
        hist->GetYaxis()->SetRangeUser(miny, maxy);
    }
    hist->GetXaxis()->SetTitle(titleX);
    hist->GetYaxis()->SetTitle(titleY);

    hist_SS ->SetLineColor(7);
    hist_SS ->SetMarkerColor(7);

    h_m_sum ->SetLineColor(1);
    h_m_1   ->SetLineColor(2);
    h_m_2   ->SetLineColor(4);
    h_m_3   ->SetLineColor(6);


    // Draw the histogram
    hist->Draw("E1X0");
    hist_SS->Draw("E1X0same");
    h_m_sum->Draw("HISTsame");
    h_m_1  ->Draw("HISTsame");
    h_m_2  ->Draw("HISTsame");
    h_m_3  ->Draw("HISTsame");


    // Optional: draw Gaussian separately
    //gaus->SetLineColor(kBlue);
    //gaus->Draw("same");

    // Add legend
    TLegend *leg_mc = new TLegend(0.6, 0.58, 0.9, 0.8);
    leg_mc->AddEntry(hist, "Data", "p");
    leg_mc->AddEntry(hist_SS, "Data SS", "p");
    leg_mc->AddEntry(h_m_sum,"MC sum","l");
    leg_mc->AddEntry(h_m_1  ,"#gammap#uparrow#rightarrowJ/#psip","l");
    leg_mc->AddEntry(h_m_2  ,"#gammaAu#rightarrowJ/#psiAu","l");
    leg_mc->AddEntry(h_m_3  ,"#gamma#gamma#rightarrowe^{#plus}e^{#minus}","l");
    //leg_mc->AddEntry(gaus, "Injected Gaussian Peak", "l");
    leg_mc->SetTextSize(0.048);
    leg_mc-> SetTextFont(42);
    leg_mc->SetFillColor(0);
    leg_mc->SetFillStyle(0);
    leg_mc->SetBorderSize(0);
    leg_mc->Draw();
    char name[100];
    //sprintf(name,"N = %d - %d(%f)", hist->Integral(hist->FindBin(2.8),hist->FindBin(3.2)),expFit->Integral(2.8,3.2), bkgCounts);
    //draw_Mark(name, 0.6, 0.65, 0.88, 0.7);
    if(Mark!="empty"){
        draw_Mark(Mark, 0.2, 0.85, 0.3, 0.9);
    }
    if(Mark2!="empty"){
        draw_Mark(Mark2, 0.2, 0.8, 0.3, 0.83);
    }


    c -> Print("./Plots/"+canvasName+".png");

}

void MC_background_pp(TH2F *h, TString canvasName, TString titleX, TString titleY, TString Mark="empty", TString Mark2="empty", float minx=0, float maxx=-10000, float miny=0, float maxy=-10000, int Nrebin=-1){
    h->GetXaxis()->SetRangeUser(2,6);
    TH1F *hist = (TH1F*)h->ProjectionY()->Clone("h_"+canvasName+"_p");

    TFile *inputFile = new TFile("/star/u/wschmidk/starana/UPC/pA_2015/MCplay/root/fit_mpthists.root");
    TH2D*     h2MCsum = (TH2D*)inputFile->Get("h2MCsum")->Clone("h2MCsum");
    TH2D*     h2mpt_1 = (TH2D*)inputFile->Get("h2mpt_1")->Clone("h2mpt_1");
    //TH2D*     h2mpt_2 = (TH2D*)inputFile->Get("h2mpt_2")->Clone("h2mpt_2");
    TH2D*     h2mpt_3 = (TH2D*)inputFile->Get("h2mpt_3")->Clone("h2mpt_3");

    TH1F *h_m_sum_tmp  = (TH1F*)h2MCsum -> ProjectionX()->Clone("h_m_sum_tmp");
    TH1F *h_m_1_tmp    = (TH1F*)h2mpt_1 -> ProjectionX()->Clone("h_m_1_tmp");
    //TH1F *h_m_2_tmp    = (TH1F*)h2mpt_2 -> ProjectionX()->Clone("h_m_2_tmp");
    TH1F *h_m_3_tmp    = (TH1F*)h2mpt_3 -> ProjectionX()->Clone("h_m_3_tmp");

    TH1F *h_m_sum  = (TH1F*) hist->Clone("h_m_sum");
    TH1F *h_m_1    = (TH1F*) hist->Clone("h_m_1");
    //TH1F *h_m_2    = (TH1F*) hist->Clone("h_m_2");
    TH1F *h_m_3    = (TH1F*) hist->Clone("h_m_3");

    h_m_sum->Reset();    
    h_m_1  ->Reset();
    //h_m_2  ->Reset();
    h_m_3  ->Reset();  
    
    for(int i=1; i<=h_m_sum_tmp->GetNbinsX(); i++){
        int bin = hist->FindBin(h_m_sum_tmp->GetBinCenter(i));
        h_m_sum->SetBinContent(bin, h_m_sum_tmp->GetBinContent(i));
        h_m_1  ->SetBinContent(bin, h_m_1_tmp  ->GetBinContent(i));
        //h_m_2  ->SetBinContent(bin, h_m_2_tmp  ->GetBinContent(i));
        h_m_3  ->SetBinContent(bin, h_m_3_tmp  ->GetBinContent(i));
    }
    fit_Jpsi_with_templates_pp(canvasName, hist, h_m_1 , h_m_3, titleX, titleY, Mark, Mark2);

    TCanvas *c = new TCanvas("c_"+canvasName,"c_"+canvasName,800,600);

    // Simulate histogram with exponential background + gaussian peak
    //TH1F *hist = new TH1F("hist", "Histogram with Exponential + Gaussian", 40, 2, 6);

    // Fill exponential background
    //TF1 *expTrue = new TF1("expTrue", "expo", 2, 6);
    //expTrue->SetParameters(3, -0.1); // amplitude, exponent
    //hist->FillRandom("expTrue", 1000);

    // Add gaussian peak
    //TF1 *gaus = new TF1("gaus", "gaus", 2.5, 4.5);
    //gaus->SetParameters(200, 3, 0.1); // amplitude, mean, sigma
    //hist->FillRandom("gaus", 400);
    if(Nrebin>0){
        hist   ->Rebin(Nrebin);
    }
    if(maxx!=-10000){
        hist->GetXaxis()->SetRangeUser(minx, maxx);
    }
    hist->GetYaxis()->SetRangeUser(hist->GetMinimum(), 1.3*hist->GetMaximum());
    if(maxy!=-10000){
        hist->GetYaxis()->SetRangeUser(miny, maxy);
    }
    hist->GetXaxis()->SetTitle(titleX);
    hist->GetYaxis()->SetTitle(titleY);

    h_m_sum->SetLineColor(1);
    h_m_1  ->SetLineColor(2);
    //h_m_2  ->SetLineColor(4);
    h_m_3  ->SetLineColor(6);


    // Draw the histogram
    hist->Draw("E1X0");
    h_m_sum->Draw("HISTsame");
    h_m_1  ->Draw("HISTsame");
    //h_m_2  ->Draw("HISTsame");
    h_m_3  ->Draw("HISTsame");


    // Optional: draw Gaussian separately
    //gaus->SetLineColor(kBlue);
    //gaus->Draw("same");

    // Add legend
    TLegend *leg_mc = new TLegend(0.6, 0.58, 0.9, 0.8);
    leg_mc->AddEntry(hist, "Data", "p");
    leg_mc->AddEntry(h_m_sum,"MC sum","l");
    leg_mc->AddEntry(h_m_1  ,"#gammap#uparrow#rightarrowJ/#psip","l");
    //leg_mc->AddEntry(h_m_2  ,"#gammaAu#rightarrowJ/#psiAu","l");
    leg_mc->AddEntry(h_m_3  ,"#gamma#gamma#rightarrowe^{#plus}e^{#minus}","l");
    //leg_mc->AddEntry(gaus, "Injected Gaussian Peak", "l");
    leg_mc->SetTextSize(0.048);
    leg_mc-> SetTextFont(42);
    leg_mc->SetFillColor(0);
    leg_mc->SetFillStyle(0);
    leg_mc->SetBorderSize(0);
    leg_mc->Draw();
    char name[100];
    //sprintf(name,"N = %d - %d(%f)", hist->Integral(hist->FindBin(2.8),hist->FindBin(3.2)),expFit->Integral(2.8,3.2), bkgCounts);
    //draw_Mark(name, 0.6, 0.65, 0.88, 0.7);
    if(Mark!="empty"){
        draw_Mark(Mark, 0.2, 0.85, 0.3, 0.9);
    }
    if(Mark2!="empty"){
        draw_Mark(Mark2, 0.2, 0.8, 0.3, 0.83);
    }


    c -> Print("./Plots/"+canvasName+".png");

}

void raw_tot_corr_asymmetry_m(TH2F *h, TH2F *h_SS, TH2F *h_L, TH2F *h_R,TH2F *h_SS_L, TH2F *h_SS_R, TString canvasName, TString titleX, TString titleY, TString Mark="empty", TString Mark2="empty", TString name_begin="pp_", bool doBkgSubtraction=true){
    char name[100];

    h->GetXaxis()->SetRangeUser(p_bin, p_bin);
    TH1F *h_p = (TH1F*)h->ProjectionY()->Clone("h_"+canvasName+"_p");

    h_SS->GetXaxis()->SetRangeUser(p_bin, p_bin);
    TH1F *h_p_SS = (TH1F*)h_SS->ProjectionY()->Clone("h_"+canvasName+"_SS_p");

    h->GetXaxis()->SetRangeUser(n_bin, n_bin);
    TH1F *h_n = (TH1F*)h->ProjectionY()->Clone("h_"+canvasName+"_n");
    h_SS->GetXaxis()->SetRangeUser(n_bin, n_bin);
    TH1F *h_n_SS = (TH1F*)h_SS->ProjectionY()->Clone("h_"+canvasName+"_SS_n");

    // Left
    h_L->GetXaxis()->SetRangeUser(p_bin, p_bin);
    TH1F *h_p_L = (TH1F*)h_L->ProjectionY()->Clone("h_"+canvasName+"_p_L");
    h_SS_L->GetXaxis()->SetRangeUser(p_bin, p_bin);
    TH1F *h_p_SS_L = (TH1F*)h_SS_L->ProjectionY()->Clone("h_"+canvasName+"_SS_p_L");

    h_L->GetXaxis()->SetRangeUser(n_bin, n_bin);
    TH1F *h_n_L = (TH1F*)h_L->ProjectionY()->Clone("h_"+canvasName+"_n_L");
    h_SS_L->GetXaxis()->SetRangeUser(n_bin, n_bin);
    TH1F *h_n_SS_L = (TH1F*)h_SS_L->ProjectionY()->Clone("h_"+canvasName+"_SS_n_L");

    // Right
    h_R->GetXaxis()->SetRangeUser(p_bin, p_bin);
    TH1F *h_p_R = (TH1F*)h_R->ProjectionY()->Clone("h_"+canvasName+"_p_R");
    h_SS_R->GetXaxis()->SetRangeUser(p_bin, p_bin);
    TH1F *h_p_SS_R = (TH1F*)h_SS_R->ProjectionY()->Clone("h_"+canvasName+"_SS_p_R");

    h_R->GetXaxis()->SetRangeUser(n_bin, n_bin);
    TH1F *h_n_R = (TH1F*)h_R->ProjectionY()->Clone("h_"+canvasName+"_n_R");
    h_SS_R->GetXaxis()->SetRangeUser(n_bin, n_bin);
    TH1F *h_n_SS_R = (TH1F*)h_SS_R->ProjectionY()->Clone("h_"+canvasName+"_SS_n_R");
    h_p ->GetXaxis()->SetRangeUser(-1,1);
    
    //h_p_L    ->GetXaxis()->SetRangeUser(-1,1);
    //h_p_SS_L ->GetXaxis()->SetRangeUser(-1,1);
    //h_n_R    ->GetXaxis()->SetRangeUser(-1,1);
    //h_n_SS_R ->GetXaxis()->SetRangeUser(-1,1);
    //h_p_R    ->GetXaxis()->SetRangeUser(-1,1);
    //h_p_SS_R ->GetXaxis()->SetRangeUser(-1,1);
    //h_n_L    ->GetXaxis()->SetRangeUser(-1,1);
    //h_n_SS_L ->GetXaxis()->SetRangeUser(-1,1);
    if(doBkgSubtraction){
        h_p_L    ->Add(h_p_SS_L,-1);
        h_n_R    ->Add(h_n_SS_R,-1);
        h_p_R    ->Add(h_p_SS_R,-1);
        h_n_L    ->Add(h_n_SS_L,-1);
    }

    h_p_L    ->GetXaxis()->SetTitle("M_{ee} [GeV/c^{2}]");
    h_p_SS_L ->GetXaxis()->SetTitle("M_{ee} [GeV/c^{2}]");
    h_n_R    ->GetXaxis()->SetTitle("M_{ee} [GeV/c^{2}]");
    h_n_SS_R ->GetXaxis()->SetTitle("M_{ee} [GeV/c^{2}]");
    h_p_R    ->GetXaxis()->SetTitle("M_{ee} [GeV/c^{2}]");
    h_p_SS_R ->GetXaxis()->SetTitle("M_{ee} [GeV/c^{2}]");
    h_n_L    ->GetXaxis()->SetTitle("M_{ee} [GeV/c^{2}]");
    h_n_SS_L ->GetXaxis()->SetTitle("M_{ee} [GeV/c^{2}]");

    h_p_L    ->GetYaxis()->SetTitle("#frac{Entries}{0.1 [GeV/c^{2}]}");
    h_p_SS_L ->GetYaxis()->SetTitle("#frac{Entries}{0.1 [GeV/c^{2}]}");
    h_n_R    ->GetYaxis()->SetTitle("#frac{Entries}{0.1 [GeV/c^{2}]}");
    h_n_SS_R ->GetYaxis()->SetTitle("#frac{Entries}{0.1 [GeV/c^{2}]}");
    h_p_R    ->GetYaxis()->SetTitle("#frac{Entries}{0.1 [GeV/c^{2}]}");
    h_p_SS_R ->GetYaxis()->SetTitle("#frac{Entries}{0.1 [GeV/c^{2}]}");
    h_n_L    ->GetYaxis()->SetTitle("#frac{Entries}{0.1 [GeV/c^{2}]}");
    h_n_SS_L ->GetYaxis()->SetTitle("#frac{Entries}{0.1 [GeV/c^{2}]}");

    double int_min = 2.8;
    double int_max = 3.2;

    // Define exponential function
    TF1 *expFit[4];
    expFit[0] = new TF1("expFit0", "([0]+1.0/[1]*x+1.0/([2]*[2])*x*x)", 1, 10);
    expFit[1] = new TF1("expFit1", "([0]+1.0/[1]*x+1.0/([2]*[2])*x*x)", 1, 10);
    expFit[2] = new TF1("expFit2", "([0]+1.0/[1]*x+1.0/([2]*[2])*x*x)", 1, 10);
    expFit[3] = new TF1("expFit3", "([0]+1.0/[1]*x+1.0/([2]*[2])*x*x)", 1, 10); 

    // Define full fit function
    TF1 *fitFunc[4];
    fitFunc[0] = new TF1("fitFunc0", totalFitPol, 2.0, 6.0, 8);
    fitFunc[1] = new TF1("fitFunc1", totalFitPol, 2.0, 6.0, 8);
    fitFunc[2] = new TF1("fitFunc2", totalFitPol, 2.0, 6.0, 8);
    fitFunc[3] = new TF1("fitFunc3", totalFitPol, 2.0, 6.0, 8);

    // Background estimate in window (using background part only)
    TF1 *bkgOnly[4]; 
    bkgOnly[0] = new TF1("bkgOnly0", "([0]+1.0/[1]*x+1.0/([2]*[2])*x*x)", 2.8, 3.2);
    bkgOnly[1] = new TF1("bkgOnly1", "([0]+1.0/[1]*x+1.0/([2]*[2])*x*x)", 2.8, 3.2);
    bkgOnly[2] = new TF1("bkgOnly2", "([0]+1.0/[1]*x+1.0/([2]*[2])*x*x)", 2.8, 3.2);
    bkgOnly[3] = new TF1("bkgOnly3", "([0]+1.0/[1]*x+1.0/([2]*[2])*x*x)", 2.8, 3.2);

    for(int i=0;i<4;i++){
        //fitFunc[i]->SetParNames("alpha", "n", "mean", "sigma", "CBnorm", "BKGnorm", "BKGslope","B");
        // Initial parameter guesses
        //fitFunc[i]->SetParameters(1e-2, 3.0, 3.1, 1.1e-3, 400, 5, -0.8);
        fitFunc[i]->SetParameters(1, 1, Mjpsi, 6.1e-2, 3000, 2.4, 5e6, 4e4);
        //fitFunc[i]->FixParameter(2, Mjpsi);

        bkgOnly[i]->SetLineColor(kGreen);
        //bkgOnly[i]->SetLineStyle(7);
    
    }
    TCanvas *c_4p = new TCanvas("c_4p_"+canvasName,"c_4p_"+canvasName,800,600);
    c_4p->Divide(2 , 2 );
    c_4p->cd(1);
    h_p_L->GetYaxis()->SetRangeUser(0, h_p_L->GetMaximum()*1.2);
    h_p_L->Draw("E1X0");
    // --- Fit the background excluding the peak region ---
    // Fit only sidebands: from 0–4.5 and 5.5–10
    h_p_L->Fit("expFit3", "REM0", "", 3.3, 6); // right side (use "+" to append to previous fit)
    //h_p_L->Fit("expFit3", "+REM0", "", 2.0, 2.7);   // left side

    // Draw the fitted exponential
    expFit[0]->SetLineColor(kRed);
    expFit[0]->SetLineStyle(7);
    expFit[0]->Draw("same");

    //fitFunc[0]->SetParameters(1.3e-6, 1.0, 3.07, 1.8e-10, 70, 17, -0.6);

    // Fit
    h_p_L->Fit(fitFunc[0], "ERM", "", 2.5, 4.0);
    //h_p_L->Fit(fitFunc[0], "ERM", "", 2.3, 4.0);

    //// Optional: inspect results
    //for (int i = 0; i < fitFunc[0]->GetNpar(); i++) {
    //    fitFunc[0]->SetParameters(i, fitFunc[0]->GetParameter(i));
    //}
    //h_p_L->Fit(fitFunc[0], "ERMW", "", 2.0, 6.0);
    //h_p_L->Fit(fitFunc[0], "ERM", "", 2.0, 4.0);

    // Get the Exp. part of the fit
    bkgOnly[0]->SetParameters(fitFunc[0]->GetParameter(5), fitFunc[0]->GetParameter(6), fitFunc[0]->GetParameter(7));
    double bkgCounts = 0;
    if(doBkgSubtraction) {
        bkgCounts = bkgOnly[0]->Integral(int_min, int_max)/h_p_L->GetBinWidth(1);
    }

    // Draw the fitted fitFunc function
    fitFunc[0]->SetLineColor(kBlue);
    fitFunc[0]->Draw("same");
    bkgOnly[0]->Draw("same");
    h_p_L->GetListOfFunctions()->Clear();
    h_p_L->Draw("E1X0same");


    // Add legend
    TLegend *leg = new TLegend(0.55, 0.65, 0.85, 0.87);
    leg->AddEntry(h_p_L, "Data", "p");
    leg->AddEntry(expFit[0], "Exp Fit", "l");
    leg->AddEntry(fitFunc[0], "Crystal Ball + Exp Bkg.", "l");
    leg->AddEntry(bkgOnly[0], "Exp Bkg.", "l");
    //leg->AddEntry(gaus, "Injected Gaussian Peak", "l");
    leg->SetTextSize(0.048);
    leg-> SetTextFont(42);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->Draw();


    sprintf(name,"N = %.2f(#color[4]{%.2f}) - #color[2]{%.2f}(#color[3]{%.2f})", h_p_L->Integral(h_p_L->FindBin(int_min),h_p_L->FindBin(int_max)),fitFunc[0]->Integral(int_min,int_max)/h_p_L->GetBinWidth(1),expFit[0]->Integral(int_min,int_max)/h_p_L->GetBinWidth(1), bkgCounts);
    draw_Mark("#uparrow Left", 0.5, 0.55, 0.75, 0.6);
    draw_Mark(name, 0.4, 0.5, 0.6, 0.55);
    double N_L_p_fit = fitFunc[0]->Integral(int_min,int_max)/h_p_L->GetBinWidth(1) - bkgCounts;
    double N_L_p_fit_bkg = fitFunc[0]->Integral(int_min,int_max)/h_p_L->GetBinWidth(1);
    double N_L_p     = h_p_L->Integral(h_p_L->FindBin(int_min),h_p_L->FindBin(int_max)) - bkgCounts; 
    double N_L_p_bkg = h_p_L->Integral(h_p_L->FindBin(int_min),h_p_L->FindBin(int_max)); 
    c_4p->cd(2);
    h_n_R->GetYaxis()->SetRangeUser(0, h_n_R->GetMaximum()*1.2);
    h_n_R->Draw("E1X0");
    // --- Fit the background excluding the peak region ---
    h_n_R->Fit("expFit2", "REM0", "", 3.3, 4.1); // right side (use "+" to append to previous fit)
    //h_n_R->Fit("expFit2", "+REM0", "", 2.1, 2.7);   // left side

    // Draw the fitted exponential
    expFit[1]->SetLineColor(kRed);
    expFit[1]->SetLineStyle(7);
    expFit[1]->Draw("same");
    h_n_R->GetListOfFunctions()->Clear();
    h_n_R->Draw("E1X0same");
    // Initial parameter guesses
    //fitFunc[1]->SetParameters(1.3e-6, 1.0, 3.07, 1.8e-10, 70, 17, -0.6);
    // Fit
    h_n_R->Fit(fitFunc[1], "ERM", "", 2.5, 4.0);
    //h_n_R->Fit(fitFunc[1], "ERM", "", 2.3, 4.0);
    // Get the Exp. part of the fit
    bkgOnly[1]->SetParameters(fitFunc[1]->GetParameter(5), fitFunc[1]->GetParameter(6), fitFunc[1]->GetParameter(7));
    if(doBkgSubtraction) {
        bkgCounts = bkgOnly[1]->Integral(int_min,int_max)/h_n_R->GetBinWidth(1);
    } else {
        bkgCounts = 0;
    }
    // Draw the fitted fitFunc function
    fitFunc[1]->SetLineColor(kBlue);
    fitFunc[1]->Draw("same");
    bkgOnly[1]->Draw("same");
    sprintf(name,"N = %.2f(#color[4]{%.2f}) - #color[2]{%.2f}(#color[3]{%.2f})", h_n_R->Integral(h_n_R->FindBin(int_min),h_n_R->FindBin(int_max)),fitFunc[1]->Integral(int_min,int_max)/h_n_R->GetBinWidth(1),expFit[1]->Integral(int_min,int_max)/h_n_R->GetBinWidth(1), bkgCounts);
    draw_Mark("#downarrow Right", 0.5, 0.55, 0.75, 0.6);
    draw_Mark(name, 0.4, 0.5, 0.6, 0.55);
    if(Mark!="empty"){
        draw_Mark(Mark, 0.45, 0.85, 0.5, 0.9);
    }
    if(Mark2!="empty"){
        draw_Mark(Mark2, 0.45, 0.8, 0.5, 0.83);
    }
    double N_R_n_fit = fitFunc[1]->Integral(int_min,int_max)/h_n_R->GetBinWidth(1) - bkgCounts;
    double N_R_n_fit_bkg = fitFunc[1]->Integral(int_min,int_max)/h_n_R->GetBinWidth(1);
    double N_R_n     = h_n_R->Integral(h_n_R->FindBin(int_min),h_n_R->FindBin(int_max)) - bkgCounts; 
    double N_R_n_bkg = h_n_R->Integral(h_n_R->FindBin(int_min),h_n_R->FindBin(int_max)); 
    c_4p->cd(3);
    h_n_L->GetYaxis()->SetRangeUser(0, h_n_L->GetMaximum()*1.2);
    h_n_L->Draw("E1X0");
    // --- Fit the background excluding the peak region ---
    h_n_L->Fit("expFit0", "REM0", ""  , 3.4, 6); // right side (use "+" to append to previous fit)
    //h_n_L->Fit("expFit0", "+REM0", "" , 1.9, 2.5);   // left side

    // Draw the fitted exponential
    expFit[3]->SetLineColor(kRed);
    expFit[3]->SetLineStyle(7);
    expFit[3]->Draw("same");
    // Initial parameter guesses
    //fitFunc[3]->SetParameters(1e-3, 1.0, 3.07, 1.7e-5, 100, 17, -0.6);
    // Fit
    h_n_L->Fit(fitFunc[3],"ERM", "", 2.5, 4.0);
    //h_n_L->Fit(fitFunc[3],"ERM", "", 2.3, 4.0);
    // Get the Exp. part of the fit
    bkgOnly[3]->SetParameters(fitFunc[3]->GetParameter(5), fitFunc[3]->GetParameter(6), fitFunc[3]->GetParameter(7));
    if(doBkgSubtraction) {
        bkgCounts = bkgOnly[3]->Integral(int_min,int_max)/h_n_R->GetBinWidth(1);
    } else {
        bkgCounts = 0;
    }
    // Draw the fitted fitFunc function
    fitFunc[3]->SetLineColor(kBlue);
    fitFunc[3]->Draw("same");
    bkgOnly[3]->Draw("same");
    h_n_L->GetListOfFunctions()->Clear();
    h_n_L->Draw("E1X0same");
    sprintf(name,"N = %.2f(#color[4]{%.2f}) - #color[2]{%.2f}(#color[3]{%.2f})", h_n_L->Integral(h_n_L->FindBin(int_min),h_n_L->FindBin(int_max)),fitFunc[3]->Integral(int_min,int_max)/h_n_L->GetBinWidth(1),expFit[3]->Integral(int_min,int_max)/h_n_L->GetBinWidth(1), bkgCounts);
    draw_Mark("#downarrow Left", 0.5, 0.55, 0.75, 0.6);
    draw_Mark(name, 0.4, 0.5, 0.6, 0.55);

    double N_L_n_fit = fitFunc[3]->Integral(int_min,int_max)/h_n_R->GetBinWidth(1) - bkgCounts;
    double N_L_n_fit_bkg = fitFunc[3]->Integral(int_min,int_max)/h_n_R->GetBinWidth(1);
    double N_L_n     = h_n_L->Integral(h_n_L->FindBin(int_min),h_n_L->FindBin(int_max)) - bkgCounts;
    double N_L_n_bkg = h_n_L->Integral(h_n_L->FindBin(int_min),h_n_L->FindBin(int_max));
    c_4p->cd(4);
    h_p_R->GetYaxis()->SetRangeUser(0, h_p_R->GetMaximum()*1.2);
    h_p_R->Draw("E1X0");
    // --- Fit the background excluding the peak region ---
    h_p_R->Fit("expFit1", "REM0", "", 3.4, 6); // right side (use "+" to append to previous fit)
    //h_p_R->Fit("expFit1", "+REM0", "", 1.9, 2.5);   // left side
    // Draw the fitted exponential
    expFit[2]->SetLineColor(kRed);
    expFit[2]->SetLineStyle(7);
    expFit[2]->Draw("same");
    // Initial parameter guesses
    //fitFunc[2]->SetParameters(1.3e-6, 1.0, 3.07, 1.8e-10, 70, 17, -0.6);
    // Fit
    h_p_R->Fit(fitFunc[2], "ERM", "", 2.5, 4.0);
    //h_p_R->Fit(fitFunc[2], "ERM", "", 2.3, 4.0);
    // Get the Exp. part of the fit
    bkgOnly[2]->SetParameters(fitFunc[2]->GetParameter(5), fitFunc[2]->GetParameter(6), fitFunc[2]->GetParameter(7));
    if(doBkgSubtraction) {
        bkgCounts = bkgOnly[2]->Integral(int_min,int_max)/h_n_R->GetBinWidth(1);
    } else {
        bkgCounts = 0;
    }
    // Draw the fitted fitFunc function
    fitFunc[2]->SetLineColor(kBlue);
    fitFunc[2]->Draw("same");
    bkgOnly[2]->Draw("same");
    h_p_R->GetListOfFunctions()->Clear();
    h_p_R->Draw("E1X0same");

    sprintf(name,"N = %.2f(#color[4]{%.2f}) - #color[2]{%.2f}(#color[3]{%.2f})", h_p_R->Integral(h_p_R->FindBin(int_min),h_p_R->FindBin(int_max)),fitFunc[2]->Integral(int_min,int_max)/h_p_R->GetBinWidth(1),expFit[2]->Integral(int_min,int_max)/h_p_R->GetBinWidth(1), bkgCounts);
    draw_Mark("#uparrow Right", 0.5, 0.55, 0.75, 0.6);
    draw_Mark(name, 0.4, 0.5, 0.6, 0.55);

    double N_R_p_fit = fitFunc[2]->Integral(int_min,int_max)/h_n_R->GetBinWidth(1) - bkgCounts;
    double N_R_p_fit_bkg = fitFunc[2]->Integral(int_min,int_max)/h_n_R->GetBinWidth(1);
    double N_R_p     = h_p_R->Integral(h_p_R->FindBin(int_min),h_p_R->FindBin(int_max)) - bkgCounts; 
    double N_R_p_bkg = h_p_R->Integral(h_p_R->FindBin(int_min),h_p_R->FindBin(int_max)); 


    c_4p->Print("./Plots/c_4p_"+canvasName+".png");
    c_4p->Print("./Plots/c_4p_"+canvasName+".pdf");










    //double P_avg = 0.601;
    //double sqrt_cosPhiLR = 0.592;
    //double correction = 1.0/(P_avg*sqrt_cosPhiLR);
    //double N_L_p = h_p_L -> Integral() - h_p_SS_L -> Integral(); 
    //double N_R_n = h_n_R -> Integral() - h_n_SS_R -> Integral(); 
    //double N_R_p = h_p_R -> Integral() - h_p_SS_R -> Integral(); 
    //double N_L_n = h_n_L -> Integral() - h_n_SS_L -> Integral(); 
    double eps_fit   = correction*( sqrt(N_L_p_fit*N_R_n_fit) - sqrt(N_R_p_fit*N_L_n_fit) )/( sqrt(N_L_p_fit*N_R_n_fit) + sqrt(N_R_p_fit*N_L_n_fit) );
    double eps   = correction*( sqrt(N_L_p*N_R_n) - sqrt(N_R_p*N_L_n) )/( sqrt(N_L_p*N_R_n) + sqrt(N_R_p*N_L_n) );
    //double eps   = correction*( sqrt(N_R_p*N_L_n) - sqrt(N_L_p*N_R_n) )/( sqrt(N_R_p*N_L_n) + sqrt(N_L_p*N_R_n) );
    //double d_eps = correction*sqrt(N_L_p*N_R_n) * sqrt(N_R_p*N_L_n)/pow(sqrt(N_L_p*N_R_n) + sqrt(N_R_p*N_L_n),2)*sqrt(1/N_L_p+1/N_R_n+1/N_R_p+1/N_L_n);
    double A_fit = sqrt(N_L_p_fit*N_R_n_fit);
    double B_fit = sqrt(N_R_p_fit*N_L_n_fit);
    double d_eps_fit = correction * A_fit * B_fit /pow( A_fit + B_fit,2) * sqrt(N_L_p_fit_bkg/pow(N_L_p_fit,2)+N_R_n_fit_bkg/pow(N_R_n_fit,2)+N_R_p_fit_bkg/pow(N_R_p_fit,2)+N_L_n_fit_bkg/pow(N_L_n_fit,2));
    double A = sqrt(N_L_p*N_R_n);
    double B = sqrt(N_R_p*N_L_n);
    double d_eps = correction * A * B /pow( A + B,2) * sqrt(N_L_p_bkg/pow(N_L_p,2)+N_R_n_bkg/pow(N_R_n,2)+N_R_p_bkg/pow(N_R_p,2)+N_L_n_bkg/pow(N_L_n,2));
    //double eps = ( sqrt(N_L_p*N_L_n) - sqrt(N_R_p*N_R_n) )/( sqrt(N_L_p*N_L_n) + sqrt(N_R_p*N_R_n) );
    //double d_eps = sqrt(N_L_p*N_L_n) * sqrt(N_R_p*N_R_n)/pow(sqrt(N_L_p*N_L_n) + sqrt(N_R_p*N_R_n),2)*sqrt(1/N_L_p+1/N_R_n+1/N_R_p+1/N_L_n);

    double Ep = 250;
    if(name_begin=="pA_")Ep = 103.732; 
    double y = y_Mean;
    double w = sqrt(2*Ep*Mjpsi*exp(-1*y));

    // Define the single point (X, Y) and its errors
    double wErr = 0.001*w;

    std::cout<< "Integral:" << " | " << N_L_p << " | " << N_R_n << " | " << N_R_p << " | " << N_L_n << " | " << eps<< " | " << d_eps << std::endl;
    std::cout<< "A:" <<  eps<< " +/- " << d_eps << " | " << canvasName << std::endl;
    std::cout<< "A fit:" <<  eps_fit<< " +/- " << d_eps_fit << " | " << canvasName << std::endl;
    std::cout<<"<y>="<<y_Mean<<std::endl;
    std::cout<<"W="<<w<<std::endl;

    // Bill S results for comparison
    // THE A_N^gamma MEASUREMENT (from ../prelim_xasym_nevt/ANsum_sb.txt)
    double ANgam = 0.0502138;
    double eANgam = 0.199776;
    // @ mean W
    double Wmean = 23.8;    
    // TGraphErrors for data
    TGraphErrors* gdat = new TGraphErrors(1);
    gdat->SetPoint(0,Wmean,ANgam);
    gdat->SetPointError(0,0.,eANgam);
    gdat->SetMarkerColor(6);
    gdat->SetLineColor(6);
    gdat->SetMarkerStyle(20); gdat->SetMarkerSize(1.5);
    gdat->SetLineWidth(2);
    //gdat->Draw(""); return;


    double d_sys = sqrt(pow(0.09,2)+pow(0.02,2)+pow(0.025,2)+pow(0.025,2)); // 9% fit variation + 2% polarization + 2% cosPhiLR
    double w_sys = 1.0; // placeholder for future W systematic uncertainty
    double d_eps_sys = fabs(eps)*d_sys; 
    TGraphErrors *graph = new TGraphErrors(1, &w, &eps, &wErr, &d_eps);
    graph->SetName("graph");
    TGraphErrors *graph_fit = new TGraphErrors(1, &w, &eps_fit, &wErr, &d_eps_fit);
    graph_fit->SetName("graph_fit");

    TGraphErrors *graph_sys = new TGraphErrors(1, &w, &eps, &w_sys, &d_eps_sys);
    graph_sys->SetName("graph_sys");

    // TGraph for curve from Jakub
    TGraph* gJW = new TGraph("dat/AsymV4_BNL.dat");
    cout << "npoints curve: " << gJW->GetN() << endl;
    gJW->SetLineWidth(2); gJW->SetLineColor(2);

    TCanvas *c = new TCanvas("c_"+canvasName,"c_"+canvasName,800,600);

    graph_sys->SetFillColorAlpha(kGray+1, 0.4);   // transparent fill
    graph_sys->SetLineColor(kGray+2);
    graph_sys->SetMarkerSize(0);

    // Set marker style and color
    graph->SetMarkerStyle(20);  // Full circle
    //graph->SetMarkerSize(1.5);
    graph->SetMarkerColor(1);
    graph->SetLineColor(1);
    graph->SetLineWidth(2);

    graph_fit->SetMarkerStyle(20);  // Full circle
    //graph_fit->SetMarkerSize(1.5);
    graph_fit->SetMarkerColor(7);
    graph_fit->SetLineColor(7);
    graph_fit->SetLineWidth(1);

    // Set axis limits
    graph->GetXaxis()->SetLimits(4, 45);  // X-axis range
    graph->GetYaxis()->SetRangeUser(-0.2, 0.5);  // Y-axis range

    graph->GetXaxis()->SetTitle(titleX);
    graph->GetYaxis()->SetTitle(titleY);
    // Draw the graph with error bars
    graph->Draw("AP");  

    // Draw horizontal line at yLine
    TLine *hLine = new TLine(4, 0, 45, 0);
    hLine->SetLineColor(1);
    hLine->SetLineWidth(3);
    hLine->SetLineStyle(9);
    hLine->Draw("same");

    graph_sys->Draw("same 2");  
    graph->Draw("same P");  
    gJW->Draw("same C");
    gdat->Draw("same P");
    graph_fit->Draw("same P");

    if(Mark!="empty"){
        draw_Mark(Mark, 0.2, 0.85, 0.3, 0.9);
    }
    if(Mark2!="empty"){
        draw_Mark(Mark2, 0.2, 0.8, 0.3, 0.83);
    }
    // Display canvas
    c->Update();

    c->Print("./Plots/"+canvasName+".png");
    c->Print("./Plots/"+canvasName+".pdf");


    // Save to ROOT file
    TFile* file = new TFile("./Files/"+canvasName+".root", "RECREATE");
    graph->Write();  // Saves the graph with the name "myGraph"
    graph_sys->Write();
    file->Close();





}

void draw_projections(TH2F *h, TH2F *h_SS, TString canvasName, TString titleX, TString titleY, TString Mark="empty", TString Mark2="empty", float minx=0, float maxx=-10000, float miny=0, float maxy=-10000, int Nrebin=-1){
    h->GetXaxis()->SetRangeUser(p_bin, p_bin);
    TH1F *h_p = (TH1F*)h->ProjectionY()->Clone("h_"+canvasName+"_p");
    h_SS->GetXaxis()->SetRangeUser(p_bin, p_bin);
    TH1F *h_p_SS = (TH1F*)h_SS->ProjectionY()->Clone("h_"+canvasName+"_SS_p");

    h->GetXaxis()->SetRangeUser(n_bin, n_bin);
    TH1F *h_n = (TH1F*)h->ProjectionY()->Clone("h_"+canvasName+"_n");
    h_SS->GetXaxis()->SetRangeUser(n_bin, n_bin);
    TH1F *h_n_SS = (TH1F*)h_SS->ProjectionY()->Clone("h_"+canvasName+"_SS_n");

    TCanvas *c = new TCanvas("c_"+canvasName,"c_"+canvasName,800,600);
    h_p->SetLineColor(2);
    h_n->SetLineColor(4);
    h_p->SetMarkerColor(2);
    h_n->SetMarkerColor(4);
    h_p_SS->SetLineColor(2);
    h_n_SS->SetLineColor(4);
    h_p_SS->SetMarkerColor(2);
    h_n_SS->SetMarkerColor(4);
    h_p_SS->SetLineStyle(7);
    h_n_SS->SetLineStyle(7);
    h_p_SS->SetMarkerStyle(24);
    h_n_SS->SetMarkerStyle(24);

    h_p->GetXaxis()->SetTitle(titleX);
    h_p->GetYaxis()->SetTitle(titleY);
 
   
    if(Nrebin>0){
        h_p   ->Rebin(Nrebin);
        h_p_SS->Rebin(Nrebin);
        h_n   ->Rebin(Nrebin);
        h_n_SS->Rebin(Nrebin);
    }

    if(maxx!=-10000){
        h_p->GetXaxis()->SetRangeUser(minx, maxx);
    }
    if(maxy!=-10000){
        h_p->GetYaxis()->SetRangeUser(miny, maxy);
    }

    h_p   ->Draw("E1");
    //axis->Draw();
    h_n   ->Draw("E1same");
    h_p_SS->Draw("E1same");
    h_n_SS->Draw("E1same");

    TLegend* leg = new TLegend(0.72, 0.7, 0.99, 0.9, NULL, "brNDC");
    leg->AddEntry(h_p   , "e^{+}e^{-} #uparrow"  ,"l");
    leg->AddEntry(h_n   , "e^{+}e^{-} #downarrow","l");
    leg->AddEntry(h_p_SS, "e^{s}e^{s} #uparrow"  ,"l");
    leg->AddEntry(h_n_SS, "e^{s}e^{s} #downarrow","l");
    leg->SetTextSize(0.048);
    leg-> SetTextFont(42);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->Draw();

    if(Mark!="empty"){
        draw_Mark(Mark, 0.2, 0.85, 0.3, 0.9);
    }
    if(Mark2!="empty"){
        draw_Mark(Mark2, 0.2, 0.8, 0.3, 0.83);
    }

    c->Print("./Plots/"+canvasName+".png");


}

void draw_tot_projections(TH2F *h, TH2F *h_SS, TString canvasName, TString titleX, TString titleY, TString Mark="empty", TString Mark2="empty", float minx=0, float maxx=-10000, float miny=0, float maxy=-10000, int Nrebin=-1){
    h->GetXaxis()->SetRangeUser(2,6);
    TH1F *h_p = (TH1F*)h->ProjectionY()->Clone("h_"+canvasName+"_p");
    h_SS->GetXaxis()->SetRangeUser(2,6);
    TH1F *h_p_SS = (TH1F*)h_SS->ProjectionY()->Clone("h_"+canvasName+"_SS_p");

    //h->GetXaxis()->SetRangeUser(5,5);
    //TH1F *h_n = (TH1F*)h->ProjectionY()->Clone("h_"+canvasName+"_n");
    //h_SS->GetXaxis()->SetRangeUser(5,5);
    //TH1F *h_n_SS = (TH1F*)h_SS->ProjectionY()->Clone("h_"+canvasName+"_SS_n");

    TCanvas *c = new TCanvas("c_"+canvasName,"c_"+canvasName,800,600);



    h_p   ->SetLineColor(2);
    h_p_SS->SetLineColor(4);
    h_p   ->SetMarkerColor(2);
    h_p_SS->SetMarkerColor(4);
    h_p   ->SetMarkerStyle(20);
    h_p_SS->SetMarkerStyle(24);   

    h_p->GetXaxis()->SetTitle(titleX);
    h_p->GetYaxis()->SetTitle(titleY);
    
    
    h_p->GetYaxis()->SetRangeUser(h_p->GetMinimum(), 1.3*h_p->GetMaximum());



    if(Nrebin>0){
        h_p   ->Rebin(Nrebin);
        h_p_SS->Rebin(Nrebin);

    }
    if(maxx!=-10000){
        h_p->GetXaxis()->SetRangeUser(minx, maxx);
    }
    if(maxy!=-10000){
        h_p->GetYaxis()->SetRangeUser(miny, maxy);
    }

    h_p   ->Draw("E1X0");
    //h_n   ->Draw("same");
    h_p_SS->Draw("E1X0same");
    //h_n_SS->Draw("same");

    TLegend* leg = new TLegend(0.72, 0.7, 0.99, 0.8, NULL, "brNDC");
    leg->AddEntry(h_p   , "e^{+}e^{-}"  ,"p");
    //leg->AddEntry(h_n   , "e^{+}e^{-} #downarrow","l");
    leg->AddEntry(h_p_SS, "e^{s}e^{s}"  ,"p");
    //leg->AddEntry(h_n_SS, "e^{s}e^{s} #downarrow","l");
    leg->SetTextSize(0.048);
    leg-> SetTextFont(42);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->Draw();
    if(Mark!="empty"){
        draw_Mark(Mark, 0.2, 0.85, 0.3, 0.9);
    }
    if(Mark2!="empty"){
        draw_Mark(Mark2, 0.2, 0.8, 0.3, 0.83);
    }
    //if(Mark2!="empty"){
    char name[100];
    sprintf(name,"Mean: %.3f", h_p->GetMean());


    draw_Mark(name, 0.72, 0.6, 0.8, 0.7);
    //}

    c->Print("./Plots/"+canvasName+".png");


}

void draw_updown_projections(TH2F *h, TH2F *h_d, TString canvasName, TString titleX, TString titleY, TString Mark="empty", TString Mark2="empty", float minx=0, float maxx=-10000, float miny=0, float maxy=-10000, int Nrebin=-1){
    //h->GetXaxis()->SetRangeUser(2,6);
    //TH1F *h_p = (TH1F*)h->ProjectionY()->Clone("h_"+canvasName+"_p");
    //h_d->GetXaxis()->SetRangeUser(2,6);
    //TH1F *h_p_d = (TH1F*)h_d->ProjectionY()->Clone("h_"+canvasName+"_d_p");
    TH1F *h_ax = new TH1F("h_ax","h_ax",400,-2,2);

    h->GetXaxis()->SetRangeUser(-2,-0.001);
    TH1F *h_l = (TH1F*)h->ProjectionX()->Clone("h_"+canvasName+"_l");
    h->GetXaxis()->SetRangeUser(0.001,2);
    TH1F *h_r = (TH1F*)h->ProjectionX()->Clone("h_"+canvasName+"_r");

    h_d->GetXaxis()->SetRangeUser(-2,-0.001);
    TH1F *h_d_l = (TH1F*)h_d->ProjectionX()->Clone("h_d_"+canvasName+"_l");
    h_d->GetXaxis()->SetRangeUser(0.001,2);
    TH1F *h_d_r = (TH1F*)h_d->ProjectionX()->Clone("h_d_"+canvasName+"_r");

    TCanvas *c = new TCanvas("c_"+canvasName,"c_"+canvasName,800,600);



    h_l   ->SetLineColor(2);
    h_r   ->SetLineColor(4);
    h_l   ->SetMarkerColor(2);
    h_r   ->SetMarkerColor(4);
    h_l   ->SetMarkerStyle(20);
    h_r   ->SetMarkerStyle(24);   

    h_d_l   ->SetLineColor(1);
    h_d_r   ->SetLineColor(6);
    h_d_l   ->SetMarkerColor(1);
    h_d_r   ->SetMarkerColor(6);
    h_d_l   ->SetMarkerStyle(20);
    h_d_r   ->SetMarkerStyle(24);   


    
    
    h_l->GetYaxis()->SetRangeUser(h_l->GetMinimum(), 1.3*h_l->GetMaximum());



    if(Nrebin>0){
        h_l   ->Rebin(Nrebin);
        h_r->Rebin(Nrebin);

    }
    if(maxx!=-10000){
        h_ax->GetXaxis()->SetRangeUser(minx, maxx);
    }
    if(maxy!=-10000){
        h_ax->GetYaxis()->SetRangeUser(miny, maxy);
    }

    //h_ax->GetXaxis()->SetRangeUser(-1.2, 1.2);
    //h_ax->GetYaxis()->SetRangeUser(0,30);

    h_ax->GetXaxis()->SetTitle(titleX);
    h_ax->GetYaxis()->SetTitle(titleY);

    h_ax    ->Draw("AXIS");
    h_l     ->Draw("HISTsame");
    h_d_l   ->Draw("HISTsame");
    h_r     ->Draw("HISTsame");
    h_d_r   ->Draw("HISTsame");

    std::cout<< "L up "<<h_l  ->GetMean()<<std::endl;
    std::cout<< "R up "<<h_r  ->GetMean()<<std::endl;
    std::cout<< "L down "<<h_d_l     ->GetMean()<<std::endl;
    std::cout<< "R down "<<h_d_r   ->GetMean()<<std::endl;
    //if(Mark2!="empty"){
    char name[100];
    sprintf(name,"up <cos#phi_{L}>: %.3f", h_l->GetMean());
    draw_Mark(name, 0.32, 0.7, 0.4, 0.75);
    sprintf(name,"up <cos#phi_{R}>: %.3f", h_r   ->GetMean());
    draw_Mark(name, 0.32, 0.65, 0.4, 0.7);
    sprintf(name,"down <cos#phi_{L}>: %.3f", h_d_l     ->GetMean());
    draw_Mark(name, 0.32, 0.6, 0.4, 0.65);
    sprintf(name,"down <cos#phi_{R}>: %.3f", h_d_r   ->GetMean());
    draw_Mark(name, 0.32, 0.55, 0.4, 0.6);
    double avg_cosphi = sqrt(
        (fabs(h_l->GetMean())+fabs(h_d_l->GetMean()))
        *(fabs(h_r->GetMean()) + fabs(h_d_r->GetMean()))
    )/2;
    sqrt_cosPhiLR = avg_cosphi;
    sprintf(name,"#sqrt{<cos#phi_{L}><cos#phi_{R}>}: %.3f", avg_cosphi);
    draw_Mark(name, 0.32, 0.47, 0.4, 0.53);

    TLegend* leg = new TLegend(0.72, 0.65, 0.99, 0.85, NULL, "brNDC");
    leg->AddEntry(h_l   , "e^{+}e^{-} L"  ,"l");
    leg->AddEntry(h_d_l   , "e^{+}e^{-} #downarrow L","l");
    leg->AddEntry(h_r, "e^{s}e^{s} R"  ,"l");
    leg->AddEntry(h_d_r, "e^{s}e^{s} #downarrow R","l");
    leg->SetTextSize(0.048);
    leg-> SetTextFont(42);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->Draw();
    if(Mark!="empty"){
        draw_Mark(Mark, 0.2, 0.85, 0.3, 0.9);
    }
    if(Mark2!="empty"){
        draw_Mark(Mark2, 0.2, 0.8, 0.3, 0.83);
    }

    c->Print("./Plots/"+canvasName+".png");


}

void draw_mee_fast(){

    //TFile *inputRootFile = new TFile("./hist_pA_piP_UPC_mkkWithTrigger_cuts_new.root"); TString name_begin = "pA_";
    TFile *inputRootFile = new TFile("./hist_pp_piP_UPC_mkkWithTrigger_cuts_new.root"); TString name_begin = "pp_";

    if(name_begin == "pA_"){
        y_Mean = 0.154;
        P_avg = 0.601;
        sqrt_cosPhiLR = 0.592;
    }else{
        y_Mean = 0.0;
        P_avg = 0.5984;
        sqrt_cosPhiLR = 0.627;
    }

    // Fill polarisation information into a map for easy access later when we want to plot vs fill number
    std::map<int, std::vector<double>> fillMap;
    if(name_begin == "pA_"){

        ifstream file("polar/dat/fill_poll.csv");
        if (!file.is_open()) {
            cerr << "Error opening file: " << filename << endl;
            return;
        }

        string line;
        getline(file, line); // skip header line

        while (getline(file, line)) {
            if (line.empty()) continue;

            stringstream ss(line);
            //string fill_str, pol_str, pol_err_str;
            string fill_str,energy_str,start_str,stop_str,pol_str,pol_err_str,dpdt_str,dpdt_err_str;

            getline(ss, fill_str, ',');
            getline(ss, energy_str, ',');
            getline(ss, start_str, ',');
            getline(ss, stop_str, ',');
            getline(ss, pol_str, ',');
            getline(ss, pol_err_str, ',');
            getline(ss, dpdt_str, ',');
            getline(ss, dpdt_err_str, ',' );

            //std::cout<< fill_str << "," << pol_str << "," << pol_err_str << std::endl;
            int fill = atoi(fill_str.c_str());
            double energy = atof(energy_str.c_str());
            int start = atoi(start_str.c_str());
            int stop = atoi(stop_str.c_str());
            double pol = atof(pol_str.c_str());
            double pol_err = atof(pol_err_str.c_str());
            double dpdt = atof(dpdt_str.c_str());
            double dpdt_err = atof(dpdt_err_str.c_str());

            std::vector<double> vals;

            vals.push_back(energy);  //0
            vals.push_back(double(start)); //1
            vals.push_back(double(stop)); //2
            vals.push_back(pol); //3
            vals.push_back(pol_err); //4
            vals.push_back(dpdt); //5
            vals.push_back(dpdt_err); //6

            fillMap[fill] = vals; // assign vector to map
            //               0       1                2            3    4        5       6
            //fillMap[fill] = {energy, (double)start, (double)stop, pol, pol_err, dpdt, dpdt_err};   

        }
        file.close();
    }

    if(name_begin == "pp_"){

        ifstream file_pp("polar/dat/fill_poll_pp.csv");
        if (!file_pp.is_open()) {
            cerr << "Error opening file_pp: " << filename << endl;
            return;
        }

        string line_pp;
        getline(file_pp, line_pp); // skip header line

        while (getline(file_pp, line_pp)) {
            if (line_pp.empty()) continue;

            stringstream ss(line_pp);
            //string fill_str, pol_str, pol_err_str;
            string fill_str,energy_str,start_str,stop_str,pol_b_str,pol_b_err_str,dpdt_b_str,dpdt_b_err_str,pol_y_str,pol_y_err_str,dpdt_y_str,dpdt_y_err_str;

            getline(ss, fill_str, ',');
            getline(ss, energy_str, ',');
            getline(ss, start_str, ',');
            getline(ss, stop_str, ',');
            getline(ss, pol_b_str, ',');
            getline(ss, pol_b_err_str, ',');
            getline(ss, dpdt_b_str, ',');
            getline(ss, dpdt_b_err_str, ',' );
            getline(ss, pol_y_str, ',');
            getline(ss, pol_y_err_str, ',');
            getline(ss, dpdt_y_str, ',');
            getline(ss, dpdt_y_err_str, ',' );

            //std::cout<< fill_str << "," << pol_b_str << "," << pol_b_err_str << std::endl;
            int fill = atoi(fill_str.c_str());
            double energy = atof(energy_str.c_str());
            int start = atoi(start_str.c_str());
            int stop = atoi(stop_str.c_str());
            double pol_b = atof(pol_b_str.c_str());
            double pol_b_err = atof(pol_b_err_str.c_str());
            double dpdt_b = atof(dpdt_b_str.c_str());
            double dpdt_b_err = atof(dpdt_b_err_str.c_str());
            double pol_y = atof(pol_y_str.c_str());
            double pol_y_err = atof(pol_y_err_str.c_str());
            double dpdt_y = atof(dpdt_y_str.c_str());
            double dpdt_y_err = atof(dpdt_y_err_str.c_str());

            std::vector<double> vals;

            vals.push_back(energy);  //0
            vals.push_back(double(start)); //1
            vals.push_back(double(stop)); //2
            vals.push_back(pol_b); //3
            vals.push_back(pol_b_err); //4
            vals.push_back(dpdt_b); //5
            vals.push_back(dpdt_b_err); //6
            vals.push_back(pol_y); //7
            vals.push_back(pol_y_err); //8
            vals.push_back(dpdt_y); //9
            vals.push_back(dpdt_y_err); //10

            fillMap[fill] = vals; // assign vector to map
            //               0       1                2            3    4        5       6
            //fillMap[fill] = {energy, (double)start, (double)stop, pol, pol_err, dpdt, dpdt_err};   

        }
        file_pp.close();
    }

    TH1F *h_evt_z_rank = NULL;
    TH1F *h_evt_z_high_rank = NULL;
    h_evt_z_rank = (TH1F*) inputRootFile->Get("h_evt_z_rank")->Clone("h_evt_z_rank_c");
    h_evt_z_high_rank = (TH1F*) inputRootFile->Get("h_evt_z_high_rank")->Clone("h_evt_z_high_rank_c");

    // Define histograms
    TH1F *h_pol     = new TH1F("h_pol","h_pol",100,0,100);
    TH1F *h_pol_err = new TH1F("h_pol_err","h_pol_err",25,0,5);
    TH1F *h_w_pol     = new TH1F("h_w_pol",    "h_w_pol",100,0,100);
    TH1F *h_w_pol_err = new TH1F("h_w_pol_err","h_w_pol_err",25,0,5);

    TH1F *h_y_OS  = new TH1F( "h_y_OS", "Pair OS y for e", 40, -2, 2 ) ; 
    TH1F *h_y_SS  = new TH1F( "h_y_SS", "Pair SS y for e", 40, -2, 2 ) ; 

    TH2F *h_dca_spin_y         = new TH2F("h_dca_spin_y",     "h_dca_spin_y",      20, -3.5, 16.5, 60, -3, 3);
    //TH2F *h_dca_spin_y_L       = new TH2F("h_dca_spin_y_L",   "h_dca_spin_y_L",    20, -3.5, 16.5, 60, -3, 3);
    //TH2F *h_dca_spin_y_R       = new TH2F("h_dca_spin_y_R",   "h_dca_spin_y_R",    20, -3.5, 16.5, 60, -3, 3);
    TH2F *h_dca_spin_y_SS      = new TH2F("h_dca_spin_y_SS",  "h_dca_spin_y_SS",   20, -3.5, 16.5, 60, -3, 3);
    //TH2F *h_dca_spin_y_SS_L    = new TH2F("h_dca_spin_y_SS_L","h_dca_spin_y_SS_L", 20, -3.5, 16.5, 60, -3, 3);
    //TH2F *h_dca_spin_y_SS_R    = new TH2F("h_dca_spin_y_SS_R","h_dca_spin_y_SS_R", 20, -3.5, 16.5, 60, -3, 3);


    TH2F *h_dca_spin_mee       = new TH2F("h_dca_spin_mee",     "h_dca_spin_mee", 20, -3.5, 16.5, 40, 2,6);
    TH2F *h_dca_spin_mee_R     = new TH2F("h_dca_spin_mee_R",   "h_dca_spin_mee_R", 20, -3.5, 16.5, 40, 2,6);
    TH2F *h_dca_spin_mee_L     = new TH2F("h_dca_spin_mee_L",   "h_dca_spin_mee_L", 20, -3.5, 16.5, 40, 2,6);
    TH2F *h_dca_spin_mee_SS    = new TH2F("h_dca_spin_mee_SS",  "h_dca_spin_mee_SS", 20, -3.5, 16.5, 40, 2,6);
    TH2F *h_dca_spin_mee_SS_R  = new TH2F("h_dca_spin_mee_SS_R","h_dca_spin_mee_SS_R", 20, -3.5, 16.5, 40, 2,6);
    TH2F *h_dca_spin_mee_SS_L  = new TH2F("h_dca_spin_mee_SS_L","h_dca_spin_mee_SS_L", 20, -3.5, 16.5, 40, 2,6);

    TH2F *h_w_dca_spin_mee       = new TH2F("h_w_dca_spin_mee",     "h_w_dca_spin_mee", 20, -3.5, 16.5, 40, 2,6);
    TH2F *h_w_dca_spin_mee_R     = new TH2F("h_w_dca_spin_mee_R",   "h_w_dca_spin_mee_R", 20, -3.5, 16.5, 40, 2,6);
    TH2F *h_w_dca_spin_mee_L     = new TH2F("h_w_dca_spin_mee_L",   "h_w_dca_spin_mee_L", 20, -3.5, 16.5, 40, 2,6);
    TH2F *h_w_dca_spin_mee_SS    = new TH2F("h_w_dca_spin_mee_SS",  "h_w_dca_spin_mee_SS", 20, -3.5, 16.5, 40, 2,6);
    TH2F *h_w_dca_spin_mee_SS_R  = new TH2F("h_w_dca_spin_mee_SS_R","h_w_dca_spin_mee_SS_R", 20, -3.5, 16.5, 40, 2,6);
    TH2F *h_w_dca_spin_mee_SS_L  = new TH2F("h_w_dca_spin_mee_SS_L","h_w_dca_spin_mee_SS_L", 20, -3.5, 16.5, 40, 2,6);

    TH2F* h_dca_cosphi_sinphi_u       = new TH2F("h_dca_cosphi_sinphi_u",   "h_dca_cosphi_sinphi_u", 44, -1.1,1.1, 44, -1.1,1.1);
    TH2F* h_dca_cosphi_sinphi_d       = new TH2F("h_dca_cosphi_sinphi_d",   "h_dca_cosphi_sinphi_d", 44, -1.1,1.1, 44, -1.1,1.1);
    TH2F* h_dca_cosphi_sinphi_u_SS       = new TH2F("h_dca_cosphi_sinphi_u_SS",   "h_dca_cosphi_sinphi_u_SS", 44, -1.1,1.1, 44, -1.1,1.1);
    TH2F* h_dca_cosphi_sinphi_d_SS       = new TH2F("h_dca_cosphi_sinphi_d_SS",   "h_dca_cosphi_sinphi_d_SS", 44, -1.1,1.1, 44, -1.1,1.1);
    //Setup and read the Tree    
    int eventT_t ;
    int runN_t ;
    int fillN_t ;
    double mee_t ;
    int spin_b_t ;
    int spin_y_t ;
    double phi_t ;
    double y_t ;
    double pt_t ;
    int ss_t ;
    double mee_weight_t;

    bool trigger_cut;
    bool vertex_cut;
    bool calo_match_cut;
    bool hits_cut;
    bool dca_cut;
    bool emc_energy_cut;
    bool chiSquare_ee_cut;
    bool chiSquare_non_ee_cut;
    bool bemc_wedge_b2b_cut;

    TTree *runT = (TTree*)inputRootFile->Get("runT");

    runT->SetBranchAddress("runN",&runN_t);
    runT->SetBranchAddress("eventT",&eventT_t);
    runT->SetBranchAddress("fillN",&fillN_t);
    runT->SetBranchAddress("mee",&mee_t);
    runT->SetBranchAddress("spin_b",&spin_b_t);
    runT->SetBranchAddress("spin_y",&spin_y_t);
    runT->SetBranchAddress("phi",&phi_t);
    runT->SetBranchAddress("y",&y_t);
    runT->SetBranchAddress("pt",&pt_t);
    runT->SetBranchAddress("ss",&ss_t);
    runT->SetBranchAddress("mee_weight",&mee_weight_t);
    runT-> SetBranchAddress("trigger_cut",&trigger_cut);
    runT-> SetBranchAddress("vertex_cut",&vertex_cut);
    runT-> SetBranchAddress("calo_match_cut",&calo_match_cut);
    runT-> SetBranchAddress("hits_cut",&hits_cut);
    runT-> SetBranchAddress("dca_cut",&dca_cut);
    runT-> SetBranchAddress("bemc_wedge_b2b_cut",&bemc_wedge_b2b_cut);
    runT-> SetBranchAddress("emc_energy_cut",&emc_energy_cut);
    runT-> SetBranchAddress("chiSquare_ee_cut",&chiSquare_ee_cut);
    runT-> SetBranchAddress("chiSquare_non_ee_cut",&chiSquare_non_ee_cut);


    Long64_t nentries_t = runT->GetEntries();
    double totalPol = 0;
    double totalErr = 0;
    int count = 0;
    for (Long64_t it = 0; it < nentries_t; it++) {
        runT->GetEntry(it);
        
        //std::cout << "Run: " << runN_t << " | Fill: " << fillN_t << std::endl;
        if (fillN_t == 0) continue; // skip if fill number is 0
        //if ( fillN_t == 19134 || fillN_t == 19159 || fillN_t == 19169)   
        //    std::cout << runN_t << std::endl;     
        //if ( fillN_t != 19134 && fillN_t != 19159 ){

            double fill_start = fillMap[fillN_t][1];
            double pol = fillMap[fillN_t][3];
            double pol_err = fillMap[fillN_t][4];
            totalPol += pol;
            totalErr += pol_err; // sum of squares for error propagation
            count++;
            double dt = double(int(fill_start)-int(eventT_t))/3600; // hours since start of fill
            //std::cout << runN_t << " " << fillN_t << " " << pol << " +/- " << pol_err << " eventT_t " << eventT_t << " fill_start "<< fill_start << " dt = " << dt << " dP/dt = " << fillMap[fillN_t][5] << " P = " <<  pol*(1+fillMap[fillN_t][5]/100*dt) << std::endl;
            double evt_pol = pol*(1+fillMap[fillN_t][5]/100*dt);
            //if(!trigger_cut || !vertex_cut   || !emc_energy_cut || !chiSquare_ee_cut || !chiSquare_non_ee_cut || !calo_match_cut || !dca_cut || !bemc_wedge_b2b_cut) continue;
            if(!trigger_cut || !vertex_cut  || !hits_cut || !emc_energy_cut || !chiSquare_ee_cut || !chiSquare_non_ee_cut || !calo_match_cut || !dca_cut || !bemc_wedge_b2b_cut) continue;// || !calo_match_cut
            if(evt_pol>50){
                if(mee_t > 2.8 && mee_t < 3.2){
                    if(ss_t==1){
                        h_y_SS->Fill(y_t);
                        h_dca_spin_y_SS->Fill(spin_b_t,y_t);
                        if(spin_b_t==3){
                           h_dca_cosphi_sinphi_u_SS         ->Fill(cos(phi_t), sin(phi_t));   
                        }
                        if(spin_b_t==5){
                           h_dca_cosphi_sinphi_d_SS         ->Fill(cos(phi_t), sin(phi_t));   
                        } 
                    }else{
                        h_y_OS->Fill(y_t);
                        h_dca_spin_y->Fill(spin_b_t,y_t);
                        if(spin_b_t==3){
                           h_dca_cosphi_sinphi_u         ->Fill(cos(phi_t), sin(phi_t));   
                        }
                        if(spin_b_t==5){
                           h_dca_cosphi_sinphi_d         ->Fill(cos(phi_t), sin(phi_t));   
                        } 
                    }
                }
                if(fabs(y_t)<=2.0 ){ 
                    if(mee_t > 2.8 && mee_t < 3.2){
                        h_pol->Fill(evt_pol);
                        h_pol_err->Fill(pol_err);
                        h_w_pol->Fill(evt_pol,mee_weight_t);
                        h_w_pol_err->Fill(pol_err,mee_weight_t);
                    }

                    if(ss_t==1){
                        h_dca_spin_mee_SS  ->Fill(spin_b_t,mee_t);


                        if(mee_weight_t>0)h_w_dca_spin_mee_SS  ->Fill(spin_b_t,mee_t,mee_weight_t);
                        if(cos(phi_t)>=0){
                            h_dca_spin_mee_SS_R->Fill(spin_b_t,mee_t);
                            if(mee_weight_t>0)h_w_dca_spin_mee_SS_R->Fill(spin_b_t,mee_t,mee_weight_t);
                        }
                        else{
                            h_dca_spin_mee_SS_L->Fill(spin_b_t,mee_t);
                            if(mee_weight_t>0)h_w_dca_spin_mee_SS_L->Fill(spin_b_t,mee_t,mee_weight_t);
                        }
                    }
                    else{
                        h_dca_spin_mee     ->Fill(spin_b_t,mee_t);
                        if(mee_weight_t>0)h_w_dca_spin_mee     ->Fill(spin_b_t,mee_t,mee_weight_t);
                        if(cos(phi_t)>=0){
                            h_dca_spin_mee_R   ->Fill(spin_b_t,mee_t);
                            if(mee_weight_t>0)h_w_dca_spin_mee_R   ->Fill(spin_b_t,mee_t,mee_weight_t);
                        }
                        else{
                            h_dca_spin_mee_L   ->Fill(spin_b_t,mee_t);
                            if(mee_weight_t>0)h_w_dca_spin_mee_L   ->Fill(spin_b_t,mee_t,mee_weight_t);
                        }
                    }
                }
            }else{
                std::cout << runN_t << " " << fillN_t << " " << pol << " +/- " << pol_err << " eventT_t " << eventT_t << " fill_start "<< fill_start << " dt = " << dt << " dP/dt = " << fillMap[fillN_t][5] << " P = " <<  pol*(1+fillMap[fillN_t][5]/100*dt) << std::endl;
            }

        //}
    }

    double avgPol = totalPol / count;
    double avgErr = totalErr / count; // standard error propagation
    std::cout << "Average Polarisation: " << avgPol << " ± " << avgErr << std::endl;
    TString Mark = "empty";
    if(name_begin == "pA_"){
        Mark = "p#uparrowAu#rightarrowe^{#plus}e^{#minus}pAu #sqrt{s_{pN}}=200GeV";
    }
    if(name_begin == "pp_"){
        Mark = "p#uparrowp#downarrow#rightarrowe^{#plus}e^{#minus}pp #sqrt{s_{pp}}=500GeV";
    }
    h_evt_z_rank ->Rebin(10);
    h_evt_z_high_rank ->Rebin(10);

    if(name_begin=="pA_") {
        draw_1d_z( h_evt_z_rank, name_begin + "z_rank", "z", "N", Mark, " ", -199, 199, 500, 10*1399, false);
    }else{
        draw_1d_z( h_evt_z_rank, name_begin + "z_rank", "z", "N", Mark, " ", -199, 199, 500, 1000*1899, false);
    }
    //draw_1d( h_evt_z_high_rank, name_begin + "z_high_rank", "z", "N", Mark, " ", -199, 199, 500, 10*1399, true);
    /*
    // Draw histograms
    char name[100];
    sprintf(name,"2.8<M_{ee}<3.2 GeV/c^{2}; Mean: %.3f", P_avg);

    if(name_begin=="pA_") {
        draw_1d( h_pol, name_begin + "pol", "P_{0} [%]", "N", Mark, name, 51, 69, 0, 99);
    }else{
        draw_1d( h_pol, name_begin + "pol", "P_{0} [%]", "N", Mark, name, 40, 80, 0, 899);
    }
    sprintf(name,"2.8<M_{ee}<3.2 GeV/c^{2}; Mean: %.3f", h_y_SS->GetMean());
    draw_1d(h_y_SS , name_begin + "y_SS", "y", "N", Mark, name, y_min, y_max, 0, 490);
    sprintf(name,"2.8<M_{ee}<3.2 GeV/c^{2}; Mean: %.3f", h_y_OS->GetMean());
    draw_1d(h_y_OS , name_begin + "y_OS", "y", "N", Mark, name, y_min, y_max, 0, 490);

    h_y_OS->Add(h_y_SS,-1); // subtract SS from OS to get signal distribution
    sprintf(name,"2.8<M_{ee}<3.2 GeV/c^{2}; Mean: %.3f", h_y_OS->GetMean());
    draw_1d(h_y_OS , name_begin + "y_OS_noSS", "y", "N", Mark, name, y_min, y_max, 0, 490);

    sprintf(name,"Mean: %.3f", h_pol_err->GetMean());    
    if(name_begin=="pA_") {
        draw_1d( h_pol_err, name_begin + "pol_err", "P_{err} [%]", "N", Mark, name, 1, 4);
    }else{
        draw_1d( h_pol_err, name_begin + "pol_err", "P_{err} [%]", "N", Mark, name, 1, 4, 0, 1899);
    }

    h_dca_cosphi_sinphi_u->Add(h_dca_cosphi_sinphi_u_SS,-1);
    h_dca_cosphi_sinphi_d->Add(h_dca_cosphi_sinphi_d_SS,-1);
    if(name_begin=="pA_") {
        draw_updown_projections(h_dca_cosphi_sinphi_u, h_dca_cosphi_sinphi_d, name_begin+"dca_cosphi_sinphi_udprojection", "cos#phi" , "#", Mark , " 2.8<M_{ee}<3.2 GeV/c^{2}", -1.2, 1.2, 0, 30);
    }else{
        draw_updown_projections(h_dca_cosphi_sinphi_u, h_dca_cosphi_sinphi_d, name_begin+"dca_cosphi_sinphi_udprojection", "cos#phi" , "#", Mark , " 2.8<M_{ee}<3.2 GeV/c^{2}", -1.2, 1.2, 0, 200);
    }


    draw_projections(h_dca_spin_mee,   h_dca_spin_mee_SS,   name_begin+"RP2E_dca_spin_mee_projection",   "M_{ee} [GeV/c^{2}]" , "#", Mark);
    //draw_projections(h_dca_spin_mee_R, h_dca_spin_mee_SS_R, name_begin+"RP2E_dca_spin_mee_R_projection", "M_{ee} [GeV/c^{2}]" , "#", Mark, "#phi>0");
    //draw_projections(h_dca_spin_mee_L, h_dca_spin_mee_SS_L, name_begin+"RP2E_dca_spin_mee_L_projection", "M_{ee} [GeV/c^{2}]" , "#", Mark, "#phi<0");

    fit_background_with_gap( h_dca_spin_mee  , name_begin + "fit_dca_spin_mee_projection", "M_{ee} [GeV/c^{2}]" , "#frac{Entries}{0.1 [GeV/c^{2}]}", Mark, "empty", mee_min, mee_max); 
    //MC_background( h_dca_spin_mee  , h_dca_spin_mee_SS  , name_begin + "MC_background", "M_{ee} [GeV/c^{2}]" , "#frac{Entries}{0.1 [GeV/c^{2}]}", Mark, "empty", mee_min, mee_max); 

    draw_tot_projections(h_dca_spin_mee,   h_dca_spin_mee_SS,   name_begin+"dca_mee_projection",   "M_{ee} [GeV/c^{2}]" , "#frac{Entries}{0.1 [GeV/c^{2}]}", Mark, "empty", mee_min, mee_max);
    draw_tot_projections(h_dca_spin_y,     h_dca_spin_y_SS,     name_begin+"dca_y_projection",     "y" ,                  "#frac{Entries}{0.1 y}" ,          Mark , " 2.8<M_{ee}<3.2 GeV/c^{2}", y_min, y_max);

    y_Mean = h_y_OS->GetMean();
    P_avg = h_pol->GetMean()/100.0; // convert from percentage to fraction
    P_avg_err = h_pol_err->GetMean()/100.0; // convert from percentage to fraction
    //sqrt_cosPhiLR = 0.630;
    correction = 1.0/(P_avg*sqrt_cosPhiLR);

    if(name_begin=="pA_") {
        raw_tot_corr_asymmetry_m(h_dca_spin_mee, h_dca_spin_mee_SS, h_dca_spin_mee_L, h_dca_spin_mee_R,h_dca_spin_mee_SS_L, h_dca_spin_mee_SS_R,  name_begin+"dca_tot_corr_assym_mee", "W_{#gammap} [GeV]", "A_{N}", Mark, "empty", name_begin);
    }else{
        raw_tot_corr_asymmetry_m(h_dca_spin_mee, h_dca_spin_mee_SS, h_dca_spin_mee_L, h_dca_spin_mee_R,h_dca_spin_mee_SS_L, h_dca_spin_mee_SS_R,  name_begin+"dca_tot_corr_assym_mee", "W_{#gammap} [GeV]", "A_{N}", Mark, "empty", name_begin, false);
    }
    std::cout << "P_avg: " << P_avg << " ± " << P_avg_err << std::endl;
    std::cout << "sqrt_cosPhiLR: " << sqrt_cosPhiLR << std::endl;
    std::cout << "Correction factor (1/(P_avg*sqrt_cosPhiLR)): " << correction << std::endl;
    */
}


// pp variations
//All cuts:
// A:-0.00439278 +/- 0.0504165 | pp_dca_tot_corr_assym_mee
// A fit:0.031257 +/- 0.0531283 | pp_dca_tot_corr_assym_mee
// P_avg: 0.598085 ± 0.0238452
// sqrt_cosPhiLR: 0.632134
// Correction factor (1/(P_avg*sqrt_cosPhiLR)): 2.64501

// no dca cut:
// A:0.00374624 +/- 0.0481091 | pp_dca_tot_corr_assym_mee
// A fit:0.0723306 +/- 0.0506604 | pp_dca_tot_corr_assym_mee
// P_avg: 0.597921 ± 0.0238702
// sqrt_cosPhiLR: 0.631697
// Correction factor (1/(P_avg*sqrt_cosPhiLR)): 2.64757

// no calo match cut:
// A:-0.00439278 +/- 0.0504165 | pp_dca_tot_corr_assym_mee
// A fit:0.031257 +/- 0.0531283 | pp_dca_tot_corr_assym_mee
// P_avg: 0.598085 ± 0.0238452
// sqrt_cosPhiLR: 0.632134
// Correction factor (1/(P_avg*sqrt_cosPhiLR)): 2.64501

// soft hits cut:
// A:-0.0228529 +/- 0.0500311 | pp_dca_tot_corr_assym_mee
// A fit:0.121831 +/- 0.0536334 | pp_dca_tot_corr_assym_mee
// P_avg: 0.59807 ± 0.0238529
// sqrt_cosPhiLR: 0.631156
// Correction factor (1/(P_avg*sqrt_cosPhiLR)): 2.64918

// soft emc energy cut:
// A:-0.0141156 +/- 0.0499129 | pp_dca_tot_corr_assym_mee
// A fit:0.0118565 +/- 0.0528778 | pp_dca_tot_corr_assym_mee
// P_avg: 0.598064 ± 0.0238409
// sqrt_cosPhiLR: 0.632668
// Correction factor (1/(P_avg*sqrt_cosPhiLR)): 2.64287

// very soft cuts:
// A:-0.0219874 +/- 0.0455354 | pp_dca_tot_corr_assym_mee
// A fit:-0.130366 +/- 0.0587361 | pp_dca_tot_corr_assym_mee
// P_avg: 0.597356 ± 0.0239047
// sqrt_cosPhiLR: 0.62942
// Correction factor (1/(P_avg*sqrt_cosPhiLR)): 2.65966