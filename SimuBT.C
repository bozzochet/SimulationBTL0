#include "TCanvas.h"
#include "TH1F.h"
#include "TF1.h"
#include "TRandom.h"
#include "TString.h"
#include "TGraphErrors.h"

double resoX = 0.01; //10 um

double SingleTrial(int events=1000000, bool kNoDraw=false);

void Single(){
  SingleTrial();
}

void Full(){
  int trials=1000;

  TH1* resoXmeas = new TH1F("resoXmeas", "reasoXmeas", 100, -0.02, 0.02);
  resoXmeas->GetXaxis()->SetTitle("(#sigma_x measured - #sigma_x true) / #sigma_x true");
  resoXmeas->GetYaxis()->SetTitle("Trials");
  
  for (int ii=0; ii<trials; ii++) {
    double meas = 0.001*SingleTrial(1000000, true);
    resoXmeas->Fill((meas-resoX)/resoX);
    printf("%d) (%f -%f)/%f = %f\n", ii, meas, resoX, resoX, (meas-resoX)/resoX);
  }

  TCanvas* c = new TCanvas("RelativeErrorOnResolution", "RelativeErrorOnResolution");
  resoXmeas->Draw();
}

double SingleTrial(int events, bool kNoDraw) {

  // do not use the downstream telescope
  bool kNoDowTel = false;

  // half distances between Up - DUT - Down
  bool kShorten = false;
  
  //------------
  
  double sigmaBeamPipeX = 2; // 2 mm
  double sigmaThetaX = 1.0;// 1 mrad

  double BeamPipeZ = 0.0;

  const int nUpTel = 6;
  double UpTelZ[nUpTel] = {1000.0, 1000.1, 1030.0, 1030.1, 1060.0, 1060.1};
  TString UpTelK[nUpTel] = {"Y", "X", "X", "Y", "Y", "X"};
  
  double DUTZ = 1500.0;
  if (kShorten) {
    DUTZ = 1250.0;
  }
  if (kNoDowTel) {
    DUTZ = UpTelZ[nUpTel-1] + 20.0;
  }

  const int nDowTel = 4;
  double DowTelZ[nDowTel] = {2000.0, 2000.1, 2060.0, 2060.1};
  if (kShorten) {
    double DowTelZShorter[nDowTel] = {1500.0, 1500.1, 1560.0, 1560.1};
    for (int ii=0; ii<nDowTel; ii++) {
      DowTelZ[ii] = DowTelZShorter[ii];
    }
  }
  TString DowTelK[nDowTel] = {"X", "Y", "X", "Y"};

  //------------
  
  TF1* mct = new TF1("mct", "[0] + [1]*x", 0.0, 2.0);
  TF1* fit = new TF1("fit", "[0] + [1]*x", 0.0, 2.0);

  TF1* smegaus = new TF1("smegaus", "[0]*exp(-0.5*((x-[1])/[2])*((x-[1])/[2]))", -50, 50);
  smegaus->SetLineColor(kBlack);
  smegaus->SetLineWidth(3);
  TF1* resgaus = new TF1("resgaus", "[0]*exp(-0.5*((x-[1])/[2])*((x-[1])/[2]))", -50, 50);
  resgaus->SetLineColor(kRed+2);
  resgaus->SetLineWidth(3);
  
  TH1* mctBeamPipeX = new TH1F("mctBeamPipeX", "mctBeamPipeX", 1000, -10.0, 10.0);
  mctBeamPipeX->GetXaxis()->SetTitle("X-coo at Beam Pipe (mm)");
  mctBeamPipeX->GetYaxis()->SetTitle("Entries");
  TH1* mctThetaX = new TH1F("mctThetaX", "mctThetaX", 1000, -5.0, 5.0);
  mctThetaX->GetXaxis()->SetTitle("Theta (mrad)");
  mctThetaX->GetYaxis()->SetTitle("Entries");
  
  TH1* mctUpTelX[nUpTel];
  for (int ii=0; ii<nUpTel; ii++) {
    mctUpTelX[ii] = NULL;
    if (UpTelK[ii] == "X") {
      mctUpTelX[ii] = new TH1F(Form("mctUpTelX_%d", ii), Form("mctUpTelX_%d", ii), 1000, -10.0, 10.0);
      mctUpTelX[ii]->GetXaxis()->SetTitle(Form("X-coo at Upstream (%d) (mm)", ii));
      mctUpTelX[ii]->GetYaxis()->SetTitle("Entries");
    }
  }

  TH1* mctDUTX = new TH1F("mctDUTX", "mctDUTX", 1000, -10.0, 10.0);
  mctDUTX->GetXaxis()->SetTitle("X-coo at DUT (mm)");
  mctDUTX->GetYaxis()->SetTitle("Entries");
  
  TH1* mctDowTelX[nDowTel];
  for (int ii=0; ii<nDowTel; ii++) {
    mctDowTelX[ii] = NULL;
    if (!kNoDowTel) {
      if (DowTelK[ii] == "X") {
	mctDowTelX[ii] = new TH1F(Form("mctDowTelX_%d", ii), Form("mctDowTelX_%d", ii), 1000, -10.0, 10.0);
	mctDowTelX[ii]->GetXaxis()->SetTitle(Form("X-coo at Downstream (%d) (mm)", ii));
	mctDowTelX[ii]->GetYaxis()->SetTitle("Entries");
      }
    }
  }
  
  TH1* fitBeamPipeX = new TH1F("fitBeamPipeX", "fitBeamPipeX", 1000, -10.0, 10.0);
  fitBeamPipeX->SetLineColor(kRed+2);
  fitBeamPipeX->SetMarkerColor(kRed+2);
  TH1* fitThetaX = new TH1F("fitThetaX", "fitThetaX", 1000, -5.0, 5.0);  
  fitThetaX->SetLineColor(kRed+2);
  fitThetaX->SetMarkerColor(kRed+2);
  
  TH1* fitDUTX = new TH1F("fitDUTX", "fitDUTX", 1000, -10.0, 10.0);
  fitDUTX->SetLineColor(kRed+2);
  fitDUTX->SetMarkerColor(kRed+2);

  double reslim = 50.0;
  if (kNoDowTel) {
    reslim = 100.0;
  }
  
  TH1* fitDUTSmeX = new TH1F("fitDUTSmeX", "fitDUTSmeX", 1000, -reslim, reslim);
  fitDUTSmeX->GetXaxis()->SetTitle("X-coo residual at DUT (#mum)");
  fitDUTSmeX->GetYaxis()->SetTitle("Entries");
  TH1* fitDUTResX = new TH1F("fitDUTResX", "fitDUTResX", 1000, -reslim, reslim);
  fitDUTResX->GetXaxis()->SetTitle("X-coo residual at DUT (#mum)");
  fitDUTResX->GetYaxis()->SetTitle("Entries");
  fitDUTResX->SetLineColor(kRed+2);
  fitDUTResX->SetMarkerColor(kRed+2);
  
  for (int ii=0; ii<events; ii++) {
    
    double beampipex = gRandom->Gaus(0, sigmaBeamPipeX);
    double thetax = gRandom->Gaus(0, sigmaThetaX);

    mct->SetParameter(0, beampipex);
    mct->SetParameter(1, TMath::Tan(0.001*thetax));

    mctBeamPipeX->Fill(mct->Eval(BeamPipeZ));
    mctThetaX->Fill(thetax);

    for (int ii=0; ii<nUpTel; ii++) {
      if (UpTelK[ii] == "X") {
	mctUpTelX[ii]->Fill(mct->Eval(UpTelZ[ii]));
      }
    }

    mctDUTX->Fill(mct->Eval(DUTZ));

    if (!kNoDowTel) {
      for (int ii=0; ii<nDowTel; ii++) {
	if (DowTelK[ii] == "X") {
	  mctDowTelX[ii]->Fill(mct->Eval(DowTelZ[ii]));
	}
      }
    }

    TGraphErrors gr;

    for (int ii=0; ii<nUpTel; ii++) {
      if (UpTelK[ii] == "X") {
	double ut = gRandom->Gaus(mct->Eval(UpTelZ[ii]), resoX);
	gr.SetPoint(gr.GetN(), UpTelZ[ii], ut);
	gr.SetPointError(gr.GetN()-1, 0.0, resoX);
      }
    }

    if (!kNoDowTel) {
      for (int ii=0; ii<nDowTel; ii++) {
	if (DowTelK[ii] == "X") {
	  double dt = gRandom->Gaus(mct->Eval(DowTelZ[ii]), resoX);
	  gr.SetPoint(gr.GetN(), DowTelZ[ii], dt);
	  gr.SetPointError(gr.GetN()-1, 0.0, resoX);
	}
      }
    }

    gr.Fit(fit, "NQ");

    fitBeamPipeX->Fill(fit->Eval(BeamPipeZ));
    fitThetaX->Fill(TMath::ATan(1000.0*fit->GetParameter(1)));

    fitDUTX->Fill(fit->Eval(DUTZ));

    fitDUTSmeX->Fill(1000.0*(fit->Eval(DUTZ) - mct->Eval(DUTZ)));
    fitDUTResX->Fill(1000.0*(gRandom->Gaus(fit->Eval(DUTZ), resoX) - mct->Eval(DUTZ)));
  }
  
  if (!kNoDraw) {
    TCanvas* cbpx = new TCanvas("BeamPipeX", "BeamPipeX");
    mctBeamPipeX->Draw();
    fitBeamPipeX->Draw("same");
    
    TCanvas* cstx = new TCanvas("ThetaX", "ThetaX");
    mctThetaX->Draw();
    fitThetaX->Draw("same");
    
    TCanvas* cutx[nUpTel];
    for (int ii=0; ii<nUpTel; ii++) {
      if (UpTelK[ii] == "X") {
	cutx[ii] = new TCanvas(Form("UpstreamTelescopeX_%d", ii), Form("UpstreamTelescopeX_%d", ii));
	mctUpTelX[ii]->Draw();
      }
    }
    
    TCanvas* cdutx = new TCanvas("DUTX", "DUTX");
    mctDUTX->Draw();
    fitDUTX->Draw("same");
    
    TCanvas* cdtx[nDowTel];
    if (!kNoDowTel) {
      for (int ii=0; ii<nDowTel; ii++) {
	if (DowTelK[ii] == "X") {
	  cdtx[ii] = new TCanvas(Form("DownstreamTelescopeX_%d", ii), Form("DownstreamTelescopeX_%d", ii));
	  mctDowTelX[ii]->Draw();
	}
      }
    }
  }

  if (!kNoDraw) {
    TCanvas* cdsx = new TCanvas("DUTSmearingX", "DUTSmearingX");
    fitDUTSmeX->Draw();
  }
  smegaus->SetParameter(0, fitDUTSmeX->GetMaximum());
  smegaus->SetParameter(1, 0.0);
  smegaus->SetParLimits(2, 0.0, 10000.0*resoX);
  smegaus->SetParameter(2, 1000.0*resoX);
  if (kNoDraw) {
    fitDUTSmeX->Fit(smegaus, "NQ", "L");
  }
  else {
    fitDUTSmeX->Fit(smegaus, "Q", "L");
  }
  
  if (!kNoDraw) {
    TCanvas* cdrx = new TCanvas("DUTResidualX", "DUTResidualX");
    fitDUTResX->GetYaxis()->SetRangeUser(0, 2.2*fitDUTResX->GetMaximum());
    fitDUTResX->Draw();
  }
  resgaus->SetParameter(0, fitDUTResX->GetMaximum());
  resgaus->SetParameter(1, 0.0);
  resgaus->SetParameter(2, 1000.0*resoX);
  if (kNoDraw) {
    fitDUTResX->Fit(resgaus, "NQ", "L");
  }
  else {
    fitDUTResX->Fit(resgaus, "Q", "L");
    fitDUTSmeX->Draw("same");
  }
      
  //--------------
  
  double sigmaXmeas = resgaus->GetParameter(2);
  double smeaXmeas = smegaus->GetParameter(2);
  
  double resoXmeas = TMath::Sqrt(TMath::Power(sigmaXmeas, 2.0) - TMath::Power(smeaXmeas, 2.0));
  if (!kNoDraw) {
    printf("DUT X) sigma measured: %f, smearing measured: %f\n", sigmaXmeas, smeaXmeas);
    printf("DUT X) resolution: %f\n", resoXmeas);
  }

  //--------------

  if (kNoDraw){
    if (mct) delete mct;
    if (fit) delete fit;
    
    if (smegaus) delete smegaus;
    if (resgaus) delete resgaus;
    
    if (mctBeamPipeX) delete mctBeamPipeX;
    if (mctThetaX) delete mctThetaX;
    
    for (int ii=0; ii<nUpTel; ii++) {
      if (mctUpTelX[ii]) delete mctUpTelX[ii];
    }
    
    if (mctDUTX) delete mctDUTX;
    
    for (int ii=0; ii<nDowTel; ii++) {
      if (mctDowTelX[ii]) delete mctDowTelX[ii];
    }
    
    if (fitBeamPipeX) delete fitBeamPipeX;
    if (fitThetaX) delete fitThetaX;
    
    if (fitDUTX) delete fitDUTX;
    
    if (fitDUTSmeX) delete fitDUTSmeX;
    if (fitDUTResX) delete fitDUTResX;
  }

  //--------------
  
  return resoXmeas;
}
