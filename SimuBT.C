#include "TCanvas.h"
#include "TH1F.h"
#include "TF1.h"
#include "TRandom.h"
#include "TString.h"
#include "TGraphErrors.h"

double resoX = 0.01; //10 um
double resoY = 0.01; //10 um
double resoU = 0.01; //10 um

double SingleTrial(int events=1000000, bool kNoDraw=false, int component=0);

void Single(){
  SingleTrial();
}

void Full(){
  int trials=100;

  TH1* resoXmeas = new TH1F("resoXmeas", "reasoXmeas", 100, -0.02, 0.02);
  resoXmeas->GetXaxis()->SetTitle("(#sigma_x measured - #sigma_x true) / #sigma_x true");
  resoXmeas->GetYaxis()->SetTitle("Trials");
  
  for (int ii=0; ii<trials; ii++) {
    double meas = 0.001*SingleTrial(100000, true);
    resoXmeas->Fill((meas-resoX)/resoX);
    printf("%d) (%f -%f)/%f = %f\n", ii, meas, resoX, resoX, (meas-resoX)/resoX);
  }

  TCanvas* c = new TCanvas("RelativeErrorOnResolution", "RelativeErrorOnResolution");
  resoXmeas->Draw();
}

double SingleTrial(int events, bool kNoDraw, int component) {

  // do not use the downstream telescope
  bool kNoDowTel = false;

  // half distances between Up - DUT - Down
  bool kShorten = false;
  
  //------------
  
  double sigmaBeamPipeX = 2; // 2 mm
  double sigmaBeamPipeY = 2; // 2 mm
  double sigmaThetaX = 1.0;// 1 mrad
  double sigmaThetaY = 1.0;// 1 mrad

  double BeamPipeZ = 0.0;

  const int nUpTel = 6;
  double UpTelZ[nUpTel] = {1000.0, 1000.1, 1030.0, 1030.1, 1060.0, 1060.1};
  TString UpTelK[nUpTel] = {"Y", "X", "X", "Y", "Y", "X"};

  const int nDUT = 2;
  double DUTZ[nDUT] = {1500.0, 1500.1};
  TString DUTK[nDUT] = {"X", "U"};
  double rotU = 45;//degrees
  
  if (kShorten) {
    DUTZ[0] = 1250.0;
    DUTZ[1] = 1250.1;
  }
  if (kNoDowTel) {
    DUTZ[0] = UpTelZ[nUpTel-1] + 20.0;
    DUTZ[1] = UpTelZ[nUpTel-1] + 20.1;
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
  
  TF1* mctX = new TF1("mctX", "[0] + [1]*x", 0.0, 2.0);
  TF1* mctY = new TF1("mctY", "[0] + [1]*x", 0.0, 2.0);
  TF1* fitX = new TF1("fitX", "[0] + [1]*x", 0.0, 2.0);
  TF1* fitY = new TF1("fitY", "[0] + [1]*x", 0.0, 2.0);

  TF1* smegaus = new TF1("smegaus", "[0]*exp(-0.5*((x-[1])/[2])*((x-[1])/[2]))", -50, 50);
  smegaus->SetLineColor(kBlack);
  smegaus->SetLineWidth(3);
  TF1* resgaus = new TF1("resgaus", "[0]*exp(-0.5*((x-[1])/[2])*((x-[1])/[2]))", -50, 50);
  resgaus->SetLineColor(kRed+2);
  resgaus->SetLineWidth(3);
  
  TH1* mctBeamPipeX = new TH1F("mctBeamPipeX", "mctBeamPipeX", 1000, -10.0, 10.0);
  mctBeamPipeX->GetXaxis()->SetTitle("X-coo at Beam Pipe (mm)");
  mctBeamPipeX->GetYaxis()->SetTitle("Entries");
  TH1* mctBeamPipeY = new TH1F("mctBeamPipeY", "mctBeamPipeY", 1000, -10.0, 10.0);
  mctBeamPipeY->GetXaxis()->SetTitle("Y-coo at Beam Pipe (mm)");
  mctBeamPipeY->GetYaxis()->SetTitle("Entries");
  TH1* mctThetaX = new TH1F("mctThetaX", "mctThetaX", 1000, -5.0, 5.0);
  mctThetaX->GetXaxis()->SetTitle("Theta (mrad)");
  mctThetaX->GetYaxis()->SetTitle("Entries");
  TH1* mctThetaY = new TH1F("mctThetaY", "mctThetaY", 1000, -5.0, 5.0);
  mctThetaY->GetXaxis()->SetTitle("Theta (mrad)");
  mctThetaY->GetYaxis()->SetTitle("Entries");
  
  TH1* mctUpTelX[nUpTel];
  for (int ii=0; ii<nUpTel; ii++) {
    mctUpTelX[ii] = NULL;
    if (UpTelK[ii] == "X") {
      mctUpTelX[ii] = new TH1F(Form("mctUpTelX_%d", ii), Form("mctUpTelX_%d", ii), 1000, -10.0, 10.0);
      mctUpTelX[ii]->GetXaxis()->SetTitle(Form("X-coo at Upstream (%d) (mm)", ii));
      mctUpTelX[ii]->GetYaxis()->SetTitle("Entries");
    }
  }

  TH1* mctDUT[nDUT];
  for (int ii=0; ii<nDUT; ii++) {
    mctDUT[ii] = NULL;
    if (DUTK[ii] == "X") {
      mctDUT[ii] = new TH1F(Form("mctDUTX_%d", ii), Form("mctDUTX_%d", ii), 1000, -10.0, 10.0);
      mctDUT[ii]->GetXaxis()->SetTitle("X-coo at DUT (mm)");
      mctDUT[ii]->GetYaxis()->SetTitle("Entries");
    }
    else if (DUTK[ii] == "U") {
      mctDUT[ii] = new TH1F(Form("mctDUTU_%d", ii), Form("mctDUTU_%d", ii), 1000, -10.0, 10.0);
      mctDUT[ii]->GetXaxis()->SetTitle("U-coo at DUT (mm)");
      mctDUT[ii]->GetYaxis()->SetTitle("Entries");
    }
  }
  
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
  TH1* fitBeamPipeY = new TH1F("fitBeamPipeY", "fitBeamPipeY", 1000, -10.0, 10.0);
  fitBeamPipeY->SetLineColor(kRed+2);
  fitBeamPipeY->SetMarkerColor(kRed+2);
  TH1* fitThetaX = new TH1F("fitThetaX", "fitThetaX", 1000, -5.0, 5.0);  
  fitThetaX->SetLineColor(kRed+2);
  fitThetaX->SetMarkerColor(kRed+2);
  TH1* fitThetaY = new TH1F("fitThetaY", "fitThetaY", 1000, -5.0, 5.0);  
  fitThetaY->SetLineColor(kRed+2);
  fitThetaY->SetMarkerColor(kRed+2);
  
  double reslim = 50.0;
  if (kNoDowTel) {
    reslim = 100.0;
  }
  
  TH1* fitDUT[nDUT];
  TH1* fitDUTSmeX[nDUT];
  TH1* fitDUTResX[nDUT];
  for (int ii=0; ii<nDUT; ii++) {
    fitDUT[ii] = NULL;
    if (DUTK[ii] == "X") {
      fitDUT[ii] = new TH1F(Form("fitDUTX_%d", ii), Form("fitDUTX_%d", ii), 1000, -10.0, 10.0);
      fitDUT[ii]->SetLineColor(kRed+2);
      fitDUT[ii]->SetMarkerColor(kRed+2);
    }
    else if (DUTK[ii] == "U") {
      fitDUT[ii] = new TH1F(Form("fitDUTU_%d", ii), Form("fitDUTU_%d", ii), 1000, -10.0, 10.0);
      fitDUT[ii]->SetLineColor(kRed+2);
      fitDUT[ii]->SetMarkerColor(kRed+2);
    }
    fitDUTSmeX[ii] = new TH1F(Form("fitDUTSmeX_%d", ii), Form("fitDUTSmeX_%d", ii), 1000, -reslim, reslim);
    fitDUTSmeX[ii]->GetXaxis()->SetTitle("X-coo residual at DUT (#mum)");
    fitDUTSmeX[ii]->GetYaxis()->SetTitle("Entries");
    fitDUTResX[ii] = new TH1F(Form("fitDUTResX_%d", ii), Form("fitDUTResX_%d", ii), 1000, -reslim, reslim);
    fitDUTResX[ii]->GetXaxis()->SetTitle("X-coo residual at DUT (#mum)");
    fitDUTResX[ii]->GetYaxis()->SetTitle("Entries");
    fitDUTResX[ii]->SetLineColor(kRed+2);
    fitDUTResX[ii]->SetMarkerColor(kRed+2);
  }
    
  for (int ii=0; ii<events; ii++) {
    
    double beampipex = gRandom->Gaus(0, sigmaBeamPipeX);
    double beampipey = gRandom->Gaus(0, sigmaBeamPipeY);
    double thetax = gRandom->Gaus(0, sigmaThetaX);
    double thetay = gRandom->Gaus(0, sigmaThetaY);
    
    mctX->SetParameter(0, beampipex);
    mctX->SetParameter(1, TMath::Tan(0.001*thetax));

    mctY->SetParameter(0, beampipey);
    mctY->SetParameter(1, TMath::Tan(0.001*thetay));
    
    mctBeamPipeX->Fill(mctX->Eval(BeamPipeZ));
    mctBeamPipeY->Fill(mctY->Eval(BeamPipeZ));
    mctThetaX->Fill(thetax);
    mctThetaY->Fill(thetay);

    for (int ii=0; ii<nUpTel; ii++) {
      if (UpTelK[ii] == "X") {
	mctUpTelX[ii]->Fill(mctX->Eval(UpTelZ[ii]));
      }
    }

    for (int ii=0; ii<nDUT; ii++) {
      if (DUTK[ii] == "X") {
	mctDUT[ii]->Fill(mctX->Eval(DUTZ[ii]));
      }
      else if (DUTK[ii] == "U") {
	double x = mctX->Eval(DUTZ[ii]);
	double y = mctY->Eval(DUTZ[ii]);
	double u = x*TMath::Cos(TMath::DegToRad()*rotU) - y*TMath::Sin(TMath::DegToRad()*rotU);
	mctDUT[ii]->Fill(u);
      }
    }

    if (!kNoDowTel) {
      for (int ii=0; ii<nDowTel; ii++) {
	if (DowTelK[ii] == "X") {
	  mctDowTelX[ii]->Fill(mctX->Eval(DowTelZ[ii]));
	}
      }
    }

    TGraphErrors grX;
    TGraphErrors grY;

    for (int ii=0; ii<nUpTel; ii++) {
      if (UpTelK[ii] == "X") {
	double ut = gRandom->Gaus(mctX->Eval(UpTelZ[ii]), resoX);
	grX.SetPoint(grX.GetN(), UpTelZ[ii], ut);
	grX.SetPointError(grX.GetN()-1, 0.0, resoX);
      }
      if (UpTelK[ii] == "Y") {
	double ut = gRandom->Gaus(mctY->Eval(UpTelZ[ii]), resoY);
	grY.SetPoint(grY.GetN(), UpTelZ[ii], ut);
	grY.SetPointError(grY.GetN()-1, 0.0, resoY);
      }
    }

    if (!kNoDowTel) {
      for (int ii=0; ii<nDowTel; ii++) {
	if (DowTelK[ii] == "X") {
	  double dt = gRandom->Gaus(mctX->Eval(DowTelZ[ii]), resoX);
	  grX.SetPoint(grX.GetN(), DowTelZ[ii], dt);
	  grX.SetPointError(grX.GetN()-1, 0.0, resoX);
	}
	else if (DowTelK[ii] == "Y") {
	  double dt = gRandom->Gaus(mctY->Eval(DowTelZ[ii]), resoY);
	  grY.SetPoint(grY.GetN(), DowTelZ[ii], dt);
	  grY.SetPointError(grY.GetN()-1, 0.0, resoY);
	}
      }
    }

    grX.Fit(fitX, "NQ");
    grY.Fit(fitY, "NQ");

    fitBeamPipeX->Fill(fitX->Eval(BeamPipeZ));
    fitBeamPipeY->Fill(fitY->Eval(BeamPipeZ));
    fitThetaX->Fill(TMath::ATan(1000.0*fitX->GetParameter(1)));
    fitThetaY->Fill(TMath::ATan(1000.0*fitY->GetParameter(1)));

    for (int ii=0; ii<nDUT; ii++) {
      if (DUTK[ii] == "X") {
	double xtruth = mctX->Eval(DUTZ[ii]);
	double xtrack = fitX->Eval(DUTZ[ii]);
	double xmeas = gRandom->Gaus(xtruth, resoX);
	fitDUT[ii]->Fill(xtrack);
	fitDUTSmeX[ii]->Fill(1000.0*(xtrack - xtruth));
	fitDUTResX[ii]->Fill(1000.0*(xmeas - xtrack));
      }
      else if (DUTK[ii] == "U") {
	double xtruth = mctX->Eval(DUTZ[ii]);
	double xtrack = fitX->Eval(DUTZ[ii]);
	double ytruth = mctY->Eval(DUTZ[ii]);
	double ytrack = fitY->Eval(DUTZ[ii]);
	double utruth = xtruth*TMath::Cos(TMath::DegToRad()*rotU) -
	  ytruth*TMath::Sin(TMath::DegToRad()*rotU);
	double utrack = xtrack*TMath::Cos(TMath::DegToRad()*rotU) -
	  ytrack*TMath::Sin(TMath::DegToRad()*rotU);
	double umeas = gRandom->Gaus(utruth, resoU);
	double xmeas = gRandom->Gaus(xtruth, resoX);
	fitDUT[ii]->Fill(utrack);
	/* fitDUTSmeX[ii]->Fill(1000.0*(xtrack - xtruth)); */
	/* fitDUTResX[ii]->Fill(1000.0*(xmeas - xtrack)); */
      }
    }

  }
  
  if (!kNoDraw) {
    TCanvas* cbpx = new TCanvas("BeamPipeX", "BeamPipeX");
    mctBeamPipeX->Draw();
    fitBeamPipeX->Draw("same");

    TCanvas* cbpy = new TCanvas("BeamPipeY", "BeamPipeY");
    mctBeamPipeY->Draw();
    fitBeamPipeY->Draw("same");
    
    TCanvas* cstx = new TCanvas("ThetaX", "ThetaX");
    mctThetaX->Draw();
    fitThetaX->Draw("same");
    
    TCanvas* csty = new TCanvas("ThetaY", "ThetaY");
    mctThetaY->Draw();
    fitThetaY->Draw("same");
    
    TCanvas* cutx[nUpTel];
    for (int ii=0; ii<nUpTel; ii++) {
      if (UpTelK[ii] == "X") {
	cutx[ii] = new TCanvas(Form("UpstreamTelescopeX_%d", ii), Form("UpstreamTelescopeX_%d", ii));
	mctUpTelX[ii]->Draw();
      }
    }
    
    TCanvas* cdut[nDUT];
    for (int ii=0; ii<nDUT; ii++) {
      if (DUTK[ii] == "X") {
	cdut[ii] = new TCanvas(Form("DUTX_%d", ii), Form("DUTX_%d", ii));
	mctDUT[ii]->Draw();
	fitDUT[ii]->Draw("same");
      }
      if (DUTK[ii] == "U") {
	cdut[ii] = new TCanvas(Form("DUTU_%d", ii), Form("DUTU_%d", ii));
	mctDUT[ii]->Draw();
	fitDUT[ii]->Draw("same");
      }
    }
      
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

  double resoXmeas[nDUT];
  
  for (int ii=0; ii<nDUT; ii++) {
    if (!kNoDraw) {
      TCanvas* cdsx = new TCanvas(Form("DUTSmearingX_%d", ii), Form("DUTSmearingX_%d", ii));
      fitDUTSmeX[ii]->Draw();
    }
    smegaus->SetParameter(0, fitDUTSmeX[ii]->GetMaximum());
    smegaus->SetParameter(1, fitDUTSmeX[ii]->GetMean());
    smegaus->SetParLimits(2, 0.0, 10.0*fitDUTSmeX[ii]->GetRMS());
    smegaus->SetParameter(2, fitDUTSmeX[ii]->GetRMS());
    if (kNoDraw) {
      fitDUTSmeX[ii]->Fit(smegaus, "NQ", "L");
    }
    else {
      fitDUTSmeX[ii]->Fit(smegaus, "Q", "L");
    }
    
    if (!kNoDraw) {
      TCanvas* cdrx = new TCanvas(Form("DUTResidualX_%d", ii), Form("DUTResidualX_%d", ii));
      fitDUTResX[ii]->GetYaxis()->SetRangeUser(0, 1.1*fitDUTSmeX[ii]->GetMaximum());
      fitDUTResX[ii]->Draw();
    }
    resgaus->SetParameter(0, fitDUTResX[ii]->GetMaximum());
    resgaus->SetParameter(1, fitDUTResX[ii]->GetMean());
    resgaus->SetParLimits(2, 0.0, 10.0*fitDUTResX[ii]->GetRMS());
    resgaus->SetParameter(2, fitDUTResX[ii]->GetRMS());
    if (kNoDraw) {
      fitDUTResX[ii]->Fit(resgaus, "NQ", "L");
    }
    else {
      fitDUTResX[ii]->Fit(resgaus, "Q", "L");
      fitDUTSmeX[ii]->Draw("same");
    }
    
    //--------------
    
    double sigmaXmeas = resgaus->GetParameter(2);
    double smeaXmeas = smegaus->GetParameter(2);
    
    resoXmeas[ii] = TMath::Sqrt(TMath::Power(sigmaXmeas, 2.0) - TMath::Power(smeaXmeas, 2.0));
    if (!kNoDraw) {
      printf("DUT X %d) sigma measured: %f, smearing measured: %f\n", ii, sigmaXmeas, smeaXmeas);
      printf("DUT X %d) resolution: %f\n", ii, resoXmeas[ii]);
    }
  }

  //--------------

  if (kNoDraw){
    if (mctX) delete mctX;
    if (mctY) delete mctY;
    if (fitX) delete fitX;
    if (fitY) delete fitY;
    
    if (smegaus) delete smegaus;
    if (resgaus) delete resgaus;
    
    if (mctBeamPipeX) delete mctBeamPipeX;
    if (mctBeamPipeY) delete mctBeamPipeY;
    if (mctThetaX) delete mctThetaX;
    if (mctThetaY) delete mctThetaY;
    
    for (int ii=0; ii<nUpTel; ii++) {
      if (mctUpTelX[ii]) delete mctUpTelX[ii];
    }

    for (int ii=0; ii<nDUT; ii++) {
      if (mctDUT[ii]) delete mctDUT[ii];
    }
    
    for (int ii=0; ii<nDowTel; ii++) {
      if (mctDowTelX[ii]) delete mctDowTelX[ii];
    }
    
    if (fitBeamPipeX) delete fitBeamPipeX;
    if (fitBeamPipeY) delete fitBeamPipeY;
    if (fitThetaX) delete fitThetaX;
    if (fitThetaY) delete fitThetaY;

    for (int ii=0; ii<nDUT; ii++) {
      if (fitDUT[ii]) delete fitDUT[ii];

      if (fitDUTSmeX[ii]) delete fitDUTSmeX[ii];
      if (fitDUTResX[ii]) delete fitDUTResX[ii];
    }
  }

  //--------------
  
  return resoXmeas[component];
}
