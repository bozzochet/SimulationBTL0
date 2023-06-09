#include "TCanvas.h"
#include "TH1F.h"
#include "TF1.h"
#include "TRandom.h"
#include "TString.h"
#include "TGraphErrors.h"

double resoX = 0.01; //10 um
double resoY = 0.01; //10 um
double resoU = 0.01; //10 um
double rotU = 45;//degrees

struct STResult {
  TString Type;
  double ExpResoX;//mm
  double MeasResoX;//mm
  double MeasResoErrX;//mm
};

void Single();
void Full();
std::vector<STResult> SingleTrial(int events=1000000, bool kNoDraw=false);

//----

/* rotation and anti-rotation:
   u = x*cosα - y*sinα;
   v = x*sinα + y*cosα;
   
   x = u*cosα + v*sinα;
   y = -u*sinα + v*cosα;
*/

/* mixed:
   x = f(u,y) = u*cosα + (x*sinα+y*cosα)*sinα;
   x = u*cosα + x*sin^2α + y*cosα*sinα;
   x (1 - sin^2α) = u*cosα + y*cosα*sinα;
   x = (u*cosα + y*cosα*sinα)/(1 - sin^2α);

   ...

   v = f(u, y) = (u*cosα + v*sinα)*sinα + y*cosα;
   v = u*cosα*sinα + v*sin^2α + y*cosα;
   v (1 - sin^2α) = u*cosα*sinα + y*cosα;
   v = (u*cosα*sinα + y*cosα)/(1 - sin^2α);
*/

/*
  x = f(u,y) --> dx = dx/du * du + dx/dy * dY;
  --> dx = du*cosα/(1 - sin^2α) + dy*sinα*cosα/(1 - sin^2α);
  = du/cosα + dy*tanα
*/

double Rotation(double x, double y, double angle, int component){

  if (component==0) { //u
    return x*TMath::Cos(TMath::DegToRad()*angle) - y*TMath::Sin(TMath::DegToRad()*angle);
  }
  else if (component==1) {//v
    return x*TMath::Sin(TMath::DegToRad()*angle) + y*TMath::Cos(TMath::DegToRad()*angle);
  }
  else {
    printf("Rotation) Wrong component (%d)...\n", component);
  }
  
  return -999.9;  
}

double MixedRotation(double a, double b, double angle, int component){

  if (component==0) { //x from u and y
    double u = a;
    double y = b;
    return (u*TMath::Cos(TMath::DegToRad()*angle) + y*TMath::Cos(TMath::DegToRad()*angle)*TMath::Sin(TMath::DegToRad()*angle))/(1.0 - TMath::Power(TMath::Sin(TMath::DegToRad()*angle), 2.0));
  }
  else if (component==3) {//v from u and y
    double u = a;
    double y = b;
    return (u*TMath::Cos(TMath::DegToRad()*angle)*TMath::Sin(TMath::DegToRad()*angle) + y*TMath::Cos(TMath::DegToRad()*angle))/(1.0 - TMath::Power(TMath::Sin(TMath::DegToRad()*angle), 2.0));
  }
  else {
    printf("MixedRotation) Wrong or not yet implemented component (%d)...\n", component);
  }
  
  return -999.9;  
}

void Single(){

  SingleTrial();
  
  return;
}

void Full(){
  //  int trials=100;
  int trials=1000;
  //  int events=10000;//10k
  int events=1000000;//1M

  std::vector<STResult> _results = SingleTrial(events, true);//first time to retrive number of "types"

  std::vector<TH1*> resoXmeas(_results.size());
  for (int rr=0; rr<((int)(_results.size())); rr++) {
    resoXmeas[rr] = new TH1F(Form("resoXmeas_%d", rr), Form("reasoXmeas_%d (%s)", rr, _results[rr].Type.Data()), 1000, -0.1, 0.1);
    resoXmeas[rr]->GetXaxis()->SetTitle("(#sigma_x measured - #sigma_x true) / #sigma_x true");
    resoXmeas[rr]->GetYaxis()->SetTitle("Trials");
  }
  
  for (int tt=0; tt<trials; tt++) {
    for (int rr=0; rr<((int)(_results.size())); rr++) {
      double exp = _results[rr].ExpResoX;
      double meas = _results[rr].MeasResoX;
      resoXmeas[rr]->Fill((meas-exp)/exp);
      printf("%d, %d) (%f -%f)/%f = %f\n", tt, rr, meas, exp, exp, (meas-exp)/exp);
    }
    if (tt!=(trials-1)) {//last time not needed...
      _results = SingleTrial(events, true);
    }
  }

  for (int rr=0; rr<((int)(_results.size())); rr++) {
    TCanvas* c = new TCanvas(Form("RelativeErrorOnResolution_%d", rr), Form("RelativeErrorOnResolution_%d (%s)", rr, _results[rr].Type.Data()));
    resoXmeas[rr]->Draw();
  }

  return;
}

std::vector<STResult> SingleTrial(int events, bool kNoDraw) {

  // do not use the downstream telescope
  bool kNoDowTel =false;

  // half distances between Up - DUT - Down
  bool kShorten = false;

  // double distances between Up - DUT - Down
  bool kExtend = false;
  
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
  
  if (kShorten) {
    DUTZ[0] = 1250.0;
    DUTZ[1] = 1250.1;
  }
  else if (kExtend) {
    DUTZ[0] = 2000.0;
    DUTZ[1] = 2000.1;
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
  else if (kExtend) {
    double DowTelZLonger[nDowTel] = {3000.0, 3000.1, 3060.0, 3060.1};
    for (int ii=0; ii<nDowTel; ii++) {
      DowTelZ[ii] = DowTelZLonger[ii];
    }
  }
  TString DowTelK[nDowTel] = {"X", "Y", "X", "Y"};

  //------------
  
  TF1* mctX = new TF1("mctX", "[0] + [1]*x", 0.0, 2.0);
  TF1* mctY = new TF1("mctY", "[0] + [1]*x", 0.0, 2.0);
  TF1* fitX = new TF1("fitX", "[0] + [1]*x", 0.0, 2.0);
  TF1* fitY = new TF1("fitY", "[0] + [1]*x", 0.0, 2.0);

  TF1* smegausX = new TF1("smegausX", "[0]*exp(-0.5*((x-[1])/[2])*((x-[1])/[2]))", -50, 50);
  smegausX->SetLineColor(kBlack);
  smegausX->SetLineWidth(3);
  TF1* smegausY = new TF1("smegausY", "[0]*exp(-0.5*((x-[1])/[2])*((x-[1])/[2]))", -50, 50);
  smegausY->SetLineColor(kBlue+1);
  smegausY->SetLineWidth(3);
  TF1* resgausX = new TF1("resgausX", "[0]*exp(-0.5*((x-[1])/[2])*((x-[1])/[2]))", -50, 50);
  resgausX->SetLineColor(kRed+2);
  resgausX->SetLineWidth(3);  
  TF1* resgausY = new TF1("resgausY", "[0]*exp(-0.5*((x-[1])/[2])*((x-[1])/[2]))", -50, 50);
  resgausY->SetLineColor(kRed+2);
  resgausY->SetLineWidth(3);
  
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
  TH1* fitDUTSmeY[nDUT];
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
    fitDUTSmeY[ii] = new TH1F(Form("fitDUTSmeY_%d", ii), Form("fitDUTSmeY_%d", ii), 1000, -reslim, reslim);
    fitDUTSmeY[ii]->GetXaxis()->SetTitle("Y-coo residual at DUT (#mum)");
    fitDUTSmeY[ii]->GetYaxis()->SetTitle("Entries");
    fitDUTSmeY[ii]->SetLineColor(kBlue+1);
    fitDUTSmeY[ii]->SetMarkerColor(kBlue+1);
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
	double u = Rotation(x, y, rotU, 0);
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

    /*
      just for double check:
      ****************************

      if y = q + mx:
      
      m = (mean(x*y) - mean(x)*mean(y)) / (mean(x^2) - mean(x)*mean(x))
      q = (mean(x^2)*mean(y) - mean(x)*mean(x*y)) / (mean(x^2) - mean(x)*mean(x)) = mean(y) - m*mean(x)

      if errors on y not all equal:
      
      mean(x) = sum(x_i/sigma_i^2)/sum(1/sigma_i^2)
      mean(y) = sum(y_i/sigma_i^2)/sum(1/sigma_i^2)
      sigma^2 = N/sum(1/sigma_i^2)
      
      S1 = sum(1/sigma_i^2)
      Sz = sum(x_i/sigma_i^2)
      Szz = sum(x_i*x_i/sigma_i^2)
      Sx = sum(y_i/sigma_i^2)
      Szx = sum(x_y*y_i/sigma_i^2)
      Dz = mean(x^2) - mean(x)*mean(x)

      U = k*Vmatrix
      Vmatrix = mean(x^2)    -mean(x)
                 -mean(x)       1


      where
      k = (sigma^2/(N*Dz))
        = (N/sum(1/sigma_i^2))/(N*Dz))
        = 1/(sum(1/sigma_i^2)*Dz))
	= 1/(S1*Dz)

      and so
      mean(x) = sum(x_i/sigma_i^2)/sum(1/sigma_i^2) = Sz/S1
      mean(x^2) = sum(x_i*x_i/sigma_i^2)/sum(1/sigma_i^2) = Szz/S1 
      
      -->
      sigma^2(q): U00 = k*mean(x^2)
      sigma^2(m): U11 = k
      mixed_term: U01 = U10 = -k*mean(x)

      -->
      sigma^2(y) = sigma^2(q) + x^2*sigma^2(m) + 2x*mixed_term
      sigma^2(y) = U00 + U11*x^2 + 2x*U01
     */
    double S1 = 0, Sz = 0, Szz = 0, Sx = 0, Szx = 0;
    
    for (int ii=0; ii<nUpTel; ii++) {
      if (UpTelK[ii] == "X") {
	double z = UpTelZ[ii];
	double x = gRandom->Gaus(mctX->Eval(z), resoX);
	grX.SetPoint(grX.GetN(), z, x);
	grX.SetPointError(grX.GetN()-1, 0.0, resoX);
	S1 += 1.0/TMath::Power(resoX, 2.0);
	Sz += z/TMath::Power(resoX, 2.0);
	Szz += z*z/TMath::Power(resoX, 2.0);
	Sx += x/TMath::Power(resoX, 2.0);
	Szx += z*x/TMath::Power(resoX, 2.0);
      }
      if (UpTelK[ii] == "Y") {
	double z = UpTelZ[ii];
	double y = gRandom->Gaus(mctY->Eval(z), resoY);
	grY.SetPoint(grY.GetN(), z, y);
	grY.SetPointError(grY.GetN()-1, 0.0, resoY);
      }
    }

    if (!kNoDowTel) {
      for (int ii=0; ii<nDowTel; ii++) {
	if (DowTelK[ii] == "X") {
	  double z = DowTelZ[ii];
	  double x = gRandom->Gaus(mctX->Eval(z), resoX);
	  grX.SetPoint(grX.GetN(), z, x);
	  grX.SetPointError(grX.GetN()-1, 0.0, resoX);
	  S1 += 1.0/TMath::Power(resoX, 2.0);
	  Sz += z/TMath::Power(resoX, 2.0);
	  Szz += z*z/TMath::Power(resoX, 2.0);
	  Sx += x/TMath::Power(resoX, 2.0);
	  Szx += z*x/TMath::Power(resoX, 2.0);
	}
	else if (DowTelK[ii] == "Y") {
	  double z = DowTelZ[ii];
	  double y = gRandom->Gaus(mctY->Eval(z), resoY);
	  grY.SetPoint(grY.GetN(), z, y);
	  grY.SetPointError(grY.GetN()-1, 0.0, resoY);
	}
      }
    }

    double Dz = S1 * Szz - Sz * Sz;
    double X0 = (Sx * Szz - Sz * Szx) / Dz;
    double mX = (S1 * Szx - Sz * Sx) / Dz;
    double X0err = TMath::Sqrt(Szz / Dz);
    double mXerr = TMath::Sqrt(S1 / Dz);

    double meanz = Sz/S1;
    double meanzsq = Szz/S1;
    double Dz_plain = meanzsq - meanz*meanz; 
    double k = 1.0/(S1*Dz_plain);
    //    printf("%f %f %f --> %E\n", meanz, meanzsq, Dz_plain, k);

    double U00 = k*meanzsq;
    double U11 = k;
    double U01 = -k*meanz;
    
    grX.Fit(fitX, "NQ");
    grY.Fit(fitY, "NQ");

    /*
    printf("%f %f\n", fitX->GetParameter(0), X0);
    printf("%f %f\n", fitX->GetParameter(1), mX);
    printf("%f %f\n", fitX->GetParError(0), X0err);
    printf("%f %f\n", fitX->GetParError(1), mXerr);
    */

    fitBeamPipeX->Fill(fitX->Eval(BeamPipeZ));
    fitBeamPipeY->Fill(fitY->Eval(BeamPipeZ));
    fitThetaX->Fill(TMath::ATan(1000.0*fitX->GetParameter(1)));
    fitThetaY->Fill(TMath::ATan(1000.0*fitY->GetParameter(1)));

    for (int ii=0; ii<nDUT; ii++) {
      if (DUTK[ii] == "X") {
	double z = DUTZ[ii];
	double xtruth = mctX->Eval(z);
	double xtrack = fitX->Eval(z);
	double xmeas = gRandom->Gaus(xtruth, resoX);
	double ytruth = mctY->Eval(z);
	double ytrack = fitY->Eval(z);
	fitDUT[ii]->Fill(xtrack);
	fitDUTSmeX[ii]->Fill(1000.0*(xtrack - xtruth));
	fitDUTResX[ii]->Fill(1000.0*(xmeas - xtrack));
	fitDUTSmeY[ii]->Fill(1000.0*(ytrack - ytruth));
	double sigmaxsq = U00 + U11*z*z + 2*z*U01;
	//	printf("%f (μm)\n", 1000.0*TMath::Sqrt(sigmaxsq));
      }
      else if (DUTK[ii] == "U") {
	double xtruth = mctX->Eval(DUTZ[ii]);
	double xtrack = fitX->Eval(DUTZ[ii]);
	double ytruth = mctY->Eval(DUTZ[ii]);
	double ytrack = fitY->Eval(DUTZ[ii]);
	double utruth = Rotation(xtruth, ytruth, rotU, 0);
	double utrack = Rotation(xtrack, ytrack, rotU, 0);
	double vtruth = Rotation(xtruth, ytruth, rotU, 1);
	double vtrack = Rotation(xtrack, ytrack, rotU, 1);
	double umeas = gRandom->Gaus(utruth, resoU);
	// we need a "v":
	//	double vmeas = vtruth; //v from truth
	//	double vmeas = vtrack; //v from track
	//	double vmeas = gRandom->Gaus(vtruth, resoU); //v from another layer
	double vmeas = MixedRotation(umeas, ytrack, rotU, 3);//from u and y
	double xmeas = Rotation(umeas, vmeas, -rotU, 0);
	double ymeas = Rotation(umeas, vmeas, -rotU, 1);
	fitDUT[ii]->Fill(utrack);
	fitDUTSmeX[ii]->Fill(1000.0*(xtrack - xtruth));
	fitDUTResX[ii]->Fill(1000.0*(xmeas - xtrack));
	fitDUTSmeY[ii]->Fill(1000.0*(ytrack - ytruth));
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

  std::vector<STResult> resoXmeas;
  
  for (int ii=0; ii<nDUT; ii++) {
    
    if (!kNoDraw) {
      TCanvas* cdsx = new TCanvas(Form("DUTSmearingX_%d", ii), Form("DUTSmearingX_%d", ii));
      fitDUTSmeX[ii]->Draw();
    }
    smegausX->SetParameter(0, fitDUTSmeX[ii]->GetMaximum());
    smegausX->SetParameter(1, fitDUTSmeX[ii]->GetMean());
    smegausX->SetParLimits(2, 0.0, 10.0*fitDUTSmeX[ii]->GetRMS());
    smegausX->SetParameter(2, fitDUTSmeX[ii]->GetRMS());
    if (kNoDraw) {
      fitDUTSmeX[ii]->Fit(smegausX, "NQ", "L");
    }
    else {
      fitDUTSmeX[ii]->Fit(smegausX, "Q", "L");
    }

    if (!kNoDraw) {
      TCanvas* cdsy = new TCanvas(Form("DUTSmearingY_%d", ii), Form("DUTSmearingY_%d", ii));
      fitDUTSmeY[ii]->Draw();
    }
    smegausY->SetParameter(0, fitDUTSmeY[ii]->GetMaximum());
    smegausY->SetParameter(1, fitDUTSmeY[ii]->GetMean());
    smegausY->SetParLimits(2, 0.0, 10.0*fitDUTSmeY[ii]->GetRMS());
    smegausY->SetParameter(2, fitDUTSmeY[ii]->GetRMS());
    if (kNoDraw) {
      fitDUTSmeY[ii]->Fit(smegausY, "NQ", "L");
    }
    else {
      fitDUTSmeY[ii]->Fit(smegausY, "Q", "L");
    }
    
    if (!kNoDraw) {
      TCanvas* cdrx = new TCanvas(Form("DUTResidualX_%d", ii), Form("DUTResidualX_%d", ii));
      fitDUTResX[ii]->GetYaxis()->SetRangeUser(0, 1.1*fitDUTSmeX[ii]->GetMaximum());
      fitDUTResX[ii]->Draw();
    }
    resgausX->SetParameter(0, fitDUTResX[ii]->GetMaximum());
    resgausX->SetParameter(1, fitDUTResX[ii]->GetMean());
    resgausX->SetParLimits(2, 0.0, 10.0*fitDUTResX[ii]->GetRMS());
    resgausX->SetParameter(2, fitDUTResX[ii]->GetRMS());
    if (kNoDraw) {
      fitDUTResX[ii]->Fit(resgausX, "NQ", "L");
    }
    else {
      fitDUTResX[ii]->Fit(resgausX, "Q", "L");
      fitDUTSmeX[ii]->Draw("same");
      if (DUTK[ii] == "U") {
	fitDUTSmeY[ii]->Draw("same");
      }
    }
    
    //--------------
    
    double sigmaXmeas = resgausX->GetParameter(2);
    double sigmaYmeas = resgausY->GetParameter(2);
    double smeaYmeas = smegausY->GetParameter(2);
    double smeaXmeas = smegausX->GetParameter(2);

    double sigmaXmeas_err = resgausX->GetParError(2);
    double sigmaYmeas_err = resgausY->GetParError(2);
    double smeaYmeas_err = smegausY->GetParError(2);
    double smeaXmeas_err = smegausX->GetParError(2);
    
    STResult _stresult;

    _stresult.Type = DUTK[ii];
    if (DUTK[ii] == "X") {
      _stresult.ExpResoX = resoX;
      double measresoX = TMath::Sqrt(TMath::Power(sigmaXmeas, 2.0)
				     - TMath::Power(smeaXmeas, 2.0)
				     );
      double measresoerrX = TMath::Sqrt(TMath::Power(sigmaXmeas*sigmaXmeas_err/measresoX, 2.0)
						 + TMath::Power(smeaXmeas*smeaXmeas_err/measresoX, 2.0)
						 );
      if (!kNoDraw) {
	printf("DUT X %d (%s)) ", ii, DUTK[ii].Data());
	printf("sigma measured: %f +- %f (μm), ", sigmaXmeas, sigmaXmeas_err);
	printf("smearing measured: %f +- %f (μm)\n", smeaXmeas, smeaXmeas_err);
	printf("--> resolution: %f +- %f (μm) - %f %%\n", measresoX, measresoerrX, 100.0*measresoerrX/measresoX);
      }
      _stresult.MeasResoX = 0.001*measresoX;
      _stresult.MeasResoErrX = 0.001*measresoerrX;
    }
    else if (DUTK[ii] == "U") {
      /*
	x = f(u,y) --> dx = dx/du * du + dx/dy * dY;
	--> dx = du*cosα/(1 - sin^2α) + dy*sinα*cosα/(1 - sin^2α);
	= du/cosα + dy*tanα
      */
      _stresult.ExpResoX = resoU/TMath::Cos(TMath::DegToRad()*rotU);
      double measresoX = TMath::Sqrt(TMath::Power(sigmaXmeas, 2.0)
				     - TMath::Power(smeaXmeas, 2.0)
				     - TMath::Power(TMath::Tan(TMath::DegToRad()*rotU)*smeaYmeas, 2.0)
				     );
      double measresoerrX = TMath::Sqrt(TMath::Power(sigmaXmeas*sigmaXmeas_err/measresoX, 2.0)
					+ TMath::Power(smeaXmeas*smeaXmeas_err/measresoX, 2.0)
					+ TMath::Power(TMath::Power(TMath::Tan(TMath::DegToRad()*rotU), 2.0)*smeaYmeas*smeaYmeas_err/measresoX, 2.0)
					);

      if (!kNoDraw) {
	printf("DUT X %d (%s)) ", ii, DUTK[ii].Data());
	printf("sigma measured: %f +- %f (μm), ", sigmaXmeas, sigmaXmeas_err);
	printf("smearing X measured: %f +- %f (μm), ", smeaXmeas, smeaXmeas_err);
	printf("smearing Y measured: %f +- %f (μm)\n", smeaYmeas, smeaYmeas_err);
	printf("--> resolution: %f +- %f (μm) - %f %%\n", measresoX, measresoerrX, 100.0*measresoerrX/measresoX);
      }
      _stresult.MeasResoX = 0.001*measresoX;
      _stresult.MeasResoErrX = 0.001*measresoerrX;
    }

    resoXmeas.push_back(_stresult);
  }

  //--------------

  if (kNoDraw){
    if (mctX) delete mctX;
    if (mctY) delete mctY;
    if (fitX) delete fitX;
    if (fitY) delete fitY;
    
    if (smegausX) delete smegausX;
    if (smegausY) delete smegausY;
    if (resgausX) delete resgausX;
    if (resgausY) delete resgausY;
    
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
      if (fitDUTSmeY[ii]) delete fitDUTSmeY[ii];
      if (fitDUTResX[ii]) delete fitDUTResX[ii];
    }
  }

  //--------------

  return resoXmeas;
}
