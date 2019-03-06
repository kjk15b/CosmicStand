#include <TPad.h>
#include <math.h>
#include <iostream>
using namespace std;

// I deleted a fair amount of lines of code. It did not seem important and seemed rather specific to what this code may have been previously used for.

  // cout << "Starting "  << endl;
  
 #define V_o 192174652.6 // assumed speed of light through a scintillator with index of refraction of ~1.56
 #define h_o 0.3175 // initial pre-defined plane separation distance (m)
 #define delta 0.1 // 10 cm from the top plane PMT

int evtdisplay_4ch(const char *fname = "azm_70cm.root", const bool interactive=true)
{
  gStyle->SetOptStat("nemruo");

  TString name;
  TString title;
  TFile *tfile = new TFile(fname,"READ");
  
  TFile *outfile = new TFile("outfile.root","RECREATE", "histograms from DRS4 waveforms" ,9);
  outfile->Cd("/");

  Int_t evt;
  const int NCH = 4;
  const int NSAMPLES = 1024;
  Float_t t[NCH][NSAMPLES];  // time
  Float_t ch[NCH][NSAMPLES]; // voltage
  TTree *tree = (TTree*)tfile->Get("t");
  tree->SetBranchAddress("evt",&evt);
  for (int ich=0; ich<NCH; ich++) {
     TString bName="t"; bName+=ich; // t0, t1, ...
     tree->SetBranchAddress(bName,&t[ich]);
     bName="ch"; bName+=ich;  // ch0, ch1, ...
     tree->SetBranchAddress(bName,&ch[ich]);
  }

  TTree *dst = new TTree("tdst", "");

  TCanvas *c_display = new TCanvas("drs4","drs4 display",1000,700);
  c_display->Divide(1,NCH);

  TGraph *gpulse[NCH];
  for (int ich=0; ich<NCH; ich++)
  {
    gpulse[ich] = new TGraph();
    name = "gch"; name += ich;
    gpulse[ich]->SetName(name);
    gpulse[ich]->SetMarkerStyle(20);
    gpulse[ich]->SetMarkerSize(0.4);
    gpulse[ich]->SetMarkerColor(ich+1);
    gpulse[ich]->SetLineColor(ich+1);
  }

  c_display->cd(1);
  
  
  
  // Comparing timing differences of detectors with peaks, using detector 0 (finger scintillator) as reference
  TH1F* hPeakTime0minus1 = new TH1F("hPeakTime0minus1","Detector 0 time minus detector 1 time; Time (ps); Counts; ",200,-10,90);
  TH1F* hPeakTime0minus2 = new TH1F("hPeakTime0minus2","Detector 0 time minus detector 2 time; Time (ps); Counts; ",200,-10,90);
  TH1F* hPeakTime0minus3 = new TH1F("hPeakTime0minus3","Detector 0 time minus detector 3 time; Time(ps); Counts;",200,-10,90);
  
  TH1F* hPeakTime2minus3 = new TH1F("hPeakTime2minus3","Detector 2 time minus detector 3 time; Time(ps); Counts;", 200,-10,90);
  
 
 // Comparing the timing differences of the finger scintillator relative to the rest of the scintillators with CFD
  TH1F* hCFDTime0minus1 = new TH1F("hCFDTime0minus1","Detector 0 time minus detector 1 time; Time (ps); Counts; ",200,-10,90);
  TH1F* hCFDTime0minus2 = new TH1F("hCFDTime0minus2","Detector 0 time minus detector 2 time; Time (ps); Counts; ",200,-10,90);
  TH1F* hCFDTime0minus3 = new TH1F("hCFDTime0minus3","Detector 0 time minus detector 3 time; Time(ps); Counts;",200,-10,90);

  // CFD timing difference for either end of scintillating plane
TH1F* hCFDTime2minus1 = new TH1F("hCFDTime2minus3","Detector 2 time minus detector 3 time; Time(ps);Counts",200,-10,90);

// Histogram of Azimuthal Gamma angle
TH1F* hGamma = new TH1F("hGamma","Azimuthal Cosmic Ray Distribution;Angle(Degrees);Counts",200,-10,90);

// Histogram of bottom plane distribution
TH1F* hPosition = new TH1F("hPosition","Scintillator Plane Distribution; Position(m); Counts", 200, -10,90)'


  TH1F* hSaturation[NCH];
  TH1F* hPeak[NCH];
  TH1F* hPeakTime[NCH];
  TH1F* hCFDTime[NCH];
  TH1F* hgamma;
  TH1F* hposition;
  for (int ich=0; ich<NCH; ich++) {
    name="hSaturationCh"; name+=ich;
    hSaturation[ich] = new TH1F(name,name,NSAMPLES,0,NSAMPLES);
    name="hPeakCh"; name+=ich;
    hPeak[ich] = new TH1F(name,name,100,-0.5,0.5);
    name="hPeakTimeCh"; name+=ich;
    hPeakTime[ich] = new TH1F(name,name,200,0,200);
    name="hCFDTimeCh"; name+=ich;
    hCFDTime[ich] = new TH1F(name,name,200,0,200);
  }
  
  hgamma = new TH1F("hgamma", "hgamma", 200,0,200); // Plots for new histograms for positional tracking and angular distributions
  hposition = new TH1F("hposition", "hposition",200,0,200);
  

  int nskipped = 0; // number of skipped, uninteresting events

  Double_t satThresh=-0.490;
  Double_t trigThresh=-0.016; // negative pulses (> this or >= this?)

  Long64_t nentries = tree->GetEntries();
  cout << "Found " << nentries << " events" << endl;
  for (int ievt=0; ievt<nentries; ievt++)
  {
    tree->GetEntry(ievt);

    // Read in data
    int trig = 0;
    int ctrig[NCH];
    Double_t cPeak[NCH];
    Double_t cPeakTime[NCH];
    Double_t cMax[NCH];
    Double_t cMaxTime[NCH];
    Double_t cfdTime[NCH];
    int cSatSamples[NCH];
    for (int ich=0; ich<NCH; ich++)
    {
      int nSatSamples=0;
      cSatSamples[ich]=0;
      ctrig[ich]=0;
      cPeak[ich]=0.0;
      cPeakTime[ich]=0.0;
      cfdTime[ich]=0.0;
      for (int isamp=0; isamp<NSAMPLES; isamp++)
      {
        gpulse[ich]->SetPoint(isamp,t[ich][isamp],ch[ich][isamp]);

        if (ch[ich][isamp] < cPeak[ich]) {
          cPeak[ich]=ch[ich][isamp];
          cPeakTime[ich]=t[ich][isamp];
        }
        if (ch[ich][isamp] > cMax[ich]) {
          cMax[ich]=ch[ich][isamp];
          cMaxTime[ich]=t[ich][isamp];
        }


        if (ch[ich][isamp] <= satThresh) {
          nSatSamples++;
          cSatSamples[ich]++;
        }
        if (ch[ich][isamp] <= trigThresh) {
          trig=1;
          ctrig[ich]=1;
        }
      }
      
      // After peak finding, perform CFD search
      for (int isamp=0; isamp<NSAMPLES-1; isamp++)
      {
        if (isamp >= NSAMPLES) continue;
        const double CFDfraction=0.25;
        if (ch[ich][isamp] < cPeak[ich]*CFDfraction && ch[ich][isamp+1] >= cPeak[ich]*CFDfraction) {
          cfdTime[ich]=t[ich][isamp];
        }
      }
      if (nSatSamples > 0) hSaturation[ich]->Fill(nSatSamples);
      hPeak[ich]->Fill(cPeak[ich]);
      hPeakTime[ich]->Fill(cPeakTime[ich]);
      hCFDTime[ich]->Fill(cfdTime[ich]);
      printf("peak[%d] = %g\n", ich, cPeak[ich]);
    }

	
	// For hPeakTime0minus1
	// Fill is done only when not saturated
    if (cSatSamples[0]==0 && cSatSamples[1]==0) { 
      hPeakTime0minus1->Fill(cPeakTime[0]-cPeakTime[1]);
      hCFDTime0minus1->Fill(cfdTime[0]-cfdTime[1]);
    }
	
	// For hPeakTime0minus2
	if (cSatSamples[0]==0 && cSatSamples[2] == 0) {
		hPeakTime0minus2->Fill(cPeak[0]-cPeak[2];
		hCFDTime0minus2->Fill(cfdTime[0]-cfdTime[2]);
	}
	
	// For hPeakTime0minus3
    if (cSatSamples[0]==0 && cSatSamples[3]==0) { 
      hPeakTime0minus3->Fill(cPeakTime[0]-cPeakTime[3]);
      hCFDTime0minus3->Fill(cfdTime[0]-cfdTime[3]);
    }
	
	// For hPeakTime2minus3
	if(cSatSamples[2]==0 && cSatSamples[3]==0)
	{
		hPeakTime2minus3->Fill(cPeakTime[2]-cPeakTime[3]);
		hCFDTime2minus3->Fill(cfdTime[2]-cfdTime[3]);
	}
	
// Insert data analysis for gamma distribution and positional resolution
	if(cSatSamples[3]==0 && cSatSamples[2]==0)
	{
		Double_t  xf = V_o * (cfdTime[2] - cfdTime[3]);
		hPosition->Fill(xf);
		xf = h_o / (xf - delta);
		Double_t gamma = atan(xf);
		gamma = (gamma * 180) / 3.1456;
		hGamma->Fill(gamma);
	}

    if (trig)
    {
      for (int ich=0; ich<NCH; ich++)
      {
        c_display->cd(ich+1);
        if (interactive) gpulse[ich]->Draw("alp");
      }
      if ( gPad!= 0 )
      {
        gPad->Modified();
        gPad->Update();
      }

      cout << "Skipped " << nskipped << " events" << endl;
      nskipped = 0;

      string junk="";
      cout << evt << " ? ";
      if (interactive) cin >> junk;
      if ( junk[0] == 'q' ) break;
      else if ( junk[0] == 'w' )
      {
        c_display->SaveAs(".png");
      }
    }
    else
    {
      nskipped++;
    }

  }

  outfile->Cd("/");
  hPeakTime0minus1->Write();
  hPeakTime0minus2->Write();
  hPeakTime0minus3->Write();
  hPeakTime2minus3->Write();
  hCFDTime0minus1->Write();
  hCFDTime0minus2->Write();
  hCFDTime0minus3->Write();
  hCFDTime2minus3->Write();
  hGamma->Write();
  hPosition->Write();

  
  for (int ich=0; ich<NCH; ich++) {
    hSaturation[ich]->Write();
    hPeak[ich]->Write();
    hPeakTime[ich]->Write();
    hCFDTime[ich]->Write();
  }

  outfile->Close();
  return 1;
}

