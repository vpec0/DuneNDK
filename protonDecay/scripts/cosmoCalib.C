#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TEventList.h"
#include "TH1D.h"

void cosmoCalib() {

  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetOptStat(0);
  gStyle->SetStatX(0.9);
  gStyle->SetStatY(0.937);
  gStyle->SetStatW(0.3);
  gStyle->SetStatH(0.15);
  gStyle->SetHistLineColor(1);
  TGaxis::SetMaxDigits(3);

  TCanvas* canv = new TCanvas("canv","",800,600);
  TLegend* leg = new TLegend(0.6,0.7,0.85,0.85);

  TFile* inFile = new TFile("/pnfs/dune/scratch/users/wallbank/v06_45_01/ana/comsic_calib/anahist.root","r");
  //TFile* inFile = new TFile("cosmo_calib.root","r");
  TTree* tree = (TTree*)inFile->Get("cosmocalib/CosmoCalib");
  if (!tree) {
    std::cout << "Input tree not found." << std::endl;
    return;
  }

  TEventList tpcList("tpcList"), stoppingList("stoppingList"), crossList("crossList"), innerCrossList("innerCrossList"), gapCrossList("gapCrossList");
  int NMuons = tree->GetEntriesFast();
  tree->Draw(">>tpcList","MuonTPC");
  int NMuonsTPC = tpcList.GetN();
  tree->Draw(">>stoppingList","MuonStopping");
  int NMuonsStopping = stoppingList.GetN();
  tree->Draw(">>crossList","MuonCross");
  int NMuonsCross = crossList.GetN();
  tree->Draw(">>innerCrossList","MuonInnerCross");
  int NMuonsInnerCross = innerCrossList.GetN();
  tree->Draw(">>gapCrossList","MuonGapCross");
  int NMuonsGapCross = gapCrossList.GetN();

  std::cout << "There were " << NMuons << " muons (not just primaries)" << std::endl
  	    << "  -- " << NMuonsTPC << " depositing in the TPC (" << (float)NMuonsTPC*100/(float)NMuons << "%)" << std::endl
  	    << "  -- " << NMuonsStopping << " stopping (" << (float)NMuonsStopping*100/(float)NMuonsTPC << "%)" << std::endl
  	    << "  -- " << NMuonsCross << " crossing APA/CPA frames (" << (float)NMuonsCross*100/(float)NMuonsTPC << "%)" << std::endl
  	    << "  -- " << NMuonsInnerCross << " crossing inner APA/CPA frames (" << (float)NMuonsInnerCross*100/(float)NMuonsTPC << "%)" << std::endl
  	    << "  -- " << NMuonsGapCross << " crossing APA/CPA gaps (" << (float)NMuonsGapCross*100/(float)NMuonsTPC << "%)" << std::endl;

  TH1D* hStoppingEnergy = new TH1D("StoppingEnergy",";True Total Energy (GeV);",100,0,50000);
  TH1D* hNonStoppingEnergy = new TH1D("NonStoppingEnergy",";True Total Energy (GeV);",100,0,50000);
  tree->Draw("MuonEnergy>>StoppingEnergy","MuonStopping");
  tree->Draw("MuonEnergy>>NonStoppingEnergy","!MuonStopping");
  canv->cd();
  canv->Clear();
  canv->SetLogy();
  hStoppingEnergy->SetLineColor(1);
  hStoppingEnergy->SetFillColor(1);
  hStoppingEnergy->SetFillStyle(3003);
  //hStoppingEnergy->Scale(1./hStoppingEnergy->GetEntries());
  hNonStoppingEnergy->SetLineColor(2);
  hNonStoppingEnergy->SetFillColor(2);
  hNonStoppingEnergy->SetFillStyle(3003);
  //hNonStoppingEnergy->Scale(1./hNonStoppingEnergy->GetEntries());
  hNonStoppingEnergy->Draw();
  hStoppingEnergy->Draw("same");
  leg->Clear();
  leg->AddEntry(hStoppingEnergy,    "Stopping",     "f");
  leg->AddEntry(hNonStoppingEnergy, "Non stopping", "f");
  leg->Draw();
  canv->SaveAs("StoppingMuonsEnergy.eps");
  canv->SetLogy(0);

  TH1D* hStoppingAngleTheta = new TH1D("StoppingAngleTheta",";True Angle Theta (deg);",100,0,180);
  TH1D* hNonStoppingAngleTheta = new TH1D("NonStoppingAngleTheta",";True Angle Theta (deg);",100,0,180);
  tree->Draw("MuonAngleTheta>>StoppingAngleTheta","MuonStopping");
  tree->Draw("MuonAngleTheta>>NonStoppingAngleTheta","!MuonStopping");
  canv->cd();
  canv->Clear();
  canv->SetLogy();
  hStoppingAngleTheta->SetLineColor(1);
  hStoppingAngleTheta->SetFillColor(1);
  hStoppingAngleTheta->SetFillStyle(3003);
  //hStoppingAngleTheta->Scale(1./hStoppingAngleTheta->GetEntries());
  hNonStoppingAngleTheta->SetLineColor(2);
  hNonStoppingAngleTheta->SetFillColor(2);
  hNonStoppingAngleTheta->SetFillStyle(3003);
  //hNonStoppingAngleTheta->Scale(1./hNonStoppingAngleTheta->GetEntries());
  hNonStoppingAngleTheta->Draw();
  hStoppingAngleTheta->Draw("same");
  leg->Clear();
  leg->AddEntry(hStoppingAngleTheta,    "Stopping",     "f");
  leg->AddEntry(hNonStoppingAngleTheta, "Non stopping", "f");
  leg->Draw();
  canv->SaveAs("StoppingMuonsAngleTheta.eps");
  canv->SetLogy(0);

  TH1D* hStoppingAnglePhi = new TH1D("StoppingAnglePhi",";True Angle Phi (deg);",100,0,180);
  TH1D* hNonStoppingAnglePhi = new TH1D("NonStoppingAnglePhi",";True Angle Phi (deg);",100,0,180);
  tree->Draw("MuonAnglePhi>>StoppingAnglePhi","MuonStopping");
  tree->Draw("MuonAnglePhi>>NonStoppingAnglePhi","!MuonStopping");
  canv->cd();
  canv->Clear();
  canv->SetLogy();
  hStoppingAnglePhi->SetLineColor(1);
  hStoppingAnglePhi->SetFillColor(1);
  hStoppingAnglePhi->SetFillStyle(3003);
  //hStoppingAnglePhi->Scale(1./hStoppingAnglePhi->GetEntries());
  hNonStoppingAnglePhi->SetLineColor(2);
  hNonStoppingAnglePhi->SetFillColor(2);
  hNonStoppingAnglePhi->SetFillStyle(3003);
  //hNonStoppingAnglePhi->Scale(1./hNonStoppingAnglePhi->GetEntries());
  hNonStoppingAnglePhi->Draw();
  hStoppingAnglePhi->Draw("same");
  leg->Clear();
  leg->AddEntry(hStoppingAnglePhi,    "Stopping",     "f");
  leg->AddEntry(hNonStoppingAnglePhi, "Non stopping", "f");
  leg->Draw();
  canv->SaveAs("StoppingMuonsAnglePhi.eps");
  canv->SetLogy(0);

  TH1D* hCrossingEnergy = new TH1D("CrossingEnergy",";True Total Energy (GeV);",100,0,50000);
  TH1D* hNonCrossingEnergy = new TH1D("NonCrossingEnergy",";True Total Energy (GeV);",100,0,50000);
  tree->Draw("MuonEnergy>>CrossingEnergy","MuonInnerCross");
  tree->Draw("MuonEnergy>>NonCrossingEnergy","!MuonInnerCross");
  canv->cd();
  canv->Clear();
  canv->SetLogy();
  hCrossingEnergy->SetLineColor(1);
  hCrossingEnergy->SetFillColor(1);
  hCrossingEnergy->SetFillStyle(3003);
  //hCrossingEnergy->Scale(1./hCrossingEnergy->GetEntries());
  hNonCrossingEnergy->SetLineColor(2);
  hNonCrossingEnergy->SetFillColor(2);
  hNonCrossingEnergy->SetFillStyle(3003);
  //hNonCrossingEnergy->Scale(1./hNonCrossingEnergy->GetEntries());
  hNonCrossingEnergy->Draw();
  hCrossingEnergy->Draw("same");
  leg->Clear();
  leg->AddEntry(hCrossingEnergy,    "Crossing",     "f");
  leg->AddEntry(hNonCrossingEnergy, "Non crossing", "f");
  leg->Draw();
  canv->SaveAs("CrossingMuonsEnergy.eps");
  canv->SetLogy(0);

  TH1D* hCrossingAngleTheta = new TH1D("CrossingAngleTheta",";True Angle Theta (deg);",100,0,180);
  TH1D* hNonCrossingAngleTheta = new TH1D("NonCrossingAngleTheta",";True Angle Theta (deg);",100,0,180);
  tree->Draw("MuonAngleTheta>>CrossingAngleTheta","MuonInnerCross");
  tree->Draw("MuonAngleTheta>>NonCrossingAngleTheta","!MuonInnerCross");
  canv->cd();
  canv->Clear();
  canv->SetLogy();
  hCrossingAngleTheta->SetLineColor(1);
  hCrossingAngleTheta->SetFillColor(1);
  hCrossingAngleTheta->SetFillStyle(3003);
  //hCrossingAngleTheta->Scale(1./hCrossingAngleTheta->GetEntries());
  hNonCrossingAngleTheta->SetLineColor(2);
  hNonCrossingAngleTheta->SetFillColor(2);
  hNonCrossingAngleTheta->SetFillStyle(3003);
  //hNonCrossingAngleTheta->Scale(1./hNonCrossingAngleTheta->GetEntries());
  hNonCrossingAngleTheta->Draw();
  hCrossingAngleTheta->Draw("same");
  leg->Clear();
  leg->AddEntry(hCrossingAngleTheta,    "Crossing",     "f");
  leg->AddEntry(hNonCrossingAngleTheta, "Non crossing", "f");
  leg->Draw();
  canv->SaveAs("CrossingMuonsAngleTheta.eps");
  canv->SetLogy(0);

  TH1D* hCrossingAnglePhi = new TH1D("CrossingAnglePhi",";True Angle Phi (deg);",100,0,180);
  TH1D* hNonCrossingAnglePhi = new TH1D("NonCrossingAnglePhi",";True Angle Phi (deg);",100,0,180);
  tree->Draw("MuonAnglePhi>>CrossingAnglePhi","MuonInnerCross");
  tree->Draw("MuonAnglePhi>>NonCrossingAnglePhi","!MuonInnerCross");
  canv->cd();
  canv->Clear();
  canv->SetLogy();
  hCrossingAnglePhi->SetLineColor(1);
  hCrossingAnglePhi->SetFillColor(1);
  hCrossingAnglePhi->SetFillStyle(3003);
  //hCrossingAnglePhi->Scale(1./hCrossingAnglePhi->GetEntries());
  hNonCrossingAnglePhi->SetLineColor(2);
  hNonCrossingAnglePhi->SetFillColor(2);
  hNonCrossingAnglePhi->SetFillStyle(3003);
  //hNonCrossingAnglePhi->Scale(1./hNonCrossingAnglePhi->GetEntries());
  hNonCrossingAnglePhi->Draw();
  hCrossingAnglePhi->Draw("same");
  leg->Clear();
  leg->AddEntry(hCrossingAnglePhi,    "Crossing",     "f");
  leg->AddEntry(hNonCrossingAnglePhi, "Non crossing", "f");
  leg->Draw();
  canv->SaveAs("CrossingMuonsAnglePhi.eps");
  canv->SetLogy(0);

}
