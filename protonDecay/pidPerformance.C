#include "TMath.h"

#include <string>
#include <map>

void pidPerformance() {

  gStyle->SetOptStat(0);
  const int kMaxParticles = 1000;
  const int kMaxOtherParticles = 10;

  TFile* inFile = new TFile("/pnfs/dune/scratch/users/wallbank/v06_45_01/reconana/pdkcosmobg/anahist.root");
  //TFile* inFile = new TFile("/dune/app/users/wallbank/larsoft-pdk/workspace/protonDecay/pdk_recon_ana.root");
  TTree* inTree = (TTree*)inFile->Get("pdkreconana/ReconPerformance");
  if (!inTree) {
    std::cout << "Input tree not found." << std::endl;
    return;
  }

  std::map<int,std::string> particlePDGs;
  particlePDGs[2212] = "Proton";
  particlePDGs[13]   = "Muon";
  particlePDGs[321]  = "Kaon";
  particlePDGs[211]  = "Pion";

  int fNKaons;
  int fNMuons;
  int fNMuons20cm;
  int fNHits2cm;

  int fNParticles;
  int fPDG[kMaxParticles];
  float fTrueEnergy[kMaxParticles];
  float fDepositedEnergy[kMaxParticles];
  int fNParticleTracks[kMaxParticles];
  float fEfficiency[kMaxParticles];
  float fPurity[kMaxParticles];
  int fReconHits[kMaxParticles][3];
  int fReconPDG[kMaxParticles][3];
  float fPIDA[kMaxParticles][3];

  int fNTracks;
  int fTrackPDG[kMaxParticles];
  float fTrackTrueEnergy[kMaxParticles];
  float fTrackDepositedEnergy[kMaxParticles];
  int fTrackReconHits[kMaxParticles][3];
  int fTrackReconPDG[kMaxParticles][3];
  float fTrackPIDA[kMaxParticles][3];
  int fOtherParticles[kMaxParticles];
  int fOtherParticlePDG[kMaxParticles][kMaxOtherParticles];
  float fOtherParticleContamination[kMaxParticles][kMaxOtherParticles];
  float fOtherParticleFraction[kMaxParticles][kMaxOtherParticles];

  inTree->SetBranchAddress("NKaons",     &fNKaons);
  inTree->SetBranchAddress("NMuons",     &fNMuons);
  inTree->SetBranchAddress("NMuons20cm", &fNMuons20cm);
  inTree->SetBranchAddress("NHits2cm",   &fNHits2cm);

  inTree->SetBranchAddress("NParticles",      &fNParticles);
  inTree->SetBranchAddress("PDG",             fPDG);
  inTree->SetBranchAddress("TrueEnergy",      fTrueEnergy);
  inTree->SetBranchAddress("DepositedEnergy", fDepositedEnergy);
  inTree->SetBranchAddress("NParticleTracks", fNParticleTracks);
  inTree->SetBranchAddress("Efficiency",      fEfficiency);
  inTree->SetBranchAddress("Purity",          fPurity);
  inTree->SetBranchAddress("ReconHits",       fReconHits);
  inTree->SetBranchAddress("ReconPDG",        fReconPDG);
  inTree->SetBranchAddress("PIDA",            fPIDA);

  inTree->SetBranchAddress("NTracks",                    &fNTracks);
  inTree->SetBranchAddress("TrackPDG",                   fTrackPDG);
  inTree->SetBranchAddress("TrackTrueEnergy",            fTrackTrueEnergy);
  inTree->SetBranchAddress("TrackDepositedEnergy",       fTrackDepositedEnergy);
  inTree->SetBranchAddress("TrackReconHits",             fTrackReconHits);
  inTree->SetBranchAddress("TrackReconPDG",              fTrackReconPDG);
  inTree->SetBranchAddress("TrackPIDA",                  fTrackPIDA);
  inTree->SetBranchAddress("OtherParticles",             fOtherParticles);
  inTree->SetBranchAddress("OtherParticlePDG",           fOtherParticlePDG);
  inTree->SetBranchAddress("OtherParticleContamination", fOtherParticleContamination);
  inTree->SetBranchAddress("OtherParticleFraction",      fOtherParticleFraction);

  TH1D* hTrueProtonEnergy  = new TH1D("TrueProtonEnergy", ";True Deposited Energy (GeV);True PID",100,0,3);
  TH1D* hTrueMuonEnergy    = new TH1D("TrueMuonEnergy",   ";True Deposited Energy (GeV);True PID",100,0,3);
  TH1D* hTrueKaonEnergy    = new TH1D("TrueKaonEnergy",   ";True Deposited Energy (GeV);True PID",100,0,3);
  TH1D* hTruePionEnergy    = new TH1D("TruePionEnergy",   ";True Deposited Energy (GeV);True PID",100,0,3);
  TH1D* hReconProtonEnergy = new TH1D("ReconProtonEnergy",";True Deposited Energy (GeV);Recon PID",100,0,3);
  TH1D* hReconMuonEnergy   = new TH1D("ReconMuonEnergy",  ";True Deposited Energy (GeV);Recon PID",100,0,3);
  TH1D* hReconKaonEnergy   = new TH1D("ReconKaonEnergy",  ";True Deposited Energy (GeV);Recon PID",100,0,3);
  TH1D* hReconPionEnergy   = new TH1D("ReconPionEnergy",  ";True Deposited Energy (GeV);Recon PID",100,0,3);

  TH2I* hParticleEnergy = new TH2I("ParticleEnergy",";True Deposited Energy (GeV);",100,0,3,5,0,5);
  hParticleEnergy->GetYaxis()->SetBinLabel(1,"Proton");
  hParticleEnergy->GetYaxis()->SetBinLabel(2,"Muon");
  hParticleEnergy->GetYaxis()->SetBinLabel(3,"Kaon");
  hParticleEnergy->GetYaxis()->SetBinLabel(4,"Pion");
  hParticleEnergy->GetYaxis()->SetBinLabel(5,"Other");

  TH1D* hProtonEfficiency  = new TH1D("ProtonEfficiency", ";Reconstruction Efficiency;",50,0,1.5);
  TH1D* hProtonPurity      = new TH1D("ProtonPurity",     ";Reconstruction Purity;",50,0,1.5);
  TH1D* hMuonEfficiency    = new TH1D("MuonEfficiency",   ";Reconstruction Efficiency;",50,0,1.5);
  TH1D* hMuonPurity        = new TH1D("MuonPurity",       ";Reconstruction Purity;",50,0,1.5);
  TH1D* hKaonEfficiency    = new TH1D("KaonEfficiency",   ";Reconstruction Efficiency;",100,0,1.5);
  TH1D* hKaonPurity        = new TH1D("KaonPurity",       ";Reconstruction Purity;",100,0,1.5);
  TH1D* hPionEfficiency    = new TH1D("PionEfficiency",   ";Reconstruction Efficiency;",50,0,1.5);
  TH1D* hPionPurity        = new TH1D("PionPurity",       ";Reconstruction Purity;",50,0,1.5);
  TH2I* hMisReconstruction = new TH2I("MisReconstruction",";True particle;Reconstructed with;",5,0,5,5,0,5);
  hMisReconstruction->GetXaxis()->SetBinLabel(1,"Proton");
  hMisReconstruction->GetXaxis()->SetBinLabel(2,"Muon");
  hMisReconstruction->GetXaxis()->SetBinLabel(3,"Kaon");
  hMisReconstruction->GetXaxis()->SetBinLabel(4,"Pion");
  hMisReconstruction->GetXaxis()->SetBinLabel(5,"Other");
  hMisReconstruction->GetYaxis()->SetBinLabel(1,"Proton");
  hMisReconstruction->GetYaxis()->SetBinLabel(2,"Muon");
  hMisReconstruction->GetYaxis()->SetBinLabel(3,"Kaon");
  hMisReconstruction->GetYaxis()->SetBinLabel(4,"Pion");
  hMisReconstruction->GetYaxis()->SetBinLabel(5,"Other");
  TH2I* hMisIdentification = new TH2I("MisIdentification",";True particle;Identified as;",5,0,5,5,0,5);
  hMisIdentification->GetXaxis()->SetBinLabel(1,"Proton");
  hMisIdentification->GetXaxis()->SetBinLabel(2,"Muon");
  hMisIdentification->GetXaxis()->SetBinLabel(3,"Kaon");
  hMisIdentification->GetXaxis()->SetBinLabel(4,"Pion");
  hMisIdentification->GetXaxis()->SetBinLabel(5,"Other");
  hMisIdentification->GetYaxis()->SetBinLabel(1,"Proton");
  hMisIdentification->GetYaxis()->SetBinLabel(2,"Muon");
  hMisIdentification->GetYaxis()->SetBinLabel(3,"Kaon");
  hMisIdentification->GetYaxis()->SetBinLabel(4,"Pion");
  hMisIdentification->GetYaxis()->SetBinLabel(5,"Other");

  TH1D* hProtonPIDA = new TH1D("ProtonPIDA",";PIDA;",100,0,40);
  TH1D* hMuonPIDA   = new TH1D("MuonPIDA",  ";PIDA;",100,0,40);
  TH1D* hKaonPIDA   = new TH1D("KaonPIDA",  ";PIDA;",100,0,40);
  TH1D* hPionPIDA   = new TH1D("PionPIDA",  ";PIDA;",100,0,40);

  TH1D* hProtonPIDAPlane = new TH1D("ProtonPIDAPlane",";PIDA;",100,0,40);
  TH1D* hMuonPIDAPlane   = new TH1D("MuonPIDAPlane",";PIDA;",100,0,40);
  TH1D* hKaonPIDAPlane   = new TH1D("KaonPIDAPlane",";PIDA;",100,0,40);
  TH1D* hPionPIDAPlane   = new TH1D("PionPIDAPlane",";PIDA;",100,0,40);

  TH1D* hProtonProtonPID = new TH1D("ProtonProtonPID",";True Deposited Energy (GeV);",100,0,3);
  TH1D* hProtonMuonPID   = new TH1D("ProtonMuonPID",  ";True Deposited Energy (GeV);",100,0,3);
  TH1D* hProtonKaonPID   = new TH1D("ProtonKaonPID",  ";True Deposited Energy (GeV);",100,0,3);
  TH1D* hProtonPionPID   = new TH1D("ProtonPionPID",  ";True Deposited Energy (GeV);",100,0,3);
  TH1D* hProtonOtherPID  = new TH1D("ProtonOtherPID", ";True Deposited Energy (GeV);",100,0,3);
  TH1D* hMuonProtonPID   = new TH1D("MuonProtonPID",  ";True Deposited Energy (GeV);",100,0,3);
  TH1D* hMuonMuonPID     = new TH1D("MuonMuonPID",    ";True Deposited Energy (GeV);",100,0,3);
  TH1D* hMuonKaonPID     = new TH1D("MuonKaonPID",    ";True Deposited Energy (GeV);",100,0,3);
  TH1D* hMuonPionPID     = new TH1D("MuonPionPID",    ";True Deposited Energy (GeV);",100,0,3);
  TH1D* hMuonOtherPID    = new TH1D("MuonOtherPID",   ";True Deposited Energy (GeV);",100,0,3);
  TH1D* hKaonProtonPID   = new TH1D("KaonProtonPID",  ";True Deposited Energy (GeV);",100,0,3);
  TH1D* hKaonMuonPID     = new TH1D("KaonMuonPID",    ";True Deposited Energy (GeV);",100,0,3);
  TH1D* hKaonKaonPID     = new TH1D("KaonKaonPID",    ";True Deposited Energy (GeV);",100,0,3);
  TH1D* hKaonPionPID     = new TH1D("KaonPionPID",    ";True Deposited Energy (GeV);",100,0,3);
  TH1D* hKaonOtherPID    = new TH1D("KaonOtherPID",   ";True Deposited Energy (GeV);",100,0,3);
  TH1D* hPionProtonPID   = new TH1D("PionProtonPID",  ";True Deposited Energy (GeV);",100,0,3);
  TH1D* hPionMuonPID     = new TH1D("PionMuonPID",    ";True Deposited Energy (GeV);",100,0,3);
  TH1D* hPionKaonPID     = new TH1D("PionKaonPID",    ";True Deposited Energy (GeV);",100,0,3);
  TH1D* hPionPionPID     = new TH1D("PionPionPID",    ";True Deposited Energy (GeV);",100,0,3);
  TH1D* hPionOtherPID    = new TH1D("PionOtherPID",   ";True Deposited Energy (GeV);",100,0,3);

  TH1D* hProtonTrackProtonPID = new TH1D("ProtonTrackProtonPID",";True Deposited Energy (GeV);",100,0,3);
  TH1D* hProtonTrackMuonPID   = new TH1D("ProtonTrackMuonPID",  ";True Deposited Energy (GeV);",100,0,3);
  TH1D* hProtonTrackKaonPID   = new TH1D("ProtonTrackKaonPID",  ";True Deposited Energy (GeV);",100,0,3);
  TH1D* hProtonTrackPionPID   = new TH1D("ProtonTrackPionPID",  ";True Deposited Energy (GeV);",100,0,3);
  TH1D* hProtonTrackOtherPID  = new TH1D("ProtonTrackOtherPID", ";True Deposited Energy (GeV);",100,0,3);
  TH1D* hMuonTrackProtonPID   = new TH1D("MuonTrackProtonPID",  ";True Deposited Energy (GeV);",100,0,3);
  TH1D* hMuonTrackMuonPID     = new TH1D("MuonTrackMuonPID",    ";True Deposited Energy (GeV);",100,0,3);
  TH1D* hMuonTrackKaonPID     = new TH1D("MuonTrackKaonPID",    ";True Deposited Energy (GeV);",100,0,3);
  TH1D* hMuonTrackPionPID     = new TH1D("MuonTrackPionPID",    ";True Deposited Energy (GeV);",100,0,3);
  TH1D* hMuonTrackOtherPID    = new TH1D("MuonTrackOtherPID",   ";True Deposited Energy (GeV);",100,0,3);
  TH1D* hKaonTrackProtonPID   = new TH1D("KaonTrackProtonPID",  ";True Deposited Energy (GeV);",100,0,3);
  TH1D* hKaonTrackMuonPID     = new TH1D("KaonTrackMuonPID",    ";True Deposited Energy (GeV);",100,0,3);
  TH1D* hKaonTrackKaonPID     = new TH1D("KaonTrackKaonPID",    ";True Deposited Energy (GeV);",100,0,3);
  TH1D* hKaonTrackPionPID     = new TH1D("KaonTrackPionPID",    ";True Deposited Energy (GeV);",100,0,3);
  TH1D* hKaonTrackOtherPID    = new TH1D("KaonTrackOtherPID",   ";True Deposited Energy (GeV);",100,0,3);
  TH1D* hPionTrackProtonPID   = new TH1D("PionTrackProtonPID",  ";True Deposited Energy (GeV);",100,0,3);
  TH1D* hPionTrackMuonPID     = new TH1D("PionTrackMuonPID",    ";True Deposited Energy (GeV);",100,0,3);
  TH1D* hPionTrackKaonPID     = new TH1D("PionTrackKaonPID",    ";True Deposited Energy (GeV);",100,0,3);
  TH1D* hPionTrackPionPID     = new TH1D("PionTrackPionPID",    ";True Deposited Energy (GeV);",100,0,3);
  TH1D* hPionTrackOtherPID    = new TH1D("PionTrackOtherPID",   ";True Deposited Energy (GeV);",100,0,3);

  for (int event = 0; event < inTree->GetEntriesFast(); ++event) {

    if (event % 1000 == 0)
      std::cout << "Processing event " << event << std::endl;

    inTree->GetEntry(event);

    for (int track = 0; track < fNTracks; ++track) {
      int best_plane = -1;
      int largest_nhits = 0;
      for (int plane = 0; plane < 3; ++plane) {
  	if (fTrackReconHits[track][plane] > largest_nhits) {
  	  best_plane = plane;
  	  largest_nhits = fTrackReconHits[track][plane];
  	}
  	if (fTrackPIDA[track][plane] == 0)
  	  continue;
  	if (fTrackPDG[track] == 2212) {
  	  hTrueProtonEnergy->Fill(fTrackDepositedEnergy[track]);
  	  hProtonPIDA->Fill(fTrackPIDA[track][plane]);
  	}
  	else if (TMath::Abs(fTrackPDG[track]) == 13) {
  	  hTrueMuonEnergy->Fill(fTrackDepositedEnergy[track]);
  	  hMuonPIDA->Fill(fTrackPIDA[track][plane]);
  	}
  	else if (fTrackPDG[track] == 321) {
  	  hTrueKaonEnergy->Fill(fTrackDepositedEnergy[track]);
  	  hKaonPIDA->Fill(fTrackPIDA[track][plane]);
  	}
  	else if (TMath::Abs(fTrackPDG[track]) == 211) {
  	  hTruePionEnergy->Fill(fTrackDepositedEnergy[track]);
	  hPionPIDA->Fill(fTrackPIDA[track][plane]);
	}
      } // plane
      if (best_plane > -1) {
  	if (fTrackPDG[track] == 2212) {
  	  hProtonPIDAPlane->Fill(fTrackPIDA[track][best_plane]);
  	  if (fTrackReconPDG[track][best_plane] == 2212) {
  	    hProtonTrackProtonPID->Fill(fTrackDepositedEnergy[track]);
	    //hMisIdentification->Fill("Proton","Proton",1);
	  }
  	  else if (TMath::Abs(fTrackReconPDG[track][best_plane]) == 13) {
  	    hProtonTrackMuonPID->Fill(fTrackDepositedEnergy[track]);
	    hMisIdentification->Fill("Proton","Muon",1);
	  }
  	  else if (fTrackReconPDG[track][best_plane] == 321) {
  	    hProtonTrackKaonPID->Fill(fTrackDepositedEnergy[track]);
	    hMisIdentification->Fill("Proton","Kaon",1);
	  }
	  else if (TMath::Abs(fTrackReconPDG[track][best_plane]) == 211) {
	    hProtonTrackPionPID->Fill(fTrackDepositedEnergy[track]);
	    hMisIdentification->Fill("Proton","Pion",1);
	  }
	  else {
  	    hProtonTrackOtherPID->Fill(fTrackDepositedEnergy[track]);
	    hMisIdentification->Fill("Proton","Other",1);
	  }
  	}
  	else if (TMath::Abs(fTrackPDG[track]) == 13) {
  	  hMuonPIDAPlane->Fill(fTrackPIDA[track][best_plane]);
  	  if (fTrackReconPDG[track][best_plane] == 2212) {
  	    hMuonTrackProtonPID->Fill(fTrackDepositedEnergy[track]);
	    hMisIdentification->Fill("Muon","Proton",1);
	  }
  	  else if (TMath::Abs(fTrackReconPDG[track][best_plane]) == 13) {
  	    hMuonTrackMuonPID->Fill(fTrackDepositedEnergy[track]);
	    //hMisIdentification->Fill("Muon","Muon",1);
	  }
  	  else if (fTrackReconPDG[track][best_plane] == 321) {
  	    hMuonTrackKaonPID->Fill(fTrackDepositedEnergy[track]);
	    hMisIdentification->Fill("Muon","Kaon",1);
	  }
	  else if (TMath::Abs(fTrackReconPDG[track][best_plane]) == 211) {
	    hMuonTrackPionPID->Fill(fTrackDepositedEnergy[track]);
	    hMisIdentification->Fill("Muon","Pion",1);
	  }
  	  else {
  	    hMuonTrackOtherPID->Fill(fTrackDepositedEnergy[track]);
	    hMisIdentification->Fill("Muon","Other",1);
	  }
  	}
  	else if (fTrackPDG[track] == 321) {
  	  hKaonPIDAPlane->Fill(fTrackPIDA[track][best_plane]);
  	  if (fTrackReconPDG[track][best_plane] == 2212) {
  	    hKaonTrackProtonPID->Fill(fTrackDepositedEnergy[track]);
	    hMisIdentification->Fill("Kaon","Proton",1);
	  }
  	  else if (TMath::Abs(fTrackReconPDG[track][best_plane]) == 13) {
  	    hKaonTrackMuonPID->Fill(fTrackDepositedEnergy[track]);
	    hMisIdentification->Fill("Kaon","Muon",1);
	  }
  	  else if (fTrackReconPDG[track][best_plane] == 321) {
  	    hKaonTrackKaonPID->Fill(fTrackDepositedEnergy[track]);
	    //hMisIdentification->Fill("Kaon","Kaon",1);
	  }
	  else if (TMath::Abs(fTrackReconPDG[track][best_plane]) == 211) {
	    hKaonTrackPionPID->Fill(fTrackDepositedEnergy[track]);
	    hMisIdentification->Fill("Kaon","Pion",1);
	  }
  	  else {
  	    hKaonTrackOtherPID->Fill(fTrackDepositedEnergy[track]);
	    hMisIdentification->Fill("Kaon","Other",1);
	  }
  	}
	else if (TMath::Abs(fTrackPDG[track]) == 211) {
  	  hPionPIDAPlane->Fill(fTrackPIDA[track][best_plane]);
  	  if (fTrackReconPDG[track][best_plane] == 2212) {
  	    hPionTrackProtonPID->Fill(fTrackDepositedEnergy[track]);
	    hMisIdentification->Fill("Pion","Proton",1);
	  }
  	  else if (TMath::Abs(fTrackReconPDG[track][best_plane]) == 13) {
  	    hPionTrackMuonPID->Fill(fTrackDepositedEnergy[track]);
	    hMisIdentification->Fill("Pion","Muon",1);
	  }
  	  else if (fTrackReconPDG[track][best_plane] == 321) {
  	    hPionTrackKaonPID->Fill(fTrackDepositedEnergy[track]);
	    hMisIdentification->Fill("Pion","Kaon",1);
	  }
	  else if (TMath::Abs(fTrackReconPDG[track][best_plane]) == 211) {
	    hPionTrackPionPID->Fill(fTrackDepositedEnergy[track]);
	    //hMisIdentification->Fill("Pion","Pion",1);
	  }
  	  else {
  	    hPionTrackOtherPID->Fill(fTrackDepositedEnergy[track]);
	    hMisIdentification->Fill("Pion","Other",1);
	  }
  	}
  	if (fTrackReconPDG[track][best_plane] == 2212)
  	  hReconProtonEnergy->Fill(fTrackDepositedEnergy[track]);
  	else if (TMath::Abs(fTrackReconPDG[track][best_plane]) == 13)
  	  hReconMuonEnergy->Fill(fTrackDepositedEnergy[track]);
  	else if (fTrackReconPDG[track][best_plane] == 321)
  	  hReconKaonEnergy->Fill(fTrackDepositedEnergy[track]);
  	else if (TMath::Abs(fTrackReconPDG[track][best_plane]) == 211)
  	  hReconPionEnergy->Fill(fTrackDepositedEnergy[track]);
      } // best_plane
      for (int otherParticle = 0; otherParticle < fOtherParticles[track]; ++otherParticle) {
      	int trackPDG = TMath::Abs(fTrackPDG[track]);
      	int otherParticlePDG = TMath::Abs(fOtherParticlePDG[track][otherParticle]);
      	if (otherParticlePDG == 2212 or otherParticlePDG == 13 or otherParticlePDG == 321 or otherParticlePDG == 211) {
      	  if (trackPDG == 2212 or trackPDG == 13 or trackPDG == 321 or trackPDG == 211)
      	    hMisReconstruction->Fill(particlePDGs[otherParticlePDG].c_str(),particlePDGs[trackPDG].c_str(),1);
      	  else
      	    hMisReconstruction->Fill(particlePDGs[otherParticlePDG].c_str(),"Other",1);
      	}
      	else {
      	  if (trackPDG == 2212 or trackPDG == 13 or trackPDG == 321 or trackPDG == 211)
      	    hMisReconstruction->Fill("Other",particlePDGs[trackPDG].c_str(),1);
      	  else
      	    hMisReconstruction->Fill("Other","Other",1);
      	}
      } // other particles
    } // track

    for (int particle = 0; particle < fNParticles; ++particle) {
      int best_plane = -1;
      int largest_nhits = 0;
      for (int plane = 0; plane < 3; ++plane) {
      	if (fReconHits[particle][plane] > largest_nhits) {
      	  best_plane = plane;
      	  largest_nhits = fReconHits[particle][plane];
      	}
      }
      if (best_plane < 0)
      	continue;
      if (fPDG[particle] == 2212) {
	hParticleEnergy->Fill(fDepositedEnergy[particle],"Proton",1);
	if (fEfficiency[particle] != 0)
	  hProtonEfficiency->Fill(fEfficiency[particle]);
	if (fPurity[particle] != 0)
	  hProtonPurity->Fill(fPurity[particle]);
    	if (fReconPDG[particle][best_plane] == 2212)
    	  hProtonProtonPID->Fill(fDepositedEnergy[particle]);
    	else if (TMath::Abs(fReconPDG[particle][best_plane]) == 13)
    	  hProtonMuonPID->Fill(fDepositedEnergy[particle]);
    	else if (fReconPDG[particle][best_plane] == 321)
    	  hProtonKaonPID->Fill(fDepositedEnergy[particle]);
	else if (TMath::Abs(fReconPDG[particle][best_plane]) == 211)
	  hProtonPionPID->Fill(fDepositedEnergy[particle]);
    	else
    	  hProtonOtherPID->Fill(fDepositedEnergy[particle]);
      }
      else if (TMath::Abs(fPDG[particle]) == 13) {
	hParticleEnergy->Fill(fDepositedEnergy[particle],"Muon",1);
	if (fEfficiency[particle] != 0)
	  hMuonEfficiency->Fill(fEfficiency[particle]);
	if (fPurity[particle] != 0)
	  hMuonPurity->Fill(fPurity[particle]);
    	if (fReconPDG[particle][best_plane] == 2212)
    	  hMuonProtonPID->Fill(fDepositedEnergy[particle]);
    	else if (TMath::Abs(fReconPDG[particle][best_plane]) == 13)
    	  hMuonMuonPID->Fill(fDepositedEnergy[particle]);
    	else if (fReconPDG[particle][best_plane] == 321)
    	  hMuonKaonPID->Fill(fDepositedEnergy[particle]);
	else if (TMath::Abs(fReconPDG[particle][best_plane]) == 211)
	  hMuonPionPID->Fill(fDepositedEnergy[particle]);
    	else
    	  hMuonOtherPID->Fill(fDepositedEnergy[particle]);
      }
      else if (fPDG[particle] == 321) {
	hParticleEnergy->Fill(fDepositedEnergy[particle],"Kaon",1);
	if (fEfficiency[particle] != 0)
	  hKaonEfficiency->Fill(fEfficiency[particle]);
	if (fPurity[particle] != 0)
	  hKaonPurity->Fill(fPurity[particle]);
    	if (fReconPDG[particle][best_plane] == 2212)
    	  hKaonProtonPID->Fill(fDepositedEnergy[particle]);
    	else if (TMath::Abs(fReconPDG[particle][best_plane]) == 13)
    	  hKaonMuonPID->Fill(fDepositedEnergy[particle]);
    	else if (fReconPDG[particle][best_plane] == 321)
    	  hKaonKaonPID->Fill(fDepositedEnergy[particle]);
	else if (TMath::Abs(fReconPDG[particle][best_plane]) == 211)
	  hKaonPionPID->Fill(fDepositedEnergy[particle]);
    	else
    	  hKaonOtherPID->Fill(fDepositedEnergy[particle]);
      }
      else if (TMath::Abs(fPDG[particle]) == 211) {
	hParticleEnergy->Fill(fDepositedEnergy[particle],"Pion",1);
  	if (fEfficiency[particle] != 0)
	  hPionEfficiency->Fill(fEfficiency[particle]);
	if (fPurity[particle] != 0)
	  hPionPurity->Fill(fPurity[particle]);
    	if (fReconPDG[particle][best_plane] == 2212)
    	  hPionProtonPID->Fill(fDepositedEnergy[particle]);
    	else if (TMath::Abs(fReconPDG[particle][best_plane]) == 13)
    	  hPionMuonPID->Fill(fDepositedEnergy[particle]);
    	else if (fReconPDG[particle][best_plane] == 321)
    	  hPionKaonPID->Fill(fDepositedEnergy[particle]);
	else if (TMath::Abs(fReconPDG[particle][best_plane]) == 211)
	  hPionPionPID->Fill(fDepositedEnergy[particle]);
    	else
    	  hPionOtherPID->Fill(fDepositedEnergy[particle]);
      }
      else
	hParticleEnergy->Fill(fDepositedEnergy[particle],"Other",1);
    } // particle

  } // event

  TCanvas* canv = new TCanvas("canv","",800,600);
  hProtonPIDA->Scale(1./hProtonPIDA->GetEntries());
  hProtonPIDA->SetLineColor(1);
  hProtonPIDA->SetFillColor(1);
  hProtonPIDA->SetFillStyle(3003);
  hMuonPIDA->Scale(1./hMuonPIDA->GetEntries());
  hMuonPIDA->SetLineColor(2);
  hMuonPIDA->SetFillColor(2);
  hMuonPIDA->SetFillStyle(3003);
  hKaonPIDA->Scale(1./hKaonPIDA->GetEntries());
  hKaonPIDA->SetLineColor(3);
  hKaonPIDA->SetFillColor(3);
  hKaonPIDA->SetFillStyle(3003);
  hPionPIDA->Scale(1./hPionPIDA->GetEntries());
  hPionPIDA->SetLineColor(4);
  hPionPIDA->SetFillColor(4);
  hPionPIDA->SetFillStyle(3003);
  hKaonPIDA->Draw("hist");
  hMuonPIDA->Draw("hist same");
  hProtonPIDA->Draw("hist same");
  hPionPIDA->Draw("hist same");
  TLegend* leg = new TLegend(0.6,0.7,0.9,0.9);
  leg->AddEntry(hProtonPIDA, "Proton", "f");
  leg->AddEntry(hMuonPIDA,   "Muon",   "f");
  leg->AddEntry(hKaonPIDA,   "Kaon",   "f");
  leg->AddEntry(hPionPIDA,   "Pion",   "f");
  leg->Draw();
  canv->SaveAs("PIDA.eps");

  canv->cd();
  canv->Clear();
  hProtonPIDAPlane->Scale(1./hProtonPIDAPlane->GetEntries());
  hProtonPIDAPlane->SetLineColor(1);
  hProtonPIDAPlane->SetFillColor(1);
  hProtonPIDAPlane->SetFillStyle(3003);
  hMuonPIDAPlane->Scale(1./hMuonPIDAPlane->GetEntries());
  hMuonPIDAPlane->SetLineColor(2);
  hMuonPIDAPlane->SetFillColor(2);
  hMuonPIDAPlane->SetFillStyle(3003);
  hKaonPIDAPlane->Scale(1./hKaonPIDAPlane->GetEntries());
  hKaonPIDAPlane->SetLineColor(3);
  hKaonPIDAPlane->SetFillColor(3);
  hKaonPIDAPlane->SetFillStyle(3003);
  hPionPIDAPlane->Scale(1./hPionPIDAPlane->GetEntries());
  hPionPIDAPlane->SetLineColor(4);
  hPionPIDAPlane->SetFillColor(4);
  hPionPIDAPlane->SetFillStyle(3003);
  hKaonPIDAPlane->Draw("hist");
  hMuonPIDAPlane->Draw("hist same");
  hProtonPIDAPlane->Draw("hist same");
  hPionPIDAPlane->Draw("hist same");
  leg->Clear();
  leg->AddEntry(hProtonPIDAPlane, "Proton", "f");
  leg->AddEntry(hMuonPIDAPlane,   "Muon",   "f");
  leg->AddEntry(hKaonPIDAPlane,   "Kaon",   "f");
  leg->AddEntry(hPionPIDAPlane,   "Pion",   "f");
  leg->Draw();
  canv->SaveAs("PIDAPlane.eps");

  canv->cd();
  canv->Clear();
  canv->SetLogy();
  hProtonProtonPID->SetLineColor(1);
  hProtonProtonPID->SetFillColor(1);
  hProtonProtonPID->SetFillStyle(3003);
  hProtonMuonPID->SetLineColor(2);
  hProtonMuonPID->SetFillColor(2);
  hProtonMuonPID->SetFillStyle(3003);
  hProtonKaonPID->SetLineColor(3);
  hProtonKaonPID->SetFillColor(3);
  hProtonKaonPID->SetFillStyle(3003);
  hProtonPionPID->SetLineColor(4);
  hProtonPionPID->SetFillColor(4);
  hProtonPionPID->SetFillStyle(3003);
  hProtonOtherPID->SetLineColor(5);
  hProtonOtherPID->SetFillColor(5);
  hProtonOtherPID->SetFillStyle(3003);
  hProtonProtonPID->Draw();
  hProtonMuonPID->Draw("same");
  hProtonKaonPID->Draw("same");
  hProtonPionPID->Draw("same");
  hProtonOtherPID->Draw("same");
  leg->Clear();
  leg->AddEntry(hProtonProtonPID, "Proton", "f");
  leg->AddEntry(hProtonMuonPID,   "Muon",   "f");
  leg->AddEntry(hProtonKaonPID,   "Kaon",   "f");
  leg->AddEntry(hProtonPionPID,   "Pion",   "f");
  leg->AddEntry(hProtonOtherPID,  "Other",  "f");
  leg->Draw();
  canv->SaveAs("ProtonPID.eps");
  canv->SetLogy(0);

  canv->cd();
  canv->Clear();
  canv->SetLogy();
  hMuonProtonPID->SetLineColor(1);
  hMuonProtonPID->SetFillColor(1);
  hMuonProtonPID->SetFillStyle(3003);
  hMuonMuonPID->SetLineColor(2);
  hMuonMuonPID->SetFillColor(2);
  hMuonMuonPID->SetFillStyle(3003);
  hMuonKaonPID->SetLineColor(3);
  hMuonKaonPID->SetFillColor(3);
  hMuonKaonPID->SetFillStyle(3003);
  hMuonPionPID->SetLineColor(4);
  hMuonPionPID->SetFillColor(4);
  hMuonPionPID->SetFillStyle(3003);
  hMuonOtherPID->SetLineColor(5);
  hMuonOtherPID->SetFillColor(5);
  hMuonOtherPID->SetFillStyle(3003);
  hMuonMuonPID->Draw();
  hMuonKaonPID->Draw("same");
  hMuonPionPID->Draw("same");
  hMuonOtherPID->Draw("same");
  hMuonProtonPID->Draw("same");
  leg->Clear();
  leg->AddEntry(hMuonProtonPID, "Proton", "f");
  leg->AddEntry(hMuonMuonPID,   "Muon",   "f");
  leg->AddEntry(hMuonKaonPID,   "Kaon",   "f");
  leg->AddEntry(hMuonPionPID,   "Pion",   "f");
  leg->AddEntry(hMuonOtherPID,  "Other",  "f");
  leg->Draw();
  canv->SaveAs("MuonPID.eps");
  canv->SetLogy(0);

  canv->cd();
  canv->Clear();
  hKaonProtonPID->SetLineColor(1);
  hKaonProtonPID->SetFillColor(1);
  hKaonProtonPID->SetFillStyle(3003);
  hKaonMuonPID->SetLineColor(2);
  hKaonMuonPID->SetFillColor(2);
  hKaonMuonPID->SetFillStyle(3003);
  hKaonKaonPID->SetLineColor(3);
  hKaonKaonPID->SetFillColor(3);
  hKaonKaonPID->SetFillStyle(3003);
  hKaonPionPID->SetLineColor(4);
  hKaonPionPID->SetFillColor(4);
  hKaonPionPID->SetFillStyle(3003);
  hKaonOtherPID->SetLineColor(5);
  hKaonOtherPID->SetFillColor(5);
  hKaonOtherPID->SetFillStyle(3003);
  hKaonProtonPID->Draw();
  hKaonMuonPID->Draw("same");
  hKaonKaonPID->Draw("same");
  hKaonPionPID->Draw("same");
  hKaonOtherPID->Draw("same");
  leg->Clear();
  leg->AddEntry(hKaonProtonPID, "Proton", "f");
  leg->AddEntry(hKaonMuonPID,   "Muon",   "f");
  leg->AddEntry(hKaonKaonPID,   "Kaon",   "f");
  leg->AddEntry(hKaonPionPID,   "Pion",   "f");
  leg->AddEntry(hKaonOtherPID,  "Other",  "f");
  leg->Draw();
  canv->SaveAs("KaonPID.eps");

  canv->cd();
  canv->Clear();
  canv->SetLogy();
  hPionProtonPID->SetLineColor(1);
  hPionProtonPID->SetFillColor(1);
  hPionProtonPID->SetFillStyle(3003);
  hPionMuonPID->SetLineColor(2);
  hPionMuonPID->SetFillColor(2);
  hPionMuonPID->SetFillStyle(3003);
  hPionKaonPID->SetLineColor(3);
  hPionKaonPID->SetFillColor(3);
  hPionKaonPID->SetFillStyle(3003);
  hPionPionPID->SetLineColor(4);
  hPionPionPID->SetFillColor(4);
  hPionPionPID->SetFillStyle(3003);
  hPionOtherPID->SetLineColor(5);
  hPionOtherPID->SetFillColor(5);
  hPionOtherPID->SetFillStyle(3003);
  hPionMuonPID->Draw();
  hPionPionPID->Draw("same");
  hPionKaonPID->Draw("same");
  hPionProtonPID->Draw("same");
  hPionOtherPID->Draw("same");
  leg->Clear();
  leg->AddEntry(hPionProtonPID, "Proton", "f");
  leg->AddEntry(hPionMuonPID,   "Muon",   "f");
  leg->AddEntry(hPionKaonPID,   "Kaon",   "f");
  leg->AddEntry(hPionPionPID,   "Pion",   "f");
  leg->AddEntry(hPionOtherPID,  "Other",  "f");
  leg->Draw();
  canv->SaveAs("PionPID.eps");
  canv->SetLogy(0);

  canv->cd();
  canv->Clear();
  canv->SetLogy();
  hProtonTrackProtonPID->SetLineColor(1);
  hProtonTrackProtonPID->SetFillColor(1);
  hProtonTrackProtonPID->SetFillStyle(3003);
  hProtonTrackMuonPID->SetLineColor(2);
  hProtonTrackMuonPID->SetFillColor(2);
  hProtonTrackMuonPID->SetFillStyle(3003);
  hProtonTrackKaonPID->SetLineColor(3);
  hProtonTrackKaonPID->SetFillColor(3);
  hProtonTrackKaonPID->SetFillStyle(3003);
  hProtonTrackPionPID->SetLineColor(4);
  hProtonTrackPionPID->SetFillColor(4);
  hProtonTrackPionPID->SetFillStyle(3003);
  hProtonTrackOtherPID->SetLineColor(5);
  hProtonTrackOtherPID->SetFillColor(5);
  hProtonTrackOtherPID->SetFillStyle(3003);
  hProtonTrackProtonPID->Draw();
  hProtonTrackMuonPID->Draw("same");
  hProtonTrackKaonPID->Draw("same");
  hProtonTrackPionPID->Draw("same");
  hProtonTrackOtherPID->Draw("same");
  leg->Clear();
  leg->AddEntry(hProtonTrackProtonPID, "Proton", "f");
  leg->AddEntry(hProtonTrackMuonPID,   "Muon",   "f");
  leg->AddEntry(hProtonTrackKaonPID,   "Kaon",   "f");
  leg->AddEntry(hProtonTrackPionPID,   "Pion",   "f");
  leg->AddEntry(hProtonTrackOtherPID,  "Other",  "f");
  leg->Draw();
  canv->SaveAs("ProtonTrackPID.eps");
  canv->SetLogy(0);

  canv->cd();
  canv->Clear();
  canv->SetLogy();
  hMuonTrackProtonPID->SetLineColor(1);
  hMuonTrackProtonPID->SetFillColor(1);
  hMuonTrackProtonPID->SetFillStyle(3003);
  hMuonTrackMuonPID->SetLineColor(2);
  hMuonTrackMuonPID->SetFillColor(2);
  hMuonTrackMuonPID->SetFillStyle(3003);
  hMuonTrackKaonPID->SetLineColor(3);
  hMuonTrackKaonPID->SetFillColor(3);
  hMuonTrackKaonPID->SetFillStyle(3003);
  hMuonTrackPionPID->SetLineColor(4);
  hMuonTrackPionPID->SetFillColor(4);
  hMuonTrackPionPID->SetFillStyle(3003);
  hMuonTrackOtherPID->SetLineColor(5);
  hMuonTrackOtherPID->SetFillColor(5);
  hMuonTrackOtherPID->SetFillStyle(3003);
  hMuonTrackMuonPID->Draw();
  hMuonTrackPionPID->Draw("same");
  hMuonTrackOtherPID->Draw("same");
  hMuonTrackKaonPID->Draw("same");
  hMuonTrackProtonPID->Draw("same");
  leg->Clear();
  leg->AddEntry(hMuonTrackProtonPID, "Proton", "f");
  leg->AddEntry(hMuonTrackMuonPID,   "Muon",   "f");
  leg->AddEntry(hMuonTrackKaonPID,   "Kaon",   "f");
  leg->AddEntry(hMuonTrackPionPID,   "Pion",   "f");
  leg->AddEntry(hMuonTrackOtherPID,  "Other",  "f");
  leg->Draw();
  canv->SaveAs("MuonTrackPID.eps");
  canv->SetLogy(0);

  canv->cd();
  canv->Clear();
  hKaonTrackProtonPID->SetLineColor(1);
  hKaonTrackProtonPID->SetFillColor(1);
  hKaonTrackProtonPID->SetFillStyle(3003);
  hKaonTrackMuonPID->SetLineColor(2);
  hKaonTrackMuonPID->SetFillColor(2);
  hKaonTrackMuonPID->SetFillStyle(3003);
  hKaonTrackKaonPID->SetLineColor(3);
  hKaonTrackKaonPID->SetFillColor(3);
  hKaonTrackKaonPID->SetFillStyle(3003);
  hKaonTrackPionPID->SetLineColor(4);
  hKaonTrackPionPID->SetFillColor(4);
  hKaonTrackPionPID->SetFillStyle(3003);
  hKaonTrackOtherPID->SetLineColor(5);
  hKaonTrackOtherPID->SetFillColor(5);
  hKaonTrackOtherPID->SetFillStyle(3003);
  hKaonTrackProtonPID->Draw();
  hKaonTrackMuonPID->Draw("same");
  hKaonTrackKaonPID->Draw("same");
  hKaonTrackPionPID->Draw("same");
  hKaonTrackOtherPID->Draw("same");
  leg->Clear();
  leg->AddEntry(hKaonTrackProtonPID, "Proton", "f");
  leg->AddEntry(hKaonTrackMuonPID,   "Muon",   "f");
  leg->AddEntry(hKaonTrackKaonPID,   "Kaon",   "f");
  leg->AddEntry(hKaonTrackPionPID,   "Pion",   "f");
  leg->AddEntry(hKaonTrackOtherPID,  "Other",  "f");
  leg->Draw();
  canv->SaveAs("KaonTrackPID.eps");

  canv->cd();
  canv->Clear();
  canv->SetLogy();
  hPionTrackProtonPID->SetLineColor(1);
  hPionTrackProtonPID->SetFillColor(1);
  hPionTrackProtonPID->SetFillStyle(3003);
  hPionTrackMuonPID->SetLineColor(2);
  hPionTrackMuonPID->SetFillColor(2);
  hPionTrackMuonPID->SetFillStyle(3003);
  hPionTrackKaonPID->SetLineColor(3);
  hPionTrackKaonPID->SetFillColor(3);
  hPionTrackKaonPID->SetFillStyle(3003);
  hPionTrackPionPID->SetLineColor(4);
  hPionTrackPionPID->SetFillColor(4);
  hPionTrackPionPID->SetFillStyle(3003);
  hPionTrackOtherPID->SetLineColor(5);
  hPionTrackOtherPID->SetFillColor(5);
  hPionTrackOtherPID->SetFillStyle(3003);
  hPionTrackMuonPID->Draw();
  hPionTrackPionPID->Draw("same");
  hPionTrackProtonPID->Draw("same");
  hPionTrackKaonPID->Draw("same");
  hPionTrackOtherPID->Draw("same");
  leg->Clear();
  leg->AddEntry(hPionTrackProtonPID, "Proton", "f");
  leg->AddEntry(hPionTrackMuonPID,   "Muon",   "f");
  leg->AddEntry(hPionTrackKaonPID,   "Kaon",   "f");
  leg->AddEntry(hPionTrackPionPID,   "Pion",   "f");
  leg->AddEntry(hPionTrackOtherPID,  "Other",  "f");
  leg->Draw();
  canv->SaveAs("PionTrackPID.eps");
  canv->SetLogy(0);

  canv->cd();
  canv->Clear();
  hMisReconstruction->Draw("colz");
  canv->SaveAs("MisReconstruction.eps");

  canv->cd();
  canv->Clear();
  canv->SetLogz();
  hMisIdentification->Draw("colz");
  canv->SaveAs("MisIdentification.eps");
  canv->SetLogz(0);

  canv->cd();
  canv->Clear();
  hProtonEfficiency->Scale(1./hProtonEfficiency->GetEntries());
  hMuonEfficiency->Scale(1./hMuonEfficiency->GetEntries());
  hKaonEfficiency->Scale(1./hKaonEfficiency->GetEntries());
  hPionEfficiency->Scale(1./hPionEfficiency->GetEntries());
  hProtonEfficiency->SetLineColor(1);
  hProtonEfficiency->SetFillColor(1);
  hProtonEfficiency->SetFillStyle(3003);
  hMuonEfficiency->SetLineColor(2);
  hMuonEfficiency->SetFillColor(2);
  hMuonEfficiency->SetFillStyle(3003);
  hKaonEfficiency->SetLineColor(3);
  hKaonEfficiency->SetFillColor(3);
  hKaonEfficiency->SetFillStyle(3003);
  hPionEfficiency->SetLineColor(4);
  hPionEfficiency->SetFillColor(4);
  hPionEfficiency->SetFillStyle(3003);
  hProtonEfficiency->Draw("hist");
  hMuonEfficiency->Draw("hist same");
  hKaonEfficiency->Draw("hist same");
  hPionEfficiency->Draw("hist same");
  TLegend* leftleg = new TLegend(0.15,0.6,0.45,0.85);
  leftleg->AddEntry(hProtonEfficiency, "Proton", "f");
  leftleg->AddEntry(hMuonEfficiency,   "Muon",   "f");
  leftleg->AddEntry(hKaonEfficiency,   "Kaon",   "f");
  leftleg->AddEntry(hPionEfficiency,   "Pion",   "f");
  leftleg->Draw();
  canv->SaveAs("ReconstructionEfficiency.eps");

  canv->cd();
  canv->Clear();
  hProtonPurity->Scale(1./hProtonPurity->GetEntries());
  hMuonPurity->Scale(1./hMuonPurity->GetEntries());
  hKaonPurity->Scale(1./hKaonPurity->GetEntries());
  hPionPurity->Scale(1./hPionPurity->GetEntries());
  hProtonPurity->SetLineColor(1);
  hProtonPurity->SetFillColor(1);
  hProtonPurity->SetFillStyle(3003);
  hMuonPurity->SetLineColor(2);
  hMuonPurity->SetFillColor(2);
  hMuonPurity->SetFillStyle(3003);
  hKaonPurity->SetLineColor(3);
  hKaonPurity->SetFillColor(3);
  hKaonPurity->SetFillStyle(3003);
  hPionPurity->SetLineColor(4);
  hPionPurity->SetFillColor(4);
  hPionPurity->SetFillStyle(3003);
  hProtonPurity->Draw("hist");
  hMuonPurity->Draw("hist same");
  hKaonPurity->Draw("hist same");
  hPionPurity->Draw("hist same");
  leftleg->Clear();
  leftleg->AddEntry(hProtonPurity, "Proton", "f");
  leftleg->AddEntry(hMuonPurity,   "Muon",   "f");
  leftleg->AddEntry(hKaonPurity,   "Kaon",   "f");
  leftleg->AddEntry(hPionPurity,   "Pion",   "f");
  leftleg->Draw();
  canv->SaveAs("ReconstructionPurity.eps");

  canv->cd();
  canv->Clear();
  canv->SetLogy();
  hTrueProtonEnergy->SetLineColor(1);
  hTrueProtonEnergy->SetFillColor(1);
  hTrueProtonEnergy->SetFillStyle(3003);
  hTrueMuonEnergy->SetLineColor(2);
  hTrueMuonEnergy->SetFillColor(2);
  hTrueMuonEnergy->SetFillStyle(3003);
  hTrueKaonEnergy->SetLineColor(3);
  hTrueKaonEnergy->SetFillColor(3);
  hTrueKaonEnergy->SetFillStyle(3003);
  hTruePionEnergy->SetLineColor(4);
  hTruePionEnergy->SetFillColor(4);
  hTruePionEnergy->SetFillStyle(3003);
  hTrueMuonEnergy->Draw();
  hTrueProtonEnergy->Draw("same");
  hTruePionEnergy->Draw("same");
  hTrueKaonEnergy->Draw("same");
  leg->Clear();
  leg->AddEntry(hTrueProtonEnergy, "Proton", "f");
  leg->AddEntry(hTrueMuonEnergy,   "Muon",   "f");
  leg->AddEntry(hTrueKaonEnergy,   "Kaon",   "f");
  leg->AddEntry(hTruePionEnergy,   "Pion",   "f");
  leg->Draw();
  canv->SaveAs("TrueEnergy.eps");
  canv->SetLogy(0);

  canv->cd();
  canv->Clear();
  canv->SetLogy();
  hReconProtonEnergy->SetLineColor(1);
  hReconProtonEnergy->SetFillColor(1);
  hReconProtonEnergy->SetFillStyle(3003);
  hReconMuonEnergy->SetLineColor(2);
  hReconMuonEnergy->SetFillColor(2);
  hReconMuonEnergy->SetFillStyle(3003);
  hReconKaonEnergy->SetLineColor(3);
  hReconKaonEnergy->SetFillColor(3);
  hReconKaonEnergy->SetFillStyle(3003);
  hReconPionEnergy->SetLineColor(4);
  hReconPionEnergy->SetFillColor(4);
  hReconPionEnergy->SetFillStyle(3003);
  hReconMuonEnergy->Draw();
  hReconPionEnergy->Draw("same");
  hReconProtonEnergy->Draw("same");
  hReconKaonEnergy->Draw("same");
  leg->Clear();
  leg->AddEntry(hReconProtonEnergy, "Proton", "f");
  leg->AddEntry(hReconMuonEnergy,   "Muon",   "f");
  leg->AddEntry(hReconKaonEnergy,   "Kaon",   "f");
  leg->AddEntry(hReconPionEnergy,   "Pion",   "f");
  leg->Draw();
  canv->SaveAs("ReconEnergy.eps");
  canv->SetLogy(0);

  canv->cd();
  canv->Clear();
  canv->SetLogz();
  hParticleEnergy->Draw("colz");
  canv->SaveAs("ParticleEnergy.eps");
  canv->SetLogz(0);

  return;

}
