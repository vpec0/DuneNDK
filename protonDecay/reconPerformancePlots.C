{

  const int kMaxParticles = 1000000;

  // TFile* inFile = new TFile("pdk_recon_ana.root");
  TFile* inFile = new TFile("/pnfs/dune/scratch/users/wallbank/v06_45_01/reconana/pdkcosmobg/anahist.root");
  TTree* inTree = (TTree*)inFile->Get("pdkreconana/ReconPerformance");

  TProfile* hKaonAll = new TProfile("KaonAll",";True Energy (MeV);Reconstruction efficiency;",100,0,100000);
  TProfile* hKaonKaon = new TProfile("KaonKaon",";True Energy (MeV);Reconstruction efficiency;",100,0,100000);
  TProfile* hKaonOneKaon = new TProfile("KaonOneKaon",";True Energy (MeV);Reconstruction efficiency;",100,0,100000);
  TProfile* hKaonNoMuon = new TProfile("KaonNoMuon",";True Energy (MeV);Reconstruction efficiency;",100,0,100000);
  TProfile* hKaonNoHits = new TProfile("KaonNoHits",";True Energy (MeV);Reconstruction efficiency;",100,0,100000);

  TProfile* hMuonAll = new TProfile("MuonAll",";True Energy (MeV);Reconstruction efficiency;",100,0,100000);
  TProfile* hMuonKaon = new TProfile("MuonKaon",";True Energy (MeV);Reconstruction efficiency;",100,0,100000);
  TProfile* hMuonOneKaon = new TProfile("MuonOneKaon",";True Energy (MeV);Reconstruction efficiency;",100,0,100000);
  TProfile* hMuonNoMuon = new TProfile("MuonNoMuon",";True Energy (MeV);Reconstruction efficiency;",100,0,100000);
  TProfile* hMuonNoHits = new TProfile("MuonNoHits",";True Energy (MeV);Reconstruction efficiency;",100,0,100000);

  TProfile* hProtonAll = new TProfile("ProtonAll",";True Energy (MeV);Reconstruction efficiency;",100,0,100000);
  TProfile* hProtonKaon = new TProfile("ProtonKaon",";True Energy (MeV);Reconstruction efficiency;",100,0,100000);
  TProfile* hProtonOneKaon = new TProfile("ProtonOneKaon",";True Energy (MeV);Reconstruction efficiency;",100,0,100000);
  TProfile* hProtonNoMuon = new TProfile("ProtonNoMuon",";True Energy (MeV);Reconstruction efficiency;",100,0,100000);
  TProfile* hProtonNoHits = new TProfile("ProtonNoHits",";True Energy (MeV);Reconstruction efficiency;",100,0,100000);

  int fNKaons;
  int fNMuons;
  int fNMuons20cm;
  int fNHits2cm;
  int fNParticles;
  int fPDG[kMaxParticles];
  double fTrueEnergy[kMaxParticles];
  int fNTracks[kMaxParticles];
  double fEfficiency[kMaxParticles];
  inTree->SetBranchAddress("NKaons",     &fNKaons);
  inTree->SetBranchAddress("NMuons",     &fNMuons);
  inTree->SetBranchAddress("NMuons20cm", &fNMuons20cm);
  inTree->SetBranchAddress("NHits2cm",   &fNHits2cm);
  inTree->SetBranchAddress("NParticles", &fNParticles);
  inTree->SetBranchAddress("PDG",        fPDG);
  inTree->SetBranchAddress("TrueEnergy", fTrueEnergy);
  inTree->SetBranchAddress("NTracks",    fNTracks);
  inTree->SetBranchAddress("Efficiency", fEfficiency);

  // look through events
  for (int event = 0; event < inTree->GetEntriesFast(); ++event) {

    if (event % 1000 == 0)
      std::cout << "Analysing event " << event << " out of " << inTree->GetEntriesFast() << std::endl;

    inTree->GetEntry(event);

    for (int particle = 0; particle < fNParticles; ++particle) {

      if (fPDG[particle] == 321) {
	hKaonAll->Fill(fTrueEnergy[particle], fEfficiency[particle]>0.7);
	if (fNKaons > 0) {
	  hKaonKaon->Fill(fTrueEnergy[particle], fEfficiency[particle]>0.7);
	  if (fNKaons == 1) {
	    hKaonOneKaon->Fill(fTrueEnergy[particle], fEfficiency[particle]>0.7);
	    if (fNMuons20cm == 0) {
	      hKaonNoMuon->Fill(fTrueEnergy[particle], fEfficiency[particle]>0.7);
	      if (fNHits2cm == 0) {
		hKaonNoHits->Fill(fTrueEnergy[particle], fEfficiency[particle]>0.7);
	      }
	    }
	  }
	}
      }

      if (fPDG[particle] == 11) {
	hMuonAll->Fill(fTrueEnergy[particle], fEfficiency[particle]>0.7);
	if (fNKaons > 0) {
	  hMuonKaon->Fill(fTrueEnergy[particle], fEfficiency[particle]>0.7);
	  if (fNKaons == 1) {
	    hMuonOneKaon->Fill(fTrueEnergy[particle], fEfficiency[particle]>0.7);
	    if (fNMuons20cm == 0) {
	      hMuonNoMuon->Fill(fTrueEnergy[particle], fEfficiency[particle]>0.7);
	      if (fNHits2cm == 0) {
		hMuonNoHits->Fill(fTrueEnergy[particle], fEfficiency[particle]>0.7);
	      }
	    }
	  }
	}
      }

      if (fPDG[particle] == 2212) {
	hProtonAll->Fill(fTrueEnergy[particle], fEfficiency[particle]>0.7);
	if (fNKaons > 0) {
	  hProtonKaon->Fill(fTrueEnergy[particle], fEfficiency[particle]>0.7);
	  if (fNKaons == 1) {
	    hProtonOneKaon->Fill(fTrueEnergy[particle], fEfficiency[particle]>0.7);
	    if (fNMuons20cm == 0) {
	      hProtonNoMuon->Fill(fTrueEnergy[particle], fEfficiency[particle]>0.7);
	      if (fNHits2cm == 0) {
		hProtonNoHits->Fill(fTrueEnergy[particle], fEfficiency[particle]>0.7);
	      }
	    }
	  }
	}
      }
    }

  }

  TCanvas* canv = new TCanvas("canv","",800,600);
  TLegend* leg;

  canv->cd();
  canv->Clear();
  leg = new TLegend(0.5,0.2,0.9,0.5);
  hKaonAll->SetLineColor(1);
  leg->AddEntry(hKaonAll,"No cuts","l");
  hKaonAll->Draw();
  hKaonKaon->SetLineColor(2);
  leg->AddEntry(hKaonKaon,"With Kaon track","l");
  hKaonKaon->Draw("same");
  hKaonOneKaon->SetLineColor(3);
  leg->AddEntry(hKaonOneKaon,"and One Kaon track in event","l");
  hKaonOneKaon->Draw("same");
  hKaonNoMuon->SetLineColor(4);
  leg->AddEntry(hKaonNoMuon,"and No Muon track >20cm","l");
  hKaonNoMuon->Draw("same");
  hKaonNoHits->SetLineColor(5);
  leg->AddEntry(hKaonNoMuon,"and No space points within 2cm of walls","l");
  hKaonNoHits->Draw("same");
  leg->Draw("same");
  canv->SaveAs("KaonReco.eps");
  delete leg;

  canv->cd();
  canv->Clear();
  leg = new TLegend(0.5,0.2,0.9,0.5);
  hMuonAll->SetLineColor(1);
  leg->AddEntry(hMuonAll,"No cuts","l");
  hMuonAll->Draw();
  hMuonKaon->SetLineColor(2);
  leg->AddEntry(hMuonKaon,"With Kaon track","l");
  hMuonKaon->Draw("same");
  hMuonOneKaon->SetLineColor(3);
  leg->AddEntry(hMuonOneKaon,"and One Kaon track in event","l");
  hMuonOneKaon->Draw("same");
  hMuonNoMuon->SetLineColor(4);
  leg->AddEntry(hKaonNoMuon,"and No Muon track >20cm","l");
  hMuonNoMuon->Draw("same");
  hMuonNoHits->SetLineColor(5);
  leg->AddEntry(hMuonNoMuon,"and No space points within 2cm of walls","l");
  hMuonNoHits->Draw("same");
  leg->Draw("same");
  canv->SaveAs("MuonReco.eps");
  delete leg;

  canv->cd();
  canv->Clear();
  leg = new TLegend(0.5,0.2,0.9,0.5);
  hProtonAll->SetLineColor(1);
  leg->AddEntry(hProtonAll,"No cuts","l");
  hProtonAll->Draw();
  hProtonKaon->SetLineColor(2);
  leg->AddEntry(hProtonKaon,"With Kaon track","l");
  hProtonKaon->Draw("same");
  hProtonOneKaon->SetLineColor(3);
  leg->AddEntry(hProtonOneKaon,"and One Kaon track in event","l");
  hProtonOneKaon->Draw("same");
  hProtonNoMuon->SetLineColor(4);
  leg->AddEntry(hProtonNoMuon,"and No Muon track >20cm","l");
  hProtonNoMuon->Draw("same");
  hProtonNoHits->SetLineColor(5);
  leg->AddEntry(hProtonNoMuon,"and No space points within 2cm of walls","l");
  hProtonNoHits->Draw("same");
  leg->Draw("same");
  canv->SaveAs("ProtonReco.eps");
  delete leg;

  return;

}
