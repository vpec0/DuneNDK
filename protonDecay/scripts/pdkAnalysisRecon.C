#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"
#include "TMath.h"

#include "TStyle.h"
#include "TAxis.h"
#include "TGaxis.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TH1D.h"


#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <sstream>

void pdkAnalysisRecon() {

  const int kMaxObjects = 3000;

  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetOptStat(0);
  gStyle->SetStatX(0.9);
  gStyle->SetStatY(0.937);
  gStyle->SetStatW(0.3);
  gStyle->SetStatH(0.15);
  gStyle->SetHistLineColor(1);
  TGaxis::SetMaxDigits(3);

  // CUTS --------------------------------

  float closeTrackDist = 5; //cm
  float outerSpacePointDist = 2; //cm
  float longTrackLength = 20; //cm
  float kaonEnergyCut = 250.; //MeV
  float kaonMassCut = 600.; //MeV

  // fiducial volume -- from M Robinson
  float edges[3][2];
  edges[0][0]=-745.744;
  edges[0][1]=745.744;
  edges[1][0]=-600;
  edges[1][1]=600;
  edges[2][0]=-1;
  edges[2][1]=5808.87;

  // from dune10kt_v1.gdml:
  /**
   * top: 627.58 (LAr) position w.r.t to Envelope -100
   * bottom: -677.58 (-100)
   * left/right: 745.744
   * front/back: 0/6159.35 ?
   * frame 10.16 cm thick
   * first vertical beam placed at -2772.60
   * last vertical beam at 3024.60
   * top of top APA beam at 581.80
   * bottom of bottom beam at -623.30
   **/

  TFile* inFile = new TFile("data/10k/ana_hists_merged.root","r");
  TTree* tree = (TTree*)inFile->Get("pdkrecon/PDKRecon");
  if (!tree) {
    std::cout << "Input tree not found." << std::endl;
    return;
  }

  int fNSpacePoints;
  float fSpacePoint[kMaxObjects][3];
  int fNTracks;
  int fNTrackHits[kMaxObjects];
  int fNTrackPlaneHits[kMaxObjects][3];
  int fTrackBestPlane[kMaxObjects];
  int fTrackPDG[kMaxObjects][3];
  float fTrackStart[kMaxObjects][3];
  float fTrackEnd[kMaxObjects][3];
  float fTrackLength[kMaxObjects];
  float fTrackEnergy[kMaxObjects][3];

  tree->SetBranchAddress("NSpacePoints",    &fNSpacePoints);
  tree->SetBranchAddress("SpacePoint",      fSpacePoint);
  tree->SetBranchAddress("NTracks",         &fNTracks);
  tree->SetBranchAddress("NTrackHits",      fNTrackHits);
  tree->SetBranchAddress("NTrackPlaneHits", fNTrackPlaneHits);
  tree->SetBranchAddress("TrackBestPlane",  fTrackBestPlane);
  tree->SetBranchAddress("TrackPDG",        fTrackPDG);
  tree->SetBranchAddress("TrackStart",      fTrackStart);
  tree->SetBranchAddress("TrackEnd",        fTrackEnd);
  tree->SetBranchAddress("TrackLength",     fTrackLength);
  tree->SetBranchAddress("TrackEnergy",     fTrackEnergy);

  TH1D* hCutAll                = new TH1D("CutAll",               ";Recon Deposited Energy (MeV);No. of events;",60,0,1000);
  TH1D* hCutNoLongMuon         = new TH1D("CutNoLongMuon",        ";Recon Deposited Energy (MeV);No. of events;",60,0,1000);
  TH1D* hCutOneKaon            = new TH1D("CutOneKaon",           ";Recon Deposited Energy (MeV);No. of events;",60,0,1000);
  TH1D* hCutNoOuterSpacePoints = new TH1D("CutNoOuterSpacePoints",";Recon Deposited Energy (MeV);No. of events;",60,0,1000);
  TH1D* hCutNoBothActivity     = new TH1D("CutNoBothActivity",    ";Recon Deposited Energy (MeV);No. of events;",60,0,1000);
  TH1D* hCutKaonEnergy         = new TH1D("CutKaonTrackEnergy",   ";Recon Deposited Energy (MeV);No. of events;",60,0,1000);
  TH1D* hCutKaonMass           = new TH1D("CutKaonMass",          ";Recon Deposited Energy (MeV);No. of events;",60,0,1000);

  for (int event = 0; event < tree->GetEntriesFast(); ++event) {

    tree->GetEntry(event);

    if (event%10000 == 0)
      std::cout << "Processing event " << event << std::endl;

    float totalDepositedEnergy = 0;
    std::vector<int> longMuonTracks, kaonTracks;
    std::map<int,std::vector<int> > closeFrontTracks, closeBackTracks;
    std::vector<int> outerSpacePoints;

    // tracks
    for (int track = 0; track < fNTracks; ++track) {
      TVector3 start = TVector3(fTrackStart[track]), end = TVector3(fTrackEnd[track]);
      totalDepositedEnergy += fTrackEnergy[track][fTrackBestPlane[track]];
      if (fTrackPDG[track][fTrackBestPlane[track]] == 321)
	kaonTracks.push_back(track);
      if (TMath::Abs(fTrackPDG[track][fTrackBestPlane[track]]) == 13 && fTrackLength[track] > longTrackLength)
	longMuonTracks.push_back(track);
      for (int otherTrack = 0; otherTrack < fNTracks; ++otherTrack) {
	if (track == otherTrack)
	  continue;
	TVector3 otherStart = TVector3(fTrackStart[otherTrack]), otherEnd = TVector3(fTrackStart[otherTrack]);
	if ((start-otherStart).Mag() < closeTrackDist || (start-otherEnd).Mag() < closeTrackDist)
	  closeFrontTracks[track].push_back(otherTrack);
	if ((end-otherStart).Mag() < closeTrackDist || (end-otherEnd).Mag() < closeTrackDist)
	  closeBackTracks[track].push_back(otherTrack);
      }
    }

    // space points
    for (int spacePoint = 0; spacePoint < fNSpacePoints; ++spacePoint) {
      const float* xyz = fSpacePoint[spacePoint];
      if ((xyz[0] >= edges[0][0] && xyz[0] < edges[0][0] + outerSpacePointDist) || (xyz[0] <= edges[0][1] && xyz[0] > edges[0][1] - outerSpacePointDist) ||
	  (xyz[1] >= edges[1][0] && xyz[1] < edges[1][0] + outerSpacePointDist) || (xyz[1] <= edges[1][1] && xyz[1] > edges[1][1] - outerSpacePointDist) ||
	  (xyz[2] >= edges[2][0] && xyz[2] < edges[2][0] + outerSpacePointDist) || (xyz[2] <= edges[2][1] && xyz[2] > edges[2][1] - outerSpacePointDist))
	outerSpacePoints.push_back(spacePoint);
    }

    //all
    hCutAll->Fill(totalDepositedEnergy);

    if (longMuonTracks.size() == 0) {
      //no long muons
      hCutNoLongMuon->Fill(totalDepositedEnergy);

      if (kaonTracks.size() == 1) {
	//one kaon track
	hCutOneKaon->Fill(totalDepositedEnergy);
	int kaon = kaonTracks[0];

	if (outerSpacePoints.size() == 0) {
	  // no hits near detector edges
	  hCutNoOuterSpacePoints->Fill(totalDepositedEnergy);

	  if (!(closeFrontTracks[kaon].size() != 0 && closeBackTracks[kaon].size() != 0)) {
	    // no activity at both ends of the kaon
	    hCutNoBothActivity->Fill(totalDepositedEnergy);
	    std::vector<int> closeTracks;
	    if (closeFrontTracks[kaon].size() || closeBackTracks[kaon].size())
	      closeTracks = closeFrontTracks[kaon].size() ? closeFrontTracks[kaon] : closeBackTracks[kaon];

	    if (fTrackEnergy[kaon][fTrackBestPlane[kaon]] <= kaonEnergyCut) {
	      // not too energetic
	      hCutKaonEnergy->Fill(totalDepositedEnergy);

	      float depositedEnergy = 0;
	      for (std::vector<int>::const_iterator closeTrackIt = closeTracks.begin(); closeTrackIt != closeTracks.end(); ++closeTrackIt) {
		depositedEnergy += fTrackEnergy[*closeTrackIt][fTrackBestPlane[*closeTrackIt]];
		for (std::vector<int>::const_iterator ccTrackIt = closeFrontTracks[*closeTrackIt].begin(); ccTrackIt != closeFrontTracks[*closeTrackIt].end(); ++ccTrackIt)
		  if (*ccTrackIt != kaon)
		    depositedEnergy += fTrackEnergy[*ccTrackIt][fTrackBestPlane[*ccTrackIt]];
		for (std::vector<int>::const_iterator ccTrackIt = closeBackTracks[*closeTrackIt].begin(); ccTrackIt != closeBackTracks[*closeTrackIt].end(); ++ccTrackIt)
		  if (*ccTrackIt != kaon)
		    depositedEnergy += fTrackEnergy[*ccTrackIt][fTrackBestPlane[*ccTrackIt]];
	      }

	      if (depositedEnergy <= kaonMassCut) {
		// not too massive
		hCutKaonMass->Fill(totalDepositedEnergy);

		// PASS ALL CUTS!
		std::cout << "Event " << event << " passes all cuts" << std::endl;
		std::cout << "  Kaon energy " << fTrackEnergy[kaon][fTrackBestPlane[kaon]] << "; following deposited energy " << depositedEnergy << std::endl;

	      }
	    }
	  }
	}
      }
    }

  } // event

  TCanvas* canv = new TCanvas("canv","",800,600);
  canv->SetLogy();
  hCutAll->SetLineWidth(2);
  hCutAll->SetLineColor(1);
  //hCutAll->SetFillColor(1);
  //hCutAll->SetFillStyle(3003);
  hCutAll->Draw();
  std::stringstream cutAll;
  cutAll << "All events (" << hCutAll->GetEntries() << ")";
  hCutNoLongMuon->SetLineWidth(2);
  hCutNoLongMuon->SetLineColor(2);
  //hCutNoLongMuon->SetFillColor(2);
  //hCutNoLongMuon->SetFillStyle(3003);
  hCutNoLongMuon->Draw("same");
  std::stringstream cutNoLongMuon;
  cutNoLongMuon << "No long muon track (>" << longTrackLength << " cm) (" << hCutNoLongMuon->GetEntries() << ")";
  hCutOneKaon->SetLineWidth(2);
  hCutOneKaon->SetLineColor(3);
  //hCutOneKaon->SetFillColor(3);
  //hCutOneKaon->SetFillStyle(3003);
  hCutOneKaon->Draw("same");
  std::stringstream cutOneKaon;
  cutOneKaon << "One kaon (" << hCutOneKaon->GetEntries() << ")";
  hCutNoOuterSpacePoints->SetLineWidth(2);
  hCutNoOuterSpacePoints->SetLineColor(4);
  //hCutNoOuterSpacePoints->SetFillColor(4);
  //hCutNoOuterSpacePoints->SetFillStyle(3003);
  hCutNoOuterSpacePoints->Draw("same");
  std::stringstream cutNoOuterSpacePoints;
  cutNoOuterSpacePoints << "No hits <" << outerSpacePointDist << " cm from walls (" << hCutNoOuterSpacePoints->GetEntries() << ")";
  hCutNoBothActivity->SetLineWidth(2);
  hCutNoBothActivity->SetLineColor(5);
  //hCutNoBothActivity->SetFillColor(5);
  //hCutNoBothActivity->SetFillStyle(3003);
  hCutNoBothActivity->Draw("same");
  std::stringstream cutNoBothActivity;
  cutNoBothActivity << "Kaons without activity at both ends (" << hCutNoBothActivity->GetEntries() << ")";
  hCutKaonEnergy->SetLineWidth(2);
  hCutKaonEnergy->SetLineColor(6);
  //hCutKaonEnergy->SetFillColor(6);
  //hCutKaonEnergy->SetFillStyle(3003);
  hCutKaonEnergy->Draw("same");
  std::stringstream cutKaonEnergy;
  cutKaonEnergy << "Kaon energy <" << kaonEnergyCut << " MeV (" << hCutKaonEnergy->GetEntries() << ")";
  hCutKaonMass->SetLineWidth(2);
  hCutKaonMass->SetLineColor(7);
  //hCutKaonMass->SetFillColor(7);
  //hCutKaonMass->SetFillStyle(3003);
  hCutKaonMass->Draw("same");
  std::stringstream cutKaonMass;
  cutKaonMass << "Kaon mass <" << kaonMassCut << " MeV (" << hCutKaonMass->GetEntries() << ")";
  TLegend* leg = new TLegend(0.45,0.45,0.87,0.87);
  leg->AddEntry(hCutAll,                cutAll.str().c_str(),                "l");
  leg->AddEntry(hCutNoLongMuon,         cutNoLongMuon.str().c_str(),         "l");
  leg->AddEntry(hCutOneKaon,            cutOneKaon.str().c_str(),            "l");
  leg->AddEntry(hCutNoOuterSpacePoints, cutNoOuterSpacePoints.str().c_str(), "l");
  leg->AddEntry(hCutNoBothActivity,     cutNoBothActivity.str().c_str(),     "l");
  leg->AddEntry(hCutKaonEnergy,         cutKaonEnergy.str().c_str(),         "l");
  leg->AddEntry(hCutKaonMass,           cutKaonMass.str().c_str(),           "l");
  leg->Draw();
  canv->SaveAs("plots/10k_1/pdk_selection/PDKReconSelection.pdf");
  canv->SetLogy(0);

  return;

}
