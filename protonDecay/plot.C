{

  //8cm gap at top and bottom
  //0.7 cm left and right
  //0cm front and back

  const int MAXKAON=1000;
  const int MAXHISTORY=100;

  TFile* inFile = new TFile("pdk_truth.root","READ");
  TTree* tree = (TTree*)inFile->Get("pdktruth/kaons");

  int nkaon;
  unsigned runid;
  unsigned eventid;
  double ekaon     [MAXKAON];
  double edaughters[MAXKAON];
  double edecay    [MAXKAON];
  double enear     [MAXKAON];
  int pdgcode      [MAXKAON];
  double muonrange;
  double etotal;
  double extent[3][2];

  tree->SetBranchAddress("run",        &runid);
  tree->SetBranchAddress("event",      &eventid);
  tree->SetBranchAddress("muonrange",  &muonrange);
  tree->SetBranchAddress("etotal",     &etotal);
  tree->SetBranchAddress("enear",      &enear);
  tree->SetBranchAddress("extent",     extent);
  tree->SetBranchAddress("nkaon",      &nkaon);
  tree->SetBranchAddress("ekaon",      ekaon);
  tree->SetBranchAddress("edaughters", edaughters);
  tree->SetBranchAddress("edecay",     edecay);
  tree->SetBranchAddress("pdgcode",    pdgcode);

  const double cmin = 0.1;
  const double cmax = 300;
  const int cnbins = 40;
  double cbins[cnbins];
  for (int i = 0; i <= cnbins; ++i)
    cbins[i] = pow(10, log10(cmin) + i*(log10(cmax)-log10(cmin))/cnbins);
  const double bmin = 1;
  const double bmax = 100000;
  const int nbins = 40;
  double bins[nbins];
  for (int i = 0; i <= nbins; ++i)
    bins[i] = pow(10, log10(bmin) + i*(log10(bmax)-log10(bmin))/nbins);

  TH2D* h2       = new TH2D("h2", ";Energy deposition from K^{+/-} decay products [MeV];K^{+/-} Energy Deposition [MeV]", 25, 0., 250., 25, 0., 200.);
  //TH1D* hall     = new TH1D("hall", "Charged Kaon Spectrum, 4x10^{8} filtered muons;MeV;/MeV", nbins, bins);
  TH1D* hall     = new TH1D("hall", ";MeV;/MeV", nbins, bins);
  TH1D* hseen    = new TH1D("hseen", "Charged Kaon Spectrum;MeV;/MeV", nbins, bins);
  TH1D* hone     = new TH1D("hone", "Charged Kaon Spectrum;MeV;/MeV", nbins, bins);
  TH1D* hnomu    = new TH1D("hnomu", "Charged Kaon Spectrum;MeV;/MeV", nbins, bins);
  TH1D* hedge    = new TH1D("hedge", "Charged Kaon Spectrum;MeV;/MeV", nbins, bins);
  TH1D* hplus    = new TH1D("hplus", "K^{+} Spectrum;MeV;/MeV", nbins, bins);
  TH1D* hminus   = new TH1D("hminus", "K^{-} Spectrum;MeV;/MeV", nbins, bins);
  //TH1D* hclosest = new TH1D("hclosest", "Closest hit to wall;cm;/cm", cnbins, cbins);
  TH1D* hclosest = new TH1D("hclosest", ";cm;/cm", cnbins, cbins);
  //TH1D* hxmin    = new TH1D("hxmin", "Closest hit to wall (active edge of cryostat);cm;/cm", cnbins, cbins);
  TH1D* hxmin    = new TH1D("hxmin", ";cm;/cm", cnbins, cbins);
  TH1D* hxmax    = new TH1D("hxmax" , ";cm;/cm", cnbins, cbins);
  TH1D* hymin    = new TH1D("hymin", ";cm;/cm", cnbins, cbins);
  TH1D* hymax    = new TH1D("hymax", ";cm;/cm", cnbins, cbins);
  TH1D* hzmin    = new TH1D("hzmin", ";cm;/cm", cnbins, cbins);
  TH1D* hzmax    = new TH1D("hzmax", ";cm;/cm", cnbins, cbins);
  //TH1D* hold     = new TH1D("hold", "Charged Kaon Spectrum;MeV;/MeV", nbins, bins);

  const int N = 10000;
  double xall[N];
  double yall[N];
  int nall = 0;
  double xmupass[N];
  double ymupass[N];
  int nmupass = 0;

  int ntotal = 0;
  int nseen = 0;
  int none = 0;
  int nnomu = 0;
  int nplus = 0;
  int nminus = 0;
  int nedge = 0;
  TDatabasePDG* db = new TDatabasePDG();
  db->ReadPDGTable();

  // fiducial volume
  double edges[3][2], edgedist[3][2];
  edges[0][0]=-745.744;
  edges[0][1]=745.744;
  edges[1][0]=-600;
  edges[1][1]=600;
  edges[2][0]=-1;
  edges[2][1]=5808.87;

  // look at events
  for (int event = 0; event < tree->GetEntriesFast(); ++event) {

    tree->GetEntry(event);
    std::cout << std::endl << "Event " << event << std::endl;

    // determine distance to detector boundary
    for (int coord = 0; coord < 3; ++coord) {
      std::cout << extent[coord][0] << ", " << extent[coord][1] << std::endl;
      edgedist[coord][0] = extent[coord][0] - edges[coord][0];
      edgedist[coord][1] = edges[coord][1] - extent[coord][1];
    }

    // get closest ide to detector boundary
    double closest=9999999.;
    for (int coord = 0; coord < 3; ++coord) {
      for (int hilo = 0; hilo < 2; ++hilo) {
	edgedist[coord][hilo] = edgedist[coord][hilo] < 0 ? 0 : edgedist[coord][hilo];
	if (edgedist[coord][hilo] < closest)
	  closest = edgedist[coord][hilo];
      }
    }

    // fill boundary edge histograms
    hxmin->Fill(edgedist[0][0]);
    hxmax->Fill(edgedist[0][1]);
    hymin->Fill(edgedist[1][0]);
    hymax->Fill(edgedist[1][1]);
    hzmin->Fill(edgedist[2][0]);
    hzmax->Fill(edgedist[2][1]);
    hclosest->Fill(closest);

    std::cout << nkaon << " charged kaons" << std::endl;
    std::cout << "total energy " << etotal << " MeV" << std::endl;
    std::cout << "muon range " << muonrange << " cm" << std::endl;
    std::cout << "distance to wall " << closest << std::endl;

    ntotal += nkaon;

    // count how many kaons how non-zero energy
    int nk = 0;
    for (int kaonIt = 0; kaonIt < nkaon; ++kaonIt) {
      if (ekaon[kaonIt])
	++nk;
    }

    // look at all kaons
    for (int kaonIt = 0; kaonIt < nkaon; ++kaonIt) {

      // energy
      double e = ekaon[kaonIt] + edaughters[kaonIt];
      hall->Fill(e);
      if (pdgcode[kaonIt] == 321) {
	hplus->Fill(e);
	++nplus;
      }
      else if (pdgcode[kaonIt] == -321) {
	hminus->Fill(e);
	++nminus;
      }

      if (ekaon[kaonIt]) {

	std::cout << "kaon " << kaonIt << ":" << std::endl;
	std::cout << "\t ekaon = " << ekaon[kaonIt] << std::endl;
	std::cout << "\t edaughters = " << edaughters[kaonIt] << std::endl;
	std::cout << "\t edecay = " << edecay[kaonIt] << std::endl;
	std::cout << "\t enear = " << enear[kaonIt] << std::endl;

	hseen->Fill(e);
	++nseen;

	// single kaon
	if (nk == 1) {

	  ++none;
	  hone->Fill(e);

	  // no muon
	  if (muonrange < 20) {

	    hnomu->Fill(e);
	    ++nnomu;

	    double ey = e;//ekaon[kaonIt]+edaughters[kaonIt];//+edecay[fk];
	    h2->Fill(edecay[kaonIt], ey);

	    double ex = etotal - ey - edecay[kaonIt];

	    // no hits within 2cm of walls
	    if (closest > 2) {

	      hedge->Fill(e);
	      ++nedge;

	      std::cout << "ALL PASSED (" << closest << " cm): " << runid << ", " << eventid;
	      //std::cout << "edge dists: " << extent[2][fk] << ", " << extent[2][1] << std::endl;
	      yall[nall] = ey;
	      xall[nall] = ex;
	      if (ex < 30 && ey < 300)
		std::cout << " ROI, Entry " << event;

	      ++nall;

	    }

	    // hits within 2cm of walls -- muon passing!
	    else {

	      ymupass[nmupass] = ey;
	      xmupass[nmupass] = ex;

	      if (ekaon[kaonIt] + edaughters[kaonIt] < 10 and edecay[kaonIt] < 10)
		std::cout << "LOW E: "  << runid << ", " << eventid;
	      else {
		std::cout << "MOST PASSED: " << runid << ", " << eventid;
		if (ex < 30 && ey < 300)
		  std::cout << " ROI, Entry " << event;
	      }

	      ++nmupass;

	    }

	    TParticlePDG* particle = db->GetParticle(pdgcode[kaonIt]);
	    std::cout << "\t" << (particle ? particle->GetName() : "UNKNOWN")
		      << " ekaon = " << ekaon[kaonIt]
		      << ", secondaries = " << edaughters[kaonIt]
		      << ", edecay = " << edecay[kaonIt]
		      << ", enear = " << enear[kaonIt]
		      << ", eother = " << (etotal-ekaon[kaonIt]-edaughters[kaonIt]-edecay[kaonIt]-enear[kaonIt]) << ":";

	  } // muon range

	} // single kaon

      } // kaon track

    } // loop over kaons

    std::cout << std::endl;

  } // event loop

  std::cout << ntotal << " charged Kaon" << std::endl;
  //  hold->Fill(55.3036+96.8055);
  //  hold->Fill(21.4718+126.812);
  //  hold->Fill(131.078+155.738);
  //  hold->Fill(97.5077+672.39);
  //  hold->Fill(103.878+670.39);
  //  hold->Fill(651.744+469.394);
  //  hold->Fill(322.957+198.038);
  //  hold->Fill(562.8484+1387.62);

  for (int bin = 0; bin < nbins; ++bin) {
    //hold  ->SetBinContent(bin+1, hold  ->GetBinContent(bin+1)/(bins[bin+1]-bins[bin]));
    hall  ->SetBinContent(bin+1, hall  ->GetBinContent(bin+1)/(bins[bin+1]-bins[bin]));
    hnomu ->SetBinContent(bin+1, hnomu ->GetBinContent(bin+1)/(bins[bin+1]-bins[bin]));
    hone  ->SetBinContent(bin+1, hone  ->GetBinContent(bin+1)/(bins[bin+1]-bins[bin]));
    hseen ->SetBinContent(bin+1, hseen ->GetBinContent(bin+1)/(bins[bin+1]-bins[bin]));
    hplus ->SetBinContent(bin+1, hplus ->GetBinContent(bin+1)/(bins[bin+1]-bins[bin]));
    hminus->SetBinContent(bin+1, hminus->GetBinContent(bin+1)/(bins[bin+1]-bins[bin]));
    hedge ->SetBinContent(bin+1, hedge ->GetBinContent(bin+1)/(bins[bin+1]-bins[bin]));
  }

  for (int bin = 0; bin < cnbins; ++bin) {
    hclosest->SetBinContent(bin+1, hclosest->GetBinContent(bin+1)/(cbins[bin+1]-cbins[bin]));
    hxmin   ->SetBinContent(bin+1, hxmin   ->GetBinContent(bin+1)/(cbins[bin+1]-cbins[bin]));
    hxmax   ->SetBinContent(bin+1, hxmax   ->GetBinContent(bin+1)/(cbins[bin+1]-cbins[bin]));
    hymin   ->SetBinContent(bin+1, hymin   ->GetBinContent(bin+1)/(cbins[bin+1]-cbins[bin]));
    hymax   ->SetBinContent(bin+1, hymax   ->GetBinContent(bin+1)/(cbins[bin+1]-cbins[bin]));
    hzmin   ->SetBinContent(bin+1, hzmin   ->GetBinContent(bin+1)/(cbins[bin+1]-cbins[bin]));
    hzmax   ->SetBinContent(bin+1, hzmax   ->GetBinContent(bin+1)/(cbins[bin+1]-cbins[bin]));
  }
  //hall->SetEntries(ntotal);

  // plot
  TCanvas* can = new TCanvas;
  can->SetLogy(kTRUE);
  can->SetLogx(kTRUE);
  can->cd();
  hall->SetStats(kFALSE);
  hall->GetXaxis()->SetNoExponent(kTRUE);
  hall->SetLineColor(kBlack);
  hall->SetLineWidth(2);
  hall->Draw();
  hseen->SetLineColor(kBlue);
  hseen->SetLineWidth(2);
  hseen->Draw("SAME");
  hone->SetLineColor(kRed);
  hone->SetLineWidth(2);
  hone->Draw("SAME");
  hnomu->SetLineColor(kBlack);
  hnomu->SetLineStyle(2);
  hnomu->SetLineWidth(2);
  hnomu->Draw("SAME");
  hedge->SetLineColor(kBlue);
  hedge->SetLineStyle(2);
  hedge->SetLineWidth(2);
  hedge->Draw("SAME");
  ostringstream tall;
  ostringstream tseen;
  ostringstream tone;
  ostringstream tnomu;
  ostringstream tedge;
  tall  << "Total depositing energy ("          << ntotal << ")"; 
  tseen << "With Kaon track ("                  << nseen  << ")";
  tone  << "and One Kaon track in event ("      << none   << ")";
  tnomu << "and No Muon track >20 cm ("         << nnomu  << ")";
  tedge << "and No hits within 2 cm of walls (" << nedge  << ")";
  TLegend* leg = new TLegend(0.50, 0.90, 0.95, 0.65);
  leg->AddEntry(hall,  tall.str().c_str(),  "l");
  leg->AddEntry(hseen, tseen.str().c_str(), "l");
  leg->AddEntry(hone,  tone.str().c_str(),  "l");
  leg->AddEntry(hnomu, tnomu.str().c_str(), "l");
  leg->AddEntry(hedge, tedge.str().c_str(), "l");
  leg->Draw();
  can->SaveAs("kaon_spectrum.eps");
  can->SaveAs("kaon_spectrum.png");
  can->SaveAs("kaon_spectrum.pdf");

  // plot
  can = new TCanvas;
  can->SetLogy(kTRUE);
  can->SetLogx(kTRUE);
  can->cd();
  hall->SetLineColor(kBlack);
  hall->Draw();
  hplus->SetLineColor(kRed);
  hplus->Draw("SAME");
  hminus->SetLineColor(kBlue);
  hminus->Draw("SAME");
  tall.str("");
  ostringstream tplus;
  ostringstream tminus;
  tall   << "charged K tree (" << ntotal << ")";
  tplus  << "K+ tree only ("   << nplus  << ")";
  tminus << "K- tree only ("   << nminus << ")";
  leg = new TLegend(0.65, 0.90, 0.95, 0.75);
  leg->AddEntry(hall, tall.str().c_str(), "l");
  leg->AddEntry(hplus, tplus.str().c_str(), "l");
  leg->AddEntry(hminus, tminus.str().c_str(), "l");
  leg->Draw();
  //  can->SaveAs("kaon_spectrum.pdf");

  // plot
  can = new TCanvas;
  can->SetLogx(kTRUE);
  can->SetLogy(kTRUE);
  can->cd();
  tmp = can->DrawFrame(0.1,0.03,300,10000);
  tmp->GetXaxis()->SetNoExponent(kTRUE);
  //tmp->SetTitle("Closest hit to wall (edge of cryostat active);cm;/cm");
  tmp->SetTitle(";cm;/cm");
  can->Modified();
  hclosest->SetLineColor(kBlack);
  hclosest->Draw("SAME");
  hxmin->SetLineColor(kRed);
  hxmin->SetLineStyle(2);
  hxmin->Draw("SAME");
  hxmax->SetLineColor(kBlue);
  hxmax->SetLineStyle(2);
  hxmax->Draw("SAME");
  hymin->SetLineColor(kRed);
  hymin->SetLineStyle(3);
  hymin->Draw("SAME");
  hymax->SetLineColor(kBlue);
  hymax->SetLineStyle(3);
  hymax->Draw("SAME");
  hzmin->SetLineColor(kRed);
  hzmin->SetLineStyle(4);
  hzmin->Draw("SAME");
  hzmax->SetLineColor(kBlue);
  hzmax->SetLineStyle(4);
  hzmax->Draw("SAME");
  leg = new TLegend(0.65, 0.90, 0.95, 0.60);
  leg->AddEntry(hclosest, "Closest Overall", "l");
  leg->AddEntry(hxmin, "To Low X Limit", "l");
  leg->AddEntry(hxmax, "To High X Limit", "l");
  leg->AddEntry(hymin, "To Low Y Limit", "l");
  leg->AddEntry(hymax, "To High Y Limit", "l");
  leg->AddEntry(hzmin, "To Low Z Limit", "l");
  leg->AddEntry(hzmax, "To High Z Limit", "l");
  leg->Draw();
  //  can->SaveAs("kaon_spectrum.pdf");

  // plot
  can = new TCanvas;
  can->cd();
  tmp = can->DrawFrame(1, 1, 100000, 10000);
  tmp->SetXTitle("Other energy deposition [MeV]");
  tmp->SetYTitle("K^{+/-} energy deposition [MeV]");
  //tmp->SetTitle("Background Candidates");
  can->Modified();
  TText *t = new TText(6.1,15.8,"Region of interest");
  t->SetTextAlign(22);
  t->SetTextColor(kBlack);
  t->SetTextFont(43);
  t->SetTextSize(26);
  t->SetTextAngle(60);
  t->Draw();
  can->Modified();
  can->SetLogx(kTRUE);
  can->SetLogy(kTRUE);
  TGraph* gmupass = new TGraph(nmupass, xmupass, ymupass);
  gmupass->SetMarkerColor(kBlue);
  gmupass->SetMarkerSize(1);
  gmupass->SetMarkerStyle(21);
  gmupass->Draw("p");
  TGraph* gall = new TGraph(nall, xall, yall);
  gall->SetMarkerColor(kRed);
  gall->SetMarkerSize(1);
  gall->SetMarkerStyle(20);
  gall->Draw("p");
  int nl = 3;
  double xl[3];
  double yl[3];
  xl[0]=0.01;
  yl[0]=300;
  xl[1]=30;
  yl[1]=300;
  xl[2]=30;
  yl[2]=0.01;
  TGraph* line=new TGraph(nl, xl, yl);
  line->SetLineColor(kBlack);
  line->SetLineStyle(2);
  //line->GetXaxis()->SetLimits(40,100000);  
  //line->GetYaxis()->SetLimits(0,1000);
  line->Draw("l");
  leg = new TLegend(0.15, 0.90, 0.65, 0.80);
  ostringstream tmupass;
  tmupass << "Passing all but fiducial and energy cuts (" << nmupass << ")";
  leg->AddEntry(gmupass, tmupass.str().c_str(), "p");
  ostringstream tallpass;
  tallpass << "Also passing fiducial cut (" << nall << ")";
  leg->AddEntry(gall, tallpass.str().c_str(), "p");
  leg->Draw();
  can->SaveAs("kaon_scatter.pdf");
  can->SaveAs("kaon_scatter.eps");
  can->SaveAs("kaon_scatter.png");

  can = new TCanvas;
  gStyle->SetStatX(0.9);
  //can->SetLogz(kTRUE);
  //h2->SetStats(kFALSE);
  h2->SetStats(kFALSE);
  h2->Draw("colz");
  can->SaveAs("kaon_2d.pdf");
  can->SaveAs("kaon_2d.png");
  can->SaveAs("kaon_2d.eps");

}
