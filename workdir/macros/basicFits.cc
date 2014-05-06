void fitMll()
{
  RooRealVar m_Z0( "m_Z0",     "Z0 mass",  
		   91.188, "GeV/c^{2}" );
  RooRealVar gamma_Z0( "gamma_Z0", "Z0 width",  
		       2.4952, "GeV/c^{2}" );

  RooRealVar m_ee( "m_ee","m_ee",50,120,"GeV/c^{2}" );
  RooRealVar cb_bias_ee( "cb_bias_ee", "bias", 
			 0.0, -1.0, 1.0,"GeV/c^{2}" ); 
  RooRealVar cb_width_ee( "cb_width_ee","width", 
			  1.31,0.9,1.8,"GeV/c^{2}" ); 
  RooRealVar cb_alpha_ee( "cb_alpha_ee","alpha", 1.1,0.6,2.0 ); 
  RooRealVar cb_power_ee( "cb_power_ee","power", 3.5, 0.5, 5.0 ); 
  RooCBShape cb_pdf_ee( "cb_pdf_ee", "A  Crystal Ball Lineshape", 
			m_ee, 
			cb_bias_ee,  cb_width_ee, 
			cb_alpha_ee, cb_power_ee );
  RooBreitWigner bw_Z0( "bw_Z0","The true Z0 lineshape (BW)", 
			m_ee , m_Z0, gamma_Z0 );
  RooNumConvPdf bw_ee("bw_ee","Convolution", m_ee, bw_Z0, cb_pdf_ee );
  bw_ee.setConvolutionWindow( cb_bias_ee, cb_width_ee, 50 );
  RooRealVar  lambda_ee( "lambda_ee",    "exponent of background", 
			 0., -100., 0.);
  RooExponential exp_ee( "pdf_bkg_ee", 
			 "Exponential background", 
			 m_ee, lambda_ee );
  RooRealVar n_Z_ee( "n_Z_ee","n_Z_ee",     0., 20000.);
  RooRealVar n_bkg_ee( "n_bkg_ee","n_bkg_ee", 0., 1000.);
  RooRealVar f_bkg_ee( "f_bkg_ee","f_bkg_ee", 0., 1.);

  //RooRealVar xp("xp","xp",91,80,100,"GeV/c^{2}");
  //RooRealVar sigp("sigp","sigp",2.5,0,10,"GeV/c^{2}");
  //RooRealVar xi("xi","xi",-0.2002,-5,5);
  //RooRealVar rho1("rho1","rho1",0.1165,0,5);
  //RooRealVar rho2("rho2","rho2",-0.168,-5,5);
  //RooBukinPdf bukin_ee("bukin", "bukin", m_ee, xp, sigp, xi, rho1, rho2 );
 
  //  RooRealVar mean("xp","xp",91,80,100,"GeV/c^{2}");
  //  RooRealVar sigL("sigL","sigL",2.5,0,10,"GeV/c^{2}");
  //  RooRealVar sigR("sigR","sigR",2.5,0,10,"GeV/c^{2}");
  //  RooBifurGauss bg_ee("bg_ee", "bg_ee", m_ee, mean, sigL, sigR );
  
  //  RooArgList listPdf(       exp_ee, bg_ee   );
  RooArgList listPdf(       exp_ee, bw_ee   );
  //RooArgList listPdf(       exp_ee, bukin_ee   );
  //  RooArgList listPdfCoef( n_Z_ee,  n_bkg_ee );
  RooArgList listPdfCoef( f_bkg_ee );
  RooAddPdf pdf( "pdf", "PDF ee", listPdf, listPdfCoef );
  RooArgSet* params = pdf.getParameters( RooDataSet() );
  params->readFromFile( "config/Zee.config", "READ", "Zee: Fit Parameters" );
  cout << "Parameters configured" << endl;
  params->selectByAttrib("READ",kTRUE)->Print("v") ;
  cout << "Parameters NOT configured" << endl;
  params->selectByAttrib("READ",kFALSE)->Print("v") ;

  //
  // fitting 
  //
  //  RooFitResult* result =  pdf.fitTo(*dataset,RooFit::FitOptions("rm0") );
  RooFitResult* result =  RooFitUtils::fit( pdf, *dataset );

  //
  // RooPlot
  //
  RooPlot* plot = m_ee.frame(140);
  dataset->plotOn(plot);
  pdf.plotOn(plot, RooFit::Components("pdf_bkg_ee"), RooFit::LineStyle(kDashed), RooFit::LineColor(kRed));
  pdf.plotOn(plot);
  plot->Draw();

  //
  // Gautier's version of RooPlot... :)
  //
  TCanvas* c = RooFitUtils::newAddPdfPlot( pdf, *dataset, m_ee );
  c->SetLogy();
  c->Print();
}

void fitMet()
{
  RooRealVar peak("peak","peak",10,0,50);
  RooRealVar tail("tail","tail",0.1,-1,5);
  RooRealVar width("width","width",10,0,50);
  RooRealVar met("met","met",0,100,"GeV/c^{2}");
  RooNovosibirsk pdf("novosibirsk", "novosibirsk", met, peak, width, tail );
  RooArgSet* params = pdf.getParameters( RooDataSet() );
  params->readFromFile( "config/Zee.config", "READ", "Zee: Fit Parameters" );
  cout << "Parameters configured" << endl;
  params->selectByAttrib("READ",kTRUE)->Print("v") ;
  cout << "Parameters NOT configured" << endl;
  params->selectByAttrib("READ",kFALSE)->Print("v") ;

  RooFitResult* result =  pdf.fitTo(*dataset,"rm0","c");
  result->Print();
  RooPlot* plot = met.frame(100);
  dataset->plotOn(plot);
  pdf.plotOn(plot);
  plot->Draw();
}

void
fitMatthieu()
{
   TH1 *MET = new TH1D("MET","MET with PU",150,0,150);
   MET->SetBinContent(1,36);
   MET->SetBinContent(2,53);
   MET->SetBinContent(3,100);
   MET->SetBinContent(4,140);
   MET->SetBinContent(5,188);
   MET->SetBinContent(6,243);
   MET->SetBinContent(7,275);
   MET->SetBinContent(8,350);
   MET->SetBinContent(9,339);
   MET->SetBinContent(10,411);
   MET->SetBinContent(11,435);
   MET->SetBinContent(12,528);
   MET->SetBinContent(13,537);
   MET->SetBinContent(14,585);
   MET->SetBinContent(15,567);
   MET->SetBinContent(16,603);
   MET->SetBinContent(17,624);
   MET->SetBinContent(18,679);
   MET->SetBinContent(19,677);
   MET->SetBinContent(20,707);
   MET->SetBinContent(21,719);
   MET->SetBinContent(22,749);
   MET->SetBinContent(23,733);
   MET->SetBinContent(24,776);
   MET->SetBinContent(25,754);
   MET->SetBinContent(26,850);
   MET->SetBinContent(27,812);
   MET->SetBinContent(28,848);
   MET->SetBinContent(29,823);
   MET->SetBinContent(30,839);
   MET->SetBinContent(31,863);
   MET->SetBinContent(32,780);
   MET->SetBinContent(33,784);
   MET->SetBinContent(34,776);
   MET->SetBinContent(35,869);
   MET->SetBinContent(36,799);
   MET->SetBinContent(37,765);
   MET->SetBinContent(38,760);
   MET->SetBinContent(39,764);
   MET->SetBinContent(40,758);
   MET->SetBinContent(41,697);
   MET->SetBinContent(42,735);
   MET->SetBinContent(43,640);
   MET->SetBinContent(44,694);
   MET->SetBinContent(45,683);
   MET->SetBinContent(46,653);
   MET->SetBinContent(47,619);
   MET->SetBinContent(48,628);
   MET->SetBinContent(49,589);
   MET->SetBinContent(50,564);
   MET->SetBinContent(51,574);
   MET->SetBinContent(52,560);
   MET->SetBinContent(53,490);
   MET->SetBinContent(54,519);
   MET->SetBinContent(55,476);
   MET->SetBinContent(56,407);
   MET->SetBinContent(57,410);
   MET->SetBinContent(58,414);
   MET->SetBinContent(59,346);
   MET->SetBinContent(60,367);
   MET->SetBinContent(61,349);
   MET->SetBinContent(62,339);
   MET->SetBinContent(63,317);
   MET->SetBinContent(64,279);
   MET->SetBinContent(65,270);
   MET->SetBinContent(66,258);
   MET->SetBinContent(67,233);
   MET->SetBinContent(68,222);
   MET->SetBinContent(69,205);
   MET->SetBinContent(70,202);
   MET->SetBinContent(71,183);
   MET->SetBinContent(72,192);
   MET->SetBinContent(73,160);
   MET->SetBinContent(74,168);
   MET->SetBinContent(75,137);
   MET->SetBinContent(76,124);
   MET->SetBinContent(77,133);
   MET->SetBinContent(78,87);
   MET->SetBinContent(79,99);
   MET->SetBinContent(80,80);
   MET->SetBinContent(81,91);
   MET->SetBinContent(82,83);
   MET->SetBinContent(83,92);
   MET->SetBinContent(84,61);
   MET->SetBinContent(85,70);
   MET->SetBinContent(86,59);
   MET->SetBinContent(87,53);
   MET->SetBinContent(88,47);
   MET->SetBinContent(89,43);
   MET->SetBinContent(90,48);
   MET->SetBinContent(91,35);
   MET->SetBinContent(92,28);
   MET->SetBinContent(93,32);
   MET->SetBinContent(94,25);
   MET->SetBinContent(95,23);
   MET->SetBinContent(96,26);
   MET->SetBinContent(97,20);
   MET->SetBinContent(98,12);
   MET->SetBinContent(99,16);
   MET->SetBinContent(100,12);
   MET->SetBinContent(101,17);
   MET->SetBinContent(102,12);
   MET->SetBinContent(103,13);
   MET->SetBinContent(104,6);
   MET->SetBinContent(105,10);
   MET->SetBinContent(106,10);
   MET->SetBinContent(107,7);
   MET->SetBinContent(108,2);
   MET->SetBinContent(109,4);
   MET->SetBinContent(110,7);
   MET->SetBinContent(111,4);
   MET->SetBinContent(112,2);
   MET->SetBinContent(113,4);
   MET->SetBinContent(114,4);
   MET->SetBinContent(115,4);
   MET->SetBinContent(116,2);
   MET->SetBinContent(118,5);
   MET->SetBinContent(120,2);
   MET->SetBinContent(124,2);
   MET->SetBinContent(125,1);
   MET->SetBinContent(128,1);
   MET->SetBinContent(129,1);
   MET->SetBinContent(131,1);
   MET->SetBinContent(132,1);
   MET->SetBinContent(134,1);
   MET->SetBinContent(145,1);

  MET->Draw();

  RooRealVar peak("peak","peak",10,0,50);
  RooRealVar tail("tail","tail",0.1,-1,5);
  RooRealVar width("width","width",10,0,50);
  RooRealVar lambda( "lambda","lambda",0, -100., 0.);

  RooRealVar met("met","met",0,200,"GeV/c^{2}");

  RooNovosibirsk sig("novosibirsk", "novosibirsk", met, peak, width, tail );
  RooExponential bkg( "exponential","exponential", met, lambda );
  RooRealVar     f( "f","f", 0, 0., 1. );
  RooAddPdf      pdf( "pdf", "pdf", bkg, sig, f );  

  Chi2 chi2( "chi2", MET, met, pdf );  
  cout << "Done. Now fit." << endl;
  RooFitResult* fitResult = chi2.fit();
  //  chi2.evaluate();
  cout << "Done." << endl;
  fitResult->Print();

  TFile* dumFile = TFile::Open("dum.root", "RECREATE" );
  chi2.writeHist( dumFile );
}
