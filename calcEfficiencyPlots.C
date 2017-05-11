

const double zMax = 80;

//for the L1L1 dataset
double vertEffL1L1(double gct, double zz, double mA){
    return 1./gct*exp((-5-zz)/gct)*exp((6.67021E-5-3.23874E-6*mA+3.00751E-8*mA*mA)*pow(zz,3)
				     +(-0.00473355+0.000229457*mA-2.38217E-6*mA*mA)*pow(zz,2)
				      +(-0.00606069-0.00183954*mA+2.63933E-5*mA*mA)*zz
				      +(-3.28232+0.101364*mA-0.000815324*mA*mA));
    }

double fitZcutL1L1(double mass){
  if (mass<0.02){
    mass = 0.02;
  }
  else if (mass>0.06){
    mass = 0.06;
  }
  return -68.72+5738*mass-126929*pow(mass,2)+883250*pow(mass,3);
 
}
//for the L1L2 dataset
double vertEffL1L2(double gct, double zz, double mA){
  return 1./gct*exp((-5-zz)/gct)*(-0.36085+0.0241208*mA-0.000276979*pow(mA,2))*exp(-pow( ( zz-(-57.5354+3.5771*mA-0.0281698*pow(mA,2)))/(2*(3.82466+0.135087*mA)) ,2) );
    }

double fitZcutL1L2(double mass){
  //return (15.86 + 1857*mass - 62743.2*mass*mass + 536084*mass*mass*mass-5.0);
  if (mass<0.02){
    mass = 0.02;
  }
  else if (mass>0.06){
    mass = 0.06;
  }
  return -38.02+4895*mass-111334*pow(mass,2)+744842*pow(mass,3);
 
}

//for the L2L2 dataset
double vertEffL2L2(double gct, double zz, double mA){
  return 1./gct*exp((-5-zz)/gct)*(-0.745739+0.0513159*mA-0.00063263*pow(mA,2))*exp(-pow( ( zz-(-37.2262+4.70208*mA-0.0441197*pow(mA,2)))/(2*(5.33717+0.104304*mA)) ,2) );
}
double fitZcutL2L2(double mass){
  return -17.76 +8563*mass-293498*mass*mass+2678060*mass*mass*mass;
}


//here's the integral part
//mA in MeV, gct in mm, zCut in mm
double integrateZ0(double gct, double mA, double zCut){
  double step = 0.05;
  double integral = 0;
  for (double iz=zCut; iz<=zMax; iz+=step){
    integral += (step/6)*(vertEffL1L1(gct,iz,mA)+4*vertEffL1L1(gct,(2*iz+step)/2,mA)+vertEffL1L1(gct,iz+step,mA));
  }
  return integral;
}

double integrateZ1(double gct, double mA, double zCut){
  double step = 0.05;
  double integral = 0;
  for (double iz=zCut; iz<=zMax; iz+=step){
    integral += (step/6)*(vertEffL1L2(gct,iz,mA)+4*vertEffL1L2(gct,(2*iz+step)/2,mA)+vertEffL1L2(gct,iz+step,mA));
  }
  return integral;
}
double integrateZ2(double gct, double mA, double zCut){
  double step = 0.05;
  double integral = 0;
  for (double iz=zCut; iz<=zMax; iz+=step){
    integral += (step/6)*(vertEffL2L2(gct,iz,mA)+4*vertEffL2L2(gct,(2*iz+step)/2,mA)+vertEffL2L2(gct,iz+step,mA));
  }
  return integral;
}
//E0 [GeV], mA [MeV]
//returns [mm]
double gct(double E0, double eps2, double mA){
  return 8*(E0/10)*(pow(0.0001,2)/eps2)*pow(100/mA,2);
}


void calcEfficiencyPlots(){
  gROOT->SetBatch(true);

  //define by layer (L1L1=0,L1L2=2,L2L2=3)
  const int ntypes = 3;
  
  //loop over mass
  const int nmass = 90;
  double mass[nmass];
  for (int ii=0;ii<nmass;ii++){
    mass[ii] = ii;
  }
  
 

  //loop over eps2
  const int neps2 = 13;
  double eps2[neps2] = {7E-10, 8E-10, 9E-10, 1E-9, 2E-9, 3E-9, 4E-9, 5E-9, 6E-9, 7E-9, 8E-9, 9E-9, 1E-8};
  int eps2cheat[neps2] = {7,8,9,1,2,3,4,5,6,7,8,9,1};
  int eps2cheat2[neps2] = {10,10,10,9,9,9,9,9,9,9,9,9,8};
  double gctf[neps2][nmass];
  const int nzpos = 115;
  int zpos[nzpos];
  for (int jj=0; jj<nzpos; jj++){
    zpos[jj] = jj-5;//mm
  }
  
  //store the integral values
  double integral[ntypes][nzpos][neps2][nmass];  

  
  for(int kk=0; kk<nzpos; kk++){
    for (int ii=0; ii<nmass; ii++){
      for (int jj=0; jj<neps2;jj++){
	gctf[jj][ii] = gct(1.056,eps2[jj],mass[ii]);
	double zz = zpos[kk];
	integral[0][kk][jj][ii] = integrateZ0(gctf[jj][ii], mass[ii], zz);
	integral[1][kk][jj][ii] = integrateZ1(gctf[jj][ii], mass[ii], zz);
	integral[2][kk][jj][ii] = integrateZ2(gctf[jj][ii], mass[ii], zz);
      }
    }
  }


  
  //make canvas
  TCanvas *canvas = new TCanvas("canvas"," Integral Plots", 700, 700);
  canvas->SetFillColor(0);
  canvas->SetBorderMode(0);
  canvas->SetBorderSize(0);
  canvas->SetFrameFillColor(0);
  canvas->SetFrameBorderMode(0);
  //canvas->SetLeftMargin(0.1);
  canvas->SetRightMargin(0.2);
  std::string pdf_file_name = "../integralResults2.pdf";
  canvas->Update();
  gStyle->SetOptStat(0);

  ///////////////////////////////////////////////
  //Point for drawing the zCut
 
  const int nPts = 90;//55;//35;
  const int nPts2 = 90;//20;
  double cmassX[nPts];
  double zcutX[ntypes][nPts];
  for (int ijk=0; ijk<nPts; ijk++){
    cmassX[ijk] = ijk;
    zcutX[0][ijk] = fitZcutL1L1(cmassX[ijk]/1000);
    zcutX[1][ijk] = fitZcutL1L2(cmassX[ijk]/1000);
    cout<<"mass:\t"<<cmassX[ijk]<<"\t zCut:\t"<<zcutX[0][ijk]<<endl;
    
    //zcutX[2][ijk] = fitZcutL2L2(cmassX[ijk]/1000);
  }
  for (int ijk=0; ijk<nPts2; ijk++){
    cmassX[ijk] = ijk;//+15+0.5*ijk;
    zcutX[2][ijk] = fitZcutL2L2(cmassX[ijk]/1000);


  }
  


  //plots to make:  // double integral[ntypes][nzpos][neps2][nmass];  
  TH2F *h2_mVz[ntypes];
  TH2F *h2_epsVz[ntypes];
  const char *l1l1 = "L1L1, #epsilon^{2}=5E-9; vertex z position [mm]; mass [MeV];full integral where zCut = z vertex position";
  const char *l1l2 = "L1L2, #epsilon^{2}=5E-9; vertex z position [mm]; mass [Mev];full integral where zCut = z vertex position";
  const char *l2l2 = "L2L2, #epsilon^{2}=5E-9; vertex z position [mm]; mass [MeV];full integral where zCut = z vertex position";
  const char *l1l1e = "L1L1, mass=35MeV; vertex z position [mm]; #epsilon^{2};full integral where zCut = z vertex position";
  const char *l1l2e = "L1L2, mass=35MeV; vertex z position [mm]; #epsilon^{2};full integral where zCut = z vertex position";
  const char *l2l2e = "L2L2, mass=35MeV; vertex z position [mm]; #epsilon^{2};full integral where zCut = z vertex position";
 
  for (int kk=0; kk<ntypes; kk++){
    const char *layer;
    const char *layere;
    if (kk==0){layer = l1l1;layere = l1l1e;}
    else if (kk==1) {layer=l1l2;layere = l1l2e;}
    else {layer = l2l2;layere = l2l2e;}

    h2_mVz[kk] = new TH2F(Form("h2_mVz[%d]",kk),layer,117, -5, 100,100,0,90);
    h2_epsVz[kk] = new TH2F(Form("h2_epsVz[%d]",kk),layere,117, -5, 100,200,1E-10,4E-8);
    h2_mVz[kk]->GetZaxis()->SetRangeUser(1E-5,1E0); 
    h2_epsVz[kk]->GetZaxis()->SetRangeUser(1E-5,1E0); 
  }

 	
  //mass vs z, weighted by efficiency
  for (int ik=0; ik<ntypes; ik++){
    for (int iz=0; iz<nzpos; iz++){
      int zbin = h2_mVz[ik]->GetXaxis()->FindBin(zpos[iz]);

      for (int ii=0; ii<nmass;ii++){
	int mbin = h2_mVz[ik]->GetYaxis()->FindBin(mass[ii]);
	h2_mVz[ik]->SetBinContent(zbin,mbin,integral[ik][iz][7][ii]);
      }
      for (int jj=0; jj<neps2;jj++){
	int ebin = h2_epsVz[ik]->GetYaxis()->FindBin(eps2[jj]);
	h2_epsVz[ik]->SetBinContent(zbin,ebin,integral[ik][iz][jj][6]);
      }
      
    }
  }

  //gct (fix epsilon) vs z, weighted by efficiency

  //gct (fix mass) vs z, weighted by efficiency
  canvas->SetLogz();
  TGraph *gx[ntypes];

  for (int ik=0; ik<ntypes; ik++){

    h2_mVz[ik]->Draw("colz");
    h2_mVz[ik]->GetZaxis()->SetTitleOffset(1.5);
    if (ik!=2){
      gx[ik] = new TGraph(nPts,zcutX[ik],cmassX);
      gx[ik]->SetLineColor(kRed);
      gx[ik]->SetLineWidth(3);
      gx[ik]->Draw("lsame");
    }
    /*
    else{
      gx[ik] = new TGraph(nPts2,zcutX[ik],cmassX);
      gx[ik]->SetLineColor(kRed);
      gx[ik]->SetLineWidth(3);
      gx[ik]->Draw("lsame");
    }
    */
    canvas->Print( (pdf_file_name + "(").c_str());

  }
  canvas->SetLogy();

  for (int ik=0; ik<ntypes; ik++){

    h2_epsVz[ik]->Draw("colz");
    h2_epsVz[ik]->GetZaxis()->SetTitleOffset(1.5);
    canvas->Print( (pdf_file_name + "(").c_str());

  }
    
  canvas->Print( (pdf_file_name + ")").c_str());


}
