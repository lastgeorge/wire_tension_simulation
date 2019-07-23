void calc_multiple(){
  gStyle->SetOptStat(0);
  Double_t r = 76e-6; // m, wire radius
  Double_t h = 4.76e-3; // m, wire plane gap
  Double_t s = 4.76e-3; // m, wire pitech
  Double_t d = 50e-6; // m, displacement
  
  // 0 --> 10
  const int n = 15;
  // move the middle wire, e.g. 5 if n = 11
  int middle_n = (n-1)/2;

  TMatrixD L_n(n,n);
  TMatrixD L_v(n,n);
  TMatrixD L_h(n,n);

  for (Int_t i=0;i!=n;i++){
    for (Int_t j=0;j!=n;j++){
      if (i==j){
	L_n(i,j) = log(2*h/r);
	L_h(i,j) = log(2*h/r);
	if (i==middle_n){
	  L_v(i,j) = log(2*(h-d)/r);
	}else{
	  L_v(i,j) = log(2*h/r);
	}
      }else{
	L_n(i,j) = log(1+4*h*h/s/s/pow(i-j,2))/2.;
	if (i==middle_n || j==middle_n){
	  L_v(i,j) = log(1+4*h*(h-d)/s/s/pow(i-j,2))/2.;
	  L_h(i,j) = log(1+4*h*h/pow((i-j)*s-d,2))/2.;
	}else{
	  L_v(i,j) = log(1+4*h*h/s/s/pow(i-j,2))/2.;
	  L_h(i,j) = log(1+4*h*h/s/s/pow(i-j,2))/2.;
	}
      }
    }
  }
  L_n.Invert();
  L_v.Invert();
  L_h.Invert();

  Double_t Q_n[n], Q_v[n], Q_h[n];
  for (Int_t i=0;i!=n;i++){
    Q_n[i] = 0;
    Q_v[i] = 0;
    Q_h[i] = 0;
  }
  for (Int_t i=0;i!=n;i++){
    for (Int_t j=0;j!=n;j++){
      Q_n[i] += L_n(i,j);
      Q_v[i] += L_v(i,j);
      Q_h[i] += L_h(i,j);
    }
  }
  TGraph *g1 = new TGraph();
  TGraph *g2 = new TGraph();
  for (Int_t i=0;i!=n;i++){
    g1->SetPoint(i, i+1, (Q_v[i]-Q_n[i])/Q_n[i] );
    g2->SetPoint(i, i+1, (Q_h[i]-Q_n[i])/Q_n[i] );
  }
  TCanvas *c1 = new TCanvas("c1","c1",1200,600);
  c1->Divide(2,1);
  c1->cd(1);
  g1->Draw("A*");
  c1->cd(2);
  g2->Draw("A*");
  g1->SetTitle("Vertical");
  g2->SetTitle("Horizontal");
  //L_v.Invert();
  //  L_v.Draw("COLZ");
  
}
