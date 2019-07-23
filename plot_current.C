void plot_current(){
  Double_t h = 4.76e-3; // wire plance spacing
  Double_t r = 76e-6; // wire radius
  Double_t s[4] = {4.5e-3,4.89e-3,4.89e-3,4.5e-3}; // wire pitch
  
  Double_t d = 10e-6; // shift in position

  const int n = 5; // number of wires in a plane would be 2*n + 1,
  const int num = 2*n+1; // number of wires in a plane
  const int tot_num = 4*num; // number of wires in total

  Double_t x[tot_num];
  Double_t y[tot_num];

  // initialize the position of all wires
  for (Int_t i=0;i!=4;i++){

    x[i*num + n] = 0; // central wire is always at zero ... 
    y[i*num + n] = (4-i) * h;
    
    for (Int_t j=0;j!=n;j++){
      x[i*num+j] = -(n-j)*s[i];
      y[i*num+j] = (4-i) * h;
      
      x[i*num+j+n+1] = (j+1)*s[i];
      y[i*num+j+n+1] = (4-i)*h;
    }

    // for (Int_t j=0;j!=num;j++){
    //   if (i==0) x[i*num+j] += s[i]/2.;
    // }
  }

  // current L, initialize ...
  TMatrixD L(tot_num, tot_num);
  for (Int_t i=0;i!=tot_num;i++){
    for (Int_t j=0;j!=tot_num;j++){
      if (i==j){
        L(i,j) = log(2*y[i]/r);
      }else{
	L(i,j) = 0.5*log(1+4*y[i]*y[j]/(pow(x[i]-x[j],2)+pow(y[i]-y[j],2)));
      }
    }
  }


  TVectorD Vp(tot_num);
  TVectorD Vp_curr(tot_num);
  TVectorD Vp_prev(tot_num);
  Double_t potential[4]={-40,0,0,0}; //G, U, V, X
  for (Int_t i=0;i!=4;i++){
    for (Int_t j=0;j!=num;j++){
      Vp(i*num+j) = potential[i];
      Vp_curr(i*num+j) = potential[i];
      Vp_prev(i*num+j) = potential[i];
    }
  }

  
  // Now simulate the motion of the wire ...
  TMatrixD L_curr = L;
  TMatrixD L_prev = L;

  Int_t plane_no = 1;   //G, U, V, X
  Int_t wire_no = n; // central wire 

  const int nticks = 1000;

  TH1F **hh = new TH1F*[tot_num];
  for (Int_t i=0;i!=tot_num;i++){
    hh[i] = new TH1F(Form("hh_%d",i), Form("hh_%d",i), nticks,0,nticks);
  }
  
  
  
  for (Int_t i=0;i!=nticks+1;i++){
    // move two wires
    Int_t n1 = plane_no * num + wire_no;
    Double_t x1_pos = x[n1] + sin(i/100.*2*3.1415926)*0;
    Double_t y1_pos = y[n1] + sin(i/100.*2*3.1415926)*d;

    Int_t n2 = plane_no*num + wire_no + 1;
    Double_t x2_pos = x[n2] + sin(i/150.*2*3.1415926)*0;
    Double_t y2_pos = y[n2] + sin(i/150.*2*3.1415926)*0;
    
    if (i==0){
      for (Int_t j=0;j!=tot_num;j++){
	
	if (j == n1){
	  L_curr(j,j) = log(2*y1_pos/r);
	  L_prev(j,j) = log(2*y1_pos/r);
	  L_curr(j,n2) = 0.5*log(1+4*y2_pos*y1_pos/(pow(x2_pos-x1_pos,2)+pow(y2_pos-y1_pos,2)));
	}else if (j==n2){
	  L_curr(j,j) = log(2*y2_pos/r);
	  L_prev(j,j) = log(2*y2_pos/r);
	  L_curr(j,n1) = 0.5*log(1+4*y2_pos*y1_pos/(pow(x2_pos-x1_pos,2)+pow(y2_pos-y1_pos,2)));
	}else{
	  L_curr(n1,j) = 0.5*log(1+4*y1_pos*y[j]/(pow(x1_pos-x[j],2)+pow(y1_pos-y[j],2)));
	  L_curr(j,n1) =  L_curr(n1,j);

	  L_prev(n1,j) = 0.5*log(1+4*y1_pos*y[j]/(pow(x1_pos-x[j],2)+pow(y1_pos-y[j],2)));
	  L_prev(j,n1) =  L_prev(n1,j);

	  L_curr(n2,j) = 0.5*log(1+4*y2_pos*y[j]/(pow(x2_pos-x[j],2)+pow(y2_pos-y[j],2)));
	  L_curr(j,n2) =  L_curr(n2,j);

	  L_prev(n2,j) = 0.5*log(1+4*y2_pos*y[j]/(pow(x2_pos-x[j],2)+pow(y2_pos-y[j],2)));
	  L_prev(j,n2) =  L_prev(n2,j);
	}
      }
    }else{
      for (Int_t j=0;j!=tot_num;j++){
	if (j == n1){
	  L_curr(j,j) = log(2*y1_pos/r);
	  L_curr(j,n2) = 0.5*log(1+4*y2_pos*y1_pos/(pow(x2_pos-x1_pos,2)+pow(y2_pos-y1_pos,2)));
	}else if (j==n2){
	  L_curr(j,j) = log(2*y2_pos/r);
	  L_curr(j,n1) = 0.5*log(1+4*y2_pos*y1_pos/(pow(x2_pos-x1_pos,2)+pow(y2_pos-y1_pos,2)));
	}else{
	  L_curr(n1,j) = 0.5*log(1+4*y1_pos*y[j]/(pow(x1_pos-x[j],2)+pow(y1_pos-y[j],2)));
	  L_curr(j,n1) =  L_curr(n1,j);

	  L_curr(n2,j) = 0.5*log(1+4*y2_pos*y[j]/(pow(x2_pos-x[j],2)+pow(y2_pos-y[j],2)));
	  L_curr(j,n2) =  L_curr(n2,j);
	}
      }
    }

    TMatrixD L_curr_inv = L_curr;
    TMatrixD L_prev_inv = L_prev;
    L_curr_inv.Invert();
    L_prev_inv.Invert();

    
    // for (Int_t j=0;j!=num;j++){
    //   Vp_curr(2*num+j) = potential[2] + sin(i/150.*2.*3.1415926);
    // }
    
    
    TVector dq = (L_prev_inv * Vp_prev - L_curr_inv * Vp_curr);
    if (i!=0){
      for (Int_t j=0;j!=tot_num;j++){
	hh[j]->SetBinContent(i,dq(j));
      }
      //      hc->SetBinContent(i,dq(plane_no * num + wire_no));
      //      hp->SetBinContent(i,dq(plane_no * num + wire_no + 1));
      //   hn->SetBinContent(i,dq(plane_no * num + wire_no - 1));
    }
    L_prev = L_curr;
    Vp_prev = Vp_curr;
  }
  
  TCanvas *c1 = new TCanvas("c1","c1",1200,600);
  c1->Divide(2,1);
  c1->cd(1);
  gStyle->SetOptStat(0);

  hh[(plane_no) * num + wire_no]->Draw();
  hh[(plane_no) * num + wire_no]->SetLineColor(1);
  hh[(plane_no) * num + wire_no+1]->Draw("same");
  hh[(plane_no) * num + wire_no+1]->SetLineColor(2);
  hh[(plane_no) * num + wire_no+2]->Draw("same");
  hh[(plane_no) * num + wire_no+2]->SetLineColor(6);
  hh[(plane_no) * num + wire_no-1]->Draw("same");
  hh[(plane_no) * num + wire_no-1]->SetLineColor(4);
  
  Double_t max = hh[plane_no * num + wire_no]->GetMaximum(); Double_t min = hh[plane_no * num + wire_no]->GetMinimum();
  if (max <  hh[plane_no * num + wire_no+1]->GetMaximum()){
    max = hh[plane_no * num + wire_no+1]->GetMaximum();
  }
  if (min > hh[plane_no * num + wire_no+1]->GetMinimum()){
    min = hh[plane_no * num + wire_no+1]->GetMinimum();
  }
  hh[plane_no * num + wire_no]->GetYaxis()->SetRangeUser(min*1.2, max*1.2);
  

  TLegend *le1 = new TLegend(0.6,0.6,0.89,0.89);
  le1->SetFillColor(0);
  le1->AddEntry(hh[plane_no * num + wire_no],"C wire","l");
  le1->AddEntry(hh[plane_no * num + wire_no+1],"C+1 wire","l");
  le1->AddEntry(hh[plane_no * num + wire_no+2],"C+2 wire","l");
  le1->AddEntry(hh[plane_no * num + wire_no-1],"C-1 wire","l");
  
  
  // hc->Draw();
  // hc->SetLineColor(1);
  // hp->Draw("same");
  // hp->SetLineColor(2);
  // hn->Draw("same");
  // hn->SetLineColor(6);
  c1->cd(2);
  TH1 *hm = hh[plane_no * num + wire_no]->FFT(0,"MAG");
  hm->Draw();
  hm->SetLineColor(1);
  TH1 *hm1 = hh[plane_no * num + wire_no+1]->FFT(0,"MAG");
  hm1->Draw("same");
  hm1->SetLineColor(2);
  TH1 *hm2 = hh[plane_no * num + wire_no+2]->FFT(0,"MAG");
  hm2->Draw("same");
  hm2->SetLineColor(6);
  TH1 *hm3 = hh[plane_no * num + wire_no-1]->FFT(0,"MAG");
  hm3->Draw("same");
  hm3->SetLineColor(4);
  hm->SetLineStyle(1);
  hm1->SetLineStyle(2);
  hm2->SetLineStyle(3);
  hm3->SetLineStyle(4);
  
  le1->Draw();

  max = hm->GetMaximum();
  min = hm->GetMinimum();
  if (max <  hm1->GetMaximum()){
    max = hm1->GetMaximum();
  }
  if (min > hm1->GetMinimum()){
    min = hm1->GetMinimum();
  }
  hm->GetYaxis()->SetRangeUser(min*0.5,max*1.2);


  //hm1->Scale(8);
  
  // TH1 *hm_c = hc->FFT(0,"MAG");
  // TH1 *hm_p = hp->FFT(0,"MAG");
  // TH1 *hm_n = hn->FFT(0,"MAG");

  // hm_c->Draw();
  // hm_p->SetLineColor(2);
  // hm_p->Draw("same");
  // hm_n->SetLineColor(6);
  // hm_n->Draw("same");
  


  
  
  
}
