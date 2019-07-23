void test1(){
  Double_t h = 4.76e-3;
  Double_t r = 76e-6;
  Double_t h3 = 3*h;
  Double_t h2 = 3*h;
  Double_t h1 = 4*h;
  Double_t s = 4.76e-3;
  Double_t d = 50e-6;

  gStyle->SetOptStat(0);
  
  TMatrixD L(3,3);
  L(0,0) = log(2*h1/r)/100.;
  L(1,1) =log(2*h2/r);
  L(2,2) = log(2*h3/r);

  L(1,2) = 0.5 * log(1+4*h2*h3/s/s);
  L(2,1) = L(1,2);

  L(0,1) = 0.5*log(1+4*h1*h2/(h1-h2)/(h1-h2))/100.; //log(2*(h1-h2)/r) + log(2*(h1+h2)/r);
  L(1,0) = L(0,1);

  L(0,2) = 0.5*log(1+4*h1*h3/(h1-h3)/(h1-h3))/100.; //log(2*(h1-h3)/r) + log(2*(h1+h3)/r);
  L(2,0) = L(0,2);

  // L.Draw("COLZ");

  TMatrixD L_p1(3,3);
  
  // horizontal move ...
  L_p1 = L;
  L_p1(1,2) = 0.5 * log(1+4*h2*h3/(s-d)/(s-d));
  L_p1(2,1) = L_p1(1,2);

  //  vertical move ...
  h2 = h2-d;
  L_p1(0,0) = log(2*h1/r)/100;
  L_p1(1,1) =log(2*h2/r);
  L_p1(2,2) = log(2*h3/r);

  L_p1(1,2) = 0.5 * log(1+4*h2*h3/s/s);
  L_p1(2,1) = L_p1(1,2);

  L_p1(0,1) = 0.5*log(1+4*h1*h2/(h1-h2)/(h1-h2))/100; //log(2*(h1-h2)/r) + log(2*(h1+h2)/r);
  L_p1(1,0) = L(0,1);

  L_p1(0,2) = 0.5*log(1+4*h1*h3/(h1-h3)/(h1-h3))/100; //log(2*(h1-h3)/r) + log(2*(h1+h3)/r);
  L_p1(2,0) = L_p1(0,2);

  
  

  
  L.Invert();
  L_p1.Invert();
  TVectorD v1(3);
  v1(0) = 0;
  v1(1) = 10;
  v1(2) = 10;
  
  TVectorD v2(3);
  v2(0) = -10;
  v2(1) = 0;
  v2(2) = 0;

  TVectorD q1 = L * v1;
  TVectorD q1_p = L_p1 * v1;
  std::cout << q1(0)-q1_p(0) << " " << q1(1)-q1_p(1) << " " << q1(2)-q1_p(2) << std::endl;
  // std::cout <<  << " " <<  << " " <<  << std::endl;

  TVectorD q2 = L * v2;
  TVectorD q2_p = L_p1 * v2;
  std::cout << q2(0)-q2_p(0) << " " << q2(1)- q2_p(1) << " " << q2(2)-q2_p(2) << std::endl;
  //std::cout <<  << " " << << " " <<  << std::endl;
}
