void test(){
  TMatrixD L(3,3);
  TMatrixD L_p(3,3);
  Double_t L_init[3][3]={{100,10,10},
			{10,30,20},
			{10,20,30}};
  for (Int_t i=0;i!=3;i++){
    for (Int_t j=0;j!=3;j++){
      L(i,j) = L_init[i][j];
      L_p(i,j) = L_init[i][j];
    }
  }
  // horizontal change ...
  L_p(1,2) -= 5;
  L_p(2,1) -= 5;
  
  // vertical charnge ...
  //L_p(0,2)-=5;
  //L_p(2,0)-=5;
  
  L.Invert();
  L_p.Invert();
  TVectorD v1(3);
  v1(0) = 0;
  v1(1) = 10;
  v1(2) = 10;
  TVectorD v2(3);
  v2(0) = -10;
  v2(1) = 0;
  v2(2) = 0;

  TVectorD q1 = L * v1;
  TVectorD q1_p = L_p * v1;
  std::cout << q1(0) << " " << q1(1) << " " << q1(2) << std::endl;
  std::cout << q1_p(0) << " " << q1_p(1) << " " << q1_p(2) << std::endl;

  TVectorD q2 = L * v2;
  TVectorD q2_p = L_p * v2;
  std::cout << q2(0) << " " << q2(1) << " " << q2(2) << std::endl;
  std::cout << q2_p(0) << " " << q2_p(1) << " " << q2_p(2) << std::endl;
  //  L.Draw("COLZ");
  
}
