void test2(){
  //  Double_t s = ;
  Double_t h = 4.5e-3;//4.76e-3; // height between plane, wire pitch ...
  Double_t r = 76e-6; // wire radius 
  Double_t d = 10e-6; // shift 


  Double_t x[28]={-3*h,-2*h,-1*h,0, h, 2*h, 3*h,   // G plane wires
		  -3*h,-2*h,-1*h,0, h, 2*h, 3*h, // U plane wires
		  -3*h,-2*h,-1*h,0, h, 2*h, 3*h, // V plane wires
		  -3*h,-2*h,-1*h,0, h, 2*h, 3*h // X plane wires
  };

  // // quick hack
  // h = 5.0e-3;
  // for (Int_t i=0;i!=7;i++){
  //   x[i] *=1.7;
  // }
  
  Double_t y[28]={4*h, 4*h, 4*h, 4*h, 4*h, 4*h, 4*h, // G plane wires
		  3*h, 3*h, 3*h, 3*h, 3*h, 3*h, 3*h, // U plane wires
		  2*h, 2*h, 2*h, 2*h, 2*h, 2*h, 2*h, // V plane wires
		  h,   h,   h,   h,   h,   h,   h // X plane wires
  };
  
  // // hack back ...
  // h = 4.5e-3;
  
  TMatrixD L(28,28);
  for (Int_t i=0;i!=28;i++){
    for (Int_t j=0;j!=28;j++){
      if (i==j){
        L(i,j) = log(2*y[i]/r);
      }else{
	L(i,j) = 0.5*log(1+4*y[i]*y[j]/(pow(x[i]-x[j],2)+pow(y[i]-y[j],2)));
      }
    }
  }

  TMatrixD L_p(28,28);
  //y[10] +=d;
  int offset = 0; //0 for U, 7 for V, 14 for X
  y[10+offset] += d;
  for (Int_t i=0;i!=28;i++){
    for (Int_t j=0;j!=28;j++){
      if (i==j){
        L_p(i,j) = log(2*y[i]/r);
      }else{
	L_p(i,j) = 0.5*log(1+4*y[i]*y[j]/(pow(x[i]-x[j],2)+pow(y[i]-y[j],2)));
      }
    }
  }

  L.Invert();
  L_p.Invert();

  TVectorD v1(28);
  TVectorD v2(28);
  for (Int_t i=0;i!=7;i++){
    v1(i) = -40;
    v1(7+i) = 0;
    v1(14+i) = 0;
    v1(21+i) = 0;
    
    v2(i) =  -40;
    v2(7+i) = -20;
    v2(14+i) = 0;
    v2(21+i) = 20;
  }

  
  TVectorD q1 = L * v1;
  TVectorD q1_p = L_p * v1;
  

  for (Int_t i=0;i!=7;i++){
    std::cout << q1(i+7+offset)-q1_p(i+7+offset) << " ";
  }
  std::cout << std::endl;

  TVectorD q2 = L * v2;
  TVectorD q2_p = L_p * v2;
   for (Int_t i=0;i!=7;i++){
    std::cout << q2(i+7+offset)-q2_p(i+7+offset) << " ";
  }
  std::cout << std::endl;
}
