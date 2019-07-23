void calc(){
  Double_t r = 76e-6; // m, wire radius
  Double_t h = 4.76e-3; // m, wire plane gap
  Double_t s = 4.76e-3; // m, wire pitech
  Double_t d = 10e-6; // m, displacement

  // Nominal case (no displacement) two wires close to each other ...
  Double_t l11 = log(2*h/r);
  Double_t dl = log(1+4*h*h/s/s)/2.;
  Double_t l22 = log(2*h/r);

  Double_t q1_n = (l22-dl)/(l11*l22-dl*dl);
  Double_t q2_n = (l11-dl)/(l11*l22-dl*dl);

  std::cout << q1_n << " " << q2_n << std::endl;

  // Vertical displacement for wire 2;
  l11 = log(2*h/r);
  l22 = log(2*(h-d)/r);
  dl = log(1+4*h*(h-d)/s/s)/2.;
  Double_t q1_v = (l22-dl)/(l11*l22-dl*dl);
  Double_t q2_v = (l11-dl)/(l11*l22-dl*dl);

  std::cout << (q1_v-q1_n)/q1_n << " " << (q2_v-q2_n)/q2_n << std::endl;

  // horizontal displacement for wire 2;
  l11 = log(2*h/r);
  l22 = log(2*h/r);
  dl = log(1+4*h*h/(s-d)/(s-d))/2.;
  Double_t q1_h = (l22-dl)/(l11*l22-dl*dl);
  Double_t q2_h = (l11-dl)/(l11*l22-dl*dl);

  std::cout << (q1_h-q1_n)/q1_n << " " << (q2_h-q2_n)/q2_n << std::endl;
}
