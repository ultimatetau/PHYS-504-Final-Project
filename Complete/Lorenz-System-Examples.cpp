#include "Lorenz.H"

int main()
{

  double tmax;
  double dt;
  std::string his0 = "LorenzHistory.txt";
  std::string his1 = "LorenzHistory1.txt";
  std::string his2 = "LorenzHistory2.txt";
  std::string his3 = "LorenzHistory3.txt";
  std::string his4 = "LorenzHistory4.txt";
  std::string his5 = "LorenzHistory5.txt";
  
  double x0,y0,z0,S,B,P;
  std::cout<<"Please enter an initial value for x(include decimal): ";
  std::cin >> x0;
  std::cout<<"Please enter an initial value for y(include decimal): ";
  std::cin >> y0;
  std::cout<<"Please enter an initial value for z(include decimal): ";
  std::cin >> z0;
  std::cout<<"Please enter a time step: ";
  std::cin >> dt;
  std::cout<<"Please enter a total time: ";
  std::cin >> tmax;
  
  std::cout<<"Please enter Sigma: ";
  std::cin >> S;

  std::cout<<"Please enter Rho: ";
  std::cin >> P;

  std::cout<<"Please enter Beta: ";
  std::cin >> B;

  std::cout<<"Test cases of the algorithm will be output as LorenzHistory files, which can be plotted using GNUplot> splot <filename> using 2:3:4"<<std::endl;
  
  auto orbit_history = integrate_rk4(x0,y0,z0,S,P,B, tmax, dt);
  write_history(orbit_history,his0);
  auto orbit_history1 = integrate_rk4(x0,y0,z0,S,P-15,B, tmax, dt);
  write_history(orbit_history1,his1);
  auto orbit_history2 = integrate_rk4(x0,y0,z0,10,28,2.66666, tmax, dt);
  write_history(orbit_history2,his2);
  auto orbit_history3 = integrate_rk4(x0+.01,y0,z0,10,28,2.66666, tmax, dt);
  write_history(orbit_history3,his3);

  auto orbit_history4 = integrate_rk4(x0,y0,z0,10,14,2.66666, tmax, dt);
  write_history(orbit_history4,his4);
  auto orbit_history5 = integrate_rk4(x0+.01,y0,z0,10,14,2.66666, tmax, dt);
  write_history(orbit_history5,his5);
  
}


