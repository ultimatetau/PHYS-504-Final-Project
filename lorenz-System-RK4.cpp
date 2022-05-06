#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>


struct LorenzState {
  
  double t;
  
  double x;
  double y;
  double z;

  double vx;
  double vy;
  double vz;

};

  LorenzState rhs(const LorenzState& state);
void write_history(const std::vector<LorenzState>& history);
std::vector<LorenzState> integrate(const double B, const double S, const double  P, const double tmax, const double dt);

const double S=10;
const double P=28;
const double B=8.0/3.0;


LorenzState rhs(const LorenzState& state)
{

  LorenzState dldt{};

  dldt.x=state.vx;
  dldt.y=state.vy;
  dldt.z=state.vz;

  dldt.x = S*(state.y-state.x);
  dldt.y = state.x*(P-state.z)-state.y;
  dldt.z = state.x*state.y-B*state.z;
  
  return dldt;

}


void write_history(const std::vector<LorenzState>& history,std::string str)
{

  std::ofstream of;

  of.open(str);
  for (auto o : history) {
      
      
      of << std::setw(12) << o.t
	 << std::setw(12) << o.x
	 << std::setw(12) << o.y
	 << std::setw(12) << o.z
	 << std::setw(12) << o.vx
	 << std::setw(12) << o.vy
	 << std::setw(12) <<o.vz
	 << std::endl;
      
    }
    of.close();

}


std::vector<LorenzState> integrate_rk4(const double a, const double b, const double c, const double tmax, const double dt)

{
  std::vector<LorenzState> o_history{};

  LorenzState state{};
  LorenzState state1{};
  LorenzState state2{};
  LorenzState state3{};
  LorenzState state4{};

  
  state.t = 0.0;
  state.x = a;
  state.y = b;
  state.z = c;



  state.vx = S*(state.y-state.x);
  state.vy = state.x*(P-state.z)-state.y;
  state.vz = state.x*state.y-B*state.z;

  o_history.push_back(state);

  while(state.t <tmax)

    {
      auto k0 = rhs(state);

        state1=state;
	
	state1.t += dt;
	state1.x += dt*k0.x;
	state1.y += dt*k0.y;
	state1.z += dt*k0.z;
	
	state1.vx+=k0.vx*dt;
	state1.vy+=k0.vy*dt;
	state1.vz+=k0.vz*dt;
	
	auto k1 = rhs(state1);

	state2 = state1;
	
	state2.t += .5*dt;
	state2.x += dt*k1.x;
	state2.y += dt*k1.y;   
	state2.z += dt*k1.z;
	
	state2.vx+=.5*k1.vx*dt;
	state2.vy+=.5*k1.vy*dt;
	state2.vz+=.5*k1.vz*dt;

	auto k2 = rhs(state2);

	state3 = state2;
	
	state3.t +=.5*dt;
	state3.x += dt*k2.x;
	state3.y += dt*k2.y;   
	state3.z += dt*k2.z;

	state3.vx+=.5*k2.vx*dt;
	state3.vy+=.5*k2.vy*dt;
	state3.vz+=.5*k2.vz*dt;

	
	auto k3 = rhs(state);
	state = state3;
	
	state.t +=dt;
	state.x += dt*k3.x;
	state.y += dt*k3.y;   
	state.z += dt*k3.z;

	state.vx+=k3.vx*dt;
	state.vy+=k3.vy*dt;
	state.vz+=k3.vz*dt;
	

	
	o_history.push_back(state);
	
      }

    return o_history;
      
      
      
      
    }


int main()
{

  double tmax;
  double dt;
  std::string his0 = "LorenzHistory.txt";
  double x0,y0,z0;
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
  
  auto orbit_history = integrate_rk4(x0,y0,z0, tmax, dt);
  write_history(orbit_history,his0);

    
}



