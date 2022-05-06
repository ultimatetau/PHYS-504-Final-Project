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
	
	state1.vx+= .25*k0.vx*dt;
	state1.vy+=.25*k0.vy*dt;
	state1.vz+=.25*k0.vz*dt;
	
	auto k1 = rhs(state1);
	
	
	state2.t += dt;
	state2.x += dt*k1.x;
	state2.y += dt*k1.y;  
	state2.z += dt*k1.z;
	
	state2.vx+=k1.vx*dt;
	state2.vy+=k1.vy*dt;
	state2.vz+=k1.vz*dt;

	auto k2 = rhs(state2);

	state.t += .5*dt;
	state.x += .5*dt*k2.x;
	state.y += .5*dt*k2.y;
	state.z += .5*dt*k2.z;
	
	state.vx+=k2.vx*dt;
	state.vy+=k2.vy*dt;
	state.vz+=k2.vz*dt;

	//auto k3 = rhs(state3);

	//state.t += dt;
	//state.x += .5*dt*k3.x;
	//state.y += .5*dt*k3.y;
	//state.z += .5*dt*k3.z;
	
	//state.vx+= .5*k3.vx*dt;
	//state.vy+=.5*k3.vy*dt;
	//state.vz+=.5*k3.vz*dt;
	
	o_history.push_back(state);
	
      }

    return o_history;
      
      
      
      
    }


int main()
{

 double tmax = 100.0;
    double dt = 0.01;
    std::string his0 = "LorenzHistory.txt";

    

    auto orbit_history = integrate_rk4(1,1,1, tmax, dt);
    write_history(orbit_history,his0);

    
}



