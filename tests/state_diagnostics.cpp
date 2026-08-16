#include <array>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using V3 = std::array<double, 3>;
using Q4 = std::array<double, 4>;
using M3 = std::array<double, 9>;

V3 add(V3 a, V3 b) { for (int i=0;i<3;++i) a[i]+=b[i]; return a; }
V3 mul(V3 a, double s) { for (double& x:a) x*=s; return a; }
double dot(V3 a, V3 b) { return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]; }
V3 cross(V3 a, V3 b) {
    return {a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2],
            a[0]*b[1]-a[1]*b[0]};
}
V3 matVec(const M3& m, V3 v) {
    return {m[0]*v[0]+m[1]*v[1]+m[2]*v[2],
            m[3]*v[0]+m[4]*v[1]+m[5]*v[2],
            m[6]*v[0]+m[7]*v[1]+m[8]*v[2]};
}
M3 bodyToWorld(Q4 q) {
    const double w=q[0],x=q[1],y=q[2],z=q[3];
    return {1-2*(y*y+z*z), 2*(x*y-w*z), 2*(x*z+w*y),
            2*(x*y+w*z), 1-2*(x*x+z*z), 2*(y*z-w*x),
            2*(x*z-w*y), 2*(y*z+w*x), 1-2*(x*x+y*y)};
}

struct Totals {
    double mass=0, kinetic=0, potential=0;
    V3 momentum{}, angular{}, force{}, torque{};
};

void includeBody(Totals& t, double mass, const M3& inertia, V3 position,
                 V3 velocity, V3 omegaBody, Q4 orientation, V3 gravity) {
    const M3 rotation=bodyToWorld(orientation);
    const V3 linearMomentum=mul(velocity,mass);
    const V3 spinWorld=matVec(rotation,matVec(inertia,omegaBody));
    t.mass+=mass;
    t.momentum=add(t.momentum,linearMomentum);
    t.angular=add(t.angular,add(cross(position,linearMomentum),spinWorld));
    t.kinetic+=0.5*mass*dot(velocity,velocity)+0.5*dot(omegaBody,matVec(inertia,omegaBody));
    t.potential-=mass*dot(gravity,position);
}

std::vector<std::string> splitCsv(const std::string& line) {
    std::vector<std::string> cells; std::stringstream in(line); std::string cell;
    while(std::getline(in,cell,',')) cells.push_back(cell);
    return cells;
}

int main(int argc,char** argv) {
    if(argc!=7){std::cerr<<"usage: state_diagnostics PARTICLES [BODY_STATE|-] [BODY_CSV|-] gx gy gz\n";return EXIT_FAILURE;}
    try {
        const V3 gravity{std::stod(argv[4]),std::stod(argv[5]),std::stod(argv[6])};
        Totals totals;
        std::ifstream particles(argv[1]);
        if(!particles) throw std::runtime_error("cannot open particle state");
        int number,status; double mass,radius; V3 pos,vel,omega; Q4 quat;
        while(particles>>number>>status>>mass>>radius
              >>pos[0]>>pos[1]>>pos[2]>>vel[0]>>vel[1]>>vel[2]
              >>quat[0]>>quat[1]>>quat[2]>>quat[3]
              >>omega[0]>>omega[1]>>omega[2]) {
            const double inertia=0.4*mass*radius*radius;
            includeBody(totals,mass,{inertia,0,0,0,inertia,0,0,0,inertia},
                        pos,vel,omega,quat,gravity);
        }
        if(std::string(argv[2])!="-") {
            std::ifstream csv(argv[3]); if(!csv) throw std::runtime_error("cannot open body CSV");
            std::string line; std::getline(csv,line); std::getline(csv,line);
            std::vector<double> masses; std::vector<M3> inertias;
            while(std::getline(csv,line)) { if(line.empty()) continue; auto c=splitCsv(line); if(c.size()<24) throw std::runtime_error("body CSV has fewer than 24 columns");
                masses.push_back(std::stod(c[1])); M3 I{}; for(int k=0;k<9;++k) I[k]=std::stod(c[15+k]); inertias.push_back(I); }
            std::ifstream bodies(argv[2]); if(!bodies) throw std::runtime_error("cannot open body state");
            for(std::size_t i=0;i<masses.size();++i){double time;V3 p,v,w,f,tau;Q4 q;
                if(!(bodies>>time>>p[0]>>p[1]>>p[2]>>v[0]>>v[1]>>v[2]
                    >>w[0]>>w[1]>>w[2]>>q[0]>>q[1]>>q[2]>>q[3]
                    >>f[0]>>f[1]>>f[2]>>tau[0]>>tau[1]>>tau[2])) throw std::runtime_error("body state row count differs from body CSV");
                includeBody(totals,masses[i],inertias[i],p,v,w,q,gravity);
                totals.force=add(totals.force,f); totals.torque=add(totals.torque,tau);
            }
        }
        std::cout<<std::scientific<<std::setprecision(12)
                 <<"mass="<<totals.mass<<"\n"
                 <<"momentum="<<totals.momentum[0]<<","<<totals.momentum[1]<<","<<totals.momentum[2]<<"\n"
                 <<"angular_momentum="<<totals.angular[0]<<","<<totals.angular[1]<<","<<totals.angular[2]<<"\n"
                 <<"kinetic="<<totals.kinetic<<"\n"
                 <<"potential="<<totals.potential<<"\n"
                 <<"mechanical="<<totals.kinetic+totals.potential<<"\n"
                 <<"body_force_sum="<<totals.force[0]<<","<<totals.force[1]<<","<<totals.force[2]<<"\n"
                 <<"body_torque_sum="<<totals.torque[0]<<","<<totals.torque[1]<<","<<totals.torque[2]<<"\n";
    }catch(const std::exception& e){std::cerr<<"state diagnostics failed: "<<e.what()<<"\n";return EXIT_FAILURE;}
    return EXIT_SUCCESS;
}
