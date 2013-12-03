#include "cell.h"
#include <complex>

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

//void CMP_half(double& half, double& p_value, double& c_value, double& n_value, double& nn_value, double& h, unsigned int gasIndex) {

//}

double expDen(double mass, double temp, double speed) {
    return exp(-mass*speed*speed/2/temp);
}

// limitters
namespace limitter{
double super_bee(double x, double y, double z) {

    if( (z-y)*(y-x) <= 0 ) return 0.0;
    else return std::max(0, std::min(2*std::abs(y-x), std::min(std::abs(z-y), std::min(std::abs(y-x), 2*std::abs(z-y)))))*sgn(z-y);
}
}

void Cell::Init() {
    for(unsigned int gasIndex=0;gasIndex<P->gas.size();gasIndex++) {
        double mT = P->gas[gasIndex]->mass*T;
        double C = 0.0;
        for(unsigned int impulseIndex;impulseIndex<P->impulse.size();impuleIndex++) {
            C += std::exp(-P->impulse[impulseIndex].lenght()/(2*mT));
        }
        C *= P->d3Impulse;
        C = 1.0/C;
        for(unsigned int impulseIndex;impulseIndex<P->impulse.size();impuleIndex++) {
            m_value[gasIndex][impulseIndex] = C*std::exp(-P->impulse[impulseIndex].lenght()/(2*mT));
        }
    }
}

void Cell::computeHalf(Dimention dim)
{
    for(unsigned int gasIndex=0;gasIndex<P->gas.size();gasIndex++) {
        double y = P->timestep/P->gas[gasIndex]->mass;
        for(unsigned int impulseIndex;impulseIndex<P->impulse.size();impuleIndex++) {
            y *= std::abs(P->impulse[impulseIndex][dim]/m_h[dim]);
            if(P->impulse[impulseIndex][dim] > 0) {
                m_half[gasIndex][impulseIndex] = m_value[gasIndex][impulseIndex] + (1-y)/2*limitter::super_bee(m_prev[dim]->m_value[gasIndex][impulseIndex],
                                                                                           m_value[gasIndex][impulseIndex],
                                                                                           m_next[dim]->m_value[gasIndex][impulseIndex]);
            } else {
                m_half[gasIndex][impulseIndex] = m_next[dim]->m_value[gasIndex][impulseIndex] - (1-y)/2*limitter::super_bee(m_value[gasIndex][impulseIndex],
                                                                                                        m_next[dim]->m_value[gasIndex][impulseIndex],
                                                                                                        m_next[dim]->m_next[dim]->m_value[gasIndex][impulseIndex]);
            }
        }
    }
}

void Cell::computeValue(Cell::Dimention dim)
{
    for(unsigned int gasIndex=0;gasIndex<P->gas.size();gasIndex++) {
        double y = P->timestep/P->gas[gasIndex]->mass;
        for(unsigned int impulseIndex;impulseIndex<P->impulse.size();impuleIndex++) {
            y *= P->impulse[impulseIndex][dim]/m_h[dim];
            m_value[gasIndex][impulseIndex] += -y*(m_half[gasIndex][impulseIndex] - m_prev[dim]->m_half[gasIndex][impulseIndex]);
        }
    }
}

