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
    double mT = P->gas[gasIndex]->mass*T;
    double C = 0.0;
    for(unsigned int impulseIndex;impulseIndex<P->impulse.size();impuleIndex++) {
        C += std::exp(-P->impulse[impulseIndex].lenght()/(2*mT));
    }
    C *= P->d3Impulse;
    C = 1.0/C;
    for(unsigned int impulseIndex;impulseIndex<P->impulse.size();impuleIndex++) {
        m_value[impulseIndex] = C*std::exp(-P->impulse[impulseIndex].lenght()/(2*mT));
    }
}

void Cell::ComputeHalf(Cell::Dimention dim)
{
    if(m_prev[dim]) {
        if(m_next[dim]) {
            if(m_next[dim]->m_next[dim]) {
                computeHalf_Normal(dim);
            } else {
                //computeHalf_PreRight(dim);
            }
        } else {
            //computeHalf_Right(dim);
        }
    } else {
        if(m_next[dim]) {
            //computeHalf_Left(dim);
        } else {
            std::cout << "This cell has no neightbourns on dim = <" << dim << ">" << std::endl;
        }
    }
}

void Cell::ComputeValue(Cell::Dimention dim)
{
    if(m_prev[dim]) {
        if(m_next[dim]) {
            if(m_next[dim]->m_next[dim]) {
                computeValue_Normal(dim);
            } else {
                //computeValue_PreRight(dim);
            }
        } else {
            //computeValue_Right(dim);
        }
    } else {
        if(m_next[dim]) {
            //computeValue_Left(dim);
        } else {
            std::cout << "This cell has no neightbourns on dim = <" << dim << ">" << std::endl;
        }
    }
}

void Cell::computeHalf_Normal(Dimention dim)
{
    double y = P->timestep/P->gas[gasIndex]->mass;
    for(unsigned int impulseIndex;impulseIndex<P->impulse.size();impuleIndex++) {
        y *= std::abs(P->impulse[impulseIndex][dim]/m_h[dim]);
        if(P->impulse[impulseIndex][dim] > 0) {
            m_half[impulseIndex] = m_value[impulseIndex] + (1-y)/2*limitter::super_bee(m_prev[dim]->m_value[impulseIndex],
                                                                                       m_value[impulseIndex],
                                                                                       m_next[dim]->m_value[impulseIndex]);
        } else {
            m_half[impulseIndex] = m_next[dim]->m_value[impulseIndex] - (1-y)/2*limitter::super_bee(m_value[impulseIndex],
                                                                                                    m_next[dim]->m_value[impulseIndex],
                                                                                                    m_next[dim]->m_next[dim]->m_value[impulseIndex]);
        }
    }
}

void Cell::computeValue_Normal(Cell::Dimention dim)
{
    double y = P->timestep/P->gas[gasIndex]->mass;
    for(unsigned int impulseIndex;impulseIndex<P->impulse.size();impuleIndex++) {
        y *= P->impulse[impulseIndex][dim]/m_h[dim];
        m_value[impulseIndex] += -y*(m_half[impulseIndex] - m_prev[dim]->m_half[impulseIndex]);
    }
}

