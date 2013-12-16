#include "cell.h"
#include <complex>

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

//void CMP_half(double& half, double& p_value, double& c_value, double& n_value, double& nn_value, double& h, unsigned int gasIndex) {

//}

inline double exp_erg(const double& mass, const double& temp, const double& imp) {
  return exp(-imp*imp/mass/2/temp);
}

// limitters
namespace limitter{
double super_bee(double x, double y, double z) {

    if( (z-y)*(y-x) <= 0 ) return 0.0;
    else return std::max(0.0, std::min(2*std::abs(y-x), std::min(std::abs(z-y), std::min(std::abs(y-x), 2*std::abs(z-y)))))*sgn(z-y);
}
}

void Cell::Init() {
    double mT = P->gas[gasIndex]->mass*T;
    double C = 0.0;
    for(unsigned int impulseIndex;impulseIndex<P->impulse->value.size();impulseIndex++) {
        C += exp_erg(P->gas[gasIndex]->mass, T, mod_vec(P->impulse->value[impulseIndex]));
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
                computeHalf_PreRight(dim);
            }
        } else {
            computeHalf_Right(dim);
        }
    } else {
        if(m_next[dim]) {
            computeHalf_Left(dim);
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
                computeValue_PreRight(dim);
            }
        } else {
            computeValue_Right(dim);
        }
    } else {
        if(m_next[dim]) {
            computeValue_Left(dim);
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

void Cell::computeHalf_Left(Cell::Dimention dim)
{
  double C1_up = 0.0;
  double C1_down = 0.0;
  double C2_up = 0.0;
  double y = P->timestep/P->gas[gasIndex]->mass;
  for(unsigned int impulseIndex;impulseIndex<P->impulse.size();impulseIndex++) {
    if(P->impulse[impulseIndex][dim] < 0) {
        y *= std::abs(P->impulse[impulseIndex][dim]/m_h[dim]);

        m_value[impulseIndex] = 2*m_next[dim]->m_next[dim]->m_value[impulseIndex] - m_next[dim]->m_value[impulseIndex];
        if (m_value[impulseIndex] < 0) m_value[impulseIndex] = 0;

        m_half[impulseIndex] = m_next[dim]->m_value[impulseIndex] - (1-y)/2*limitter::super_bee(m_value[impulseIndex],
                                                                                                m_next[dim]->m_value[impulseIndex],
                                                                                                m_next[dim]->m_next[dim]->m_value[impulseIndex]);

        C1_up += abs(P->impulse[impulseIndex][dim]*m_half[impulseIndex]);
        C2_up += abs(P->impulse[impulseIndex][dim]*(m_value[impulseIndex] + m_next[dim]->m_value[impulseIndex])/2);
    } else {
        C1_down += abs(P->impulse[impulseIndex][dim]*exp_erg(P->gas[gasIndex]->mass, T, P->impulse[impulseIndex][dim]));
    }
  }

  for(unsigned int impulseIndex;impulseIndex<P->impulse.size();impulseIndex++) {
    if(P->impulse[impulseIndex][dim] > 0) {
        m_half[impulseIndex] = C1_up/C1_down*exp_erg(P->gas[gasIndex]->mass, T, P->impulse[impulseIndex]);
        m_value[impulseIndex] = 2*C2_up/C1_down*exp_erg(P->gas[gasIndex]->mass, T, P->impulse[impulseIndex] - m_next[dim]->m_value[impulseIndex]);
        if(m_value[impulseIndex] < 0.0) m_value[impulseIndex] = 0.0;
    }
  }
}

void Cell::computeValue_Left(Cell::Dimention dim)
{
}

void Cell::computeHalf_Right(Cell::Dimention dim)
{
}

void Cell::computeValue_Right(Cell::Dimention dim)
{
}

void Cell::computeHalf_PreRight(Cell::Dimention dim)
{
}

void Cell::computeValue_PreRight(Cell::Dimention dim)
{
}

