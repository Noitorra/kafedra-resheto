#include "cell.h"
#include <complex>

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

//void CMP_half(double& half, double& p_value, double& c_value, double& n_value, double& nn_value, double& h, unsigned int gasIndex) {

//}
inline double exp_erg(const double& mass, const double& temp, std::vector<double>& imp) {
	return exp(-mod_vec(imp)/mass/2/temp);
}

// limitters
namespace limitter{
double super_bee(double x, double y, double z) {
  //return 0.0;
    if( (z-y)*(y-x) <= 0 ) return 0.0;
    else return std::max(0.0, std::min(2*std::abs(y-x), std::min(std::abs(z-y), std::min(std::abs(y-x), 2*std::abs(z-y)))))*sgn(z-y);
}
}

Cell::Cell() {
  m_prev.push_back(NULL);
  m_prev.push_back(NULL);
  m_next.push_back(NULL);
  m_next.push_back(NULL);

  m_half.resize(P->gas.size(), NULL);
  m_value.resize(P->gas.size(), NULL);
}

void Cell::Init() {
  for(int gasIndex=0;gasIndex<P->gas.size();gasIndex++) {
    double C = 0.0;
    for(unsigned int impulseIndex=0;impulseIndex<P->impulse->value.size();impulseIndex++) {
        C += exp_erg(P->gas[gasIndex]->mass, T, P->impulse->value[impulseIndex]);
    }

    C *= P->impulse->d3P;
    C = 1.0/C;

    // Allocating space for values and half's
    m_half[gasIndex] = new double[P->impulse->value.size()];
    m_value[gasIndex] = new double[P->impulse->value.size()];

          for(unsigned int impulseIndex=0;impulseIndex<P->impulse->value.size();impulseIndex++) {
                  m_value[gasIndex][impulseIndex] = C*exp_erg(P->gas[gasIndex]->mass, T, P->impulse->value[impulseIndex]);
      //m_half[impulseIndex] = 0.0;
    }
  }
}

void Cell::InitBase() {
  m_half = new double[P->impulse->value.size()];
  m_value = new double[P->impulse->value.size()];
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
      //std::cout << "This cell has no neightbourns on dim = <" << dim << ">" << std::endl;
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
        computeValue_Normal(dim);
      }
    } else {
      computeValue_Right(dim);
    }
  } else {
    if(m_next[dim]) {
      computeValue_Left(dim);
    } else {
      //std::cout << "This cell has no neightbourns on dim = <" << dim << ">" << std::endl;
    }
  }
  //for(unsigned int impulseIndex=0;impulseIndex<P->impulse->value.size();impulseIndex++) {
  //  if(abs(m_value[impulseIndex]) > 1.0) {
  //    std::cout << m_value[impulseIndex] << " on " << impulseIndex << std::endl; 
  //  }
  //}
}

void Cell::computeHalf_Normal(Dimention dim)
{
  for(unsigned int impulseIndex=0;impulseIndex<P->impulse->value.size();impulseIndex++) {
		double y = P->timestep/P->gas[gasIndex]->mass*std::abs(P->impulse->value[impulseIndex][dim]/m_h[dim]);
    if(P->impulse->value[impulseIndex][dim] > 0) {
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
    for(unsigned int impulseIndex=0;impulseIndex<P->impulse->value.size();impulseIndex++) {
		    double y = P->timestep/P->gas[gasIndex]->mass*P->impulse->value[impulseIndex][dim]/m_h[dim];
        m_value[impulseIndex] = m_value[impulseIndex] - y*(m_half[impulseIndex] - m_prev[dim]->m_half[impulseIndex]);
    }
}

void Cell::computeHalf_Left(Cell::Dimention dim)
{
  double C1_up = 0.0;
  double C1_down = 0.0;
  double C2_up = 0.0;
  for(unsigned int impulseIndex=0;impulseIndex<P->impulse->value.size();impulseIndex++) {
    if(P->impulse->value[impulseIndex][dim] < 0) {
		    double y = P->timestep/P->gas[gasIndex]->mass*std::abs(P->impulse->value[impulseIndex][dim]/m_h[dim]);

        m_value[impulseIndex] = 2*m_next[dim]->m_value[impulseIndex] - m_next[dim]->m_next[dim]->m_value[impulseIndex];
        if (m_value[impulseIndex] < 0) m_value[impulseIndex] = 0;

        m_half[impulseIndex] = m_next[dim]->m_value[impulseIndex] - (1-y)/2*limitter::super_bee(m_value[impulseIndex],
                                                                                                m_next[dim]->m_value[impulseIndex],
                                                                                                m_next[dim]->m_next[dim]->m_value[impulseIndex]);

        C1_up += abs(P->impulse->value[impulseIndex][dim]*m_half[impulseIndex]);
        C2_up += abs(P->impulse->value[impulseIndex][dim]*(m_value[impulseIndex] + m_next[dim]->m_value[impulseIndex])/2);
    } else {
        C1_down += abs(P->impulse->value[impulseIndex][dim]*exp_erg(P->gas[gasIndex]->mass, T, P->impulse->value[impulseIndex]));
    }
  }

  for(unsigned int impulseIndex=0;impulseIndex<P->impulse->value.size();impulseIndex++) {
    if(P->impulse->value[impulseIndex][dim] > 0) {
        m_half[impulseIndex] = C1_up/C1_down*exp_erg(P->gas[gasIndex]->mass, T, P->impulse->value[impulseIndex]);
		    m_value[impulseIndex] = 2*C2_up/C1_down*exp_erg(P->gas[gasIndex]->mass, T, P->impulse->value[impulseIndex]) - m_next[dim]->m_value[impulseIndex];
        if(m_value[impulseIndex] < 0.0) m_value[impulseIndex] = 0.0;
    }
  }
}

void Cell::computeValue_Left(Cell::Dimention dim)
{
}

void Cell::computeHalf_Right(Cell::Dimention dim)
{
	double C1_up = 0.0;
  double C1_down = 0.0;
  double C2_up = 0.0;
  for(unsigned int impulseIndex=0;impulseIndex<P->impulse->value.size();impulseIndex++) {
    if(P->impulse->value[impulseIndex][dim] > 0) {
		    double y = P->timestep/P->gas[gasIndex]->mass*std::abs(P->impulse->value[impulseIndex][dim]/m_h[dim]);

        m_value[impulseIndex] = 2*m_prev[dim]->m_value[impulseIndex] - m_prev[dim]->m_prev[dim]->m_value[impulseIndex];
        if (m_value[impulseIndex] < 0) m_value[impulseIndex] = 0;

        m_prev[dim]->m_half[impulseIndex] = m_prev[dim]->m_value[impulseIndex] + (1-y)/2*limitter::super_bee(m_prev[dim]->m_prev[dim]->m_value[impulseIndex],
                                                                                                             m_prev[dim]->m_value[impulseIndex],
                                                                                                             m_value[impulseIndex]);

        C1_up += abs(P->impulse->value[impulseIndex][dim]*m_prev[dim]->m_half[impulseIndex]);
        C2_up += abs(P->impulse->value[impulseIndex][dim]*(m_value[impulseIndex] + m_prev[dim]->m_value[impulseIndex])/2);
    } else {
        C1_down += abs(P->impulse->value[impulseIndex][dim]*exp_erg(P->gas[gasIndex]->mass, T, P->impulse->value[impulseIndex]));
    }
  }

  for(unsigned int impulseIndex=0;impulseIndex<P->impulse->value.size();impulseIndex++) {
    if(P->impulse->value[impulseIndex][dim] < 0) {
        m_prev[dim]->m_half[impulseIndex] = C1_up/C1_down*exp_erg(P->gas[gasIndex]->mass, T, P->impulse->value[impulseIndex]);
		    m_value[impulseIndex] = 2*C2_up/C1_down*exp_erg(P->gas[gasIndex]->mass, T, P->impulse->value[impulseIndex]) - m_prev[dim]->m_value[impulseIndex];
        if(m_value[impulseIndex] < 0.0) m_value[impulseIndex] = 0.0;
    }
  }
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

// Macroparams
double Cell::getTemperature() {
  return T;
}

double Cell::getDensity() {
  if(!m_value) return 0.0;

  double density = 0.0;
  for(unsigned int impulseIndex=0;impulseIndex<P->impulse->value.size();impulseIndex++) {
    density += m_value[impulseIndex];
  }
  density *= P->impulse->d3P;
  return density;
}
