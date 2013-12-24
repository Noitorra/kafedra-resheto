#ifndef CELL_H
#define CELL_H

#include "param.h"

inline double exp_erg(const double& mass, const double& temp, const std::vector<double>& imp);

class Cell {
public:
    enum Dimention { X = 0,
                     Y = 1,
                     Z = 2 };
	// start var
    double T;

    std::vector<double> m_h;

    std::vector<Cell*> m_next;
    std::vector<Cell*> m_prev;

    std::vector<double*> m_value;
    std::vector<double*> m_half;
public:
    Cell();
    void Init();
    void InitBase();
    double getTemperature(int gas);
    double getDensity(int gas);

    void ComputeHalf(Dimention dim);
    void ComputeValue(Dimention dim);
    // Normal
    void computeHalf_Normal(Dimention dim);
    void computeValue_Normal(Dimention dim);
    // Left
    void computeHalf_Left(Dimention dim);
    void computeValue_Left(Dimention dim);
    // Right
    void computeHalf_Right(Dimention dim);
    void computeValue_Right(Dimention dim);
    // PreRight
    void computeHalf_PreRight(Dimention dim);
    void computeValue_PreRight(Dimention dim);
};

#endif // CELL_H
