#ifndef CELL_H
#define CELL_H

#include "param.h"

class Cell {
public:
    enum Dimention { X = 0,
                     Y = 0,
                     Z = 0 };
public:
    Cell() {}
    void Init();

    void ComputeHalf(Dimention dim);
    void ComputeValue(Dimention dim);

    void computeHalf_Normal(Dimention dim);
    void computeValue_Normal(Dimention dim);
    // start var
    double T;
    unsigned int gasIndex;
    std::vector<double> m_h;

    std::vector<Cell*> m_next;
    std::vector<Cell*> m_prev;

    double* m_value;
    double* m_half;
};

#endif // CELL_H
