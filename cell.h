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
    void computeHalf(Dimention dim);
    void computeValue(Dimention dim);
    // start var
    double T;
    std::vector<double> m_h;

    std::vector<Cell*> m_next;
    std::vector<Cell*> m_prev;
    std::vector<double*> m_value;
    std::vector<double*> m_half;
};

#endif // CELL_H
