#ifndef SPEED_H
#define SPEED_H

#include "header.h"

class Gas {
public:
    Gas(double getMass)
        :mass(getMass)
    {}
    double mass;

};

class Param
{
public:
    Param() {}
    static Param* instance() {
        static Param *s_Param = new Param();
        return s_Param;
    }

    void Init() {
        n = 20;
        double* spd= new double[n];
        for(unsigned int i=0;i<n;i++) {
            spd[i] = cut*(2.0*i/(n-1)-1);
        }
        for(unsigned int i1=0;i1<n;i1++)
           for(unsigned int i2=0;i2<n;i2++)
               for(unsigned int i3=0;i3<n;i3++) {
                   osg::Vec3d vec_speed(spd[i1],
                                        spd[i2],
                                        spd[i3]);
                   if (vec_speed.length() < cut) {
                       impulse.push_back(vec_speed);
                   }
               }
    }

    // Speed
    unsigned int n;
    double cut;
    std::vector<osg::Vec3d> impulse;
    double d3Impulse;
    double dImpulse;
    // Global Program Params
    std::vector<Gas*> gas;
    double timestep;
};

#define P Param::instance()

#endif // SPEED_H
