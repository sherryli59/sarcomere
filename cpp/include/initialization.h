#ifndef INITIALIZATION_H
#define INITIALIZATION_H

#include "sarcomere.h"

namespace initialization {
    void partial_fix(Sarcomere& s, int& n_fixed);
    void cb(Sarcomere& s);
    void set_myosin_direction_x_noise(Sarcomere& s, double noise_std);
    void cb_off_angle(Sarcomere& s);
    void am_off_angle(Sarcomere& s);
    void sarcomeric_structure_tight(Sarcomere& s);
    void sarcomeric_structure(Sarcomere& s);
}

#endif // INITIALIZATION_H
