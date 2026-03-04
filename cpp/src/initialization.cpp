#include "initialization.h"
#include <vector>
#include <cmath>
#include <cstdio>
#include <gsl/gsl_randist.h>

// Moved initialization/setup methods from sarcomere.cpp

void Sarcomere::set_bundling_parameters(double max_strength, int ramp_steps) {
    k_bundle_max = (max_strength > 0.0) ? max_strength : 0.0;
    bundle_ramp_steps = (ramp_steps > 0) ? ramp_steps : 0;
}

void Sarcomere::partial_fix(int& n_fixed){
    std::vector<vector> myosin_positions;
    myosin_positions = {
        {0, -2, 0}, {0, -1, 0}, {0, 0, 0}, {0, 1, 0}, {0, 2, 0},
        {0, -2.5, 0}, {0, -1.5, 0}, {0, -0.5, 0}, {0, 0.5, 0}, {0, 1.5, 0}
    };
    for (int i = 0; i < n_fixed; i++){
        myosin.center[i].x = myosin_positions[i][0];
        myosin.center[i].y = myosin_positions[i][1];
        myosin.center[i].z = myosin_positions[i][2];
    }
    for (int i = 0; i < myosin.n; i++){
        myosin.direction[i] = {1, 0, 0};
    }
    myosin.update_endpoints();
    update_system();
}

void Sarcomere::cb(){
    std::vector<vector> actin_positions = {
        {0.5, -0.015, 0}, {-0.5, 0.015, 0}
    };
    for (int i = 0; i < actin_positions.size(); i++){
        actin.center[i].x = actin_positions[i][0];
        actin.center[i].y = actin_positions[i][1];
        actin.center[i].z = actin_positions[i][2];
    }
    actin.direction[0] = {1, 0, 0};
    actin.direction[1] = {-1, 0, 0};
    
    std::vector<vector> myosin_positions = {
        {-1.33, 0.045, 0}, {-1.33, -0.015, 0}, {1.33, 0.015, 0}, {1.33, -0.045, 0}
    };
    for (int i = 0; i < myosin_positions.size(); i++){
        myosin.center[i].x = myosin_positions[i][0];
        myosin.center[i].y = myosin_positions[i][1];
        myosin.center[i].z = myosin_positions[i][2];
        myosin.direction[i] = {1, 0, 0};
    }
}

void Sarcomere::set_myosin_direction_x_noise(double noise_std){
    const double sigma = (noise_std > 0.0) ? noise_std : 0.0;
    for (int i = 0; i < myosin.n; ++i) {
        double dx = 1.0;
        double dy = 0.0;
        double dz = 0.0;
        if (rng != nullptr && sigma > 0.0) {
            dx += gsl_ran_gaussian(rng, sigma);
            dy = gsl_ran_gaussian(rng, sigma);
            dz = gsl_ran_gaussian(rng, sigma);
        }
        vec dir{dx, dy, dz};
        double norm = dir.norm();
        if (norm < EPS) {
            dir = {1.0, 0.0, 0.0};
        } else {
            dir = dir / norm;
        }
        myosin.direction[i] = dir;
    }
    myosin.update_endpoints();
}

void Sarcomere::cb_off_angle(){
    if (actin.n < 2 || myosin.n < 4) {
        return;
    }

    constexpr double ANG2_DEG     = 150.0;
    constexpr double Y_ANCH_LEFT  = -0.015;
    constexpr double Y_ANCH_RIGHT =  0.015;

    const double angle_rad = ANG2_DEG * M_PI / 180.0;
    vec dir0{1.0, 0.0, 0.0};
    vec dir1{std::cos(angle_rad), std::sin(angle_rad), 0.0};
    dir1.normalize();
    printf("dir1: (%f, %f, %f)\n", dir1.x, dir1.y, dir1.z);

    actin.direction[0] = dir0;
    actin.direction[1] = dir1;

    const double half_len = 0.5 * actin.length;
    vec anchor0{0.0, Y_ANCH_LEFT, 0.0};
    vec anchor1{0.0, Y_ANCH_RIGHT, 0.0};

    vec center0 = anchor0 + dir0 * half_len;
    vec center1 = anchor1 + dir1 * half_len;
    printf("Actin 1 center: (%f, %f, %f)\n", center1.x, center1.y, center1.z);
    actin.center[0] = center0;
    actin.center[1] = center1;
    actin.update_endpoints();

    double y0 = actin.right_end_y[0];
    double y1 = actin.right_end_y[1];
    printf("Actin 0 right endpoint y: %f\n", y0);
    printf("Actin 1 right endpoint y: %f\n", y1);
    printf("Actin 0 left endpoint y: %f\n", actin.left_end_y[0]);
    printf("Actin 1 left endpoint y: %f\n", actin.left_end_y[1]);
    std::vector<vec> myosin_positions = {
        { 1.33, y0 - 0.03, 0.0 },
        { 1.33, y0 + 0.03, 0.0 },
        {-1.33, y1 - 0.03, 0.0 },
        {-1.33, y1 + 0.03, 0.0 }
    };


    for (size_t i = 0; i < myosin_positions.size() && i < static_cast<size_t>(myosin.n); ++i) {
        myosin.center[static_cast<int>(i)] = myosin_positions[i];
        myosin.direction[static_cast<int>(i)] = vec{1.0, 0.0, 0.0};
    }
    myosin.update_endpoints();
}

void Sarcomere::am_off_angle(){
    std::vector<vector> actin_positions = {
       {-0.5, 0.015, 0}
    };
    for (int i = 0; i < actin_positions.size(); i++){
        actin.center[i].x = actin_positions[i][0];
        actin.center[i].y = actin_positions[i][1];
        actin.center[i].z = actin_positions[i][2];
    }

    actin.direction[0] = {-0.9397, 0.3420, 0};

    std::vector<vector> myosin_positions = {
        {-1.33, 0, 0}
    };
    for (int i = 0; i < myosin_positions.size(); i++){
        myosin.center[i].x = myosin_positions[i][0];
        myosin.center[i].y = myosin_positions[i][1];
        myosin.center[i].z = myosin_positions[i][2];
        myosin.direction[i] = {1, 0, 0};
    }
}

void Sarcomere::sarcomeric_structure_tight(){
    box[0] = 5.32;
    box[1] = 5.32;
    box[2] = 5.32;

    std::vector<vector> myosin_positions = {
        {-1, -0.03, 0.0}, {1, -0.06, 0.0},
        {-1,  0.03, 0.0}, {1,  0, 0.0},
        {-1,  0.09, 0.0}, {1,  0.06, 0.0}
    };
    for (int i = 0; i < myosin_positions.size(); i++){
        myosin.center[i].x   = myosin_positions[i][0];
        myosin.center[i].y   = myosin_positions[i][1];
        myosin.center[i].z   = myosin_positions[i][2];
        myosin.direction[i] = {1, 0, 0};
    }
    myosin.update_endpoints();

    std::vector<vector> actin_positions = {
        {-2.16, -0.06, 0.0}, 
        {-2.16, 0.0, 0.0}, 
        {-2.16,  0.06, 0.0}, 
        {0.1, -0.03, 0.0}, 
        {0.1, 0.03, 0.0},
        {0.1,  0.09, 0.0}
    };
    for (int i = 0; i < actin_positions.size(); i++){
        printf("Setting actin %d position to (%f, %f, %f)\n", i, actin_positions[i][0], actin_positions[i][1], actin_positions[i][2]);
        actin.center[i].x   = actin_positions[i][0];
        actin.center[i].y   = actin_positions[i][1];
        actin.center[i].z   = actin_positions[i][2];
        actin.direction[i] = {1, 0, 0};
    }
    int n = actin_positions.size();
    actin_positions = {
        {-0.1, -0.06, 0.0},
        {-0.1,  0.00, 0.0},
        {-0.1,  0.06, 0.0},
        { 2.16, -0.03, 0.0},
        { 2.16,  0.03, 0.0},
        { 2.16,  0.09, 0.0}
    };
    for (int i = 0; i < actin_positions.size(); i++){
        printf("Setting actin %d position to (%f, %f, %f)\n", i+n, actin_positions[i][0], actin_positions[i][1], actin_positions[i][2]);
        actin.center[i+n].x   = actin_positions[i][0];
        actin.center[i+n].y   = actin_positions[i][1];
        actin.center[i+n].z   = actin_positions[i][2];
        actin.direction[i+n] = {-1, 0, 0};
    }
    actin.update_endpoints();

    update_system();
}

void Sarcomere::sarcomeric_structure(){
    box[0] = 5.32;
    box[1] = 5.32;
    box[2] = 5.32;

    std::vector<vector> myosin_positions = {
        {-1.33, -0.03, 0.0}, {1.33, -0.06, 0.0},
        {-1.33,  0.03, 0.0}, {1.33,  0, 0.0},
        {-1.33,  0.09, 0.0}, {1.33,  0.06, 0.0}
    };
    for (int i = 0; i < myosin_positions.size(); i++){
        myosin.center[i].x   = myosin_positions[i][0];
        myosin.center[i].y   = myosin_positions[i][1];
        myosin.center[i].z   = myosin_positions[i][2];
        myosin.direction[i] = {1, 0, 0};
    }
    myosin.update_endpoints();

    std::vector<vector> actin_positions = {
        {-2.16, -0.06, 0.0}, 
        {-2.16, 0.0, 0.0}, 
        {-2.16,  0.06, 0.0}, 
        {0.5, -0.03, 0.0}, 
        {0.5, 0.03, 0.0},
        {0.5,  0.09, 0.0}
    };
    for (int i = 0; i < actin_positions.size(); i++){
        printf("Setting actin %d position to (%f, %f, %f)\n", i, actin_positions[i][0], actin_positions[i][1], actin_positions[i][2]);
        actin.center[i].x   = actin_positions[i][0];
        actin.center[i].y   = actin_positions[i][1];
        actin.center[i].z   = actin_positions[i][2];
        actin.direction[i] = {1, 0, 0};
    }

    int n = actin_positions.size();
    actin_positions = {
        {-0.5, -0.06, 0.0},
        {-0.5,  0.00, 0.0},
        {-0.5,  0.06, 0.0},
        { 2.16, -0.03, 0.0},
        { 2.16,  0.03, 0.0},
        { 2.16,  0.09, 0.0}
    };
    for (int i = 0; i < actin_positions.size(); i++){
        printf("Setting actin %d position to (%f, %f, %f)\n", i+n, actin_positions[i][0], actin_positions[i][1], actin_positions[i][2]);
        actin.center[i+n].x   = actin_positions[i][0];
        actin.center[i+n].y   = actin_positions[i][1];
        actin.center[i+n].z   = actin_positions[i][2];
        actin.direction[i+n] = {-1, 0, 0};
    }
    actin.update_endpoints();

    update_system();
}
#include "initialization.h"
#include <vector>
#include <cstdio>
#include <cmath>

namespace initialization {

void partial_fix(Sarcomere& s, int& n_fixed) {
    std::vector<vector> myosin_positions;
    myosin_positions = {
        {0, -2, 0}, {0, -1, 0}, {0, 0, 0}, {0, 1, 0}, {0, 2, 0},
        {0, -2.5, 0}, {0, -1.5, 0}, {0, -0.5, 0}, {0, 0.5, 0}, {0, 1.5, 0}
    };
    for (int i = 0; i < n_fixed; i++){
        s.myosin.center[i].x = myosin_positions[i][0];
        s.myosin.center[i].y = myosin_positions[i][1];
        s.myosin.center[i].z = myosin_positions[i][2];
    }
    for (int i = 0; i < s.myosin.n; i++){
        s.myosin.direction[i] = {1, 0, 0};
    }
    s.myosin.update_endpoints();
    s.update_system();
}

void cb(Sarcomere& s) {
    std::vector<vector> actin_positions = {
        {0.5, -0.015, 0}, {-0.5, 0.015, 0}
    };
    for (int i = 0; i < actin_positions.size(); i++){
        s.actin.center[i].x = actin_positions[i][0];
        s.actin.center[i].y = actin_positions[i][1];
        s.actin.center[i].z = actin_positions[i][2];
    }
    s.actin.direction[0] = {1, 0, 0};
    s.actin.direction[1] = {-1, 0, 0};

    std::vector<vector> myosin_positions = {
        {-1.33, 0.045, 0}, {-1.33, -0.015, 0}, {1.33, 0.015, 0}, {1.33, -0.045, 0}
    };
    for (int i = 0; i < myosin_positions.size(); i++){
        s.myosin.center[i].x = myosin_positions[i][0];
        s.myosin.center[i].y = myosin_positions[i][1];
        s.myosin.center[i].z = myosin_positions[i][2];
        s.myosin.direction[i] = {1, 0, 0};
    }
}

void set_myosin_direction_x_noise(Sarcomere& s, double noise_std) {
    const double sigma = (noise_std > 0.0) ? noise_std : 0.0;
    for (int i = 0; i < s.myosin.n; ++i) {
        double dx = 1.0;
        double dy = 0.0;
        double dz = 0.0;
        if (s.rng != nullptr && sigma > 0.0) {
            dx += gsl_ran_gaussian(s.rng, sigma);
            dy = gsl_ran_gaussian(s.rng, sigma);
            dz = gsl_ran_gaussian(s.rng, sigma);
        }
        vec dir{dx, dy, dz};
        double norm = dir.norm();
        if (norm < EPS) {
            dir = {1.0, 0.0, 0.0};
        } else {
            dir = dir / norm;
        }
        s.myosin.direction[i] = dir;
    }
    s.myosin.update_endpoints();
}

void cb_off_angle(Sarcomere& s) {
    if (s.actin.n < 2 || s.myosin.n < 4) {
        return;
    }
    constexpr double ANG2_DEG     = 150.0;
    constexpr double Y_ANCH_LEFT  = -0.015;
    constexpr double Y_ANCH_RIGHT =  0.015;

    const double angle_rad = ANG2_DEG * M_PI / 180.0;
    vec dir0{1.0, 0.0, 0.0};
    vec dir1{std::cos(angle_rad), std::sin(angle_rad), 0.0};
    dir1.normalize();
    printf("dir1: (%f, %f, %f)\n", dir1.x, dir1.y, dir1.z);

    s.actin.direction[0] = dir0;
    s.actin.direction[1] = dir1;

    const double half_len = 0.5 * s.actin.length;
    vec anchor0{0.0, Y_ANCH_LEFT, 0.0};
    vec anchor1{0.0, Y_ANCH_RIGHT, 0.0};

    vec center0 = anchor0 + dir0 * half_len;
    vec center1 = anchor1 + dir1 * half_len;
    printf("Actin 1 center: (%f, %f, %f)\n", center1.x, center1.y, center1.z);
    s.actin.center[0] = center0;
    s.actin.center[1] = center1;
    s.actin.update_endpoints();

    double y0 = s.actin.right_end_y[0];
    double y1 = s.actin.right_end_y[1];
    printf("Actin 0 right endpoint y: %f\n", y0);
    printf("Actin 1 right endpoint y: %f\n", y1);
    printf("Actin 0 left endpoint y: %f\n", s.actin.left_end_y[0]);
    printf("Actin 1 left endpoint y: %f\n", s.actin.left_end_y[1]);
    std::vector<vec> myosin_positions = {
        { 1.33, y0 - 0.03, 0.0 },
        { 1.33, y0 + 0.03, 0.0 },
        {-1.33, y1 - 0.03, 0.0 },
        {-1.33, y1 + 0.03, 0.0 }
    };

    for (size_t i = 0; i < myosin_positions.size() && i < static_cast<size_t>(s.myosin.n); ++i) {
        s.myosin.center[static_cast<int>(i)] = myosin_positions[i];
        s.myosin.direction[static_cast<int>(i)] = vec{1.0, 0.0, 0.0};
    }
    s.myosin.update_endpoints();
}

void am_off_angle(Sarcomere& s) {
    std::vector<vector> actin_positions = {
       {-0.5, 0.015, 0}
    };
    for (int i = 0; i < actin_positions.size(); i++){
        s.actin.center[i].x = actin_positions[i][0];
        s.actin.center[i].y = actin_positions[i][1];
        s.actin.center[i].z = actin_positions[i][2];
    }
    s.actin.direction[0] = {-0.9397, 0.3420, 0};

    std::vector<vector> myosin_positions = {
        {-1.33, 0, 0}
    };
    for (int i = 0; i < myosin_positions.size(); i++){
        s.myosin.center[i].x = myosin_positions[i][0];
        s.myosin.center[i].y = myosin_positions[i][1];
        s.myosin.center[i].z = myosin_positions[i][2];
        s.myosin.direction[i] = {1, 0, 0};
    }
}

void sarcomeric_structure_tight(Sarcomere& s) {
    s.box[0] = 5.32;
    s.box[1] = 5.32;
    s.box[2] = 5.32;

    std::vector<vector> myosin_positions = {
        {-1, -0.03, 0.0}, {1, -0.06, 0.0},
        {-1,  0.03, 0.0}, {1,  0, 0.0},
        {-1,  0.09, 0.0}, {1,  0.06, 0.0}
    };
    for (int i = 0; i < myosin_positions.size(); i++){
        s.myosin.center[i].x   = myosin_positions[i][0];
        s.myosin.center[i].y   = myosin_positions[i][1];
        s.myosin.center[i].z   = myosin_positions[i][2];
        s.myosin.direction[i] = {1, 0, 0};
    }
    s.myosin.update_endpoints();

    std::vector<vector> actin_positions = {
        {-2.16, -0.06, 0.0}, 
        {-2.16, 0.0, 0.0}, 
        {-2.16,  0.06, 0.0}, 
        {0.1, -0.03, 0.0}, 
        {0.1, 0.03, 0.0},
        {0.1,  0.09, 0.0}
    };
    for (int i = 0; i < actin_positions.size(); i++){
        printf("Setting actin %d position to (%f, %f, %f)\n", i, actin_positions[i][0], actin_positions[i][1], actin_positions[i][2]);
        s.actin.center[i].x   = actin_positions[i][0];
        s.actin.center[i].y   = actin_positions[i][1];
        s.actin.center[i].z   = actin_positions[i][2];
        s.actin.direction[i] = {1, 0, 0};
    }
    int n = actin_positions.size();
    actin_positions = {
        {-0.1, -0.06, 0.0},
        {-0.1,  0.00, 0.0},
        {-0.1,  0.06, 0.0},
        { 2.16, -0.03, 0.0},
        { 2.16,  0.03, 0.0},
        { 2.16,  0.09, 0.0}
    };
    for (int i = 0; i < actin_positions.size(); i++){
        printf("Setting actin %d position to (%f, %f, %f)\n", i+n, actin_positions[i][0], actin_positions[i][1], actin_positions[i][2]);
        s.actin.center[i+n].x   = actin_positions[i][0];
        s.actin.center[i+n].y   = actin_positions[i][1];
        s.actin.center[i+n].z   = actin_positions[i][2];
        s.actin.direction[i+n] = {-1, 0, 0};
    }
    s.actin.update_endpoints();
    s.update_system();
}

void sarcomeric_structure(Sarcomere& s) {
    s.box[0] = 5.32;
    s.box[1] = 5.32;
    s.box[2] = 5.32;

    std::vector<vector> myosin_positions = {
        {-1.33, -0.03, 0.0}, {1.33, -0.06, 0.0},
        {-1.33,  0.03, 0.0}, {1.33,  0, 0.0},
        {-1.33,  0.09, 0.0}, {1.33,  0.06, 0.0}
    };
    for (int i = 0; i < myosin_positions.size(); i++){
        s.myosin.center[i].x   = myosin_positions[i][0];
        s.myosin.center[i].y   = myosin_positions[i][1];
        s.myosin.center[i].z   = myosin_positions[i][2];
        s.myosin.direction[i] = {1, 0, 0};
    }
    s.myosin.update_endpoints();

    std::vector<vector> actin_positions = {
        {-2.16, -0.06, 0.0}, 
        {-2.16, 0.0, 0.0}, 
        {-2.16,  0.06, 0.0}, 
        {0.5, -0.03, 0.0}, 
        {0.5, 0.03, 0.0},
        {0.5,  0.09, 0.0}
    };
    for (int i = 0; i < actin_positions.size(); i++){
        printf("Setting actin %d position to (%f, %f, %f)\n", i, actin_positions[i][0], actin_positions[i][1], actin_positions[i][2]);
        s.actin.center[i].x   = actin_positions[i][0];
        s.actin.center[i].y   = actin_positions[i][1];
        s.actin.center[i].z   = actin_positions[i][2];
        s.actin.direction[i] = {1, 0, 0};
    }

    int n = actin_positions.size();
    actin_positions = {
        {-0.5, -0.06, 0.0}, 
        {-0.5,  0.00, 0.0}, 
        {-0.5,  0.06, 0.0}, 
        { 2.16, -0.03, 0.0}, 
        { 2.16,  0.03, 0.0}, 
        { 2.16,  0.09, 0.0}
    };
    for (int i = 0; i < actin_positions.size(); i++){
        printf("Setting actin %d position to (%f, %f, %f)\n", i+n, actin_positions[i][0], actin_positions[i][1], actin_positions[i][2]);
        s.actin.center[i+n].x   = actin_positions[i][0];
        s.actin.center[i+n].y   = actin_positions[i][1];
        s.actin.center[i+n].z   = actin_positions[i][2];
        s.actin.direction[i+n] = {-1, 0, 0};
    }
    s.actin.update_endpoints();
    s.update_system();
}

} // namespace initialization
