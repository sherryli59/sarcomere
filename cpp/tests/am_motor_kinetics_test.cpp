#include "am_motor_kinetics.h"

#include <cassert>
#include <cmath>

int main() {
    using namespace am_motor_kinetics;
    const double q_eq = thin_activation_calcium_equilibrium(
        1.0, 20.0, 0.5, 3.0, 20.0);
    assert(std::abs(q_eq - 8.0 / 17.0) < 1e-12);

    assert(std::abs(weak_attachment_gate(0.0, 0.02) - 0.02) < 1e-12);
    assert(std::abs(weak_attachment_gate(1.0, 0.02) - 1.0) < 1e-12);
    assert(weak_attachment_gate(0.0, 0.02) > 0.0);
    assert(std::abs(passive_crossbridge_force(100.0, 0.01) - 1.0) < 1e-12);
    assert(std::abs(passive_crossbridge_force(100.0, -0.01) + 1.0) < 1e-12);

    assert(transition_probability(0.0, 1e-3) == 0.0);
    assert(std::abs(transition_probability(100.0, 1e-3) -
                    (1.0 - std::exp(-0.1))) < 1e-14);
    const auto initial_on_count = [](double fraction) {
        int count = 0;
        for (int head = 0; head < 9; ++head) {
            count += initially_on_head(head, 9, fraction) ? 1 : 0;
        }
        return count;
    };
    assert(initial_on_count(0.1) == 1);
    assert(initial_on_count(0.3) == 3 && initial_on_count(0.5) == 5);

    assert(activation(0.2, 0.0, 0.02, 0.2, 0.1) == 1.0);
    const double a0 = activation(0.0, 1.0, 0.02, 0.2, 0.1);
    const double amid = activation(0.05, 1.0, 0.02, 0.2, 0.1);
    assert(std::abs(a0 - 0.1) < 1e-12);
    assert(amid > a0 && amid <= 1.0);
    assert(activation(0.030, 0.050, 0.002, 0.008, 0.0, 0.5) == 0.0);
    assert(activation(0.030, 0.050, 0.002, 0.008, 0.1, 0.5) == 0.1);
    assert(activation(0.004, 0.050, 0.002, 0.008, 0.0, 0.5) > 0.9);
    assert(activation(0.054, 0.050, 0.002, 0.008, 0.0, 0.5) > 0.9);

    assert(std::abs(hill_gate(1.0, 1.0, 2.0) - 0.5) < 1e-12);
    const double q_calcium = thin_activation_step(
        0.0, 1.0, 0.0, 0.01,
        100.0, 0.5, 2.0, 0.0, 0.2, 2.0, 10.0);
    const double q_calcium_strong = thin_activation_step(
        0.0, 1.0, 1.0, 0.01,
        100.0, 0.5, 2.0, 100.0, 0.2, 2.0, 10.0);
    assert(q_calcium > 0.0 && q_calcium < 1.0);
    assert(q_calcium_strong > q_calcium);
    const double q_no_calcium_strong = thin_activation_step(
        1.0, 0.0, 1.0, 0.1,
        100.0, 0.5, 2.0, 100.0, 0.2, 2.0, 10.0);
    assert(q_no_calcium_strong < 1.0);
    double q_decay = 1.0;
    for (int step = 0; step < 1000; ++step) {
        q_decay = thin_activation_step(
            q_decay, 0.0, 0.0, 0.001,
            100.0, 0.5, 2.0, 100.0, 0.2, 2.0, 10.0);
    }
    assert(q_decay < 1e-4);

    const double unloaded = load_slowed_rate(80.0, 0.0, 2.0);
    const double loaded = load_slowed_rate(80.0, 2.0, 2.0);
    assert(std::abs(unloaded - 80.0) < 1e-12);
    assert(loaded < unloaded);
    assert(std::abs(loaded - 80.0 / std::exp(1.0)) < 1e-12);

    assert(thick_recruitment_rate(2.0, 0.0, 1.0) == 2.0);
    assert(thick_recruitment_rate(2.0, 1.0, 1.0) > 2.0);
    assert(thick_folding_rate(10.0, 1.0, 1.0) < 10.0);

    assert(active_crossbridge_force(10.0, -0.2, 3.0) == 0.0);
    assert(active_crossbridge_force(10.0, 0.2, 3.0) == 2.0);
    assert(active_crossbridge_force(10.0, 0.5, 3.0) == 3.0);
    assert(active_crossbridge_force(10.0, 0.5, 0.0) == 5.0);
    assert(working_stroke_velocity(5.0, 0.0, 4.0) == 5.0);
    assert(working_stroke_velocity(5.0, 2.0, 4.0) == 2.5);
    assert(working_stroke_velocity(5.0, 4.0, 4.0) == 0.0);
    assert(working_stroke_velocity(5.0, 8.0, 4.0) == 0.0);
    assert(working_stroke_velocity(5.0, 2.0, 0.0) == 5.0);

    assert(std::abs(strain_accelerated_rate(10.0, 5.0, 0.1) -
                    10.0 * std::exp(0.5)) < 1e-12);
    assert(strain_accelerated_rate(10.0, 5.0, 0.0) == 10.0);
    assert(strain_accelerated_rate(10.0, 0.0, 0.1) == 10.0);
    assert(std::abs(strain_accelerated_rate(10.0, 1e6, 1.0, 2.0) -
                    10.0 * std::exp(2.0)) < 1e-12);
    const double split_axial = anisotropic_strain_accelerated_rate(
        10.0, 5.0, 2.0, 0.1, 7.0, 0.0);
    assert(std::abs(split_axial - 10.0 * std::exp(1.0)) < 1e-12);
    const double split_both = anisotropic_strain_accelerated_rate(
        10.0, 5.0, 2.0, 0.1, 3.0, 0.2);
    assert(std::abs(split_both - 10.0 * std::exp(4.0)) < 1e-12);
    assert(anisotropic_strain_accelerated_rate(
               10.0, 5.0, 2.0, 0.0, 7.0, 0.0) == 10.0);
    assert(std::abs(anisotropic_strain_accelerated_rate(
                        10.0, 5.0, 2.0, 1.0, 7.0, 1.0, 2.0) -
                    10.0 * std::exp(2.0)) < 1e-12);
    assert(occupancy_scaled_per_head_stiffness(600.0, 0, 6) == 0.0);
    assert(occupancy_scaled_per_head_stiffness(600.0, 1, 6) == 100.0);
    assert(occupancy_scaled_per_head_stiffness(600.0, 3, 6) == 100.0);
    assert(occupancy_scaled_per_head_stiffness(600.0, 6, 6) == 100.0);
    assert(occupancy_scaled_per_head_stiffness(600.0, 12, 6) == 50.0);
    assert(overlap_limited_neighbor_capacity(141, 24, 0.0, 0.675) == 0);
    assert(overlap_limited_neighbor_capacity(141, 24, 0.01, 0.675) == 1);
    assert(overlap_limited_neighbor_capacity(141, 24, 0.3375, 0.675) == 12);
    assert(overlap_limited_neighbor_capacity(141, 24, 0.675, 0.675) == 24);
    assert(overlap_limited_neighbor_capacity(12, 24, 0.675, 0.675) == 12);
    assert(overlap_limited_neighbor_capacity(12, 0, 0.3375, 0.675) == 6);
    assert(weak_slots_after_strong(12, 0) == 12);
    assert(weak_slots_after_strong(12, 7) == 5);
    assert(weak_slots_after_strong(12, 12) == 0);
    assert(weak_slots_after_strong(12, 20) == 0);
    assert(std::abs(load_slowed_rate_with_floor(80.0, 4.5, 1.0, 0.2) -
                    80.0 * (0.2 + 0.8 * std::exp(-4.5))) < 1e-12);
    assert(load_slowed_rate_with_floor(80.0, 1e6, 1.0, 0.2) == 16.0);
    assert(load_slowed_rate_with_floor(80.0, 1e6, 1.0, 0.8) == 64.0);

    double q_cycle = 0.0;
    double q_diastolic_max = 0.0;
    for (int step = 1; step <= 200; ++step) {
        const double time = step * 0.001;
        const double calcium = activation(
            time, 0.050, 0.002, 0.008, 0.0, 0.5);
        q_cycle = thin_activation_step(
            q_cycle, calcium, 0.0, 0.001,
            100.0, 0.5, 3.0, 100.0, 0.02, 2.0, 100.0);
        if (std::fmod(time, 0.050) > 0.049 - 1e-12) {
            q_diastolic_max = std::max(q_diastolic_max, q_cycle);
        }
    }
    assert(q_diastolic_max < 0.05);
    return 0;
}
