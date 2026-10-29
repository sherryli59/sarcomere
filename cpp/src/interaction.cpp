#include "interaction.h"
#include "geometry.h"
#include "rod_mobility.h"
#include <cstdio>
#include <cstdlib>  // For exit()
#include <cmath>
#include <deque>
#include <limits>
#include <algorithm>
#include <utility>
#include <vector>

//---------------------------------------------------------------------
// Non-Template Function Definitions
//---------------------------------------------------------------------




real am_energy1(const ArrayXreal& center1, const double& length1, const ArrayXreal& dir1,
    const ArrayXreal& center2, const double& length2, const ArrayXreal& dir2,
    const std::vector<double>& box, const std::array<bool,3>& periodic, const double k_am, const double kappa_am,
    const double cutoff, const double optimal)
{
    // Compute endpoints for both filaments
    ArrayXreal a = center1 - 0.5 * length1 * dir1;
    ArrayXreal b = center1 + 0.5 * length1 * dir1;
    ArrayXreal c = center2 - 0.5 * length2 * dir2;
    ArrayXreal d = center2 + 0.5 * length2 * dir2;

    // Convert to raw arrays for segment-segment distance
    real a_arr[3] = {a[0], a[1], a[2]};
    real b_arr[3] = {b[0], b[1], b[2]};
    real c_arr[3] = {c[0], c[1], c[2]};
    real d_arr[3] = {d[0], d[1], d[2]};
    real dist = segment_segment_distance(a_arr, b_arr, c_arr, d_arr, box, periodic);
    // Angle and strength
    real dot_val = dir1[0]*dir2[0] + dir1[1]*dir2[1] + dir1[2]*dir2[2];
    real strength = abs(dot_val);
    dot_val = std::max(real(-1), std::min(real(1), dot_val));  // clamp
    real angle_energy = 0.5 * kappa_am * (1.0 - dot_val * dot_val); 
    if (dist > cutoff) {
        return angle_energy;
    }
    else {
        const real extension = dist - optimal;
        real energy = 0.5 * k_am * strength * extension * extension + angle_energy;
        // real repulsion_threshold = real(0.9) * optimal;
        // if (dist < repulsion_threshold) {
        //     real rep_offset = (repulsion_threshold - dist)/optimal;
        //     energy += 0.5 * 5 * k_am * rep_offset * rep_offset;
        // }
        return energy;
    }
}

real am_energy(const ArrayXreal& dir1, const ArrayXreal& dir2, const double kappa_am)
{
    real dot_val = dir1[0]*dir2[0] + dir1[1]*dir2[1] + dir1[2]*dir2[2];
    dot_val = std::max(real(-1), std::min(real(1), dot_val));  
    real strength = abs(dot_val);
    real angle_energy = 0.5 * kappa_am * (1.0 - dot_val * dot_val); 
    return angle_energy;
}


real aa_energy(const ArrayXreal& center1, const double& length1, 
    const ArrayXreal& dir1,
    const ArrayXreal& center2, const double& length2, 
    const ArrayXreal& dir2,
    const std::vector<double>& box, const std::array<bool,3>& periodic, 
    const double k_aa, const double kappa_aa,
    const double cutoff, const double optimal,
    const double antiparallel_overlap)
{
    // Compute endpoints for filament 1 in 3D.
    ArrayXreal a = center1 - 0.5 * length1 * dir1;
    ArrayXreal b = center1 + 0.5 * length1 * dir1;

    // Compute endpoints for filament 2 in 3D.
    ArrayXreal c = center2 - 0.5 * length2 * dir2;
    ArrayXreal d = center2 + 0.5 * length2 * dir2;

    // Convert to raw arrays for segment_segment_distance
    real a_arr[3] = {a[0], a[1], a[2]};
    real b_arr[3] = {b[0], b[1], b[2]};
    real c_arr[3] = {c[0], c[1], c[2]};
    real d_arr[3] = {d[0], d[1], d[2]};

    real raw_dist = segment_segment_distance(a_arr, b_arr, c_arr, d_arr, box, periodic);
    real dot_val = dir1[0]*dir2[0] + dir1[1]*dir2[1] + dir1[2]*dir2[2];
    dot_val = std::max(real(-1), std::min(real(1), dot_val));  // clamp for safety
    real angle_energy = 0.5 * kappa_aa * (1.0 - dot_val * dot_val);
    // `cutoff` is an association/capture criterion, not a rupture rule.  Once
    // an AA bond exists its mechanical spring remains active until KMC breaks
    // it; otherwise a strained bond silently loses all restoring force.
    const real extension = raw_dist - optimal;

    // Preserve the working coarse-grained geometry: ACTN2 registers the
    // barbed-end overlap zone (`left_end`), not filament centers.  Parallel
    // neighbors register their barbed ends, whereas antiparallel neighbors
    // interdigitate by the prescribed Z-disc overlap.  The polarity weight is
    // smooth, and using both filament projections makes an AA edge invariant
    // to pair ordering.
    real dx = c_arr[0] - a_arr[0];
    real dy = c_arr[1] - a_arr[1];
    real dz = c_arr[2] - a_arr[2];
    real box_arr[3] = {box.size() > 0 ? box[0] : 0.0,
                       box.size() > 1 ? box[1] : 0.0,
                       box.size() > 2 ? box[2] : 0.0};
    if (periodic[0] && box_arr[0] != real(0.0)) {
        dx -= box_arr[0] * smooth_round(dx / box_arr[0]);
    }
    if (periodic[1] && box_arr[1] != real(0.0)) {
        dy -= box_arr[1] * smooth_round(dy / box_arr[1]);
    }
    if (periodic[2] && box_arr[2] != real(0.0)) {
        dz -= box_arr[2] * smooth_round(dz / box_arr[2]);
    }
    real ds1 = dx*dir1[0] + dy*dir1[1] + dz*dir1[2];
    real ds2 = dx*dir2[0] + dy*dir2[1] + dz*dir2[2];
    const real polarity_weight = 0.5 * (1.0 - dot_val);
    const real target_overlap = std::max(0.0, antiparallel_overlap) * polarity_weight;
    const real axial_error1 = ds1 - target_overlap;
    const real axial_error2 = ds2 + target_overlap;
    real U_ax = 0.25 * k_aa *
        (axial_error1 * axial_error1 + axial_error2 * axial_error2);
    real energy = 0.5 * k_aa * extension * extension + angle_energy + U_ax;
    return energy;
}


std::vector<double> compute_aa_force_and_energy(Filament& actin,
                                                int& actin1_index, int& actin2_index,
                                                const std::vector<double>& box,
                                                const double k_aa, const double kappa_aa,
                                                const double cutoff, const double optimal,
                                                const double antiparallel_overlap)
{
    // Canonicalize an unordered AA edge before differentiating.  Besides
    // avoiding ambiguous closest-point choices for exactly parallel,
    // overlapping rods, this guarantees that asking for (i,j) or (j,i)
    // describes the same physical edge.
    const bool reversed = actin1_index > actin2_index;
    const int canonical_first = reversed ? actin2_index : actin1_index;
    const int canonical_second = reversed ? actin1_index : actin2_index;

    // Construct 3D centers.
    ArrayXreal center1(3);
    center1 << actin.center[canonical_first].x, actin.center[canonical_first].y, actin.center[canonical_first].z;
    ArrayXreal dir1(3);
    dir1 << actin.direction[canonical_first].x, actin.direction[canonical_first].y, actin.direction[canonical_first].z;

    ArrayXreal center2(3);
    center2 << actin.center[canonical_second].x, actin.center[canonical_second].y, actin.center[canonical_second].z;
    ArrayXreal dir2(3);
    dir2 << actin.direction[canonical_second].x, actin.direction[canonical_second].y, actin.direction[canonical_second].z;

    real energy;
    const VectorXd gradient_result = -gradient(
        aa_energy, wrt(center1, dir1, dir2),
        at(center1, actin.length, dir1,
           center2, actin.length, dir2,
           box, actin.periodic_axes, k_aa, kappa_aa, cutoff, optimal,
           antiparallel_overlap), energy);
    std::vector<double> forces(9, 0.0);
    for (int component = 0; component < 3; ++component) {
        forces[component] = (reversed ? -1.0 : 1.0) *
                            gradient_result[component];
        forces[3 + component] = gradient_result[
            (reversed ? 6 : 3) + component];
        forces[6 + component] = gradient_result[
            (reversed ? 3 : 6) + component];
    }
    return forces;

}


namespace {
// log(cosh(x)) without overflow: |x| + log1p(exp(-2|x|)) - log 2.
inline double stable_log_cosh(double x) {
    const double a = std::abs(x);
    return a + std::log1p(std::exp(-2.0 * a)) - 0.6931471805599453;
}
// Saturating radial law shared by the AA link and the AM tether. Linear with
// slope `k` for small extension, force bounded by `f_max`. Setting the
// saturation scale to f_max/k keeps the small-displacement stiffness exactly
// `k`, so a bond still holds as tightly as before and only its behaviour under
// a large pull changes. This is a genuine potential, unlike a displacement
// clamp, so detailed balance and the equilibrium distribution are preserved.
inline vec mul_vec(const vec& v, double s) {
    return vec{v.x * s, v.y * s, v.z * s};
}

inline void saturating_radial(double k, double f_max, double extension,
                              double& force, double& energy) {
    if (!(f_max > 0.0) || !(k > 0.0)) {
        force = k * extension;
        energy = 0.5 * k * extension * extension;
        return;
    }
    const double scale = f_max / k;
    force = f_max * std::tanh(extension / scale);
    energy = f_max * scale * stable_log_cosh(extension / scale);
}
}  // namespace

LocalAALinkResult compute_local_aa_link(
    const Filament& actin, int first, int second,
    double s_first, double s_second,
    const std::vector<double>& box,
    double radial_stiffness, double perpendicular_stiffness,
    double alignment_stiffness,
    double preferred_centerline,
    double antiparallel_overlap,
    double axial_stiffness_ratio,
    double axial_saturation_error,
    double axial_force_max,
    double radial_force_max) {
    LocalAALinkResult result;
    const vec& direction_first = actin.direction[first];
    const vec& direction_second = actin.direction[second];
    const vec center_first = actin.center[first];
    const vec center_second = actin.center[second];
    result.connector = center_second + s_second * direction_second -
                       center_first - s_first * direction_first;
    result.connector.pbc_wrap(box, actin.periodic_axes);
    result.length = result.connector.norm();

    const double radial_k = std::max(0.0, radial_stiffness);
    const double perpendicular_k = std::max(0.0, perpendicular_stiffness);
    const double alignment_k = std::max(0.0, alignment_stiffness);
    const double rest_length = std::max(0.0, preferred_centerline);
    const double projection_first = result.connector.dot(direction_first);
    const double projection_second = result.connector.dot(direction_second);
    const double direction_dot = std::max(
        -1.0, std::min(1.0, direction_first.dot(direction_second)));
    result.direction_dot = direction_dot;
    const double extension = result.length - rest_length;

    double radial_force_unused = 0.0;
    double radial_energy_term = 0.0;
    saturating_radial(radial_k, radial_force_max, extension,
                      radial_force_unused, radial_energy_term);
    result.energy = radial_energy_term +
                    0.5 * perpendicular_k *
                        (projection_first * projection_first +
                         projection_second * projection_second) +
                    0.5 * alignment_k *
                        (1.0 - direction_dot * direction_dot);

    // Exact gradient of the same energy previously differentiated with
    // autodiff. The minimum-image shift is constant except exactly at a box
    // boundary, just as smooth_round was in the autodiff implementation.
    vec connector_gradient {0.0, 0.0, 0.0};
    if (radial_k > 0.0 && result.length > 1e-14) {
        double radial_force = 0.0;
        double radial_energy = 0.0;
        saturating_radial(radial_k, radial_force_max, extension,
                          radial_force, radial_energy);
        connector_gradient += (radial_force / result.length) * result.connector;
    }
    connector_gradient += perpendicular_k *
        (projection_first * direction_first +
         projection_second * direction_second);
    // Status-2-only barbed-end registration. The caller disables this term for
    // status 1 by passing axial_stiffness_ratio=0. The same bounded energy is
    // differentiated for translations and orientations, so it cannot create
    // the force-only, unbounded prestress of the previous implementation.
    vec axial_force_first {0.0, 0.0, 0.0};
    vec axial_grad_dir_first {0.0, 0.0, 0.0};
    vec axial_grad_dir_second {0.0, 0.0, 0.0};
    const double antiparallel_ov = std::max(0.0, antiparallel_overlap);
    if (antiparallel_ov > 0.0) {
        const double link_length = actin.length;
        vec end_first = center_first - 0.5 * link_length * direction_first;
        vec end_second = center_second - 0.5 * link_length * direction_second;
        vec axial_connector = end_second - end_first;
        axial_connector.pbc_wrap(box, actin.periodic_axes);
        const double ds_first = axial_connector.dot(direction_first);
        const double ds_second = axial_connector.dot(direction_second);
        const double polarity_weight = 0.5 * (1.0 - direction_dot);
        const double target_overlap = antiparallel_ov * polarity_weight;
        const double axial_error_first = ds_first - target_overlap;
        const double axial_error_second = ds_second + target_overlap;
        result.registry_error = std::sqrt(0.5 *
            (axial_error_first * axial_error_first +
             axial_error_second * axial_error_second));
        const double axial_k = std::max(0.0, axial_stiffness_ratio) * radial_k;
        // Small-error stiffness, shared by both branches so that a registered
        // pair sees identical mechanics either way.
        const double small_error_k = 0.5 * axial_k;
        double grad_first = 0.0;
        double grad_second = 0.0;
        if (axial_force_max > 0.0 && small_error_k > 0.0) {
            // Saturating register: U = F*s*[log cosh(e1/s) + log cosh(e2/s)]
            // with s = F/small_error_k. dU/de -> small_error_k*e as e->0 and
            // -> +/-F as |e| grows, so a pair that starts far out of register
            // keeps being pulled toward the target instead of feeling nothing.
            const double scale = axial_force_max / small_error_k;
            auto log_cosh = [](double x) {
                const double magnitude = std::abs(x);
                // log(cosh x) -> |x| - log 2 to machine precision well before
                // cosh(x) itself would overflow.
                return magnitude > 20.0 ? magnitude - std::log(2.0)
                                        : std::log(std::cosh(x));
            };
            grad_first = axial_force_max * std::tanh(axial_error_first / scale);
            grad_second = axial_force_max * std::tanh(axial_error_second / scale);
            result.energy += axial_force_max * scale *
                (log_cosh(axial_error_first / scale) +
                 log_cosh(axial_error_second / scale));
        } else {
            const double sigma = std::max(1e-12, axial_saturation_error);
            const double attenuation_first = std::exp(
                -0.5 * axial_error_first * axial_error_first / (sigma * sigma));
            const double attenuation_second = std::exp(
                -0.5 * axial_error_second * axial_error_second / (sigma * sigma));
            grad_first = small_error_k * axial_error_first * attenuation_first;
            grad_second = small_error_k * axial_error_second * attenuation_second;
            result.energy += 0.5 * axial_k * sigma * sigma *
                (2.0 - attenuation_first - attenuation_second);
        }
        axial_force_first = grad_first * direction_first +
                            grad_second * direction_second;
        axial_grad_dir_first =
            grad_first *
                (axial_connector + (0.5 * link_length) * direction_first +
                 (0.5 * antiparallel_ov) * direction_second) +
            grad_second *
                ((0.5 * (link_length - antiparallel_ov)) * direction_second);
        axial_grad_dir_second =
            grad_first *
                ((-0.5 * (link_length - antiparallel_ov)) * direction_first) +
            grad_second *
                (axial_connector - (0.5 * link_length) * direction_second -
                 (0.5 * antiparallel_ov) * direction_first);
    }
    result.force_on_first = connector_gradient + axial_force_first;

    const vec gradient_direction_first =
        -s_first * connector_gradient +
        perpendicular_k * projection_first * result.connector -
        alignment_k * direction_dot * direction_second +
        axial_grad_dir_first;
    const vec gradient_direction_second =
        s_second * connector_gradient +
        perpendicular_k * projection_second * result.connector -
        alignment_k * direction_dot * direction_first +
        axial_grad_dir_second;
    result.torque_on_first = -gradient_direction_first;
    result.torque_on_second = -gradient_direction_second;
    // The same four products, kept separately so a run can say WHICH ACTN2
    // stiffness supplies the torque.  Signs match the negation above.
    result.torque_first_lever = s_first * connector_gradient;
    result.torque_first_perp =
        -(perpendicular_k * projection_first) * result.connector;
    result.torque_first_align =
        (alignment_k * direction_dot) * direction_second;
    result.torque_first_axial = -1.0 * axial_grad_dir_first;
    result.torque_second_lever = -s_second * connector_gradient;
    result.torque_second_perp =
        -(perpendicular_k * projection_second) * result.connector;
    result.torque_second_align =
        (alignment_k * direction_dot) * direction_first;
    result.torque_second_axial = -1.0 * axial_grad_dir_second;
    return result;
}

std::vector<double> compute_am_force_and_energy(Filament& actin, Myosin& myosin,
                                                int& actin_index, int& myosin_index,
                                                const std::vector<double>& box,
                                                const double k_am, const double kappa_am,
                                                const double cutoff, const double optimal)
{
    // Define 3D center positions
    ArrayXreal center1(3);
    center1[0] = actin.center[actin_index].x;
    center1[1] = actin.center[actin_index].y;
    center1[2] = actin.center[actin_index].z;
    ArrayXreal dir1(3);
    dir1[0] = actin.direction[actin_index].x;
    dir1[1] = actin.direction[actin_index].y;
    dir1[2] = actin.direction[actin_index].z;
    ArrayXreal center2(3);
    center2[0] = myosin.center[myosin_index].x;
    center2[1] = myosin.center[myosin_index].y;
    center2[2] = myosin.center[myosin_index].z;
    ArrayXreal dir2(3);
    dir2[0] = myosin.direction[myosin_index].x;
    dir2[1] = myosin.direction[myosin_index].y;
    dir2[2] = myosin.direction[myosin_index].z;
    real u;
    std::vector<double> forces;
    VectorXd forces_3; // 3D version of force vector
    if (k_am > 1e-6) {
        forces_3 = -gradient(am_energy1, wrt(center1, dir1, dir2),
                                at(center1, actin.length, dir1,
                                   center2, myosin.length, dir2, box, actin.periodic_axes,
                                   k_am, kappa_am, cutoff, optimal), u);
        forces.resize(forces_3.size());
        VectorXd::Map(&forces[0], forces_3.size()) = forces_3;
        // The passive AM interaction is radial proximity only.  The former
        // axial "endstop" generated propulsion-like forces even with v_am=0,
        // duplicating motor mechanics and ejecting a correctly overlapped thin
        // filament. Axial force belongs exclusively to explicit bound heads.
    }
    else {
        forces_3 = -gradient(am_energy, wrt(dir1, dir2),
                                at(dir1, dir2, kappa_am), u);
        // Prepend three zeros for the force vector (x, y, z)
        forces.resize(forces_3.size() + 3);
        forces[0] = 0;
        forces[1] = 0;
        forces[2] = 0;
        VectorXd::Map(&forces[3], forces_3.size()) = forces_3;
    }
    return forces;
}

AMTetherResult compute_am_proximity_tether(const Filament& actin,
                                           const Myosin& myosin,
                                           int actin_index,
                                           int myosin_index,
                                           const std::vector<double>& box,
                                           double stiffness,
                                           double optimal,
                                           bool directional) {
    AMTetherResult result;
    vec myosin_start = myosin.left_end[myosin_index];
    vec myosin_end = myosin.right_end[myosin_index];
    if (directional) {
        // Restrict the attachment point to the productive, head-bearing outer
        // third selected by actin polarity.  This prevents the moment arm from
        // drifting into the bare zone or the wrong half of a bipolar filament.
        const vec full_axis = myosin_end - myosin_start;
        if (actin.direction[actin_index].dot(myosin.direction[myosin_index]) > 0.0) {
            myosin_end = myosin_start +
                          geometry::am_head_bearing_fraction * full_axis;
        } else {
            myosin_start = myosin_end -
                            geometry::am_head_bearing_fraction * full_axis;
        }
    }
    const auto closest = geometry::segment_segment_distance_w_normal(
        actin.left_end[actin_index], actin.right_end[actin_index],
        myosin_start, myosin_end,
        box, actin.periodic_axes);
    result.actin_attachment = closest.second.start;
    result.myosin_attachment = closest.second.end;

    vec axis_a = actin.direction[actin_index];
    vec axis_m = myosin.direction[myosin_index];
    if (axis_a.dot(axis_m) < 0.0) axis_m = -axis_m;
    vec common_axis = axis_a + axis_m;
    if (common_axis.norm() <= 1e-12) common_axis = axis_a;
    common_axis.normalize();

    vec transverse = closest.second.normal -
                     common_axis * closest.second.normal.dot(common_axis);
    const double distance = transverse.norm();
    result.transverse_distance = distance;
    if (!(stiffness > 0.0) || !(optimal > 0.0) || distance <= 1e-12) {
        return result;
    }
    const vec unit = transverse / distance;
    const double magnitude = -stiffness * (distance - optimal) /
                             (optimal * optimal);
    result.force_on_actin = magnitude * unit;
    result.force_on_myosin = -result.force_on_actin;
    return result;
}

AMTetherResult compute_am_head_tether(const Filament& actin,
                                      const Myosin& myosin,
                                      int actin_index,
                                      int myosin_index,
                                      const std::vector<double>& box,
                                      double myosin_anchor_axial,
                                      double stiffness,
                                      double optimal,
                                      double myosin_anchor_azimuth,
                                      double myosin_radial_offset,
                                      const vec* material_e1,
                                      const vec* material_e2,
                                      double actin_radius,
                                      double radial_slack,
    double tether_force_max) {
    AMTetherResult result;
    vec axis_a = actin.direction[actin_index];
    vec axis_m = myosin.direction[myosin_index];
    axis_a.normalize();
    axis_m.normalize();

    const auto axis_geometry = geometry::segment_segment_distance_w_normal(
        actin.left_end[actin_index], actin.right_end[actin_index],
        myosin.left_end[myosin_index], myosin.right_end[myosin_index],
        box, actin.periodic_axes);
    result.centerline_distance = axis_geometry.first;
    result.surface_gap = result.centerline_distance -
                         (myosin.radius + std::max(0.0, actin_radius));

    vec transverse_1;
    vec transverse_2;
    if (material_e1 != nullptr && material_e2 != nullptr) {
        transverse_1 = *material_e1 - axis_m * material_e1->dot(axis_m);
        if (transverse_1.norm() <= 1e-12) transverse_1 = *material_e2;
        transverse_1.normalize();
        transverse_2 = axis_m.cross(transverse_1);
        transverse_2.normalize();
        if (transverse_2.dot(*material_e2) < 0.0) transverse_2 = -transverse_2;
    } else {
        vec reference = std::abs(axis_m.z) < 0.9
            ? vec{0.0, 0.0, 1.0} : vec{0.0, 1.0, 0.0};
        transverse_1 = axis_m.cross(reference);
        if (transverse_1.norm() <= 1e-12) transverse_1 = {1.0, 0.0, 0.0};
        transverse_1.normalize();
        transverse_2 = axis_m.cross(transverse_1);
        transverse_2.normalize();
    }
    const vec radial_anchor = myosin_radial_offset *
        (std::cos(myosin_anchor_azimuth) * transverse_1 +
         std::sin(myosin_anchor_azimuth) * transverse_2);
    result.myosin_attachment = static_cast<vec>(myosin.center[myosin_index]) +
                                myosin_anchor_axial * axis_m + radial_anchor;
    vec center_disp = static_cast<vec>(actin.center[actin_index]) -
                      result.myosin_attachment;
    center_disp.pbc_wrap(box, actin.periodic_axes);
    const vec actin_center_image = result.myosin_attachment + center_disp;
    result.unbounded_actin_coordinate =
        (result.myosin_attachment - actin_center_image).dot(axis_a);
    const double actin_half_length = 0.5 * actin.length;
    result.projection_inside_actin =
        result.unbounded_actin_coordinate >= -actin_half_length &&
        result.unbounded_actin_coordinate <= actin_half_length;
    const double actin_coordinate = std::clamp(
        result.unbounded_actin_coordinate,
        -actin_half_length, actin_half_length);
    result.actin_coordinate = actin_coordinate;
    result.actin_attachment = actin_center_image + actin_coordinate * axis_a;

    if (axis_a.dot(axis_m) < 0.0) axis_m = -axis_m;
    vec common_axis = axis_a + axis_m;
    if (common_axis.norm() <= 1e-12) common_axis = axis_a;
    common_axis.normalize();
    vec connector = result.actin_attachment - result.myosin_attachment;
    connector.pbc_wrap(box, actin.periodic_axes);
    result.head_connector_distance = connector.norm();
    result.head_reach_surface_distance =
        result.head_connector_distance - std::max(0.0, actin_radius);
    vec transverse = connector - common_axis * connector.dot(common_axis);
    result.transverse_distance = transverse.norm();
    result.head_surface_distance = result.transverse_distance -
                                   std::max(0.0, actin_radius);
    if (!(stiffness > 0.0) || !(optimal > 0.0) ||
        result.centerline_distance <= 1e-12) return result;

    const vec unit_normal = axis_geometry.second.normal /
                            result.centerline_distance;
    const double extension = result.centerline_distance - optimal;
    const double slack = std::max(0.0, radial_slack);
    double loaded_extension = 0.0;
    if (extension > slack) loaded_extension = extension - slack;
    if (extension < -slack) loaded_extension = extension + slack;
    double tether_force = 0.0;
    double tether_energy = 0.0;
    saturating_radial(stiffness, tether_force_max, loaded_extension,
                      tether_force, tether_energy);
    const double magnitude = -tether_force;
    result.force_on_actin = magnitude * unit_normal;
    result.force_on_myosin = -result.force_on_actin;
    return result;
}

AMTwoSpringHeadResult compute_am_two_spring_head(
    const Filament& actin, const Myosin& myosin,
    int actin_index, int myosin_index, const std::vector<double>& box,
    double actin_coordinate, double myosin_anchor_axial,
    double myosin_anchor_azimuth, double myosin_radial_offset,
    const vec* material_e1, const vec* material_e2,
    double radial_preferred, double radial_stiffness, double radial_slack,
    double axial_stiffness, double stroke_progress, double stall_force,
    bool strong, bool bidirectional_strong,
    double radial_force_max,
    AmRadialMode radial_mode) {
    AMTwoSpringHeadResult result;
    vec axis_a = actin.direction[actin_index];
    vec axis_m = myosin.direction[myosin_index];
    if (axis_a.norm() <= 1e-12 || axis_m.norm() <= 1e-12) return result;
    axis_a.normalize();
    axis_m.normalize();

    vec e1;
    vec e2;
    if (material_e1 != nullptr && material_e2 != nullptr) {
        e1 = *material_e1 - axis_m * material_e1->dot(axis_m);
        if (e1.norm() <= 1e-12) e1 = *material_e2;
        e1.normalize();
        e2 = axis_m.cross(e1);
        e2.normalize();
        if (e2.dot(*material_e2) < 0.0) e2 = -1.0 * e2;
    } else {
        const vec reference = std::abs(axis_m.z) < 0.9
            ? vec{0.0, 0.0, 1.0} : vec{0.0, 1.0, 0.0};
        e1 = axis_m.cross(reference);
        if (e1.norm() <= 1e-12) e1 = {1.0, 0.0, 0.0};
        e1.normalize();
        e2 = axis_m.cross(e1);
        e2.normalize();
    }
    const vec radial_anchor = std::max(0.0, myosin_radial_offset) *
        (std::cos(myosin_anchor_azimuth) * e1 +
         std::sin(myosin_anchor_azimuth) * e2);
    result.myosin_attachment = static_cast<vec>(myosin.center[myosin_index]) +
        myosin_anchor_axial * axis_m + radial_anchor;

    vec center_disp = static_cast<vec>(actin.center[actin_index]) -
                      result.myosin_attachment;
    center_disp.pbc_wrap(box, actin.periodic_axes);
    const vec actin_center_image = result.myosin_attachment + center_disp;
    const double bounded_coordinate = std::clamp(
        actin_coordinate, -0.5 * actin.length, 0.5 * actin.length);
    result.actin_attachment = actin_center_image + bounded_coordinate * axis_a;
    result.connector = result.actin_attachment - result.myosin_attachment;
    result.connector.pbc_wrap(box, actin.periodic_axes);

    const double z = result.connector.dot(axis_a);
    const vec transverse = result.connector - z * axis_a;
    const double rho = transverse.norm();
    result.longitudinal_extension = z;
    result.transverse_distance = rho;
    result.head_reach_surface_distance = result.connector.norm() -
                                         std::max(0.0, radial_preferred);

    vec gradient_h {0.0, 0.0, 0.0};
    vec explicit_gradient_axis {0.0, 0.0, 0.0};

    // What the radial spring penalises, and the direction it pulls.
    // kAmRadialLine reproduces the shipped behaviour exactly.
    double q = rho;
    vec q_hat = (rho > 1e-12) ? mul_vec(transverse, 1.0 / rho)
                              : vec{0.0, 0.0, 0.0};
    bool explicit_axis_term = true;
    if (radial_mode == AmRadialMode::kAmRadialConnector) {
        q = result.connector.norm();
        q_hat = (q > 1e-12) ? mul_vec(result.connector, 1.0 / q)
                            : vec{0.0, 0.0, 0.0};
        explicit_axis_term = false;   // |h| depends on u_a only via s*grad
    } else if (radial_mode == AmRadialMode::kAmRadialAxis) {
        vec from_centre = result.actin_attachment -
                          static_cast<vec>(myosin.center[myosin_index]);
        from_centre.pbc_wrap(box, actin.periodic_axes);
        const vec perp = from_centre - axis_m * from_centre.dot(axis_m);
        q = perp.norm();
        q_hat = (q > 1e-12) ? mul_vec(perp, 1.0 / q) : vec{0.0, 0.0, 0.0};
        explicit_axis_term = false;   // a projection of P; s*grad covers u_a
    }
    const double radial_extension = q - std::max(0.0, radial_preferred);
    const double slack = std::max(0.0, radial_slack);
    double loaded_radial_extension = 0.0;
    if (radial_extension > slack) {
        loaded_radial_extension = radial_extension - slack;
    } else if (radial_extension < -slack) {
        loaded_radial_extension = radial_extension + slack;
    }
    if (radial_stiffness > 0.0 && q > 1e-12 &&
        loaded_radial_extension != 0.0) {
        // Saturating when radial_force_max > 0: f_max*tanh(ext/scale) with
        // scale = f_max/k.  The force is bounded while log cosh grows LINEARLY,
        // so the well still confines and there is no finite-energy escape.
        // radial_force_max <= 0 keeps the unbounded linear spring, which is what
        // the legacy azimuth-ON arm uses.
        double dU_drho = 0.0;
        double radial_energy_term = 0.0;
        saturating_radial(radial_stiffness, radial_force_max,
                          loaded_radial_extension, dU_drho, radial_energy_term);
        if (radial_mode == AmRadialMode::kAmRadialLine) {
            // Written exactly as before the mode switch existed, so the
            // default path is BITWISE identical -- (dU/rho)*t and t*(dU/rho)
            // do not round the same way, and this is the shipped behaviour.
            gradient_h += (dU_drho / rho) * transverse;
            // rho^2 = h.h - (h.u)^2.  Explicit orientation derivative, needed
            // in addition to the material-point lever arm.
            explicit_gradient_axis += -(dU_drho * z / rho) * result.connector;
        } else {
            gradient_h += mul_vec(q_hat, dU_drho);
        }
        result.radial_energy = radial_energy_term;
        result.energy += result.radial_energy;
    }

    double dU_dz = 0.0;
    if (strong) {
        result.productive_strain = stroke_progress - z;
        // The attached head is an elastic element whose rest coordinate is
        // advanced one-way by the power stroke.  Tension-only rectification
        // makes an overshooting thin filament feel exactly zero force, so
        // nothing opposes penetration past the M-line.  In bidirectional mode
        // the head is a signed spring about stroke_progress: it pulls while the
        // actin lags and pushes back once the actin overshoots, matching the
        // weak-head branch below and real cross-bridge elasticity.
        const bool axial_active = axial_stiffness > 0.0 &&
            (bidirectional_strong ? result.productive_strain != 0.0
                                  : result.productive_strain > 0.0);
        if (axial_active) {
            const double uncapped_force =
                axial_stiffness * result.productive_strain;
            const double magnitude = std::abs(uncapped_force);
            result.active_force = stall_force > 0.0
                ? std::copysign(std::min(magnitude, stall_force), uncapped_force)
                : uncapped_force;
            dU_dz = -result.active_force;
            const double strain_magnitude = std::abs(result.productive_strain);
            if (stall_force > 0.0 && magnitude > stall_force) {
                result.axial_energy = stall_force * strain_magnitude -
                    0.5 * stall_force * stall_force / axial_stiffness;
            } else {
                result.axial_energy = 0.5 * axial_stiffness *
                    result.productive_strain * result.productive_strain;
            }
            result.energy += result.axial_energy;
        }
    } else {
        result.productive_strain = -z;
        if (axial_stiffness > 0.0) {
            dU_dz = axial_stiffness * z;
            result.axial_energy = 0.5 * axial_stiffness * z * z;
            result.energy += result.axial_energy;
        }
    }
    gradient_h += dU_dz * axis_a;
    explicit_gradient_axis += dU_dz * result.connector;

    result.force_on_actin = -1.0 * gradient_h;
    result.force_on_myosin = gradient_h;
    result.axial_force_on_actin = -dU_dz * axis_a;
    result.radial_force_on_actin =
        result.force_on_actin - result.axial_force_on_actin;

    const vec actin_axis_gradient =
        bounded_coordinate * gradient_h + explicit_gradient_axis;
    result.physical_torque_on_actin =
        -1.0 * axis_a.cross(actin_axis_gradient);
    vec myosin_lever = result.myosin_attachment -
                       static_cast<vec>(myosin.center[myosin_index]);
    myosin_lever.pbc_wrap(box, actin.periodic_axes);
    result.physical_torque_on_myosin = myosin_lever.cross(gradient_h);
    return result;
}

double am_head_azimuth_offset(
    const Filament& actin, const Myosin& myosin,
    int actin_index, int myosin_index, const std::vector<double>& box,
    double myosin_anchor_axial, double crown_azimuth,
    const vec& material_e1, const vec& material_e2) {
    const double kNaN = std::numeric_limits<double>::quiet_NaN();
    vec axis_m = myosin.direction[myosin_index];
    vec axis_a = actin.direction[actin_index];
    if (axis_m.norm() <= 1e-12 || axis_a.norm() <= 1e-12) return kNaN;
    axis_m.normalize();
    axis_a.normalize();
    vec e1 = material_e1 - axis_m * material_e1.dot(axis_m);
    if (e1.norm() <= 1e-12) return kNaN;
    e1.normalize();
    vec e2 = axis_m.cross(e1);
    e2.normalize();
    if (e2.dot(material_e2) < 0.0) e2 = -1.0 * e2;

    const vec axial_anchor = static_cast<vec>(myosin.center[myosin_index]) +
                             myosin_anchor_axial * axis_m;
    vec center_disp = static_cast<vec>(actin.center[actin_index]) - axial_anchor;
    center_disp.pbc_wrap(box, actin.periodic_axes);
    const vec actin_center_image = axial_anchor + center_disp;
    const double actin_coordinate = std::clamp(
        (axial_anchor - actin_center_image).dot(axis_a),
        -0.5 * actin.length, 0.5 * actin.length);
    const vec closest_actin = actin_center_image + actin_coordinate * axis_a;
    vec radial = closest_actin - axial_anchor;
    radial -= axis_m * radial.dot(axis_m);
    if (radial.norm() <= 1e-12) return kNaN;
    radial.normalize();
    const double target = std::atan2(radial.dot(e2), radial.dot(e1));
    return std::remainder(target - crown_azimuth, 2.0 * M_PI);
}

bool am_head_azimuth_admits(
    const Filament& actin, const Myosin& myosin,
    int actin_index, int myosin_index, const std::vector<double>& box,
    double myosin_anchor_axial, double crown_azimuth, double half_width,
    const vec& material_e1, const vec& material_e2) {
    if (!(half_width > 0.0)) return true;
    if (half_width >= M_PI) return true;
    const double offset = am_head_azimuth_offset(
        actin, myosin, actin_index, myosin_index, box,
        myosin_anchor_axial, crown_azimuth, material_e1, material_e2);
    // Degenerate geometry carries no azimuth information, so it must not be
    // allowed to silently forbid every head.
    if (!std::isfinite(offset)) return true;
    return std::abs(offset) <= half_width;
}

double select_am_binding_azimuth(
    const Filament& actin, const Myosin& myosin,
    int actin_index, int myosin_index, const std::vector<double>& box,
    double myosin_anchor_axial, double crown_azimuth, double max_swing,
    const vec& material_e1, const vec& material_e2) {
    vec axis_m = myosin.direction[myosin_index];
    vec axis_a = actin.direction[actin_index];
    if (axis_m.norm() <= 1e-12 || axis_a.norm() <= 1e-12) {
        return crown_azimuth;
    }
    axis_m.normalize();
    axis_a.normalize();
    vec e1 = material_e1 - axis_m * material_e1.dot(axis_m);
    if (e1.norm() <= 1e-12) return crown_azimuth;
    e1.normalize();
    vec e2 = axis_m.cross(e1);
    e2.normalize();
    if (e2.dot(material_e2) < 0.0) e2 = -1.0 * e2;

    const vec axial_anchor = static_cast<vec>(myosin.center[myosin_index]) +
                             myosin_anchor_axial * axis_m;
    vec center_disp = static_cast<vec>(actin.center[actin_index]) - axial_anchor;
    center_disp.pbc_wrap(box, actin.periodic_axes);
    const vec actin_center_image = axial_anchor + center_disp;
    const double actin_coordinate = std::clamp(
        (axial_anchor - actin_center_image).dot(axis_a),
        -0.5 * actin.length, 0.5 * actin.length);
    const vec closest_actin = actin_center_image + actin_coordinate * axis_a;
    vec radial = closest_actin - axial_anchor;
    radial -= axis_m * radial.dot(axis_m);
    if (radial.norm() <= 1e-12) return crown_azimuth;
    radial.normalize();
    const double target = std::atan2(radial.dot(e2), radial.dot(e1));
    const double delta = std::remainder(target - crown_azimuth, 2.0 * M_PI);
    const double bounded = std::clamp(delta, -std::max(0.0, max_swing),
                                     std::max(0.0, max_swing));
    return crown_azimuth + bounded;
}

namespace {

vec regularized_steric_normal(
    vec first_axis, vec second_axis, vec center_displacement,
    int first_id, int second_id) {
    constexpr double eps = 1e-12;
    first_axis.normalize();
    second_axis.normalize();

    // A true nonparallel intersection has zero distance and therefore no
    // unique signed distance gradient.  The normal to both axes is the local
    // capsule-separation direction; its sign is selected by the resolved
    // center displacement when possible.  Swapping the rods reverses both the
    // cross product and center displacement, preserving pair antisymmetry.
    vec normal = first_axis.cross(second_axis);
    if (normal.norm() > eps) {
        normal.normalize();
        const double projection = normal.dot(center_displacement);
        if (projection < -eps) normal = -normal;
        return normal;
    }

    // Parallel/coincident rods use their transverse center separation.  If
    // that also vanishes, choose a deterministic pair-specific transverse
    // direction so exact coincidence still experiences a reciprocal force.
    if (first_axis.dot(second_axis) < 0.0) second_axis = -second_axis;
    vec common_axis = first_axis + second_axis;
    if (common_axis.norm() <= eps) common_axis = first_axis;
    common_axis.normalize();
    normal = center_displacement -
             common_axis * center_displacement.dot(common_axis);
    if (normal.norm() > eps) {
        normal.normalize();
        return normal;
    }

    const std::uint64_t low =
        static_cast<std::uint64_t>(std::min(first_id, second_id) + 1);
    const std::uint64_t high =
        static_cast<std::uint64_t>(std::max(first_id, second_id) + 1);
    const std::uint64_t hash =
        (low * 0x9e3779b97f4a7c15ULL) ^
        (high * 0xbf58476d1ce4e5b9ULL);
    const std::array<vec, 3> references{
        vec{1.0, 0.0, 0.0}, vec{0.0, 1.0, 0.0}, vec{0.0, 0.0, 1.0}};
    vec reference = references[static_cast<size_t>(hash % 3ULL)];
    if (std::abs(reference.dot(common_axis)) > 0.9) {
        reference = references[static_cast<size_t>((hash + 1ULL) % 3ULL)];
    }
    normal = reference - common_axis * reference.dot(common_axis);
    normal.normalize();
    if (first_id > second_id) normal = -normal;
    return normal;
}

vec steric_contact_orientation_drive(
    const vec& center, vec axis, const vec& contact, const vec& force,
    const std::vector<double>& box, const std::array<bool,3>& periodic,
    const vec& other_axis) {
    axis.normalize();
    vec normalized_other = other_axis;
    normalized_other.normalize();
    // The closest point is non-unique along parallel overlapping rods.  Do
    // not manufacture a torque from whichever endpoint the distance solver
    // happened to choose; their reciprocal transverse force is sufficient.
    if (axis.cross(normalized_other).norm() <= 1e-8) {
        return {0.0, 0.0, 0.0};
    }
    vec lever = contact - center;
    lever.pbc_wrap(box, periodic);
    const vec physical_torque = lever.cross(force);
    return physical_torque.cross(axis);
}

}  // namespace

vec compute_actin_myosin_radial_repulsion(
    const Filament& actin, const Myosin& myosin,
    int actin_index, int myosin_index, const std::vector<double>& box,
    double excluded_centerline_distance, double stiffness) {
    return compute_actin_myosin_radial_repulsion_result(
        actin, myosin, actin_index, myosin_index, box,
        excluded_centerline_distance, stiffness).force_on_first;
}

RepulsionResult compute_actin_myosin_radial_repulsion_result(
    const Filament& actin, const Myosin& myosin,
    int actin_index, int myosin_index, const std::vector<double>& box,
    double excluded_centerline_distance, double stiffness) {
    RepulsionResult result{};
    if (!(excluded_centerline_distance > 0.0) || !(stiffness > 0.0)) {
        return result;
    }
    const auto geometry_result = geometry::segment_segment_distance_w_normal(
        actin.left_end[actin_index], actin.right_end[actin_index],
        myosin.left_end[myosin_index], myosin.right_end[myosin_index],
        box, actin.periodic_axes);
    const double distance = geometry_result.first;
    result.segment_distance = distance;
    if (distance >= excluded_centerline_distance) return result;
    vec normal = geometry_result.second.normal;
    if (normal.norm() <= 1e-12) {
        vec center_displacement =
            static_cast<vec>(actin.center[actin_index]) -
            static_cast<vec>(myosin.center[myosin_index]);
        center_displacement.pbc_wrap(box, actin.periodic_axes);
        normal = regularized_steric_normal(
            actin.direction[actin_index], myosin.direction[myosin_index],
            center_displacement, actin_index, actin.n + myosin_index);
    }
    normal.normalize();
    const vec force =
        stiffness * (excluded_centerline_distance - distance) * normal;
    result.applied = true;
    result.force_on_first = force;
    result.force_on_second = -force;
    result.torque_on_first = steric_contact_orientation_drive(
        actin.center[actin_index], actin.direction[actin_index],
        geometry_result.second.start, force, box, actin.periodic_axes,
        myosin.direction[myosin_index]);
    result.torque_on_second = steric_contact_orientation_drive(
        myosin.center[myosin_index], myosin.direction[myosin_index],
        geometry_result.second.end, -force, box, myosin.periodic_axes,
        actin.direction[actin_index]);
    return result;
}

ACTN2PatchGeometry compute_actn2_patch_geometry(
    const Filament& actin, int actin_i, int actin_j, double s_i, double s_j,
    const std::vector<double>& box) {
    ACTN2PatchGeometry patch{};
    if (actin_i < 0 || actin_j < 0 || actin_i >= actin.n ||
        actin_j >= actin.n) {
        return patch;
    }
    const vec center_i = actin.center[actin_i];
    const vec center_j = actin.center[actin_j];
    patch.first = center_i + s_i * actin.direction[actin_i];
    vec second_raw = center_j + s_j * actin.direction[actin_j];
    vec connector = second_raw - patch.first;
    connector.pbc_wrap(box, actin.periodic_axes);
    patch.second = patch.first + connector;
    patch.midpoint = patch.first + 0.5 * connector;
    patch.half_length = 0.5 * connector.norm();
    return patch;
}

namespace {

double periodic_center_distance_squared(
    const vec& first, const vec& second, const std::vector<double>& box,
    const std::array<bool, 3>& periodic_axes) {
    vec displacement = first - second;
    displacement.pbc_wrap(box, periodic_axes);
    return displacement.norm_squared();
}

}  // namespace

bool actn2_patch_may_contact_myosin_tip(
    const ACTN2PatchGeometry& patch, const Myosin& myosin, int myosin_index,
    const std::vector<double>& box, double axial_margin,
    double transverse_cutoff) {
    if (myosin_index < 0 || myosin_index >= myosin.n || axial_margin < 0.0 ||
        transverse_cutoff < 0.0) {
        return false;
    }
    const double axial_reach = 0.5 * myosin.length + axial_margin;
    const double reach = std::hypot(axial_reach, transverse_cutoff) +
        patch.half_length;
    return periodic_center_distance_squared(
        patch.midpoint, myosin.center[myosin_index], box,
        myosin.periodic_axes) <= reach * reach;
}

bool actn2_patch_may_contact_myosin_body(
    const ACTN2PatchGeometry& patch, const Myosin& myosin, int myosin_index,
    const std::vector<double>& box, double outer_contact_distance) {
    if (myosin_index < 0 || myosin_index >= myosin.n ||
        outer_contact_distance < 0.0) {
        return false;
    }
    const double reach = patch.half_length + 0.5 * myosin.length +
        outer_contact_distance;
    return periodic_center_distance_squared(
        patch.midpoint, myosin.center[myosin_index], box,
        myosin.periodic_axes) <= reach * reach;
}

ACTN2MyosinContactResult compute_actn2_myosin_contact(
    const Filament& actin, const Myosin& myosin,
    int actin_i, int actin_j, double s_i, double s_j, int myosin_index,
    const std::vector<double>& box, double hard_core_distance,
    double epsilon, double force_cap, double shell_width,
    double shell_energy) {
    ACTN2MyosinContactResult result{};
    if (actin_i < 0 || actin_j < 0 || actin_i >= actin.n ||
        actin_j >= actin.n || myosin_index < 0 ||
        myosin_index >= myosin.n || !(hard_core_distance > 0.0) ||
        !std::isfinite(hard_core_distance) || epsilon < 0.0 ||
        !std::isfinite(epsilon) || shell_width < 0.0 ||
        !std::isfinite(shell_width) || shell_energy < 0.0 ||
        !std::isfinite(shell_energy) ||
        (!(epsilon > 0.0) && !(shell_width > 0.0 && shell_energy > 0.0))) {
        return result;
    }

    const vec center_i = actin.center[actin_i];
    const vec center_j = actin.center[actin_j];
    vec first = center_i + s_i * actin.direction[actin_i];
    vec second_raw = center_j + s_j * actin.direction[actin_j];
    vec connector = second_raw - first;
    connector.pbc_wrap(box, actin.periodic_axes);
    vec second = first + connector;
    const auto geometry_result = geometry::segment_segment_distance_w_normal(
        first, second, myosin.left_end[myosin_index],
        myosin.right_end[myosin_index], box, actin.periodic_axes);
    result.distance = geometry_result.first;
    result.penetration = std::max(0.0, hard_core_distance - result.distance);
    const bool shell_enabled = shell_width > 0.0 && shell_energy > 0.0;
    const double shell_outer_distance = hard_core_distance +
        (shell_enabled ? shell_width : 0.0);
    result.shell_penetration = shell_enabled
        ? std::max(0.0, shell_outer_distance - result.distance)
        : 0.0;
    if (!(result.distance < shell_outer_distance)) return result;

    vec normal = geometry_result.second.normal;
    if (normal.norm() <= 1e-12) {
        vec midpoint = 0.5 * (first + second);
        vec displacement = midpoint - myosin.center[myosin_index];
        displacement.pbc_wrap(box, actin.periodic_axes);
        normal = regularized_steric_normal(
            connector, myosin.direction[myosin_index], displacement,
            actin_i, actin.n + myosin_index);
    }
    normal.normalize();

    // The soft corona is a conservative harmonic shoulder whose energy at
    // the hard-core surface is shell_energy.  Continuing the shoulder inside
    // the core keeps both energy and force continuous.  It represents the
    // crowded myosin-head layer, not a larger geometric rupture radius.
    double magnitude = 0.0;
    if (shell_enabled) {
        const double shell_k = 2.0 * shell_energy /
            (shell_width * shell_width);
        magnitude += shell_k * result.shell_penetration;
        result.energy += 0.5 * shell_k * result.shell_penetration *
            result.shell_penetration;
    }

    // The original WCA hard core is unchanged. Choose sigma so its cutoff is
    // exactly hard_core_distance, regularize deep overlap, and add its shifted
    // potential and force only inside that surface.
    if (epsilon > 0.0 && result.distance < hard_core_distance) {
        const double sigma = hard_core_distance /
            std::pow(2.0, 1.0 / 6.0);
        const double r = std::max(0.05 * sigma, result.distance);
        const double sr = sigma / r;
        const double sr2 = sr * sr;
        const double sr6 = sr2 * sr2 * sr2;
        const double sr12 = sr6 * sr6;
        magnitude += std::max(
            0.0, 24.0 * epsilon / r * (2.0 * sr12 - sr6));
        result.energy += std::max(
            0.0, 4.0 * epsilon * (sr12 - sr6) + epsilon);
    }
    if (force_cap > 0.0 && std::isfinite(force_cap)) {
        magnitude = std::min(magnitude, force_cap);
    }
    if (!(magnitude > 0.0) || !std::isfinite(magnitude)) return result;

    double fraction_j = 0.5;
    const double connector_length_sq = connector.norm_squared();
    if (connector_length_sq > 1e-16) {
        vec from_first = geometry_result.second.start - first;
        from_first.pbc_wrap(box, actin.periodic_axes);
        fraction_j = std::clamp(
            from_first.dot(connector) / connector_length_sq, 0.0, 1.0);
    }
    result.arm_weight_j = fraction_j;
    result.arm_weight_i = 1.0 - fraction_j;
    const vec connector_force = magnitude * normal;
    result.force_on_actin_i = result.arm_weight_i * connector_force;
    result.force_on_actin_j = result.arm_weight_j * connector_force;
    result.force_on_myosin = -connector_force;
    result.torque_on_actin_i = steric_contact_orientation_drive(
        actin.center[actin_i], actin.direction[actin_i], first,
        result.force_on_actin_i, box, actin.periodic_axes,
        actin.direction[actin_j]);
    result.torque_on_actin_j = steric_contact_orientation_drive(
        actin.center[actin_j], actin.direction[actin_j], second,
        result.force_on_actin_j, box, actin.periodic_axes,
        actin.direction[actin_i]);
    result.torque_on_myosin = steric_contact_orientation_drive(
        myosin.center[myosin_index], myosin.direction[myosin_index],
        geometry_result.second.end, result.force_on_myosin, box,
        myosin.periodic_axes, connector);
    result.force = magnitude;
    result.applied = true;
    return result;
}

ZDiscMyosinTipGap zdisc_myosin_tip_gap(
    const Filament& actin, const Myosin& myosin,
    int actin_i, int actin_j, double s_i, double s_j, int myosin_index,
    int end_sign, const std::vector<double>& box, double onset,
    double transverse_cutoff) {
    // Shared full-connector geometry for soft force and hard-contact queries.
    ZDiscMyosinTipGap result{};
    if (actin_i < 0 || actin_j < 0 || actin_i >= actin.n ||
        actin_j >= actin.n || myosin_index < 0 || myosin_index >= myosin.n ||
        (end_sign != -1 && end_sign != 1) || !(onset > 0.0) ||
        !(transverse_cutoff > 0.0)) {
        return result;
    }
    const vec center_i = actin.center[actin_i];
    const vec center_j = actin.center[actin_j];
    const vec first = center_i + s_i * actin.direction[actin_i];
    const vec second_raw = center_j + s_j * actin.direction[actin_j];
    vec connector = second_raw - first;
    connector.pbc_wrap(box, actin.periodic_axes);
    vec axis = myosin.direction[myosin_index];
    if (axis.norm() <= 1e-12) return result;
    axis.normalize();
    const vec outward = static_cast<double>(end_sign) * axis;
    const double half_length = 0.5 * myosin.length;
    if (!(half_length > 0.0)) return result;
    vec midpoint = first + 0.5 * connector - myosin.center[myosin_index];
    midpoint.pbc_wrap(box, myosin.periodic_axes);
    // A non-nearest periodic myosin image can intersect a long or tilted
    // connector. Enumerate every image within the conservative sphere.
    const double reach = std::hypot(half_length + onset, transverse_cutoff) +
        0.5 * connector.norm();
    const double coord[3]{midpoint.x, midpoint.y, midpoint.z};
    double period[3]{0.0, 0.0, 0.0};
    int lower[3]{0, 0, 0}, upper[3]{0, 0, 0};
    for (int k = 0; k < 3; ++k) {
        if (myosin.periodic_axes[k] &&
            k < static_cast<int>(box.size()) && box[k] > 0.0) {
            period[k] = box[k];
            lower[k] = static_cast<int>(std::ceil((-reach - coord[k]) / box[k]));
            upper[k] = static_cast<int>(std::floor((reach - coord[k]) / box[k]));
            if (lower[k] > upper[k]) return result;
        } else if (std::abs(coord[k]) > reach) {
            return result;
        }
    }
    for (int ix = lower[0]; ix <= upper[0]; ++ix)
    for (int iy = lower[1]; iy <= upper[1]; ++iy)
    for (int iz = lower[2]; iz <= upper[2]; ++iz) {
        const vec image_midpoint = midpoint +
            vec{ix * period[0], iy * period[1], iz * period[2]};
        if (image_midpoint.norm_squared() > reach * reach) continue;
        const vec start = image_midpoint - 0.5 * connector;
        const double z0 = start.dot(outward);
        const double dz = connector.dot(outward);
        double lo = 0.0, hi = 1.0;
        // Same-point axial eligibility: this tip's half, up to its onset.
        if (dz == 0.0) {
            if (z0 < 0.0 || z0 >= half_length + onset) continue;
        } else {
            double entry = -z0 / dz;
            double exit = (half_length + onset - z0) / dz;
            if (entry > exit) std::swap(entry, exit);
            lo = std::max(lo, entry);
            hi = std::min(hi, exit);
            if (lo > hi) continue;
        }
        const vec radial_start = start - z0 * outward;
        const vec radial_delta = connector - dz * outward;
        const double a = radial_delta.norm_squared();
        const double minimum = a > 0.0
            ? -radial_start.dot(radial_delta) / a
            : (dz < 0.0 ? hi : lo);
        double t = std::clamp(minimum, lo, hi);
        double rho = (radial_start + t * radial_delta).norm();
        double gap = z0 + t * dz - half_length;
        if (rho >= transverse_cutoff) continue;
        // If the radial minimum is exactly at the onset plane, choose an
        // interior point of the simultaneously eligible radial/axial interval.
        // Otherwise an eligible connector could falsely request zero force.
        if (gap >= onset && dz != 0.0 && lo < hi) {
            double inside_lo = lo, inside_hi = hi;
            if (a > 0.0) {
                const double minimum_sq =
                    (radial_start + minimum * radial_delta).norm_squared();
                const double span = std::sqrt(std::max(
                    0.0, (transverse_cutoff * transverse_cutoff - minimum_sq) / a));
                inside_lo = std::max(lo, minimum - span);
                inside_hi = std::min(hi, minimum + span);
            }
            t = 0.5 * (inside_lo + inside_hi);
            rho = (radial_start + t * radial_delta).norm();
            gap = z0 + t * dz - half_length;
        }
        if (!(rho < transverse_cutoff && gap < onset)) continue;
        // Prefer radial proximity; break equal-radius ties by penetration.
        if (!result.eligible || rho < result.transverse_distance ||
            (rho == result.transverse_distance && gap < result.axial_gap)) {
            result.eligible = true;
            result.axial_gap = gap;
            result.transverse_distance = rho;
        }
    }
    return result;
}

ZDiscMyosinTipContactResult compute_zdisc_myosin_tip_contact(
    const Filament& actin, const Myosin& myosin,
    int actin_i, int actin_j, double s_i, double s_j, int myosin_index,
    int end_sign, const std::vector<double>& box, double onset,
    double stiffness, double force_max, double transverse_flat,
    double transverse_cutoff) {
    ZDiscMyosinTipContactResult result{};
    if (actin_i < 0 || actin_j < 0 || actin_i >= actin.n ||
        actin_j >= actin.n || myosin_index < 0 || myosin_index >= myosin.n ||
        (end_sign != -1 && end_sign != 1) || !(onset > 0.0) ||
        !(stiffness > 0.0) || !(force_max > 0.0) ||
        transverse_flat < 0.0 || !(transverse_cutoff > transverse_flat)) {
        return result;
    }

    const vec center_i = actin.center[actin_i];
    const vec center_j = actin.center[actin_j];
    const vec first = center_i + s_i * actin.direction[actin_i];
    const vec second_raw = center_j + s_j * actin.direction[actin_j];
    vec connector = second_raw - first;
    connector.pbc_wrap(box, actin.periodic_axes);
    const vec second = first + connector;

    vec axis = myosin.direction[myosin_index];
    if (axis.norm() <= 1e-12) return result;
    axis.normalize();
    const vec outward = static_cast<double>(end_sign) * axis;
    const vec tip = end_sign < 0 ? myosin.left_end[myosin_index]
                                 : myosin.right_end[myosin_index];
    const auto geometry = zdisc_myosin_tip_gap(
        actin, myosin, actin_i, actin_j, s_i, s_j, myosin_index,
        end_sign, box, onset, transverse_cutoff);
    if (!geometry.eligible) return result;
    result.eligible_half = true;
    result.axial_gap = geometry.axial_gap;
    result.transverse_distance = geometry.transverse_distance;

    if (result.transverse_distance <= transverse_flat) {
        result.radial_gate = 1.0;
    } else {
        const double q = (result.transverse_distance - transverse_flat) /
            (transverse_cutoff - transverse_flat);
        result.radial_gate = (1.0 - q) * (1.0 - q) * (1.0 + 2.0 * q);
    }
    result.compression = onset - result.axial_gap;
    const double x = stiffness * result.compression / force_max;
    result.force = result.radial_gate * force_max * std::tanh(x);
    const double ax = std::abs(x);
    const double log_cosh = ax + std::log1p(std::exp(-2.0 * ax)) -
        std::log(2.0);
    result.energy = result.radial_gate * force_max * force_max / stiffness *
        log_cosh;
    if (!(result.force > 0.0) || !std::isfinite(result.force)) return result;

    const vec patch_force = result.force * outward;
    result.force_on_actin_i = 0.5 * patch_force;
    result.force_on_actin_j = 0.5 * patch_force;
    result.force_on_myosin = -patch_force;
    result.torque_on_actin_i = steric_contact_orientation_drive(
        actin.center[actin_i], actin.direction[actin_i], first,
        result.force_on_actin_i, box, actin.periodic_axes,
        actin.direction[actin_j]);
    result.torque_on_actin_j = steric_contact_orientation_drive(
        actin.center[actin_j], actin.direction[actin_j], second,
        result.force_on_actin_j, box, actin.periodic_axes,
        actin.direction[actin_i]);
    result.torque_on_myosin = steric_contact_orientation_drive(
        myosin.center[myosin_index], myosin.direction[myosin_index], tip,
        result.force_on_myosin, box, myosin.periodic_axes, connector);
    result.applied = true;
    return result;
}

namespace {

// Gauss-Legendre rule on [0, 1]; weights sum to 1.
struct UnitGaussRule {
    std::vector<double> nodes;
    std::vector<double> weights;
};

UnitGaussRule make_unit_gauss_rule(int n) {
    UnitGaussRule rule;
    rule.nodes.resize(static_cast<size_t>(n));
    rule.weights.resize(static_cast<size_t>(n));
    const double pi = std::acos(-1.0);
    for (int i = 0; i < n; ++i) {
        // Newton iteration on P_n from the usual Chebyshev-like guess.
        double x = std::cos(pi * (i + 0.75) / (n + 0.5));
        double derivative = 1.0;
        for (int iteration = 0; iteration < 100; ++iteration) {
            double previous = 1.0, current = x;
            for (int k = 2; k <= n; ++k) {
                const double next =
                    ((2.0 * k - 1.0) * x * current - (k - 1.0) * previous) / k;
                previous = current;
                current = next;
            }
            derivative = n == 1 ? 1.0
                : n * (x * current - previous) / (x * x - 1.0);
            const double step = current / derivative;
            x -= step;
            if (std::abs(step) < 1e-15) break;
        }
        rule.nodes[static_cast<size_t>(i)] = 0.5 * (1.0 - x);
        rule.weights[static_cast<size_t>(i)] =
            1.0 / ((1.0 - x * x) * derivative * derivative);
    }
    return rule;
}

const UnitGaussRule& unit_gauss_rule(int n) {
    static const UnitGaussRule production_rule = make_unit_gauss_rule(16);
    if (n == 16) return production_rule;
    // Other orders are for refinement tests; a deque keeps references stable.
    thread_local std::deque<std::pair<int, UnitGaussRule>> cache;
    for (const auto& entry : cache) {
        if (entry.first == n) return entry.second;
    }
    cache.emplace_back(n, make_unit_gauss_rule(n));
    return cache.back().second;
}

}  // namespace

ZDiscAxialStopContactResult compute_zdisc_axial_stop_contact(
    const Filament& actin, const Myosin& myosin,
    int actin_i, int actin_j, double s_i, double s_j, int myosin_index,
    int end_sign, const std::vector<double>& box, double onset,
    double stiffness, double force_max, double gate_radius,
    double gate_cutoff, double epsilon, int quadrature_nodes) {
    ZDiscAxialStopContactResult result{};
    if (actin_i < 0 || actin_j < 0 || actin_i >= actin.n ||
        actin_j >= actin.n || myosin_index < 0 || myosin_index >= myosin.n ||
        (end_sign != -1 && end_sign != 1) || !(onset > 0.0) ||
        !(stiffness > 0.0) || !(force_max > 0.0) || !(gate_radius >= 0.0) ||
        !(gate_cutoff > gate_radius) || !(epsilon > 0.0) ||
        quadrature_nodes < 1) {
        return result;
    }
    vec axis = myosin.direction[myosin_index];
    if (axis.norm() <= 1e-12) return result;
    axis.normalize();
    const vec outward = static_cast<double>(end_sign) * axis;
    const double half_length = 0.5 * myosin.length;
    const double contact_half = half_length + onset;  // H
    if (!(half_length > 0.0) || !(contact_half > epsilon)) return result;
    // sqrt(z^2 + epsilon^2) < H exactly when z < axial_support.
    const double axial_support =
        std::sqrt(contact_half * contact_half - epsilon * epsilon);
    const double saturation_depth = force_max / stiffness;

    const vec center_i = actin.center[actin_i];
    const vec center_j = actin.center[actin_j];
    const vec first = center_i + s_i * actin.direction[actin_i];
    const vec second_raw = center_j + s_j * actin.direction[actin_j];
    vec connector = second_raw - first;
    connector.pbc_wrap(box, actin.periodic_axes);
    const vec second = first + connector;
    vec midpoint = first + 0.5 * connector - myosin.center[myosin_index];
    midpoint.pbc_wrap(box, myosin.periodic_axes);

    // Every periodic myosin image within the conservative sphere, as in
    // zdisc_myosin_tip_gap(); distinct images are distinct contacts and add.
    const double reach = std::hypot(contact_half, gate_cutoff) +
        0.5 * connector.norm();
    const double coord[3]{midpoint.x, midpoint.y, midpoint.z};
    double period[3]{0.0, 0.0, 0.0};
    int lower[3]{0, 0, 0}, upper[3]{0, 0, 0};
    for (int k = 0; k < 3; ++k) {
        if (myosin.periodic_axes[k] &&
            k < static_cast<int>(box.size()) && box[k] > 0.0) {
            period[k] = box[k];
            lower[k] = static_cast<int>(std::ceil((-reach - coord[k]) / box[k]));
            upper[k] = static_cast<int>(std::floor((reach - coord[k]) / box[k]));
            if (lower[k] > upper[k]) return result;
        } else if (std::abs(coord[k]) > reach) {
            return result;
        }
    }

    const auto& rule = unit_gauss_rule(quadrature_nodes);
    const double dz = connector.dot(outward);
    const vec radial_delta = connector - dz * outward;
    const double radial_a = radial_delta.norm_squared();
    vec force_i{0.0, 0.0, 0.0}, force_j{0.0, 0.0, 0.0};
    vec torque_m{0.0, 0.0, 0.0};
    double energy = 0.0;
    bool have_deepest = false;
    for (int ix = lower[0]; ix <= upper[0]; ++ix)
    for (int iy = lower[1]; iy <= upper[1]; ++iy)
    for (int iz = lower[2]; iz <= upper[2]; ++iz) {
        const vec image_midpoint = midpoint +
            vec{ix * period[0], iy * period[1], iz * period[2]};
        if (image_midpoint.norm_squared() > reach * reach) continue;
        // Connector relative to this image's myosin center.
        const vec start = image_midpoint - 0.5 * connector;
        const double z0 = start.dot(outward);
        // Support on this half: 0 <= z(t) < axial_support, rho(t) < cutoff.
        double lo = 0.0, hi = 1.0;
        if (dz == 0.0) {
            if (z0 < 0.0 || z0 >= axial_support) continue;
        } else {
            double entry = -z0 / dz;
            double exit = (axial_support - z0) / dz;
            if (entry > exit) std::swap(entry, exit);
            lo = std::max(lo, entry);
            hi = std::min(hi, exit);
        }
        const vec radial_start = start - z0 * outward;
        const double radial_b = radial_start.dot(radial_delta);
        const double radial_c =
            radial_start.norm_squared() - gate_cutoff * gate_cutoff;
        if (radial_a > 0.0) {
            const double discriminant =
                radial_b * radial_b - radial_a * radial_c;
            if (!(discriminant > 0.0)) continue;
            const double root = std::sqrt(discriminant);
            lo = std::max(lo, (-radial_b - root) / radial_a);
            hi = std::min(hi, (-radial_b + root) / radial_a);
        } else if (!(radial_c < 0.0)) {
            continue;
        }
        // A zero-length connector keeps [0, 1]: its point is either inside
        // the support or rejected above.
        if (!(hi > lo)) continue;
        // Split where the integrand is not smooth, so each piece converges
        // spectrally: Phi changes from quadratic to linear (a kink in the
        // force), and the gate leaves 1 at gate_radius.
        double breaks[5]{lo, 0.0, 0.0, 0.0, 0.0};
        int n_breaks = 1;
        auto add_break = [&](double t) {
            if (t > lo && t < hi) breaks[n_breaks++] = t;
        };
        const double saturation_eps = contact_half - saturation_depth;
        if (dz != 0.0 && saturation_eps > epsilon) {
            const double z_saturation = std::sqrt(
                saturation_eps * saturation_eps - epsilon * epsilon);
            add_break((z_saturation - z0) / dz);
        }
        if (radial_a > 0.0 && gate_radius > 0.0) {
            const double inner_c = radial_start.norm_squared() -
                gate_radius * gate_radius;
            const double discriminant =
                radial_b * radial_b - radial_a * inner_c;
            if (discriminant > 0.0) {
                const double root = std::sqrt(discriminant);
                add_break((-radial_b - root) / radial_a);
                add_break((-radial_b + root) / radial_a);
            }
        }
        breaks[n_breaks++] = hi;
        std::sort(breaks + 1, breaks + n_breaks - 1);
        for (int piece = 0; piece + 1 < n_breaks; ++piece) {
            const double piece_lo = breaks[piece];
            const double span = breaks[piece + 1] - piece_lo;
            if (!(span > 0.0)) continue;
            for (size_t k = 0; k < rule.nodes.size(); ++k) {
                const double t = piece_lo + span * rule.nodes[k];
                const double weight = span * rule.weights[k];
                const vec x = start + t * connector;
                const double z = x.dot(outward);
                if (z < 0.0) continue;
                const double z_eps = std::sqrt(z * z + epsilon * epsilon);
                const double depth = contact_half - z_eps;
                if (!(depth > 0.0)) continue;
                const double rho = (x - z * outward).norm();
                if (!(rho < gate_cutoff)) continue;
                double gate = 1.0;
                if (rho > gate_radius) {
                    const double q =
                        (rho - gate_radius) / (gate_cutoff - gate_radius);
                    gate = (1.0 - q) * (1.0 - q) * (1.0 + 2.0 * q);
                }
                const double phi = depth <= saturation_depth
                    ? 0.5 * stiffness * depth * depth
                    : force_max * (depth - 0.5 * saturation_depth);
                const double slope = std::min(stiffness * depth, force_max);
                // Axial part of -grad u only; no radial component.
                const vec f = (gate * slope * z / z_eps) * outward;
                energy += weight * gate * phi;
                force_i += (weight * (1.0 - t)) * f;
                force_j += (weight * t) * f;
                torque_m += (-weight) * x.cross(f);
                const double gap = z - half_length;
                if (!have_deepest || gap < result.deepest_gap) {
                    have_deepest = true;
                    result.deepest_gap = gap;
                    result.deepest_rho = rho;
                    result.deepest_gate = gate;
                }
            }
        }
    }
    if (!(energy > 0.0) || !std::isfinite(energy)) return result;

    result.energy = energy;
    result.force_on_actin_i = force_i;
    result.force_on_actin_j = force_j;
    result.force_on_myosin = -(force_i + force_j);
    result.axial_force = result.force_on_myosin.norm();
    // Forces act at the material attachment sites of the two actins.
    vec lever_i = first - center_i;
    lever_i.pbc_wrap(box, actin.periodic_axes);
    vec lever_j = second - center_j;
    lever_j.pbc_wrap(box, actin.periodic_axes);
    vec axis_i = actin.direction[actin_i];
    vec axis_j = actin.direction[actin_j];
    axis_i.normalize();
    axis_j.normalize();
    result.physical_torque_on_actin_i = lever_i.cross(force_i);
    result.physical_torque_on_actin_j = lever_j.cross(force_j);
    result.physical_torque_on_myosin = torque_m;
    result.torque_on_actin_i = result.physical_torque_on_actin_i.cross(axis_i);
    result.torque_on_actin_j = result.physical_torque_on_actin_j.cross(axis_j);
    result.torque_on_myosin = torque_m.cross(axis);
    result.applied = true;
    return result;
}

BarbedTipMyosinEndstopResult compute_barbed_tip_myosin_endstop(
    const Filament& actin, const Myosin& myosin,
    int actin_index, int myosin_index, const std::vector<double>& box,
    double capsule_radius, double onset, double stiffness,
    double force_max, bool slab_exclusion, double transverse_cutoff) {
    BarbedTipMyosinEndstopResult result{};
    if (actin_index < 0 || actin_index >= actin.n || myosin_index < 0 ||
        myosin_index >= myosin.n || !(capsule_radius > 0.0) ||
        !(onset > 0.0) || !(stiffness > 0.0) || !(force_max > 0.0)) {
        return result;
    }
    vec axis = myosin.direction[myosin_index];
    if (axis.norm() <= 1e-12) return result;
    axis.normalize();
    const vec barbed_tip = actin.left_end[actin_index];
    vec relative = barbed_tip - myosin.center[myosin_index];
    relative.pbc_wrap(box, myosin.periodic_axes);
    result.axial_coordinate = relative.dot(axis);
    result.end_sign = result.axial_coordinate < 0.0 ? -1 : 1;
    const vec transverse = relative - result.axial_coordinate * axis;
    result.transverse_distance = transverse.norm();
    // SLAB MODE.  The A-band a barbed tip must stay out of is a slab, and a
    // slab has no sides: it can only be entered through one of its two
    // faces, so the onset shell covers every approach and penetration
    // always grows continuously from zero.  Two things change.
    //
    // (1) The transverse HARD CUTOFF at capsule_radius becomes a wide
    //     locality bound.  In production capsule_radius (0.080) is exactly
    //     the resting A-M centerline (myosin_radius 0.050 + actin_radius
    //     0.015 + am_surface_gap 0.015), so every bound actin sits ON the
    //     gate boundary and dithers across it ~1000 times/s, and any tip
    //     arriving from outside it meets no opposition at all until the
    //     gate flips -- MEASURED: 3441 production cases went from gate-shut
    //     to past the kill depth in one frame, median 347.6 nm deep, 100%
    //     already fatal the first time the code could see them, with the
    //     force at 0.0% of the preceding 600 us.
    // (2) The spherocylinder cap extension is dropped, so the exclusion
    //     boundary is the flat A-band face.  cap_extension makes the axial
    //     boundary a function of the transverse distance, which gives this
    //     purely axial force a nonzero curl -- F_z depends on r while F_r
    //     is identically zero, so cycling a tip in and out at two radii
    //     does net work.  With the boundary at 0.5*length the force is a
    //     function of the axial coordinate alone and is a true gradient.
    const double eligibility_radius =
        (slab_exclusion && transverse_cutoff > 0.0) ? transverse_cutoff
                                                    : capsule_radius;
    if (result.transverse_distance >= eligibility_radius) return result;
    result.radially_eligible = true;
    const double cap_extension = slab_exclusion ? 0.0 : std::sqrt(std::max(
        0.0, capsule_radius * capsule_radius -
                 result.transverse_distance * result.transverse_distance));
    result.capsule_axial_boundary = 0.5 * myosin.length + cap_extension;
    result.signed_clearance = std::abs(result.axial_coordinate) -
                              result.capsule_axial_boundary;
    result.penetration = std::max(0.0, -result.signed_clearance);
    if (result.signed_clearance >= onset) return result;
    result.compression = onset - result.signed_clearance;
    const double x = stiffness * result.compression / force_max;
    result.force = force_max * std::tanh(x);
    if (!(result.force > 0.0) || !std::isfinite(result.force)) return result;
    const vec outward = static_cast<double>(result.end_sign) * axis;
    result.force_on_actin = result.force * outward;
    result.force_on_myosin = -result.force_on_actin;
    result.torque_on_actin = steric_contact_orientation_drive(
        actin.center[actin_index], actin.direction[actin_index], barbed_tip,
        result.force_on_actin, box, actin.periodic_axes, outward);
    const vec myosin_center = myosin.center[myosin_index];
    const vec myosin_contact = myosin_center +
        static_cast<double>(result.end_sign) * 0.5 * myosin.length * axis;
    result.torque_on_myosin = steric_contact_orientation_drive(
        myosin.center[myosin_index], myosin.direction[myosin_index],
        myosin_contact, result.force_on_myosin, box,
        myosin.periodic_axes, actin.direction[actin_index]);
    result.applied = true;
    return result;
}

AMOffCenterTorqueResult compute_am_offcenter_torque(
    const vec& center,
    const vec& axis,
    const vec& attachment,
    const vec& force,
    const std::vector<double>& box,
    const std::array<bool,3>& periodic) {
    AMOffCenterTorqueResult result;
    result.lever_arm = attachment - center;
    result.lever_arm.pbc_wrap(box, periodic);
    result.physical_torque = result.lever_arm.cross(force);

    vec unit_axis = axis;
    if (unit_axis.norm() <= 1e-12) return result;
    unit_axis.normalize();
    result.orientation_drive = result.physical_torque.cross(unit_axis);
    return result;
}

double project_nematic_alignment_drive(
    const vec& moving_axis, const vec& target_axis,
    const vec& orientation_drive) {
    vec moving = moving_axis;
    vec target = target_axis;
    if (moving.norm() <= 1e-12 || target.norm() <= 1e-12) return 0.0;
    moving.normalize();
    target.normalize();
    if (moving.dot(target) < 0.0) target = -1.0 * target;
    vec tangent = target - moving * target.dot(moving);
    if (tangent.norm() <= 1e-12) return 0.0;
    tangent.normalize();
    return orientation_drive.dot(tangent);
}

AMNematicTetherResult compute_am_nematic_flat_bottom_tether(
    const vec& actin_axis, const vec& myosin_axis,
    double stiffness, double slack_angle) {
    AMNematicTetherResult result;
    vec axis_a = actin_axis;
    vec axis_m = myosin_axis;
    if (!(stiffness > 0.0) ||
        axis_a.norm() <= 1e-12 || axis_m.norm() <= 1e-12) {
        return result;
    }
    axis_a.normalize();
    axis_m.normalize();
    const double dot = std::clamp(axis_a.dot(axis_m), -1.0, 1.0);
    const double sine = std::sqrt(std::max(0.0, 1.0 - dot * dot));
    result.angle = std::asin(std::clamp(sine, 0.0, 1.0));
    const double bounded_slack = std::clamp(
        slack_angle, 0.0, 0.5 * M_PI);
    const double slack_sine = std::sin(bounded_slack);
    if (sine <= slack_sine || sine <= 1e-12) return result;

    const double excess = sine - slack_sine;
    result.energy = 0.5 * stiffness * excess * excess;
    // Negative constrained gradient of the nematic flat-bottom potential.
    const double coefficient = stiffness * excess * dot / sine;
    result.drive_on_actin =
        coefficient * (axis_m - dot * axis_a);
    result.drive_on_myosin =
        coefficient * (axis_a - dot * axis_m);
    result.active = true;
    return result;
}

RepulsionResult compute_myosin_repulsion(const Filament& actin,
                                         const Myosin& myosin,
                                         int i,
                                         int j,
                                         const std::vector<double>& box,
                                         int fix_myosin,
                                         const utils::MoleculeConnection& actinIndicesPerMyosin,
                                         double stiffness,
                                         double exclusion_distance)
{
    RepulsionResult result{};
    // exclusion_distance <= 0 preserves the historical hard-core-only range.
    const double cutoff = (exclusion_distance > 0.0)
                              ? exclusion_distance
                              : 2.0 * myosin.radius;
    const double EPS = 1e-9;

    vec center_displacement = myosin.center[i] - myosin.center[j];
    center_displacement.pbc_wrap(box, myosin.periodic_axes);

    auto geom = geometry::segment_segment_distance_w_normal(
        myosin.left_end[i], myosin.right_end[i],
        myosin.left_end[j], myosin.right_end[j],
        box, myosin.periodic_axes);

    const double distance = geom.first;
    result.segment_distance = distance;
    if (distance >= cutoff) {
        return result;
    }

    vec normal_vector = geom.second.normal;
    double norm = normal_vector.norm();
    if (norm <= EPS) {
        normal_vector = regularized_steric_normal(
            myosin.direction[i], myosin.direction[j], center_displacement,
            i, j);
    } else {
        normal_vector = normal_vector / norm;
    }

    double overlap = (cutoff - distance)/cutoff;
    if (overlap <= 0.0) {
        return result;
    }

    double effective_stiffness = (std::isfinite(stiffness) && stiffness > 0.0) ? stiffness : 0.0;
    if (effective_stiffness <= 0.0) {
        return result;
    }
    double force_mag = effective_stiffness * overlap;
    if (!std::isfinite(force_mag) || force_mag < 0.0) {
        force_mag = 0.0;
    }

    vec force_vec = force_mag * normal_vector;

    result.applied = true;

    if (i < fix_myosin) {
        result.force_on_second -= force_vec * 2.0;
    } else if (j < fix_myosin) {
        result.force_on_first += force_vec * 2.0;
    } else {
        // Biochemical occupancy does not make a thick filament immobile.
        // Excluded volume must remain reciprocal unless the caller explicitly
        // fixes a filament; otherwise binding injects net momentum and lets
        // unlike-status rods interpenetrate.
        result.force_on_first += force_vec;
        result.force_on_second -= force_vec;
    }
    result.torque_on_first = steric_contact_orientation_drive(
        myosin.center[i], myosin.direction[i], geom.second.start,
        result.force_on_first, box, myosin.periodic_axes,
        myosin.direction[j]);
    result.torque_on_second = steric_contact_orientation_drive(
        myosin.center[j], myosin.direction[j], geom.second.end,
        result.force_on_second, box, myosin.periodic_axes,
        myosin.direction[i]);
    // if (overlap > 0.1) {
    //     printf("Myosin-myosin repulsion between myosins %d and %d: distance=%f, overlap=%f, force_mag=%f,force_vec=(%f, %f, %f)\n",
    //            i, j, distance, overlap, force_mag,
    //             force_vec.x, force_vec.y, force_vec.z);
    // }
    return result;
}

ContactProjectionResult project_myosin_hard_contact(
    Myosin& myosin,
    int i,
    int j,
    const std::vector<double>& box,
    int fix_myosin,
    double target_distance,
    double translation_weight,
    double rotation_weight)
{
    ContactProjectionResult result{};
    if (i == j || i < 0 || j < 0 || i >= myosin.n || j >= myosin.n ||
        !std::isfinite(target_distance) || target_distance <= 0.0) {
        return result;
    }

    const auto geom = geometry::segment_segment_distance_w_normal(
        myosin.left_end[i], myosin.right_end[i],
        myosin.left_end[j], myosin.right_end[j],
        box, myosin.periodic_axes);
    result.distance_before = geom.first;
    if (!std::isfinite(result.distance_before) ||
        result.distance_before >= target_distance) {
        return result;
    }

    const bool first_fixed = i < fix_myosin;
    const bool second_fixed = j < fix_myosin;
    if (first_fixed && second_fixed) return result;

    vec center_displacement = myosin.center[i] - myosin.center[j];
    center_displacement.pbc_wrap(box, myosin.periodic_axes);
    vec normal = geom.second.normal;
    if (normal.norm() <= 1e-9) {
        normal = regularized_steric_normal(
            myosin.direction[i], myosin.direction[j],
            center_displacement, i, j);
    } else {
        normal.normalize();
    }

    // A tiny clearance avoids immediately re-entering the soft potential
    // through roundoff.  This is 50 pm for a 50 nm thick-filament diameter.
    const double clearance = std::max(1e-12, target_distance * 1e-6);
    const double gap = target_distance - result.distance_before + clearance;
    const double wt = (std::isfinite(translation_weight) &&
                       translation_weight > 0.0)
        ? translation_weight : 1.0;
    const double wr = (std::isfinite(rotation_weight) &&
                       rotation_weight > 0.0)
        ? rotation_weight : 0.0;
    const vec first_orientation_gradient = steric_contact_orientation_drive(
        myosin.center[i], myosin.direction[i], geom.second.start,
        normal, box, myosin.periodic_axes, myosin.direction[j]);
    const vec second_orientation_gradient = steric_contact_orientation_drive(
        myosin.center[j], myosin.direction[j], geom.second.end,
        -normal, box, myosin.periodic_axes, myosin.direction[i]);
    double denominator = 0.0;
    if (!first_fixed) {
        denominator += wt + wr * first_orientation_gradient.dot(
            first_orientation_gradient);
    }
    if (!second_fixed) {
        denominator += wt + wr * second_orientation_gradient.dot(
            second_orientation_gradient);
    }
    if (!(denominator > 0.0) || !std::isfinite(denominator)) return result;
    const double multiplier = gap / denominator;

    vec first_delta = wt * multiplier * normal;
    vec second_delta = -wt * multiplier * normal;
    if (!first_fixed) {
        double dx = first_delta.x, dy = first_delta.y, dz = first_delta.z;
        myosin.displace(i, dx, dy, dz);
        if (wr > 0.0 && first_orientation_gradient.norm() > 0.0) {
            vec axis = myosin.direction[i];
            axis.normalize();
            myosin.direction[i] = rod_mobility::sphere_exp(
                axis, wr * multiplier * first_orientation_gradient);
            myosin.update_endpoints(i);
        }
    }
    if (!second_fixed) {
        double dx = second_delta.x, dy = second_delta.y, dz = second_delta.z;
        myosin.displace(j, dx, dy, dz);
        if (wr > 0.0 && second_orientation_gradient.norm() > 0.0) {
            vec axis = myosin.direction[j];
            axis.normalize();
            myosin.direction[j] = rod_mobility::sphere_exp(
                axis, wr * multiplier * second_orientation_gradient);
            myosin.update_endpoints(j);
        }
    }
    result.applied = true;
    result.correction = gap;
    return result;
}

RepulsionResult compute_actin_repulsion(const Filament& actin,
                                        int i,
                                        int j,
                                        const std::vector<double>& box,
                                        double crosslinker_length,
                                        double stiffness)
{
    RepulsionResult result{};
    const double EPS = 1e-9;

    vec center_displacement = actin.center[i] - actin.center[j];
    center_displacement.pbc_wrap(box, actin.periodic_axes);
    const double center_distance = center_displacement.norm();
    if (center_distance > crosslinker_length + actin.length) {
        return result;
    }

    auto geom = geometry::segment_segment_distance_w_normal(
        actin.left_end[i], actin.right_end[i],
        actin.left_end[j], actin.right_end[j],
        box, actin.periodic_axes);

    const double distance = geom.first;
    result.segment_distance = distance;
    if (distance >= crosslinker_length) {
        return result;
    }

    vec normal_vector = geom.second.normal;
    double norm = normal_vector.norm();
    if (norm <= EPS) {
        normal_vector = regularized_steric_normal(
            actin.direction[i], actin.direction[j], center_displacement,
            i, j);
    } else {
        normal_vector = normal_vector / norm;
    }

    double overlap = crosslinker_length - distance;
    if (overlap <= 0.0) {
        return result;
    }

    double force_mag = stiffness * overlap;
    if (!std::isfinite(force_mag) || force_mag < 0.0) {
        force_mag = 0.0;
    }
    vec force_vec = force_mag * normal_vector;

    // Excluded volume is independent of biochemical catch state.  Always
    // apply an equal-and-opposite pair force; suppressing the reaction on a
    // status-2 filament injects momentum and lets unlike-status rods overlap.
    result.force_on_first += force_vec;
    result.force_on_second -= force_vec;
    result.torque_on_first = steric_contact_orientation_drive(
        actin.center[i], actin.direction[i], geom.second.start,
        result.force_on_first, box, actin.periodic_axes,
        actin.direction[j]);
    result.torque_on_second = steric_contact_orientation_drive(
        actin.center[j], actin.direction[j], geom.second.end,
        result.force_on_second, box, actin.periodic_axes,
        actin.direction[i]);
    result.applied = true;
    return result;
}

RepulsionResult compute_barbed_end_lateral_repulsion(
    const Filament& actin,
    int i,
    int j,
    const std::vector<double>& box,
    double range,
    double stiffness,
    double force_cap,
    double axial_window)
{
    RepulsionResult result{};
    if (i == j || i < 0 || j < 0 || i >= actin.n || j >= actin.n) return result;
    if (!(range > 0.0) || !(stiffness > 0.0) ||
        !std::isfinite(range) || !std::isfinite(stiffness) ||
        !(axial_window > 0.0) || !std::isfinite(axial_window)) {
        return result;
    }

    vec ui = actin.direction[i];
    vec uj = actin.direction[j];
    if (!(ui.norm() > 0.0) || !(uj.norm() > 0.0)) return result;
    ui.normalize();
    uj.normalize();

    // Filament::direction points from the barbed end to the pointed end, so
    // the material barbed-end cap is center - L/2 * direction.
    const vec barbed_i = static_cast<vec>(actin.center[i]) -
                         0.5 * actin.length * ui;
    const vec barbed_j = static_cast<vec>(actin.center[j]) -
                         0.5 * actin.length * uj;
    vec separation = barbed_j - barbed_i;
    separation.pbc_wrap(box, actin.periodic_axes);

    // Common axis of the pair, independent of relative polarity.  Antiparallel
    // partners are the physically interesting case for ACTN2, and negating uj
    // when the pair opposes keeps the axis well defined for both.
    vec axis = ui + (ui.dot(uj) >= 0.0 ? uj : -1.0 * uj);
    if (!(axis.norm() > 1e-12)) return result;
    axis.normalize();

    const double axial_offset = separation.dot(axis);
    result.axial_gap = axial_offset;
    if (std::abs(axial_offset) > axial_window) return result;

    const vec transverse = separation - axial_offset * axis;
    const double distance = transverse.norm();
    result.transverse_distance = distance;
    result.segment_distance = distance;
    if (!(distance < range)) return result;
    // Two barbed ends sharing a lattice site have no defined push direction.
    // Leave those to the isotropic capsule sterics rather than inventing one.
    const double minimum_distance = 1e-4 * range;
    if (!(distance > minimum_distance)) return result;

    const double overlap = range - distance;
    result.contact_penetration = overlap;
    double force_magnitude = stiffness * overlap;
    if (force_cap > 0.0 && std::isfinite(force_cap)) {
        force_magnitude = std::min(force_magnitude, force_cap);
    }
    if (!(force_magnitude > 0.0) || !std::isfinite(force_magnitude)) {
        return result;
    }

    // transverse points from i toward j, so pushing apart drives i backwards
    // along it.  The pair is reciprocal and injects no net momentum.
    const vec direction = (1.0 / distance) * transverse;
    const vec force_on_i = -force_magnitude * direction;
    const vec force_on_j = -1.0 * force_on_i;
    result.force_on_first = force_on_i;
    result.force_on_second = force_on_j;
    result.torque_on_first = steric_contact_orientation_drive(
        actin.center[i], ui, barbed_i, force_on_i,
        box, actin.periodic_axes, uj);
    result.torque_on_second = steric_contact_orientation_drive(
        actin.center[j], uj, barbed_j, force_on_j,
        box, actin.periodic_axes, ui);
    result.applied = true;
    return result;
}

RepulsionResult compute_opposing_pointed_end_contact(
    const Filament& actin,
    int i,
    int j,
    const std::vector<double>& box,
    double contact_diameter,
    double epsilon,
    double force_cap,
    double antiparallel_dot_max)
{
    RepulsionResult result{};
    if (i == j || i < 0 || j < 0 || i >= actin.n || j >= actin.n ||
        !(contact_diameter > 0.0) || !(epsilon > 0.0) ||
        !std::isfinite(contact_diameter) || !std::isfinite(epsilon)) {
        return result;
    }

    vec ui = actin.direction[i];
    vec uj = actin.direction[j];
    if (!(ui.norm() > 0.0) || !(uj.norm() > 0.0)) return result;
    ui.normalize();
    uj.normalize();
    const double direction_dot = ui.dot(uj);
    if (direction_dot > antiparallel_dot_max) return result;

    // Filament::direction points from the barbed end to the pointed end, so
    // center + L/2 direction is the material pointed-end cap.
    const vec pointed_i = static_cast<vec>(actin.center[i]) +
                          0.5 * actin.length * ui;
    const vec pointed_j = static_cast<vec>(actin.center[j]) +
                          0.5 * actin.length * uj;
    vec separation = pointed_j - pointed_i;
    separation.pbc_wrap(box, actin.periodic_axes);

    // e points from i toward its pointed end and, for an antiparallel pair,
    // toward j's pointed end before collision.  Preserve this signed order:
    // after an attempted crossing axial_gap becomes negative instead of being
    // mistaken for two separated endpoint spheres.
    vec e = ui - uj;
    if (!(e.norm() > 1e-12)) return result;
    e.normalize();
    const double axial_gap = separation.dot(e);
    const vec transverse = separation - axial_gap * e;
    const double transverse_distance = transverse.norm();
    result.axial_gap = axial_gap;
    result.transverse_distance = transverse_distance;
    result.segment_distance = std::sqrt(
        axial_gap * axial_gap +
        transverse_distance * transverse_distance);

    // Two circular end caps on distinct lattice tracks do not collide.  For
    // caps in the same lane, their transverse offset leaves this much allowed
    // axial clearance at first contact.
    if (transverse_distance >= contact_diameter) return result;
    const double allowed_gap = std::sqrt(std::max(
        0.0, contact_diameter * contact_diameter -
             transverse_distance * transverse_distance));
    const double penetration = allowed_gap - axial_gap;
    result.contact_penetration = std::max(0.0, penetration);
    if (!(penetration > 0.0)) return result;

    // Shift the WCA radial coordinate so the force is exactly zero at first
    // contact (penetration=0) and rises steeply only under compression.
    const double wca_cutoff = std::pow(2.0, 1.0 / 6.0) * contact_diameter;
    const double minimum_r = 0.05 * contact_diameter;
    const double r = std::max(minimum_r, wca_cutoff - penetration);
    const double sr = contact_diameter / r;
    const double sr2 = sr * sr;
    const double sr6 = sr2 * sr2 * sr2;
    const double sr12 = sr6 * sr6;
    double force_magnitude =
        24.0 * epsilon / r * (2.0 * sr12 - sr6);
    force_magnitude = std::max(0.0, force_magnitude);
    if (force_cap > 0.0 && std::isfinite(force_cap)) {
        force_magnitude = std::min(force_magnitude, force_cap);
    }
    if (!(force_magnitude > 0.0) || !std::isfinite(force_magnitude)) {
        return result;
    }

    // Increase the signed gap.  The pair is reciprocal and therefore cannot
    // inject net momentum.  Applying the forces at the material caps also
    // supplies the correct off-axis rotation drive when the rods are skewed.
    const vec force_on_i = -force_magnitude * e;
    const vec force_on_j = -force_on_i;
    result.force_on_first = force_on_i;
    result.force_on_second = force_on_j;
    result.torque_on_first = steric_contact_orientation_drive(
        actin.center[i], ui, pointed_i, force_on_i,
        box, actin.periodic_axes, uj);
    result.torque_on_second = steric_contact_orientation_drive(
        actin.center[j], uj, pointed_j, force_on_j,
        box, actin.periodic_axes, ui);
    result.applied = true;
    return result;
}

vec compute_actin_myosin_repulsion(const Filament& actin,
                                   const Myosin& myosin,
                                   int act_idx,
                                   int myo_idx,
                                   const std::vector<double>& box,
                                   double radius,
                                   double stiffness)
{
    // const double EPS = 1e-9;
    // const double max_force_limit = (max_force_cap > 0.0 && std::isfinite(max_force_cap))
    //                                    ? max_force_cap
    //                                    : std::numeric_limits<double>::infinity();

    auto geom = geometry::segment_segment_distance_w_normal(
        actin.left_end[act_idx], actin.right_end[act_idx],
        myosin.left_end[myo_idx], myosin.right_end[myo_idx],
        box, actin.periodic_axes);

    const double distance = geom.first;
    if (distance >= radius) {
        return {0.0, 0.0, 0.0};
    }

    // auto normal_it = geom.second.find("normal");
    // vec normal_vector = normal_it->second;
    // double norm = normal_vector.norm();
    // vec dir;
    // if (norm <= EPS) {
    //     vec center_displacement = actin.center[act_idx] - myosin.center[myo_idx];
    //     center_displacement.pbc_wrap(box, actin.periodic_axes);
    //     double center_norm = center_displacement.norm();
    //     dir = center_displacement / center_norm;
    // } else {
    //     dir = normal_vector / norm;
    // }

    // double overlap = radius - distance;

    // double magnitude = stiffness * overlap;
    
    // if (magnitude > max_force_limit) {
    //     magnitude = max_force_limit;
    // }
    // vec repulsive_force = magnitude * dir;


    vec u = myosin.direction[myo_idx];
    // Signed axial coordinate of the actin *left* tip relative to myosin center
    vec tip_disp = actin.left_end[act_idx] - myosin.center[myo_idx];
    tip_disp.pbc_wrap(box, actin.periodic_axes);
    double s  = tip_disp.dot(u);                         // +s toward "right" tip, −s toward "left"
    double Lm = 0.5 * myosin.length;                     // half-length of myosin
    vec endstop_force = {0.0, 0.0, 0.0};
    if (std::fabs(s) <= Lm) {
        double dist_to_end = Lm - std::fabs(s);       // ∈ [0, Lm]
        double mag_mid = stiffness * dist_to_end/Lm;         // increase toward center
        //endstop_force = mag_mid * actin.direction[act_idx]; // pushes outward along actin direction
        endstop_force = mag_mid * (s >= 0 ? u : -u);      //force on actin pushes outward along myosin axis
        if (dist_to_end/Lm > 0.1){
}
    }
    return endstop_force; //force on actin
}

bool apply_myosin_repulsion(const Filament& actin,
                            const Myosin& myosin,
                            int i,
                            int j,
                            const std::vector<double>& box,
                            int fix_myosin,
                            const utils::MoleculeConnection& actinIndicesPerMyosin,
                            double stiffness,
                            vec& force_on_first,
                            vec& force_on_second,
                            double& segment_distance,
                            double exclusion_distance)
{
    auto result = compute_myosin_repulsion(
        actin,
        myosin,
        i,
        j,
        box,
        fix_myosin,
        actinIndicesPerMyosin,
        stiffness,
        exclusion_distance);
    segment_distance = result.segment_distance;
    if (!result.applied) {
        return false;
    }
    force_on_first += result.force_on_first;
    force_on_second += result.force_on_second;
    return true;
}

bool apply_actin_repulsion(const Filament& actin,
                           int i,
                           int j,
                           const std::vector<double>& box,
                           double crosslinker_length,
                           double stiffness,
                           vec& force_on_first,
                           vec& force_on_second)
{
    auto result = compute_actin_repulsion(
        actin,
        i,
        j,
        box,
        crosslinker_length,
        stiffness);
    if (!result.applied) {
        return false;
    }
    force_on_first += result.force_on_first;
    force_on_second += result.force_on_second;
    return true;
}

bool apply_actin_myosin_repulsion(const Filament& actin,
                                  const Myosin& myosin,
                                  int act_idx,
                                  int myo_idx,
                                  const std::vector<double>& box,
                                  double radius,
                                  double stiffness,
                                  vec& force_on_actin,
                                  vec& force_on_myosin)
{
    constexpr double EPS_FORCE = 1e-9;
    vec repulsive_force = compute_actin_myosin_repulsion(
        actin,
        myosin,
        act_idx,
        myo_idx,
        box,
        radius,
        stiffness);
    if (repulsive_force.norm() <= EPS_FORCE) {
        return false;
    }
    force_on_actin += repulsive_force;
    force_on_myosin -= repulsive_force;
    return true;
}
