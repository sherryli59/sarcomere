#include "interaction.h"
#include "geometry.h"
#include <cstdio>
#include <cstdlib>  // For exit()
#include <cmath>
#include <limits>
#include <algorithm>

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
    double axial_force_max) {
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

    result.energy = 0.5 * radial_k * extension * extension +
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
        connector_gradient +=
            (radial_k * extension / result.length) * result.connector;
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
                                      double radial_slack) {
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
    const double actin_coordinate = std::clamp(
        (result.myosin_attachment - actin_center_image).dot(axis_a),
        -0.5 * actin.length, 0.5 * actin.length);
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
    const double magnitude = -stiffness * loaded_extension;
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
    bool strong, bool bidirectional_strong) {
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

    const double radial_extension = rho - std::max(0.0, radial_preferred);
    const double slack = std::max(0.0, radial_slack);
    double loaded_radial_extension = 0.0;
    if (radial_extension > slack) {
        loaded_radial_extension = radial_extension - slack;
    } else if (radial_extension < -slack) {
        loaded_radial_extension = radial_extension + slack;
    }
    if (radial_stiffness > 0.0 && rho > 1e-12 &&
        loaded_radial_extension != 0.0) {
        const double dU_drho = radial_stiffness * loaded_radial_extension;
        gradient_h += (dU_drho / rho) * transverse;
        // rho^2 = h.h - (h.u)^2.  This explicit orientation derivative is
        // required in addition to the material-point lever arm.
        explicit_gradient_axis += -(dU_drho * z / rho) * result.connector;
        result.radial_energy = 0.5 * radial_stiffness *
                               loaded_radial_extension * loaded_radial_extension;
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

vec compute_actin_myosin_radial_repulsion(
    const Filament& actin, const Myosin& myosin,
    int actin_index, int myosin_index, const std::vector<double>& box,
    double excluded_centerline_distance, double stiffness) {
    if (!(excluded_centerline_distance > 0.0) || !(stiffness > 0.0)) {
        return {0.0, 0.0, 0.0};
    }
    const auto geometry_result = geometry::segment_segment_distance_w_normal(
        actin.left_end[actin_index], actin.right_end[actin_index],
        myosin.left_end[myosin_index], myosin.right_end[myosin_index],
        box, actin.periodic_axes);
    const double distance = geometry_result.first;
    if (distance >= excluded_centerline_distance) return {0.0, 0.0, 0.0};
    vec normal = geometry_result.second.normal;
    if (normal.norm() <= 1e-12) {
        normal = static_cast<vec>(actin.center[actin_index]) -
                 static_cast<vec>(myosin.center[myosin_index]);
        normal.pbc_wrap(box, actin.periodic_axes);
        vec axis = myosin.direction[myosin_index];
        axis.normalize();
        normal -= axis * normal.dot(axis);
    }
    if (normal.norm() <= 1e-12) normal = {0.0, 1.0, 0.0};
    normal.normalize();
    return stiffness * (excluded_centerline_distance - distance) * normal;
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
                                         double stiffness)
{
    RepulsionResult result{};
    const double cutoff = 2.0 * myosin.radius;
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
        const double denom = std::max(center_displacement.norm(), EPS);
        if (denom <= EPS) {
            return result;
        }
        normal_vector = center_displacement / denom;
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
        return result;
    }
    if (j < fix_myosin) {
        result.force_on_first += force_vec * 2.0;
        return result;
    }

    int status_i = 0;
    const auto& connections_i = actinIndicesPerMyosin.getConnections(i);
    for (int act_idx : connections_i) {
        status_i = std::max(status_i, actin.cb_status[act_idx]);
    }

    int status_j = 0;
    const auto& connections_j = actinIndicesPerMyosin.getConnections(j);
    for (int act_idx : connections_j) {
        status_j = std::max(status_j, actin.cb_status[act_idx]);
    }

    if (status_i > 0 && status_j == 0) {
        result.force_on_second -= force_vec * 2.0;
    } else if (status_i == 0 && status_j > 0) {
        result.force_on_first += force_vec * 2.0;
    } else {
        result.force_on_first += force_vec;
        result.force_on_second -= force_vec;
    }
    // if (overlap > 0.1) {
    //     printf("Myosin-myosin repulsion between myosins %d and %d: distance=%f, overlap=%f, force_mag=%f,force_vec=(%f, %f, %f)\n",
    //            i, j, distance, overlap, force_mag,
    //             force_vec.x, force_vec.y, force_vec.z);
    // }
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
    if (distance >= crosslinker_length) {
        return result;
    }

    vec normal_vector = geom.second.normal;
    double norm = normal_vector.norm();
    if (norm <= EPS) {
        // Parallel coincident rods have a real excluded-volume overlap but no
        // unique distance gradient.  First use the transverse part of their
        // center separation.  If that also vanishes, choose a deterministic
        // pair-specific transverse direction.  The index sign makes swapping
        // (i,j) reverse the regularized normal, preserving pair symmetry.
        vec axis_i = actin.direction[i];
        vec axis_j = actin.direction[j];
        axis_i.normalize();
        axis_j.normalize();
        if (axis_i.dot(axis_j) < 0.0) axis_j = -axis_j;
        vec common_axis = axis_i + axis_j;
        if (common_axis.norm() <= EPS) common_axis = axis_i;
        common_axis.normalize();

        normal_vector =
            center_displacement - common_axis * center_displacement.dot(common_axis);
        double fallback_norm = normal_vector.norm();
        if (fallback_norm <= EPS) {
            const std::uint64_t low =
                static_cast<std::uint64_t>(std::min(i, j) + 1);
            const std::uint64_t high =
                static_cast<std::uint64_t>(std::max(i, j) + 1);
            const std::uint64_t hash =
                (low * 0x9e3779b97f4a7c15ULL) ^
                (high * 0xbf58476d1ce4e5b9ULL);
            const std::array<vec, 3> references{
                vec{1.0, 0.0, 0.0},
                vec{0.0, 1.0, 0.0},
                vec{0.0, 0.0, 1.0}};
            vec reference = references[static_cast<size_t>(hash % 3ULL)];
            if (std::abs(reference.dot(common_axis)) > 0.9) {
                reference = references[
                    static_cast<size_t>((hash + 1ULL) % 3ULL)];
            }
            normal_vector =
                reference - common_axis * reference.dot(common_axis);
            normal_vector.normalize();
            if (i > j) normal_vector = -normal_vector;
        } else {
            normal_vector = normal_vector / fallback_norm;
        }
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
                            double& segment_distance)
{
    auto result = compute_myosin_repulsion(
        actin,
        myosin,
        i,
        j,
        box,
        fix_myosin,
        actinIndicesPerMyosin,
        stiffness);
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
