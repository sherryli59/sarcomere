#include "langevin.h"
#include "sarcomere.h"
#include "components.h"
#include "lattice_spacing.h"
#include "cxxopts.hpp"
#include <array>
#include <vector>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <omp.h>

int main(int argc, char* argv[]){
    // Complete record of the run configuration, filled immediately after
    // option parsing and handed to the model before the HDF5 file is created.
    std::string cli_invocation_capture;
    std::string cli_defaults_capture;
    int nsteps;
    int seed;
    double dt;
    double beta;
    double actin_diff_coeff_trans;
    double actin_diff_coeff_rot;
    double max_actin_displacement;
    double max_myosin_displacement;
    double myosin_diff_coeff_trans;
    double myosin_diff_coeff_rot;
    bool physical_mobility;
    double temperature_kelvin;
    double thermal_noise_scale;
    double solvent_viscosity;
    double actin_hydrodynamic_diameter;
    double myosin_hydrodynamic_diameter;
    double actin_crowding_parallel;
    double actin_crowding_perpendicular;
    double actin_crowding_tumble;
    double actin_crowding_roll;
    double myosin_crowding_parallel;
    double myosin_crowding_perpendicular;
    double myosin_crowding_tumble;
    double myosin_crowding_roll;
    double max_displacement;
    double max_actin_rotation;
    double max_myosin_rotation;
    int save_every;
    double k_on;
    double k_off;
    double base_lifetime;
    double directional_base_lifetime;
    double lifetime_coeff;
    double k_aa;
    double kappa_aa;
    double k_am;
    double k_mm;
    int myomesin_sites_per_filament;
    double myomesin_k_on;
    double myomesin_k_off;
    double myomesin_axial_capture;
    double myomesin_radial_rest;
    double myomesin_radial_capture;
    double myomesin_radial_k;
    double myomesin_alignment_cos;
    double myomesin_bell_distance;
    double myomesin_force_max;
    double kappa_am;
    double v_am;
    int n_actins;
    int n_myosins;
    double Lx;
    double Ly;
    double Lz;
    double actin_length;
    double myosin_length;
    double myosin_radius;
    double actin_direction_noise;
    double myosin_direction_noise;
    double myosin_direction_tilt_deg;
    double actin_direction_polar_fraction;
    double myosin_direction_polar_fraction;
    double am_cutoff;
    double am_optimal;
    double actin_radius;
    double am_surface_gap;
    double am_steric_k;
    double actin_steric_k;
    double myosin_steric_k;
    double myosin_exclusion_distance;
    bool am_stroke_detailed_balance;
    double am_stroke_base_rate;
    double am_stroke_theta;
    double am_stroke_free_energy;
    bool aa_double_overlap_contact;
    double aa_double_overlap_contact_diameter;
    double aa_double_overlap_contact_epsilon;
    double aa_double_overlap_contact_force_cap;
    double six_ring_angle_noise_deg;
    double six_ring_radial_noise;
    double six_ring_axis_noise_deg;
    bool six_ring_prebind_aa;
    bool prebundle_assembly;
    double three_membered_transverse_jitter;
    double three_membered_axial_jitter;
    double three_membered_tilt_jitter;
    double aa_cutoff;
    double aa_optimal;
    double aa_spacing_ratio;
    double aa_preferred_centerline;
    double aa_capture_margin;
    double aa_antiparallel_overlap;
    double init_aa_overlap;
    double minimal_catch_tilt_deg;
    double minimal_catch_myosin_tilt_deg;
    double minimal_catch_myosin_inward_shift;
    double aa_barbed_zone;
    double aa_axial_register_force_max;
    double aa_radial_force_max;
    double am_tether_force_max;
    double am_tether_occupancy_floor;
    double am_tether_occupancy_floor_reach;
    double probe_tilt_deg;
    double probe_confine_k;
    double aa_axial_register_ratio;
    bool aa_register_all_statuses;
    bool aa_perpendicular_all_statuses;
    double aa_status2_support_window;
    double aa_motor_pair_support_window;
    double aa_motor_pair_commit_time;
    double aa_motor_pair_switch_ratio;
    double aa_motor_pair_exchange_time;
    double aa_motor_pair_dropout_time;
    double aa_status2_registry_on;
    double aa_status2_registry_off;
    double aa_register_saturation_error;
    bool am_bidirectional_crossbridge;
    bool am_crossover_pull_through;
    double am_overlap_target_ratio;
    bool aa_local_multilink;
    int aa_prebundle_links_per_pair;
    double aa_site_exclusion;
    double aa_overlap_capture;
    double aa_perpendicular_k;
    double aa_status1_radial_k;
    double aa_status1_site_diffusion;
    bool aa_status2_locked_until_detach;
    double aa_status2_unlock_rate;
    double aa_status2_unlock_force_scale;
    double aa_status2_unlock_registry_scale;
    bool aa_arm_state_model;
    double aa_arm_open_rate_weak;
    double aa_arm_open_rate_mature;
    double aa_arm_bell_distance;
    double aa_tether_off_rate;
    double aa_tether_rebind_rate;
    double aa_tether_capture_distance;
    double aa_arm_rebind_theta;
    double aa_actn2_contact_radius;
    double aa_myosin_contact_epsilon;
    double aa_myosin_contact_force_cap;
    double aa_myosin_contact_shell_width;
    double aa_myosin_contact_shell_energy;
    bool aa_zdisc_tip_contact;
    double aa_zdisc_tip_onset;
    double aa_zdisc_tip_stiffness;
    double aa_zdisc_tip_force_max;
    bool aa_zdisc_tip_contact_requires_link;
    bool aa_zdisc_tip_hard_projection;
    double aa_barbed_lateral_range;
    double aa_barbed_lateral_k;
    double aa_barbed_lateral_force_max;
    double aa_barbed_lateral_axial_window;
    double aa_zdisc_tip_transverse_flat;
    double aa_zdisc_tip_transverse_cutoff;
    double aa_zdisc_tip_stop_stiffness;
    double aa_zdisc_tip_stop_force_max;
    double aa_zdisc_tip_stop_gate_radius;
    double aa_zdisc_tip_stop_gate_cutoff;
    double aa_zdisc_tip_stop_epsilon;
    double aa_occlusion_rupture_penetration;
    double aa_occlusion_capture_radius;
    bool aa_overlap_occlusion_veto;
    double aa_occlusion_rupture_dwell;
    double aa_occlusion_rupture_rate;
    bool aa_barbed_tip_endstop_slab;
    double aa_barbed_tip_endstop_transverse_cutoff;
    double aa_association_max_direction_dot;
    double aa_link_rupture_length;
    int aa_stiffness_reference_links;
    int aa_total_link_pool;
    bool aa_neighbor_identity;
    double aa_neighbor_slot_memory;
    double aa_status2_slot_memory;
    double tau_rec;
    double titin_k;
    double titin_rest_length_option;
    std::string titin_mode_name;
    int titin_links_per_half;
    double titin_capture_distance;
    double titin_k_on;
    double titin_k_off;
    bool titin_strain_gated_kinetics;
    double titin_transition_theta;
    double titin_rebind_memory;
    double titin_rebind_rate_multiplier;
    bool titin_require_live_actn2;
    double titin_k_off_unsupported;
    bool titin_rebind_requires_actn2;
    int fixed_lattice_preequil_steps;
    double k_bundle_max;
    int bundle_ramp_steps;
    bool resume;
    int resume_frame;
    bool directional;
    bool stabilize_status1;
    std::string aa_status_mode_name;
    actn2_load_kinetics::StatusMode parsed_aa_status_mode;
    double aa_status2_force_on;
    double aa_status2_force_off;
    bool dynamic_am;
    bool am_offcenter_torque;
    bool am_spatial_heads;
    bool sterics_enabled;
    int am_heads_per_half;
    int am_spokes_per_crown;
    bool am_head_azimuth;
    bool am_azimuth_gate;
    double am_azimuth_gate_half_width_deg;
    bool torque_breakdown;
    bool legacy_actin_torque_cap;
    double am_head_radial_force_max;
    int am_max_heads_per_actin_half;
    int am_sites_per_actin;
    int am_max_heads_per_actin;
    double am_k_on;
    double am_weak_binding_floor;
    double am_k_weak_off;
    double am_k_powerstroke;
    double am_k_strong_off;
    double am_force_scale;
    double am_motor_k;
    double am_power_stroke;
    double am_stall_force;
    double am_unloaded_velocity;
    double am_weak_stiffness_ratio;
    double am_tether_k;
    double am_tether_slack;
    double am_head_transverse_k;
    double am_alignment_k;
    double am_alignment_slack_deg;
    int am_tether_saturation_heads;
    int am_head_radial_mode;
    double am_head_radial_k_override;
    double am_head_radial_rest_override;
    int am_alignment_saturation_heads;
    int am_motor_saturation_heads;
    bool am_binding_azimuth_selection;
    double am_binding_azimuth_half_width_deg;
    double am_rupture_distance;
    bool am_strain_dependent_detachment;
    double am_strain_detachment_theta;
    double am_transverse_strain_detachment_theta;
    double am_strain_detachment_max_exponent;
    double am_strong_off_floor_fraction;
    double am_cycle_period;
    double am_activation_rise;
    double am_activation_decay;
    double am_activation_floor;
    double am_activation_active_fraction;
    bool am_thin_activation;
    bool am_thin_feedback;
    double am_q_calcium_open_rate;
    double am_q_calcium_half;
    double am_q_calcium_hill;
    double am_q_strong_feedback_rate;
    double am_q_strong_half;
    double am_q_strong_hill;
    double am_q_close_rate;
    bool am_initialize_q_calcium_equilibrium;
    bool am_thick_feedback;
    double am_thick_recruit_rate;
    double am_thick_fold_rate;
    double am_thick_tension_scale;
    double am_thick_initial_on_fraction;
    double am_thick_force_clamp;
    bool actin_turnover;
    double actin_turnover_rate;
    double actin_turnover_start_time;
    int actin_turnover_max_placement_attempts;
    double actin_turnover_clearance_ratio;
    double am_status2_memory;
    bool deterministic;
    bool isometric_actin_axial;
    int n_fixed_myosins;
    int max_myosin_bonds;
    int max_aa_bonds_per_actin;
    bool aa_energy_gated_association;
    double aa_transition_state_theta;
    bool aa_passive_detailed_balance;
    double aa_encounter_k;
    double aa_encounter_kappa;
    int dimension;
    int pre_phase_steps;
    int interaction_ramp_steps;
    std::array<bool,3> periodic_axes{true, true, true};

    auto parse_pbc_mask = [](const std::string& mask) -> std::array<bool,3> {
        if (mask.size() != 3) {
            throw std::invalid_argument("pbc_mask must include exactly three characters (e.g., 110).");
        }
        std::array<bool,3> axes{};
        for (size_t idx = 0; idx < 3; ++idx) {
            char c = mask[idx];
            if (c == '1' || c == 'T' || c == 't' || c == 'y' || c == 'Y') {
                axes[idx] = true;
            } else if (c == '0' || c == 'F' || c == 'f' || c == 'n' || c == 'N') {
                axes[idx] = false;
            } else {
                throw std::invalid_argument("Invalid character in pbc_mask. Use only 0/1 (or T/F).");
            }
        }
        return axes;
    };

    std::string filename;
    std::string init_struc;
    std::string initial_geometry_file;
    int initial_geometry_frame;
    bool beta_was_explicit = false;
    try {
        cxxopts::Options options("sarcomere", "simulate sarcomere assembly using Monte Carlo");

        options.add_options()
            ("nsteps", "Number of steps", cxxopts::value<int>(nsteps)->default_value("500000"))
            ("seed", "Seed value", cxxopts::value<int>(seed)->default_value("0"))
            ("dt", "Time step", cxxopts::value<double>(dt)->default_value("0.00001"))
            ("beta", "[UNUSED 2026-09-17: escape hatch, kept; superseded by --temperature_kelvin] Beta value", cxxopts::value<double>(beta)->default_value("241.0"))
            ("actin_diff_coeff_trans", "Actin translational diffusion coefficient", cxxopts::value<double>(actin_diff_coeff_trans)->default_value("0.2"))
            ("actin_diff_coeff_rot", "Actin rotational diffusion coefficient", cxxopts::value<double>(actin_diff_coeff_rot)->default_value("0"))
            ("myosin_diff_coeff_trans", "Myosin translational diffusion coefficient", cxxopts::value<double>(myosin_diff_coeff_trans)->default_value("0.05"))
            ("myosin_diff_coeff_rot", "Myosin rotational diffusion coefficient", cxxopts::value<double>(myosin_diff_coeff_rot)->default_value("0.05"))
            ("physical_mobility", "Use temperature-derived anisotropic rigid-rod mobility (production default; pass false for legacy scalar diffusion)", cxxopts::value<bool>(physical_mobility)->default_value("true")->implicit_value("true"))
            ("temperature_kelvin", "Absolute temperature for thermal mobility", cxxopts::value<double>(temperature_kelvin)->default_value("310"))
            ("thermal_noise_scale", "Multiplier on Brownian translation and rotation; zero gives deterministic mechanics", cxxopts::value<double>(thermal_noise_scale)->default_value("1"))
            ("solvent_viscosity", "Reference solvent viscosity in Pa s", cxxopts::value<double>(solvent_viscosity)->default_value("0.001"))
            ("actin_hydrodynamic_diameter", "Thin-filament hydrodynamic diameter (um)", cxxopts::value<double>(actin_hydrodynamic_diameter)->default_value("0.010"))
            ("myosin_hydrodynamic_diameter", "Thick-filament hydrodynamic diameter (um)", cxxopts::value<double>(myosin_hydrodynamic_diameter)->default_value("0.035"))
            ("actin_crowding_parallel", "Actin parallel drag multiplier", cxxopts::value<double>(actin_crowding_parallel)->default_value("10"))
            ("actin_crowding_perpendicular", "Actin perpendicular drag multiplier", cxxopts::value<double>(actin_crowding_perpendicular)->default_value("10"))
            ("actin_crowding_tumble", "Actin tumbling drag multiplier", cxxopts::value<double>(actin_crowding_tumble)->default_value("160"))
            ("actin_crowding_roll", "Actin axial-roll drag multiplier", cxxopts::value<double>(actin_crowding_roll)->default_value("160"))
            ("myosin_crowding_parallel", "Myosin parallel drag multiplier", cxxopts::value<double>(myosin_crowding_parallel)->default_value("20"))
            ("myosin_crowding_perpendicular", "Myosin perpendicular drag multiplier", cxxopts::value<double>(myosin_crowding_perpendicular)->default_value("20"))
            ("myosin_crowding_tumble", "Myosin tumbling drag multiplier", cxxopts::value<double>(myosin_crowding_tumble)->default_value("64"))
            ("myosin_crowding_roll", "Myosin axial-roll drag multiplier", cxxopts::value<double>(myosin_crowding_roll)->default_value("64"))
            ("max_actin_displacement", "Maximum deterministic displacement per step from forces for actin",
             cxxopts::value<double>(max_actin_displacement)->default_value("0.005"))
            ("max_myosin_displacement", "Maximum deterministic displacement per step from forces for myosin",
             cxxopts::value<double>(max_myosin_displacement)->default_value("0.005"))
            ("max_actin_rotation", "Maximum rotational magnitude per step for actin",
             cxxopts::value<double>(max_actin_rotation)->default_value("0.005"))
            ("max_myosin_rotation", "Maximum rotational magnitude per step for myosin",
             cxxopts::value<double>(max_myosin_rotation)->default_value("0.005"))
            ("save_every", "Save every", cxxopts::value<int>(save_every)->default_value("200"))
            ("k_on", "k_on", cxxopts::value<double>(k_on)->default_value("1000"))
            ("k_off", "k_off", cxxopts::value<double>(k_off)->default_value("1"))
            ("base_lifetime", "Load-independent lifetime for status-1 and status-2 bonds", cxxopts::value<double>(base_lifetime)->default_value("0.1"))
            ("directional_base_lifetime", "Additional constant lifetime for status-2 bonds", cxxopts::value<double>(directional_base_lifetime)->default_value("0.05"))
            ("lifetime_coeff", "Load-dependent lifetime coefficient", cxxopts::value<double>(lifetime_coeff)->default_value("0.4"))
            ("stabilize_status1", "Also apply lifetime_coeff*f_load to status-1 bonds (original kinetics)", cxxopts::value<bool>(stabilize_status1)->default_value("false")->implicit_value("true"))
            ("aa_status_mode", "AA status classifier: motor_transmitted_register uses a persistent, source-resolved myosin pair and the actual ACTN2 connector reaction; productive_register retains the aggregate-vector rule", cxxopts::value<std::string>(aa_status_mode_name)->default_value("productive_register"))
            ("aa_status2_force_on", "Load required for status 1->2", cxxopts::value<double>(aa_status2_force_on)->default_value("1.0"))
            ("aa_status2_force_off", "Load below which status 2->1", cxxopts::value<double>(aa_status2_force_off)->default_value("0.5"))
            ("k_aa", "k_aa", cxxopts::value<double>(k_aa)->default_value("10"))
            ("kappa_aa", "ACTN2 angular stiffness (pN*um) in U = 0.5*kappa_aa*(1-cos^2 theta); divided by aa_stiffness_reference_links for the per-link value. Default lowered 100 -> 3 on 2026-09-18: 100 pN*um = 23,364 kT/rad^2 froze the Z-pair at 0.53 deg rms and a 33x cut costs zero delivered rotational drive", cxxopts::value<double>(kappa_aa)->default_value("3"))
            ("k_am", "k_am", cxxopts::value<double>(k_am)->default_value("10"))
            ("k_mm", "Effective axial stiffness (pN/um) of one occupied sparse coarse-grained myomesin module; 0 disables all modules", cxxopts::value<double>(k_mm)->default_value("0"))
            ("myomesin_sites_per_filament", "Finite local-module edge valency of each thick filament", cxxopts::value<int>(myomesin_sites_per_filament)->default_value("2"))
            ("myomesin_k_on", "Local coarse-grained myomesin-module association rate (1/s)", cxxopts::value<double>(myomesin_k_on)->default_value("5000"))
            ("myomesin_k_off", "Zero-load coarse-grained myomesin-module dissociation rate (1/s)", cxxopts::value<double>(myomesin_k_off)->default_value("1"))
            ("myomesin_axial_capture", "Maximum bare-zone axial mismatch for new module capture (um)", cxxopts::value<double>(myomesin_axial_capture)->default_value("0.25"))
            ("myomesin_radial_rest", "Nonzero transverse rest spacing of a myomesin module (um)", cxxopts::value<double>(myomesin_radial_rest)->default_value("0.11"))
            ("myomesin_radial_capture", "Allowed deviation from myomesin radial rest spacing at capture (um)", cxxopts::value<double>(myomesin_radial_capture)->default_value("0.04"))
            ("myomesin_radial_k", "Independent transverse stiffness of one coarse-grained module (pN/um); 0 uses k_mm", cxxopts::value<double>(myomesin_radial_k)->default_value("300"))
            ("myomesin_alignment_cos", "Minimum absolute thick-axis dot product for new bridge capture", cxxopts::value<double>(myomesin_alignment_cos)->default_value("0.94"))
            ("myomesin_bell_distance", "Bell-model transition distance for load-accelerated module release (um); 0 gives load-independent turnover", cxxopts::value<double>(myomesin_bell_distance)->default_value("0.0005"))
            ("myomesin_force_max", "Maximum resultant force carried by one coarse-grained myomesin module (pN); scale with represented molecular occupancy; 0 disables the cap", cxxopts::value<double>(myomesin_force_max)->default_value("40"))
            ("kappa_am", "kappa_am", cxxopts::value<double>(kappa_am)->default_value("100"))
            ("v_am", "[UNUSED 2026-09-17: not set by any of 177 sbatch files] v_am", cxxopts::value<double>(v_am)->default_value("10"))
            ("dynamic_am", "Use persistent detached/weak/strong actin-myosin kinetics", cxxopts::value<bool>(dynamic_am)->default_value("true")->implicit_value("true"))
            ("am_offcenter_torque", "Apply explicit AM forces at head locations to rotate filaments", cxxopts::value<bool>(am_offcenter_torque)->default_value("true")->implicit_value("true"))
            ("am_spatial_heads", "Use azimuthally resolved thick-filament head anchors", cxxopts::value<bool>(am_spatial_heads)->default_value("true")->implicit_value("true"))
            ("sterics", "Enable filament steric repulsion", cxxopts::value<bool>(sterics_enabled)->default_value("true")->implicit_value("true"))
            ("am_spokes_per_crown", "Head sites per axial crown; 3 is the physiological layout (default), 1 gives one head per axial station and so 3x the crowns at fixed am_heads_per_half", cxxopts::value<int>(am_spokes_per_crown)->default_value("3"))
            ("legacy_actin_torque_cap", "Restore the pre-2026-09-18 actin torque clamp, which budgeted the FULL gradient magnitude rather than the component that can rotate the filament. Actin has no roll channel, so the axial component is inert; budgeting it made the clamp fire on 98.8% of entries in the canonical wrap rig and discard 50.2% of the real rotational drive while the effective torque sat 4.7x below the cap. Set this ONLY to reproduce runs made before the fix", cxxopts::value<bool>(legacy_actin_torque_cap)->default_value("false")->implicit_value("true"))
            ("torque_breakdown", "Record per-source ACTIN TORQUE contributions as tb_<source>_{x,y,z} plus tb_total_precap. Off by default and bitwise inert. The ACTN2 link torque is split by which stiffness supplies it: lever, aa_perpendicular_k, kappa_aa (align) and the axial register", cxxopts::value<bool>(torque_breakdown)->default_value("false")->implicit_value("true"))
            ("am_head_azimuth", "Radial channel mode. false (default) anchors the head on the myosin AXIS and measures its radial spring against am_preferred_centerline, so the azimuth cannot strain a head by geometry alone and the aggregate pair tether is folded per-head. true restores the legacy surface-anchored behaviour with the separate pair tether", cxxopts::value<bool>(am_head_azimuth)->default_value("false"))
            ("am_azimuth_gate", "Azimuth-selection gate: a head may only bind an actin that its FIXED crown spoke points at, within am_azimuth_gate_half_width_deg.  A binding permission, not a geometry change, so it composes with am_head_azimuth=false.  Default false is bitwise inert", cxxopts::value<bool>(am_azimuth_gate)->default_value("false")->implicit_value("true"))
            ("am_azimuth_gate_half_width_deg", "Half-width of the azimuth gate in degrees.  <=0 selects 180/am_spokes_per_crown, which partitions the circle so exactly one spoke per crown admits any actin.  >=180 admits everything", cxxopts::value<double>(am_azimuth_gate_half_width_deg)->default_value("0"))
            ("am_head_radial_force_max", "Saturating cap on the per-head radial force (pN). Negative (default) DERIVES it from the integrator displacement clamp as max_actin_displacement/(mu_perp*dt*am_tether_saturation_heads), so a full complement of heads exactly reaches the clamp and the cap tracks dt. 0 keeps the unbounded linear spring; >0 is taken literally", cxxopts::value<double>(am_head_radial_force_max)->default_value("-1"))
            ("am_heads_per_half", "Conserved myosin head slots on each thick-filament half (three sites per crown; 141 gives ~14.3-nm crown spacing at the default length)", cxxopts::value<int>(am_heads_per_half)->default_value("141"))
            ("am_max_heads_per_actin_half", "Full-overlap weak+strong attachment capacity for one actin on one myosin half; scales with productive overlap, strong heads have priority, 0 disables", cxxopts::value<int>(am_max_heads_per_actin_half)->default_value("24"))
            ("am_sites_per_actin", "Material actin binding-site lattice shared by all neighboring myosins (~5.5-nm spacing at 182 sites/um)", cxxopts::value<int>(am_sites_per_actin)->default_value("182"))
            ("am_max_heads_per_actin", "Global simultaneous weak+strong head occupancy across every myosin touching one actin", cxxopts::value<int>(am_max_heads_per_actin)->default_value("24"))
            ("am_k_on", "Weak AM attachment rate (1/time)", cxxopts::value<double>(am_k_on)->default_value("100"))
            ("am_weak_binding_floor", "Basal weak AM attachment fraction at q=0 (weak binding only)", cxxopts::value<double>(am_weak_binding_floor)->default_value("0.1"))
            ("am_k_weak_off", "Weak AM detachment rate (1/time)", cxxopts::value<double>(am_k_weak_off)->default_value("200"))
            ("am_k_powerstroke", "Zero-load weak-to-strong rate (1/time)", cxxopts::value<double>(am_k_powerstroke)->default_value("200"))
            ("am_k_strong_off", "Zero-load strong AM detachment rate (1/time)", cxxopts::value<double>(am_k_strong_off)->default_value("160"))
            ("am_force_scale", "Force scale for load-slowed AM rates", cxxopts::value<double>(am_force_scale)->default_value("1"))
            ("am_motor_k", "Axial motor spring constant (pN/um)", cxxopts::value<double>(am_motor_k)->default_value("1000"))
            ("am_power_stroke", "Axial power-stroke distance", cxxopts::value<double>(am_power_stroke)->default_value("0.005"))
            ("am_stall_force", "Maximum productive force per strong AM cross-bridge (pN)", cxxopts::value<double>(am_stall_force)->default_value("4.5"))
            ("am_unloaded_velocity", "Zero-load speed of the finite AM working stroke (um/s)", cxxopts::value<double>(am_unloaded_velocity)->default_value("5"))
            ("am_weak_stiffness_ratio", "Weak passive axial stiffness / strong motor stiffness", cxxopts::value<double>(am_weak_stiffness_ratio)->default_value("0.1"))
            ("am_tether_k", "Pair-level aggregate flat-bottom AM containment stiffness outside the slack band", cxxopts::value<double>(am_tether_k)->default_value("3000"))
            ("am_tether_slack", "Half-width of the force-free AM tether distance band", cxxopts::value<double>(am_tether_slack)->default_value("0.002"))
            ("am_head_transverse_k", "Passive transverse stiffness of each bound explicit AM head (pN/um; zero disables)", cxxopts::value<double>(am_head_transverse_k)->default_value("0"))
            ("am_alignment_k", "Weak occupied-pair AM nematic alignment stiffness (pN um/rad^2 near alignment)", cxxopts::value<double>(am_alignment_k)->default_value("30"))
            ("am_alignment_slack_deg", "Force-free AM axial misalignment half-width in degrees", cxxopts::value<double>(am_alignment_slack_deg)->default_value("2"))
            ("am_tether_saturation_heads", "Occupied heads needed to reach the pair-level AM tether stiffness", cxxopts::value<int>(am_tether_saturation_heads)->default_value("6"))
            ("am_head_radial_k", "Override the AM head radial stiffness (pN/um). The composite default is am_head_transverse_k + am_tether_k/am_tether_saturation_heads = 600, which with the derived cap gives a saturation length lambda = f_max/k of only 3.2 nm -- far below the tens of nm of radial difference a tilted actin produces, so the two-head aligning differential saturates to zero. Lowering k alone raises lambda at UNCHANGED peak pull. <=0 (default) keeps the composite and is bitwise inert.", cxxopts::value<double>(am_head_radial_k_override)->default_value("0"))
            ("am_head_radial_mode", "What the AM head radial spring penalises. 0 = line (DEFAULT, shipped): the perpendicular from the anchor to the ACTIN LINE, identical for every point on that line, so the bound site effectively slides and a tilt lands entirely on it. 1 = connector: the full length |h| from the bound material site. 2 = axis: distance from the bound material site to the MYOSIN AXIS, keeping radial and axial separate so the power stroke is untouched. Standalone two-filament escape barrier from 36 deg: 9.9 kT for mode 0, 0.9 for mode 1, 0.0 for mode 2.", cxxopts::value<int>(am_head_radial_mode)->default_value("0"))
            ("am_head_radial_rest", "Override the AM head radial rest length (um). <=0 (default) inherits am_preferred_centerline and is bitwise inert.", cxxopts::value<double>(am_head_radial_rest_override)->default_value("0"))
            ("am_alignment_saturation_heads", "Occupancy count for the AM ALIGNING TORQUE alone. am_alignment_k is multiplied by min(1, bound heads on the pair / this), so with the default it inherits am_tether_saturation_heads=6 and the delivered aligning constant scales with head number -- MEASURED 0.6-0.7 at hph27 vs 1.2 at hph54 against a nominal 3, which confounds every head-count comparison. Set to 1 to deliver the full nominal am_alignment_k to any engaged pair, head-count independent. 0 (default) inherits the tether count and is bitwise inert.", cxxopts::value<int>(am_alignment_saturation_heads)->default_value("0"))
            ("am_motor_saturation_heads", "Axial counterpart of am_tether_saturation_heads: scales each head's axial stiffness AND stall cap by 1/min(bound heads on that pair, this value), so a pair's TOTAL sliding force saturates instead of growing with head count. 0 (default) disables it and is bitwise inert", cxxopts::value<int>(am_motor_saturation_heads)->default_value("0"))
            ("am_binding_azimuth_selection", "Let a head swivel toward reachable actin at binding", cxxopts::value<bool>(am_binding_azimuth_selection)->default_value("true")->implicit_value("true"))
            ("am_binding_azimuth_half_width_deg", "Maximum binding-time swivel away from the fixed crown azimuth", cxxopts::value<double>(am_binding_azimuth_half_width_deg)->default_value("30"))
            ("am_rupture_distance", "Emergency AM contour extension beyond the relaxed connector", cxxopts::value<double>(am_rupture_distance)->default_value("0.1"))
            ("am_strain_dependent_detachment", "Accelerate AM detachment continuously with nonproductive elastic strain", cxxopts::value<bool>(am_strain_dependent_detachment)->default_value("true")->implicit_value("true"))
            ("am_strain_detachment_theta", "Fraction of negative strong-head axial strain energy coupled to the detachment barrier", cxxopts::value<double>(am_strain_detachment_theta)->default_value("0.1"))
            ("am_transverse_strain_detachment_theta", "Fraction of transverse bound-head spring energy coupled to the detachment barrier", cxxopts::value<double>(am_transverse_strain_detachment_theta)->default_value("0"))
            ("am_strain_detachment_max_exponent", "Numerical cap on the strain-accelerated detachment exponent", cxxopts::value<double>(am_strain_detachment_max_exponent)->default_value("20"))
            ("am_strong_off_floor_fraction", "Minimum strong-head detachment rate as a fraction of the zero-load rate", cxxopts::value<double>(am_strong_off_floor_fraction)->default_value("1"))
            ("am_cycle_period", "Heartbeat period; <=0 gives tonic activation", cxxopts::value<double>(am_cycle_period)->default_value("0"))
            ("am_activation_rise", "Activation transient rise time", cxxopts::value<double>(am_activation_rise)->default_value("0.02"))
            ("am_activation_decay", "Activation transient decay time", cxxopts::value<double>(am_activation_decay)->default_value("0.2"))
            ("am_activation_floor", "Minimum thin-filament activation [0,1]", cxxopts::value<double>(am_activation_floor)->default_value("0"))
            ("am_activation_active_fraction", "Fraction of each heartbeat occupied by the calcium transient; the remaining phase is exactly at am_activation_floor", cxxopts::value<double>(am_activation_active_fraction)->default_value("0.5"))
            ("am_thin_activation", "Evolve per-actin thin-filament activation q_i", cxxopts::value<bool>(am_thin_activation)->default_value("true")->implicit_value("true"))
            ("am_thin_feedback", "Allow strong heads to increase q_i", cxxopts::value<bool>(am_thin_feedback)->default_value("true")->implicit_value("true"))
            ("am_q_calcium_open_rate", "Calcium-driven thin opening rate", cxxopts::value<double>(am_q_calcium_open_rate)->default_value("20"))
            ("am_q_calcium_half", "Half-activation calcium level", cxxopts::value<double>(am_q_calcium_half)->default_value("0.5"))
            ("am_q_calcium_hill", "Calcium Hill exponent", cxxopts::value<double>(am_q_calcium_hill)->default_value("3"))
            ("am_q_strong_feedback_rate", "Strong-head thin-feedback opening rate", cxxopts::value<double>(am_q_strong_feedback_rate)->default_value("100"))
            ("am_q_strong_half", "Strong-head fraction at half feedback", cxxopts::value<double>(am_q_strong_half)->default_value("0.02"))
            ("am_q_strong_hill", "Strong-head feedback Hill exponent", cxxopts::value<double>(am_q_strong_hill)->default_value("2"))
            ("am_q_close_rate", "Thin-filament closing rate", cxxopts::value<double>(am_q_close_rate)->default_value("20"))
            ("am_initialize_q_calcium_equilibrium", "Initialize q_i at the calcium-only stationary state while leaving AM heads, ACTN2, and titin unbound", cxxopts::value<bool>(am_initialize_q_calcium_equilibrium)->default_value("true")->implicit_value("true"))
            ("am_thick_feedback", "Enable tension-dependent OFF/SRX-to-ON recruitment", cxxopts::value<bool>(am_thick_feedback)->default_value("false")->implicit_value("true"))
            ("am_thick_recruit_rate", "Zero-tension OFF/SRX-to-ON rate", cxxopts::value<double>(am_thick_recruit_rate)->default_value("2"))
            ("am_thick_fold_rate", "Zero-tension detached-ON-to-OFF/SRX rate", cxxopts::value<double>(am_thick_fold_rate)->default_value("20"))
            ("am_thick_tension_scale", "Half-filament tension scale for thick recruitment", cxxopts::value<double>(am_thick_tension_scale)->default_value("3"))
            ("am_thick_initial_on_fraction", "Initial detached-ON fraction per thick-filament half", cxxopts::value<double>(am_thick_initial_on_fraction)->default_value("0.1"))
            ("am_thick_force_clamp", "[UNUSED 2026-09-17: not set by any of 177 sbatch files] Diagnostic half-filament tension clamp; negative uses mechanical tension", cxxopts::value<double>(am_thick_force_clamp)->default_value("-1"))
            ("actin_turnover", "Enable fixed-population whole-actin renewal", cxxopts::value<bool>(actin_turnover)->default_value("false")->implicit_value("true"))
            ("actin_turnover_rate", "Per-actin Poisson renewal rate (1/s)", cxxopts::value<double>(actin_turnover_rate)->default_value("0"))
            ("actin_turnover_start_time", "Simulation time before actin renewal begins (s)", cxxopts::value<double>(actin_turnover_start_time)->default_value("0"))
            ("actin_turnover_max_placement_attempts", "Maximum sterically rejected placement proposals per renewal event", cxxopts::value<int>(actin_turnover_max_placement_attempts)->default_value("5000"))
            ("actin_turnover_clearance_ratio", "Fraction of the nominal hard-core contact distance a renewed actin must clear; the soft force-capped sterics let the equilibrium population sit inside the nominal distance, so 1.0 rejects nearly every proposal", cxxopts::value<double>(actin_turnover_clearance_ratio)->default_value("0.8"))
            ("am_status2_memory", "Time that an established status-2 AA bond retains status 2 after mechanical support is lost", cxxopts::value<double>(am_status2_memory)->default_value("0"))
            ("n_actins", "Number of actins", cxxopts::value<int>(n_actins)->default_value("50"))
            ("n_myosins", "Number of myosins", cxxopts::value<int>(n_myosins)->default_value("4"))
            ("Lx", "Lx", cxxopts::value<double>(Lx)->default_value("10"))
            ("Ly", "Ly", cxxopts::value<double>(Ly)->default_value("10"))
            ("Lz", "Lz", cxxopts::value<double>(Lz)->default_value("10"))
            ("pbc_mask", "Periodicity mask as three digits (e.g., 110 => periodic in x,y)", cxxopts::value<std::string>()->default_value("111"))
            ("periodic_mask", "Deprecated alias for --pbc_mask", cxxopts::value<std::string>())
            ("actin_length", "Actin length", cxxopts::value<double>(actin_length)->default_value("1"))
            ("myosin_length", "Myosin length", cxxopts::value<double>(myosin_length)->default_value("1.5"))
            ("myosin_radius", "Myosin radius", cxxopts::value<double>(myosin_radius)->default_value("0.025"))
            ("actin_direction_noise", "Std dev for aligning actin directions to +x; negative leaves isotropic initialization", cxxopts::value<double>(actin_direction_noise)->default_value("-1"))
            ("myosin_direction_noise", "Std dev for aligning myosin directions to +x during initialization", cxxopts::value<double>(myosin_direction_noise)->default_value("1"))
            ("actin_direction_polar_fraction", "Fraction of actins seeded along +x rather than -x when actin_direction_noise >= 0. 1 (default) is the historical POLAR start and is bitwise identical to it. Nematic order S is sign-blind, so 0.5 gives the SAME order with mixed polarity -- which is what ACTN2 needs, since it only links pairs with direction_dot <= aa_association_max_direction_dot. A fully polar start has ~0 antiparallel pairs and cannot build a Z-disc", cxxopts::value<double>(actin_direction_polar_fraction)->default_value("1"))
            ("myosin_direction_polar_fraction", "Same, for myosin", cxxopts::value<double>(myosin_direction_polar_fraction)->default_value("1"))
            ("myosin_direction_tilt_deg", "Mean in-plane tilt of initialized myosin directions from +x", cxxopts::value<double>(myosin_direction_tilt_deg)->default_value("0"))
            ("am_cutoff", "cutoff for am interaction range", cxxopts::value<double>(am_cutoff)->default_value("0.05"))
            ("am_optimal", "optimal distance for am interaction", cxxopts::value<double>(am_optimal)->default_value("0.03"))
            ("actin_radius", "Effective thin-filament radius", cxxopts::value<double>(actin_radius)->default_value("0.015"))
            ("am_surface_gap", "Preferred surface gap between actin and myosin", cxxopts::value<double>(am_surface_gap)->default_value("0.015"))
            ("am_steric_k", "Actin-myosin radial hard-core stiffness", cxxopts::value<double>(am_steric_k)->default_value("1000"))
            ("actin_steric_k", "Actin-actin hard-core penalty independent of AA mechanics", cxxopts::value<double>(actin_steric_k)->default_value("10000"))
            ("myosin_steric_k", "Myosin-myosin hard-core stiffness independent of AA mechanics", cxxopts::value<double>(myosin_steric_k)->default_value("1000"))
            ("am_stroke_detailed_balance", "[UNUSED 2026-09-17: alternative detailed-balance stroke model, never exercised] Discrete working stroke with strain-dependent forward/reverse rates obeying detailed balance, instead of the prescribed force-velocity ramp", cxxopts::value<bool>(am_stroke_detailed_balance)->default_value("false")->implicit_value("true"))
            ("am_stroke_base_rate", "[UNUSED 2026-09-17: alternative detailed-balance stroke model, never exercised] Attempt rate for the discrete working stroke (1/s)", cxxopts::value<double>(am_stroke_base_rate)->default_value("1000"))
            ("am_stroke_theta", "[UNUSED 2026-09-17: alternative detailed-balance stroke model, never exercised] Bell split of the stroke barrier (0..1)", cxxopts::value<double>(am_stroke_theta)->default_value("0.5"))
            ("am_stroke_free_energy", "[UNUSED 2026-09-17: alternative detailed-balance stroke model, never exercised] Chemical free energy driving one stroke (pN.um); <0 derives it from am_stall_force", cxxopts::value<double>(am_stroke_free_energy)->default_value("-1"))
            ("myosin_exclusion_distance", "Range of the myosin-myosin repulsion (0 = 2*myosin_radius hard core only; ~0.0953 holds the thick-filament lattice at the spacing that seats an actin in the trigonal interstice)", cxxopts::value<double>(myosin_exclusion_distance)->default_value("0"))
            ("aa_double_overlap_contact", "Enable local antiparallel pointed-end unilateral contact (sandbox; default off)", cxxopts::value<bool>(aa_double_overlap_contact)->default_value("false")->implicit_value("true"))
            ("aa_double_overlap_contact_diameter", "Decorated pointed-end contact diameter (um)", cxxopts::value<double>(aa_double_overlap_contact_diameter)->default_value("0.012"))
            ("aa_double_overlap_contact_epsilon", "Pointed-end WCA energy scale (pN um)", cxxopts::value<double>(aa_double_overlap_contact_epsilon)->default_value("0.00428"))
            ("aa_double_overlap_contact_force_cap", "Optional pointed-end compression/load-sharing force ceiling (pN; <=0 uncapped)", cxxopts::value<double>(aa_double_overlap_contact_force_cap)->default_value("0"))
            ("aa_spacing_ratio", "Preferred AA/AM centerline-spacing ratio, used only when aa_preferred_centerline <= 0. Default 0.92 -> 0.9875 (= 79/80) on 2026-09-18 so the derived path agrees with the direct default rather than restoring the 73.6 nm inconsistency", cxxopts::value<double>(aa_spacing_ratio)->default_value("0.9875"))
            ("aa_preferred_centerline", "Direct AA centerline target (um); <=0 derives it from aa_spacing_ratio. Default changed -1 (which derived 0.92*80 = 73.6 nm) -> 0.079 on 2026-09-18: aa_barbed_lateral_range is a one-sided wall at 80 nm whose 2000 pN/um is 6x the ACTN2 per-link 333 pN/um, so it parks the pair at 79 nm and a 73.6 nm rest leaves every link 5.4 nm strained", cxxopts::value<double>(aa_preferred_centerline)->default_value("0.079"))
            ("aa_capture_margin", "AA capture distance beyond its preferred centerline", cxxopts::value<double>(aa_capture_margin)->default_value("0.02"))
            ("aa_antiparallel_overlap", "Preferred axial overlap of antiparallel actin barbed-end regions", cxxopts::value<double>(aa_antiparallel_overlap)->default_value("0.1"))
            ("init_aa_overlap", "[UNUSED 2026-09-17: not set by any of 177 sbatch files] Initial Z-disc interdigitation for two_sarcomere_unit; <=0 starts perfectly registered at aa_antiparallel_overlap", cxxopts::value<double>(init_aa_overlap)->default_value("-1"))
            ("minimal_catch_tilt_deg", "Initial thin-filament misalignment for minimal_catch_* rigs (degrees), applied as a rigid rotation of each prebound ACTN2 pair about its barbed-end midpoint. Moves the thin filaments out of motor reach; prefer --minimal_catch_myosin_tilt_deg", cxxopts::value<double>(minimal_catch_tilt_deg)->default_value("0"))
            ("minimal_catch_myosin_tilt_deg", "Initial THICK-filament misalignment for minimal_catch_* rigs (degrees), about the actin-myosin offset direction through the filament centre, which preserves the capture separation exactly and never touches the Z-disc", cxxopts::value<double>(minimal_catch_myosin_tilt_deg)->default_value("0"))
            ("minimal_catch_myosin_inward_shift", "Axial inward shift (um) of BOTH thick filaments toward the central Z-disc, for minimal_catch_balanced* AND two_sarcomere_unit. Default -1e6 = use the variant value (0 for balanced, 0.55 for balanced_buried), so this flag is bitwise inert unless set. Clearance of the central barbed tip is 0.1975-shift um at Lx=4.", cxxopts::value<double>(minimal_catch_myosin_inward_shift)->default_value("-1e6"))
            ("aa_barbed_zone", "Z-disc width: ACTN2 binds only within this distance (um) of a barbed end; <=0 binds anywhere in the mutual overlap", cxxopts::value<double>(aa_barbed_zone)->default_value("-1"))
            ("aa_axial_register_force_max", "Ceiling (pN) for a saturating status-2 axial register force that does not decay far from the target; <=0 keeps the Gaussian-attenuated form", cxxopts::value<double>(aa_axial_register_force_max)->default_value("36"))
            ("aa_radial_force_max", "[UNUSED 2026-09-17: not set by any of 177 sbatch files] Saturating cap on the ACTN2 radial force (pN); <=0 keeps the unbounded quadratic. Small-displacement stiffness is unchanged; bounding the force is what permits a larger dt without a displacement clamp", cxxopts::value<double>(aa_radial_force_max)->default_value("0"))
            ("am_tether_force_max", "[UNUSED 2026-09-17: not set by any of 177 sbatch files] Saturating cap on the AM head tether force (pN); <=0 keeps the unbounded quadratic", cxxopts::value<double>(am_tether_force_max)->default_value("0"))
            ("am_tether_occupancy_floor", "Floor on the radial pair tether occupancy multiplier; 0 keeps the bound-head gate", cxxopts::value<double>(am_tether_occupancy_floor)->default_value("0"))
            ("am_tether_occupancy_floor_reach", "Reach past am_preferred_centerline within which the occupancy floor applies (um)", cxxopts::value<double>(am_tether_occupancy_floor_reach)->default_value("0.04"))
            ("probe_tilt_deg", "Initial actin-myosin misalignment for single_am_tilted (degrees)", cxxopts::value<double>(probe_tilt_deg)->default_value("20"))
            ("probe_confine_k", "Harmonic trap (pN/um) on each actin centre toward its initial position; torque-free; 0 disables", cxxopts::value<double>(probe_confine_k)->default_value("0"))
            ("aa_axial_register_ratio", "Status-2 axial register stiffness as a fraction of radial AA stiffness", cxxopts::value<double>(aa_axial_register_ratio)->default_value("1.0"))
            ("aa_register_all_statuses", "Let every attached local ACTN2 link carry axial register mechanics while status 1/2 still controls turnover lifetime", cxxopts::value<bool>(aa_register_all_statuses)->default_value("false")->implicit_value("true"))
            ("aa_perpendicular_all_statuses", "Whether status-1 ACTN2 links also carry the perpendicular (axial-sliding) constraint. Default true preserves historical behaviour; set false to make status 1 a sliding bond of radial + alignment only, so crosslinked filaments can still slide to search for register. NOTE: aa_status_mode=arm_load_productive gates status-2 promotion on exactly this load, so false makes promotion impossible in that mode", cxxopts::value<bool>(aa_perpendicular_all_statuses)->default_value("true")->implicit_value("true"))
            ("aa_status2_registry_on", "Maximum barbed-end register error for status 1 to status 2 capture", cxxopts::value<double>(aa_status2_registry_on)->default_value("0.1"))
            ("aa_status2_support_window", "Motor-timescale window (s) over which productive opposed load is averaged for status-2 FORMATION only; <=0 requires instantaneous bilateral loading", cxxopts::value<double>(aa_status2_support_window)->default_value("0"))
            ("aa_motor_pair_support_window", "EMA window (s) for one selected motor pair's ACTN2-transmitted support in motor_transmitted_register", cxxopts::value<double>(aa_motor_pair_support_window)->default_value("0.005"))
            ("aa_motor_pair_commit_time", "Minimum continuous candidate-identity residence time (s) before status-2 maturation", cxxopts::value<double>(aa_motor_pair_commit_time)->default_value("0.010"))
            ("aa_motor_pair_switch_ratio", "Alternative/current transmitted-support ratio required to exchange the weak-state motor pair (>=1)", cxxopts::value<double>(aa_motor_pair_switch_ratio)->default_value("1.5"))
            ("aa_motor_pair_exchange_time", "Minimum continuous time (s) that the same stronger challenger must remain best before weak-state partner exchange", cxxopts::value<double>(aa_motor_pair_exchange_time)->default_value("0.001"))
            ("aa_motor_pair_dropout_time", "Maximum time (s) a weak selected motor pair may provide zero transmitted support before its identity and commitment clock are cleared", cxxopts::value<double>(aa_motor_pair_dropout_time)->default_value("0.005"))
            ("aa_status2_registry_off", "Maximum barbed-end register error for status-2 retention", cxxopts::value<double>(aa_status2_registry_off)->default_value("0.15"))
            ("aa_register_saturation_error", "Error scale that bounds the status-2 register energy and force", cxxopts::value<double>(aa_register_saturation_error)->default_value("0.05"))
            ("am_bidirectional_crossbridge", "Strong heads act as bidirectional elastic springs (restoring force on overshoot) instead of tension-only", cxxopts::value<bool>(am_bidirectional_crossbridge)->default_value("false")->implicit_value("true"))
            ("aa_zdisc_tip_hard_projection", "Run the hard position projection that displaces a penetrating ACTN2 patch back out of a myosin; false leaves penetration to be resolved by actn2_myosin_collision instead", cxxopts::value<bool>(aa_zdisc_tip_hard_projection)->default_value("true")->implicit_value("true"))
            ("aa_zdisc_tip_contact_requires_link", "Apply the Z-disc barbed-tip endstop only to actins carrying a live ACTN2 link; a bare barbed tip then interdigitates freely. Default true (2026-10-25): the capsule excludes the Z-disc, not a free thin filament", cxxopts::value<bool>(aa_zdisc_tip_contact_requires_link)->default_value("true")->implicit_value("true"))
            ("am_crossover_pull_through", "Once an actin slides past the bare zone into the opposite head-bearing half, let that half power stroke too; direction follows actin polarity so it keeps pulling the same way until the actin exits", cxxopts::value<bool>(am_crossover_pull_through)->default_value("true")->implicit_value("true"))
            ("am_overlap_target_ratio", "Legacy sliding throttle: engaged-overlap fraction at which motor sliding stalls (<=0 disables)", cxxopts::value<double>(am_overlap_target_ratio)->default_value("0"))
            ("aa_local_multilink", "Use sparse local multi-link ACTN2 mechanics (production default; pass false for the legacy pair edge)", cxxopts::value<bool>(aa_local_multilink)->default_value("true")->implicit_value("true"))
            ("aa_prebundle_links_per_pair", "Maximum transient local ACTN2 links per actin pair", cxxopts::value<int>(aa_prebundle_links_per_pair)->default_value("1"))
            ("aa_site_exclusion", "Minimum contour-coordinate separation between local ACTN2 attachments", cxxopts::value<double>(aa_site_exclusion)->default_value("0.0192"))
            ("aa_overlap_capture", "Reference overlap used to calibrate overlap-extensive local ACTN2 association", cxxopts::value<double>(aa_overlap_capture)->default_value("0.1"))
            ("aa_perpendicular_k", "Local ACTN2 perpendicularity stiffness", cxxopts::value<double>(aa_perpendicular_k)->default_value("50"))
            ("aa_status1_radial_k", "Per-molecule radial stiffness (pN/um) for weak/status-1 ACTN2 links; <=0 uses the locked per-link stiffness", cxxopts::value<double>(aa_status1_radial_k)->default_value("-1"))
            ("aa_status1_site_diffusion", "One-dimensional diffusion coefficient (um^2/s) of each status-1 ACTN2 binding domain along its actin barbed-end zone", cxxopts::value<double>(aa_status1_site_diffusion)->default_value("0.010"))
            ("aa_status2_locked_until_detach", "Keep stereospecific status-2 ACTN2 material sites locked until molecular detachment or contour rupture", cxxopts::value<bool>(aa_status2_locked_until_detach)->default_value("false")->implicit_value("true"))
            ("aa_status2_unlock_rate", "[UNUSED 2026-09-17: status-2 unlock path, never exercised] Zero-support mature-to-weak ACTN2 conformational rate (1/s); 0 disables this transition", cxxopts::value<double>(aa_status2_unlock_rate)->default_value("0"))
            ("aa_status2_unlock_force_scale", "[UNUSED 2026-09-17: status-2 unlock path, never exercised] Motor-support scale (pN) suppressing mature-to-weak ACTN2 relaxation", cxxopts::value<double>(aa_status2_unlock_force_scale)->default_value("0.1"))
            ("aa_status2_unlock_registry_scale", "[UNUSED 2026-09-17: status-2 unlock path, never exercised] Registry-error scale (um) accelerating mature-to-weak ACTN2 relaxation outside the retention band", cxxopts::value<double>(aa_status2_unlock_registry_scale)->default_value("0.02"))
            ("aa_arm_state_model", "Use explicit two-arm ACTN2 opening, single-arm tethering, and local rebinding; disables continuous status-1 site diffusion and permanent status-2 locking", cxxopts::value<bool>(aa_arm_state_model)->default_value("false")->implicit_value("true"))
            ("aa_arm_open_rate_weak", "[UNUSED 2026-09-17: two-arm ACTN2 state model, disabled by aa_arm_state_model=false in all 32 uses] Zero-force opening rate (1/s) of each arm in a weak two-arm ACTN2 bridge", cxxopts::value<double>(aa_arm_open_rate_weak)->default_value("20"))
            ("aa_arm_open_rate_mature", "[UNUSED 2026-09-17: two-arm ACTN2 state model, disabled by aa_arm_state_model=false in all 32 uses] Zero-force opening rate (1/s) of each arm in a mature two-arm ACTN2 bridge", cxxopts::value<double>(aa_arm_open_rate_mature)->default_value("1"))
            ("aa_arm_bell_distance", "[UNUSED 2026-09-17: two-arm ACTN2 state model, disabled by aa_arm_state_model=false in all 32 uses] Bell distance (um) for tensile/contact-assisted opening of one ACTN2 arm", cxxopts::value<double>(aa_arm_bell_distance)->default_value("0.001"))
            ("aa_tether_off_rate", "[UNUSED 2026-09-17: superseded ACTN2 tether model] Zero-force rate (1/s) for loss of the remaining arm of a singly tethered ACTN2", cxxopts::value<double>(aa_tether_off_rate)->default_value("5"))
            ("aa_tether_rebind_rate", "[UNUSED 2026-09-17: superseded ACTN2 tether model] Attempt rate (1/s) for a tethered ACTN2 free arm to close onto a local exposed actin site", cxxopts::value<double>(aa_tether_rebind_rate)->default_value("100"))
            ("aa_tether_capture_distance", "[UNUSED 2026-09-17: superseded ACTN2 tether model] Maximum connector length (um) considered for tethered ACTN2 free-arm rebinding", cxxopts::value<double>(aa_tether_capture_distance)->default_value("0.10"))
            ("aa_arm_rebind_theta", "[UNUSED 2026-09-17: two-arm ACTN2 state model, disabled by aa_arm_state_model=false in all 32 uses] Transition-state fraction multiplying candidate bridge energy in the tethered-arm rebinding rate", cxxopts::value<double>(aa_arm_rebind_theta)->default_value("0.5"))
            ("aa_actn2_contact_radius", "Effective ACTN2 connector radius (um) used for finite contact with a thick filament", cxxopts::value<double>(aa_actn2_contact_radius)->default_value("0.005"))
            ("aa_myosin_contact_epsilon", "WCA energy scale (pN um) for ACTN2-connector/thick-filament contact", cxxopts::value<double>(aa_myosin_contact_epsilon)->default_value("0.00428"))
            ("aa_myosin_contact_force_cap", "Coarse-grained force cap (pN) for one ACTN2/thick-filament contact; <=0 is uncapped", cxxopts::value<double>(aa_myosin_contact_force_cap)->default_value("40"))
            ("aa_myosin_contact_shell_width", "Width (um) of the optional soft ACTN2/thick-filament crowded-head corona; 0 disables it and preserves the original WCA contact", cxxopts::value<double>(aa_myosin_contact_shell_width)->default_value("0.015"))
            ("aa_myosin_contact_shell_energy", "Corona energy (pN um) at the physical hard-core surface; 0 disables the soft corona", cxxopts::value<double>(aa_myosin_contact_shell_energy)->default_value("0.0214"))
            ("aa_zdisc_tip_contact", "Enable head-independent axial exclusion and nonpenetration between myosin tips and bilateral ACTN2/Z-disc patches", cxxopts::value<bool>(aa_zdisc_tip_contact)->default_value("true")->implicit_value("true"))
            ("aa_zdisc_tip_onset", "Axial myosin-tip/Z-patch gap (um) where the end-stop begins", cxxopts::value<double>(aa_zdisc_tip_onset)->default_value("0.050"))
            ("aa_zdisc_tip_stiffness", "Small-compression axial stiffness (pN/um) of the barbed-tip end-stop; the ACTN2 connector force uses aa_zdisc_tip_stop_stiffness", cxxopts::value<double>(aa_zdisc_tip_stiffness)->default_value("1000"))
            ("aa_zdisc_tip_force_max", "Saturating force (pN) of the barbed-tip end-stop; the ACTN2 connector force uses aa_zdisc_tip_stop_force_max", cxxopts::value<double>(aa_zdisc_tip_force_max)->default_value("20"))
            ("aa_zdisc_tip_transverse_flat", "Transverse radius (um) of full strength for the geometric tip-gap queries (exclusion audit, hard projection)", cxxopts::value<double>(aa_zdisc_tip_transverse_flat)->default_value("0.100"))
            ("aa_zdisc_tip_transverse_cutoff", "Transverse cutoff (um) for the geometric tip-gap queries (exclusion audit, hard projection)", cxxopts::value<double>(aa_zdisc_tip_transverse_cutoff)->default_value("0.130"))
            ("aa_zdisc_tip_stop_stiffness", "Axial-stop connector force: stiffness (pN/um) of the quadratic part of Phi", cxxopts::value<double>(aa_zdisc_tip_stop_stiffness)->default_value("1000"))
            ("aa_zdisc_tip_stop_force_max", "Axial-stop connector force: saturation slope (pN) of Phi, the per-connector axial force bound", cxxopts::value<double>(aa_zdisc_tip_stop_force_max)->default_value("30"))
            ("aa_zdisc_tip_stop_gate_radius", "Axial-stop connector force: radius (um) inside which the radial gate is 1", cxxopts::value<double>(aa_zdisc_tip_stop_gate_radius)->default_value("0.130"))
            ("aa_zdisc_tip_stop_gate_cutoff", "Axial-stop connector force: radius (um) where the smoothstep gate reaches 0", cxxopts::value<double>(aa_zdisc_tip_stop_gate_cutoff)->default_value("0.140"))
            ("aa_zdisc_tip_stop_epsilon", "Axial-stop connector force: axial smoothing length (um) in sqrt(z^2 + epsilon^2)", cxxopts::value<double>(aa_zdisc_tip_stop_epsilon)->default_value("0.005"))
            ("aa_barbed_lateral_range", "Target transverse barbed-end spacing (um) enforced between ALL barbed-end pairs; 0 disables. ACTN2 alone contracts a crosslinked cluster to the actin hard core, so this supplies the lateral pressure a Z-disc face needs to tile", cxxopts::value<double>(aa_barbed_lateral_range)->default_value("0"))
            ("aa_barbed_lateral_k", "Stiffness (pN/um) of the barbed-end lateral repulsion; 0 disables", cxxopts::value<double>(aa_barbed_lateral_k)->default_value("0"))
            ("aa_barbed_lateral_force_max", "Saturating force (pN) for one barbed-end lateral contact; <=0 is uncapped", cxxopts::value<double>(aa_barbed_lateral_force_max)->default_value("0"))
            ("aa_barbed_lateral_axial_window", "Axial gate (um): barbed ends separated by more than this along the filament axis belong to different Z-discs and do not repel", cxxopts::value<double>(aa_barbed_lateral_axial_window)->default_value("0.250"))
            ("aa_occlusion_rupture_penetration", "Deep ACTN2/thick-filament penetration (um) causing immediate geometric rupture", cxxopts::value<double>(aa_occlusion_rupture_penetration)->default_value("0.010"))
            ("aa_overlap_occlusion_veto", "Enable reason-3 ACTN2 rupture only when either linked actin barbed tip penetrates the myosin capture capsule by aa_occlusion_rupture_penetration; AM reach/contour overlap alone never ruptures", cxxopts::value<bool>(aa_overlap_occlusion_veto)->default_value("true")->implicit_value("true"))
            ("aa_occlusion_rupture_dwell", "Grace period (SECONDS) that barbed-tip burial must persist before the reason-3 ACTN2 deletion fires. 0 (default) keeps the original probability-1 immediate deletion and is bitwise inert. Converted internally to steps via dt, so the physics does not change with dt. Measured burial episodes last a median 20 us / max 320 us in a registered lattice, so 1e-4 to 5e-4 spares 92-100% of transients while still punishing sustained burial. The per-link counter resets on resume.", cxxopts::value<double>(aa_occlusion_rupture_dwell)->default_value("0"))
            ("aa_occlusion_rupture_rate", "Zero-depth hazard rate (1/s) for reason-3 barbed-tip burial. 0 (default) keeps the original probability-1 deletion and is bitwise inert. Nonzero makes burial a RATE like every other ACTN2 death channel: k(p) = rate * exp((p - aa_occlusion_rupture_penetration)/aa_occlusion_rupture_penetration), summed over the two barbed tips. A soft repulsion cannot defend a probability-1 threshold, because a spring in a thermal bath holds a distribution of depths and the first excursion past the line is fatal; a rate charges a transient dip a small probability while still killing sustained deep burial.", cxxopts::value<double>(aa_occlusion_rupture_rate)->default_value("0"))
            ("aa_barbed_tip_endstop_slab", "Treat the A-band barbed-tip exclusion as a slab rather than a transversely gated capsule end cap. Drops the hard cutoff at aa_occlusion_capture_radius (which in production equals the resting A-M centerline, so every bound actin sits on it and a tip arriving from outside it meets no force at all until it is already past the kill depth) and drops the spherocylinder cap extension (which makes the axial boundary depend on the transverse distance and so gives this purely axial force a nonzero curl). false (default) is bitwise inert.", cxxopts::value<bool>(aa_barbed_tip_endstop_slab)->default_value("false")->implicit_value("true"))
            ("aa_barbed_tip_endstop_transverse_cutoff", "Transverse locality bound (um) for the slab-mode barbed-tip endstop. 0 (default) falls back to aa_occlusion_capture_radius, i.e. the legacy gate. Only meaningful with aa_barbed_tip_endstop_slab.", cxxopts::value<double>(aa_barbed_tip_endstop_transverse_cutoff)->default_value("0"))
            ("aa_occlusion_capture_radius", "Steric centerline radius (um) used by ACTN2 contour exclusion, the barbed-tip axial end-stop, and reason-3 penetration; default 0.080 excludes the separate 0.050-um motor-head reach, while <=0 restores the legacy full A-M capture centerline", cxxopts::value<double>(aa_occlusion_capture_radius)->default_value("0.080"))
            ("aa_association_max_direction_dot", "Maximum actin director dot product for a new local ACTN2 association; >=1 disables", cxxopts::value<double>(aa_association_max_direction_dot)->default_value("-0.5"))
            ("aa_link_rupture_length", "Finite explicit ACTN2 contour length in um; <=0 disables forced contour rupture", cxxopts::value<double>(aa_link_rupture_length)->default_value("0"))
            ("aa_stiffness_reference_links", "Reference explicit-link count used to convert legacy k_aa and kappa_aa pair stiffnesses into fixed per-ACTN2 stiffnesses", cxxopts::value<int>(aa_stiffness_reference_links)->default_value("3"))
            ("aa_total_link_pool", "Maximum total local ACTN2 links; <=0 uses a reservoir", cxxopts::value<int>(aa_total_link_pool)->default_value("0"))
            ("aa_neighbor_identity", "Retain AA lattice-neighbor slots across molecular turnover", cxxopts::value<bool>(aa_neighbor_identity)->default_value("true")->implicit_value("true"))
            ("aa_neighbor_slot_memory", "Status-1 neighbor-slot memory after the last molecular link detaches (s)", cxxopts::value<double>(aa_neighbor_slot_memory)->default_value("0.02"))
            ("aa_status2_slot_memory", "Status-2-protected neighbor-slot memory after the last molecular link detaches (s)", cxxopts::value<double>(aa_status2_slot_memory)->default_value("0.20"))
            ("six_ring_angle_noise_deg", "Initial angular disorder for six_actin_ring", cxxopts::value<double>(six_ring_angle_noise_deg)->default_value("20"))
            ("six_ring_radial_noise", "Initial radial disorder for six_actin_ring", cxxopts::value<double>(six_ring_radial_noise)->default_value("0.01"))
            ("six_ring_axis_noise_deg", "[UNUSED 2026-09-17: not set by any of 177 sbatch files] Initial actin-axis tilt disorder for six_actin_ring", cxxopts::value<double>(six_ring_axis_noise_deg)->default_value("0"))
            ("six_ring_prebind_aa", "Initialize the six intended nearest-neighbor AA edges before KMC", cxxopts::value<bool>(six_ring_prebind_aa)->default_value("false")->implicit_value("true"))
            ("prebundle_assembly", "[UNUSED 2026-09-17: not set by any of 177 sbatch files] Seed AA links on initially-overlapping actin pairs (both polarities) at t=0", cxxopts::value<bool>(prebundle_assembly)->default_value("false")->implicit_value("true"))
            ("three_membered_transverse_jitter", "[UNUSED 2026-09-17: three-membered init structure, never built] Gaussian transverse position jitter for the unassembled three-membered diagnostic", cxxopts::value<double>(three_membered_transverse_jitter)->default_value("0"))
            ("three_membered_axial_jitter", "[UNUSED 2026-09-17: three-membered init structure, never built] Gaussian axial position jitter for the unassembled three-membered diagnostic", cxxopts::value<double>(three_membered_axial_jitter)->default_value("0"))
            ("three_membered_tilt_jitter", "[UNUSED 2026-09-17: three-membered init structure, never built] Gaussian y/z director jitter for the unassembled three-membered diagnostic", cxxopts::value<double>(three_membered_tilt_jitter)->default_value("0"))
            // ("aa_cutoff", "cutoff for aa interaction range", cxxopts::value<double>(aa_cutoff)->default_value("0.05"))
            // ("aa_optimal", "optimal distance for aa interaction", cxxopts::value<double>(aa_optimal)->default_value("0.03"))
            ("tau_rec", "AA reassociation refractory time after KMC detachment", cxxopts::value<double>(tau_rec)->default_value("0.05"))
            ("titin_k", "Titin effective spring constant", cxxopts::value<double>(titin_k)->default_value("0.0"))
            ("titin_rest_length", "Titin spring rest length (um); <=0 uses the historical length-derived value", cxxopts::value<double>(titin_rest_length_option)->default_value("-1"))
            ("titin_mode", "Titin topology: legacy_crosslinked, balanced_owning, or explicit_dynamic", cxxopts::value<std::string>(titin_mode_name)->default_value("legacy_crosslinked"))
            ("titin_links_per_half", "Explicit titin material domains per myosin half", cxxopts::value<int>(titin_links_per_half)->default_value("3"))
            ("titin_capture_distance", "Maximum distance (um) from a thick-tip titin domain to an actin barbed-end domain", cxxopts::value<double>(titin_capture_distance)->default_value("0.70"))
            ("titin_k_on", "Zero-strain explicit titin association rate (1/s)", cxxopts::value<double>(titin_k_on)->default_value("5000"))
            ("titin_k_off", "Zero-strain explicit titin dissociation rate (1/s)", cxxopts::value<double>(titin_k_off)->default_value("0.5"))
            ("titin_strain_gated_kinetics", "Apply elastic-strain gating to explicit titin dissociation; association is distance-limited and does not prepay spring energy", cxxopts::value<bool>(titin_strain_gated_kinetics)->default_value("true")->implicit_value("true"))
            ("titin_transition_theta", "Transition-state fraction [0,1] controlling the explicit titin strain-dependent off-rate", cxxopts::value<double>(titin_transition_theta)->default_value("0.5"))
            ("titin_rebind_memory", "Time (s) during which a detached explicit titin domain may locally recapture the same actin", cxxopts::value<double>(titin_rebind_memory)->default_value("0.020"))
            ("titin_rebind_rate_multiplier", "Association-rate multiplier for same-domain recapture during titin_rebind_memory", cxxopts::value<double>(titin_rebind_rate_multiplier)->default_value("100"))
            ("titin_k_off_unsupported", "Explicit titin dissociation rate (1/s) used while the bound actin carries no live ACTN2 link. Models the serial actin-ACTN2-titin anchor, whose composite lifetime is set by its weakest element. <=0 keeps a single global k_off", cxxopts::value<double>(titin_k_off_unsupported)->default_value("0"))
            ("titin_rebind_requires_actn2", "Apply the live-ACTN2 requirement to same-partner recapture as well as to new associations. Required for the serial anchor to have any effect: otherwise every unsupported release is immediately undone by an exempt, rate-multiplied recapture", cxxopts::value<bool>(titin_rebind_requires_actn2)->default_value("false")->implicit_value("true"))
            ("titin_require_live_actn2", "Require a live ACTN2 A-A link (any status) on an actin barbed end before a NEW explicit titin association may form; already bound titins are never released by ACTN2 turnover and same-partner recapture inside titin_rebind_memory is exempt", cxxopts::value<bool>(titin_require_live_actn2)->default_value("false")->implicit_value("true"))
            ("fixed_lattice_preequil_steps", "Advance chemistry and binding kinetics for this many steps at fixed filament geometry before release", cxxopts::value<int>(fixed_lattice_preequil_steps)->default_value("0"))
            ("k_bundle_max", "Maximum transverse bundling spring constant", cxxopts::value<double>(k_bundle_max)->default_value("0.0"))
            ("bundle_ramp_steps", "Number of steps used to ramp bundling strength", cxxopts::value<int>(bundle_ramp_steps)->default_value("0"))
            ("resume", "Resume", cxxopts::value<bool>(resume)->default_value("false"))
            ("resume_frame", "[UNUSED 2026-09-17: escape hatch, kept] Frame index to load when resuming (0-based; default loads latest)", cxxopts::value<int>(resume_frame)->default_value("-1"))
            ("directional", "Directional", cxxopts::value<bool>(directional)->default_value("true"))
            ("deterministic", "[UNUSED 2026-09-17: escape hatch, kept] Force deterministic single-threaded execution", cxxopts::value<bool>(deterministic)->default_value("false")->implicit_value("true"))
            ("isometric_actin_axial", "Constrain actin displacement along its own axis (diagnostic isometric mechanics)", cxxopts::value<bool>(isometric_actin_axial)->default_value("false")->implicit_value("true"))
            ("n_fixed_myosins", "Number of fixed myosins", cxxopts::value<int>(n_fixed_myosins)->default_value("0"))
            ("filename", "Filename", cxxopts::value<std::string>(filename)->default_value("data/traj.h5"))
            ("initial_structure", "Type of initial structure", cxxopts::value<std::string>(init_struc)->default_value("random"))
            ("initial_geometry_file", "[UNUSED 2026-09-17: escape hatch, kept] HDF5 trajectory supplying initial filament centers and directions", cxxopts::value<std::string>(initial_geometry_file)->default_value(""))
            ("initial_geometry_frame", "[UNUSED 2026-09-17: escape hatch, kept] Zero-based source frame for geometry import; -1 selects the latest committed frame", cxxopts::value<int>(initial_geometry_frame)->default_value("-1"))
            ("pre_phase_steps", "Number of sterics-only pre-phase steps (used by myosin_x_noise and default init paths)", cxxopts::value<int>(pre_phase_steps)->default_value("500"))
            ("interaction_ramp_steps", "Steps used to smooth-ramp mature AM/AA stiffness after pre-relaxation", cxxopts::value<int>(interaction_ramp_steps)->default_value("0"))
            ("max_myosin_bonds", "Maximum actin bonds per myosin",cxxopts::value<int>(max_myosin_bonds)->default_value("10"))
            ("max_aa_bonds_per_actin", "Maximum retained AA neighbor slots per actin; <=0 is unlimited", cxxopts::value<int>(max_aa_bonds_per_actin)->default_value("4"))
            ("aa_energy_gated_association", "Apply transition-state elastic-energy gating to AA association", cxxopts::value<bool>(aa_energy_gated_association)->default_value("false")->implicit_value("true"))
            ("aa_transition_state_theta", "Fraction [0,1] of AA elastic cost placed in the association barrier (provisional)", cxxopts::value<double>(aa_transition_state_theta)->default_value("0.1"))
            ("aa_passive_detailed_balance", "Passive control: apply complementary elastic factor to AA off-rate", cxxopts::value<bool>(aa_passive_detailed_balance)->default_value("false")->implicit_value("true"))
            ("aa_encounter_k", "Soft radial stiffness used only for AA association gating", cxxopts::value<double>(aa_encounter_k)->default_value("10"))
            ("aa_encounter_kappa", "Angular stiffness used only for AA association gating", cxxopts::value<double>(aa_encounter_kappa)->default_value("0"))
            ("dimension", "Simulation dimensionality (2 or 3)", cxxopts::value<int>(dimension)->default_value("3"))
            ("h, help", "Print usage");

        auto result = options.parse(argc, argv);

        // Capture the full configuration here: `options` does not outlive this
        // try block, and 84 of 259 options have no attribute of their own.
        for (int i = 0; i < argc; ++i) {
            if (i > 0) cli_invocation_capture += ' ';
            cli_invocation_capture += argv[i];
        }
        for (const auto& group : options.groups()) {
            for (const auto& option : options.group_help(group).options) {
                if (option.l.empty()) continue;
                cli_defaults_capture += option.l.front();
                cli_defaults_capture += '=';
                cli_defaults_capture += option.has_default
                    ? option.default_value : std::string("<required>");
                cli_defaults_capture += ';';
            }
        }
        beta_was_explicit = result.count("beta") > 0;
        if (!result.unmatched().empty()) {
            std::string message = "Unrecognized positional/option tokens:";
            for (const auto& token : result.unmatched()) message += " " + token;
            message += ". For boolean values use --flag=true or --flag=false.";
            throw std::invalid_argument(message);
        }

        if (result.count("help")) {
            std::cout << options.help() << std::endl;
            return 0;
        }

        std::string mask_value = result["pbc_mask"].as<std::string>();
        if (result.count("periodic_mask")) {
            mask_value = result["periodic_mask"].as<std::string>();
        }
        periodic_axes = parse_pbc_mask(mask_value);

        if (deterministic) {
            omp_set_dynamic(0);
            omp_set_num_threads(1);
            omp_set_schedule(omp_sched_static, 0);
            std::cout << "Deterministic mode enabled: forcing single-threaded static scheduling.\n";
        } else {
            omp_set_schedule(omp_sched_dynamic, 1);
        }
    } catch (const std::exception& e)
	  {
		    std::cerr << "Error parsing options: " << e.what() << std::endl;
		    return 1;
	  }

    try {
        parsed_aa_status_mode =
            actn2_load_kinetics::parse_status_mode(aa_status_mode_name);
    } catch (const std::exception& e) {
        std::cerr << "Invalid aa_status_mode: " << e.what() << std::endl;
        return 1;
    }


    if (dimension != 2 && dimension != 3) {
        std::cerr << "dimension must be 2 or 3" << std::endl;
        return 1;
    }
    bool is3D = (dimension == 3);
    if (!(thermal_noise_scale >= 0.0) || !std::isfinite(thermal_noise_scale)) {
        std::cerr << "thermal_noise_scale must be finite and non-negative" << std::endl;
        return 1;
    }
    myosin_direction_noise = std::max(0.0, myosin_direction_noise);
    if (!std::isfinite(actin_direction_noise) ||
        !std::isfinite(myosin_direction_tilt_deg)) {
        std::cerr << "Orientation initialization parameters must be finite" << std::endl;
        return 1;
    }
    pre_phase_steps = std::max(0, pre_phase_steps);
    interaction_ramp_steps = std::max(0, interaction_ramp_steps);
    if (initial_geometry_frame < -1) {
        std::cerr << "initial_geometry_frame must be -1 or non-negative"
                  << std::endl;
        return 1;
    }
    if (init_struc == "geometry_h5") {
        if (resume) {
            std::cerr << "geometry_h5 is a fresh molecular initialization, not "
                         "a checkpoint resume; set resume=false"
                      << std::endl;
            return 1;
        }
        if (initial_geometry_file.empty()) {
            std::cerr << "initial_structure=geometry_h5 requires "
                         "initial_geometry_file"
                      << std::endl;
            return 1;
        }
    } else if (!initial_geometry_file.empty()) {
        std::cerr << "initial_geometry_file is only used with "
                     "initial_structure=geometry_h5"
                  << std::endl;
        return 1;
    }
    if (!(myosin_exclusion_distance >= 0.0) ||
        !std::isfinite(myosin_exclusion_distance)) {
        std::cerr << "myosin_exclusion_distance must be finite and non-negative"
                  << std::endl;
        return 1;
    }
    if (k_mm < 0.0 || !std::isfinite(k_mm) ||
        myomesin_sites_per_filament < 1 ||
        myomesin_sites_per_filament > 16 ||
        myomesin_k_on < 0.0 || myomesin_k_off < 0.0 ||
        !(myomesin_axial_capture > 0.0) ||
        !(myomesin_radial_rest > 0.0) ||
        myomesin_radial_capture < 0.0 || myomesin_radial_k < 0.0 ||
        myomesin_alignment_cos < 0.0 || myomesin_alignment_cos > 1.0 ||
        myomesin_bell_distance < 0.0 || myomesin_force_max < 0.0 ||
        !std::isfinite(myomesin_k_on) ||
        !std::isfinite(myomesin_k_off) ||
        !std::isfinite(myomesin_axial_capture) ||
        !std::isfinite(myomesin_radial_rest) ||
        !std::isfinite(myomesin_radial_capture) ||
        !std::isfinite(myomesin_radial_k) ||
        !std::isfinite(myomesin_alignment_cos) ||
        !std::isfinite(myomesin_bell_distance) ||
        !std::isfinite(myomesin_force_max)) {
        std::cerr << "Invalid sparse myomesin module parameters" << std::endl;
        return 1;
    }
    if (!(actin_steric_k >= 0.0) || !(myosin_steric_k >= 0.0) ||
        !(aa_encounter_k >= 0.0) || !(aa_encounter_kappa >= 0.0) ||
        !std::isfinite(actin_steric_k) || !std::isfinite(myosin_steric_k) ||
        !std::isfinite(aa_encounter_k) || !std::isfinite(aa_encounter_kappa)) {
        std::cerr << "Steric and AA encounter stiffnesses must be finite and non-negative" << std::endl;
        return 1;
    }
    if (!(aa_double_overlap_contact_diameter >= 0.0) ||
        !(aa_double_overlap_contact_epsilon >= 0.0) ||
        !(aa_double_overlap_contact_force_cap >= 0.0) ||
        !std::isfinite(aa_double_overlap_contact_diameter) ||
        !std::isfinite(aa_double_overlap_contact_epsilon) ||
        !std::isfinite(aa_double_overlap_contact_force_cap)) {
        std::cerr << "Double-overlap contact parameters must be finite and non-negative"
                  << std::endl;
        return 1;
    }
    if (!(myosin_radius >= 0.0) || !(actin_radius >= 0.0) ||
        !(am_surface_gap >= 0.0) || !std::isfinite(myosin_radius) ||
        !std::isfinite(actin_radius) || !std::isfinite(am_surface_gap)) {
        std::cerr << "Filament radii and am_surface_gap must be finite and non-negative" << std::endl;
        return 1;
    }
    if (!(am_motor_k >= 0.0) || !std::isfinite(am_motor_k) ||
        !(am_stall_force >= 0.0) || !std::isfinite(am_stall_force) ||
        !(am_power_stroke >= 0.0) || !std::isfinite(am_power_stroke) ||
        !(am_unloaded_velocity >= 0.0) ||
        !std::isfinite(am_unloaded_velocity)) {
        std::cerr << "AM motor stiffness, stall force, power stroke, and "
                     "unloaded velocity must "
                     "be finite and non-negative" << std::endl;
        return 1;
    }
    if (am_max_heads_per_actin_half < 0 ||
        am_max_heads_per_actin_half > 256) {
        std::cerr << "am_max_heads_per_actin_half must be between 0 and 256"
                  << std::endl;
        return 1;
    }
    if (am_sites_per_actin < 1 || am_sites_per_actin > 512) {
        std::cerr << "am_sites_per_actin must be between 1 and 512"
                  << std::endl;
        return 1;
    }
    if (am_max_heads_per_actin < 1 || am_max_heads_per_actin > 512) {
        std::cerr << "am_max_heads_per_actin must be between 1 and 512"
                  << std::endl;
        return 1;
    }
    if (!(aa_spacing_ratio > 0.0) || !std::isfinite(aa_spacing_ratio)) {
        std::cerr << "aa_spacing_ratio must be finite and positive" << std::endl;
        return 1;
    }
    if (!(aa_capture_margin > 0.0) || !std::isfinite(aa_capture_margin)) {
        std::cerr << "aa_capture_margin must be finite and positive" << std::endl;
        return 1;
    }
    if (aa_prebundle_links_per_pair < 1 || aa_prebundle_links_per_pair > 16 ||
        aa_site_exclusion < 0.0 || !std::isfinite(aa_site_exclusion) ||
        aa_overlap_capture <= 0.0 || !std::isfinite(aa_overlap_capture) ||
        aa_status2_registry_on < 0.0 ||
        !std::isfinite(aa_status2_registry_on) ||
        aa_status2_registry_off < aa_status2_registry_on ||
        !std::isfinite(aa_status2_registry_off) ||
        aa_register_saturation_error <= 0.0 ||
        !std::isfinite(aa_register_saturation_error) ||
        aa_perpendicular_k < 0.0 || !std::isfinite(aa_perpendicular_k) ||
        !std::isfinite(aa_status1_radial_k) || aa_status1_radial_k == 0.0 ||
        aa_status1_site_diffusion < 0.0 ||
        !std::isfinite(aa_status1_site_diffusion) ||
        aa_status2_unlock_rate < 0.0 ||
        !std::isfinite(aa_status2_unlock_rate) ||
        aa_status2_unlock_force_scale <= 0.0 ||
        !std::isfinite(aa_status2_unlock_force_scale) ||
        aa_status2_unlock_registry_scale <= 0.0 ||
        !std::isfinite(aa_status2_unlock_registry_scale) ||
        aa_arm_open_rate_weak < 0.0 || !std::isfinite(aa_arm_open_rate_weak) ||
        aa_arm_open_rate_mature < 0.0 || !std::isfinite(aa_arm_open_rate_mature) ||
        aa_arm_bell_distance < 0.0 || !std::isfinite(aa_arm_bell_distance) ||
        aa_tether_off_rate < 0.0 || !std::isfinite(aa_tether_off_rate) ||
        aa_tether_rebind_rate < 0.0 || !std::isfinite(aa_tether_rebind_rate) ||
        aa_tether_capture_distance <= 0.0 || !std::isfinite(aa_tether_capture_distance) ||
        aa_arm_rebind_theta < 0.0 || aa_arm_rebind_theta > 1.0 ||
        !std::isfinite(aa_arm_rebind_theta) ||
        aa_actn2_contact_radius < 0.0 || !std::isfinite(aa_actn2_contact_radius) ||
        aa_myosin_contact_epsilon < 0.0 || !std::isfinite(aa_myosin_contact_epsilon) ||
        !std::isfinite(aa_myosin_contact_force_cap) ||
        aa_myosin_contact_shell_width < 0.0 ||
        !std::isfinite(aa_myosin_contact_shell_width) ||
        aa_myosin_contact_shell_energy < 0.0 ||
        !std::isfinite(aa_myosin_contact_shell_energy) ||
        ((aa_myosin_contact_shell_width > 0.0) !=
         (aa_myosin_contact_shell_energy > 0.0)) ||
        aa_zdisc_tip_onset <= 0.0 || !std::isfinite(aa_zdisc_tip_onset) ||
        aa_zdisc_tip_stiffness <= 0.0 ||
        !std::isfinite(aa_zdisc_tip_stiffness) ||
        aa_zdisc_tip_force_max <= 0.0 ||
        !std::isfinite(aa_zdisc_tip_force_max) ||
        aa_zdisc_tip_transverse_flat < 0.0 ||
        !std::isfinite(aa_zdisc_tip_transverse_flat) ||
        aa_zdisc_tip_transverse_cutoff <= aa_zdisc_tip_transverse_flat ||
        !std::isfinite(aa_zdisc_tip_transverse_cutoff) ||
        !(aa_zdisc_tip_stop_stiffness > 0.0) ||
        !std::isfinite(aa_zdisc_tip_stop_stiffness) ||
        !(aa_zdisc_tip_stop_force_max > 0.0) ||
        !std::isfinite(aa_zdisc_tip_stop_force_max) ||
        !(aa_zdisc_tip_stop_gate_radius >= 0.0) ||
        !(aa_zdisc_tip_stop_gate_cutoff > aa_zdisc_tip_stop_gate_radius) ||
        !std::isfinite(aa_zdisc_tip_stop_gate_cutoff) ||
        !(aa_zdisc_tip_stop_epsilon > 0.0) ||
        !(aa_zdisc_tip_stop_epsilon < aa_zdisc_tip_onset + 0.5 * myosin_length) ||
        aa_barbed_lateral_range < 0.0 ||
        !std::isfinite(aa_barbed_lateral_range) ||
        aa_barbed_lateral_k < 0.0 || !std::isfinite(aa_barbed_lateral_k) ||
        !std::isfinite(aa_barbed_lateral_force_max) ||
        aa_barbed_lateral_axial_window <= 0.0 ||
        !std::isfinite(aa_barbed_lateral_axial_window) ||
        aa_occlusion_rupture_penetration < 0.0 ||
        !std::isfinite(aa_occlusion_rupture_penetration) ||
        aa_occlusion_rupture_dwell < 0.0 ||
        !std::isfinite(aa_occlusion_rupture_dwell) ||
        aa_occlusion_rupture_rate < 0.0 ||
        !std::isfinite(aa_occlusion_rupture_rate) ||
        aa_barbed_tip_endstop_transverse_cutoff < 0.0 ||
        !std::isfinite(aa_barbed_tip_endstop_transverse_cutoff) ||
        !std::isfinite(aa_occlusion_capture_radius) ||
        aa_motor_pair_support_window < 0.0 ||
        !std::isfinite(aa_motor_pair_support_window) ||
        aa_motor_pair_commit_time < 0.0 ||
        !std::isfinite(aa_motor_pair_commit_time) ||
        aa_motor_pair_switch_ratio < 1.0 ||
        !std::isfinite(aa_motor_pair_switch_ratio) ||
        aa_motor_pair_exchange_time < 0.0 ||
        !std::isfinite(aa_motor_pair_exchange_time) ||
        aa_motor_pair_dropout_time < 0.0 ||
        !std::isfinite(aa_motor_pair_dropout_time) ||
        aa_association_max_direction_dot < -1.0 ||
        !std::isfinite(aa_association_max_direction_dot) ||
        aa_link_rupture_length < 0.0 ||
        !std::isfinite(aa_link_rupture_length) ||
        aa_total_link_pool < 0 || aa_stiffness_reference_links < 1 ||
        aa_stiffness_reference_links > 16) {
        std::cerr << "Invalid local AA multi-link parameters" << std::endl;
        return 1;
    }
    if (aa_arm_state_model && !aa_local_multilink) {
        std::cerr << "aa_arm_state_model requires aa_local_multilink=true"
                  << std::endl;
        return 1;
    }
    if (aa_status2_locked_until_detach && aa_status2_unlock_rate > 0.0) {
        std::cerr << "aa_status2_locked_until_detach=true conflicts with a "
                     "positive aa_status2_unlock_rate" << std::endl;
        return 1;
    }
    if (parsed_aa_status_mode ==
            actn2_load_kinetics::StatusMode::motor_transmitted_register &&
        (!dynamic_am || !aa_local_multilink)) {
        std::cerr << "motor_transmitted_register requires dynamic_am=true "
                     "and aa_local_multilink=true" << std::endl;
        return 1;
    }
    if (!std::isfinite(aa_preferred_centerline)) {
        std::cerr << "aa_preferred_centerline must be finite" << std::endl;
        return 1;
    }
    if (three_membered_transverse_jitter < 0.0 ||
        three_membered_axial_jitter < 0.0 ||
        !std::isfinite(minimal_catch_tilt_deg) ||
        std::abs(minimal_catch_tilt_deg) > 90.0 ||
        !std::isfinite(minimal_catch_myosin_tilt_deg) ||
        std::abs(minimal_catch_myosin_tilt_deg) > 90.0 ||
        !std::isfinite(minimal_catch_myosin_inward_shift) ||
        (minimal_catch_myosin_inward_shift > -1e5 &&
         std::abs(minimal_catch_myosin_inward_shift) > 5.0) ||
        !std::isfinite(am_head_radial_force_max) ||
        three_membered_tilt_jitter < 0.0 ||
        !std::isfinite(three_membered_transverse_jitter) ||
        !std::isfinite(three_membered_axial_jitter) ||
        !std::isfinite(three_membered_tilt_jitter)) {
        std::cerr << "Three-membered lattice jitters must be finite and non-negative"
                  << std::endl;
        return 1;
    }
    if (!std::isfinite(aa_transition_state_theta) ||
        aa_transition_state_theta < 0.0 || aa_transition_state_theta > 1.0) {
        std::cerr << "aa_transition_state_theta must be between 0 and 1" << std::endl;
        return 1;
    }
    if (!std::isfinite(am_activation_active_fraction) ||
        am_activation_active_fraction <= 0.0 ||
        am_activation_active_fraction > 1.0) {
        std::cerr << "am_activation_active_fraction must be in (0, 1]"
                  << std::endl;
        return 1;
    }
    if (!std::isfinite(actin_turnover_rate) || actin_turnover_rate < 0.0 ||
        !std::isfinite(actin_turnover_start_time) ||
        actin_turnover_start_time < 0.0 ||
        actin_turnover_max_placement_attempts <= 0 ||
        !std::isfinite(actin_turnover_clearance_ratio) ||
        actin_turnover_clearance_ratio <= 0.0 ||
        actin_turnover_clearance_ratio > 1.0) {
        std::cerr << "Actin turnover parameters must be finite and non-negative, "
                     "with positive placement attempts" << std::endl;
        return 1;
    }
    if (aa_passive_detailed_balance &&
        (lifetime_coeff != 0.0 || directional_base_lifetime != 0.0 ||
         tau_rec != 0.0)) {
        std::cerr << "aa_passive_detailed_balance requires lifetime_coeff=0, "
                     "directional_base_lifetime=0, and bond_recovery_time=0"
                  << std::endl;
        return 1;
    }

    if (!resume && resume_frame >= 0) {
        std::cout << "Warning: --resume_frame ignored because --resume was not set.\n";
        resume_frame = -1;
    }
    if (resume_frame < -1) {
        std::cout << "Warning: --resume_frame must be -1 or non-negative. Using latest frame instead.\n";
        resume_frame = -1;
    }

    gsl_rng * rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng,seed);
    std::vector <double> box(3);
    box[0] = Lx;
    box[1] = Ly;
    box[2] = Lz;
    rod_mobility::Diffusion actin_mobility;
    rod_mobility::Diffusion myosin_mobility;
    if (physical_mobility) {
        try {
            const double temperature_beta = rod_mobility::beta(temperature_kelvin);
            if (!beta_was_explicit) {
                beta = temperature_beta;
            } else {
                std::cerr << "Warning: explicit --beta=" << beta
                          << " overrides temperature-derived beta="
                          << temperature_beta << ". This run is a diagnostic "
                             "thermal override.\n";
            }
            actin_mobility = rod_mobility::finite_cylinder(
                actin_length, actin_hydrodynamic_diameter,
                temperature_kelvin, solvent_viscosity,
                actin_crowding_parallel, actin_crowding_perpendicular,
                actin_crowding_tumble, actin_crowding_roll);
            myosin_mobility = rod_mobility::finite_cylinder(
                myosin_length, myosin_hydrodynamic_diameter,
                temperature_kelvin, solvent_viscosity,
                myosin_crowding_parallel, myosin_crowding_perpendicular,
                myosin_crowding_tumble, myosin_crowding_roll);
        } catch (const std::exception& error) {
            std::cerr << "Invalid physical mobility configuration: "
                      << error.what() << std::endl;
            gsl_rng_free(rng);
            return 1;
        }
        std::cout << "Physical rod mobility at T=" << temperature_kelvin
                  << " K, beta=" << beta << "\n"
                  << "  actin D_parallel=" << actin_mobility.parallel
                  << ", D_perp=" << actin_mobility.perpendicular
                  << ", D_tumble=" << actin_mobility.tumble
                  << ", D_roll=" << actin_mobility.roll << "\n"
                  << "  myosin D_parallel=" << myosin_mobility.parallel
                  << ", D_perp=" << myosin_mobility.perpendicular
                  << ", D_tumble=" << myosin_mobility.tumble
                  << ", D_roll=" << myosin_mobility.roll << '\n';
    }
    double diff_coeff_ratio = myosin_diff_coeff_trans > 0.0
        ? actin_diff_coeff_trans / myosin_diff_coeff_trans
        : std::numeric_limits<double>::infinity();
    double titin_rest_length = titin_rest_length_option > 0.0
        ? titin_rest_length_option
        : 0.5 * myosin_length + (2.0 / 3.0) * actin_length;
    Sarcomere::TitinMode titin_mode;
    if (titin_mode_name == "legacy_crosslinked" ||
        titin_mode_name == "legacy") {
        titin_mode = Sarcomere::TitinMode::legacy_crosslinked;
    } else if (titin_mode_name == "balanced_owning" ||
               titin_mode_name == "balanced") {
        titin_mode = Sarcomere::TitinMode::balanced_owning;
    } else if (titin_mode_name == "explicit_dynamic" ||
               titin_mode_name == "explicit") {
        titin_mode = Sarcomere::TitinMode::explicit_dynamic;
    } else {
        std::cerr << "Unknown --titin_mode=" << titin_mode_name
                  << "; expected legacy_crosslinked, balanced_owning, or "
                     "explicit_dynamic"
                  << std::endl;
        gsl_rng_free(rng);
        return 1;
    }
    if (fixed_lattice_preequil_steps < 0) {
        std::cerr << "fixed_lattice_preequil_steps must be non-negative"
                  << std::endl;
        gsl_rng_free(rng);
        return 1;
    }
    if (titin_mode == Sarcomere::TitinMode::balanced_owning &&
        (n_actins != 4 || n_myosins != 2)) {
        std::cerr << "balanced_owning titin requires --n_actins=4 and "
                     "--n_myosins=2" << std::endl;
        gsl_rng_free(rng);
        return 1;
    }
    if (titin_mode == Sarcomere::TitinMode::explicit_dynamic &&
        (!(titin_k > 0.0) || !(titin_rest_length_option > 0.0) ||
         titin_links_per_half <= 0 || !(titin_capture_distance > 0.0) ||
         titin_k_on < 0.0 || titin_k_off < 0.0 ||
         titin_transition_theta < 0.0 || titin_transition_theta > 1.0 ||
         titin_rebind_memory < 0.0 || !std::isfinite(titin_rebind_memory) ||
         titin_rebind_rate_multiplier < 0.0 ||
         !std::isfinite(titin_rebind_rate_multiplier) ||
         !std::isfinite(titin_k_off_unsupported))) {
        std::cerr << "explicit_dynamic titin requires positive --titin_k, "
                     "--titin_rest_length, --titin_links_per_half, and "
                     "--titin_capture_distance, non-negative rates, and a "
                     "transition theta in [0,1]" << std::endl;
        gsl_rng_free(rng);
        return 1;
    }
    // The gate reads the sparse local ACTN2 link index. Without that model
    // there are no A-A links at all, so every ordinary association would be
    // rejected and the run would silently produce zero bound titins.
    if (titin_mode == Sarcomere::TitinMode::explicit_dynamic &&
        (titin_require_live_actn2 || titin_k_off_unsupported > 0.0) &&
        !aa_local_multilink) {
        std::cerr << "--titin_require_live_actn2 and "
                     "--titin_k_off_unsupported require "
                     "--aa_local_multilink=true" << std::endl;
        gsl_rng_free(rng);
        return 1;
    }
    auto compute_max_force = [&](double max_disp, double diff_coeff) {
        if (max_disp <= 0.0 || diff_coeff <= 0.0 || !std::isfinite(beta) || !std::isfinite(dt)) {
            return std::numeric_limits<double>::infinity();
        }
        double val = max_disp / (beta * diff_coeff * dt);
        if (!(val > 0.0) || !std::isfinite(val)) {
            return std::numeric_limits<double>::infinity();
        }
        return val;
    };
    auto compute_max_torque = [&](double max_rotation, double diff_coeff_rot) {
        if (max_rotation <= 0.0 || diff_coeff_rot <= 0.0 || !std::isfinite(beta) || !std::isfinite(dt)) {
            return std::numeric_limits<double>::infinity();
        }
        double val = max_rotation / (beta * diff_coeff_rot * dt);
        if (!(val > 0.0) || !std::isfinite(val)) {
            return std::numeric_limits<double>::infinity();
        }
        return val;
    };
    const double actin_force_diffusion = physical_mobility
        ? std::max(actin_mobility.parallel, actin_mobility.perpendicular)
        : actin_diff_coeff_trans;
    const double myosin_force_diffusion = physical_mobility
        ? std::max(myosin_mobility.parallel, myosin_mobility.perpendicular)
        : myosin_diff_coeff_trans;
    const double actin_rotation_diffusion = physical_mobility
        ? actin_mobility.tumble
        : actin_diff_coeff_rot;
    const double myosin_rotation_diffusion = physical_mobility
        ? myosin_mobility.tumble
        : myosin_diff_coeff_rot;
    double max_actin_force = compute_max_force(max_actin_displacement, actin_force_diffusion);
    double max_myosin_force = compute_max_force(max_myosin_displacement, myosin_force_diffusion);
    double max_actin_torque = compute_max_torque(max_actin_rotation, actin_rotation_diffusion);
    double max_myosin_torque = compute_max_torque(max_myosin_rotation, myosin_rotation_diffusion);
    if (std::isfinite(max_actin_torque)) {
    } else {
    }
    if (std::isfinite(max_myosin_torque)) {
    } else {
    }
    const double am_preferred_centerline =
        lattice_spacing::am_preferred_centerline(
            myosin_radius, actin_radius, am_surface_gap);
    aa_optimal = aa_preferred_centerline > 0.0
        ? aa_preferred_centerline
        : lattice_spacing::aa_preferred_centerline(
              am_preferred_centerline, aa_spacing_ratio);
    aa_spacing_ratio = lattice_spacing::effective_ratio(
        aa_optimal, am_preferred_centerline);
    aa_cutoff = aa_optimal + aa_capture_margin;
    std::cout << "Lattice spacing: d_AM*=" << am_preferred_centerline
              << ", d_AA*=" << aa_optimal
              << ", rho=" << aa_spacing_ratio
              << ", AA capture cutoff=" << aa_cutoff << '\n';
    Sarcomere model(n_actins, n_myosins, box, actin_length, myosin_length,
                        myosin_radius, am_cutoff, am_optimal, aa_cutoff, aa_optimal,
                        k_on, k_off,
                        base_lifetime, directional_base_lifetime, lifetime_coeff, diff_coeff_ratio,
                          k_aa, kappa_aa, k_am, kappa_am, k_mm, v_am,
                        filename,rng, seed, n_fixed_myosins, dt, tau_rec,
                        titin_k, titin_rest_length,
                        directional, max_myosin_bonds, max_actin_force, max_myosin_force,
                        max_actin_torque, max_myosin_torque);
    model.stabilize_status1 = stabilize_status1;
    model.set_aa_status_kinetics(
        parsed_aa_status_mode,
        aa_status2_force_on, aa_status2_force_off);
    model.set_am_kinetics(
        dynamic_am, am_k_on, am_k_weak_off, am_k_powerstroke,
        am_k_strong_off, am_force_scale, am_motor_k, am_power_stroke,
        am_stall_force,
        am_unloaded_velocity,
        am_weak_stiffness_ratio, am_rupture_distance, am_cycle_period,
        am_activation_rise, am_activation_decay, am_activation_floor,
        am_status2_memory, am_activation_active_fraction);
    model.set_am_proximity_tether(
        am_tether_k, am_tether_slack, am_binding_azimuth_selection,
        am_binding_azimuth_half_width_deg, am_alignment_k,
        am_alignment_slack_deg);
    model.set_am_head_transverse_tether(am_head_transverse_k);
    model.set_am_retention_kinetics(
        am_strain_dependent_detachment, am_strain_detachment_theta,
        am_transverse_strain_detachment_theta,
        am_strain_detachment_max_exponent, am_tether_saturation_heads,
        am_strong_off_floor_fraction);
    model.am_offcenter_torque = am_offcenter_torque;
    model.am_spatial_heads = am_spatial_heads;
    model.sterics_enabled = sterics_enabled;
    model.set_am_radial_geometry(actin_radius, am_surface_gap, am_steric_k);
    model.set_steric_stiffness(actin_steric_k, myosin_steric_k, myosin_exclusion_distance);
    model.configure_myomesin_bridges(
        myomesin_sites_per_filament, myomesin_k_on, myomesin_k_off,
        myomesin_axial_capture, myomesin_radial_rest,
        myomesin_radial_capture, myomesin_radial_k,
        myomesin_alignment_cos, myomesin_bell_distance,
        myomesin_force_max);
    model.am_stroke_detailed_balance = am_stroke_detailed_balance;
    model.am_stroke_base_rate = am_stroke_base_rate;
    model.am_stroke_theta = am_stroke_theta;
    model.am_stroke_free_energy = am_stroke_free_energy;
    model.configure_double_overlap_contact(
        aa_double_overlap_contact,
        aa_double_overlap_contact_diameter,
        aa_double_overlap_contact_epsilon,
        aa_double_overlap_contact_force_cap);
    model.set_aa_spacing_geometry(
        aa_spacing_ratio, am_preferred_centerline, aa_capture_margin);
    model.set_aa_axial_overlap(
        aa_antiparallel_overlap, aa_axial_register_ratio,
        aa_status2_registry_on, aa_status2_registry_off,
        aa_register_saturation_error);
    model.am_bidirectional_crossbridge = am_bidirectional_crossbridge;
    model.am_crossover_pull_through = am_crossover_pull_through;
    model.am_overlap_target_ratio = am_overlap_target_ratio;
    model.aa_axial_register_force_max = aa_axial_register_force_max;
    model.aa_radial_force_max = aa_radial_force_max;
    model.am_tether_force_max = am_tether_force_max;
    model.am_tether_occupancy_floor = am_tether_occupancy_floor;
    model.am_tether_occupancy_floor_reach = am_tether_occupancy_floor_reach;
    model.probe_tilt_deg = probe_tilt_deg;
    model.probe_confine_k = probe_confine_k;
    model.aa_register_all_statuses = aa_register_all_statuses;
    model.aa_perpendicular_all_statuses = aa_perpendicular_all_statuses;
    model.aa_barbed_zone = aa_barbed_zone;
    model.aa_status1_radial_k = aa_status1_radial_k;
    model.aa_status1_site_diffusion = aa_status1_site_diffusion;
    model.aa_status2_locked_until_detach = aa_status2_locked_until_detach;
    model.aa_status2_unlock_rate = aa_status2_unlock_rate;
    model.aa_status2_unlock_force_scale = aa_status2_unlock_force_scale;
    model.aa_status2_unlock_registry_scale =
        aa_status2_unlock_registry_scale;
    model.aa_arm_state_model = aa_arm_state_model;
    model.aa_arm_open_rate_weak = aa_arm_open_rate_weak;
    model.aa_arm_open_rate_mature = aa_arm_open_rate_mature;
    model.aa_arm_bell_distance = aa_arm_bell_distance;
    model.aa_tether_off_rate = aa_tether_off_rate;
    model.aa_tether_rebind_rate = aa_tether_rebind_rate;
    model.aa_tether_capture_distance = aa_tether_capture_distance;
    model.aa_arm_rebind_theta = aa_arm_rebind_theta;
    model.aa_actn2_contact_radius = aa_actn2_contact_radius;
    model.aa_myosin_contact_epsilon = aa_myosin_contact_epsilon;
    model.aa_myosin_contact_force_cap = aa_myosin_contact_force_cap;
    model.aa_myosin_contact_shell_width = aa_myosin_contact_shell_width;
    model.aa_myosin_contact_shell_energy = aa_myosin_contact_shell_energy;
    model.aa_zdisc_tip_contact = aa_zdisc_tip_contact;
    model.aa_zdisc_tip_onset = aa_zdisc_tip_onset;
    model.aa_zdisc_tip_stiffness = aa_zdisc_tip_stiffness;
    model.aa_zdisc_tip_force_max = aa_zdisc_tip_force_max;
    model.aa_zdisc_tip_contact_requires_link = aa_zdisc_tip_contact_requires_link;
    model.aa_zdisc_tip_hard_projection = aa_zdisc_tip_hard_projection;
    model.aa_zdisc_tip_transverse_flat = aa_zdisc_tip_transverse_flat;
    model.aa_zdisc_tip_transverse_cutoff = aa_zdisc_tip_transverse_cutoff;
    model.aa_zdisc_tip_stop_stiffness = aa_zdisc_tip_stop_stiffness;
    model.aa_zdisc_tip_stop_force_max = aa_zdisc_tip_stop_force_max;
    model.aa_zdisc_tip_stop_gate_radius = aa_zdisc_tip_stop_gate_radius;
    model.aa_zdisc_tip_stop_gate_cutoff = aa_zdisc_tip_stop_gate_cutoff;
    model.aa_zdisc_tip_stop_epsilon = aa_zdisc_tip_stop_epsilon;
    model.aa_barbed_lateral_range = aa_barbed_lateral_range;
    model.aa_barbed_lateral_k = aa_barbed_lateral_k;
    model.aa_barbed_lateral_force_max = aa_barbed_lateral_force_max;
    model.aa_barbed_lateral_axial_window = aa_barbed_lateral_axial_window;
    model.aa_occlusion_rupture_penetration = aa_occlusion_rupture_penetration;
    model.aa_occlusion_capture_radius = aa_occlusion_capture_radius;
    model.aa_overlap_occlusion_veto = aa_overlap_occlusion_veto;
    model.aa_occlusion_rupture_dwell_steps =
        (aa_occlusion_rupture_dwell > 0.0 && dt > 0.0)
            ? static_cast<size_t>(std::llround(aa_occlusion_rupture_dwell / dt))
            : static_cast<size_t>(0);
    model.aa_occlusion_rupture_rate = aa_occlusion_rupture_rate;
    model.aa_barbed_tip_endstop_slab = aa_barbed_tip_endstop_slab;
    model.aa_barbed_tip_endstop_transverse_cutoff =
        aa_barbed_tip_endstop_transverse_cutoff;
    if (model.aa_zdisc_tip_contact) {
        std::cout
            << "ACTN2 tip projection traversal: spatial_incidence "
               "(default; exact conservative broad phase)"
            << std::endl;
    }
    model.titin_mode = titin_mode;
    if (titin_mode == Sarcomere::TitinMode::explicit_dynamic) {
        model.configure_explicit_titin(
            titin_links_per_half, titin_capture_distance,
            titin_k_on, titin_k_off, titin_strain_gated_kinetics,
            titin_transition_theta, titin_rebind_memory,
            titin_rebind_rate_multiplier, titin_require_live_actn2,
            titin_k_off_unsupported, titin_rebind_requires_actn2);
    }
    model.release_step_origin = static_cast<size_t>(
        fixed_lattice_preequil_steps);
    model.set_local_aa_model(
        aa_local_multilink, aa_prebundle_links_per_pair,
        aa_site_exclusion, aa_overlap_capture, aa_perpendicular_k,
        aa_total_link_pool, aa_stiffness_reference_links);
    model.set_aa_neighbor_identity(
        aa_neighbor_identity, aa_neighbor_slot_memory,
        aa_status2_slot_memory);
    model.set_physical_mobility(
        temperature_kelvin, beta, actin_mobility, myosin_mobility);
    model.legacy_actin_D_trans = actin_diff_coeff_trans;
    model.legacy_actin_D_rot = actin_diff_coeff_rot;
    model.legacy_myosin_D_trans = myosin_diff_coeff_trans;
    model.legacy_myosin_D_rot = myosin_diff_coeff_rot;
    model.thermal_noise_scale_metadata = thermal_noise_scale;
    model.six_ring_angle_noise_deg = std::max(0.0, six_ring_angle_noise_deg);
    model.six_ring_radial_noise = std::max(0.0, six_ring_radial_noise);
    model.six_ring_axis_noise_deg = std::max(0.0, six_ring_axis_noise_deg);
    model.max_aa_bonds_per_actin = max_aa_bonds_per_actin;
    model.aa_status2_support_window = aa_status2_support_window;
    model.aa_motor_pair_support_window = aa_motor_pair_support_window;
    model.aa_motor_pair_commit_time = aa_motor_pair_commit_time;
    model.aa_motor_pair_switch_ratio = aa_motor_pair_switch_ratio;
    model.aa_motor_pair_exchange_time = aa_motor_pair_exchange_time;
    model.aa_motor_pair_dropout_time = aa_motor_pair_dropout_time;
    model.aa_energy_gated_association = aa_energy_gated_association;
    model.aa_transition_state_theta = aa_transition_state_theta;
    model.aa_passive_detailed_balance = aa_passive_detailed_balance;
    model.aa_encounter_k = aa_encounter_k;
    model.aa_encounter_kappa = aa_encounter_kappa;
    model.aa_association_max_direction_dot =
        aa_association_max_direction_dot;
    model.aa_link_rupture_length = aa_link_rupture_length;
    model.set_interaction_ramp(0, static_cast<size_t>(interaction_ramp_steps));
    model.isometric_actin_axial = isometric_actin_axial;
    if (am_spokes_per_crown < 1 || am_spokes_per_crown > 16) {
        std::cerr << "am_spokes_per_crown must be between 1 and 16" << std::endl;
        gsl_rng_free(rng);
        return 1;
    }
    // Must precede set_am_head_count: that call builds the head pools.
    model.am_spokes_per_crown = am_spokes_per_crown;
    if (am_head_radial_mode < 0 || am_head_radial_mode > 2) {
        std::cerr << "ERROR: --am_head_radial_mode must be 0, 1 or 2\n";
        return 1;
    }
    if (am_head_radial_k_override < 0.0 ||
        !std::isfinite(am_head_radial_k_override) ||
        am_head_radial_rest_override < 0.0 ||
        !std::isfinite(am_head_radial_rest_override)) {
        std::cerr << "ERROR: --am_head_radial_k / --am_head_radial_rest must be finite and >= 0\n";
        return 1;
    }
    if (am_alignment_saturation_heads < 0 ||
        am_alignment_saturation_heads > 64) {
        std::cerr << "ERROR: --am_alignment_saturation_heads must be in [0, 64]\n";
        return 1;
    }
    if (am_motor_saturation_heads < 0 || am_motor_saturation_heads > 64) {
        std::cerr << "ERROR: --am_motor_saturation_heads must be in [0, 64]\n";
        return 1;
    }
    model.am_motor_saturation_heads = am_motor_saturation_heads;
    model.am_alignment_saturation_heads = am_alignment_saturation_heads;
    model.am_head_radial_mode = am_head_radial_mode;
    model.am_head_radial_k_override = am_head_radial_k_override;
    model.am_head_radial_rest_override = am_head_radial_rest_override;
    model.am_head_azimuth = am_head_azimuth;
    model.am_azimuth_gate = am_azimuth_gate;
    {
        // <=0 means "partition the circle": with n spokes 2*pi/n apart, a
        // half-width of pi/n makes exactly one spoke admit any given actin.
        const double deg = am_azimuth_gate_half_width_deg > 0.0
            ? am_azimuth_gate_half_width_deg
            : 180.0 / static_cast<double>(std::max(1, am_spokes_per_crown));
        model.am_azimuth_gate_half_width = deg * M_PI / 180.0;
    }
    if (am_azimuth_gate && !am_spatial_heads) {
        std::cerr << "ERROR: --am_azimuth_gate requires --am_spatial_heads, "
                     "because without spatial heads there is no crown azimuth "
                     "or myosin material frame to gate on.\n";
        return 1;
    }
    // Derive the radial cap from the displacement clamp unless given explicitly.
    // mu_perp = beta * D_perp is the mobility the clamp actually acts through,
    // so F_clamp = max_actin_displacement/(mu_perp*dt) is the largest force the
    // integrator can express on one actin.  Divide by the occupancy scale so a
    // full complement of bound heads reaches that ceiling and no more.
    if (am_head_radial_force_max < 0.0) {
        const double d_perp = physical_mobility
            ? actin_mobility.perpendicular : actin_diff_coeff_trans;
        const double mu_perp = beta * d_perp;
        const int n_sat = std::max(1, am_tether_saturation_heads);
        if (mu_perp > 0.0 && dt > 0.0 && max_actin_displacement > 0.0) {
            am_head_radial_force_max = max_actin_displacement /
                (mu_perp * dt * static_cast<double>(n_sat));
            std::cout << "am_head_radial_force_max derived from the "
                         "displacement clamp: " << am_head_radial_force_max
                      << " pN per head (clamp "
                      << max_actin_displacement / (mu_perp * dt)
                      << " pN total, / " << n_sat << " heads)" << std::endl;
        } else {
            am_head_radial_force_max = 0.0;
            std::cout << "am_head_radial_force_max: cannot derive a clamp-based "
                         "cap (mobility, dt or max_actin_displacement is zero); "
                         "falling back to the unbounded linear spring"
                      << std::endl;
        }
    }
    model.am_head_radial_force_max = am_head_radial_force_max;
    model.set_am_head_count(am_heads_per_half);
    model.set_am_neighbor_head_cap(am_max_heads_per_actin_half);
    model.set_am_site_count(am_sites_per_actin);
    model.set_am_global_head_cap(am_max_heads_per_actin);
    model.set_am_weak_binding_floor(am_weak_binding_floor);
    model.set_am_thin_activation(
        am_thin_activation, am_thin_feedback,
        am_q_calcium_open_rate, am_q_calcium_half, am_q_calcium_hill,
        am_q_strong_feedback_rate, am_q_strong_half, am_q_strong_hill,
        am_q_close_rate);
    if (am_initialize_q_calcium_equilibrium) {
        model.initialize_am_thin_activation_calcium_equilibrium();
    }
    model.set_am_thick_regulation(
        am_thick_feedback, am_thick_recruit_rate, am_thick_fold_rate,
        am_thick_tension_scale, am_thick_initial_on_fraction,
        am_thick_force_clamp);
    model.configure_actin_turnover(
        actin_turnover, actin_turnover_rate, actin_turnover_start_time,
        actin_turnover_max_placement_attempts,
        actin_turnover_clearance_ratio);
    model.set_bundling_parameters(k_bundle_max, bundle_ramp_steps);
    model.set_periodicity(periodic_axes);
    if (!resume && init_struc == "geometry_h5") {
        try {
            const int loaded_frame = model.load_initial_geometry(
                initial_geometry_file, initial_geometry_frame);
            std::cout << "Imported filament geometry from "
                      << initial_geometry_file << " frame " << loaded_frame
                      << "; molecular states initialized fresh" << std::endl;
        } catch (const std::exception& error) {
            std::cerr << "Initial geometry import failed: " << error.what()
                      << std::endl;
            return 1;
        }
    }
    if (!is3D) {
        for (int i = 0; i < n_actins; ++i) {
            model.actin.center[i].z = 0;
            model.actin.direction[i].z = 0;
            model.actin.force[i].z = 0;
            model.actin.velocity[i].z = 0;
            model.actin.torque[i].z = 0;
        }
        for (int i = 0; i < n_myosins; ++i) {
            model.myosin.center[i].z = 0;
            model.myosin.direction[i].z = 0;
            model.myosin.force[i].z = 0;
            model.myosin.velocity[i].z = 0;
            model.myosin.torque[i].z = 0;
        }
        model.actin.update_endpoints();
        model.myosin.update_endpoints();
    }

    // Configuration audit.  Each of these was a silent trap found on
    // 2026-09-17; they are reported here so a run cannot quietly inherit one.
    if (titin_mode == Sarcomere::TitinMode::explicit_dynamic &&
        titin_require_live_actn2 && !(titin_k_off_unsupported > 0.0)) {
        std::cerr << "WARNING: --titin_require_live_actn2 is set but "
                     "--titin_k_off_unsupported <= 0, which per its own help "
                     "text keeps the single global --titin_k_off="
                  << titin_k_off << "/s.  An unsupported titin therefore has a "
                     "mean lifetime of "
                  << (titin_k_off > 0.0 ? 1.0 / titin_k_off : 0.0)
                  << " s and effectively never releases, so it will pin the "
                     "thin filament after the Z-disc dies." << std::endl;
    }
    if (max_actin_rotation > 0.005 && !legacy_actin_torque_cap) {
        std::cerr << "WARNING: --max_actin_rotation = " << max_actin_rotation
                  << " exceeds the production 0.005.  Values near 0.020 were "
                     "adopted to dodge actin torque-clamp saturation, which "
                     "was a BUG: the clamp budgeted the full torque gradient "
                     "while the integrator uses only the component "
                     "perpendicular to the filament axis, and actin has no "
                     "roll channel.  That is fixed (torque cap events went "
                     "3,996,032 -> 0 per seed; measured torq_sat is 0.0000 at "
                     "0.005), so the inflated cap buys nothing and only "
                     "loosens the integrator.  Use 0.005 unless you have a "
                     "reason unrelated to clamp saturation."
                  << std::endl;
    }
    // Only meaningful in the legacy arm: azimuth-OFF folds the pair tether
    // into the per-head radial spring and skips the aggregate block entirely
    // (sarcomere.cpp, guarded by `if (am_head_azimuth)`), so the recorded
    // am_tether_force is identically 0 there.  Warning unconditionally pointed
    // at a dead code path and implied production ran an uncapped AM force when
    // it is in fact capped by am_head_radial_force_max.
    if (am_head_azimuth && am_tether_k > 0.0 && !(am_tether_force_max > 0.0)) {
        std::cerr << "WARNING: --am_tether_force_max <= 0 selects the unbounded "
                     "linear branch of the AM tether; its force grows without "
                     "limit and is bounded only by the integrator displacement "
                     "clamp.  (Legacy arm only: with --am_head_azimuth=false "
                     "the aggregate tether is not evaluated and the per-head "
                     "radial force saturates at am_head_radial_force_max.)"
                  << std::endl;
    }
    if (am_head_transverse_k > 0.0 &&
        !(am_transverse_strain_detachment_theta > 0.0)) {
        std::cerr << "WARNING: --am_transverse_strain_detachment_theta = 0, so "
                     "transverse over-stretch of a bound head does not "
                     "accelerate its detachment.  A radially stretched head "
                     "pulls indefinitely instead of letting go." << std::endl;
    }
    if (!am_spatial_heads) {
        std::cerr << "ERROR: --am_spatial_heads=false is rejected.  That path "
                     "places the head anchor on the myosin AXIS (radial offset "
                     "0) while still passing a rest length of actin_radius + "
                     "am_surface_gap = 30 nm, so every bound head starts "
                     "stretched by about 50 nm and pulls inward at roughly 5 "
                     "pN.  Supply the matching rest length before enabling it."
                  << std::endl;
        gsl_rng_free(rng);
        return 1;
    }
    if (am_tether_occupancy_floor > 0.0) {
        std::cerr << "ERROR: --am_tether_occupancy_floor has no effect and is "
                     "rejected.  The block applying it (sarcomere.cpp ~6775) "
                     "sits inside a loop that skips any pair with "
                     "am_bonds[i][j] != 1, i.e. any pair with no bound head, so "
                     "the floor can only raise occupancy for pairs that already "
                     "have one.  Its stated purpose -- tethering a pair that has "
                     "lost its heads -- is unreachable.  Wire it before "
                     "sweeping it." << std::endl;
        gsl_rng_free(rng);
        return 1;
    }

    // Captured above, while the cxxopts Options object was still in scope.
    model.cli_invocation = cli_invocation_capture;
    model.cli_defaults = cli_defaults_capture;

    // Must be set BEFORE the Langevin constructor: it calls model.new_file(),
    // which writes every meta attribute.  Setting these only in the
    // initial_structure dispatch below persisted 0.0 for every tilt, leaving the
    // runs correct but unreproducible from their own metadata.  Langevin holds
    // Sarcomere by reference, so this is the same object.
    model.minimal_catch_tilt_deg = minimal_catch_tilt_deg;
    model.minimal_catch_myosin_tilt_deg = minimal_catch_myosin_tilt_deg;
    model.minimal_catch_myosin_inward_shift = minimal_catch_myosin_inward_shift;
    // Features must be registered BEFORE the Langevin constructor creates the
    // HDF5 file, or the run aborts on a missing dataset.
    model.torque_breakdown = torque_breakdown;
    model.legacy_actin_torque_cap = legacy_actin_torque_cap;
    model.actin_direction_polar_fraction = actin_direction_polar_fraction;
    model.myosin_direction_polar_fraction = myosin_direction_polar_fraction;
    model._tb_setup();

    Langevin sim(model, beta, dt, actin_diff_coeff_trans,actin_diff_coeff_rot, myosin_diff_coeff_trans, myosin_diff_coeff_rot, save_every, resume, is3D,
                    max_actin_displacement, max_myosin_displacement,
                    max_actin_rotation, max_myosin_rotation, resume_frame);
    sim.set_rod_mobility(actin_mobility, myosin_mobility);
    sim.set_isometric_actin_axial(isometric_actin_axial);
    sim.set_thermal_noise_scale(thermal_noise_scale);
    if (!resume){
        if (init_struc == "sarcomere") {
        sim.model.sarcomeric_structure_tight();}

        else if (init_struc == "partial"){
            sim.model.partial_fix(n_fixed_myosins);
        }
        else if (init_struc == "cb"){
            sim.model.cb();
        }
        else if (init_struc == "minimal_catch") {
            sim.model.minimal_catch();
        }
        else if (init_struc.rfind("minimal_catch_", 0) == 0) {
            sim.model.init_aa_overlap = init_aa_overlap;
            sim.model.minimal_catch_tilt_deg = minimal_catch_tilt_deg;
            sim.model.minimal_catch_myosin_tilt_deg =
                minimal_catch_myosin_tilt_deg;
            sim.model.minimal_catch_myosin_inward_shift =
                minimal_catch_myosin_inward_shift;
            sim.model.minimal_catch_variant(
                init_struc.substr(std::string("minimal_catch_").size()));
        }
        else if (init_struc == "aa_bundle_common") {
            sim.model.aa_bundle_common();
        }
        else if (init_struc == "aa_bundle_opposed") {
            sim.model.aa_bundle_opposed();
        }
        else if (init_struc == "cb_off_angle"){
            sim.model.cb_off_angle();
        }
        else if (init_struc == "single_am"){
            sim.model.single_am();
        }
        else if (init_struc == "single_am_tilted"){
            sim.model.single_am_tilted();
        }
        else if (init_struc == "six_actin_ring"){
            sim.model.six_actin_ring();
            if (six_ring_prebind_aa) sim.model.prebind_six_actin_ring();
        }
        else if (init_struc == "bilateral_six_actin_ring"){
            sim.model.bilateral_six_actin_ring();
            if (six_ring_prebind_aa) {
                sim.model.prebind_bilateral_six_actin_ring();
            }
        }
        else if (init_struc == "force_balanced_six_actin_ring"){
            sim.model.force_balanced_six_actin_cell();
            if (six_ring_prebind_aa) {
                sim.model.prebind_force_balanced_six_actin_cell();
            }
            sim.model.seed_force_balanced_six_actin_heads();
        }
        else if (init_struc == "shared_thick_thin_lattice"){
            sim.model.shared_thick_thin_lattice();
            if (six_ring_prebind_aa) {
                sim.model.prebind_shared_thick_thin_lattice(
                    sim.model.aa_prebundle_links_per_pair);
            }
            sim.model.seed_shared_thick_thin_heads();
        }
        else if (init_struc == "two_sarcomere_unit"){
            sim.model.init_aa_overlap = init_aa_overlap;
            sim.model.two_sarcomere_unit();
        }
        else if (init_struc == "three_membered_lattice"){
            sim.model.init_aa_overlap = init_aa_overlap;
            sim.model.three_membered_lattice(
                six_ring_prebind_aa, three_membered_transverse_jitter,
                three_membered_axial_jitter, three_membered_tilt_jitter);
        }
        else if (init_struc == "minimal_sarcomere"){
            sim.model.minimal_sarcomere();
        }
        else if (init_struc == "myosin_x_noise"){
            if (actin_direction_noise >= 0.0) {
                sim.model.set_actin_direction_x_noise(
                    actin_direction_noise, actin_direction_polar_fraction);
            }
            sim.model.set_myosin_direction_x_noise(
                myosin_direction_noise, myosin_direction_polar_fraction,
                myosin_direction_tilt_deg);
            sim.volume_exclusion(pre_phase_steps, rng, n_fixed_myosins, false);
            if (prebundle_assembly) sim.model.prebind_assembly_actin_links();
            sim.model.current_step = 0;
        }
        else if (init_struc == "geometry_h5") {
            // A zero-step pre-phase preserves the imported pose exactly;
            // positive values request sterics-only relaxation before chemistry.
            if (pre_phase_steps > 0) {
                sim.volume_exclusion(
                    pre_phase_steps, rng, n_fixed_myosins, false);
            }
            sim.model.current_step = 0;
        }
        else{
            sim.volume_exclusion(pre_phase_steps, rng, n_fixed_myosins, false);
            sim.model.current_step = 0;}
        sim.model.reset_myosin_material_frames();
        if (fixed_lattice_preequil_steps > 0) {
            std::cout << "Fixed-lattice chemistry pre-equilibration: "
                      << fixed_lattice_preequil_steps << " steps ("
                      << fixed_lattice_preequil_steps * dt * 1.0e3
                      << " ms)" << std::endl;
            for (int step = 0; step < fixed_lattice_preequil_steps; ++step) {
                sim.model.update_system();
            }
            sim.model.reset_release_diagnostics();
            std::cout << "Released at KMC step "
                      << sim.model.current_step << std::endl;
        }
    }
    sim.run_langevin(nsteps, rng, n_fixed_myosins);
    const auto& tip_profile = sim.actn2_tip_projection_profile();
    if (tip_profile.calls > 0) {
        std::cout << "actn2_tip_projection profile: traversal=spatial_incidence"
                  << ", calls=" << tip_profile.calls
                  << ", active_links=" << tip_profile.active_links
                  << ", unique_actin_pairs=" << tip_profile.unique_actin_pairs
                  << ", constraints=" << tip_profile.constraints
                  << ", broadphase_tests=" << tip_profile.broadphase_tests
                  << ", broadphase_candidates="
                  << tip_profile.broadphase_candidates
                  << ", solve_evaluations=" << tip_profile.solve_evaluations
                  << ", audit_evaluations=" << tip_profile.audit_evaluations
                  << ", corrections=" << tip_profile.corrections
                  << ", sweeps=" << tip_profile.sweeps
                  << ", exhausted_calls=" << tip_profile.exhausted_calls
                  << ", audited_calls=" << tip_profile.audited_calls
                  << ", certified_without_audit_calls="
                  << tip_profile.certified_without_audit_calls
                  << ", build_wall_s=" << tip_profile.build_seconds
                  << ", solve_wall_s=" << tip_profile.solve_seconds
                  << ", audit_wall_s=" << tip_profile.audit_seconds
                  << ", total_wall_s=" << tip_profile.total_seconds
                  << std::endl;
    }
    if (parsed_aa_status_mode ==
        actn2_load_kinetics::StatusMode::motor_transmitted_register) {
        const auto& profile = sim.model;
        std::cout << "motor_transmitted_register profile: steps="
                  << profile.aa_motor_profile_steps
                  << ", local_update_wall_s="
                  << profile.aa_motor_local_update_seconds
                  << ", source_build_thread_s="
                  << profile.aa_motor_support_build_seconds
                  << ", candidate_scan_thread_s="
                  << profile.aa_motor_candidate_seconds
                  << ", candidate_pairs="
                  << profile.aa_motor_candidate_pairs_evaluated << std::endl;
    }
    gsl_rng_free(rng);
    return 0;
}
