#include "langevin.h"
#include <cstdio>
#include <cmath>
#include <vector>
#include <algorithm>
#include <chrono>
#include <cstdint>
#include <unordered_map>
#include <set>

namespace {
// A/B switch for the langevin work-reduction gates.  All three are exact,
// so this only exists to make them measurable and bitwise-checkable against
// the code they replaced.
const bool tip_opt_gates_enabled = [] {
    const char* v = std::getenv("SARC_OPT_GATES");
    return !(v && std::string(v) == "0");
}();
}  // namespace


using vec = utils::vec;
// Define the global mutex.
std::mutex save_mutex;

//---------------------------------------------------------------------
// Constructor: Initializes the Langevin object and, if resume is true,
// loads the previous state and prints details. Otherwise, creates a new file.
//---------------------------------------------------------------------
Langevin::Langevin(Sarcomere& model0, double& beta0, double& dt0, double& D0_actin_trans,
    double& D0_actin_rot, double& D0_myosin_trans, double& D0_myosin_rot, int& save_every0, bool& resume, bool is3D,
    double max_actin_displacement0, double max_myosin_displacement0,
    double max_actin_rotation0, double max_myosin_rotation0,
    int resume_frame)
        : model(model0), beta(beta0), dt(dt0), D_actin_trans(D0_actin_trans),
        D_actin_rot(D0_actin_rot),D_myosin_trans(D0_myosin_trans), D_myosin_rot(D0_myosin_rot),
        save_every(save_every0), is3D(is3D),
        max_actin_displacement(max_actin_displacement0), max_myosin_displacement(max_myosin_displacement0),
        max_actin_rotation(max_actin_rotation0), max_myosin_rotation(max_myosin_rotation0),
        skip_initial_save(false),
        loaded_frame_index(-1)
{
    if (resume) {
        int n_frames;
        int requested_frame = resume_frame;
        int frame_idx = model.load_state(n_frames, requested_frame);
        loaded_frame_index = frame_idx;
        // Frames are written at step multiples of save_every before integrating that step.
        // Resume should continue from that same step index.
        start_step = static_cast<int>(model.current_step);
        skip_initial_save = true;
        if (requested_frame >= 0 && frame_idx != requested_frame) {
        }
        // Print myosin details.
        // for (int i = 0; i < model.myosin.n; i++) {
        //     printf("Myosin %d: %f %f %f\n", i, model.myosin.center[i].x, model.myosin.center[i].y, model.myosin.center[i].z);
        //     printf("Myosin endpoints: (%f %f %f), (%f %f %f)\n",
        //            model.myosin.left_end[i].x, model.myosin.left_end[i].y, model.myosin.left_end[i].z,
        //            model.myosin.right_end[i].x, model.myosin.right_end[i].y, model.myosin.right_end[i].z);
        // }
        // // Print actin details.
        // for (int i = 0; i < model.actin.n; i++) {
        //     printf("Actin %d: %f %f %f\n", i, model.actin.center[i].x, model.actin.center[i].y, model.actin.center[i].z);
        //     printf("Actin endpoints: (%f %f %f), (%f %f %f)\n",
        //            model.actin.left_end[i].x, model.actin.left_end[i].y, model.actin.left_end[i].z,
        //            model.actin.right_end[i].x, model.actin.right_end[i].y, model.actin.right_end[i].z);
        // }
    } else {
        start_step = 0;
        skip_initial_save = false;
        loaded_frame_index = -1;
        model.new_file();
    }
    // If resume mode active, interpret nsteps argument as total target
    // steps (i.e., inclusive of previously-run steps) when running.
    interpret_nsteps_as_total_when_resuming = resume;
}

//---------------------------------------------------------------------
// Destructor
//---------------------------------------------------------------------
Langevin::~Langevin() {
    // No dynamic resources need explicit deallocation.
}

void Langevin::set_rod_mobility(
    const rod_mobility::Diffusion& actin_mobility,
    const rod_mobility::Diffusion& myosin_mobility) {
    if (actin_mobility.enabled != myosin_mobility.enabled) {
        throw std::invalid_argument(
            "actin and myosin tensor mobility must be enabled together");
    }
    actin_rod_mobility = actin_mobility;
    myosin_rod_mobility = myosin_mobility;
    physical_rod_mobility = actin_mobility.enabled;
}

//---------------------------------------------------------------------
// run_langevin: Runs the simulation for nsteps, periodically saving the
// state and taking sample steps.
//---------------------------------------------------------------------
void Langevin::run_langevin(int nsteps, gsl_rng* rng, int& fix_myosin) {
    using ProfileClock = std::chrono::steady_clock;
    reset_runtime_profile();
    reset_actn2_tip_projection_profile();
    model.system_profile_steps = 0;
    model.system_neighbor_seconds = 0.0;
    model.system_parallel_seconds = 0.0;
    model.system_finalize_seconds = 0.0;
    model.system_total_seconds = 0.0;
    model.parallel_section_seconds.assign(
        model.parallel_section_seconds.size(), {});
    model.aa_phase_seconds.assign(model.aa_phase_seconds.size(), {});
    const auto loop_start = ProfileClock::now();
    auto print_profile = [&](double elapsed_seconds) {
        std::cout << "runtime_profile cumulative: steps="
                  << run_profile.steps
                  << ", elapsed_wall_s=" << elapsed_seconds
                  << ", update_system_wall_s="
                  << run_profile.update_system_seconds
                  << ", neighbor_wall_s=" << model.system_neighbor_seconds
                  << ", parallel_update_wall_s="
                  << model.system_parallel_seconds
                  << ", finalize_wall_s=" << model.system_finalize_seconds
                  << ", sample_wall_s=" << run_profile.sample_total_seconds
                  << ", noise_wall_s=" << run_profile.noise_seconds
                  << ", myosin_motion_wall_s="
                  << run_profile.myosin_motion_seconds
                  << ", actin_motion_wall_s="
                  << run_profile.actin_motion_seconds
                  << ", myosin_contact_wall_s="
                  << run_profile.myosin_contact_seconds
                  << ", actn2_contact_wall_s="
                  << run_profile.actn2_contact_seconds
                  << ", save_state_wall_s="
                  << run_profile.save_state_seconds
                  << ", save_resume_wall_s="
                  << run_profile.save_resume_seconds << std::endl;
        // Per-section breakdown of parallel_update_wall_s.  Timed per
        // thread without extra barriers, so each entry includes that
        // thread's wait at the section's implicit barrier.  "mean" is
        // the honest attribution; "max" approximates the critical path.
        if (!model.parallel_section_seconds.empty()) {
            const size_t threads = model.parallel_section_seconds.size();
            std::cout << "parallel_section_profile: threads=" << threads;
            for (int k = 0; k < Sarcomere::kProfileSectionCount; ++k) {
                double total = 0.0;
                double worst = 0.0;
                for (const auto& per_thread : model.parallel_section_seconds) {
                    total += per_thread[k];
                    worst = std::max(worst, per_thread[k]);
                }
                std::cout << ", " << Sarcomere::profile_section_name(k)
                          << "_mean_s=" << total / static_cast<double>(threads)
                          << ", " << Sarcomere::profile_section_name(k)
                          << "_max_s=" << worst;
            }
            std::cout << std::endl;
        }
        // Breakdown of actn2_chemistry.  That section scales 4.2x on 8 threads
        // where the others reach 6.8-7.9x, so this is where the serial regions
        // are; a phase whose mean approaches its max is running on one thread.
        if (!model.aa_phase_seconds.empty()) {
            const size_t threads = model.aa_phase_seconds.size();
            std::cout << "aa_phase_profile: threads=" << threads;
            for (int k = 0; k < Sarcomere::kAaPhaseCount; ++k) {
                double total = 0.0;
                double worst = 0.0;
                for (const auto& per_thread : model.aa_phase_seconds) {
                    total += per_thread[k];
                    worst = std::max(worst, per_thread[k]);
                }
                std::cout << ", " << Sarcomere::aa_phase_name(k)
                          << "_mean_s=" << total / static_cast<double>(threads)
                          << ", " << Sarcomere::aa_phase_name(k)
                          << "_max_s=" << worst;
            }
            std::cout << std::endl;
        }
    };

    double start = 0.0;
    double end = 0.0;
    bool skip_first_save = skip_initial_save;
    int end_step;
    if (interpret_nsteps_as_total_when_resuming) {
        // Here `nsteps` is interpreted as the absolute target total number of steps.
        end_step = nsteps;
    } else {
        end_step = start_step + nsteps;
    }
    for (int step = start_step; step < end_step; ++step) {
        bool should_save = (step % save_every == 0);
        if (skip_first_save && step == start_step) {
            should_save = false;
        }
        if (should_save) {
            std::cout << "Step " << step << std::endl;
            const auto save_start = ProfileClock::now();
            model.save_state();
            run_profile.save_state_seconds += std::chrono::duration<double>(
                ProfileClock::now() - save_start).count();
            start = omp_get_wtime();
        }
        const auto update_start = ProfileClock::now();
        model.update_system();
        run_profile.update_system_seconds += std::chrono::duration<double>(
            ProfileClock::now() - update_start).count();
        sample_step(dt, rng, fix_myosin);
        ++run_profile.steps;
        if (should_save) {
            const auto resume_start = ProfileClock::now();
            model.save_resume_snapshot();
            run_profile.save_resume_seconds += std::chrono::duration<double>(
                ProfileClock::now() - resume_start).count();
            end = omp_get_wtime();
            std::cout << "Step " << step << " took " << (end - start)
                      << " seconds" << std::endl;
            print_profile(std::chrono::duration<double>(
                ProfileClock::now() - loop_start).count());
        }
        if (skip_first_save && step == start_step) {
            skip_first_save = false;
        }
    }
    if (end_step > start_step && (end_step % save_every == 0)) {
        std::cout << "Step " << end_step << std::endl;
        const auto save_start = ProfileClock::now();
        model.save_state();
        run_profile.save_state_seconds += std::chrono::duration<double>(
            ProfileClock::now() - save_start).count();
    }
    skip_initial_save = false;
    start_step = end_step;
    const auto resume_start = ProfileClock::now();
    model.save_resume_snapshot();
    run_profile.save_resume_seconds += std::chrono::duration<double>(
        ProfileClock::now() - resume_start).count();
    run_profile.loop_total_seconds = std::chrono::duration<double>(
        ProfileClock::now() - loop_start).count();
    print_profile(run_profile.loop_total_seconds);
}

void Langevin::volume_exclusion(int nsteps, gsl_rng* rng, int& fix_myosin,
                                bool save_frames) {
    double start, end;
    bool skip_first_save = skip_initial_save;
    int end_step;
    if (interpret_nsteps_as_total_when_resuming) {
        end_step = nsteps;
    } else {
        end_step = start_step + nsteps;
    }
    for (int step = start_step; step < end_step; ++step) {
        bool should_save = save_frames && (step % save_every == 0);
        if (skip_first_save && step == start_step) {
            should_save = false;
        }
        if (should_save) {
            std::cout << "Step " << step << std::endl;
            // Optionally protect saving with the mutex:
            // std::lock_guard<std::mutex> lock(save_mutex);
            start = omp_get_wtime();
            model.save_state();
        }
        model.update_system_sterics_only();
        sample_step(dt, rng, fix_myosin);
        if (should_save) {
            model.save_resume_snapshot();
            end = omp_get_wtime();
            std::cout << "Step " << step << " took " << (end - start) << " seconds" << std::endl;
        }
        if (skip_first_save && step == start_step) {
            skip_first_save = false;
        }
    }
    if (save_frames && end_step > start_step && (end_step % save_every == 0)) {
        std::cout << "Step " << end_step << std::endl;
        model.save_state();
    }
    skip_initial_save = false;
    start_step = end_step;
    if (save_frames) model.save_resume_snapshot();
}

void Langevin::resolve_myosin_hard_contacts(int fix_myosin, int max_sweeps) {
    if (!model.sterics_enabled || model.myosin.n < 2 || max_sweeps <= 0) {
        return;
    }

    // Refresh the cached positions.  The Verlet skin preserves the pair set
    // between rebuilds; needs_rebuild() triggers before a filament can move
    // across that skin without entering the list.
    model.neighbor_list.set_species_positions(
        model.actin.center_x, model.actin.center_y, model.actin.center_z,
        model.myosin.center_x, model.myosin.center_y, model.myosin.center_z);
    if (model.neighbor_list.needs_rebuild()) {
        model.neighbor_list.rebuild_neighbor_list();
    }

    auto& contact_events =
        model.myosin["cumulative_hard_contact_events"];
    auto& contact_correction =
        model.myosin["cumulative_hard_contact_correction"];
    const double target = 2.0 * model.myosin.radius;

    // Sequential Gauss-Seidel sweeps are deterministic and conserve the pair
    // midpoint whenever both rods are mobile.  Repeated sweeps are needed
    // because resolving one contact in a dense cluster can close another.
    for (int sweep = 0; sweep < max_sweeps; ++sweep) {
        int corrected = 0;
        for (int i = 0; i < model.myosin.n; ++i) {
            for (const auto& [neighbor_index, neighbor_type] :
                 model.neighbor_list.get_neighbors(i + model.actin.n)) {
                if (neighbor_type != ParticleType::Myosin) continue;
                const int j = neighbor_index - model.actin.n;
                if (i >= j) continue;
                // myosin hard-contact gate
                if (tip_opt_gates_enabled &&
                    !geometry::segments_may_interact(
                        model.myosin.left_end[i], model.myosin.right_end[i],
                        model.myosin.left_end[j], model.myosin.right_end[j],
                        target, model.box, model.myosin.periodic_axes)) {
                    continue;
                }
                const vec old_direction_i = model.myosin.direction[i];
                const vec old_direction_j = model.myosin.direction[j];
                const double translation_weight = physical_rod_mobility
                    ? myosin_rod_mobility.perpendicular
                    : std::max(D_myosin_trans, 1e-12);
                const double rotation_weight = physical_rod_mobility
                    ? myosin_rod_mobility.tumble
                    : std::max(D_myosin_rot, 0.0);
                const auto projection = project_myosin_hard_contact(
                    model.myosin, i, j, model.box, fix_myosin, target,
                    translation_weight, rotation_weight);
                if (!projection.applied) continue;
                if (i >= fix_myosin &&
                    (static_cast<vec>(model.myosin.direction[i]) -
                     old_direction_i).norm() > 0.0) {
                    model.transport_myosin_material_frame(
                        i, old_direction_i, model.myosin.direction[i], 0.0);
                }
                if (j >= fix_myosin &&
                    (static_cast<vec>(model.myosin.direction[j]) -
                     old_direction_j).norm() > 0.0) {
                    model.transport_myosin_material_frame(
                        j, old_direction_j, model.myosin.direction[j], 0.0);
                }
                ++corrected;
                contact_events[i] += 1.0;
                contact_events[j] += 1.0;
                contact_correction[i] += 0.5 * projection.correction;
                contact_correction[j] += 0.5 * projection.correction;
            }
        }
        if (corrected == 0) break;
    }
}

void Langevin::resolve_actn2_tip_hard_contacts(
    int fix_myosin, int max_sweeps, double tolerance,
    Actn2TipProjectionTraversal traversal) {
    if (!model.aa_zdisc_tip_contact || !model.aa_zdisc_tip_hard_projection ||
        max_sweeps <= 0 ||
        model.aa_active_link_slots.empty() || model.myosin.n <= 0) {
        return;
    }
    using Clock = std::chrono::steady_clock;
    const auto total_start = Clock::now();
    ++actn2_tip_profile.calls;
    tolerance = std::max(1e-12, tolerance);
    // Keep a tiny positive numerical buffer. A later projection for a shared
    // actin or myosin can perturb an already-satisfied constraint; correcting
    // to exactly the audit tolerance made dense coupled contacts occasionally
    // miss that tolerance by roundoff after the final sweep.
    const double clearance = 100.0 * tolerance;
    auto& actin_events =
        model.actin["cumulative_aa_zdisc_tip_projection_events"];
    auto& actin_correction =
        model.actin["cumulative_aa_zdisc_tip_projection_correction"];
    auto& myosin_events =
        model.myosin["cumulative_aa_zdisc_tip_projection_events"];
    auto& myosin_correction =
        model.myosin["cumulative_aa_zdisc_tip_projection_correction"];

    auto mobility_response = [&](const vec& gradient, const vec& axis,
                                 bool actin_species) {
        vec response;
        if (physical_rod_mobility) {
            response = rod_mobility::apply_diffusion_tensor(
                gradient, axis,
                actin_species ? actin_rod_mobility
                              : myosin_rod_mobility);
        } else {
            response = (actin_species ? std::max(D_actin_trans, 1e-12)
                                      : std::max(D_myosin_trans, 1e-12)) *
                       gradient;
        }
        if (actin_species && constrain_actin_axial) {
            response -= axis * response.dot(axis);
        }
        return response;
    };

    struct LinkRecord {
        size_t slot = 0;
        int actin_i = -1;
        int actin_j = -1;
        size_t pair_index = 0;
    };
    struct ConstraintKey {
        size_t link_order = 0;
        int myosin = -1;
        int end_sign = 0;

        bool operator<(const ConstraintKey& other) const {
            if (link_order != other.link_order) {
                return link_order < other.link_order;
            }
            if (myosin != other.myosin) return myosin < other.myosin;
            return end_sign < other.end_sign;
        }
    };

    auto pair_key = [](int first, int second) {
        const auto low = static_cast<std::uint32_t>(std::min(first, second));
        const auto high = static_cast<std::uint32_t>(std::max(first, second));
        return (static_cast<std::uint64_t>(low) << 32) |
               static_cast<std::uint64_t>(high);
    };

    // Cache only the frozen neighbor-list union for each unique ACTN2-linked
    // actin pair. The spatial traversal never materializes the much larger
    // Cartesian product of links, candidate myosins, and two tip signs.
    const auto build_start = Clock::now();
    std::unordered_map<std::uint64_t, size_t> pair_indices;
    pair_indices.reserve(model.aa_active_link_slots.size());
    // Reused across calls: the inner vectors keep their capacity, so the
    // per-pair heap allocation disappears too.
    auto& pair_candidates = tip_pair_candidates;
    size_t pair_used = 0;
    if (tip_pair_union_stamp.size() !=
        static_cast<size_t>(model.myosin.n)) {
        tip_pair_union_stamp.assign(static_cast<size_t>(model.myosin.n), 0);
        tip_pair_union_token = 0;
    }
    std::vector<LinkRecord> links;
    links.reserve(model.aa_active_link_slots.size());
    std::vector<std::vector<size_t>> links_by_actin(
        static_cast<size_t>(model.actin.n));
    // Compact reverse incidence avoids scanning every active link whenever a
    // movable myosin is projected. It stores one link index per frozen
    // link/myosin pair, not two full constraint records.
    std::vector<std::vector<size_t>> links_by_myosin(
        static_cast<size_t>(model.myosin.n));
    std::uint64_t potential_constraints = 0;
    for (size_t slot : model.aa_active_link_slots) {
        if (slot >= model.aa_links.size()) continue;
        const auto& link = model.aa_links[slot];
        if (!link.active ||
            (model.aa_arm_state_model &&
             !actn2_arm_kinetics::is_bridge(link.arm_state))) {
            continue;
        }
        const auto key = pair_key(link.actin_i, link.actin_j);
        auto [pair_it, inserted] = pair_indices.try_emplace(key, pair_used);
        if (inserted) {
            if (pair_used == pair_candidates.size()) {
                pair_candidates.emplace_back();
            }
            auto& candidates = pair_candidates[pair_used];
            candidates.clear();
            if (!tip_opt_gates_enabled) {
                candidates = model.actin_neighbors_by_species[
                    static_cast<size_t>(link.actin_i)].second;
                const auto& other = model.actin_neighbors_by_species[
                    static_cast<size_t>(link.actin_j)].second;
                candidates.insert(candidates.end(), other.begin(), other.end());
                std::sort(candidates.begin(), candidates.end());
                candidates.erase(
                    std::unique(candidates.begin(), candidates.end()),
                    candidates.end());
                ++pair_used;
            } else {
            // Stamped union.  Marking both lists and then scanning the myosin
            // range yields the ascending, duplicate-free sequence that
            // sort+unique produced -- identical contents in identical order,
            // which the Gauss-Seidel sweep below depends on -- in O(n_i + n_j
            // + n_myosin) instead of O(n log n), with no allocation.
            ++tip_pair_union_token;
            const auto token = tip_pair_union_token;
            for (int side = 0; side < 2; ++side) {
                for (int myosin_index : model.actin_neighbors_by_species[
                        static_cast<size_t>(side == 0 ? link.actin_i
                                                      : link.actin_j)].second) {
                    if (myosin_index >= 0 && myosin_index < model.myosin.n) {
                        tip_pair_union_stamp[
                            static_cast<size_t>(myosin_index)] = token;
                    }
                }
            }
            for (int myosin_index = 0; myosin_index < model.myosin.n;
                 ++myosin_index) {
                if (tip_pair_union_stamp[static_cast<size_t>(myosin_index)] ==
                    token) {
                    candidates.push_back(myosin_index);
                }
            }
            ++pair_used;
            }
        }
        const size_t link_order = links.size();
        links.push_back(
            {slot, link.actin_i, link.actin_j, pair_it->second});
        links_by_actin[static_cast<size_t>(link.actin_i)].push_back(
            link_order);
        if (link.actin_j != link.actin_i) {
            links_by_actin[static_cast<size_t>(link.actin_j)].push_back(
                link_order);
        }
        for (int myosin_index :
             pair_candidates[pair_it->second]) {
            links_by_myosin[static_cast<size_t>(myosin_index)].push_back(
                link_order);
        }
        potential_constraints +=
            2ULL * pair_candidates[pair_it->second].size();
    }
    actn2_tip_profile.active_links += links.size();
    actn2_tip_profile.unique_actin_pairs += pair_used;
    actn2_tip_profile.constraints += potential_constraints;
    actn2_tip_profile.build_seconds += std::chrono::duration<double>(
        Clock::now() - build_start).count();

    auto patch_for_link = [&](size_t link_order) {
        const auto& record = links[link_order];
        const auto& link = model.aa_links[record.slot];
        return compute_actn2_patch_geometry(
            model.actin, record.actin_i, record.actin_j,
            link.s_i, link.s_j, model.box);
    };

    auto broadphase_candidate = [&](const ACTN2PatchGeometry& patch,
                                    int myosin_index) {
        // The tip sign does not affect the conservative enclosing sphere, so
        // one geometric test covers the two canonical tip constraints.
        actn2_tip_profile.broadphase_tests += 2;
        const bool candidate = actn2_patch_may_contact_myosin_tip(
            patch, model.myosin, myosin_index, model.box,
            model.aa_zdisc_tip_onset,
            model.aa_zdisc_tip_transverse_cutoff);
        if (candidate) actn2_tip_profile.broadphase_candidates += 2;
        return candidate;
    };

    auto evaluate_and_project = [&](const ConstraintKey& key) {
        const auto& record = links[key.link_order];
        const auto& link = model.aa_links[record.slot];
        // Cheap predicate first.  `applied && axial_gap < -tolerance` is
        // equivalent to `eligible && axial_gap < -tolerance`: axial_gap below
        // -tolerance is also below onset, so compression and hence tanh(x) are
        // positive, and radial_gate is positive whenever transverse_distance is
        // strictly inside the cutoff -- which `eligible` already requires.  So
        // the force is positive and finite exactly when the gap test passes,
        // and none of the smoothstep/tanh/energy/torque work is needed to
        // decide.  It fires ~2 times per step out of ~53,000 evaluations.
        const auto gap = zdisc_myosin_tip_gap(
            model.actin, model.myosin, record.actin_i, record.actin_j,
            link.s_i, link.s_j, key.myosin, key.end_sign, model.box,
            model.aa_zdisc_tip_onset,
            model.aa_zdisc_tip_transverse_cutoff);
        if (tip_opt_gates_enabled &&
            (!gap.eligible || gap.axial_gap >= -tolerance)) {
            return false;
        }
        const auto contact = compute_zdisc_myosin_tip_contact(
            model.actin, model.myosin, record.actin_i,
            record.actin_j, link.s_i, link.s_j, key.myosin,
            key.end_sign, model.box, model.aa_zdisc_tip_onset,
            model.aa_zdisc_tip_stiffness, model.aa_zdisc_tip_force_max,
            model.aa_zdisc_tip_transverse_flat,
            model.aa_zdisc_tip_transverse_cutoff);
        if (!contact.applied || contact.axial_gap >= -tolerance) {
            return false;
        }
        vec outward = model.myosin.direction[key.myosin];
        outward.normalize();
        outward = static_cast<double>(key.end_sign) * outward;
        vec axis_i = model.actin.direction[record.actin_i];
        vec axis_j = model.actin.direction[record.actin_j];
        vec axis_m = model.myosin.direction[key.myosin];
        axis_i.normalize();
        axis_j.normalize();
        axis_m.normalize();
        const vec gradient_i = 0.5 * outward;
        const vec gradient_j = 0.5 * outward;
        const vec gradient_m = -outward;
        const vec response_i = mobility_response(gradient_i, axis_i, true);
        const vec response_j = mobility_response(gradient_j, axis_j, true);
        const bool myosin_fixed = key.myosin < fix_myosin;
        const vec response_m = myosin_fixed
            ? vec{0.0, 0.0, 0.0}
            : mobility_response(gradient_m, axis_m, false);
        const double denominator =
            gradient_i.dot(response_i) + gradient_j.dot(response_j) +
            gradient_m.dot(response_m);
        if (!(denominator > 0.0) || !std::isfinite(denominator)) {
            throw std::runtime_error(
                "ACTN2 tip hard contact has zero mobility");
        }
        const double required = -contact.axial_gap + clearance;
        const double multiplier = required / denominator;
        const vec delta_i = multiplier * response_i;
        const vec delta_j = multiplier * response_j;
        const vec delta_m = multiplier * response_m;
        auto displace = [](Filament& filament, int index,
                           const vec& delta) {
            double dx = delta.x, dy = delta.y, dz = delta.z;
            filament.displace(index, dx, dy, dz);
        };
        displace(model.actin, record.actin_i, delta_i);
        displace(model.actin, record.actin_j, delta_j);
        if (!myosin_fixed) displace(model.myosin, key.myosin, delta_m);
        actin_events[record.actin_i] += 1.0;
        actin_events[record.actin_j] += 1.0;
        myosin_events[key.myosin] += 1.0;
        actin_correction[record.actin_i] += 0.5 * required;
        actin_correction[record.actin_j] += 0.5 * required;
        myosin_correction[key.myosin] += required;
        return true;
    };

    auto audit_constraint = [&](const ConstraintKey& key, double& worst_gap,
                                std::uint64_t& worst_link,
                                int& worst_myosin, int& worst_sign) {
        ++actn2_tip_profile.audit_evaluations;
        const auto& record = links[key.link_order];
        const auto& link = model.aa_links[record.slot];
        const auto contact = compute_zdisc_myosin_tip_contact(
            model.actin, model.myosin, record.actin_i,
            record.actin_j, link.s_i, link.s_j, key.myosin,
            key.end_sign, model.box, model.aa_zdisc_tip_onset,
            model.aa_zdisc_tip_stiffness, model.aa_zdisc_tip_force_max,
            model.aa_zdisc_tip_transverse_flat,
            model.aa_zdisc_tip_transverse_cutoff);
        if (contact.applied && contact.axial_gap < worst_gap) {
            worst_gap = contact.axial_gap;
            worst_link = link.id;
            worst_myosin = key.myosin;
            worst_sign = key.end_sign;
        }
    };

    // Mobility-weighted projected Gauss-Seidel is the overdamped hard-contact
    // limit. The spatial traversal starts only with contacts inside an exact
    // conservative bound. Whenever a correction moves a body, all constraints
    // that body could newly activate are rediscovered dynamically.
    const auto solve_start = Clock::now();
    // Initial discovery is already produced in canonical order, so keep it in
    // a compact vector. Only constraints activated by a correction need the
    // logarithmic ordered-set insertion path.
    std::vector<ConstraintKey> current_constraints;
    std::set<ConstraintKey> current_injected;
    std::set<ConstraintKey> next_constraints;
    if (traversal != Actn2TipProjectionTraversal::full_scan) {
        // Initial discovery is a pure read-only filter -- ~1.1M sphere tests
        // per step at production density, and it was the single largest serial
        // block in the model.  Each link writes its survivors to its own
        // bucket and the buckets are concatenated in link order afterwards, so
        // the constraint sequence handed to the sequential Gauss-Seidel sweep
        // below is byte-for-byte the one the serial loop produced.
        if (traversal == Actn2TipProjectionTraversal::spatial_incidence &&
            tip_opt_gates_enabled) {
            if (tip_discovery_scratch.size() < links.size()) {
                tip_discovery_scratch.resize(links.size());
            }
            std::uint64_t tests = 0, hits = 0;
            #pragma omp parallel for schedule(dynamic, 8) \
                reduction(+ : tests, hits)
            for (long long lo = 0; lo < static_cast<long long>(links.size());
                 ++lo) {
                auto& keep = tip_discovery_scratch[static_cast<size_t>(lo)];
                keep.clear();
                const auto& candidates =
                    pair_candidates[links[static_cast<size_t>(lo)].pair_index];
                const auto patch = patch_for_link(static_cast<size_t>(lo));
                for (int myosin_index : candidates) {
                    tests += 2;
                    if (!actn2_patch_may_contact_myosin_tip(
                            patch, model.myosin, myosin_index, model.box,
                            model.aa_zdisc_tip_onset,
                            model.aa_zdisc_tip_transverse_cutoff)) {
                        continue;
                    }
                    hits += 2;
                    keep.push_back(myosin_index);
                }
            }
            actn2_tip_profile.broadphase_tests += tests;
            actn2_tip_profile.broadphase_candidates += hits;
            for (size_t link_order = 0; link_order < links.size();
                 ++link_order) {
                for (int myosin_index : tip_discovery_scratch[link_order]) {
                    current_constraints.push_back(
                        {link_order, myosin_index, -1});
                    current_constraints.push_back(
                        {link_order, myosin_index, 1});
                }
            }
        }
        for (size_t link_order = 0; link_order < links.size();
             ++link_order) {
            const auto& candidates =
                pair_candidates[links[link_order].pair_index];
            if (traversal ==
                Actn2TipProjectionTraversal::spatial_incidence) {
                if (tip_opt_gates_enabled) {
                    (void)candidates;   // handled by the parallel pass above
                    continue;
                }
                const auto patch = patch_for_link(link_order);
                for (int myosin_index : candidates) {
                    if (!broadphase_candidate(patch, myosin_index)) continue;
                    current_constraints.push_back(
                        {link_order, myosin_index, -1});
                    current_constraints.push_back(
                        {link_order, myosin_index, 1});
                }
            } else {
                for (int myosin_index : candidates) {
                    current_constraints.push_back(
                        {link_order, myosin_index, -1});
                    current_constraints.push_back(
                        {link_order, myosin_index, 1});
                }
            }
        }
    }

    std::vector<std::uint32_t> discovery_generation(links.size(), 0);
    std::uint32_t generation = 0;
    int last_corrected = 0;
    for (int sweep = 0; sweep < max_sweeps; ++sweep) {
        int corrected = 0;
        ++actn2_tip_profile.sweeps;

        auto process_correction = [&](const ConstraintKey& corrected_key) {
            const auto& corrected_record = links[corrected_key.link_order];
            ++generation;
            if (generation == 0) {
                std::fill(discovery_generation.begin(),
                          discovery_generation.end(), 0);
                generation = 1;
            }
            auto enqueue = [&](const ConstraintKey& key) {
                if (corrected_key < key) {
                    current_injected.insert(key);
                } else {
                    next_constraints.insert(key);
                }
            };
            auto discover_entire_link = [&](size_t link_order) {
                if (discovery_generation[link_order] == generation) return;
                discovery_generation[link_order] = generation;
                const auto& candidates =
                    pair_candidates[links[link_order].pair_index];
                if (traversal ==
                    Actn2TipProjectionTraversal::spatial_incidence) {
                    const auto patch = patch_for_link(link_order);
                    for (int myosin_index : candidates) {
                        if (!broadphase_candidate(patch, myosin_index)) {
                            continue;
                        }
                        enqueue({link_order, myosin_index, -1});
                        enqueue({link_order, myosin_index, 1});
                    }
                } else {
                    for (int myosin_index : candidates) {
                        enqueue({link_order, myosin_index, -1});
                        enqueue({link_order, myosin_index, 1});
                    }
                }
            };
            for (size_t link_order :
                 links_by_actin[
                     static_cast<size_t>(corrected_record.actin_i)]) {
                discover_entire_link(link_order);
            }
            for (size_t link_order :
                 links_by_actin[
                     static_cast<size_t>(corrected_record.actin_j)]) {
                discover_entire_link(link_order);
            }
            if (corrected_key.myosin < fix_myosin) return;

            // A moved myosin can activate a contact with any incident link.
            // Links already scanned because an endpoint actin moved need no
            // second pass.
            for (size_t link_order :
                 links_by_myosin[
                     static_cast<size_t>(corrected_key.myosin)]) {
                if (discovery_generation[link_order] == generation) continue;
                if (traversal ==
                    Actn2TipProjectionTraversal::spatial_incidence) {
                    const auto patch = patch_for_link(link_order);
                    if (!broadphase_candidate(
                            patch, corrected_key.myosin)) {
                        continue;
                    }
                }
                enqueue({link_order, corrected_key.myosin, -1});
                enqueue({link_order, corrected_key.myosin, 1});
            }
        };

        if (traversal == Actn2TipProjectionTraversal::full_scan) {
            for (size_t link_order = 0; link_order < links.size();
                 ++link_order) {
                const auto& candidates =
                    pair_candidates[links[link_order].pair_index];
                for (int myosin_index : candidates) {
                    for (int end_sign : {-1, 1}) {
                        ++actn2_tip_profile.solve_evaluations;
                        if (!evaluate_and_project(
                                {link_order, myosin_index, end_sign})) {
                            continue;
                        }
                        ++corrected;
                        ++actn2_tip_profile.corrections;
                    }
                }
            }
        } else {
            size_t current_index = 0;
            while (current_index < current_constraints.size() ||
                   !current_injected.empty()) {
                ConstraintKey key;
                if (current_index == current_constraints.size()) {
                    const auto iterator = current_injected.begin();
                    key = *iterator;
                    current_injected.erase(iterator);
                } else if (current_injected.empty() ||
                           current_constraints[current_index] <
                               *current_injected.begin()) {
                    key = current_constraints[current_index++];
                } else if (*current_injected.begin() <
                           current_constraints[current_index]) {
                    const auto iterator = current_injected.begin();
                    key = *iterator;
                    current_injected.erase(iterator);
                } else {
                    // The dynamically rediscovered key was already pending in
                    // the compact initial vector; evaluate it only once.
                    key = current_constraints[current_index++];
                    current_injected.erase(current_injected.begin());
                }
                ++actn2_tip_profile.solve_evaluations;
                if (!evaluate_and_project(key)) continue;
                ++corrected;
                ++actn2_tip_profile.corrections;
                process_correction(key);
            }
        }

        last_corrected = corrected;
        if (corrected == 0) break;
        if (traversal != Actn2TipProjectionTraversal::full_scan) {
            current_constraints.assign(next_constraints.begin(),
                                       next_constraints.end());
            next_constraints.clear();
            current_injected.clear();
        }
    }
    if (last_corrected != 0) ++actn2_tip_profile.exhausted_calls;
    actn2_tip_profile.solve_seconds += std::chrono::duration<double>(
        Clock::now() - solve_start).count();

    // A zero-correction incidence sweep is a complete certificate: the initial
    // broad phase is conservative, and every constraint sharing a subsequently
    // moved body is rediscovered. Only an exhausted solve needs the diagnostic
    // full audit. Full-scan mode retains the audit as the parity oracle.
    if (traversal != Actn2TipProjectionTraversal::full_scan &&
        last_corrected == 0) {
        ++actn2_tip_profile.certified_without_audit_calls;
        actn2_tip_profile.total_seconds += std::chrono::duration<double>(
            Clock::now() - total_start).count();
        return;
    }

    ++actn2_tip_profile.audited_calls;
    const auto audit_start = Clock::now();
    double worst_gap = 0.0;
    std::uint64_t worst_link = 0;
    int worst_myosin = -1;
    int worst_sign = 0;
    for (size_t link_order = 0; link_order < links.size(); ++link_order) {
        const auto& candidates =
            pair_candidates[links[link_order].pair_index];
        for (int myosin_index : candidates) {
            audit_constraint({link_order, myosin_index, -1}, worst_gap,
                             worst_link, worst_myosin, worst_sign);
            audit_constraint({link_order, myosin_index, 1}, worst_gap,
                             worst_link, worst_myosin, worst_sign);
        }
    }
    actn2_tip_profile.audit_seconds += std::chrono::duration<double>(
        Clock::now() - audit_start).count();

    constexpr double audit_tolerance = 1e-6;  // um = 1 pm
    if (worst_gap < -audit_tolerance) {
        std::fprintf(
            stderr,
            "WARNING: ACTN2 tip nonpenetration audit exceeded: "
            "link=%llu myosin=%d tip=%d gap_um=%.17g gap_pm=%.9g "
            "audit_tolerance_pm=%.9g; continuing\n",
            static_cast<unsigned long long>(worst_link), worst_myosin,
            worst_sign, worst_gap, 1e6 * worst_gap,
            1e6 * audit_tolerance);
    }
    actn2_tip_profile.total_seconds += std::chrono::duration<double>(
        Clock::now() - total_start).count();
}
//---------------------------------------------------------------------
// sample_step: Performs a single Langevin dynamics step by updating the 
// system, generating noise, and displacing myosin and actin particles.
//---------------------------------------------------------------------
void Langevin::sample_step(double& dt, gsl_rng* rng, int& fix_myosin) {
    using ProfileClock = std::chrono::steady_clock;
    const auto sample_start = ProfileClock::now();
    // Generate noise for both myosin and actin particles.
    int n_randns = (model.myosin.n + model.actin.n) * 6;
    std::vector<double> noise(n_randns);
    for (int i = 0; i < n_randns; i++) {
        noise[i] = thermal_noise_scale * gsl_ran_gaussian(rng, 1.0);
    }

    int n_acc_randns = model.myosin.n + model.actin.n;
    std::vector<double> acc_rand(n_acc_randns);
    for (int i = 0; i < n_acc_randns; i++) {
        acc_rand[i] = gsl_rng_uniform(rng);
    }
    const auto noise_end = ProfileClock::now();
    run_profile.noise_seconds += std::chrono::duration<double>(
        noise_end - sample_start).count();

    const double wall_eps = 1e-9;
    auto reflect_segment = [&](Filament& filament, int idx, vec& delta_pos, vec& delta_dir) {
        if (model.is_periodic[0] && model.is_periodic[1] && model.is_periodic[2]) {
            return;
        }
        int guard = 0;
        while (guard < 6) {
            vec center_current = filament.center[idx];
            vec direction_current = filament.direction[idx];
            vec center_prop = center_current + delta_pos;
            vec direction_prop = direction_current + delta_dir;
            if (!is3D) {
                direction_prop.z = 0.0;
            }
            double prop_norm = direction_prop.norm();
            if (prop_norm <= 1e-12) {
                direction_prop = direction_current;
            } else {
                direction_prop = direction_prop / prop_norm;
            }
            vec dir_norm = direction_prop;
            dir_norm.normalize();
            vec left_prop = center_prop - 0.5 * filament.length * dir_norm;
            vec right_prop = center_prop + 0.5 * filament.length * dir_norm;
            bool reflected = false;
            for (int axis = 0; axis < 3; ++axis) {
                if (model.is_periodic[axis]) {
                    continue;
                }
                if (axis >= static_cast<int>(model.box.size()) || model.box[axis] <= 0.0) {
                    continue;
                }
                double half = 0.5 * model.box[axis];
                double lower = -half;
                double upper = half;
                auto component = [&](const vec& v) -> double {
                    return axis == 0 ? v.x : (axis == 1 ? v.y : v.z);
                };
                double left_c = component(left_prop);
                double right_c = component(right_prop);
                if (left_c < lower - wall_eps || right_c < lower - wall_eps ||
                    left_c > upper + wall_eps || right_c > upper + wall_eps) {
                    if (axis == 0) {
                        delta_pos.x *= -1.0;
                        delta_dir.x *= -1.0;
                    } else if (axis == 1) {
                        delta_pos.y *= -1.0;
                        delta_dir.y *= -1.0;
                    } else {
                        delta_pos.z *= -1.0;
                        delta_dir.z *= -1.0;
                    }
                    reflected = true;
                }
            }
            if (!reflected) {
                break;
            }
            ++guard;
        }
    };

    auto clamp_center = [&](Filament& filament, int idx) {
        if (model.is_periodic[0] && model.is_periodic[1] && model.is_periodic[2]) {
            return;
        }
        vec center = filament.center[idx];
        vec dir = filament.direction[idx];
        for (int axis = 0; axis < 3; ++axis) {
            if (model.is_periodic[axis]) {
                continue;
            }
            if (axis >= static_cast<int>(model.box.size()) || model.box[axis] <= 0.0) {
                continue;
            }
            double half = 0.5 * model.box[axis];
            double lower = -half;
            double upper = half;
            double dir_component = axis == 0 ? dir.x : (axis == 1 ? dir.y : dir.z);
            double half_span = 0.5 * filament.length * std::abs(dir_component);
            double min_center = lower + half_span + wall_eps;
            double max_center = upper - half_span - wall_eps;
            double* center_component = (axis == 0) ? &center.x : (axis == 1 ? &center.y : &center.z);
            if (min_center > max_center) {
                *center_component = std::clamp(*center_component, lower + wall_eps, upper - wall_eps);
            } else {
                *center_component = std::clamp(*center_component, min_center, max_center);
            }
        }
        filament.center[idx] = center;
        filament.update_endpoints(idx);
    };

    int offset = model.myosin.n * 6;
    double D = D_myosin_trans;
    double D_rot = D_myosin_rot;
    // Update myosin particles.
    const double displacement_slack = 1.3;
    const auto myosin_motion_start = ProfileClock::now();
    for (int i = fix_myosin; i < model.myosin.n; i++) {
        if (!is3D) {
            model.myosin.force[i].z = 0;
            model.myosin.velocity[i].z = 0;
            model.myosin.torque[i].z = 0;
        }
        vec translation_noise{noise[i * 6], noise[i * 6 + 1],
                              is3D ? noise[i * 6 + 2] : 0.0};
        vec delta_pos;
        vec myosin_axis = model.myosin.direction[i];
        myosin_axis.normalize();
        if (physical_rod_mobility) {
            delta_pos = beta * dt * rod_mobility::apply_diffusion_tensor(
                model.myosin.force[i], myosin_axis, myosin_rod_mobility);
            delta_pos += dt * model.myosin.velocity[i];
            delta_pos += std::sqrt(2.0 * dt) *
                rod_mobility::apply_sqrt_diffusion_tensor(
                    translation_noise, myosin_axis, myosin_rod_mobility);
        } else {
            delta_pos = beta * D * dt * model.myosin.force[i] +
                        dt * model.myosin.velocity[i] +
                        std::sqrt(2.0 * D * dt) * translation_noise;
        }
        double disp_sq = delta_pos.x * delta_pos.x + delta_pos.y * delta_pos.y +
                         (is3D ? delta_pos.z * delta_pos.z : 0.0);
        double disp_limit = max_myosin_displacement;
        if (disp_limit > 0.0 && std::isfinite(disp_limit)) {
            double disp_mag = std::sqrt(disp_sq);
            double allowed = displacement_slack * disp_limit;
            if (disp_mag > allowed && disp_mag > 1e-12) {
                double scale = allowed / disp_mag;
                delta_pos.x *= scale;
                delta_pos.y *= scale;
                if (is3D) {
                    delta_pos.z *= scale;
                }
            }
        }
        vec rot_noise={noise[i * 6 + 3], noise[i * 6 + 4], is3D ? noise[i * 6 + 5] : 0.0};
        vec delta_u;
        double axial_roll = 0.0;
        if (physical_rod_mobility) {
            if (is3D) {
                const vec tangent_drive = rod_mobility::project_perpendicular(
                    model.myosin.torque[i], myosin_axis);
                const vec tangent_noise = rod_mobility::project_perpendicular(
                    rot_noise, myosin_axis);
                delta_u = beta * myosin_rod_mobility.tumble * dt * tangent_drive +
                    std::sqrt(2.0 * myosin_rod_mobility.tumble * dt) * tangent_noise;
                axial_roll = beta * myosin_rod_mobility.roll * dt *
                                 model.myosin.torque[i].dot(myosin_axis) +
                             std::sqrt(2.0 * myosin_rod_mobility.roll * dt) *
                                 rot_noise.dot(myosin_axis);
            } else {
                const vec tangent{-myosin_axis.y, myosin_axis.x, 0.0};
                const double scalar = beta * myosin_rod_mobility.tumble * dt *
                                          model.myosin.torque[i].dot(tangent) +
                                      std::sqrt(2.0 * myosin_rod_mobility.tumble * dt) *
                                          rot_noise.dot(tangent);
                delta_u = scalar * tangent;
            }
        } else {
            delta_u = std::sqrt(2 * D_rot * dt) * rot_noise +
                      dt * model.myosin.torque[i] * D_rot * beta;
            axial_roll = delta_u.dot(myosin_axis);
        }
        if (!is3D) {
            delta_pos.z = 0.0;
            delta_u.z = 0.0;
        }
        double rot_limit = max_myosin_rotation;
        const double rotation_scale = rod_mobility::angular_increment_scale(
            delta_u, rot_limit, displacement_slack);
        if (rotation_scale < 1.0) {
            delta_u = rotation_scale * delta_u;
            // A large tangent step and material-frame roll belong to the same
            // realized Brownian/mechanical update. Keep their time increment
            // coherent whenever the director trust radius is activated.
            axial_roll *= rotation_scale;
            model.myosin["cumulative_rotation_cap_events"][i] += 1.0;
        }
        reflect_segment(model.myosin, i, delta_pos, delta_u);

        double dx = delta_pos.x;
        double dy = delta_pos.y;
        double dz = delta_pos.z;
        model.myosin.displace(i, dx, dy, dz);
        //print myosin force and displacement
        // printf("Myosin %d force: (%f, %f, %f), displacement: (%f, %f, %f)\n",
        //        i, model.myosin.force[i].x, model.myosin.force[i].y, model.myosin.force[i].z,
        //        dx, dy, dz);
        const vec old_dir = model.myosin.direction[i];
        vec new_dir = physical_rod_mobility
            ? rod_mobility::sphere_exp(myosin_axis, delta_u)
            : old_dir + delta_u;
        if (!is3D) {
            new_dir.z = 0.0;
            new_dir.normalize();
        } else {
            new_dir.normalize();
        }
        model.myosin.direction[i] = new_dir;
        model.transport_myosin_material_frame(
            i, old_dir, new_dir, axial_roll);
        model.myosin.update_endpoints(i);
        clamp_center(model.myosin, i);
    }
    const auto actin_motion_start = ProfileClock::now();
    run_profile.myosin_motion_seconds += std::chrono::duration<double>(
        actin_motion_start - myosin_motion_start).count();
    // Update actin particles.
    for (int i = 0; i < model.actin.n; i++) {
        // if (model.actin.cb_status[i] > 1){
        //     D = D_myosin_trans;
        //     D_rot = D_myosin_rot;
        // }
        // else{
        //     D = D_actin_trans;
        //     D_rot = D_actin_rot;
        // }
        D = D_actin_trans;
        D_rot = D_actin_rot;
        if (!is3D) {
            model.actin.force[i].z = 0;
            model.actin.velocity[i].z = 0;
            model.actin.torque[i].z = 0;
        }
        vec translation_noise{noise[offset + i * 6], noise[offset + i * 6 + 1],
                              is3D ? noise[offset + i * 6 + 2] : 0.0};
        vec delta_pos;
        vec actin_axis = model.actin.direction[i];
        actin_axis.normalize();
        if (physical_rod_mobility) {
            delta_pos = beta * dt * rod_mobility::apply_diffusion_tensor(
                model.actin.force[i], actin_axis, actin_rod_mobility);
            delta_pos += dt * model.actin.velocity[i];
            delta_pos += std::sqrt(2.0 * dt) *
                rod_mobility::apply_sqrt_diffusion_tensor(
                    translation_noise, actin_axis, actin_rod_mobility);
        } else {
            delta_pos = beta * D * dt * model.actin.force[i] +
                        dt * model.actin.velocity[i] +
                        std::sqrt(2.0 * D * dt) * translation_noise;
        }
        if (constrain_actin_axial) {
            delta_pos -= actin_axis * delta_pos.dot(actin_axis);
        }
        double disp_sq = delta_pos.x * delta_pos.x + delta_pos.y * delta_pos.y +
                         (is3D ? delta_pos.z * delta_pos.z : 0.0);
        double disp_limit = max_actin_displacement;
        if (disp_limit > 0.0 && std::isfinite(disp_limit)) {
            double disp_mag = std::sqrt(disp_sq);
            double allowed = displacement_slack * disp_limit;
            if (disp_mag > allowed && disp_mag > 1e-12) {
                double scale = allowed / disp_mag;
                delta_pos.x *= scale;
                delta_pos.y *= scale;
                if (is3D) {
                    delta_pos.z *= scale;
                }
            }
        }
        vec rot_noise={noise[offset + i * 6 + 3], noise[offset + i * 6 + 4], is3D ? noise[offset + i * 6 + 5] : 0.0};
        vec delta_u;
        if (physical_rod_mobility) {
            if (is3D) {
                delta_u = beta * actin_rod_mobility.tumble * dt *
                              rod_mobility::project_perpendicular(
                                  model.actin.torque[i], actin_axis) +
                          std::sqrt(2.0 * actin_rod_mobility.tumble * dt) *
                              rod_mobility::project_perpendicular(rot_noise, actin_axis);
            } else {
                const vec tangent{-actin_axis.y, actin_axis.x, 0.0};
                const double scalar = beta * actin_rod_mobility.tumble * dt *
                                          model.actin.torque[i].dot(tangent) +
                                      std::sqrt(2.0 * actin_rod_mobility.tumble * dt) *
                                          rot_noise.dot(tangent);
                delta_u = scalar * tangent;
            }
        } else {
            delta_u = std::sqrt(2 * D_rot * dt) * rot_noise +
                      dt * model.actin.torque[i] * D_rot * beta;
        }
        if (!is3D) {
            delta_pos.z = 0.0;
            delta_u.z = 0.0;
        }
        double rot_limit = max_actin_rotation;
        const double rotation_scale = rod_mobility::angular_increment_scale(
            delta_u, rot_limit, displacement_slack);
        if (rotation_scale < 1.0) {
            delta_u = rotation_scale * delta_u;
            model.actin["cumulative_rotation_cap_events"][i] += 1.0;
        }
        reflect_segment(model.actin, i, delta_pos, delta_u);

        double dx = delta_pos.x;
        double dy = delta_pos.y;
        double dz = delta_pos.z;
        model.actin.displace(i, dx, dy, dz);
        vec new_dir_act = physical_rod_mobility
            ? rod_mobility::sphere_exp(actin_axis, delta_u)
            : static_cast<vec>(model.actin.direction[i]) + delta_u;
        if (!is3D) {
            new_dir_act.z = 0.0;
            new_dir_act.normalize();
        } else {
            new_dir_act.normalize();
        }
        model.actin.direction[i] = new_dir_act;
        model.actin.update_endpoints(i);
        clamp_center(model.actin, i);
    }
    const auto myosin_contact_start = ProfileClock::now();
    run_profile.actin_motion_seconds += std::chrono::duration<double>(
        myosin_contact_start - actin_motion_start).count();
    resolve_myosin_hard_contacts(fix_myosin);
    const auto actn2_contact_start = ProfileClock::now();
    run_profile.myosin_contact_seconds += std::chrono::duration<double>(
        actn2_contact_start - myosin_contact_start).count();
    // Production always uses the exact conservative spatial broad phase.
    // Slower traversals remain callable only by parity/performance fixtures.
    resolve_actn2_tip_hard_contacts(
        fix_myosin, 64, 1e-9,
        Actn2TipProjectionTraversal::spatial_incidence);
    const auto sample_end = ProfileClock::now();
    run_profile.actn2_contact_seconds += std::chrono::duration<double>(
        sample_end - actn2_contact_start).count();
    run_profile.sample_total_seconds += std::chrono::duration<double>(
        sample_end - sample_start).count();
}
