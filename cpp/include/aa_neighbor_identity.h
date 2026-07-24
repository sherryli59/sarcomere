#ifndef AA_NEIGHBOR_IDENTITY_H
#define AA_NEIGHBOR_IDENTITY_H

#include <algorithm>
#include <cstdint>
#include <limits>
#include <unordered_map>
#include <vector>

namespace aa_neighbor_identity {

using Step = std::uint64_t;

inline std::uint64_t pair_key(int first, int second) {
    const auto low = static_cast<std::uint32_t>(std::min(first, second));
    const auto high = static_cast<std::uint32_t>(std::max(first, second));
    return (static_cast<std::uint64_t>(low) << 32) |
           static_cast<std::uint64_t>(high);
}

inline Step saturating_add(Step first, Step second) {
    const Step maximum = std::numeric_limits<Step>::max();
    return second > maximum - first ? maximum : first + second;
}

struct Slot {
    int first = -1;
    int second = -1;
    Step last_active_step = 0;
    Step protected_until = 0;
    bool pinned = false;
};

class Registry {
public:
    void clear() { slots_.clear(); std::fill(degree_.begin(), degree_.end(), 0); }

    bool contains(int first, int second) const {
        return slots_.count(pair_key(first, second)) != 0;
    }

    int degree(int actin_index) const {
        return (actin_index >= 0 &&
                actin_index < static_cast<int>(degree_.size()))
            ? degree_[static_cast<std::size_t>(actin_index)] : 0;
    }

    const Slot* find(int first, int second) const {
        const auto iterator = slots_.find(pair_key(first, second));
        return iterator == slots_.end() ? nullptr : &iterator->second;
    }

    bool reserve(int first, int second, Step current_step,
                 int maximum_degree, bool pinned = false) {
        if (first == second || first < 0 || second < 0) return false;
        if (first > second) std::swap(first, second);
        const auto key = pair_key(first, second);
        auto iterator = slots_.find(key);
        if (iterator != slots_.end()) {
            iterator->second.pinned = iterator->second.pinned || pinned;
            return true;
        }
        if (maximum_degree > 0 &&
            (degree(first) >= maximum_degree ||
             degree(second) >= maximum_degree)) {
            return false;
        }
        slots_.emplace(key, Slot{
            first, second, current_step, current_step, pinned});
        _bump_degree(first);
        _bump_degree(second);
        return true;
    }

    void observe(int first, int second, int status, Step current_step,
                 Step status2_memory_steps, bool pinned = false) {
        if (!reserve(first, second, current_step, 0, pinned)) return;
        Slot& slot = slots_.at(pair_key(first, second));
        slot.last_active_step = current_step;
        if (status >= 2) {
            slot.protected_until = std::max(
                slot.protected_until,
                saturating_add(current_step, status2_memory_steps));
        }
    }

    void expire_unseen(Step current_step, Step status1_memory_steps) {
        std::vector<std::uint64_t> expired;
        expired.reserve(slots_.size());
        for (const auto& item : slots_) {
            const Slot& slot = item.second;
            if (slot.pinned || slot.last_active_step == current_step) continue;
            const Step ordinary_until =
                saturating_add(slot.last_active_step, status1_memory_steps);
            if (current_step > std::max(ordinary_until,
                                        slot.protected_until)) {
                expired.push_back(item.first);
                if (slot.first >= 0 &&
                    slot.first < static_cast<int>(degree_.size()))
                    --degree_[static_cast<std::size_t>(slot.first)];
                if (slot.second >= 0 &&
                    slot.second < static_cast<int>(degree_.size()))
                    --degree_[static_cast<std::size_t>(slot.second)];
            }
        }
        for (const auto key : expired) slots_.erase(key);
    }

    std::size_t size() const { return slots_.size(); }

private:
    void _bump_degree(int idx) {
        if (idx < 0) return;
        if (idx >= static_cast<int>(degree_.size()))
            degree_.resize(static_cast<std::size_t>(idx) + 1, 0);
        ++degree_[static_cast<std::size_t>(idx)];
    }
    std::unordered_map<std::uint64_t, Slot> slots_;
    std::vector<int> degree_;
};

}  // namespace aa_neighbor_identity

#endif
