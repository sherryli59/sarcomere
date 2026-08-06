#ifndef COUNTER_RNG_H
#define COUNTER_RNG_H

#include <cstdint>

namespace counter_rng {

// Stable stochastic namespaces. Adding a draw to one object or subsystem must
// not shift any other object's random sequence.
enum class Stream : std::uint64_t {
    am_head = 0x243f6a8885a308d3ULL,
    aa_link_detachment = 0x13198a2e03707344ULL,
    aa_pair_association = 0xa4093822299f31d0ULL,
    legacy_aa_pair = 0x082efa98ec4e6c89ULL,
    explicit_titin_binding = 0x452821e638d01377ULL,
    explicit_titin_release = 0xbe5466cf34e90c6cULL,
};

inline std::uint64_t mix64(std::uint64_t value) {
    value += 0x9e3779b97f4a7c15ULL;
    value = (value ^ (value >> 30)) * 0xbf58476d1ce4e5b9ULL;
    value = (value ^ (value >> 27)) * 0x94d049bb133111ebULL;
    return value ^ (value >> 31);
}

inline std::uint64_t combine(std::uint64_t state, std::uint64_t value) {
    return mix64(state ^ mix64(value));
}

// A small counter-based stream keyed by simulation seed, timestep, stochastic
// namespace, and stable object identifiers. It has no mutable global or
// thread-owned state: scheduling changes cannot alter which values an object
// receives. The local draw counter permits variable-length rejection sampling
// without perturbing other objects.
class CounterRng {
public:
    CounterRng(
        std::uint64_t seed,
        std::uint64_t step,
        Stream stream,
        std::uint64_t object_first,
        std::uint64_t object_second = 0) {
        base_ = mix64(seed);
        base_ = combine(base_, step);
        base_ = combine(base_, static_cast<std::uint64_t>(stream));
        base_ = combine(base_, object_first);
        base_ = combine(base_, object_second);
    }

    std::uint64_t next_u64() {
        return mix64(base_ + counter_++ * 0x9e3779b97f4a7c15ULL);
    }

    double uniform() {
        return static_cast<double>(next_u64() >> 11) * 0x1.0p-53;
    }

    double uniform_pos() {
        return (static_cast<double>(next_u64() >> 11) + 0.5) * 0x1.0p-53;
    }

private:
    std::uint64_t base_ = 0;
    std::uint64_t counter_ = 0;
};

inline std::uint64_t canonical_pair_key(int first, int second) {
    const auto low = static_cast<std::uint32_t>(
        first < second ? first : second);
    const auto high = static_cast<std::uint32_t>(
        first < second ? second : first);
    return (static_cast<std::uint64_t>(low) << 32) |
           static_cast<std::uint64_t>(high);
}

}  // namespace counter_rng

#endif
