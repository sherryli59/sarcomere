#include "counter_rng.h"

#include <cassert>
#include <cmath>
#include <cstdint>

int main() {
    using counter_rng::CounterRng;
    using counter_rng::Stream;

    CounterRng first(17, 29, Stream::am_head, 3, 11);
    CounterRng repeated(17, 29, Stream::am_head, 3, 11);
    CounterRng other_object(17, 29, Stream::am_head, 3, 12);
    CounterRng other_stream(17, 29, Stream::aa_link_detachment, 3, 11);
    bool object_differs = false;
    bool stream_differs = false;
    for (int draw = 0; draw < 32; ++draw) {
        const auto value = first.next_u64();
        assert(value == repeated.next_u64());
        object_differs = object_differs || value != other_object.next_u64();
        stream_differs = stream_differs || value != other_stream.next_u64();
    }
    assert(object_differs);
    assert(stream_differs);

    // Basic deterministic distribution smoke test. This is not a replacement
    // for TestU01, but catches broken conversion, short cycles, and severe
    // correlation before the generator reaches the simulator.
    constexpr int samples = 200000;
    CounterRng stats(1234, 0, Stream::aa_pair_association, 99);
    double sum = 0.0;
    double sum_squared = 0.0;
    for (int i = 0; i < samples; ++i) {
        const double value = stats.uniform();
        assert(value >= 0.0 && value < 1.0);
        sum += value;
        sum_squared += value * value;
    }
    const double mean = sum / samples;
    const double variance = sum_squared / samples - mean * mean;
    assert(std::abs(mean - 0.5) < 0.003);
    assert(std::abs(variance - 1.0 / 12.0) < 0.0015);

    CounterRng positive(1234, 0, Stream::aa_pair_association, 99);
    for (int i = 0; i < samples; ++i) {
        const double value = positive.uniform_pos();
        assert(value > 0.0 && value < 1.0);
    }

    assert(
        counter_rng::canonical_pair_key(7, 19) ==
        counter_rng::canonical_pair_key(19, 7));
    assert(
        counter_rng::canonical_pair_key(7, 19) !=
        counter_rng::canonical_pair_key(7, 20));
    return 0;
}
