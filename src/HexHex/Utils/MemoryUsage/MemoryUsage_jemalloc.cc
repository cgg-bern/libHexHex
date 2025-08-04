#include <HexHex/Utils/MemoryUsage/MemoryUsage.hh>
#include <jemalloc/jemalloc.h>
#include <stdexcept>
#include <iostream>

namespace HexHex::MemoryStats {

size_t peak_get() {
    size_t sz;
    static uint64_t epoch = 1;
    sz = sizeof(epoch);
    mallctl("epoch", &epoch, &sz, &epoch, sz);
    size_t peak = 0;
    if (0 != mallctl("thread.peak.read", &peak, &sz, NULL, 0)) {
        throw std::runtime_error("mallctl failure");
    }
    return peak;
}

void peak_reset() {
    if (0 != mallctl("thread.peak.reset", NULL, 0, NULL, 0)) {
        throw std::runtime_error("mallctl failure");
    }
}

}
