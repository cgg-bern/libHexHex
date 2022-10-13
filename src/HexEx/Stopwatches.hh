#pragma once

#include <libTimekeeper/StopWatch.hh>

namespace HexEx {
    using Timekeeper::ScopedStopWatch;
}

namespace HexEx::sw {
    using HSW = Timekeeper::HierarchicalStopWatch;

    inline HSW root("hexex");
    inline HSW read_file("read file", root);
    inline HSW write_file("write file", root);
    inline HSW extract("extract", root);
    inline HSW extractHVertices("extractHVertices()", extract);
}
