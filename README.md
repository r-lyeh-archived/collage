Collage
=======

- Collage is a lightweight C++ library to diff and patch data.
- Collage provides interface to bsdiff/bspatch libraries (for now).
- Collage is tiny. Single header and source files.
- Collage is stand-alone. All dependencies are embedded.
- Collage is cross-platform.
- Collage is MIT licensed.

### sample
```c++
#include <iostream>
#include <cassert>
#include "collage.hpp"

int main() {
    std::string source = "hello world";
    std::string target = "hello cruel \x1 world";

    std::string patch = collage::diff( source, target );
    assert( !patch.empty() );

    std::string patched = collage::patch( source, patch );
    assert( !patched.empty() );
    assert( target == patched );

    std::cout << "('" << source << "' << " << patch.size() << "-bytes patch) == '" << patched << "'" << std::endl;
    std::cout << "All ok." << std::endl;
}
```

### licenses
- [Collage](https://github.com/r-lyeh/collage), MIT licensed.
- [bsdiff](https://github.com/mendsley/bsdiff), by Colin Percival and Matthew Endsley, BSD2 licensed.
