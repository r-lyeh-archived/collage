Collage <a href="https://travis-ci.org/r-lyeh/collage"><img src="https://api.travis-ci.org/r-lyeh/collage.svg?branch=master" align="right" /></a>
=======

- Collage is a lightweight C++ library to diff and patch data.
- Collage provides interface to bsdiff/bspatch libraries (for now).
- Collage is tiny. Single header and source files.
- Collage is stand-alone. All dependencies are embedded.
- Collage is cross-platform.
- Collage is zlib/libpng licensed.

## sample
```c++
#include <iostream>
#include <cassert>
#include "collage.hpp"

int main() {
    std::string source = "hello world and thanks";
    std::string target = "hello cruel \x1 world. thanks for the fish.";

    std::string patch = collage::diff( source, target );
    assert( !patch.empty() );

    std::string patched = collage::patch( source, patch );
    assert( !patched.empty() );
    assert( target == patched );

    std::cout << "'" << source << "' + " << patch.size() << "-bytes patch == '" << patched << "'" << std::endl;
    std::cout << "All ok." << std::endl;
}
```

## possible output
```
'hello world and thanks' + 46-bytes patch == 'hello cruel ☺ world. thanks for the fish.'
All ok.
```

## licenses
- [collage](https://github.com/r-lyeh/collage), zlib/libpng licensed.
- [bsdiff](https://github.com/mendsley/bsdiff), by Colin Percival and Matthew Endsley, BSD2 licensed.
- [sais-lite](https://github.com/davehughes/sais), by Yuta Mori, MIT licensed.
