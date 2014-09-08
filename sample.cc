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
