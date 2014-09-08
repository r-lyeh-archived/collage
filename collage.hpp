#pragma once
#include <string>

namespace collage {
    // available libraries 
    enum { BSDIFF = 0 };

    // api
    std::string diff( const std::string &from, const std::string &to, unsigned library = BSDIFF );
    std::string patch( const std::string &from, const std::string &diff );
}
