#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>

#define new new_
#include "deps/bsdiff/bsdiff.c"
#include "deps/bsdiff/bspatch.c"
#undef new

#include <string>
#include <utility>

#include "collage.hpp"

namespace collage {

    namespace {
        /* variable length encoding */
        std::string vlebit( size_t i ) {
            std::string out;
            do {
                out += (unsigned char)( 0x80 | (i & 0x7f));
                i >>= 7;
            } while( i > 0 );
            *out.rbegin() ^= 0x80;
            return out;
        }
        size_t vlebit( const char *&i ) {
            size_t out = 0, j = -7;
            do {
                out |= ((size_t(*i) & 0x7f) << (j += 7) );
            } while( size_t(*i++) & 0x80 );
            return out;
        }

        /* wrappers */
        static int bs_write(struct bsdiff_stream* stream, const void* buffer, int size) {
            return *((std::string*)stream->opaque) += std::string((char *)buffer,size), 0;
        }

        static int bs_read(const struct bspatch_stream * stream, void* buffer, int size) {
            std::pair<const char *,const char *> &opaque = *( (std::pair<const char *,const char *> *)stream->opaque );
            if( opaque.first + size <= opaque.second ) {
                memcpy( buffer, opaque.first, size );
                opaque.first += size;
                return 0;
            } else {
                return -1;
            }
        }
    }

    bool diff( std::string &result, const char *from0, const char *from1, const char *to0, const char *to1, unsigned Q ) {
        result = std::string();
        switch( Q ) {
            default:
            case BSDIFF: {
                struct bsdiff_stream diff_stream;
                diff_stream.opaque = (void *)&result;
                diff_stream.malloc = malloc;
                diff_stream.free = free;
                diff_stream.write = bs_write;
                if( 0 == bsdiff( (const uint8_t *)from0, from1 - from0, (const uint8_t *)to0, to1 - to0, &diff_stream ) ) {
                    return true;
                }
            }
        }
        return false;
    }

    /* result must be resized in advance */
    bool patch( std::string &result, const char *from0, const char *from1, const char *diff0, const char *diff1, unsigned Q ) {
        switch( Q ) {
            default:
            case BSDIFF: {
                std::pair<const char *,const char *> pair( diff0, diff1 );
                struct bspatch_stream patch_stream;
                patch_stream.opaque = (void *)&pair;
                patch_stream.read = bs_read;
                if( 0 == bspatch( (const uint8_t *)from0, from1 - from0, (uint8_t *)(&result[0]), result.size(), &patch_stream ) ) {
                    return true;
                }                
            }           
        }
        return false;
    }

    std::string diff( const std::string &from, const std::string &to, unsigned Q ) {
        std::string result;
        if( diff( result, from.c_str(), from.c_str() + from.size(), to.c_str(), to.c_str() + to.size(), Q ) ) {
            return std::string() + char(Q & 0x7F) + vlebit(to.size()) + result;
        } else {
            return std::string();
        }
    }

    std::string patch( const std::string &from, const std::string &diff ) {
        const char *diff8 = diff.c_str();
        unsigned Q = *diff8++;
        std::string result( vlebit(diff8), '\0' );
        if( patch( result, from.c_str(), from.c_str() + from.size(), diff8, diff8 + (diff.size() - (diff8 - diff.c_str()) ), Q ) ) {
            return result;
        } else {
            return std::string();
        }
    }
}
