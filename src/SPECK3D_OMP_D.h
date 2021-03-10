//
// This is a class that performs SPECK3D decompression, and also utilizes OpenMP
// to achieve parallelization: input to this class is supposed to be smaller chunks
// of a bigger volume, and each chunk is decompressed individually before returning
// back the big volume.
//

#ifndef SPECK3D_OMP_D_H
#define SPECK3D_OMP_D_H

#include "SPECK3D_Decompressor.h"

using speck::RTNType;

class SPECK3D_OMP_D
{
public:

    // Parse and separate data segments for each chunk.
    // Individual decompressors will then hold segments of the bitstream.
    auto use_bitstream( const void*, size_t ) -> RTNType;

    auto set_bpp(float) -> RTNType;
    void set_num_threads( size_t );

    auto decompress() -> RTNType;

    // Give up the copy that this class holds
    auto release_data_volume() -> speck::smart_buffer_d;

    // Make a copy and then return the copy as a speck::smart_pointer
    template<typename T>
    auto get_data_volume() const -> std::pair<std::unique_ptr<T[]>, size_t>;


private:

    size_t      m_dim_x       = 0;  // Dimension of the entire volume
    size_t      m_dim_y       = 0;  // Dimension of the entire volume
    size_t      m_dim_z       = 0;  // Dimension of the entire volume
    size_t      m_chunk_x     = 0;  // Preferred dimension for a chunk.
    size_t      m_chunk_y     = 0;  // Preferred dimension for a chunk.
    size_t      m_chunk_z     = 0;  // Preferred dimension for a chunk.
    size_t      m_num_threads = 1; // number of theads to use in OpenMP sections

    float       m_bpp         = 0.0;
    
    const size_t m_header_magic = 26; // header size would be this number + num_chunks * 4

    std::vector<SPECK3D_Decompressor>   m_decompressors;
    speck::smart_buffer_d               m_vol_buf;

    //
    // Private methods
    //
    auto m_parse_header( const void*, size_t ) -> std::vector<size_t>;

};

#endif








