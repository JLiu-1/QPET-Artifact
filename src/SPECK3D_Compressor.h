//
// This is a the class object that is supposed to be used by most users, because
// it provides easy-to-use APIs.
// Functionality wise, it does not bring anything new though.
// 

#ifndef SPECK3D_COMPRESSOR_H
#define SPECK3D_COMPRESSOR_H


#include "CDF97.h"
#include "SPECK3D.h"

using speck::RTNType;

class SPECK3D_Compressor {

public:
    // Constructor
    SPECK3D_Compressor( size_t x, size_t y, size_t z );

    // Accept incoming data: copy from a raw memory block
    template< typename T >
    auto copy_data( const T* p, size_t len ) -> RTNType;

    // Accept incoming data by taking ownership of the memory block
    auto take_data( speck::buffer_type_d buf, size_t len ) -> RTNType;

    // Accept incoming data by reading a file from disk.
    // The file is supposed to contain blocks of 32-bit floating values, with
    // X direction varying fastest, and Z direction varying slowest.
    auto read_floats( const char* filename ) -> RTNType;

#ifdef QZ_TERM
    void set_qz_level( int32_t );
#else
    auto set_bpp( float ) -> RTNType;
#endif

    auto compress() -> RTNType;

    // Provide a copy of the encoded bitstream to the caller.
    auto get_encoded_bitstream() const -> std::pair<speck::buffer_type_uint8, size_t>;
    auto write_bitstream( const char* filename ) const -> RTNType;

private:
    const size_t m_dim_x, m_dim_y, m_dim_z, m_total_vals;
    speck::buffer_type_d  m_val_buf;

    speck::CDF97    m_cdf;
    speck::SPECK3D  m_encoder;

#ifdef QZ_TERM
    int32_t m_qz_lev = 0;
#else
    float m_bpp = 0.0;
#endif
};


#endif