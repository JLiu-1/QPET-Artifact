#include <algorithm> // std::all_of()
#include <cassert>
#include <cstring>   // std::memcpy()
#include <cmath>     // std::sqrt()
#include <numeric>   // std::accumulate()
#include <type_traits>   // std::is_same_v()
#include <vector>   // std::is_same_v()

#include "Conditioner.h"

speck::Conditioner::Conditioner(bool sub_mean)
{
    m_s_mean = sub_mean;
}

void speck::Conditioner::toggle_all_false()
{
    m_s_mean = false;
    m_d_rms  = false;
}

void speck::Conditioner::toggle_subtract_mean( bool b )
{
    m_s_mean = b;
}

void speck::Conditioner::toggle_divide_by_rms( bool b )
{
    m_d_rms = b;
}

auto speck::Conditioner::get_meta_size() const -> size_t
{
    return m_meta_size;
}


auto speck::Conditioner::m_calc_mean( const vecd_type& buf ) const -> double
{
    assert( buf.size() % m_num_strides == 0 );

    m_stride_buf.resize( m_num_strides );
    const size_t stride_size = buf.size() / m_num_strides;

    for (size_t s = 0; s < m_num_strides; s++) {
        auto begin = buf.begin() + stride_size * s;
        auto end   = begin + stride_size;
        m_stride_buf[s] = std::accumulate( begin, end, double{0.0} ) / double(stride_size);
    }

    double sum = std::accumulate( m_stride_buf.begin(), m_stride_buf.end(), double{0.0} );

    return (sum / double(m_stride_buf.size()));
}


auto speck::Conditioner::m_calc_rms( const vecd_type& buf ) const -> double
{
    assert( buf.size() % m_num_strides == 0 );

    m_stride_buf.resize( m_num_strides );
    const size_t stride_size = buf.size() / m_num_strides;

    for (size_t s = 0; s < m_num_strides; s++) {
        auto begin = buf.begin() + stride_size * s;
        auto end   = begin + stride_size;
        m_stride_buf[s]  = std::accumulate( begin, end, double{0.0}, 
                           [](auto a, auto b){return a + b * b;} );

        m_stride_buf[s] /= double(stride_size);
    }

    double rms = std::accumulate( m_stride_buf.begin(), m_stride_buf.end(), double{0.0} );
    rms /= double(m_stride_buf.size());
    rms  = std::sqrt(rms);

    return rms;
}


void speck::Conditioner::m_adjust_strides( size_t len ) const
{
    // Let's say 2048 is a starting point
    m_num_strides = 2048;
    if( len % m_num_strides == 0 )
        return;

    size_t num = 0;

    // First, try to increase till 2^14 = 16,384
    for( num = m_num_strides; num <= 16'384; num++ ) {
        if( len % num == 0 )
            break;
    }

    if( len % num == 0 ) {
        m_num_strides = num;
        return;
    }

    // Second, try to decrease till 1, at which point it must work.
    for( num = m_num_strides; num > 0; num-- ) {
        if( len % num == 0 )
            break;
    }

    m_num_strides = num;
}


auto speck::Conditioner::condition( vecd_type& buf ) const -> std::pair<RTNType, meta_type>
{
    meta_type meta;
    meta.fill( 0 );
    double mean = 0.0;
    double rms  = 1.0;

    // If divide_by_rms is requested, need to make sure rms is non-zero
    if( m_d_rms ) {
        if( std::all_of( buf.begin(), buf.end(), [](auto v){return v == 0.0;} ) )
            return {RTNType::Error, meta};
    }

    m_adjust_strides( buf.size() );

    // Perform subtract mean first
    if( m_s_mean ) {
        mean = m_calc_mean( buf );
        auto minus_mean = [mean](auto& v){v -= mean;};
        std::for_each( buf.begin(), buf.end(), minus_mean );
    }

    // Perform divide_by_rms second
    if( m_d_rms ) {
        rms = m_calc_rms( buf );
        auto div_rms = [rms](auto& v){v /= rms;};
        std::for_each( buf.begin(), buf.end(), div_rms );
    }

    // pack meta
    bool b8[8] = {m_s_mean, m_d_rms, false, false, false, false, false, false};
    speck::pack_8_booleans( meta[0], b8 );
    size_t pos = 1;
    std::memcpy( meta.data() + pos, &mean, sizeof(mean) );
    pos += sizeof(mean);
    std::memcpy( meta.data() + pos, &rms, sizeof(rms) );
    pos += sizeof(rms);
    assert( pos == m_meta_size );

    return {RTNType::Good, meta};
}


auto speck::Conditioner::inverse_condition( vecd_type& buf, const uint8_t* meta ) const -> RTNType
{
    // unpack meta
    bool   b8[8];
    speck::unpack_8_booleans( b8, meta[0] );
    double mean = 0.0;
    double rms  = 1.0;
    size_t pos  = 1;
    std::memcpy( &mean, meta + pos, sizeof(mean) );
    pos += sizeof(mean);
    std::memcpy( &rms, meta + pos, sizeof(rms) );
    pos += sizeof(rms);
    assert( pos == m_meta_size );

    m_adjust_strides( buf.size() );

    // Perform inverse of divide_by_rms, which is multiply by rms
    if( b8[1] ) {
        auto mul_rms = [rms](auto& v){v *= rms;};
        std::for_each( buf.begin(), buf.end(), mul_rms );
    }

    // Perform inverse of subtract_mean, which is add mean.
    if( b8[0] ) {
        auto plus_mean = [mean](auto& v){v += mean;};
        std::for_each( buf.begin(), buf.end(), plus_mean );
    }

    return RTNType::Good;
}




