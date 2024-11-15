#ifndef ZSTD_LOSSLESS
#define ZSTD_LOSSLESS

//
// This is the base class of 1D, 2D, and 3D integer SPECK implementations.
//


#include "Bitmask.h"
#include "Bitstream.h"
#include "zstd.h"
#include <vector>
namespace zstd {

class Lossless_zstd {

 public:
  // Constructor and destructor
  Lossless_zstd(){}
  
  template<class T>
  void encode(const std::vector<T> & data){

    size_t dataLength = data.size()*sizeof(T);

    size_t estimatedCompressedSize = (dataLength < 100 ? 200 : size_t(dataLength * 1.2));
    uint8_t *compressBytes = new uint8_t[estimatedCompressedSize];
    uint8_t *compressBytesPos = compressBytes;
    memcpy(compressBytesPos,&dataLength,sizeof(dataLength));
    compressBytesPos+=sizeof(dataLength);

    byte_size = ZSTD_compress(compressBytesPos, estimatedCompressedSize, data.data(), dataLength,
                              compression_level);
    byte_size+= sizeof(dataLength);
    m_bit_buffer.parse_bitstream(compressBytes, byte_size);
    delete []compressBytes;


  }
  void decode();

  // Input
  void use_bitstream(const void* p, size_t len){

    const auto* const p = static_cast<const uint8_t*>(p);
    std::memcpy(&byte_size, p, sizeof(byte_size));
    m_bit_buffer.parse_bitstream(p +sizeof(byte_size), byte_size);

  }
  size_t encoded_bitstream_len(){
    return byte_size+sizeof(byte_size);
  }
  // Output
  //auto encoded_bitstream_len() const -> size_t;
  void append_encoded_bitstream(vec8_type& buf) const{
    const auto app_size = this->encoded_bitstream_len();
    const auto orig_size = buffer.size();
    buffer.resize(orig_size + app_size);
    auto* const ptr = buffer.data() + orig_size;

    // Step 2: fill header
    size_t pos = 0;
    std::memcpy(ptr + pos, &byte_size, sizeof(byte_size));
    pos += sizeof(byte_size);

    // Step 3: assemble the right amount of bits into bytes.
    // See discussion on the number of bits to pack in function `encoded_bitstream_len()`.
    m_bit_buffer.write_bitstream(ptr + pos, byte_size);
  }


 protected:

  size_t byte_size = 0;
  Bitstream m_bit_buffer;
  compression_level = 3;
};

};  // namespace zstd

#endif
