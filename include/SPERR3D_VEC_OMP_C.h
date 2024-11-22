//
// This is a class that performs SPERR3D compression, and also utilizes OpenMP
// to achieve parallelization: the input volume is divided into smaller chunks
// and then they're processed individually.
//

#ifndef SPERR3D_VEC_OMP_C_H
#define SPERR3D_VEC_OMP_C_H

#include "SPECK3D_VEC_FLT.h"

namespace sperr {

class SPERR3D_OMP_C {
 public:
  // If 0 is passed in, the maximal number of threads will be used.
  void set_num_threads(size_t);

  // Note on `chunk_dims`: it's a preferred value, but when the volume dimension is not
  //    divisible by chunk dimensions, the actual chunk dimension will change.
  void set_dims_and_chunks(dims_type vol_dims, dims_type chunk_dims);

  void set_psnr(double);
  void set_tolerance(double);
  void set_bitrate(double);
  void set_qoi_id(int);
  void set_qoi_string(std::string);
  void set_qoi_tol(double);
  void set_qoi_block_size(int);
  void set_qoi_k(double);
#ifdef EXPERIMENTING
  void set_direct_q(double);
#endif

  // Apply compression on a volume pointed to by `buf`.
  template <typename T>
  auto compress(const T* buf1,const T* buf2,const T* buf3, size_t buf_len) -> RTNType;

  // Output: produce a vector containing the encoded bitstream.
  auto get_encoded_bitstream() const -> std::array<vec8_type,3>;

 private:
  bool m_orig_is_float = true;  // The original input precision is saved in header.
  CompMode m_mode = CompMode::Unknown;
  double m_quality = 0.0;
  std::array<double,3> m_quality_v = {0.0,0.0,0.0};
  dims_type m_dims = {0, 0, 0};        // Dimension of the entire volume
  dims_type m_chunk_dims = {0, 0, 0};  // Preferred dimensions for a chunk
  std::array<std::vector<vec8_type>,3> m_encoded_streams;
  int qoi_id = 0;
  std::string qoi_string = "x^2 + y^2 + z^2";
  double qoi_tol = 0.0;
  //double block_qoi_tol = 0.0;
  int qoi_block_size = 1;
  double qoi_k = 2;


#ifdef USE_OMP
  size_t m_num_threads = 1;

  // It turns out that the object of `SPECK3D_FLT` is not copy-constructible, so it's
  //    a little difficult to work with a container (std::vector<>), so we ask the
  //    container to store pointers (which are trivially constructible) instead.
  //
  std::vector<std::unique_ptr<SPECK3D_VEC_FLT>> m_compressors;
#else
  // This single instance of compressor doesn't need to be allocated on the heap;
  // rather, it's just to keep consistency with the USE_OMP case.
  std::unique_ptr<SPECK3D_VEC_FLT> m_compressor;
#endif

  // The eventual header size would be this magic number + num_chunks * 4
  const size_t m_header_magic_nchunks = 20;
  const size_t m_header_magic_1chunk = 14;

  //
  // Private methods
  //
  auto m_generate_header() const -> vec8_type;

  // Gather a chunk from a bigger volume.
  // If the requested chunk lives outside of the volume, whole or part,
  //    this function returns an empty vector.
  template <typename T>
  auto m_gather_chunk(const T* vol, dims_type vol_dim, std::array<size_t, 6> chunk) -> vecd_type;

  auto m_sample_center(std::array<vecd_type,3> chunk,std::array<size_t, 3> chunk_dim, std::array<size_t,3>sample_dim) -> std::array<vecd_type,3>;

  double estimate_rate_Hoeffdin(size_t n, size_t N, double q, double k);
};

}  // End of namespace sperr

#endif