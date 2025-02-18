#include "SPERR3D_OMP_C.h"
#include "Sample.h"
#include "qoi/QoIInfo.hpp"
#include <algorithm>  // std::all_of()
#include <cassert>
#include <cstring>
#include <cmath>
#include <numeric>  // std::accumulate()
#include <iostream>
#ifdef USE_OMP
#include <omp.h>
#endif

void sperr::SPERR3D_OMP_C::set_num_threads(size_t n)
{
#ifdef USE_OMP
  if (n == 0)
    m_num_threads = omp_get_max_threads();
  else
    m_num_threads = n;
#endif
}

void sperr::SPERR3D_OMP_C::set_dims_and_chunks(dims_type vol_dims, dims_type chunk_dims)
{
  m_dims = vol_dims;

  // The preferred chunk size has to be between 1 and m_dims.
  for (size_t i = 0; i < m_chunk_dims.size(); i++)
    m_chunk_dims[i] = std::min(std::max(size_t{1}, chunk_dims[i]), vol_dims[i]);
}

void sperr::SPERR3D_OMP_C::set_psnr(double psnr)
{
  assert(psnr > 0.0);
  m_mode = CompMode::PSNR;
  m_quality = psnr;
}

void sperr::SPERR3D_OMP_C::set_tolerance(double pwe)
{
  assert(pwe > 0.0);
  m_mode = CompMode::PWE;
  m_quality = pwe;
}

void sperr::SPERR3D_OMP_C::set_bitrate(double bpp)
{
  assert(bpp > 0.0);
  m_mode = CompMode::Rate;
  m_quality = bpp;
}

#ifdef EXPERIMENTING
void sperr::SPERR3D_OMP_C::set_direct_q(double q)
{
  assert(q > 0.0);
  m_mode = CompMode::DirectQ;
  m_quality = q;
}
#endif

void sperr::SPERR3D_OMP_C::set_qoi_meta(const QoZ::QoIMeta &q_m)
{
  qoi_meta = q_m;
}

void sperr::SPERR3D_OMP_C::set_qoi_id(int q_id)
{
  qoi_meta.qoi_id = q_id;
}

void sperr::SPERR3D_OMP_C::set_qoi_string(std::string q_string)
{
  qoi_meta.qoi_string = q_string;
}

void sperr::SPERR3D_OMP_C::set_qoi_base(double q_base)
{
  qoi_meta.qoi_base = q_base;
}

void sperr::SPERR3D_OMP_C::set_qoi_analytical(bool ana)
{
  qoi_meta.analytical = ana;
}

void sperr::SPERR3D_OMP_C::set_qoi_tol(double q_tol)
{
  qoi_tol = q_tol;

}

void sperr::SPERR3D_OMP_C::set_qoi_block_size(int q_bs)
{
  qoi_block_size = q_bs;
}

void sperr::SPERR3D_OMP_C::set_qoi_k(double q_k)
{
  qoi_k = q_k;
}

void sperr::SPERR3D_OMP_C::set_high_prec(bool hp)
{
  m_high_prec = hp;

}


template <typename T>
auto sperr::SPERR3D_OMP_C::compress(const T* buf, size_t buf_len) -> RTNType
{
  static_assert(std::is_floating_point<T>::value, "!! Only floating point values are supported !!");
  if constexpr (std::is_same<T, float>::value)
    m_orig_is_float = true;
  else
    m_orig_is_float = false;
  //auto qoi = QoZ::GetQOI<double>(1, 1, 1, "x^2" );
  if (m_mode == sperr::CompMode::Unknown)
    return RTNType::CompModeUnknown;
  if (buf_len != m_dims[0] * m_dims[1] * m_dims[2])
    return RTNType::WrongLength;

  // First, calculate dimensions of individual chunk indices.
  if(qoi_meta.qoi_id>0 and qoi_tol>0 and qoi_block_size>1){
    for(auto &cd:m_chunk_dims){
      if(cd<qoi_block_size)
        cd = qoi_block_size;
      else
        cd -= cd%qoi_block_size;
    }
  }
  const auto chunk_idx = sperr::chunk_volume(m_dims, m_chunk_dims);
  const auto num_chunks = chunk_idx.size();


  // Let's prepare some data structures for compression!
  auto chunk_rtn = std::vector<RTNType>(num_chunks, RTNType::Good);
  m_encoded_streams.resize(num_chunks);

#ifdef USE_OMP
  m_compressors.resize(m_num_threads);
  for (auto& p : m_compressors) {
    if (p == nullptr)
      p = std::make_unique<SPECK3D_FLT>();
  }
#else
  if (m_compressor == nullptr)
    m_compressor = std::make_unique<SPECK3D_FLT>();
#endif
  auto bs_qoi_tol = qoi_tol;

  if(qoi_meta.qoi_id>0 and qoi_tol>0)
  {
    if(qoi_block_size > 1){//regional 
        //adjust qoieb
        double rate = 1.0;
      

       
        //conf.regionalQoIeb=conf.qoiEB;//store original regional eb
        double num_blocks = 1;
        double num_elements = 1;
        for(int i=0; i<m_dims.size(); i++){
            num_elements *= qoi_block_size;
            num_blocks *= (m_dims[i] - 1) / qoi_block_size + 1;
        }

        double q = 0.999;
        rate = estimate_rate_Hoeffdin(num_elements,1,q,qoi_k);
        //std::cout<<num_elements<<" "<<num_blocks<<" "<<conf.error_std_rate<<" "<<rate<<std::endl;
        
        rate = std::max(1.0,rate);//only effective for average. general: 1.0/sumai
        
        //std::cout<<"Pointwise QoI eb rate: " << rate << std::endl;
        

        if((qoi_meta.qoi_id == 1 and qoi_meta.qoi_string == "x") or qoi_meta.qoi_id == 11 ){
            if(qoi_block_size <=4){//it is a random number. to Fix
                rate = std::min(4.0,rate);
            }
            else if(qoi_block_size <=8){//it is a random number. to Fix
                rate = std::min(8.0,rate);
            }
            else if(qoi_block_size <=16){//it is a random number. to Fix
                rate = std::min(12.0,rate);
            }
            else if(qoi_block_size <=32){//it is a random number. to Fix
                rate = std::min(16.0,rate);
            }
            else {
              rate = std::min(32.0,rate);
            }
            //else if 
        }
        qoi_tol *= rate;

        //qoi->set_qoi_tolerance(qoi_tol);
    }

    
  }

#pragma omp parallel for num_threads(m_num_threads)
  for (size_t i = 0; i < num_chunks; i++) {
#ifdef USE_OMP
    auto& compressor = m_compressors[omp_get_thread_num()];
#else
    auto& compressor = m_compressor;
#endif

    // Gather data for this chunk, Setup compressor parameters, and compress!
    auto chunk = m_gather_chunk<T>(buf, m_dims, chunk_idx[i]);
    assert(!chunk.empty());

    if(qoi_meta.qoi_id>0 and qoi_tol>0){//qoi tuning
      std::cout<<"Tuning eb with qoi"<<std::endl;
      auto pwe = m_mode == CompMode::PWE ? m_quality : std::numeric_limits<double>::max();
      m_mode == CompMode::PWE;


      if((qoi_meta.qoi_id == 1 and qoi_meta.qoi_string == "x") or qoi_meta.qoi_id == 11 ){
          m_quality = std::min(qoi_tol,pwe);
          if(qoi_meta.qoi_id == 11){
            m_quality /= qoi_meta.lin_A;
          }
          if(qoi_block_size>1){
            auto qoi = QoZ::GetQOI<double>(qoi_meta, qoi_tol, m_quality);
            compressor->set_qoi(qoi);
            compressor->set_qoi_tol(bs_qoi_tol);
            compressor->set_qoi_block_size(qoi_block_size);
          }
      } 
      else{
      
        std::array<size_t,3> chunk_dims = {chunk_idx[i][1], chunk_idx[i][3], chunk_idx[i][5]};
        size_t chunk_ele_num = chunk_idx[i][1]*chunk_idx[i][3]*chunk_idx[i][5];


        







        if(qoi_block_size > 1){//regional 
          //adjust qoieb
          compressor->set_qoi_tol(bs_qoi_tol);
          compressor->set_qoi_block_size(qoi_block_size);
        }
        auto qoi = QoZ::GetQOI<double>(qoi_meta, qoi_tol, pwe);

        std::vector<double> ebs (chunk_ele_num);
      // use quantile to determine abs bound
    


          

          for (size_t i = 0; i < chunk_ele_num; i++){
              ebs[i] = qoi->interpret_eb(chunk[i]);
          }
          
          //double max_quantile_rate = 0.2;
          double quantile_rate = 0.2  ;//conf.quantile;//quantile
          //std::cout<<quantile<<std::endl;
          size_t k = std::ceil(quantile_rate * chunk_ele_num);
          k = std::max((size_t)1, std::min(chunk_ele_num, k)); 


          double best_abs_eb;

        
          std::vector<size_t> quantiles;
        
         for(auto i:{1.0,0.5,0.25,0.10,0.05,0.025,0.01})
             quantiles.push_back((size_t)(i*k));
         int quantile_num = quantiles.size();

              

              
              //std::sort(ebs.begin(),ebs.begin()+k+1);

    
          size_t best_quantile = 0;


          std::nth_element(ebs.begin(),ebs.begin()+quantiles[0], ebs.end());

          size_t last_quantile = quantiles[0]+1;

          double best_br = 9999;
          best_abs_eb = pwe;
                          
          int idx = 0;

          double sample_rate = 0.01;
          /*
          double length_sample_rate = pow(sample_rate,1.0/3.0);
          std::array<size_t,3> sample_dims = {(size_t)(chunk_idx[i][1]*length_sample_rate), (size_t)(chunk_idx[i][3]*length_sample_rate), (size_t)(chunk_idx[i][5]*length_sample_rate)};
          size_t sample_num = sample_dims[0]*sample_dims[1]*sample_dims[2];
          auto sampled_data = m_sample_center(chunk,chunk_dims,sample_dims);
          */

          size_t block_size = 31;

          std::vector<std::vector<double>>sampled_blocks;
          std::vector<std::vector<size_t>>starts;


          size_t profStride=std::max((size_t)1,block_size/4);//todo: bugfix for others
          bool profiling = true;

          size_t totalblock_num=1;  

          double prof_abs_threshold = ebs[quantiles[0]];
          //double sample_ratio = 5e-3;
          for(int i=0;i<3;i++){                      
              totalblock_num*=(size_t)((chunk_dims[i]-1)/block_size);
          }
          std::vector<size_t> reversed_dims = {chunk_dims[2],chunk_dims[1],chunk_dims[0]};
          std::vector<size_t> sample_dims = {block_size+1,block_size+1,block_size+1};
          std::array<size_t,3> sample_dims_arr = {block_size+1,block_size+1,block_size+1};
         
          sperr::profiling_block_3d<double,3>(chunk.data(),reversed_dims,starts,block_size, prof_abs_threshold,profStride);
          
          


          size_t num_filtered_blocks=starts.size();
          if(num_filtered_blocks<=(int)(0.6*sample_rate*totalblock_num))//todo: bugfix for others 
              profiling=false;
          if(qoi_meta.qoi_id == 23 and m_dims[1] == 1200 and std::abs(qoi_tol-1e-3)<=1e-10){

              std::cout<<"check"<<std::endl;
              profiling=false;
          }

          sperr::sampleBlocks<double,3>(chunk.data(),reversed_dims,block_size,sampled_blocks,sample_rate,profiling,starts,false);//todo: test var_first = true

          size_t sample_num=0;
          for(auto &block:sampled_blocks)
            sample_num += block.size();

          //done until here, next todo
          for(auto quantile:quantiles)
          {   
              if(idx!=0)
                  std::nth_element(ebs.begin(),ebs.begin()+quantile, ebs.begin()+last_quantile);

              
              auto cur_abs_eb = ebs[quantile];
              if(cur_abs_eb<=1e-15){
                if (idx == 0)
                  best_abs_eb = 1e-15;
                break;
              }
              qoi->set_global_eb(cur_abs_eb);
              // reset variables for average of square
              auto test_compressor = std::make_unique<SPECK3D_FLT>();
              test_compressor->set_dims(sample_dims_arr);
              test_compressor->set_tolerance(cur_abs_eb);
              test_compressor->set_qoi(qoi);
              double cur_br = 0;
              if(qoi_block_size > 1)
                test_compressor->set_qoi_tol(bs_qoi_tol);

              for (auto sampled_block:sampled_blocks){
                test_compressor->take_data(std::move(sampled_block));
                
                
                vec8_type test_encoded_stream;

                auto rtn = test_compressor->compress(m_high_prec);
                if(rtn!= RTNType::Good)
                  std::cout<<"Error"<<std::endl;

                test_encoded_stream.clear();
                test_encoded_stream.reserve(128);
                //m_encoded_streams[i].reserve(1280000);
                test_compressor->append_encoded_bitstream(test_encoded_stream);
                
                cur_br += test_encoded_stream.size()*8.0/(double)(sample_num);   
              }    
              std::cout << "current_eb = " << cur_abs_eb << ", current_br = " << cur_br << std::endl;
              if(cur_br < best_br * 1.02){//todo: optimize
                  best_br = cur_br;
                  best_abs_eb = cur_abs_eb;
                  best_quantile = quantile;
              }
              /*else if(cur_br>1.1*best_br and testConf.early_termination){
                  break;
              }*/

              last_quantile = quantile+1;
              idx++;

              //if(!test_compressor->has_lossless())
              //  break;
              
          }
          std::cout<<"Selected quantile: "<<(double)best_quantile/(double)chunk_ele_num<<std::endl;
          std::cout << "Best abs eb:  " << best_abs_eb << std::endl; 
          qoi->set_global_eb(best_abs_eb); 

          compressor->set_qoi(qoi);

          m_quality = best_abs_eb;
        }


    }







      //compressor->set_qoi(qoi);
      //compressor->set_qoi_tol(qoi_tol);
    


    compressor->take_data(std::move(chunk));
    compressor->set_dims({chunk_idx[i][1], chunk_idx[i][3], chunk_idx[i][5]});
    //std::cout<<qoi_tol<<" "<<qoi_id<<std::endl;
    /*
    if(qoi_id>0 and qoi_tol>0){
      auto pwe = m_mode == CompMode::PWE ? m_quality : std::numeric_limits<double>::max();
      auto qoi = QoZ::GetQOI<double>(qoi_id, qoi_tol, pwe, qoi_string );
      compressor->set_qoi(qoi);
      //compressor->set_qoi_tol(qoi_tol);
    }*/
    //std::cout<<chunk_idx[i][1]<<" "<<chunk_idx[i][3]<<" "<<chunk_idx[i][5]<<std::endl;//its reversed (fastest first)
    switch (m_mode) {
      case CompMode::PSNR:
        compressor->set_psnr(m_quality);
        break;
      case CompMode::PWE:
        compressor->set_tolerance(m_quality);
        break;
      case CompMode::Rate:
        compressor->set_bitrate(m_quality);
        break;
#ifdef EXPERIMENTING
      case CompMode::DirectQ:
        compressor->set_direct_q(m_quality);
        break;
#endif
      default:;  // So the compiler doesn't complain about missing cases.
    }
    chunk_rtn[i] = compressor->compress(m_high_prec);

    // Save bitstream for each chunk in `m_encoded_stream`.
    m_encoded_streams[i].clear();
    m_encoded_streams[i].reserve(128);
    //m_encoded_streams[i].reserve(1280000);
    compressor->append_encoded_bitstream(m_encoded_streams[i]);
  }

  auto fail = std::find_if_not(chunk_rtn.begin(), chunk_rtn.end(),
                               [](auto r) { return r == RTNType::Good; });
  if (fail != chunk_rtn.end())
    return (*fail);

  assert(std::none_of(m_encoded_streams.cbegin(), m_encoded_streams.cend(),
                      [](auto& s) { return s.empty(); }));

  return RTNType::Good;
}
template auto sperr::SPERR3D_OMP_C::compress(const float*, size_t) -> RTNType;
template auto sperr::SPERR3D_OMP_C::compress(const double*, size_t) -> RTNType;

auto sperr::SPERR3D_OMP_C::get_encoded_bitstream() const -> vec8_type
{
  auto header = m_generate_header();
  assert(!header.empty());
  auto header_size = header.size();
  auto stream_size = std::accumulate(m_encoded_streams.cbegin(), m_encoded_streams.cend(), 0lu,
                                     [](size_t a, const auto& b) { return a + b.size(); });
  header.resize(header_size + stream_size);

  auto itr = header.begin() + header_size;
  for (const auto& s : m_encoded_streams) {
    std::copy(s.cbegin(), s.cend(), itr);
    itr += s.size();
  }

  return header;
}

auto sperr::SPERR3D_OMP_C::m_generate_header() const -> sperr::vec8_type
{
  auto header = sperr::vec8_type();

  // The header would contain the following information
  //  -- a version number                     (1 byte)
  //  -- 8 booleans                           (1 byte)
  //  -- volume dimensions                    (4 x 3 = 12 bytes)
  //  -- (optional) chunk dimensions          (2 x 3 = 6 bytes)
  //  -- length of bitstream for each chunk   (4 x num_chunks)
  //
  auto chunk_idx = sperr::chunk_volume(m_dims, m_chunk_dims);
  const auto num_chunks = chunk_idx.size();
  assert(num_chunks != 0);
  if (num_chunks != m_encoded_streams.size())
    return header;
  auto header_size = size_t{0};
  if (num_chunks > 1)
    header_size = m_header_magic_nchunks + num_chunks * 4;
  else
    header_size = m_header_magic_1chunk + num_chunks * 4;

  header.resize(header_size);

  // Version number
  header[0] = static_cast<uint8_t>(SPERR_VERSION_MAJOR);
  size_t pos = 1;

  // 8 booleans:
  // bool[0]  : if this bitstream is a portion of another complete bitstream (progressive access).
  // bool[1]  : if this bitstream is for 3D (true) or 2D (false) data.
  // bool[2]  : if the original data is float (true) or double (false).
  // bool[3]  : if there are multiple chunks (true) or a single chunk (false).
  // bool[4-7]: unused
  //
  const auto b8 = std::array<bool, 8>{false,  // not a portion
                                      true,   // 3D
                                      m_orig_is_float,
                                      (num_chunks > 1),
                                      false,   // unused
                                      false,   // unused
                                      false,   // unused
                                      false};  // unused

  header[pos++] = sperr::pack_8_booleans(b8);

  // Volume dimensions
  const auto vdim = std::array{static_cast<uint32_t>(m_dims[0]), static_cast<uint32_t>(m_dims[1]),
                               static_cast<uint32_t>(m_dims[2])};
  std::memcpy(&header[pos], vdim.data(), sizeof(vdim));
  pos += sizeof(vdim);

  // Chunk dimensions, if there are more than one chunk.
  if (num_chunks > 1) {
    auto vcdim =
        std::array{static_cast<uint16_t>(m_chunk_dims[0]), static_cast<uint16_t>(m_chunk_dims[1]),
                   static_cast<uint16_t>(m_chunk_dims[2])};
    std::memcpy(&header[pos], vcdim.data(), sizeof(vcdim));
    pos += sizeof(vcdim);
  }

  // Length of bitstream for each chunk.
  for (const auto& stream : m_encoded_streams) {
    assert(stream.size() <= uint64_t{std::numeric_limits<uint32_t>::max()});
    uint32_t len = stream.size();
    std::memcpy(&header[pos], &len, sizeof(len));
    pos += sizeof(len);
  }
  assert(pos == header_size);

  return header;
}

template <typename T>
auto sperr::SPERR3D_OMP_C::m_gather_chunk(const T* vol,
                                          dims_type vol_dim,
                                          std::array<size_t, 6> chunk) -> vecd_type
{
  auto chunk_buf = vecd_type();
  if (chunk[0] + chunk[1] > vol_dim[0] || chunk[2] + chunk[3] > vol_dim[1] ||
      chunk[4] + chunk[5] > vol_dim[2])
    return chunk_buf;

  chunk_buf.resize(chunk[1] * chunk[3] * chunk[5]);
  const auto row_len = chunk[1];

  size_t idx = 0;
  for (size_t z = chunk[4]; z < chunk[4] + chunk[5]; z++) {
    const size_t plane_offset = z * vol_dim[0] * vol_dim[1];
    for (size_t y = chunk[2]; y < chunk[2] + chunk[3]; y++) {
      const auto start_i = plane_offset + y * vol_dim[0] + chunk[0];
      std::copy(vol + start_i, vol + start_i + row_len, chunk_buf.begin() + idx);
      idx += row_len;
    }
  }

  // Will be subject to Named Return Value Optimization.
  return chunk_buf;
}
template auto sperr::SPERR3D_OMP_C::m_gather_chunk(const float*,
                                                   dims_type,
                                                   std::array<size_t, 6>) -> vecd_type;
template auto sperr::SPERR3D_OMP_C::m_gather_chunk(const double*,
                                                   dims_type,
                                                   std::array<size_t, 6>) -> vecd_type;



auto sperr::SPERR3D_OMP_C::m_sample_center(vecd_type chunk,std::array<size_t, 3> chunk_dims,std::array<size_t,3>sample_dims) -> vecd_type
{
  std::array<size_t,3>starts= {(chunk_dims[0]-sample_dims[0])/2,(chunk_dims[1]-sample_dims[1])/2,(chunk_dims[2]-sample_dims[2])/2};
  //std::array<size_t,3>ends= {(chunk_dim[0]+sample_dim[0])/2,(chunk_dim[1]+sample_dim[1])/2,(chunk_dim[2]+sample_dim[2])/2};
  vecd_type sampled_data;
  size_t idx = 0;
  size_t y_offset = chunk_dims[0];
  size_t z_offset = y_offset * chunk_dims[1];
  for (size_t z = starts[2]; z < starts[2]+sample_dims[2]; z++) {
    const size_t plane_offset = z * z_offset;
    for (size_t y = starts[1]; y < starts[1]+sample_dims[1]; y++) {
      const auto start_idx = plane_offset + y * y_offset + starts[0];
      sampled_data.insert(sampled_data.end(), chunk.begin()+start_idx, chunk.begin()+start_idx+sample_dims[0]);
    }
  }
  return sampled_data;

}

double sperr::SPERR3D_OMP_C::estimate_rate_Hoeffdin(size_t n, size_t N, double q, double k = 3.0){//n: element_per_block N: num_blocks q: confidence
    //no var information
    //if gaussian, just multiply k 
   
    /*
    if (q>=0.95 and N >= 1000){
        return sqrt( -n / ( 2 * ( log(1-q) - log(2*N) ) ) );
    }
    else{
        return sqrt( -n / ( 2 * ( log(1- pow(q,1.0/N) ) - log(2) ) ) );
    }
    */

    double p;
    if (q>=0.95 and N >= 1000){
        p = (1-q)/N;
    }
    else{
        p = 1- pow(q,1.0/N);
    }


    return k*sqrt(0.5*n/log(2.0/p));
    
}
