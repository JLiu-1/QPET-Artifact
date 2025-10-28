#include "SPERR3D_VEC_OMP_C.h"
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

void sperr::SPERR3D_VEC_OMP_C::set_num_threads(size_t n)
{
#ifdef USE_OMP
  if (n == 0)
    m_num_threads = omp_get_max_threads();
  else
    m_num_threads = n;
#endif
}

void sperr::SPERR3D_VEC_OMP_C::set_dims_and_chunks(dims_type vol_dims, dims_type chunk_dims)
{
  m_dims = vol_dims;

  // The preferred chunk size has to be between 1 and m_dims.
  for (size_t i = 0; i < m_chunk_dims.size(); i++)
    m_chunk_dims[i] = std::min(std::max(size_t{1}, chunk_dims[i]), vol_dims[i]);
}

void sperr::SPERR3D_VEC_OMP_C::set_psnr(std::vector<double> psnr)
{
  
  m_mode = CompMode::PSNR;
  for(auto i:{0,1,2})
    m_quality[i] = psnr[i];
}

void sperr::SPERR3D_VEC_OMP_C::set_tolerance(std::vector<double> pwe)
{
  assert(pwe > 0.0);
  m_mode = CompMode::PWE;
  for(auto i:{0,1,2})
    m_quality[i] = pwe[i];
}

void sperr::SPERR3D_VEC_OMP_C::set_bitrate(std::vector<double> bpp)
{
  assert(bpp > 0.0);
  m_mode = CompMode::Rate;
  for(auto i:{0,1,2})
    m_quality[i] = bpp[i];
}

#ifdef EXPERIMENTING
void sperr::SPERR3D_VEC_OMP_C::set_direct_q(std::vector<double> q)
{
  assert(q > 0.0);
  m_mode = CompMode::DirectQ;
  for(auto i:{0,1,2})
    m_quality[i] = q[i];
}
#endif

void sperr::SPERR3D_VEC_OMP_C::set_qoi_id(int q_id)
{
  qoi_id = q_id;
}

void sperr::SPERR3D_VEC_OMP_C::set_qoi_string(std::string q_string)
{
  qoi_string = q_string;
}

void sperr::SPERR3D_VEC_OMP_C::set_qoi_tol(double q_tol)
{
  qoi_tol = q_tol;

}

void sperr::SPERR3D_VEC_OMP_C::set_qoi_block_size(int q_bs)
{
  qoi_block_size = q_bs;
}

void sperr::SPERR3D_VEC_OMP_C::set_qoi_k(double q_k)
{
  qoi_k = q_k;
}


template <typename T>
auto sperr::SPERR3D_VEC_OMP_C::compress(const T* buf1, const T* buf2, const T* buf3, size_t buf_len) -> RTNType
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
  const auto chunk_idx = sperr::chunk_volume(m_dims, m_chunk_dims);
  const auto num_chunks = chunk_idx.size();

  // Let's prepare some data structures for compression!
  auto chunk_rtn = std::vector<RTNType>(num_chunks, RTNType::Good);
  for(auto i : {0,1,2}){
    m_encoded_streams[i].resize(num_chunks);
  

  #ifdef USE_OMP
    m_compressors[i].resize(m_num_threads);
    for (auto& p : m_compressors[i]) {
      if (p == nullptr)
        p = std::make_shared<SPECK3D_FLT>();
    }
  #else
    if (m_compressor[i] == nullptr)
      m_compressor[i] = std::make_shared<SPECK3D_FLT>();
  #endif
  }

#pragma omp parallel for num_threads(m_num_threads)
  for (size_t i = 0; i < num_chunks; i++) {
#ifdef USE_OMP
    std::array<std::shared_ptr<SPECK3D_FLT>,3> compressor = {m_compressors[0][omp_get_thread_num()],m_compressors[1][omp_get_thread_num()],m_compressors[2][omp_get_thread_num()]};
#else
    std::array<std::shared_ptr<SPECK3D_FLT>,3> compressor = m_compressor;
#endif

    // Gather data for this chunk, Setup compressor parameters, and compress!
    std::array<vecd_type,3> chunk = {m_gather_chunk<T>(buf1, m_dims, chunk_idx[i]),m_gather_chunk<T>(buf2, m_dims, chunk_idx[i]),m_gather_chunk<T>(buf3, m_dims, chunk_idx[i])};
    size_t chunk_ele_num = chunk_idx[i][1]*chunk_idx[i][3]*chunk_idx[i][5];
    assert(!chunk.empty() and !chunk[0].empty());

    if(qoi_id>0 and qoi_tol>0){//qoi tuning
      std::cout<<"Tuning eb with qoi"<<std::endl;
      auto pwe = m_mode == CompMode::PWE ? m_quality : std::array<double,3>{std::numeric_limits<double>::max(),std::numeric_limits<double>::max(),std::numeric_limits<double>::max()};
      m_mode == CompMode::PWE;
      
      std::array<size_t,3> chunk_dims = {chunk_idx[i][1], chunk_idx[i][3], chunk_idx[i][5]};
      
      double sample_rate = 0.01;
      double length_sample_rate = pow(sample_rate,1.0/3.0);
      std::array<size_t,3> sample_dims = {std::max((size_t)1,(size_t)(chunk_idx[i][1]*length_sample_rate)), std::max((size_t)1,(size_t)(chunk_idx[i][3]*length_sample_rate)), std::max((size_t)1,(size_t)(chunk_idx[i][5]*length_sample_rate))};
      size_t sample_num = sample_dims[0]*sample_dims[1]*sample_dims[2];

      auto sampled_data = m_sample_center(chunk,chunk_dims,sample_dims);
      /*
      if(qoi_block_size > 1){//regional 
        //adjust qoieb
        //compressor->set_qoi_tol( qoi_tol);
        //compressor->set_qoi_block_size(qoi_block_size);
        double rate = 1.0;
      

       
        //conf.regionalQoIeb=conf.qoiEB;//store original regional eb
        double num_blocks = 1;
        double num_elements = 1;
        for(int i=0; i<m_dims.size(); i++){
            num_elements *= qoi_block_size;
            num_blocks *= (m_dims[i] - 1) / qoi_block_size + 1;
        }

        double q = 0.999999;
        rate = estimate_rate_Hoeffdin(num_elements,num_blocks,q,qoi_k);
        //std::cout<<num_elements<<" "<<num_blocks<<" "<<conf.error_std_rate<<" "<<rate<<std::endl;
        
        rate = std::max(1.0,rate);//only effective for average. general: 1.0/sumai
        
        //std::cout<<"Pointwise QoI eb rate: " << rate << std::endl;
        qoi_tol *= rate;
        //qoi->set_qoi_tolerance(qoi_tol);
      } 
      */
      m_qoi = QoZ::GetQOI<double>(qoi_id, qoi_tol, pwe, qoi_string);

      std::array<std::vector<double>,3> ebs;
      for(auto i:{0,1,2})
        ebs[i].resize(chunk_ele_num);
    // use quantile to determine abs bound
  


        

        for (size_t i = 0; i < chunk_ele_num; i++){
          auto cur_ebs = m_qoi->interpret_eb(chunk[0][i],chunk[1][i],chunk[2][i]);
          for(auto j:{0,1,2})
            ebs[j][i]=cur_ebs[j];
        }
        
        //double max_quantile_rate = 0.2;
        double quantile_rate = 0.2  ;//conf.quantile;//quantile
        //std::cout<<quantile<<std::endl;
        size_t k = std::ceil(quantile_rate * chunk_ele_num);
        k = std::max((size_t)1, std::min(chunk_ele_num, k)); 

        std::vector<size_t> quantiles;
      
       for(auto i:{1.0,0.5,0.25,0.10,0.05,0.025,0.01})
           quantiles.push_back((size_t)(i*k));
       int quantile_num = quantiles.size();




       std::array<double,3> best_abs_eb=pwe;

      
        

            

            
            //std::sort(ebs.begin(),ebs.begin()+k+1);


       
  
        size_t best_quantile = 0;

         for(auto i:{0,1,2})
            std::nth_element(ebs[i].begin(),ebs[i].begin()+quantiles[0], ebs[i].end());

        size_t last_quantile = quantiles[0]+1;

        double best_br = 9999;
        //best_abs_eb = pwe;
                        
        int idx = 0;
        for(auto quantile:quantiles)
        {   
            if(idx!=0){
              for(auto i:{0,1,2})
                std::nth_element(ebs[i].begin(),ebs[i].begin()+quantile, ebs[i].begin()+last_quantile);
            }

            
            std::array<double,3> cur_abs_eb = {ebs[0][quantile],ebs[1][quantile],ebs[2][quantile]};
            //qoi->set_global_eb(cur_abs_eb);
            // reset variables for average of square
            std::array<std::unique_ptr<SPECK3D_FLT>,3> test_compressor = {std::make_unique<SPECK3D_FLT>(),std::make_unique<SPECK3D_FLT>(),std::make_unique<SPECK3D_FLT>()};
            auto sampled_copy = sampled_data;
            for(auto i:{0,1,2}){
              test_compressor[i]->take_data(std::move(sampled_copy[i]));
              test_compressor[i]->set_dims(sample_dims);
              test_compressor[i]->set_tolerance(cur_abs_eb[i]);
            }
            //test_compressor->set_qoi(qoi);


            std::array<vec8_type,3> test_encoded_stream;
            for(auto i:{0,1,2}){
              auto rtn = test_compressor[i]->compress();
              if(rtn!= RTNType::Good)
                std::cout<<"Error"<<std::endl;
            }

            std::array< vecd_type,3> offsets;
            for(auto i:{0,1,2}){
              offsets[i].resize(sample_num,0);
            }
            std::array<const vecd_type *,3>sampled_ori = {&test_compressor[0]->view_orig_data(),&test_compressor[1]->view_orig_data(),&test_compressor[2]->view_orig_data()}; 
            std::array<double,3>test_means = {test_compressor[0]->get_mean(),test_compressor[1]->get_mean(),test_compressor[2]->get_mean()};
            //std::cout<<test_means[0]<<" "<<test_means[1]<<" "<<test_means[2]<<std::endl;
            std::array<const vecd_type *,3>sampled_dec = {&test_compressor[0]->view_decoded_data(),&test_compressor[1]->view_decoded_data(),&test_compressor[2]->view_decoded_data()}; 
            

            bool outlier = false;
            size_t outlier_num = 0;
            //std::cout<<m_qoi->get_qoi_tolerance()<<std::endl;
            for(size_t i = 0; i < sample_num ; i++){
              //std::cout<<sampled_data[0][i]<<" "<<(*sampled_dec[0])[i]<<std::endl;
              if(!m_qoi->check_compliance((*sampled_ori[0])[i],(*sampled_ori[1])[i],(*sampled_ori[2])[i],
                                                   (*sampled_dec[0])[i]+test_means[0],(*sampled_dec[1])[i]+test_means[1],(*sampled_dec[2])[i]+test_means[2]) ){
                outlier = true;
                outlier_num++;
                for(auto j:{0,1,2})
                  offsets[j][i] = (*sampled_ori[j])[i] - ( (*sampled_dec[j])[i] + test_means[j] );
              }
            }

            if(outlier){
              for(auto i:{0,1,2})
                test_compressor[i]->zstd_encode(offsets[i]);
            }
            //std::cout<<outlier_num<<" "<<sample_num<<std::endl;

            //todo here

            //
            for(auto i:{0,1,2}){
              test_encoded_stream[i].clear();
              test_encoded_stream[i].reserve(128);
              test_compressor[i]->append_encoded_bitstream(test_encoded_stream[i]);
            }
            //m_encoded_streams[i].reserve(1280000);
            
            double cur_br = 0;
            for(auto i:{0,1,2})
              cur_br += test_encoded_stream[i].size()*8.0/(double)(3*sample_num);       
            std::cout << "current_eb = " << cur_abs_eb[0] <<" "<< cur_abs_eb[1]<<" "<< cur_abs_eb[2] << ", current_br = " << cur_br << std::endl;
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
            if(!outlier)
              break;
            
        }
        std::cout<<"Selected quantile: "<<(double)best_quantile/(double)chunk_ele_num<<std::endl;
        std::cout << "Best abs eb:  " << best_abs_eb[0] <<" "<< best_abs_eb[1]<<" "<< best_abs_eb[2] << std::endl; 
        //qoi->set_global_eb(best_abs_eb); 

        //compressor->set_qoi(qoi);

        m_quality = best_abs_eb;
        


    }
    







      //compressor->set_qoi(qoi);
      //compressor->set_qoi_tol(qoi_tol);
    


    
    //std::cout<<qoi_tol<<" "<<qoi_id<<std::endl;
    /*
    if(qoi_id>0 and qoi_tol>0){
      auto pwe = m_mode == CompMode::PWE ? m_quality : std::numeric_limits<double>::max();
      auto qoi = QoZ::GetQOI<double>(qoi_id, qoi_tol, pwe, qoi_string );
      compressor->set_qoi(qoi);
      //compressor->set_qoi_tol(qoi_tol);
    }*/
    //std::cout<<chunk_idx[i][1]<<" "<<chunk_idx[i][3]<<" "<<chunk_idx[i][5]<<std::endl;//its reversed (fastest first)
    for(auto j:{0,1,2}){
      compressor[j]->take_data(std::move(chunk[j]));
      compressor[j]->set_dims({chunk_idx[i][1], chunk_idx[i][3], chunk_idx[i][5]});
      switch (m_mode) {
        case CompMode::PSNR:
          compressor[j]->set_psnr(m_quality[j]);
          break;
        case CompMode::PWE:
          compressor[j]->set_tolerance(m_quality[j]);
          break;
        case CompMode::Rate:
          compressor[j]->set_bitrate(m_quality[j]);
          break;
  #ifdef EXPERIMENTING
        case CompMode::DirectQ:
          compressor[j]->set_direct_q(m_quality[j]);
          break;
  #endif
        default:;  // So the compiler doesn't complain about missing cases.
      }
      chunk_rtn[i] = compressor[j]->compress();
    }
    if(qoi_id>0 and qoi_tol>0){
      std::array<const vecd_type *,3>orig_data = {&compressor[0]->view_orig_data(),&compressor[1]->view_orig_data(),&compressor[2]->view_orig_data()};
      std::array<const vecd_type *,3>dec_data = {&compressor[0]->view_decoded_data(),&compressor[1]->view_decoded_data(),&compressor[2]->view_decoded_data()}; 
      std::array<double,3>means = {compressor[0]->get_mean(),compressor[1]->get_mean(),compressor[2]->get_mean()};
      std::array< vecd_type,3> offsets;
      for(auto j:{0,1,2}){
        offsets[j].resize(chunk_ele_num,0);
      }
      bool outlier = false;
      //std::cout<<chunk_ele_num<<std::endl;
      //std::cout<<m_qoi->get_qoi_tolerance()<<std::endl;
      for(size_t k = 0; k < chunk_ele_num ; k++){
        if(!m_qoi->check_compliance((*orig_data[0])[k],(*orig_data[1])[k],(*orig_data[2])[k],
                                             (*dec_data[0])[k]+means[0],(*dec_data[1])[k]+means[1],(*dec_data[2])[k]+means[2]) ){
          outlier = true;
          for(auto j:{0,1,2})
            offsets[j][k] = (*orig_data[j])[k] - ( (*dec_data[j])[k] + means[j]);
        }
      }
      if(outlier){
        for(auto j:{0,1,2})
          compressor[j]->zstd_encode(offsets[j]);
      }
      

    }

    // Save bitstream for each chunk in `m_encoded_stream`.
    for(auto j:{0,1,2}){
      m_encoded_streams[j][i].clear();
      m_encoded_streams[j][i].reserve(128);
      //m_encoded_streams[i].reserve(1280000);
      compressor[j]->append_encoded_bitstream(m_encoded_streams[j][i]);
    }
  }

  auto fail = std::find_if_not(chunk_rtn.begin(), chunk_rtn.end(),
                               [](auto r) { return r == RTNType::Good; });
  if (fail != chunk_rtn.end())
    return (*fail);
  for(auto j:{0,1,2})
    assert(std::none_of(m_encoded_streams[j].cbegin(), m_encoded_streams[j].cend(),
                      [](auto& s) { return s.empty(); }));

  return RTNType::Good;
}
template auto sperr::SPERR3D_VEC_OMP_C::compress(const float*,const float*,const float*, size_t) -> RTNType;
template auto sperr::SPERR3D_VEC_OMP_C::compress(const double*,const double*,const double*, size_t) -> RTNType;

auto sperr::SPERR3D_VEC_OMP_C::get_encoded_bitstream() const -> std::array<vec8_type,3>
{
  std::array<vec8_type,3> bitstream;
  for(auto i:{0,1,2}){
    auto header = m_generate_header(i);
    assert(!header.empty());
    auto header_size = header.size();
    auto stream_size = std::accumulate(m_encoded_streams[i].cbegin(), m_encoded_streams[i].cend(), 0lu,
                                       [](size_t a, const auto& b) { return a + b.size(); });
    header.resize(header_size + stream_size);

    auto itr = header.begin() + header_size;
    for (const auto& s : m_encoded_streams[i]) {
      std::copy(s.cbegin(), s.cend(), itr);
      itr += s.size();
    }
    bitstream[i] = header;
  }

  return bitstream;
}

auto sperr::SPERR3D_VEC_OMP_C::m_generate_header(int i) const -> sperr::vec8_type
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
  if (num_chunks != m_encoded_streams[i].size())
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
  for (const auto& stream : m_encoded_streams[i]) {
    assert(stream.size() <= uint64_t{std::numeric_limits<uint32_t>::max()});
    uint32_t len = stream.size();
    std::memcpy(&header[pos], &len, sizeof(len));
    pos += sizeof(len);
  }
  assert(pos == header_size);

  return header;
}

template <typename T>
auto sperr::SPERR3D_VEC_OMP_C::m_gather_chunk(const T* vol,
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
template auto sperr::SPERR3D_VEC_OMP_C::m_gather_chunk(const float*,
                                                   dims_type,
                                                   std::array<size_t, 6>) -> vecd_type;
template auto sperr::SPERR3D_VEC_OMP_C::m_gather_chunk(const double*,
                                                   dims_type,
                                                   std::array<size_t, 6>) -> vecd_type;



auto sperr::SPERR3D_VEC_OMP_C::m_sample_center(std::array<vecd_type,3> chunk,std::array<size_t, 3> chunk_dims,std::array<size_t,3>sample_dims) -> std::array<vecd_type,3>
{
  std::array<size_t,3>starts= {(chunk_dims[0]-sample_dims[0])/2,(chunk_dims[1]-sample_dims[1])/2,(chunk_dims[2]-sample_dims[2])/2};
  //std::array<size_t,3>ends= {(chunk_dim[0]+sample_dim[0])/2,(chunk_dim[1]+sample_dim[1])/2,(chunk_dim[2]+sample_dim[2])/2};
  std::array<vecd_type,3> sampled_data;
  for(auto i:{0,1,2}){
    size_t idx = 0;
    size_t y_offset = chunk_dims[0];
    size_t z_offset = y_offset * chunk_dims[1];
    for (size_t z = starts[2]; z < starts[2]+sample_dims[2]; z++) {
      const size_t plane_offset = z * z_offset;
      for (size_t y = starts[1]; y < starts[1]+sample_dims[1]; y++) {
        const auto start_idx = plane_offset + y * y_offset + starts[0];
        sampled_data[i].insert(sampled_data[i].end(), chunk[i].begin()+start_idx, chunk[i].begin()+start_idx+sample_dims[0]);
      }
    }
  }
  return sampled_data;

}

double sperr::SPERR3D_VEC_OMP_C::estimate_rate_Hoeffdin(size_t n, size_t N, double q, double k = 2.0){//n: element_per_block N: num_blocks q: confidence
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
