#ifndef SZ3_QOI_INTERFACE
#define SZ3_QOI_INTERFACE
/**
 * Interface for some specific quantities of interest (QoIs)
 */

#include<vector>
#include<string>
namespace QoZ {


    namespace concepts {

        template<class T>
        class QoIInterface {
        public:

           // using Range = multi_dimensional_range<T, N>;
            //using iterator = typename multi_dimensional_range<T, N>::iterator;

            virtual ~QoIInterface() = default;

            virtual std::array<T,3> interpret_eb(T x, T y, T z) const = 0;

            // interpret eb with iterator (Lorenzo)
            //virtual T interpret_eb(const iterator &iter) const = 0;

            // interpret eb with data pointer (Interpolation)
            //virtual T interpret_eb(const T * data, ptrdiff_t offset) = 0;

            //virtual void update_tolerance(T data, T dec_data) = 0;

            virtual bool check_compliance(T x, T y, T z, T dec_x, T dec_y, T dec_z) const = 0;

            //virtual void precompress_block(const std::shared_ptr<Range> &range) = 0;

            //virtual void postcompress_block() = 0;

            //virtual void pre_compute(const T * data) = 0;

            virtual void print() = 0;

            // for interpolation compressors
          //  virtual T get_global_eb() const = 0;

            //virtual void set_global_eb(T eb1, T eb2, T eb3) = 0;

            virtual void init() = 0;

            virtual double eval(T x,T y, T z) const = 0;

            virtual void set_qoi_tolerance(double tol) = 0;

            virtual double get_qoi_tolerance() = 0;

            virtual std::string get_expression() const = 0;

            virtual void set_dims(const std::vector<size_t>& new_dims) = 0;


            int id = 0;
        };

    }

}

#endif
