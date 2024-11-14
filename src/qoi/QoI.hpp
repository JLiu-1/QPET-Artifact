#ifndef SZ3_QOI_INTERFACE
#define SZ3_QOI_INTERFACE
/**
 * Interface for some specific quantities of interest (QoIs)
 * Created by Xin Liang on 12/06/2021.
 */


namespace QoZ {


    namespace concepts {

        template<class T>
        class QoIInterface {
        public:

            virtual ~QoIInterface() = default;

            virtual T interpret_eb(T data) const = 0;

            // interpret eb with iterator (Lorenzo)

            // interpret eb with data pointer (Interpolation)
            virtual T interpret_eb(const T * data, ptrdiff_t offset) = 0;

            virtual void update_tolerance(T data, T dec_data) = 0;

            virtual bool check_compliance(T data, T dec_data, bool verbose=false) const = 0;

            virtual void precompress_block() = 0;

            virtual void postcompress_block() = 0;

            virtual void pre_compute(const T * data) = 0;

            virtual void print() = 0;

            // for interpolation compressors
            virtual T get_global_eb() const = 0;

            virtual void set_global_eb(T eb) = 0;

            virtual void set_qoi_tolerance(double tol) = 0;

            virtual void init() = 0;

            virtual double eval(T val) const = 0;

            virtual std::string get_expression() const = 0;

            virtual void set_dims(const std::vector<size_t>& new_dims) = 0;


            int id = 0;
        };

    }

}

#endif