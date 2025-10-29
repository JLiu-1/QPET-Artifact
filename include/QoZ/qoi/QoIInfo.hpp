#ifndef SZ3_QOI_INFO
#define SZ3_QOI_INFO

#include "QoI.hpp"
#include "XSquare.hpp"
#include "XCubic.hpp"
#include "XPower.hpp"
#include "XSin.hpp"
#include "XExp.hpp"
#include "XLin.hpp"
#include "XSqrt.hpp"
#include "XAbs.hpp"
#include "XTanh.hpp"
#include "XReciprocal.hpp"
#include "XComposite.hpp"
#include "LogX.hpp"
#include "XLogX.hpp"
#include "RegionalAverage.hpp"
#include "RegionalAverageOfSquare.hpp"
#include "Isoline.hpp"
#include "MultiQoIs.hpp"
#include "FX.hpp"
#include "FX_abs.hpp"
#include "FX_New.hpp"
#include "FX_P.hpp"
#include "RegionalFX.hpp"
#include <vector>

namespace QoZ {

    template<class T, QoZ::uint N >
    std::shared_ptr<concepts::QoIInterface<T, N>> GetQOI(const Config &conf){
        auto qoi_id = conf.qoi;
        auto qoiLogBase = conf.qoiLogBase;
        auto qoi_string = conf.qoi_string;
        auto qoi_string_2 = conf.qoi_string_2;
        if(qoi_id == 14 and conf.QoIdispatch){
            if(qoi_string == "x^2" or qoi_string == "x**2")
                qoi_id = 1;
            else if(qoi_string == "logx" or qoi_string == "log(x)" or qoi_string == "Logx" or qoi_string == "Log(x)" or qoi_string == "lnx" or qoi_string == "ln(x)" or qoi_string == "Lnx" or qoi_string == "Ln(x)"){
                qoi_id = 2;
                qoiLogBase = std::exp(1.0);
            }
            else if(qoi_string == "log2(x)"  or qoi_string == "Log2(x)"  or qoi_string == "Log(x,2)"  or qoi_string == "log(x,2)"){
                qoi_id = 2;
                qoiLogBase = 2.0;
            }

            else if (qoi_string == "x^3" or qoi_string == "x**3")
                qoi_id = 9;
            else if (qoi_string == "x^0.5" or qoi_string == "x**0.5" or qoi_string == "x**1/2" or qoi_string == "x^1/2" or qoi_string == "sqrt(x)" or qoi_string == "Sqrt(x)")
                qoi_id = 10;
            //temp do not dispath lin beacuse it is preprocessed in compressor. Todo: dispatch lin.
            else if (qoi_string == "2^x" or qoi_string == "2**x"){
                qoi_id = 12;
                qoiLogBase = 2.0;
            }
            else if (qoi_string == "e^x" or qoi_string == "e**x" or qoi_string == "E**x" or qoi_string == "E^x"){
                qoi_id = 12;
                qoiLogBase = std::exp(1.0);
            }
            else if (qoi_string == "1/x" or qoi_string == "x**-1" or qoi_string == "x^-1")
                qoi_id = 13;
            //todo: dispatch x^n
            else if (qoi_string == "|x|")
                qoi_id = 19;
            else if (qoi_string == "sinx" or qoi_string == "sin(x)" or qoi_string == "Sinx" or qoi_string == "Sin(x)")
                qoi_id = 22;
            else if (qoi_string == "sin10x" or qoi_string == "sin(10x)" or qoi_string == "Sin10x" or qoi_string == "Sin(10x)"){
                qoi_id = 20;
                qoi_string = "11 10 0 22";
            }
            else if (qoi_string == "tanhx" or qoi_string == "tanh(x)" or qoi_string == "Tanhx" or qoi_string == "Tanh(x)")
                qoi_id = 23;

        }
        switch(qoi_id){
            case 1:
                return std::make_shared<QoZ::QoI_X_Square<T, N>>(conf.qoiEB, conf.absErrorBound);
            case 2:{
                if(conf.analytical)
                    return std::make_shared<QoZ::QoI_Log_X<T, N>>(conf.qoiEB, conf.absErrorBound, qoiLogBase);
                else
                    return std::make_shared<QoZ::QoI_Log_X_Approx<T, N>>(conf.qoiEB, conf.absErrorBound, qoiLogBase);
            }
            case 3:{
                if(!conf.lorenzo && !conf.lorenzo2) return std::make_shared<QoZ::QoI_RegionalAverageOfSquareInterp<T, N>>(conf.qoiEB, conf.absErrorBound, conf.qoiRegionSize, conf.dims);
                else return std::make_shared<QoZ::QoI_RegionalAverageOfSquare<T, N>>(conf.qoiEB, conf.absErrorBound);
                // return std::make_shared<SZ::QoI_RegionalAverage<T, N>>(conf.qoiEB, conf.absErrorBound);
            }
            case 4:{
            	std::vector<T> values;
            	for(int i=0; i<conf.isovalues.size(); i++){
            		values.push_back(conf.isovalues[i]);
            	}
                return std::make_shared<QoZ::QoI_Isoline<T, N>>(conf.dims, values, conf.absErrorBound);            	
            }
            case 5:{
            	// x^2 + log x
            	std::vector<std::shared_ptr<concepts::QoIInterface< T, N>>> qois;
            	qois.push_back(std::make_shared<QoZ::QoI_X_Square<T, N>>(conf.qoiEBs[0], conf.absErrorBound));
            	qois.push_back(std::make_shared<QoZ::QoI_Log_X<T, N>>(conf.qoiEBs[1], conf.absErrorBound));
                return std::make_shared<QoZ::QoI_MultiQoIs<T, N>>(qois);            	
            }
            case 6:{
            	// x^2 + isoline
            	std::vector<std::shared_ptr<concepts::QoIInterface< T, N>>> qois;
            	qois.push_back(std::make_shared<QoZ::QoI_X_Square<T, N>>(conf.qoiEBs[0], conf.absErrorBound));
                std::vector<T> values;
                for(int i=0; i<conf.isovalues.size(); i++){
                    values.push_back(conf.isovalues[i]);
                }
                qois.push_back(std::make_shared<QoZ::QoI_Isoline<T, N>>(conf.dims, values, conf.absErrorBound));
                return std::make_shared<QoZ::QoI_MultiQoIs<T, N>>(qois);            	
            }
            case 7:{
                // log x + isoline
                std::vector<std::shared_ptr<concepts::QoIInterface< T, N>>> qois;
                qois.push_back(std::make_shared<QoZ::QoI_Log_X<T, N>>(conf.qoiEBs[0], conf.absErrorBound));
                std::vector<T> values;
                for(int i=0; i<conf.isovalues.size(); i++){
                    values.push_back(conf.isovalues[i]);
                }
                qois.push_back(std::make_shared<QoZ::QoI_Isoline<T, N>>(conf.dims, values, conf.absErrorBound));
                return std::make_shared<QoZ::QoI_MultiQoIs<T, N>>(qois);             
            }
            case 8:{
            	// x^2 + log x + isoline
            	std::vector<std::shared_ptr<concepts::QoIInterface< T, N>>> qois;
            	qois.push_back(std::make_shared<QoZ::QoI_X_Square<T, N>>(conf.qoiEBs[0], conf.absErrorBound));
                qois.push_back(std::make_shared<QoZ::QoI_Log_X<T, N>>(conf.qoiEBs[1], conf.absErrorBound));
            	std::vector<T> values;
            	for(int i=0; i<conf.isovalues.size(); i++){
            		values.push_back(conf.isovalues[i]);
            	}
                qois.push_back(std::make_shared<QoZ::QoI_Isoline<T, N>>(conf.dims, values, conf.absErrorBound));
                return std::make_shared<QoZ::QoI_MultiQoIs<T, N>>(qois);            	
            }
            case 9:{
                if(conf.analytical)
                    return std::make_shared<QoZ::QoI_X_Cubic<T, N>>(conf.qoiEB, conf.absErrorBound);
                else
                    return std::make_shared<QoZ::QoI_X_Cubic_Approx<T, N>>(conf.qoiEB, conf.absErrorBound);
            }       
            case 10:{
                //return std::make_shared<QoZ::QoI_X_Sin<T, N>>(conf.qoiEB, conf.absErrorBound);
                if(conf.analytical)
                    return std::make_shared<QoZ::QoI_X_Sqrt<T, N>>(conf.qoiEB, conf.absErrorBound);
                else
                    return std::make_shared<QoZ::QoI_X_Sqrt_Approx<T, N>>(conf.qoiEB, conf.absErrorBound);
            }
            case 11:
                return std::make_shared<QoZ::QoI_X_Lin<T, N>>(conf.qoiEB, conf.absErrorBound, conf.qoi_lin_A, conf.qoi_lin_B);
            case 12:{
                if(conf.analytical)
                    return std::make_shared<QoZ::QoI_X_Exp<T, N>>(conf.qoiEB, conf.absErrorBound,qoiLogBase);
                else
                    return std::make_shared<QoZ::QoI_X_Exp_Approx<T, N>>(conf.qoiEB, conf.absErrorBound,qoiLogBase);
            }
            case 13:{
                //return std::make_shared<QoZ::QoI_XLog_X<T, N>>(conf.qoiEB, conf.absErrorBound);
                if(conf.analytical)
                    return std::make_shared<QoZ::QoI_X_Recip<T, N>>(conf.qoiEB, conf.absErrorBound);
                else
                    return std::make_shared<QoZ::QoI_X_Recip_Approx<T, N>>(conf.qoiEB, conf.absErrorBound);
            }
            case 14:
                return std::make_shared<QoZ::QoI_FX<T, N>>(conf.qoiEB, conf.absErrorBound, qoi_string, conf.isolated, conf.threshold);
            case 15:
                return std::make_shared<QoZ::QoI_FX_P<T, N>>(conf.qoiEB, conf.absErrorBound, qoi_string, qoi_string_2, conf.threshold, conf.isolated);

            case 16:{
                return std::make_shared<QoZ::QoI_RegionalFX<T, N>>(conf.qoiEB, conf.absErrorBound, conf.qoiRegionSize, conf.dims, qoi_string, conf.error_std_rate, conf.isolated, conf.threshold);
                
                // return std::make_shared<QoZ::QoI_RegionalAverage<T, N>>(conf.qoiEB, conf.absErrorBound);
            }
            case 17:
                return std::make_shared<QoZ::QoI_FX_ABS<T, N>>(conf.qoiEB, conf.absErrorBound, qoi_string, conf.isolated, conf.threshold);
            case 18:{
                if(conf.analytical)
                    return std::make_shared<QoZ::QoI_X_Power<T, N>>(conf.qoiEB, conf.absErrorBound,qoiLogBase);
                else
                    return std::make_shared<QoZ::QoI_X_Power_Approx<T, N>>(conf.qoiEB, conf.absErrorBound,qoiLogBase);

            }
            case 19:
                return std::make_shared<QoZ::QoI_X_Abs<T, N>>(conf.qoiEB, conf.absErrorBound);
            case 20:
                return std::make_shared<QoZ::QoI_X_Composite<T, N>>(conf.qoiEB, conf.absErrorBound, qoi_string, conf.analytical);
            case 21:
                return std::make_shared<QoZ::QoI_FX_New<T, N>>(conf.qoiEB, conf.absErrorBound, qoi_string, conf.isolated, conf.threshold);
            case 22:{
                if(conf.analytical)
                    return std::make_shared<QoZ::QoI_X_Sin<T, N>>(conf.qoiEB, conf.absErrorBound);
                else
                    return std::make_shared<QoZ::QoI_X_Sin_Approx<T, N>>(conf.qoiEB, conf.absErrorBound);
            }
            case 23:{
                //std::cout<<conf.analytical<<std::endl;
                if(conf.analytical)
                    return std::make_shared<QoZ::QoI_X_Tanh<T, N>>(conf.qoiEB, conf.absErrorBound);
                else
                    return std::make_shared<QoZ::QoI_X_Tanh_Approx<T, N>>(conf.qoiEB, conf.absErrorBound);
            }
        }
        return NULL;
    }

}
#endif
