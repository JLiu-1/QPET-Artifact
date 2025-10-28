#ifndef SZ_QOI_SE_HPP
#define SZ_QOI_SE_HPP

#include <algorithm>
#include <cmath>
#include <functional>
#include <symengine/expression.h>
#include <symengine/parser.h>
#include <symengine/symbol.h>
#include <symengine/derivative.h>
#include <symengine/eval.h> 
#include <symengine/solve.h>
#include <symengine/functions.h>
#include <set>

using SymEngine::Expression;
using SymEngine::Symbol;
using SymEngine::symbol;
using SymEngine::parse;
using SymEngine::diff;
using SymEngine::Constant;
using SymEngine::RealDouble;
using SymEngine::Integer;
using SymEngine::Rational;
using SymEngine::evalf;
using SymEngine::map_basic_basic;
using SymEngine::down_cast;
using SymEngine::RCP;
using SymEngine::Basic;
using SymEngine::real_double;
using SymEngine::eval_double;

using SymEngine::E;
using SymEngine::eq;
using SymEngine::solve;
using SymEngine::rcp_static_cast;
using SymEngine::Mul;
using SymEngine::Pow;
using SymEngine::Log;
using SymEngine::sqrt;
using SymEngine::is_a;
using SymEngine::FiniteSet;


inline bool is_number(const Basic & expr){
    return is_a<const RealDouble>(expr) or SymEngine::is_a<const Integer>(expr) or SymEngine::is_a<const Rational>(expr);
}


inline std::function<double(double, double, double)> convert_expression_to_function(const Basic &expr, 
                                                              const RCP<const Symbol> &x, 
                                                              const RCP<const Symbol> &y, 
                                                              const RCP<const Symbol> &z) {
            
        if (is_a<const SymEngine::Symbol>(expr) && expr.__eq__(*x)) {
            return [](double x_value, double, double) { /*std::cout<<"x="<<x_value<<std::endl;*/return x_value; };
        }
      
        else if (is_a<const SymEngine::Symbol>(expr) && expr.__eq__(*y)) {
            return [](double, double y_value, double) { /*std::cout<<"y="<<y_value<<std::endl;*/return y_value; };
        }
       
        else if (is_a<const SymEngine::Symbol>(expr) && expr.__eq__(*z)) {
            return [](double, double, double z_value) { /*std::cout<<"z="<<z_value<<std::endl;*/return z_value; };
        }
      
        else if (is_number(expr)) {
            double constant_value = eval_double(expr);
            return [constant_value](double, double, double) { /*std::cout<<"c="<<constant_value<<std::endl;*/return constant_value; };
        }

        // E
        else if ( eq(expr,*E)) {
            double e = std::exp(1);
            return [e](double, double, double) { return e; };
        }

        else if ( is_a<SymEngine::Abs>(expr)) {

            auto expr_arg = Expression(expr.get_args()[0]);

            if(is_number(expr_arg)){
                double constant_value = std::abs(eval_double(expr_arg));
                return [constant_value](double,double,double) { return constant_value; };
            }
            else if (is_a<const SymEngine::Symbol>(expr_arg) && expr.__eq__(*x)) {
                return [](double x_value,double,double) { return std::abs(x_value); };
            }
            else if (is_a<const SymEngine::Symbol>(expr_arg) && expr.__eq__(*y)) {
                return [](double ,double y_value,double) { return std::abs(y_value); };
            }
            else if (is_a<const SymEngine::Symbol>(expr_arg) && expr.__eq__(*z)) {
                return [](double ,double,double z_value) { return std::abs(z_value); };
            }

            auto arg = convert_expression_to_function(expr_arg, x,y,z);
            return [arg](double x_value,double y_value, double z_value) {
                return std::abs(arg(x_value,y_value,z_value));
            };

        }

       
        else if (is_a<SymEngine::Add>(expr)) {
            auto args = expr.get_args();

            if(args.size()==2){
                auto expr_arg1 = Expression(args[0]);
                auto expr_arg2 = Expression(args[1]);
                if(is_number(expr_arg1)){
                    if(is_number(expr_arg2)){
                        double constant_value = eval_double(expr_arg1)+eval_double(expr_arg2);
                        return [constant_value](double, double, double) { return constant_value; };
                    }
                    else if(is_a<const SymEngine::Symbol>(expr_arg2)){


                        double constant_value = eval_double(expr_arg1);
                        if(eq(expr_arg2,*x) )
                            return [constant_value](double x_value, double,double ) {
                                return x_value+constant_value;
                            };
                        if(eq(expr_arg2,*y))
                            return [constant_value](double , double y_value,double ) {
                                return y_value+constant_value;
                            };
                        if(eq(expr_arg2,*z))
                            return [constant_value](double , double,double z_value) {
                                return z_value+constant_value;
                            };
                    }
                    else{
                        double constant_value = eval_double(expr_arg1);
                        auto f = convert_expression_to_function(expr_arg2, x,y,z);
                        return [f,constant_value](double x_value,double y_value, double z_value) {
                            return f(x_value,y_value,z_value)+constant_value;
                        };
                    }
                }
                else if (is_a<const SymEngine::Symbol>(expr_arg1)){
                    if(eq(expr_arg1,*x)){
                        if(is_number(expr_arg2)){
                            double constant_value = eval_double(expr_arg2);
                            return [constant_value](double x_value, double, double) { return x_value+constant_value; };
                        }
                        else if(is_a<const SymEngine::Symbol>(expr_arg2)){

                            if(eq(expr_arg2,*x))
                            return [](double x_value, double,double ) {
                                return x_value+x_value;
                            };
                            if(eq(expr_arg2,*y))
                                return [](double x_value, double y_value,double ) {
                                    return x_value+y_value;
                                };
                            if(eq(expr_arg2,*z))
                                return [](double x_value, double,double z_value) {
                                    return x_value+z_value;
                                };
                        }
                        else{
                            auto f = convert_expression_to_function(expr_arg2, x,y,z);
                            return [f](double x_value,double y_value, double z_value) {
                                return f(x_value,y_value,z_value)+x_value;
                            };
                        }
                    }
                    if(eq(expr_arg1,*y)){
                        if(is_number(expr_arg2)){
                            double constant_value = eval_double(expr_arg2);
                            return [constant_value](double , double y_value, double) { return y_value+constant_value; };
                        }
                        else if(is_a<const SymEngine::Symbol>(expr_arg2)){

                            if(eq(expr_arg2,*x))
                            return [](double x_value, double y_value,double ) {
                                return x_value+y_value;
                            };
                            if(eq(expr_arg2,*y))
                                return [](double , double y_value,double ) {
                                    return y_value+y_value;
                                };
                            if(eq(expr_arg2,*z))
                                return [](double , double y_value,double z_value) {
                                    return y_value+z_value;
                                };
                        }
                        else{
                            auto f = convert_expression_to_function(expr_arg2, x,y,z);
                            return [f](double x_value,double y_value, double z_value) {
                                return f(x_value,y_value,z_value)+y_value;
                            };
                        }
                    }
                    if(eq(expr_arg1,*z)){
                        if(is_number(expr_arg2)){
                            double constant_value = eval_double(expr_arg2);
                            return [constant_value](double , double, double z_value) { return z_value+constant_value; };
                        }
                        else if(is_a<const SymEngine::Symbol>(expr_arg2)){

                            if(eq(expr_arg2,*x))
                            return [](double x_value, double,double z_value) {
                                return x_value+z_value;
                            };
                            if(eq(expr_arg2,*y))
                                return [](double, double y_value,double z_value) {
                                    return y_value+z_value;
                                };
                            if(eq(expr_arg2,*z))
                                return [](double , double,double z_value) {
                                    return z_value+z_value;
                                };
                        }
                        else{
                            auto f = convert_expression_to_function(expr_arg2, x,y,z);
                            return [f](double x_value,double y_value, double z_value) {
                                return f(x_value,y_value,z_value)+z_value;
                            };
                        }
                    }
                }

                else{
                    auto f = convert_expression_to_function(expr_arg1, x,y,z);
                    if(is_number(expr_arg2)){

                        double constant_value = eval_double(expr_arg2);
                        return [f,constant_value](double x_value,double y_value, double z_value) { return f(x_value,y_value,z_value)+constant_value; };
                    }
                    else if(is_a<const SymEngine::Symbol>(expr_arg2)){

                        if(eq(expr_arg2,*x))
                            return [f](double x_value,double y_value, double z_value) {
                                return f(x_value,y_value,z_value)+x_value;
                            };
                        if(eq(expr_arg2,*y))
                            return [f](double x_value,double y_value, double z_value) {
                                return f(x_value,y_value,z_value)+y_value;
                            };;
                        if(eq(expr_arg2,*z))
                            return [f](double x_value,double y_value, double z_value) {
                                return f(x_value,y_value,z_value)+z_value;
                            };


                        
                    }
                    else{
                        auto f2 = convert_expression_to_function(expr_arg2, x,y,z);
                        return [f,f2](double x_value,double y_value, double z_value) {
                            return f(x_value,y_value,z_value)+f2(x_value,y_value,z_value);
                        };
                    }

                }
            }




           // auto first = convert_expression_to_function(Expression(args[0]), x, y, z);

          
            std::vector<std::function<double(double, double, double)> > fs;
            double constant_value = 0;
            for (size_t i = 0; i < args.size(); ++i) {

                auto expr_arg = Expression(args[i]);
                    if(is_number(expr_arg)){
                        constant_value += eval_double(expr_arg);
                    }
                else

                    fs.push_back(convert_expression_to_function(expr_arg, x, y, z));
            }

           // auto first = convert_expression_to_function(Expression(args[0]), x, y, z);

            return [fs,constant_value](double x_value, double y_value, double z_value) {
                double result = constant_value;
                for (auto &fnc:fs) {
                    result += fnc(x_value, y_value, z_value);
                }
                return result;
            };
        }
        else if (is_a<SymEngine::Mul>(expr)) {
            auto args = expr.get_args();

            if(args.size()==2){
                auto expr_arg1 = Expression(args[0]);
                auto expr_arg2 = Expression(args[1]);
                if(is_number(expr_arg1)){
                    if(is_number(expr_arg2)){
                        double constant_value = eval_double(expr_arg1)*eval_double(expr_arg2);
                        return [constant_value](double, double, double) { return constant_value; };
                    }
                    else if(is_a<const SymEngine::Symbol>(expr_arg2)){


                        double constant_value = eval_double(expr_arg1);
                        if(eq(expr_arg2,*x))
                            return [constant_value](double x_value, double,double ) {
                                return x_value*constant_value;
                            };
                        if(eq(expr_arg2,*y))
                            return [constant_value](double , double y_value,double ) {
                                return y_value*constant_value;
                            };
                        if(eq(expr_arg2,*z))
                            return [constant_value](double , double,double z_value) {
                                return z_value*constant_value;
                            };
                    }
                    else{
                        double constant_value = eval_double(expr_arg1);
                        auto f = convert_expression_to_function(expr_arg2, x,y,z);
                        return [f,constant_value](double x_value,double y_value, double z_value) {
                            return f(x_value,y_value,z_value)*constant_value;
                        };
                    }
                }
                else if (is_a<const SymEngine::Symbol>(expr_arg1)){
                    if(eq(expr_arg1,*x)){
                        if(is_number(expr_arg2)){
                            double constant_value = eval_double(expr_arg2);
                            return [constant_value](double x_value, double, double) { return x_value*constant_value; };
                        }
                        else if(is_a<const SymEngine::Symbol>(expr_arg2)){

                            if(eq(expr_arg2,*x))
                            return [](double x_value, double,double ) {
                                return x_value*x_value;
                            };
                            if(eq(expr_arg2,*y))
                                return [](double x_value, double y_value,double ) {
                                    return x_value*y_value;
                                };
                            else
                                return [](double x_value, double,double z_value) {
                                    return x_value*z_value;
                                };
                        }
                        else{
                            auto f = convert_expression_to_function(expr_arg2, x,y,z);
                            return [f](double x_value,double y_value, double z_value) {
                                return f(x_value,y_value,z_value)*x_value;
                            };
                        }
                    }
                    if(eq(expr_arg1,*y)){
                        if(is_number(expr_arg2)){
                            double constant_value = eval_double(expr_arg2);
                            return [constant_value](double , double y_value, double) { return y_value*constant_value; };
                        }
                        else if(is_a<const SymEngine::Symbol>(expr_arg2)){

                            if(eq(expr_arg2,*x))
                            return [](double x_value, double y_value,double ) {
                                return x_value*y_value;
                            };
                            if(eq(expr_arg2,*y))
                                return [](double , double y_value,double ) {
                                    return y_value*y_value;
                                };
                            else
                                return [](double , double y_value,double z_value) {
                                    return y_value*z_value;
                                };
                        }
                        else{
                            auto f = convert_expression_to_function(expr_arg2, x,y,z);
                            return [f](double x_value,double y_value, double z_value) {
                                return f(x_value,y_value,z_value)*y_value;
                            };
                        }
                    }
                    else{//z
                        if(is_number(expr_arg2)){
                            double constant_value = eval_double(expr_arg2);
                            return [constant_value](double , double, double z_value) { return z_value*constant_value; };
                        }
                        else if(is_a<const SymEngine::Symbol>(expr_arg2)){

                            if(eq(expr_arg2,*x))
                            return [](double x_value, double,double z_value) {
                                return x_value*z_value;
                            };
                            if(eq(expr_arg2,*y))
                                return [](double, double y_value,double z_value) {
                                    return y_value*z_value;
                                };
                            else
                                return [](double , double,double z_value) {
                                    return z_value*z_value;
                                };
                        }
                        else{
                            auto f = convert_expression_to_function(expr_arg2, x,y,z);
                            return [f](double x_value,double y_value, double z_value) {
                                return f(x_value,y_value,z_value)*z_value;
                            };
                        }
                    }
                }

                else{
                    auto f = convert_expression_to_function(expr_arg1, x,y,z);
                    if(is_number(expr_arg2)){

                        double constant_value = eval_double(expr_arg2);
                        return [f,constant_value](double x_value,double y_value, double z_value) { return f(x_value,y_value,z_value)*constant_value; };
                    }
                    else if(is_a<const SymEngine::Symbol>(expr_arg2)){

                        if(eq(expr_arg2,*x))
                            return [f](double x_value,double y_value, double z_value) {
                                return f(x_value,y_value,z_value)*x_value;
                            };
                        if(eq(expr_arg2,*y))
                            return [f](double x_value,double y_value, double z_value) {
                                return f(x_value,y_value,z_value)*y_value;
                            };
                        else
                            return [f](double x_value,double y_value, double z_value) {
                                return f(x_value,y_value,z_value)*z_value;
                            };


                        
                    }
                    else{
                        auto f2 = convert_expression_to_function(expr_arg2, x,y,z);
                        return [f,f2](double x_value,double y_value, double z_value) {
                            return f(x_value,y_value,z_value)*f2(x_value,y_value,z_value);
                        };
                    }

                }
            }


           
            std::vector<std::function<double(double, double, double)> > fs;
            double constant_value = 1.0;
            for (size_t i = 0; i < args.size(); ++i) {

                auto expr_arg = Expression(args[i]);
                    if(is_number(expr_arg)){
                        constant_value *= eval_double(expr_arg);
                    }
                else
                    fs.push_back(convert_expression_to_function(expr_arg, x, y, z));
            }

           // auto first = convert_expression_to_function(Expression(args[0]), x, y, z);

            return [fs,constant_value](double x_value, double y_value, double z_value) {
                double result = constant_value;
                for (auto &fnc:fs) {
                    result *= fnc(x_value, y_value, z_value);
                }
                return result;
            };
        }
       
        else if (is_a<SymEngine::Pow>(expr)) {
            auto args = expr.get_args();

            auto expr_arg1 = Expression(args[0]);
            auto expr_arg2 = Expression(args[1]);
            if(is_number(expr_arg1)){
                if(is_number(expr_arg2)){
                    double constant_value = std::pow(eval_double(expr_arg1),eval_double(expr_arg2));
                    return [constant_value](double,double,double) { return constant_value; };
                }
                else if(is_a<const SymEngine::Symbol>(expr_arg2)){
                    double constant_value = eval_double(expr_arg1);

                    if(eq(expr_arg2,*x))
                        return [constant_value](double x_value,double , double ) {
                            return std::pow(constant_value,x_value);
                        };
                    if(eq(expr_arg2,*y))
                        return [constant_value](double,double y_value, double ) {
                            return std::pow(constant_value,y_value);
                        };
                    else
                        return [constant_value](double,double , double z_value) {
                            return std::pow(constant_value,z_value);
                        };

                }
                else{
                    double constant_value = eval_double(expr_arg1);
                    auto f = convert_expression_to_function(expr_arg2, x,y,z);
                    return [f,constant_value](double x_value, double y_value, double z_value) {
                        return std::pow(constant_value,f(x_value, y_value, z_value));
                    };
                }
            }
            else if (is_a<const SymEngine::Symbol>(expr_arg1)){
                if(eq(expr_arg1,*x)){
                    if(is_number(expr_arg2)){
                        double constant_value = eval_double(expr_arg2);
                        if (constant_value == 1.0 )
                            return [](double x_value,double,double){return x_value;};
                        if (constant_value == 2.0 )
                            return [](double x_value,double,double){return x_value*x_value;};
                        else if (constant_value == 3.0 )
                            return [](double x_value,double,double){return x_value*x_value*x_value;};
                        else if (constant_value == 4.0 )
                            return [](double x_value,double,double){auto a = x_value*x_value; return a*a;};
                        else if (constant_value == 0.5 )
                            return [](double x_value,double,double){return sqrt(x_value);};
                        else
                            return [constant_value](double x_value,double,double) { return std::pow(x_value,constant_value); };
                    }
                    else if(is_a<const SymEngine::Symbol>(expr_arg2)){
                        if(eq(expr_arg2,*x))
                            return [](double x_value,double,double) {
                                return std::pow(x_value,x_value);
                            };
                        if(eq(expr_arg2,*y))
                            return [](double x_value,double y_value,double) {
                                return std::pow(x_value,y_value);
                            };

                        else//z
                            return [](double x_value,double,double z_value) {
                                return std::pow(x_value,z_value);
                            };
                    }
                    else{
                        auto f = convert_expression_to_function(expr_arg2, x,y,z);
                        return [f](double x_value,double y_value,double z_value) {
                            return std::pow(x_value,f(x_value,y_value,z_value));
                        };
                    }
                }
                if(eq(expr_arg1,*y)){
                    if(is_number(expr_arg2)){
                        double constant_value = eval_double(expr_arg2);
                        if (constant_value == 1.0 )
                            return [](double ,double y_value,double){return y_value;};
                        if (constant_value == 2.0 )
                            return [](double ,double y_value,double){return y_value*y_value;};
                        else if (constant_value == 3.0 )
                            return [](double ,double y_value,double){return y_value*y_value*y_value;};
                        else if (constant_value == 4.0 )
                            return [](double ,double y_value,double){auto a = y_value*y_value; return a*a;};
                        else if (constant_value == 0.5 )
                            return [](double ,double y_value,double){return sqrt(y_value);};
                        else
                            return [constant_value](double,double y_value,double) { return std::pow(y_value,constant_value); };
                    }
                    else if(is_a<const SymEngine::Symbol>(expr_arg2)){
                        if(eq(expr_arg2,*x))
                            return [](double x_value,double y_value,double) {
                                return std::pow(y_value,x_value);
                            };
                        if(eq(expr_arg2,*y))
                            return [](double ,double y_value,double) {
                                return std::pow(y_value,y_value);
                            };

                        else//z
                            return [](double ,double y_value,double z_value) {
                                return std::pow(y_value,z_value);
                            };
                    }
                    else{
                        auto f = convert_expression_to_function(expr_arg2, x,y,z);
                        return [f](double x_value,double y_value,double z_value) {
                            return std::pow(y_value,f(x_value,y_value,z_value));
                        };
                    }
                }
                else{//z
                    if(is_number(expr_arg2)){
                        double constant_value = eval_double(expr_arg2);
                        if (constant_value == 1.0 )
                            return [](double ,double ,double z_value){return z_value;};
                        if (constant_value == 2.0 )
                            return [](double ,double ,double z_value){return z_value*z_value;};
                        else if (constant_value == 3.0 )
                            return [](double ,double ,double z_value){return z_value*z_value*z_value;};
                        else if (constant_value == 4.0 )
                            return [](double ,double,double z_value){auto a = z_value*z_value; return a*a;};
                        else if (constant_value == 0.5 )
                            return [](double ,double,double z_value){return sqrt(z_value);};
                        else
                            return [constant_value](double,double ,double z_value) { return std::pow(z_value,constant_value); };
                    }
                    else if(is_a<const SymEngine::Symbol>(expr_arg2)){
                        if(eq(expr_arg2,*x))
                            return [](double x_value,double,double z_value) {
                                return std::pow(z_value,x_value);
                            };
                        if(eq(expr_arg2,*y))
                            return [](double ,double y_value,double z_value) {
                                return std::pow(z_value,y_value);
                            };

                        else//z
                            return [](double ,double,double z_value) {
                                return std::pow(z_value,z_value);
                            };
                    }
                    else{
                        auto f = convert_expression_to_function(expr_arg2, x,y,z);
                        return [f](double x_value,double y_value,double z_value) {
                            return std::pow(z_value,f(x_value,y_value,z_value));
                        };
                    }

                }
            }
            else{
                auto f = convert_expression_to_function(expr_arg1, x,y,z);
                if(is_number(expr_arg2)){

                    double constant_value = eval_double(expr_arg2);
                    if (constant_value == 1.0 )
                        return [f](double x_value,double y_value,double z_value){return f(x_value,y_value,z_value);};
                    if (constant_value == 2.0 )
                        return [f](double x_value,double y_value,double z_value){auto a = f(x_value,y_value,z_value); return a*a;};
                    else if (constant_value == 3.0 )
                        return [f](double x_value,double y_value,double z_value){auto a = f(x_value,y_value,z_value);return a*a*a;};
                    else if (constant_value == 4.0 )
                        return [f](double x_value,double y_value,double z_value){auto a = f(x_value,y_value,z_value);a=a*a; return a*a;};
                    else if (constant_value == 0.5 )
                        return [f](double x_value,double y_value,double z_value){auto a = f(x_value,y_value,z_value);return sqrt(a);};
                    else
                        return [f,constant_value](double x_value,double y_value,double z_value) { return std::pow(f(x_value,y_value,z_value),constant_value); };                
                }
                else if(is_a<const SymEngine::Symbol>(expr_arg2)){
                    if(eq(expr_arg2,*x))
                        return [f](double x_value,double y_value,double z_value) {
                            return std::pow(f(x_value,y_value,z_value),x_value);
                        };
                    if(eq(expr_arg2,*y))
                        return [f](double x_value,double y_value,double z_value) {
                            return std::pow(f(x_value,y_value,z_value),y_value);
                        };

                    else//z
                        return [f](double x_value,double y_value,double z_value) {
                            return std::pow(f(x_value,y_value,z_value),z_value);
                        };

                }
                else{
                    auto f2 = convert_expression_to_function(expr_arg2, x,y,z);
                    return [f,f2](double x_value,double y_value,double z_value) {
                        return std::pow(f(x_value,y_value,z_value),f2(x_value,y_value,z_value));
                    };
                }

            }


            /*
            auto base = convert_expression_to_function(Expression(args[0]), x, y, z);
            auto exponent = convert_expression_to_function(Expression(args[1]), x, y, z);
            return [base, exponent](double x_value, double y_value, double z_value) {
                //std::cout<<"pow"<<std::endl;
                return std::pow(base(x_value, y_value, z_value), exponent(x_value, y_value, z_value));
            };
            */
        }
      
        else if (is_a<SymEngine::Sin>(expr)) {
            auto expr_arg = Expression(expr.get_args()[0]);

            if(is_number(expr_arg)){
                double constant_value = std::sin(eval_double(expr_arg));
                return [constant_value](double,double,double) { return constant_value; };
            }
            else if (is_a<const SymEngine::Symbol>(expr_arg)) {
                if(eq(expr_arg,*x))
                    return [](double x_value,double,double) { return std::sin(x_value); };
                if(eq(expr_arg,*y))
                    return [](double ,double y_value,double) { return std::sin(y_value); };
                else
                    return [](double ,double,double z_value) { return std::sin(z_value); };
            }
            auto arg = convert_expression_to_function(expr_arg, x,y,z);
            return [arg](double x_value, double y_value, double z_value) {
                return std::sin(arg(x_value, y_value, z_value));
            };
        }
        
        else if (is_a<SymEngine::Cos>(expr)) {
            auto expr_arg = Expression(expr.get_args()[0]);

            if(is_number(expr_arg)){
                double constant_value = std::cos(eval_double(expr_arg));
                return [constant_value](double,double,double) { return constant_value; };
            }
            else if (is_a<const SymEngine::Symbol>(expr_arg)) {
                if(eq(expr_arg,*x))
                    return [](double x_value,double,double) { return std::cos(x_value); };
                if(eq(expr_arg,*y))
                    return [](double ,double y_value,double) { return std::cos(y_value); };
                else
                    return [](double ,double,double z_value) { return std::cos(z_value); };
            }
            auto arg = convert_expression_to_function(expr_arg, x,y,z);
            return [arg](double x_value, double y_value, double z_value) {
                return std::cos(arg(x_value, y_value, z_value));
            };
        }
        
        else if (is_a<SymEngine::Log>(expr)) {
            auto args = expr.get_args();
            auto arg = convert_expression_to_function(Expression(args[0]), x, y, z);
            if (args.size() == 2) { // base log

                auto expr_arg1 = Expression(args[0]);
                auto expr_arg2 = Expression(args[1]);
                if(is_number(expr_arg1)){
                    if(is_number(expr_arg2)){
                        double constant_value = std::log(eval_double(expr_arg1))/std::log(eval_double(expr_arg2));
                        return [constant_value](double,double,double) { return constant_value; };
                    }
                    else if(is_a<const SymEngine::Symbol>(expr_arg2)){
                        double constant_value = std::log(eval_double(expr_arg1));
                        if(eq(expr_arg2,*x))
                            return [constant_value](double x_value,double y_value,double z_value) {
                                return constant_value/std::log(x_value);
                            };
                        if(eq(expr_arg2,*y))
                            return [constant_value](double x_value,double y_value,double z_value) {
                                return constant_value/std::log(y_value);
                            };
                        else
                            return [constant_value](double x_value,double y_value,double z_value) {
                                return constant_value/std::log(z_value);
                            };
                    }
                    else{
                        double constant_value = eval_double(expr_arg1);
                        auto f = convert_expression_to_function(expr_arg2, x,y,z);
                        return [f,constant_value](double x_value,double y_value,double z_value) {
                            return std::log(constant_value)/std::log(f(x_value,y_value,z_value));
                        };
                    }
                }
                else if (is_a<const SymEngine::Symbol>(expr_arg1)){
                    if(eq(expr_arg1,*x)){
                        if(is_number(expr_arg2)){
                            double constant_value = 1.0/std::log(eval_double(expr_arg2));
                            return [constant_value](double x_value,double y_value,double z_value) { return std::log(x_value)*constant_value; };
                        }
                        else if(is_a<const SymEngine::Symbol>(expr_arg2)){
                            if(eq(expr_arg2,*x))
                                return [](double x_value,double y_value,double z_value) {
                                    return 1.0;
                                };
                            if(eq(expr_arg2,*y))
                                return [](double x_value,double y_value,double z_value) {
                                    return std::log(x_value)/std::log(y_value);
                                };
                            else
                                return [](double x_value,double y_value,double z_value) {
                                    return std::log(x_value)/std::log(z_value);
                                };
                        }
                        else{
                            auto f = convert_expression_to_function(expr_arg2, x,y,z);
                            return [f](double x_value,double y_value,double z_value) {
                                return std::log(x_value)/std::log(f(x_value,y_value,z_value));
                            };
                        }
                    }
                    if(eq(expr_arg1,*y)){
                        if(is_number(expr_arg2)){
                            double constant_value = 1.0/std::log(eval_double(expr_arg2));
                            return [constant_value](double x_value,double y_value,double z_value) { return std::log(y_value)*constant_value; };
                        }
                        else if(is_a<const SymEngine::Symbol>(expr_arg2)){
                            if(eq(expr_arg2,*x))
                                return [](double x_value,double y_value,double z_value) {
                                    return std::log(y_value)/std::log(x_value);
                                };
                            if(eq(expr_arg2,*y))
                                return [](double x_value,double y_value,double z_value) {
                                    return 1.0;
                                };
                            else
                                return [](double x_value,double y_value,double z_value) {
                                    return std::log(y_value)/std::log(z_value);
                                };
                        }
                        else{
                            auto f = convert_expression_to_function(expr_arg2, x,y,z);
                            return [f](double x_value,double y_value,double z_value) {
                                return std::log(y_value)/std::log(f(x_value,y_value,z_value));
                            };
                        }

                    }
                    else{//z
                        if(is_number(expr_arg2)){
                            double constant_value = 1.0/std::log(eval_double(expr_arg2));
                            return [constant_value](double x_value,double y_value,double z_value) { return std::log(z_value)*constant_value; };
                        }
                        else if(is_a<const SymEngine::Symbol>(expr_arg2)){
                            if(eq(expr_arg2,*x))
                                return [](double x_value,double y_value,double z_value) {
                                    return std::log(z_value)/std::log(x_value);
                                };
                            if(eq(expr_arg2,*y))
                                return [](double x_value,double y_value,double z_value) {
                                    return std::log(z_value)/std::log(y_value);
                                };
                            else
                                return [](double x_value,double y_value,double z_value) {
                                    return 1.0;
                                };
                        }
                        else{
                            auto f = convert_expression_to_function(expr_arg2, x,y,z);
                            return [f](double x_value,double y_value,double z_value) {
                                return std::log(z_value)/std::log(f(x_value,y_value,z_value));
                            };
                        }

                    }
                }
                else{
                    auto f = convert_expression_to_function(expr_arg1, x,y,z);
                    if(is_number(expr_arg2)){

                        double constant_value = 1.0/std::log(eval_double(expr_arg2));
                        return [f,constant_value](double x_value,double y_value,double z_value) { return std::log(f(x_value,y_value,z_value))*constant_value; };
                    }
                    else if(is_a<const SymEngine::Symbol>(expr_arg2)){
                        if(eq(expr_arg2,*x))
                            return [f](double x_value,double y_value,double z_value) {
                                return std::log(f(x_value,y_value,z_value))/std::log(x_value);
                            };
                        if(eq(expr_arg2,*y))
                            return [f](double x_value,double y_value,double z_value) {
                                return std::log(f(x_value,y_value,z_value))/std::log(y_value);
                            };
                        else
                            return [f](double x_value,double y_value,double z_value) {
                                return std::log(f(x_value,y_value,z_value))/std::log(z_value);
                            };
                    }
                    else{
                        auto f2 = convert_expression_to_function(expr_arg2, x,y,z);
                        return [f,f2](double x_value,double y_value,double z_value) {
                            return std::log(f(x_value,y_value,z_value))/std::log(f2(x_value,y_value,z_value));
                        };
                    }

                }

            } else { // ln

                auto expr_arg = Expression(args[0]);

                if(is_number(expr_arg)){
                    double constant_value = std::log(eval_double(expr_arg));
                    return [constant_value](double,double,double) { return constant_value; };
                }
                else if (is_a<const SymEngine::Symbol>(expr_arg)) {
                    if(eq(expr_arg,*x))
                        return [](double x_value,double y_value,double z_value) { return std::log(x_value); };
                    if(eq(expr_arg,*y))
                        return [](double x_value,double y_value,double z_value) { return std::log(y_value); };
                    else
                        return [](double x_value,double y_value,double z_value) { return std::log(z_value); };
                }
                auto arg = convert_expression_to_function(expr_arg, x,y,z);
                return [arg](double x_value, double y_value, double z_value) {
                    return std::log(arg(x_value, y_value, z_value));
                };
            }
        }

        throw std::runtime_error("Unsupported expression type");
    }
    /*
    std::set<double> find_singularities(const Expression& expr, const RCP<const Symbol> &x) {
        std::set<double> singularities;

        if (is_a<Mul>(expr)) {
            auto mul_expr = rcp_static_cast<const Mul>(expr.get_basic());
            for (auto arg : mul_expr->get_args()) {
                if (is_a<Pow>(*arg)) {
                    auto pow_expr = rcp_static_cast<const Pow>(arg);
                    if (pow_expr->get_exp()->__eq__(*SymEngine::minus_one)) {
                        auto denominator = pow_expr->get_base();
                        auto solutions = solve(denominator, x);
                        if (is_a<FiniteSet>(*solutions)) {
                            auto finite_set_casted = rcp_static_cast<const FiniteSet>(solutions);
                            auto elements = finite_set_casted->get_container();
                            for (auto sol : elements) {

                                if (is_a<const Integer>(*sol))
                                    singularities.insert(SymEngine::rcp_static_cast<const SymEngine::Integer>(sol)->as_int());
                                else
                                    singularities.insert(SymEngine::rcp_static_cast<const SymEngine::RealDouble>(sol)->as_double());
                            }
                        }
                    }
                } else {
                    auto sub_singularities = find_singularities(Expression(arg), x);
                    singularities.insert(sub_singularities.begin(), sub_singularities.end());
                }
            }
        }
        
        if (is_a<Log>(expr)) {
            auto log_expr = rcp_static_cast<const Log>(expr.get_basic());
            auto log_argument = log_expr->get_args()[0];
            auto solutions = solve(log_argument, x);
            if (is_a<FiniteSet>(*solutions)) {
                auto finite_set_casted = rcp_static_cast<const FiniteSet>(solutions);
                auto elements = finite_set_casted->get_container();
                for (auto sol : elements) {
                    if (is_a<const Integer>(*sol))
                        singularities.insert(SymEngine::rcp_static_cast<const SymEngine::Integer>(sol)->as_int());
                    else
                        singularities.insert(SymEngine::rcp_static_cast<const SymEngine::RealDouble>(sol)->as_double());
                }
            }
        }

       
        if (is_a<Pow>(expr)) {
            auto pow_expr = rcp_static_cast<const Pow>(expr.get_basic());
            auto exponent = pow_expr->get_exp();
            //std::cout<<*exponent<<std::endl;

            if (is_a<const RealDouble>(*exponent) or is_a<const Integer>(*exponent)) {
                double exp_val;
                if (is_a<const Integer>(*exponent))
                    exp_val=SymEngine::rcp_static_cast<const SymEngine::Integer>(exponent)->as_int();
                else
                    exp_val=SymEngine::rcp_static_cast<const SymEngine::RealDouble>(exponent)->as_double();
        
                if (exp_val < 0 || (exp_val < 2 && std::floor(exp_val) != exp_val)) {
                    auto base = pow_expr->get_base();
                    auto solutions = solve(base, x);
                    if (is_a<FiniteSet>(*solutions)) {
                        auto finite_set_casted = rcp_static_cast<const FiniteSet>(solutions);
                        auto elements = finite_set_casted->get_container();
                        for (auto sol : elements) {
                            if (is_a<const Integer>(*sol))
                                singularities.insert(SymEngine::rcp_static_cast<const SymEngine::Integer>(sol)->as_int());
                            else
                                singularities.insert(SymEngine::rcp_static_cast<const SymEngine::RealDouble>(sol)->as_double());
                        }
                    }
                }
            }
        }

        if (is_a<SymEngine::Add>(expr)) {
            auto add_expr = rcp_static_cast<const SymEngine::Add>(expr.get_basic());
            for (auto arg : add_expr->get_args()) {
                auto sub_singularities = find_singularities(Expression(arg), x);
                singularities.insert(sub_singularities.begin(), sub_singularities.end());
            }
        }

        if (is_a<SymEngine::Mul>(expr) && !is_a<Pow>(expr)) {
            auto mul_expr = rcp_static_cast<const Mul>(expr.get_basic());
            for (auto arg : mul_expr->get_args()) {
                auto sub_singularities = find_singularities(Expression(arg), x);
                singularities.insert(sub_singularities.begin(), sub_singularities.end());
            }
        }

        return singularities;
    }*/

/*
template<class T>
inline double evaluate(const Expression & func, T val) {
            
    return (double)func.subs({{x,real_double(val)}}); 

} 
*/
#endif

