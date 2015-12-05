#ifndef ENVIRONMENT_KALMAN_FILTER_MEASUREMENT_MODEL_H
#define ENVIRONMENT_KALMAN_FILTER_MEASUREMENT_MODEL_H

#include <kalman/MeasurementModel.hpp>
#include "filter/system_model.h"

/**
 * @brief Measurement vector measuring the acceleration in x- and y-direction and the turn rate
 *
 * @param T Numeric scalar type
 */
template<typename T>
class Measurement : public Kalman::Vector<T, 12>
{
public:
    // time step from the last state to this state
    T dt = 0.01; //seconds

    KALMAN_VECTOR(Measurement, T, 12)

    //! angular velocity IMU
    static constexpr size_t OMEGA = 0;

    //! Delta S hall
    static constexpr size_t DS_HALL = 1;

    //! delta x mouse sensor A
    static constexpr size_t DX_MOUSEA = 2;

    //! delta y mouse sensor A
    static constexpr size_t DY_MOUSEA = 3;

    //! delta x mouse sensor B
    static constexpr size_t DX_MOUSEB = 4;

    //! delta y mouse sensor B
    static constexpr size_t DY_MOUSEB = 5;

    //! Y0 street model kalman tobi
    static constexpr size_t Y0 = 6;

    //! psi0 = angle between road and car x axis: see tobis street environment kalman
    static constexpr size_t PSI0 = 7;

    //! Kappa 1: tobis kalman
    static constexpr size_t KAPPA1 = 8;
    //! length of line to point where kappa 1 is
    const T l1 = 0.2;

    //! Kappa 2: tobis kalman
    static constexpr size_t KAPPA2 = 9;
    //! length of line to point where kappa 2 is
    const T l2 = 2*0.2;

    //! Kappa 3: tobis kalman
    static constexpr size_t KAPPA3 = 10;
    //! length of line to point where kappa 3 is
    const T l3 = 3*0.2;

    //! Kappa 4: tobis kalman
    static constexpr size_t KAPPA4 = 11;
    //! length of line to point where kappa 4 is
    const T l4 = 4*0.2;

    //positions of the mouse sensors A and B
    const T xA = 0.13; //meters
    const T yA = 0.05;
    const T xB = -0.13;
    const T yB = -0.05;


    T w()                 const { return (*this)[ OMEGA ]; }
    T ds_hall()           const { return (*this)[ DS_HALL ]; }
    T dx_mouse_A()        const { return (*this)[ DX_MOUSEA ]; }
    T dy_mouse_A()        const { return (*this)[ DY_MOUSEA ]; }
    T dx_mouse_B()        const { return (*this)[ DX_MOUSEB ]; }
    T dy_mouse_B()        const { return (*this)[ DY_MOUSEB ]; }
    T y0()                const { return (*this)[ Y0 ]; }
    T psi()               const { return (*this)[ PSI0 ]; }
    T kappa1()            const { return (*this)[ KAPPA1 ]; }
    T kappa2()            const { return (*this)[ KAPPA2 ]; }
    T kappa3()            const { return (*this)[ KAPPA3 ]; }
    T kappa4()            const { return (*this)[ KAPPA4 ]; }

    /*T len1()            const { return this->l1; }
    T len2()            const { return this->l2; }
    T len3()            const { return this->l3; }
    T len4()            const { return this->l4; }*/


    T& w()                  { return (*this)[ OMEGA ]; }
    T& ds_hall()            { return (*this)[ DS_HALL ]; }
    T& dx_mouse_A()         { return (*this)[ DX_MOUSEA ]; }
    T& dy_mouse_A()         { return (*this)[ DY_MOUSEA ]; }
    T& dx_mouse_B()         { return (*this)[ DX_MOUSEB ]; }
    T& dy_mouse_B()         { return (*this)[ DY_MOUSEB ]; }
    T& y0()                 { return (*this)[ Y0 ]; }
    T& psi()                { return (*this)[ PSI0 ]; }
    T& kappa1()             { return (*this)[ KAPPA1 ]; }
    T& kappa2()             { return (*this)[ KAPPA2 ]; }
    T& kappa3()             { return (*this)[ KAPPA3 ]; }
    T& kappa4()             { return (*this)[ KAPPA4 ]; }



};

/**
 * @brief Measurement model
 *
 *
 * @param T Numeric scalar type
 * @param CovarianceBase Class template to determine the covariance representation
 *                       (as covariance matrix (StandardBase) or as lower-triangular
 *                       coveriace square root (SquareRootBase))
 */
template<typename T, template<class> class CovarianceBase = Kalman::StandardBase>
class MeasurementModel : public Kalman::MeasurementModel<State<T>, Measurement<T>, CovarianceBase>
{
public:
    //! State type shortcut definition
    typedef State<T> S;

    //! Measurement type shortcut definition
    typedef Measurement<T> M;

    /**
     * @brief Definition of (possibly non-linear) measurement function
     *
     * This function maps the system state to the measurement that is expected
     * to be received from the sensor assuming the system is currently in the
     * estimated state.
     *
     * @param [in] x The system state in current time-step
     * @returns The (predicted) sensor measurement for the system state
     */
    M h(const S& x) const
    {
        M measurement;


        measurement.w() = x.w();
        measurement.ds_hall() = x.u() - 0.5*x.a()*pow(dt, 2);
        //very simple formulas. maybe better use more advanced CTRA model
        measurement.dx_mouse_A() = (cos(x.alpha())*x.u() - yA*x.w()) * dt;
        measurement.dy_mouse_A() = (sin(x.alpha())*x.u() + xA*x.w()) * dt;
        measurement.dx_mouse_B() = (cos(x.alpha())*x.u() - yB*x.w()) * dt;
        measurement.dy_mouse_B() = (sin(x.alpha())*x.u() + xB*x.w()) * dt;
        T cosAngleDiff = cos(x.phi()-x.theta());
        if (abs(cosAngleDiff) < 0.01)
        {
            //throw error, bad case
            measurement.y0 = 0;
        }
        else
        {
            measurement.y0 = x.d() / cosAngleDiff;
        }
        measurement.psi() = x.phi() - x.theta(); //angle between street and car x axis
        T tanAngleDiff = tan(x.phi() - x.theta());
        if (abs(tanAngleDiff) > 0.2)
        {
            //throw error. something went wrong
            measurement.kappa1 = l1*x.c1() + x.c0();
            measurement.kappa2 = l2*x.c2() + x.c0();
            measurement.kappa3 = l3*x.c3() + x.c0();
            measurement.kappa4 = l4*x.c4() + x.c0();
        }else
        {
            measurement.kappa1 = (tanAngleDiff*x.d() + l1*x.c1()) + x.c0();
            measurement.kappa2 = (tanAngleDiff*x.d() + l2*x.c2()) + x.c0();
            measurement.kappa3 = (tanAngleDiff*x.d() + l3*x.c3()) + x.c0();
            measurement.kappa4 = (tanAngleDiff*x.d() + l4*x.c4()) + x.c0();
        }
        return measurement;
    }

};


#endif
