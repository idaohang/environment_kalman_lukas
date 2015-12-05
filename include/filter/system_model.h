#ifndef ENVIRONMENT_KALMAN_FILTER_SYSTEM_MODEL_H
#define ENVIRONMENT_KALMAN_FILTER_SYSTEM_MODEL_H

#include <kalman/SystemModel.hpp>


/**
 * @brief System state vector-type for CTRA motion model
 *
 * @param T Numeric scalar type
 */
template<typename T>
class State : public Kalman::Vector<T, 11>
{
public:
    KALMAN_VECTOR(State, T, 11)

    // x = [u_0,a,s_0,l_0,phi_0,c_0,c_1,theta_0,omega,alpha_0,alphaDot]

    //! Velocity (scalar direction given by angle)
    static constexpr size_t U = 0;
    //! Acceleration (same direction as the velocity)
    static constexpr size_t A = 1;
    //! Frenet space: longitudinal distance from the start
    static constexpr size_t S = 2;
    //! Frenet space: lateral distance from the center line
    static constexpr size_t D = 3;
    //! Angle between the tangent of the center line and some fixed global COSY
    static constexpr size_t PHI = 4;
    //! Curvature of the center line
    static constexpr size_t C0 = 5;
    //! Change of curvature of the center line
    static constexpr size_t C1 = 6;
    //! Angle between the x axis of the car and the same fixed global COSY
    static constexpr size_t THETA = 7;
    //! Angular velocity of the car
    static constexpr size_t OMEGA = 8;
    //! Angle between the car x axis and the vector of the velocity of the car
    static constexpr size_t ALPHA = 9;
    //! Angular velocity of the vector above (relative to the car)
    static constexpr size_t ALPHADOT = 10;

    T u()       const { return (*this)[ U ]; }
    T a()       const { return (*this)[ A ]; }
    T s()       const { return (*this)[ S ]; }
    T d()       const { return (*this)[ D ]; }
    T phi()     const { return (*this)[ PHI ]; }
    T c0()      const { return (*this)[ C0 ]; }
    T c1()      const { return (*this)[ C1 ]; }
    T theta()   const { return (*this)[ THETA ]; }
    T w()       const { return (*this)[ OMEGA ]; }
    T alpha()       const { return (*this)[ ALPHA ]; }
    T alphadot()      const { return (*this)[ ALPHADOT ]; }

    T &u()        { return (*this)[ U ]; }
    T &a()        { return (*this)[ A ]; }
    T &s()        { return (*this)[ S ]; }
    T &d()        { return (*this)[ D ]; }
    T &phi()      { return (*this)[ PHI ]; }
    T &c0()       { return (*this)[ C0 ]; }
    T &c1()       { return (*this)[ C1 ]; }
    T &theta()    { return (*this)[ THETA ]; }
    T &w()        { return (*this)[ OMEGA ]; }
    T &alpha()        { return (*this)[ ALPHA ]; }
    T &alphadot()       { return (*this)[ ALPHADOT ]; }
};

/**
 * @brief System control-input vector-type for CTRA motion model
 *
 *
 * @param T Numeric scalar type
 */
template<typename T>
class Control : public Kalman::Vector<T, 1>
{
public:
    KALMAN_VECTOR(Control, T, 1)

    //! time since the filter was last called
    static constexpr size_t DT = 0;


    T dt()      const { return (*this)[ DT ]; }


    T& dt()     { return (*this)[ DT ]; }
};

/**
 * @brief System model CALC (Constant Acceleration Linear Curvature)
 *
 * This is the system model defining how our car moves from one
 * time-step to the next, i.e. how the system state evolves over time.
 *
 * @param T Numeric scalar type
 * @param CovarianceBase Class template to determine the covariance representation
 *                       (as covariance matrix (StandardBase) or as lower-triangular
 *                       coveriace square root (SquareRootBase))
 */
template<typename T, template<class> class CovarianceBase = Kalman::StandardBase>
class SystemModel : public Kalman::SystemModel<State<T>, Control<T>, CovarianceBase>
{
public:
    //! State type shortcut definition
    typedef State<T> S;

    //! Control type shortcut definition
    typedef Control<T> C;

    /**
     * @brief Definition of (non-linear) state transition function
     *
     * This function defines how the system state is propagated through time,
     * i.e. it defines in which state \f$\hat{x}_{k+1}\f$ is system is expected to
     * be in time-step \f$k+1\f$ given the current state \f$x_k\f$ in step \f$k\f$ and
     * the system control input \f$u\f$.
     *
     * @param [in] x The system state in current time-step
     * @param [in] u The control vector input
     * @returns The (predicted) system state in the next time-step
     */
    S f(const S& x, const C& u) const
    {
        //! Predicted state vector after transition
        S x_;

        //
        T dt = u.dt();
        T ds_preCalc = ds(x, u);
        T dd_preCalc = dd(x, u);

        x_.u() = x.u() + x.a()*dt; //new velocity
        x_.a() = x.a(); //constant acceleration
        x_.s() = x.s() + ds_preCalc; //s_new = s_old + ds
        x_.d() = x.d() + dd_preCalc;
        x_.phi() = x.phi() + x.c0*ds_preCalc + 0.5*x.c1()*pow(ds_preCalc, 2); //2 times integrated
        x_.c0() = x.c0() + x.c1() * ds_preCalc; //linear curvature
        x_.c1() = x.c1(); //linear curvature
        x_.theta() = x.theta() + x.w() * dt;
        x_.w() = x.w();
        x_.alpha() = x.alpha() + x.alphadot()*dt;
        x_.alphadot() = x.alphadot();


        // Return transitioned state vector
        return x_;
    }

private:

    /*
     * @brief: calculates the difference in s direction between the current and the next state
     */
    T ds(const S& x, const C& u) const
    {
        T g0 = x.theta() + x.alpha() - x.phi(); //theta_0 + alpha_0 - phi_0
        T O = x.w + x.alphadot() - x.u()*x.c0(); //omega + alphaDot - u_0*c_0
        T gT = g0 + O*u.dt(); //gamma(t+dt)
        T u0 = x.u();
        T uT = x.u() + x.a()*u.dt(); //u_0 + a*T

        T sinG0 = sin(g0);
        T cosG0 = cos(g0);

        T sinGT = sin(gT);
        T cosGT = cos(gT);

        T O2 = O^2;


        T ds = 0;

        if(abs(x.c0()*x.d()) > 0.4)
        {
            //throw error
            return ds;
        }

        // real computation
        // cases for o approx. 0
        if (abs(O) < 0.001)
        {
            //limit case
            ds = (1/( 1+x.c0()*x.d() ))  * cosG0 * (x.u()*u.dt() + 0.5*x.a()*pow(u.dt(),2));
        }
        else
        {
            ds = (1/( 1+x.c0()*x.d() )) * (1/O2) * (O*(uT*sinGT - u0*sinG0)  + x.a()*(cosGT - cosG0));
        }

        return ds;
    }

    /*
     * @brief: calculates the difference in s direction between the current and the next state
     */
    T dd(const S& x, const C& u) const
    {
        T g0 = x.theta() + x.alpha() - x.phi(); //theta_0 + alpha_0 - phi_0
        T O = x.w + x.alphadot() - x.u()*x.c0(); //omega + alphaDot - u_0*c_0
        T gT = g0 + O*u.dt(); //gamma(t+dt)
        T u0 = x.u();
        T uT = x.u() + x.a()*u.dt(); //u_0 + a*T

        T sinG0 = sin(g0);
        T cosG0 = cos(g0);

        T sinGT = sin(gT);
        T cosGT = cos(gT);

        T O2 = O^2;


        T dd = 0;


        // real computation
        // cases for o approx. 0
        if (abs(O) < 0.001)
        {
            //limit case
            dd = sinG0 * (x.u()*u.dt() + 0.5*x.a()*pow(u.dt(),2));
        }
        else
        {
            dd = (1/O2) * (O*(u0*cosG0 - uT*cosGT)  +  x.a()*(sinGT-sinG0));
        }

        return dd;

    }


};

#endif
