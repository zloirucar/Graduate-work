#ifndef _VIBRO_TRANSPORT_H_
#define _VIBRO_TRANSPORT_H_

#include "ode_num_int/OdeRhs.h"
#include "ode_num_int/util/m_const.h"

using namespace ctm;
using namespace math;

template< class VD >
class VibroTransport :
    public OdeRhs< VD >,
    public FactoryMixin< VibroTransport<VD>, OdeRhs<VD> >
    {
    public:
        typedef VectorTemplate< VD > V;
        typedef typename V::value_type real_type;
        typedef OptionalParameters::Parameters Parameters;

        explicit VibroTransport() : m_discreteState(F)
            {}

        virtual unsigned int secondOrderVarCount() const {
            return 2;
            }

        virtual unsigned int firstOrderVarCount() const {
            return 0;
            }

        virtual void rhs( V& dst, real_type time, const V& x ) const
            {
            this->odeRhsPreObservers( time, x, this );
			dst.resize(4);
            switch (m_discreteState) {
                case F:
					dst[0] = x[2];
					dst[1] = x[3];
					dst[2] = -g*sin(m_alf) + sqr(m_ow)*A*sin(m_ow*time);
					dst[3] = -g*cos(m_alf) + sqr(m_ow)*B*sin(m_ow*time+m_Eps);
                    break;
                case SL:
                case SR:
                    {
                    auto N = normalReaction(time, x);
                    auto R = -sgn(x[2]) * slidingFrictionForceNorm(N);
                    dst[0] = x[2];
					dst[1] = x[3];
					dst[2] = -g*sin(m_alf) + R/m_mass + sqr(m_ow)*A*sin(m_ow*time);
                    dst[3] = 0;
				    break;
                    }
                case K:
					dst[0] = x[2];
					dst[1] = x[3];
					dst[2] = 0;
					dst[3] = 0;
                    break;
                }
            this->odeRhsPostObservers( time, x, dst, this );
            }
		

        virtual unsigned int zeroFuncCount() const
            {
            return 2;
            }

        virtual void zeroFunctions( V& dst, real_type time, const V& x ) const
            {
			dst.resize(2);
			switch (m_discreteState) {
                case F:
                    dst[0] = x[1];
                    dst[1] = 1;
                    break;
                case SL:
                    dst[0] = normalReaction(time, x);
                    dst[1] = -x[2];
                    break;
                case SR:
                    dst[0] = normalReaction(time, x);
                    dst[1] = x[2];
                    break;
                case K: {
                    auto N = normalReaction(time, x);
                    dst[0] = N;
                    dst[1] = staticFrictionForceNorm(N) - fabs(stickingFrictionForce(time));
                    break;
                    }
                }
            }


			virtual std::vector<unsigned int> zeroFuncFlags() const {
                return std::vector<unsigned int>(2, OdeRhs<VD>::PlusMinus | OdeRhs<VD>::RecomputeAfterSwitch);
            }


            virtual void switchPhaseState(const int* transitions, real_type time, V& x) {
				switch (m_discreteState) {
                case F: {
					ASSERT(transitions[0] != 0);
                    auto fallingDown = false;
					if (x[3] < 0) {
						x_method(x); // удар по касательной
						x[3] *= -y_method(); //  нормальный удар
						x[1] = 0;
                        fallingDown = true;
                    }
                    if ((!fallingDown   ||   fabs(x[3]) < min_speed)   &&   normalReaction(time, x) > 0) {
						x[1] = 0;
						x[3] = 0;
                        m_discreteState = x[2] >= 0? SR: SL;
					}
					break;
                    }
                case SL:
                case SR:
                    if (transitions[0] != 0)
                        m_discreteState = F;
                    else {
                        ASSERT(transitions[1] != 0);
                        auto N = normalReaction(time, x);
                        if (fabs(stickingFrictionForce(time)) < staticFrictionForceNorm(N))
                            m_discreteState = K;
                        else
                            m_discreteState = m_discreteState == SL? SR: SL;
                    }
					break;
				case K:
                    if (transitions[0] != 0)
                        m_discreteState = F;
                    else {
                        ASSERT(transitions[1] != 0);
                        m_discreteState = stickingFrictionForce(time) >= 0? SL: SR;
                    }
                    break;
				default:
					ASSERT(false); // Unreachable
				}
			}

			virtual std::string describeZeroFunction(unsigned int /*index*/) const {
				return std::string();
			}

            enum DiscreteState { F, SL, SR, K };
			DiscreteState discreteState() const {
				return m_discreteState;
			}

			void setDiscreteState(DiscreteState discreteState) {
				m_discreteState = discreteState;
			}

			void computeDiscreteState(double time, const V& state)
			{
                if (state[1] > 0)
                    m_discreteState = F;
                else if (state[1] < 0)
                    throw cxx::exception("Invalid state: x[1] must be nonnegative");
                else {
                    auto N = normalReaction(time, state);
                    if (N <= 0)
                        m_discreteState = F;
                    else {
                        if (state[2] == 0) {
                            auto R = slidingFrictionForceNorm(N);
                            auto R0 = stickingFrictionForce(time);
                            m_discreteState = fabs(R0) < R ?   K :   R0 >= 0? SR: SL;
                            }
                        else
                            m_discreteState = state[2] >= 0? SR: SL;
                        }
                    }
            }

        real_type alpha() const {
            return m_alf;
            }
        void setAlpha(real_type alpha) {
            m_alf = alpha;
            }

        real_type amp_A() const {
            return A;
            }
        void setA(real_type amp_A) {
            A = amp_A;
            }

        real_type ow() const {
            return m_ow;
            }
        void setow(real_type ow) {
            m_ow = ow;
            }

	private:

		DiscreteState m_discreteState;

		//Объявление переменных
		const real_type g = 9.81;
        real_type m_alf = M_PI / 6;        // угол наклона поверхности
        real_type m_ow = 100;           // частота колебаний
		const real_type m_f = 0.7;             // трение тела о лоток
		const real_type m_f0 = 1;             // трение покоя тела о лоток
		const real_type m_mass = 1;                 // масса т.т.
        real_type A = 0;      // Амплитуда по х
        const real_type B = 0.0009;     // Амлитуда по y
        const real_type m_Eps = 0.35*M_PI/180;         // смещение фазы
		const real_type min_speed = 1e-4;           // минимальная скорость для отскока
		
		
		template< class T >
        static T sqr(T x) {
			return x*x;
        }
		


        template< class T >
        static T sgn(T x) {
            return x > 0? 1: x<0? -1: 0;
        }


        real_type normalReaction(real_type time, const V& x) const
        {
            return m_mass*g*cos(m_alf) - m_mass*sqr(m_ow)*B*sin(m_ow*time + m_Eps);
        }

        real_type slidingFrictionForceNorm(real_type normalReaction) const {
            return m_f * fabs(normalReaction);
        }

		real_type staticFrictionForceNorm(real_type normalReaction) const {
			return m_f0 * fabs(normalReaction);
		}

        real_type stickingFrictionForce(real_type time) const {
            return m_mass*g*sin(m_alf) - m_mass*sqr(m_ow)*A*sin(m_ow*time);
        }

        real_type y_method() const {
            real_type R=0.3;
            // Здесь будет метод вычисления коэфициента восстановления
            // TODO
            return R;

        }

        void x_method(V& x) const {
            real_type lam;
            lam = m_f; //Иногда принмают как коэф трения скольжения
            if (fabs(x[2]) < fabs(m_f*(y_method()*x[3] - x[3]) / lam)) {
                x[2] *= (1 - lam);
            }
            else {
                x[2] = x[2] - m_f*(y_method()*x[3] - x[3])*Sgn(x[2]);
            }
        }

	};




		

#endif // _VIBRO_TRANSPORT_H_
