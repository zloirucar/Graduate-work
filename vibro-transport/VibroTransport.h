#ifndef _VIBRO_TRANSPORT_H_
#define _VIBRO_TRANSPORT_H_

#include "ode_num_int/OdeRhs.h"

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
			real_type N;
			real_type R,R0;
            switch (m_discreteState) {
                case F1:
					dst[0] = x[2];
					dst[1] = x[3];
					dst[2] = -g*sin(m_alf) + sqr(m_ow)*A*sin(m_ow*time);
					dst[3] = -g*cos(m_alf) + sqr(m_ow)*B*sin(m_ow*time+m_Eps);
                    break;
				case F:
					dst[0] = x[2];
					dst[1] = x[3];
					dst[2] = -g*sin(m_alf) + sqr(m_ow)*A*sin(m_ow*time);
					dst[3] = -g*cos(m_alf) + sqr(m_ow)*B*sin(m_ow*time + m_Eps);
					break;
                case S:
                    {
                    computeReactions_S(N, R, time, x);
					dst[0] = x[2];
					dst[1] = x[3];
					dst[2] = -g*sin(m_alf) + R/m_mass + sqr(m_ow)*A*sin(m_ow*time);
					dst[3] = 0;
				    break;
                    }
                case K:
					computeReactions_K(N, R,R0, time, x);
					dst[0] = x[2];
					dst[1] = x[3];
					dst[2] = 0;
					dst[3] = 0;

                    break;
                    
			}
            this->odeRhsPostObservers( time, x, dst, this );
		};

		

        virtual unsigned int zeroFuncCount() const
            {
            return 3;
            }

        virtual void zeroFunctions( V& dst, real_type time, const V& x ) const
            {
			real_type N, R, R0;
			dst.resize(3);
			switch (m_discreteState) {
			case F:
				dst[0] = x[1];
				dst[1] = 1;
				dst[2] = 1;
				break;
			case F1:
				dst[0] = 1;
				dst[1] = x[1];
				dst[2] = min_speed - x[3];
				break;
			case S:
				dst[0] = N;
				dst[1] = x[2];
				dst[2] = 1;
				break;
			case K:
				computeReactions_K(N, R, R0, time, x);
				dst[0] = N;
				dst[1] = fabs(R) - fabs(R0);
				dst[2] = 1;
				break;
			};
			
			};


			virtual std::vector<unsigned int> zeroFuncFlags() const {
				return std::vector<unsigned int>{
					OdeRhs<VD>::Discontinuous | OdeRhs<VD>::PlusMinus,
					OdeRhs<VD>::Discontinuous | OdeRhs<VD>::BothDirections,
					OdeRhs<VD>::Discontinuous | OdeRhs<VD>::PlusMinus

				};
			}


			virtual void switchPhaseState(int* transitions, real_type time, V& x) {
				real_type N, R, R0;
				switch (m_discreteState) {
				case F1:
					if (transitions[0] < 0) {
						transitions[0] = 1;
						m_discreteState = F;
					};

					if (transitions[1] < 0) {
						transitions[1] = x[2] > 0 ?   1 :   x[2] < 0 ?   -1 :   0;
						x[1] = 0;
						x[3] = 0;
						m_discreteState = S;
					};
					break;
				case F:
					ASSERT(transitions[0] != 0);
					if (x[3] < 0) {
						x_method(x); // удар по касательной
						x[3] *= -y_method(); //  нормальный удар
						x[1] = 0;
						transitions[0] = 1;
					}
					else {
						transitions[0] = 1;
						m_discreteState = S;
					};
					if (fabs(x[3]) < min_speed) {
						x[1] = 0;
						x[3] = 0;
						m_discreteState = S;
					}
					break;
				case S:
						if (transitions[1] < 0) {
						x[2] = 0;
						x[1] = 0;
						m_discreteState = K;
						transitions[1] = 1;
					};
					if (transitions[0] < 0) {
						x[1] = 0;
						m_discreteState = F1;
						transitions[0] = 1;
					};
					break;
				case K:
					if (transitions[1] < 0) {
						x[3] = 0;
						m_discreteState = S;
						transitions[1] = 1;
					};
					if (transitions[0] < 0) {
						x[1] = 0;
						x[3] = 0;
						m_discreteState = F1;
						transitions[0] = 1;
					}

					break;
				default:
					ASSERT(false); // Unreachable
				}

			}

			virtual std::string describeZeroFunction(unsigned int /*index*/) const {
				return std::string();
			}

			enum DiscreteState { F, S, K, F1 };
			DiscreteState discreteState() const {
				return m_discreteState;
			}

			void setDiscreteState(DiscreteState discreteState) {
				m_discreteState = discreteState;
			}

			void computeDiscreteState(double time, const V& state)
			{
				real_type N, R, R0;
				computeReactions_K(N, R, R0, time, state);
				if (state[1] == 0 && state[2] == 0) {
					if (N < 0) {
						m_discreteState = F1;
					}
					else m_discreteState = K;
				} 

				else if (state[1] == 0 && state[2] != 0) {
					if (N < 0) {
						m_discreteState = F1;
					} else m_discreteState = S; 
				}
		
				else m_discreteState = F;

				}
			



	private:

		DiscreteState m_discreteState;

		//Объявление переменных
		const real_type g = 9.8;
		const real_type m_alf = 0.0785;// угол наклона поверхности
		const real_type m_ow = 100; // частота колебаний
		const real_type m_f = 0.4; // трение тела о лоток
		const real_type m_mass = 1; //масса т.т.
		const real_type A = 0.0003;//Амплитуда по х
		const real_type B = 0.0009;//Амлитуда по y
		const real_type m_Eps = -0.5;//смещение фазы
		const real_type min_speed = 0.05; // минимальная скорость для отскока
		
		
		template< class T >
		static T sqr(T x) {
			return x*x;
		};


		void computeReactions_S(real_type& N, real_type& R, real_type time, const V& x) const
		{
			N = m_mass*g*cos(m_alf) - m_mass*sqr(m_ow)*B*sin(m_ow*time + m_Eps);
			if (x[2] > 0) {
				R = -m_f*N;
			}
			else {
				R = m_f*N;
			}
		}

        void computeReactions_K(real_type& N, real_type& R,real_type& R0, real_type time, const V& x) const
            {
			N = m_mass*g*cos(m_alf) - m_mass*sqr(m_ow)*B*sin(m_ow*time + m_Eps);
			
			if (x[2] > 0) {
				R = -m_f*N/0.7;
			}
			else {
				R = m_f*N/0.7;
			};

			R0 = m_mass*g*sin(m_alf) - m_mass*sqr(m_ow)*A*sin(m_ow*time);
		};
   

	real_type y_method() {
		real_type R=0.7;
		// Здесь будет метод вычисления коэфициента восстановления
		// TODO
		return R;

	};

	void x_method(V& x)  {
		
		real_type lam;
		lam = m_f; //Иногда принмают как коэф трения скольжения
		
			x[2] *= (1 - lam);
		
		};
		};




		

#endif // _VIBRO_TRANSPORT_H_
