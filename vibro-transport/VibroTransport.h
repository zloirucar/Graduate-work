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
            return 1;
            }

        virtual unsigned int firstOrderVarCount() const {
            return 3;
            }

        virtual void rhs( V& dst, real_type time, const V& x ) const
            {
            this->odeRhsPreObservers( time, x, this );
			dst.resize(4);
            switch (m_discreteState) {
                case F:
					dst[0] = x[2];
					dst[1] = x[3];
					dst[2] = -g*sin(m_alf)+sqr(m_ow)*A*sin(m_ow*time);
					dst[3] = -g*cos(m_alf)+sqr(m_ow)*B*sin(m_ow*time+m_Eps);
                    //ASSERT(false); // TODO
                    break;
                case S:
                    {
                    real_type N, R;
                    computeReactions_S(N, R, time, x);
					dst[0] = x[2];
					dst[1] = x[3];
					dst[2] = -g*sin(m_alf) + R/m_mass + sqr(m_ow)*A*sin(m_ow*time);
					dst[3] = -g*cos(m_alf) + N/m_mass + sqr(m_ow)*B*sin(m_ow*time + m_Eps);
					//ASSERT(false); // TODO
                    break;
                    }
                case K:
					real_type N, R,R0;
					computeReactions_K(N, R,R0, time, x);
					dst[0] = x[2];
					dst[1] = x[3];
					dst[2] = -g*sin(m_alf) + R + sqr(m_ow)*A*sin(m_ow*time);
					dst[3] = -g*cos(m_alf) + N +sqr(m_ow)*B*sin(m_ow*time + m_Eps);
                    //ASSERT(false); // TODO
                    break;
                default:
                    ASSERT(false);
                }
            this->odeRhsPostObservers( time, x, dst, this );
		};

		

        virtual unsigned int zeroFuncCount() const
            {
            return 2;
            }

        virtual void zeroFunctions( V& dst, real_type time, const V& x ) const
            {
			real_type N, R,R0;
			dst.resize(2);
			switch (m_discreteState) {
			case F:
				dst[0] = x[1];
				dst[1] = 0;
				break;
			case S:
				
				computeReactions_S(N, R, time, x);
				dst[0] = N;
				dst[1] = x[2];
			case K:
				computeReactions_K(N, R, R0, time, x);
				dst[0] = N;
				dst[1] = R;

			};

			};

		virtual std::vector<unsigned int> zeroFuncFlags() const {
			return std::vector<unsigned int>(1, OdeRhs<VD>::BothDirections);
		}


        virtual void switchPhaseState( const int* /*transitions*/, real_type /*time*/, V& /*x*/ ) {
			ASSERT(false); //TODO
			
            }
        virtual std::string describeZeroFunction( unsigned int /*index*/ ) const {
            return std::string();
            }

        enum DiscreteState { F, S, K };
        DiscreteState discreteState() const {
            return m_discreteState;
            }
        void setDiscreteState(DiscreteState discreteState) {
            m_discreteState = discreteState;
            }

		void computeDiscreteState(double time, const V& state)
		{
			real_type N, R,R0;
			switch (m_discreteState) {
			case F:
				if (state[3] < 0.0000150) {
					m_discreteState = S;
				}
				else m_discreteState = F;
			case S:
				computeReactions_S(N, R, time, state);
					if (N < 0) {
						m_discreteState = F;
					}
					else if (state[2] < pow(1, -10)) {
						m_discreteState = K;
					}
					else m_discreteState = S;
			case K:
				computeReactions_K(N, R, R0, time, state);
				if (N < 0) {
					m_discreteState = F;
				}
				else if (R > R0) {
					m_discreteState = S;
				}
				else m_discreteState = K;

			};
		};
          

    private:

        DiscreteState m_discreteState;

		//Объявление переменных
		const real_type g = 9.8;
		const real_type m_alf = 10; // угол наклона поверхности
		const real_type m_ow = 1; // частота колебаний
		const real_type m_f = 0.8; // трение тела о лоток
		const real_type m_mass = 1; //масса т.т.
		const real_type A = 2;//Амплитуда по х
		const real_type B = 4;//Амлитуда по y
		const real_type m_Eps = 0.5;//смещение фазы

        template< class T >
        static T sqr( T x ) {
            return x*x;
            }


        void computeReactions_S(real_type& N, real_type& R, real_type time, const V& x) const
            {
			N = m_mass*g*cos(m_alf) - m_mass*sqr(m_ow)*B*sin(m_ow*time + m_Eps);
			R = m_mass*g*sin(m_alf) - m_mass*sqr(m_ow)*A*sin(m_ow*time);

           // ASSERT(m_discreteState == S);
           // ASSERT(false); // TODO
            }
        void computeReactions_K(real_type& N, real_type& R,real_type& R0, real_type time, const V& x) const
            {
			N = m_mass*g*cos(m_alf) - m_mass*sqr(m_ow)*B*cos(m_ow*time + m_Eps);
			R = m_mass*g*sin(m_alf) - m_mass*sqr(m_ow)*A*sin(m_ow*time);
			R0 = m_mass*g*sin(m_alf);//TODO
            //ASSERT(m_discreteState == K);
            //ASSERT(false); // TODO
            }
    };

#endif // _VIBRO_TRANSPORT_H_
