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
            switch (m_discreteState) {
                case F:
                    ASSERT(false); // TODO
                    break;
                case S:
                    {
                    real_type N, R;
                    computeReactions_S(N, R, time, x);
                    ASSERT(false); // TODO
                    break;
                    }
                case K:
                    ASSERT(false); // TODO
                    break;
                default:
                    ASSERT(false);
                }

			dst.resize(4);
			auto v1 = x[2];
			auto v2 = x[3];
			auto v = sqrt(v1*v1 + v2*v2);
            auto m_F1 = 0;// TODO -g*sin(m_alf) + R + pow(m_ow, 2)*A*sin(m_ow*time); //*************************Где нужно высчитывать R и N?
            auto m_F2 = 0;// TODO -g*cos(m_alf) + N + pow(m_ow, 2)*B*cos(m_ow*time + m_Eps);
			dst[0] = v1; //************************************************************это x с точкой?
			dst[1] = v2; // ***********************************************************это y с точкой?
			dst[2] = m_F1 / m_mass;
			dst[3] =m_F2/m_mass;
            this->odeRhsPostObservers( time, x, dst, this );
		};

		

        virtual unsigned int zeroFuncCount() const
            {
            return 2;
            }

        virtual void zeroFunctions( V& dst, real_type /*time*/, const V& x ) const
            {
			/*if (m_discreteState == F) {
				i1 = pow(pow(x[0], 2) + pow(x[1], 2), 0.5); // индикатор для отслеживания расстояния до поверхности
			};

			if (m_discreteState == S)
			{
				i1 = N; // индкатор нормальной реакции
				i2 = x[2]; //индикатор тангенциальной скорости
			};
			if (m_discreteState == K) {
				i1 = N; // индикатор нормальной реакции
				i2 = R; // индикатор горизонтальной реакции
			}*/
            }

        virtual void switchPhaseState( const int* /*transitions*/, real_type /*time*/, V& /*x*/ ) {
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

            ASSERT(false);// TODO

			/*if (R==0 && N==0) {
				m_discreteState = F;
			};

			if (state[3]!= 0 &&   ) {
				m_discreteState = S;
			};

			if () {
				m_discreteState = K;
			};
			*/

            
            }

    private:

        DiscreteState m_discreteState;
		double i1, i2; //************************************какого типа должны быть индикторы?

		//Объявление переменных
		const real_type g = 9.8;
		//const real_type N = m*g*cos(alf) - pow(ow, 2)*B*cos(ow*time + Eps);
		//const real_type R = m*g*sin(alf) - pow(ow, 2)*A*sin(ow*time);;
		const real_type m_alf = 10; // угол наклона поверхности
		const real_type m_ow = 1; // частота колебаний
		const real_type m_f = 0.8; // трение тела о лоток
		const real_type m_mass = 1; //масса т.т.
		const real_type A = 2;//Амплиткда по х
		const real_type B = 4;//Амлитуда по y
		const real_type m_Eps = 0.5;//смещение фазы

        template< class T >
        static T sqr( T x ) {
            return x*x;
            }
        void computeReactions_S(real_type& N, real_type& R, real_type time, const V& x) const
            {
            ASSERT(m_discreteState == S);
            ASSERT(false); // TODO
            }
        void computeReactions_K(real_type& N, real_type& R, real_type time, const V& x) const
            {
            ASSERT(m_discreteState == K);
            ASSERT(false); // TODO
            }
    };

#endif // _VIBRO_TRANSPORT_H_
