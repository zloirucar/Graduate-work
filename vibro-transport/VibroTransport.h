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

        explicit VibroTransport()
            {}

        virtual unsigned int secondOrderVarCount() const {
            return 1;
            }

        virtual unsigned int firstOrderVarCount() const {
            return 0;
            }

        virtual void rhs( V& dst, real_type time, const V& x ) const
            {
            const real_type g = 9.8;
            const real_type c = 1000;
            this->odeRhsPreObservers( time, x, this );
            dst.resize( 2 );
            dst[0] = x[1];
            dst[1] = -g;
            if (x[0] < 0)
                dst[1] -= c*x[0];
            this->odeRhsPostObservers( time, x, dst, this );
            }

        virtual unsigned int zeroFuncCount() const
            {
            return 1;
            }

        virtual void zeroFunctions( V& dst, real_type /*time*/, const V& x ) const
            {
            dst.resize( 1 );
            dst[0] = x[0];
            }

        virtual void switchPhaseState( const int* /*transitions*/ ) {
            }
        virtual std::string describeZeroFunction( unsigned int /*index*/ ) const {
            return std::string();
            }
    };

#endif // _VIBRO_TRANSPORT_H_
