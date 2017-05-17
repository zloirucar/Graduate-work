#include "ode_num_int/OdeSolverConfiguration.h"
#include "ode_num_int/OdeSolverRK4.h"
#include "reg.h"
#include "VibroTransport.h"

#include <iostream>

//namespace ctm {
//CTM_DECL_IMPLEMENTATION_TEMPLATE_TRAITS( VibroTransport, "VibroTransport" )
//}

int main()
    {
    using namespace std;
    using namespace ctm;
    using namespace math;

    registerTypes();

    typedef VectorData<double> VD;
    typedef VectorTemplate<VD> V;

    VibroTransport<VD>::Registrator vibroTransportRegistrator( "VibroTransport" );

    try {
        OdeSolverConfiguration<VD> cfg;
        cfg.setValue("rhs", "VibroTransport");
        cfg.setValue("output_con", "con_solution");
        cfg.setValue("solver", "dopri_56");
        cfg.setValue("solver.stepsizectl.tolerance", 1e-7);
        cfg.setValue("time", 10);
        V x0( 4 );
        x0[0] = 0;
		x0[1] = 0;
		x0[2] = 1;
		x0[3] = 1;
        auto sc = cfg.apply( set<unsigned int>(), 0, x0 );
        auto vt = dynamic_cast<VibroTransport<VD>*>(sc.solver()->odeRhs().get());
        ASSERT(vt);
        vt->setDiscreteState(VibroTransport<VD>::F);
        vt->computeDiscreteState( 0, x0 );
        solveOde( &cfg, &sc );
        return 0;
    }
    catch( const std::exception& e ) {
        cerr << "ERROR: " << e.what() << endl;
        return 1;
        }
    }
