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
        cfg.setValue("solver.stepsizectl.tolerance", 1e-5);
        cfg.setValue("time", 10);
        cfg.setValue("output_timing", 0.1);
        V x0( 4 );
        x0[0] = 0;
		x0[1] = 0;
		x0[2] = 0;
		x0[3] = 0;
        auto sc = cfg.apply( set<unsigned int>(), 0, x0 );
        auto vt = dynamic_cast<VibroTransport<VD>*>(sc.solver()->odeRhs().get());
        ASSERT(vt);
        vt->setDiscreteState(VibroTransport<VD>::F);
        vt->computeDiscreteState( 0, x0 );
        auto stateName = [&vt]() -> string {
            const char *stateNames[] = { "F", "SL", "SR", "K", "F1" };
            return stateNames[vt->discreteState()];
        };
        cout << "INITIAL DISCRETE STATE: " << stateName() << endl;
//        sc.solver()->odeSolverPostObservers.add([&](const OdeSolverPostObserverArg<VD>& arg) {
//            if (!(arg.stepAccepted() && arg.stepTruncated()))
//                return;
//            static unsigned int counter = 0;
//            ++counter;
//            cout << counter << ", t=" << arg.solver()->initialTime() << ", TRANSITION: " << arg.izfTrunc() << ":" << arg.transitionType() << " -> " << stateName() << endl;
//        });
        solveOde( &cfg, &sc );
        return 0;
    }
    catch( const std::exception& e ) {
        cerr << "ERROR: " << e.what() << endl;
        return 1;
        }
    }
