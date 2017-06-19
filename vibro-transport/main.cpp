#include "ode_num_int/OdeSolverConfiguration.h"
#include "ode_num_int/OdeSolverRK4.h"
#include "reg.h"
#include "VibroTransport.h"

#include <iostream>
#include <sstream>

#include <QApplication>
#include <QImage>

//namespace ctm {
//CTM_DECL_IMPLEMENTATION_TEMPLATE_TRAITS( VibroTransport, "VibroTransport" )
//}

class Average
    {
    public:
        Average() : m_v(0), m_n(0) {}

        void operator<<(double x) {
            m_v += x;
            ++m_n;
            }
        double average() const {
            return m_v / m_n;
            }

        int n() const {
            return m_n;
            }

        std::string toString() const {
            std::ostringstream s;
            s << "average: " << average() << ", n: " << m_n;
            return s.str();
            }

    private:
        double m_v;
        int m_n;
    };

int main(int argc, char *argv[])
    {
    using namespace std;
    using namespace ctm;
    using namespace math;

    registerTypes();

    typedef VectorData<double> VD;
    typedef VectorTemplate<VD> V;

    VibroTransport<VD>::Registrator vibroTransportRegistrator( "VibroTransport" );

    try {
        const double t0 = 1;
        const double t1 = t0 + 1.5;
        const double dt = 0.5;
        const double T = t1 + dt;

        OdeSolverConfiguration<VD> cfg;
        cfg.setValue("rhs", "VibroTransport");
//        cfg.setValue("output_con", "con_solution");
//		cfg.setValue("output_con.file_name", "myfile.txt");
        cfg.setValue("solver", "dopri_56");
        cfg.setValue("solver.stepsizectl.tolerance", 1e-5);
        cfg.setValue("solver.h_init", 0.001);
//        cfg.setValue("time", 60);
//        cfg.setValue("output_timing", 0.1);
        cfg.setValue("time", T);
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
            const char *stateNames[] = { "F", "SL", "SR", "K" };
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

        Average v1, v2;
        sc.solver()->odeSolverPostObservers.add([&](const OdeSolverPostObserverArg<VD>& arg) {
            if (!arg.stepAccepted())
                return;
            auto t = arg.solver()->initialTime();
            if (t < t0)
                return;
            auto v = arg.solver()->initialState()[2];
            if (t < t0+dt)
                v1 << v;
            else if (t < t1)
                return;
            else if (t < t1 + dt)
                v2 << v;
        });


        for (auto alpha = 0.0; alpha<M_PI/2 + 1e-4; alpha += M_PI/12) {
            cout << "---- alpha = " << alpha << endl;
            vt->setAlpha(alpha);
            sc.solver()->setInitialState(0, x0);
            vt->computeDiscreteState(0, x0);
            v1 = Average();
            v2 = Average();
            solveOde( &cfg, &sc );
            cout << "v1: " << v1.toString() << endl
                 << "v2: " << v2.toString() << endl;
            if (fabs(v1.average())*2 < -v2.average() && v2.average() < -4)
                cout << "It's falling!!!" << endl;
            cout << "----" << endl;
            }

        QApplication a(argc, argv);
        QImage img(800, 600, QImage::Format_ARGB32);

        return 0;
    }
    catch( const std::exception& e ) {
        cerr << "ERROR: " << e.what() << endl;
        return 1;
        }
    }
