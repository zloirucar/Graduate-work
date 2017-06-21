#include "ode_num_int/OdeSolverConfiguration.h"
#include "ode_num_int/OdeSolverRK4.h"
#include "reg.h"
#include "VibroTransport.h"

#include <iostream>
#include <sstream>

#include <QApplication>
#include <QImage>
#include <QColor>
#include <QTime>
#include <QPainter>
#include <qmath.h>

//namespace ctm {
//CTM_DECL_IMPLEMENTATION_TEMPLATE_TRAITS( VibroTransport, "VibroTransport" )
//}

// Computes average from values passed in using operator <<
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

// Returns indexed values in the specified range (see operator[])
class IndexedRange
    {
    public:
        IndexedRange(double from, double to, int count) :
            m_from(from),
            m_to(to),
            m_count(count)
            {
            ASSERT(count > 1);
            }
        double operator[](int index) const {
            ASSERT(index >= 0   &&   index < m_count);
            double t = static_cast<double>(index) / (m_count-1);
            return m_from*(1-t) + m_to*t;
            }

    private:
        double m_from;
        double m_to;
        int m_count;
    };

// Cancels ODE solver by throwing an exception if it takes too long
class DelayBasedOdeSolverCancel
    {
    public:
        typedef VectorData<double> VD;
        DelayBasedOdeSolverCancel(OdeSolver<VD> *solver, int maxDelayMsec) :
            m_maxDelayMsec(maxDelayMsec)
        {
            solver->odeSolverPostObservers.add([this](const OdeSolverPostObserverArg<VD>& arg) {
                if (--m_count)
                    return;
                m_count = Granularity;
                if(m_time.elapsed() > m_maxDelayMsec)
                    throw cxx::exception("Triggered delay-based ODE solver cancel");
                });
            reset();
        }
        void reset() {
            m_count = Granularity;
            m_time.start();
            }

    private:
        enum { Granularity = 10000 };
        int m_maxDelayMsec;
        int m_count;
        QTime m_time;
    };

class VtParameterSetter
    {
    public:
        typedef VectorData<double> VD;
        typedef VibroTransport<VD> VT;
        virtual void operator()(VT& vt, double p1, double p2) const = 0;
    };

class VtParameterSetter_ba_beta : public VtParameterSetter
    {
    public:
        explicit VtParameterSetter_ba_beta(double a) : m_a(a) {}

        void operator()(VT& vt, double p1, double p2) const override {
            const auto& ba = p1;
            const auto& beta = p2;
            const auto& a = m_a;
            auto b = a*ba;
            auto sb = sin(beta);
            auto cb = cos(beta);
            auto r = sqrt(sqr(a*cb) + sqr(b*sb));
            auto spsi = a*cb/r;
            auto cpsi = b*sb/r;
            auto A =     -a*spsi*cb - b*cpsi*sb;
            auto Bcphi = -a*spsi*sb + b*cpsi*cb;
            auto Bsphi =  a*cpsi*sb + b*spsi*cb;
            auto B = sqrt(sqr(Bcphi) + sqr(Bsphi));
            auto phi = atan2(Bsphi, Bcphi);
            vt.setA(A);
            vt.setB(B);
            vt.setEps(phi);
            }

    private:
        double m_a;
    };

class VtParameterSetter_alpha_B : public VtParameterSetter
    {
    public:
        explicit VtParameterSetter_alpha_B(double A) : m_A(A) {}

        void operator()(VT& vt, double p1, double p2) const override {
            const auto& alpha = p1;
            const auto& B = p2;
            const auto& A = m_A;
            vt.setAlpha(alpha);
            vt.setA(A);
            vt.setB(B);
            }

    private:
        double m_A;
    };

class VtIndexedParameterSetter
    {
    public:
        typedef VectorData<double> VD;
        typedef VibroTransport<VD> VT;

        VtIndexedParameterSetter(
                VT& vt, VtParameterSetter& ps,
                const IndexedRange& r1, const IndexedRange& r2) :
            m_vt(vt), m_ps(ps), m_r1(r1), m_r2(r2) {}

        void operator()(int i, int j) {
            m_ps(m_vt, m_r1[i], m_r2[j]);
            }

    private:
        VibroTransport<VD>& m_vt;
        VtParameterSetter& m_ps;
        IndexedRange m_r1;
        IndexedRange m_r2;
    };

void saveScaleImage(double min, double max)
    {
    int H = 50, W = 300;
    IndexedRange rc(0, 5./6, W);
    QImage img(W, H, QImage::Format_RGB32);
    img.fill(Qt::white);
        {
        QPainter p(&img);
        auto h = p.boundingRect(QRect(0,0,W,H), "X").height() + 5;
        for (int x=0; x<W; ++x)
            p.fillRect(QRect(x,0,1,H-h), QColor::fromHsvF(rc[x], 1, 1));
        if (min < 0 && max > 0) {
            int x = static_cast<int>((0-min)/(max-min)*W);
            p.fillRect(QRect(x,0,1,H-h+2), Qt::black);
            p.drawText(QRectF(x-20,H-h+5,40,h-5), Qt::AlignCenter, "0");
            }

        p.drawText(QRectF(0,H-h+5,W,h-5), Qt::AlignLeft | Qt::AlignVCenter, QString::number(min));
        p.drawText(QRectF(0,H-h+5,W,h-5), Qt::AlignRight | Qt::AlignVCenter, QString::number(max));
        }
    img.save("scale.png");
    }

int main(int argc, char *argv[])
    {
    QApplication app(argc, argv);
    setlocale( LC_NUMERIC, "C" );
//    saveScaleImage(-3.99932, 0.157367);
//    return 0;

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
//        cfg.setValue("output_con", "stats_full");
//		cfg.setValue("output_con.file_name", "myfile.txt");
        cfg.setValue("solver", "dopri_56");
        cfg.setValue("solver.stepsizectl.tolerance", 1e-5);
        cfg.setValue("solver.stepsizectl.max_threshold", 1e-3);
        cfg.setValue("solver.h_init", 0.0001);
        // cfg.setValue("solver.h_min", 1e-7);
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


//        for (auto alpha = 0.0; alpha<M_PI/2 + 1e-4; alpha += M_PI/12) {
//            cout << "---- alpha = " << alpha << endl;
//            vt->setAlpha(alpha);
//            sc.solver()->setInitialState(0, x0);
//            vt->computeDiscreteState(0, x0);
//            v1 = Average();
//            v2 = Average();
//            solveOde( &cfg, &sc );
//            cout << "v1: " << v1.toString() << endl
//                 << "v2: " << v2.toString() << endl;
//            if (fabs(v1.average())*2 < -v2.average() && v2.average() < -4)
//                cout << "It's falling!!!" << endl;
//            cout << "----" << endl;
//            }

        int ind = 100;
        vector<vector<double>> V(ind, vector<double>(ind));
        enum Status { Ok, SolverFailed, Accelerating };
        vector<vector<Status>> status(ind, vector<Status>(ind));
        double max = 0;
        double min = 0;
        bool hasMinMax = false;

        vt->setow(16*2*M_PI);

        // Diagram x = beta, y = ba
        vt->setAlpha(10*M_PI/180);
        VtParameterSetter_ba_beta vtParSetter(
                    0.0011 // Большая полуось эллипса, по которому движется лоток
                    );
        VtIndexedParameterSetter vtIndexedParSetter(
                    *vt, vtParSetter,
                    IndexedRange(-1, 1, ind),           // Отношение малой полуоси эллипса, по которому движется лоток, к его большой полуоси
                    IndexedRange(-M_PI/2, M_PI/2, ind)  // Угол между большой полуосью эллипса и плоскостью лотка
                    );

//        // Diagram y = B, y = alpha
//        vt->setAlpha(10*M_PI/180);
//        vt->setEps(-M_PI/2);
//        VtParameterSetter_alpha_B vtParSetter(
//                    0.002                           // A
//                    );
//        VtIndexedParameterSetter vtIndexedParSetter(
//                    *vt, vtParSetter,
//                    IndexedRange(0, M_PI/2, ind),   // alpha
//                    IndexedRange(0, 0.002, ind)     // B
//                    );

        // Skip points that compute for more than one second
        DelayBasedOdeSolverCancel solverDelayCancel(sc.solver().get(), 1000);

        for(int i=0; i<ind;i++){
            cout << "line " << (i+1) << "/" << ind << ":\t";
            for(int j=0;j<ind;j++) {
                vtIndexedParSetter(i, j);
                sc.solver()->setInitialState(0,x0);
                vt->computeDiscreteState(0, x0);
                v1 = Average();
                v2 = Average();
                try {
                    solverDelayCancel.reset();
                    solveOde( &cfg, &sc );
                    cout << "*";
                    auto v = V[i][j] = v2.average();
                    if (fabs(v1.average())*2 < -v && v < -4)
                        status[i][j] = Accelerating;
                    else {
                        status[i][j] = Ok;
                        if (hasMinMax) {
                            if (v>max)
                                max = v;
                            else
                                if (v<min) min = v;
                            }
                        else {
                            min = max = v;
                            hasMinMax = true;
                            }
                        }
                }
                catch(const std::exception&) {
                    cout << ".";
                    status[i][j] = SolverFailed;
                    }
                cout << std::flush;
            }
            cout << endl;
        }


          QImage img(ind, ind, QImage::Format_ARGB32);
          for(int i=0; i<ind;i++){
              for(int j=0;j<ind;j++) {
                // rgb= qRgb(255/(abs(max)+abs(min))*(V[i][j]+min),0,0);
                  QRgb rgb(0);
                  switch (status[i][j]) {
                      case Ok: {
                          auto v = (V[i][j] - min) / (max - min);
                          rgb = QColor::fromHsvF(v*5./6, 1, 1).rgb();
                          break;
                      }
                      case SolverFailed:
                          rgb = 0xff000000;
                          break;
                      case Accelerating:
                          rgb = 0xffffffff;
                          break;
                      default:
                          ASSERT(false);
                          break;
                }
                img.setPixel(j, ind-1-i, rgb);
            }
          }
          img.save ("picture.png");
          saveScaleImage(min, max);
          cout << "==== min: " << min << ", max: " << max << endl;

        return 0;
    }
    catch( const std::exception& e ) {
        cerr << "ERROR: " << e.what() << endl;
        return 1;
        }
    }
