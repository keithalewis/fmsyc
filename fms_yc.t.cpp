// fms_yc.t.cpp - Test yield curve code
#include <cassert>
#include <vector>
#include "fms_pwflat.h"

template<class T>
void test_fms_pwflat()
{
    using namespace fms::pwflat;

    T eps = std::numeric_limits<T>::epsilon();
    std::vector<T> t{ 1,2,3 }, f{ T(.1),T(.2),T(.3) };
    std::vector<T> t_2{ 1 }, f_2{ T(.1) };

    { // strictly_increasing
        assert(strictly_increasing(3, t.data()));
        assert(strictly_increasing(3, f.data()));
        T f2 = f[2];
        f[2] = -1;
        assert(!strictly_increasing(3, f.data()));
        f[2] = f2;
        assert(strictly_increasing(3, f.data()));
    }
    { // forward
      //0, 0, null, null, null
        assert(isnan(value<int, T>(0, 0, nullptr, nullptr)));
        //1, 0, null, null, null
        assert(isnan(value<int, T>(1, 0, nullptr, nullptr)));
        //-1, 0, null, null, null
        assert(isnan(value<int, T>(-1, 0, nullptr, nullptr)));
        //-1, 0, null, null, 0.2
        assert(isnan(value<int, T>(-1, 0, nullptr, nullptr, T(0.2))));

        int u;
        u = 1;
        T x{ T(0.2) }, x_;
        //1, 0, null, null, 0.2
        x_ = fms::pwflat::value<int, T>(u, 0, nullptr, nullptr, x);
        assert(x_ == x);

        T u_[] = { T(-1), T(0), T(0.5), T(1), T(1.5) };
        T a_[] = { T(0), T(0.1), T(0.1), T(0.1), T(0.2) };

        for (int i = 0; i < 5; i++) {
            if (i == 0 || i == 4) {
                assert(isnan(value<T, T>(u_[i], t_2.size(), t_2.data(), f_2.data())));
            }
            else {
                x_ = fms::pwflat::value<T, T>(u_[i], t_2.size(), t_2.data(), f_2.data());
                assert(x_ == a_[i]);
            }
        }

        for (int i = 0; i < 5; i++) {
            if (i == 0) {
                assert(isnan(value<T, T>(u_[i], t_2.size(), t_2.data(), f_2.data(), T(0.2))));
            }
            else {
                x_ = fms::pwflat::value<T, T>(u_[i], t_2.size(), t_2.data(), f_2.data(), T(0.2));
                assert(x_ == a_[i]);
            }
        }

        for (int i = 0; i < 3; ++i)
            assert(f[i] == value(t[i], t.size(), t.data(), f.data()));
    }
    { // integral
        T u;
        u = -1;
        assert(isnan(integral(u, t.size(), t.data(), f.data())));
        u = 4;
        assert(isnan(integral(u, t.size(), t.data(), f.data())));
        u = 0;
        assert(0 == integral(u, t.size(), t.data(), f.data()));
        u = 0.5;
        assert(T(.1)*T(.5) == integral(u, t.size(), t.data(), f.data()));
        u = 1;
        assert(T(.1) == integral(u, t.size(), t.data(), f.data()));
        u = 1.5;
        assert(T(.1) + T(.2)*T(.5) == integral(u, t.size(), t.data(), f.data()));
        u = 2.5;
        assert(T(.1) + T(.2) + T(.3)*T(.5) == integral(u, t.size(), t.data(), f.data()));
        u = 3;
        assert(fabs(T(.1) + T(.2) + T(.3) - integral(u, t.size(), t.data(), f.data())) < 2*eps);
        //		assert (.1 + .2 + .3 != .6); 
    }
    { // discount
        T u_[] = { T(-.5), T(0), T(.5), T(1), T(1.5), T(2), T(2.5), T(3), T(3.5) };
        T f_[] = { T(0), T(0), T(.05), T(.1), T(.2), T(.3), T(.45), T(.6), T(.7) };
        for (int i = 0; i < 9; i++) {
            if (i == 0 || i == 8) {
                assert(isnan(discount(u_[i], t.size(), t.data(), f.data())));
            }
            else {
                assert(fabs(exp(-f_[i]) - discount(u_[i], t.size(), t.data(), f.data())) < 2*eps);
            }
        }

        for (int i = 0; i < sizeof(u_)/sizeof(u_[0]); i++) {
            if (i == 0) {
                assert(isnan(discount(u_[i], t.size(), t.data(), f.data(), T(0.2))));
            }
            else {
                assert(fabs(exp(-f_[i]) - discount(u_[i], t.size(), t.data(), f.data(), T(0.2))) < 2*eps);
            }
        }
    }
    { // spot
        T u_[] = { T(-.5), T(0), T(.5), T(1), T(1.5), T(2), T(2.5), T(3), T(3.5) };
        T f_[] = { T(.1), T(.1), T(.1), T(.1), T(.2 / 1.5), T(.3 / 2), T(.45 / 2.5), T(.6 / 3), T(.7 / 3.5) };
        for (int i = 0; i < 9; i++) {
            if (i == 8) {
                assert(isnan(spot(u_[i], t.size(), t.data(), f.data())));
            }
            else {
                assert(fabs(f_[i] - spot(u_[i], t.size(), t.data(), f.data())) < 2*eps);
            }
        }

        for (int i = 0; i < sizeof(u_) / sizeof(u_[0]); i++) {
            assert(fabs(f_[i] - spot(u_[i], t.size(), t.data(), f.data(), T(0.2))) < 2*eps);
        }
    }
    { // present_value
        T u_[] = { T(0), T(1), T(2), T(3), T(4) };
        T d_[] = { T(0),
            discount(u_[1], t.size(), t.data(), f.data(), T(0.2)),
            discount(u_[2], t.size(), t.data(), f.data(), T(0.2)),
            discount(u_[3], t.size(), t.data(), f.data(), T(0.2)),
            discount(u_[4], t.size(), t.data(), f.data(), T(0.2))
        };
        T c_[] = { T(0), T(1), T(2), T(3), T(4) };

        //assert(isnan(present_value(1, u_, c_, t.size(), t.data(), f.data())));
        //assert(isnan(present_value(1, u_, c_, t.size(), t.data(), f.data(), 0.2)));

        T sum = 0;
        for (int i = 0; i < 5; i++) {
            sum += c_[i] * d_[i];
            if (i == 4) {
                T tmp = present_value<T, T>(i + 1, u_, c_, t.size(), t.data(), f.data(), T(0.2));
                assert(tmp == tmp);
                assert(fabs(sum - present_value(i + 1, u_, c_, t.size(), t.data(), f.data(), T(0.2))) < 2*eps);
                assert(isnan(present_value(i + 1, u_, c_, t.size(), t.data(), f.data())));
            }
            else {
                T tmp = present_value<T, T>(i + 1, u_, c_, t.size(), t.data(), f.data(), T(0.2));
                assert(tmp == tmp);
                assert(fabs(sum - present_value(i + 1, u_, c_, t.size(), t.data(), f.data(), T(0.2))) < 2*eps);
                assert(fabs(sum - present_value(i + 1, u_, c_, t.size(), t.data(), f.data())) < 1e-10);
            }
        }

    }
}
int main()
{
    test_fms_pwflat<float>();
    test_fms_pwflat<double>();

    return 0;
}
