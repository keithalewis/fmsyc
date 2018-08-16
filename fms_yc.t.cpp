// fms_yc.t.cpp - Test yield curve code
#include <cassert>
#include "fms_pwflat.h"

int main()
{
    void test_fms_pwflat();
    test_fms_pwflat();

    return 0;
}

void test_fms_pwflat()
{
    using namespace fms::pwflat;

    std::vector<double> t{ 1,2,3 }, f{ .1,.2,.3 };
    std::vector<double> t_2{ 1 }, f_2{ .1 };

    { // monotonic
        assert(monotonic(std::begin(t), std::end(t)));
        assert(monotonic(std::begin(f), std::end(f)));
        double f2 = f[2];
        f[2] = -1;
        assert(!monotonic(std::begin(f), std::end(f)));
        f[2] = f2;
        assert(!monotonic(std::rbegin(f), std::rend(f)));
    }
    { // forward
      //0, 0, null, null, null
        assert(isnan(value<int, double>(0, 0, nullptr, nullptr)));
        //1, 0, null, null, null
        assert(isnan(value<int, double>(1, 0, nullptr, nullptr)));
        //-1, 0, null, null, null
        assert(isnan(value<int, double>(-1, 0, nullptr, nullptr)));
        //-1, 0, null, null, 0.2
        assert(isnan(value<int, double>(-1, 0, nullptr, nullptr, 0.2)));

        int u;
        u = 1;
        double x{ 0.2 }, x_;
        //1, 0, null, null, 0.2
        x_ = fms::pwflat::value<int, double>(u, 0, nullptr, nullptr, x);
        assert(x_ == x);

        double u_[] = { -1, 0, 0.5, 1, 1.5 };
        double a_[] = { 0, 0.1, 0.1, 0.1, 0.2 };

        for (int i = 0; i < 5; i++) {
            if (i == 0 || i == 4) {
                assert(isnan(value<double, double>(u_[i], t_2.size(), t_2.data(), f_2.data())));
            }
            else {
                x_ = fms::pwflat::value<double, double>(u_[i], t_2.size(), t_2.data(), f_2.data());
                assert(x_ == a_[i]);
            }
        }

        for (int i = 0; i < 5; i++) {
            if (i == 0) {
                assert(isnan(value<double, double>(u_[i], t_2.size(), t_2.data(), f_2.data(), 0.2)));
            }
            else {
                x_ = fms::pwflat::value<double, double>(u_[i], t_2.size(), t_2.data(), f_2.data(), 0.2);
                assert(x_ == a_[i]);
            }
        }

        for (int i = 0; i < 3; ++i)
            assert(f[i] == value(t[i], t.size(), t.data(), f.data()));
    }
    { // integral
        double u;
        u = -1;
        assert(isnan(integral(u, t.size(), t.data(), f.data())));
        u = 4;
        assert(isnan(integral(u, t.size(), t.data(), f.data())));
        u = 0;
        assert(0 == integral(u, t.size(), t.data(), f.data()));
        u = 0.5;
        assert(.1*.5 == integral(u, t.size(), t.data(), f.data()));
        u = 1;
        assert(.1 == integral(u, t.size(), t.data(), f.data()));
        u = 1.5;
        assert(.1 + .2*.5 == integral(u, t.size(), t.data(), f.data()));
        u = 2.5;
        assert(.1 + .2 + .3*.5 == integral(u, t.size(), t.data(), f.data()));
        u = 3;
        assert(fabs(.1 + .2 + .3 - integral(u, t.size(), t.data(), f.data())) < 1e-10);
        //		assert (.1 + .2 + .3 != .6); 
    }
    { // discount
        double u_[] = { -.5, 0, .5, 1, 1.5, 2, 2.5, 3, 3.5 };
        double f_[] = { 0, 0, .05, .1, .2, .3, .45, .6, .7 };
        for (int i = 0; i < 9; i++) {
            if (i == 0 || i == 8) {
                assert(isnan(discount(u_[i], t.size(), t.data(), f.data())));
            }
            else {
                assert(fabs(exp(-f_[i]) - discount(u_[i], t.size(), t.data(), f.data())) < 1e-10);
            }
        }

        for (int i = 0; i < 9; i++) {
            if (i == 0) {
                assert(isnan(discount(u_[i], t.size(), t.data(), f.data(), 0.2)));
            }
            else {
                assert(fabs(exp(-f_[i]) - discount(u_[i], t.size(), t.data(), f.data(), 0.2)) < 1e-10);
            }
        }
    }
    { // spot
        double u_[] = { -.5, 0, .5, 1, 1.5, 2, 2.5, 3, 3.5 };
        double f_[] = { .1, .1, .1, .1, .2 / 1.5, .3 / 2, .45 / 2.5, .6 / 3, .7 / 3.5 };
        for (int i = 0; i < 9; i++) {
            if (i == 8) {
                assert(isnan(spot(u_[i], t.size(), t.data(), f.data())));
            }
            else {
                assert(fabs(f_[i] - spot(u_[i], t.size(), t.data(), f.data())) < 1e-10);
            }
        }

        for (int i = 0; i < 9; i++) {
            assert(fabs(f_[i] - spot(u_[i], t.size(), t.data(), f.data(), 0.2)) < 1e-10);
        }
    }
    { // present_value
        double u_[] = { 0, 1, 2, 3, 4 };
        double d_[] = { 0,
            discount(u_[1], t.size(), t.data(), f.data(), 0.2),
            discount(u_[2], t.size(), t.data(), f.data(), 0.2),
            discount(u_[3], t.size(), t.data(), f.data(), 0.2),
            discount(u_[4], t.size(), t.data(), f.data(), 0.2)
        };
        double c_[] = { 0, 1, 2, 3, 4 };

        //assert(isnan(present_value(1, u_, c_, t.size(), t.data(), f.data())));
        //assert(isnan(present_value(1, u_, c_, t.size(), t.data(), f.data(), 0.2)));

        double sum = 0;
        for (int i = 0; i < 5; i++) {
            sum += c_[i] * d_[i];
            if (i == 4) {
                double tmp = present_value<double, double>(i + 1, u_, c_, t.size(), t.data(), f.data(), 0.2);
                assert(tmp == tmp);
                assert(fabs(sum - present_value(i + 1, u_, c_, t.size(), t.data(), f.data(), 0.2)) < 1e-10);
                assert(isnan(present_value(i + 1, u_, c_, t.size(), t.data(), f.data())));
            }
            else {
                double tmp = present_value<double, double>(i + 1, u_, c_, t.size(), t.data(), f.data(), 0.2);
                assert(tmp == tmp);
                assert(fabs(sum - present_value(i + 1, u_, c_, t.size(), t.data(), f.data(), 0.2)) < 1e-10);
                assert(fabs(sum - present_value(i + 1, u_, c_, t.size(), t.data(), f.data())) < 1e-10);
            }
        }

    }
}
