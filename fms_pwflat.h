// fms_pwflat.h - piecewise flat curve
/*
    f(t) = f[i] if t[i-1] < t <= t[i]
         = _f   if t > t[n-1]
    and undefined if t < 0

    |                                   _f
    |        f[1]             f[n-1] (--------
    | f[0] (----- ...       (------]
    [------]      ... ------]
    |
    0-----t[0]--- ... ---t[n-2]---t[n-1]
*/
#pragma once
#include <cmath>     // exp
#include <algorithm> // adjacent_find
#include <limits>    // quiet_Nan()
#include <numeric>   // upper/lower_bound
#include <vector>
#include <gsl/span>

namespace fms {
    namespace pwflat {

        template<class X>
        constexpr X NaN() { return std::numeric_limits<X>::quiet_NaN(); }

        // strictly increasing values
        template<class I>
        inline bool monotonic(I b, I e) noexcept
        {
            using T = typename std::iterator_traits<I>::value_type;

            return e == std::adjacent_find(b, e, [](const T& t0, const T&t1) { return t0 >= t1; });
        }
        template<class T>
        inline bool monotonic(size_t n, const T* t) noexcept
        {
            return monotonic(t, t + n);
        }

        // piecewise flat curve
        // return f[i] if t[i-1] < u <= t[i], _f if u > t[n-1]
        // assumes t[i] monotonically increasing
        template<class T, class F>
        inline F value(const T& u, size_t n, const T* t, const F* f, const F& _f = NaN<F>())
        {
            if (u < 0)
                return NaN<F>();
            if (n == 0)
                return _f;

            auto ti = std::lower_bound(t, t + n, u);

            return ti == t + n ? _f : f[ti - t];
        }

        // int_0^u f(t) dt
        template<class T, class F>
        inline F integral(const T& u, size_t n, const T* t, const F* f, const F& _f = std::numeric_limits<F>::quiet_NaN())
        {
            if (u < 0)
                return NaN<F>();

            F I{ 0 };
            T t_{ 0 };

            size_t i;
            for (i = 0; i < n && t[i] <= u; ++i) {
                I += f[i] * (t[i] - t_);
                t_ = t[i];
            }
            I += (n == 0 || u > t[n - 1] ? _f : f[i]) *(u - t_);

            return I;
        }

        // discount D(u) = exp(-int_0^u f(t) dt)
        template<class T, class F>
        inline F discount(const T& u, size_t n, const T* t, const F* f, const F& _f = std::numeric_limits<F>::quiet_NaN())
        {
            return exp(-integral(u, n, t, f, _f));
        }

        // spot r(u) = (int_0^u f(t) dt)/u
        template<class T, class F>
        inline F spot(const T& u, size_t n, const T* t, const F* f, const F& _f = std::numeric_limits<F>::quiet_NaN())
        {
            return u <= t[0] ? f[0] : integral(u, n, t, f, _f) / u;
        }

        // NVI class
        // https://en.wikibooks.org/wiki/More_C%2B%2B_Idioms/Non-Virtual_Interface
        template<class T = double, class F = double>
        class interface {
            size_t n;
            const T* t;
            const F* f;
        public:
            typedef T time_type;
            typedef F rate_type;
            interface(size_t n = 0, const T* t = 0, const F* f = 0)
                : n(n), t(t), f(f)
            { }
            virtual ~interface() { }
            size_t   size() const { return _size(); }
            const T* time() const { return _time(); }
            const F* rate() const { return _rate(); }
            T value(T u, const F& _f = NaN<F>()) const
            {
                return pwflat::value(u, size(), time(), rate(), _f);
            }
            T operator()(T u, const F& _f = NaN<F>()) const
            {
                return value(u, _f);
            }
            T integral(T u, const F& _f = NaN<F>()) const
            {
                return integral(u, size(), time(), rate(), _f);
            }
            T spot(T u, const F& _f = NaN<F>()) const
            {
                return spot(u, size(), time(), rate(), _f);
            }
        private:
            // override in base class
            virtual size_t   _size() const
            {
                return n;
            }
            virtual const T* _time() const
            {
                return t;
            }
            virtual const F* _rate() const
            {
                return f;
            }
        };

        template<class T = double, class F = double>
        class curve : public interface<T, F> {
            std::vector<T> t_;
            std::vector<F> f_;
        public:
            curve()
            { }
            curve(size_t n, const T* t, const F* f)
                : t_(t, t + n), f_(f, f + n)
            { }
            curve(const std::vector<T>& t, const std::vector<F>& f)
                : t_(t), f_(f)
            {
                ensure(t_.size() == f_.size());
            }
            curve& push_back(const T& _t, const F& _f)
            {
                ensure(t_.size() == 0 || _t > t_.back());

                t_.push_back(_t);
                f_.push_back(_f);

                return *this;
            }
            curve& push_back(const std::pair<T, F>& p)
            {
                return push_back(p.first, p.second);
            }
        private:
            size_t _size() const override
            {
                return t_.size();
            }
            const T* _time() const override
            {
                return t_.data();
            }
            const F* _rate() const override
            {
                return f_.data();
            }
        };

        // value of instrument having cash flow c[i] at time u[i]
        template<class T, class F>
        inline F present_value(size_t m, const T* u, const F* c, size_t n, const T* t, const F* f, const F& _f = NaN<F>())
        {
            F p{ 0 };

            for (size_t i = 0; i < m; ++i)
                p += c[i] * pwflat::discount(u[i], n, t, f, _f);

            return p;
        }

        // derivative of present value wrt parallel shift of forward curve
        template<class T, class F>
        inline F duration(size_t m, const T* u, const F* c, size_t n, const T* t, const F* f, const F& _f = NaN())
        {
            F d{ 0 };

            for (size_t i = 0; i < m; ++i) {
                d -= u[i] * c[i] * pwflat::discount(u[i], n, t, f, _f);
            }

            return d;
        }

        // derivative of present value wrt parallel shift of forward curve after last curve time
        template<class T, class F>
        inline F partial_duration(size_t m, const T* u, const F* c, size_t n, const T* t, const F* f, const F& _f = NaN())
        {
            F d{ 0 };

            // first cash flow past end of forward curve
            size_t i0 = (n == 0) ? 0 : std::lower_bound(u, u + m, t[n - 1]) - u;
            double t0 = (n == 0) ? 0 : t[n - 1];
            for (size_t i = i0; i < m; ++i) {
                d -= (u[i] - t0)*c[i] * pwflat::discount(u[i], n, t, f, _f);
            }

            return d;
        }
    } // pwflat

} // fms
