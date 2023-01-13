#include <iostream>
#include <vector>

class CubicSpline{
public:

    CubicSpline() = default;
    ~CubicSpline() = default;

    double u(const std::vector<double>& f, const std::vector<double>& h, int i, int j)
    {
        return (f[j] - f[i]) / h[j];
    }

    double u(const std::vector<double>& f, const std::vector<double>& h, int i, int j, int k)
    {
        return (u(f, h, j, k) - u(f, h, i, j)) / (h[k] + h[j]);
    }

    void interpolate(const std::vector<double>& x, const std::vector<double>& f)
    {
        std::cout << x[0] << " " << x[1] << std::endl;

        std::vector<double> a, b, c, d, h;

        //For testing
        a.reserve(5);
        b.reserve(5);
        c.reserve(5);
        d.reserve(5);
        h.reserve(5);

        int N_1 = (int) x.size();
        int N = N_1 - 1;

        std::cout << "N_1: " << N_1 << std::endl;

        h[0] = 0;
        for (int i = 1; i < N_1; ++i)
            h[i] = x[i] - x[i - 1];

        a[0] = 0;
        b[0] = 2;
        c[0] = h[2] / (h[1] + h[2]);
        d[0] = 6 * u(f, h, 0, 1, 2);

        for (int i = 1; i < N - 2; ++i)
        {
            a[i] = h[i] / (h[i] + h[i + 1]);
            b[i] = 2;
            c[i] = h[i + 1] / (h[i] + h[i + 1]);
            d[i] = u(f, h, i - 1, i, i + 1);
        }
        a[N - 2] = h[N - 1] / (h[N - 1] + h[N]);
        b[N - 2] = 2;
        c[N - 2] = 0;
        d[N - 2] = 6 * u(f, h, N - 2, N - 1, N);

        c = solveTriangonalSlae(a, b, c, d);
        c.insert(c.begin(), 0);
        c[N] = 0;
        b[0] = 0;
        a[0] = 0;
        d[0] = 0;

        for (int i = 1; i < N_1; ++i)
        {
            b[i] = (c[i] / 3 + c[i - 1] / 6) * h[i] + u(f, h, i - 1, i);
            d[i] = (c[i] - c[i - 1]) / h[i];
            a[i] = f[i];
            //тупой вывод для двух точек
            std::cout << "\na: " << a[i] << " b: " << b[i] << " c: " << c[i] << " d: " << d[i];
        }
    }

    std::vector<double> solveTriangonalSlae(const std::vector<double>& a, const std::vector<double>& b, const std::vector<double>& c, const std::vector<double>& d)
    {
        std::vector<double> x, p ,q;
        p[0] = 0;
        q[0] = 0;
        int M_1 = b.size();
        int M = M_1 - 1;

        //For testing
        x.reserve(5);
        p.reserve(5);
        q.reserve(5);

        if (a.size() == M_1 && c.size() == M_1)
        {
            for (int i = 0; i < M; ++i)
            {
                p[i + 1] = - c[i] / (a[i] * p[i] + b[i]);
                q[i + 1] = (d[i] - a[i] * q[i]) / (a[i] * p[i] + b[i]);
                x[i] = 0;
            }
            x[M] = (d[M] - a[M] * q[M]) / (p[M] * a[M] + b[M]);
            for(unsigned int i = M - 1; i >= 0; --i)
            {
                x[i] = x[i + 1] * p[i + 1] + q[i + 1];
            }
        }
        else
            std::cout << "Length of vectors 'a', 'b' and 'c' aren't equal\n";
        return x;
    }

};

int main() {
    //По итогу ошибка в синтаксисе обращения к вектору, дома перепишу

    std::vector<double> x, f;
    x.reserve(2);
    f.reserve(2);
    x[0] = 1;
    x[1] = 2;
    f[0] = 1;
    f[1] = 2;
    std::cout << "N_1 start " << x.size() << std::endl;

    CubicSpline c;
    c.interpolate(x, f);

    return 0;
}
