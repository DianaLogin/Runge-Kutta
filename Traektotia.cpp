#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>
#include <numbers>

double pi = std::numbers::pi_v<double>;

struct Y
{
    double x, y, v, o;

    Y(double x = 0.0, double y = 0.0, double v = 0.0, double theta = 0.0)
        : x(x), y(y), v(v), o(theta) {
    }

    Y operator+(const Y& other) const
    {
        return Y(x + other.x, y + other.y, v + other.v, o + other.o);
    }

    Y operator*(double scalar) const
    {
        return Y(x * scalar, y * scalar, v * scalar, o * scalar);
    }
};

const double g = 9.81;
const double C = 0.25;
const double rho = 1.29;
const double S = 0.35;
const double T = 8.0;
const double m0 = 50.0;
const double m_fuel = 15.0;
const double burn_time = 3.0;
const double fuel_burn_m = m_fuel / burn_time;

Y f(double t, const Y& y, double m)
{
    double F_soprot = (C * rho * S * y.v * y.v) / 2.0;
    double T_tyagi = (t < burn_time) ? T : 0.0;
    double V_sgorania = (t < burn_time) ? -fuel_burn_m : 0.0;

    return Y(
        y.v * cos(y.o), // x'
        y.v * sin(y.o), // y'
        (T_tyagi - F_soprot) / m - (V_sgorania / m) * y.v - g * sin(y.o), // v'
        -g * cos(y.o) / y.v // o'
    );
}

void Runge_Kutta(double a, double b, int n, std::vector<Y>& y, std::vector<double>& time)
{
    const double h = (b - a) / n;

    for (int i = 0; i < n; ++i) {
        const double t = a + i * h;

        // Обновление массы только на текущем шаге
        double m1 = (t < burn_time) ? m0 - fuel_burn_m * t : m0 - m_fuel;

        const Y k1 = f(t, y[i], m1);
        const Y k2 = f(t + (2.0 / 3.0) * h, y[i] + k1 * (2.0 / 3.0 * h), m1);

        y[i + 1] = y[i] + (k1 + k2 * 3.0) * (h / 4.0);
        time[i + 1] = t + h;

        if (y[i + 1].y < 0) break; // проверка на падение ракеты
    }
}

int main() {
    setlocale(LC_ALL, "Russian");

    const double a = 0.0;
    const double b = 100.0; // время полета с запасом
    const int n = 1000;
    std::vector<double> start_ygol = { 10, 20, 30, 40, 50, 60, 70, 80 };
    std::vector<double> Ls;

    std::cout << std::fixed << std::setprecision(10);

    for (double theta0 : start_ygol)
    {
        std::vector<Y> y(n + 1);
        std::vector<double> time(n + 1);
        y[0] = Y(0.0, 0.0, 60.0, theta0 * pi / 180.0);
        time[0] = 0.0;

        Runge_Kutta(a, b, n, y, time);

        double L = 0.0;
        int i_max = 0;
        for (int i = 0; i <= n; ++i) {
            if (y[i].y < 0) break;
            L = y[i].x;
            i_max = i;
        }
        Ls.push_back(L);

        std::cout << "\n\n// theta0 = " << theta0 << " deg\n";

        std::cout << "t = [";
        for (int i = 0; i <= i_max; i += 10)
        {
            std::cout << time[i] << ", ";
        }
        std::cout << "]\n";

        std::cout << "x = [";
        for (int i = 0; i <= i_max; i += 10) {
            std::cout << y[i].x << ", ";
        }
        std::cout << "]\n";

        std::cout << "y = [";
        for (int i = 0; i <= i_max; i += 10) {
            std::cout << y[i].y << ", ";
        }
        std::cout << "]\n";

        std::cout << "v = [";
        for (int i = 0; i <= i_max; i += 10)
        {
            std::cout << y[i].v << ", ";
        }
        std::cout << "]\n";

        std::cout << "theta = [";
        for (int i = 0; i <= i_max; i += 10)
        {
            std::cout << y[i].o * (180/pi) << ", ";
        }
        std::cout << "]\n";

        std::cout << "L = " << L << "\n";
        std::cout << "T_flight = " << time[i_max] << "\n";
    }

    std::cout << "\n\nL(θ0):\n";
    for (size_t i = 0; i < Ls.size(); ++i)
    {
        std::cout << Ls[i] << ", ";
    }

    return 0;
}
