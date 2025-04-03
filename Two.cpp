#include <iostream>
#include <cmath>
#include <vector>
#include <locale>
#include <algorithm>

// ��������� ��� �������� ���������� ������� ������ ���������� ��������
struct �
{
    double y1;
    double y2;

    �(double y1 = 0.0, double y2 = 0.0) : y1(y1), y2(y2) {}


    � operator+(const �& other) const
    {
        return �(y1 + other.y1, y2 + other.y2);
    }

    � operator*(double scalar) const
    {
        return �(y1 * scalar, y2 * scalar);
    }
};

// ���������������� ��������� ������� dy/dt
� f(double t, const �& y)
{
    double a = 5.0, b = 10.0;

    return �(
        -a * y.y1 - b * y.y2 + (a + b - 1) * exp(-t), // y'1
        b * y.y1 - a * y.y2 + (a + b - 1) * exp(-t)   // y'2
    );
}

// ����� �����-����� 2-�� �������
void Runge_Kutta(double a, double b, int n, std::vector<�>& Yt)
{
    double h = (b - a) / n;
    Yt.resize(n + 1);

    for (int i = 0; i < n; ++i)
    {
        double t = a + i * h;

        � k1 = f(t, Yt[i]);
        � k2 = f(t + (2.0 / 3.0) * h, Yt[i] + k1 * (2.0 / 3.0 * h));

        Yt[i + 1] = Yt[i] + (k1 + k2 * 3.0) * (h / 4.0);
    }
}

int main() {
    setlocale(LC_ALL, "Russian");

    // ��������� ��������������
    double a = 0.0, b = 4.0;
    std::vector<int> n_values = { 10, 20, 40, 80, 160, 320 };

    std::cout << "n\th\tmax_error\terror/h^2\n";

    for (int n : n_values) {
        double h = (b - a) / n;
        std::vector<�> y(n + 1);

        // ��������� �������
        y[0] = �(1.0, 1.0);

        // ��������������
        Runge_Kutta(a, b, n, y);

        // ���������� ������
        double max_error = 0.0;
        for (int i = 0; i <= n; ++i) {
            double t = a + i * h;

            double error1 = fabs(y[i].y1 - exp(-t));
            double error2 = fabs(y[i].y2 - exp(-t));

            max_error = std::max({ max_error, error1, error2 });
        }

        std::cout << n << "\t" << h << "\t" << max_error
            << "\t" << max_error / (h * h) << "\n";
    }

    return 0;
}