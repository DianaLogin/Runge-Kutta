#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>

struct Y 
{
    double y1;
    double y2;

    Y(double y1 = 1.0, double y2 = 1.0) : y1(y1), y2(y2) {}

    Y operator+(const Y& other) const 
    {
        return Y(y1 + other.y1, y2 + other.y2);
    }

    Y operator*(double scalar) const
    {
        return Y(y1 * scalar, y2 * scalar);
    }
};

Y f(double t, const Y& y) 
{
    const double a = 5.0;
    const double b = 10.0;
    return Y(
        -a * y.y1 - b * y.y2 + (a + b - 1) * exp(-t),
        -b * y.y1 - a * y.y2 + (a + b - 1) * exp(-t)
            );
}

void Runge_Kutta(double a, double b, int n, std::vector<Y>& y) 
{
    const double h = (b - a) / n;

    for (int i = 0; i < n; ++i) {
        const double t = a + i * h;

        const Y k1 = f(t, y[i]);
        const Y k2 = f(t + (2.0 / 3.0) * h, y[i] + k1 * ((2.0 / 3.0) * h));

        y[i + 1] = y[i] + (k1 + k2 * 3.0) * (h / 4.0);
    }
}

int main() {
    setlocale(LC_ALL, "Russian");

    const double a = 0.0;
    const double b = 4.0;
    const std::vector<int> n_values = { 100, 200, 400, 800, 1600, 3200, 6400 };
    std::vector<double> max_errors; // Максимальные ошибки для каждого n
    std::vector<double> h_squared; // Степень шагов
    std::vector<double> err_div_h2; // Вектор для хранения err/h^2

    std::cout << std::fixed << std::setprecision(10);
    std::cout << "n\t\th\t\t\tmax_error\t\t\terror/h^2\n";
    std::cout << "======================================================================\n";

    for (int n : n_values)
    {
        double h = (b - a) / n;
        std::vector<Y> y(n + 1);
        y[0] = Y(1.0, 1.0); // Начальные условия y1(0)=1, y2(0)=1, потому что при т=0 будет e^-0 = 1
        Runge_Kutta(a, b, n, y);


        std::vector<double> err; 
        for (int i = 0; i <= n; ++i) {
            double t_i = a + i * h;
            err.push_back(fabs(y[i].y1 - exp(-t_i)) / fabs(exp(-t_i)));
        }

        // Находим максимальную ошибку для текущего n
        double max_err = (*std::max_element(err.begin(), err.end()));
        max_errors.push_back(max_err);
        h_squared.push_back(h * h);

        // Вычисление err/h^2
        double h2 = pow(h, 2);
        double err_h2 = max_err / h2;
        err_div_h2.push_back(err_h2);

        std::cout << h << ", ";
    }

    // Вывод данных для построения графиков
    std::cout << "\n\n max_errors\n";
    for (int i = 0; i < max_errors.size(); ++i)
    {
        std::cout << max_errors[i] << ", ";
    }

    std::cout << "\n\n err_over_h2\n";

    for (int i = 0; i < err_div_h2.size(); ++i)
    {
        std::cout << err_div_h2[i] << ", ";
    }


    return 0;
}