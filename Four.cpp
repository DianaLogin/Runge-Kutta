#include <iostream>
#include <cmath>
#include <vector>
#include <locale>
#include <algorithm>

// Дифференциальное уравнение dy/dx = f(x, y)
double f(double x, double y) 
{
    return y; // Пример: dy/dx = y
}

// Метод Рунге-Кутты 4-го порядка 
void Runge_Kutta(const double& a, const double& b, const int& n, std::vector<double>& y)
{
    double h = (b - a) / n; // шаг по иксам


    double m_1, m_2, m_3, m_4;

    for (int i = 0; i < n; ++i)
    {
        double x_i = a + i * h;

         m_1 = f(x_i, y[i]);
         m_2 = f(x_i + 0.5 * h, y[i] + 0.5 * h * m_1);
         m_3 = f(x_i + 0.5 * h, y[i] + 0.5 * h * m_2);
         m_4 = f(x_i + h, y[i] + h * m_3);

        y[i + 1] = y[i] + (h / 6) * (m_1 + 2 * m_2 + 2 * m_3 + m_4);
    }
}


int mainn()
{
    setlocale(LC_ALL, "Russian");
    double a = 0.0;
    double b = 1.0;
    std::vector<int> n_values = { 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024 }; 
    std::vector<double> max_errors; // Максимальные ошибки для каждого n
    std::vector<double> h_squared; // Степень шагов
    std::vector<double> err_over_h4; // Вектор для хранения err/h^4
   
    for (int n : n_values)
    {
        double h = (b - a) / n;
        std::vector<double> y(n + 1);
        y[0] = 1.0; // Начальное значение 
        Runge_Kutta(a, b, n, y);

        // Вычисление погрешности err = max |e^x - rk| \ |e^x|
        std::vector<double> err;
        for (int i = 0; i <= n; ++i) 
        {
            double x_i = a + i * h;
            err.push_back(fabs((exp(x_i) - y[i]) / fabs(exp(x_i))));
        }

        // Находим максимальную ошибку для текущего n
        double max_err = (*std::max_element(err.begin(), err.end()));
        max_errors.push_back(max_err);
        h_squared.push_back(h * h * h * h);

        // Вычисление err/h^4
        double h4 = pow(h, 4);
        double err_h4 = max_err / h4;
        err_over_h4.push_back(err_h4);
    }

    for (int i = 0; i < max_errors.size(); ++i)
    {
        std::cout << max_errors[i] << ", ";
    }

    std::cout  << "\n\n";

    for (int i = 0; i < err_over_h4.size(); ++i)
    {
        std::cout << err_over_h4[i] << ", ";
    }

    return 0;
}


//int main()
//{
//   
//    std::cout << "Приблизительные значения:\n";
//    for (int i = 0; i <= n; ++i)
//    {
//        double x_i = a + i * h;
//        std::cout << "y(" << x_i << ") = " << y[i] << "\n";
//    }
//
//    std::vector<double> err;
//    std::cout << "Точные значения:\n";
//    for (int i = 0; i <= n; ++i)
//    {
//        double x_i = a + i * h;
//        std::cout << "exp(" << x_i << ") = " << exp(x_i) << "\n";
//    }
//
//    return 0;
//}
