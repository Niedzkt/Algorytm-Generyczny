#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <algorithm>
#include <fstream>

using namespace std;

class Osobnik
{
public:
    double K;
    double T;
    double Xi;
    double J;

    double odpowiedz_skokowa(double t) const
    {
        double sqrt_term = sqrt(1 - pow(Xi, 2));
        double angle = atan(sqrt_term / Xi);
        double sine_arg = (sqrt_term / T) * t;
        double sine_val = sin(angle + sine_arg);
        double exp_val = exp((-Xi / T) * t);
        double fraction = exp_val / sqrt_term;
        double result = 1 - fraction * sine_val;
        return K * result;
    }

    double odpowiedz_impulsowa(double t) const
    {
        double sqrt_term = sqrt(1 - pow(Xi, 2));
        double exp_val = exp((-Xi / T) * t);
        double impulse = (K / T * sqrt_term) * exp_val * sin((sqrt_term / T) * t);
        return impulse;
    }
};

class AlgorytmGenetyczny
{
public:
    std::vector<Osobnik> populacja;
    int rozmiar_populacji;
    int liczba_iteracji;
    double prwdp_krzyzowania;
    double prwdp_mutacji;
    double K_min, K_max;
    double T_min, T_max;
    double Xi_min, Xi_max;
    std::vector<double> h;
    std::vector<double> g;

    AlgorytmGenetyczny(int rozmiar_populacji, int liczba_iteracji, double prwdp_krzyzowania, double prwdp_mutacji,
                       double K_min, double K_max, double T_min, double T_max, double Xi_min, double Xi_max,
                       const std::vector<double> &h, const std::vector<double> &g)
    {
        this->rozmiar_populacji = rozmiar_populacji;
        this->liczba_iteracji = liczba_iteracji;
        this->prwdp_krzyzowania = prwdp_krzyzowania;
        this->prwdp_mutacji = prwdp_mutacji;
        this->K_min = K_min;
        this->K_max = K_max;
        this->T_min = T_min;
        this->T_max = T_max;
        this->Xi_min = Xi_min;
        this->Xi_max = Xi_max;
        this->h = h;
        this->g = g;
    }

    double ocena_jakosci_osobnika(const Osobnik &osobnik, int kroki_czasowe)
    {
        double J = 0;
        for (int i = 0; i < kroki_czasowe; ++i)
        {
            double t = i * 0.1;
            double model_h = osobnik.odpowiedz_skokowa(t);
            double model_g = osobnik.odpowiedz_impulsowa(t);
            J += pow(h[i] - model_h, 2) + pow(g[i] - model_g, 2);
        }
        return J;
    }

    std::vector<Osobnik> selekcja_turniejowa()
    {
        std::vector<Osobnik> wybrana_populacja;
        int rozmiar_turnieju = 3;
        for (size_t i = 0; i < populacja.size(); ++i)
        {
            Osobnik najlepszy = populacja[rand() % populacja.size()];
            for (int j = 1; j < rozmiar_turnieju; ++j)
            {
                Osobnik zawodnik = populacja[rand() % populacja.size()];
                if (zawodnik.J < najlepszy.J)
                {
                    najlepszy = zawodnik;
                }
            }
            wybrana_populacja.push_back(najlepszy);
        }
        return wybrana_populacja;
    }

    std::vector<Osobnik> krzyzowanie(const std::vector<Osobnik> &populacja)
    {
        std::vector<Osobnik> nowa_populacja;
        for (size_t i = 0; i < populacja.size() - 1; i += 2)
        {
            Osobnik rodzic1 = populacja[i];
            Osobnik rodzic2 = populacja[i + 1];
            if ((rand() / double(RAND_MAX)) < prwdp_krzyzowania)
            {
                std::swap(rodzic1.T, rodzic2.T);
                std::swap(rodzic1.Xi, rodzic2.Xi);
            }
            nowa_populacja.push_back(rodzic1);
            nowa_populacja.push_back(rodzic2);
        }
        if (populacja.size() % 2 != 0)
        {
            nowa_populacja.push_back(populacja.back());
        }
        return nowa_populacja;
    }

    void mutacja(std::vector<Osobnik> &populacja)
    {
        for (auto &osobnik : populacja)
        {
            if ((rand() / double(RAND_MAX)) < prwdp_mutacji)
            {
                osobnik.K = K_min + (K_max - K_min) * (rand() / double(RAND_MAX));
            }
            if ((rand() / double(RAND_MAX)) < prwdp_mutacji)
            {
                osobnik.T = T_min + (T_max - T_min) * (rand() / double(RAND_MAX));
            }
            if ((rand() / double(RAND_MAX)) < prwdp_mutacji)
            {
                osobnik.Xi = Xi_min + (Xi_max - Xi_min) * (rand() / double(RAND_MAX));
            }
        }
    }

    void algorytm_genetyczny()
    {
        std::ofstream file("output.csv", std::ios::app); // Open the file for appending

        for (int i = 0; i < rozmiar_populacji; ++i)
        {
            Osobnik osobnik;
            osobnik.K = K_min + (K_max - K_min) * (rand() / double(RAND_MAX));
            osobnik.T = T_min + (T_max - T_min) * (rand() / double(RAND_MAX));
            osobnik.Xi = Xi_min + (Xi_max - Xi_min) * (rand() / double(RAND_MAX));
            osobnik.J = ocena_jakosci_osobnika(osobnik, liczba_iteracji);
            populacja.push_back(osobnik);
        }
        Osobnik najlepszy_osobnik = *std::min_element(populacja.begin(), populacja.end(), [](const Osobnik &a, const Osobnik &b)
                                                      { return a.J < b.J; });
        for (int iter = 0; iter < liczba_iteracji; ++iter)
        {
            std::vector<Osobnik> wybrana_populacja = selekcja_turniejowa();
            std::vector<Osobnik> nowa_populacja = krzyzowanie(wybrana_populacja);
            mutacja(nowa_populacja);
            for (auto &osobnik : nowa_populacja)
            {
                osobnik.J = ocena_jakosci_osobnika(osobnik, liczba_iteracji);
            }
            populacja = nowa_populacja;
            Osobnik najlepszy_nowy = *std::min_element(populacja.begin(), populacja.end(), [](const Osobnik &a, const Osobnik &b)
                                                       { return a.J < b.J; });
            if (najlepszy_nowy.J < najlepszy_osobnik.J)
            {
                najlepszy_osobnik = najlepszy_nowy;
            }
        }
        std::cout << "\nNajlepszy osobnik:\n"
                  << najlepszy_osobnik.K << "\n"
                  << najlepszy_osobnik.T << "\n"
                  << najlepszy_osobnik.Xi << "\n"
                  << najlepszy_osobnik.J << std::endl;

        // Write the best individual to the file
        file << najlepszy_osobnik.K << "," << najlepszy_osobnik.T << "," << najlepszy_osobnik.Xi << "," << najlepszy_osobnik.J << std::endl;

        file.close(); // Close the file
    }
};

int main()
{
    srand(time(0));
    int rozmiar_populacji = 70;
    int liczba_iteracji;
    double prwdp_krzyzowania = 0.5;
    double prwdp_mutacji =0.5;
    double K_min, K_max;
    double T_min, T_max;
    double Xi_min, Xi_max;
    std::cout << "Podaj liczbe iteracji: ";
    std::cin >> liczba_iteracji;
    std::cout << "Podaj rozmiar populacji: ";
    std::cin >> rozmiar_populacji;
    std::cout << "Podaj prawdopodobienstwo krzyzowania (0-1): ";
    std::cin >> prwdp_krzyzowania;
    std::cout << "Podaj prawdopodobienstwo mutacji (0-1): ";
    std::cin >> prwdp_mutacji;
    std::cout << "Podaj K_min: ";
    std::cin >> K_min;
    std::cout << "Podaj K_max: ";
    std::cin >> K_max;
    std::cout << "Podaj T_min: ";
    std::cin >> T_min;
    std::cout << "Podaj T_max: ";
    std::cin >> T_max;
    do
    {
        std::cout << "Podaj Xi_min (0-1): ";
        std::cin >> Xi_min;
        if (Xi_min < 0 || Xi_min > 1)
        {
            std::cout << "Xi_min musi byc w przedziale od 0 do 1." << std::endl;
        }
    } while (Xi_min < 0 || Xi_min > 1);

    do
    {
        std::cout << "Podaj Xi_max (0-1): ";
        std::cin >> Xi_max;
        if (Xi_max < 0 || Xi_max > 1)
        {
            std::cout << "Xi_max musi byc w przedziale od 0 do 1." << std::endl;
        }
    } while (Xi_max < 0 || Xi_max > 1);

    std::vector<double> h(liczba_iteracji, 1.0);
    std::vector<double> g(liczba_iteracji, 1.0);
    AlgorytmGenetyczny algorytm(rozmiar_populacji, liczba_iteracji, prwdp_krzyzowania, prwdp_mutacji,
                                K_min, K_max, T_min, T_max, Xi_min, Xi_max, h, g);
    std::cout << "Rozpoczynam algorytm genetyczny..." << std::endl;
    for(int i = 0; i < 10; ++i)
    {
        std::cout << "iteracja:" << i << std::endl;
        algorytm.algorytm_genetyczny();
    }
    return 0;
}