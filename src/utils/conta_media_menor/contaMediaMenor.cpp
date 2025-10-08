#include <iostream>
#include <fstream>
#include <vector>
#include <string>

// Função para calcular a média de um conjunto de números
double calcularMedia(const std::vector<double>& conjunto) {
    double soma = 0.0;
    for (double numero : conjunto) {
        soma += numero;
    }
    return soma / conjunto.size();
}

int main() {
    int contadorConjunto1Menor = 0;
    int contadorConjunto2Menor = 0;
    std::string ins[] = {"A","B", "C", "D", "E", "F", "G", "H", "I" ,"J"};
    for (int i = 0; i < 10; ++i) {
        std::string moead_s = "..\\Analises\\"+ins[i]+"\\esp_ad_moead.out";
        std::string nsga_s = "..\\Analises\\"+ins[i]+"\\esp_ad_nsga2.out";
        std::ifstream moead(moead_s);
        std::ifstream nsga(nsga_s);

        

        std::vector<double> conjunto1;
        std::vector<double> conjunto2;
        std::string linha;

        while (std::getline(moead, linha)) {
            if (linha.empty()==false) {
                double numero = std::stod(linha);
                conjunto1.push_back(numero);
           //     std::cout<<numero<<std::endl;
            }
        }
       // std::cout<<std::endl;
         while (std::getline(nsga, linha)) {
            if (linha.empty()==false) {
                double numero = std::stod(linha);
                conjunto2.push_back(numero);
         //       std::cout<<numero<<std::endl;
            }
        }


        moead.close();
        nsga.close();
       // std::cout<<std::endl;
        double mediaConjunto1 = calcularMedia(conjunto1);
        double mediaConjunto2 = calcularMedia(conjunto2);

        //std::cout<<mediaConjunto1<<std::endl;
        //std::cout<<mediaConjunto2<<std::endl;
        if (mediaConjunto1 < mediaConjunto2) {

            std::cout << "Arquivo "<<i<<" MOEA"<<std::endl;
            contadorConjunto1Menor++;
        } else if (mediaConjunto2 < mediaConjunto1) {
            contadorConjunto2Menor++;

            std::cout << "Arquivo "<<i<<" NSGA"<<std::endl;
        }
    }

    std::cout << "Quantidade de arquivos onde o conjunto 1 teve média menor que o conjunto 2: " << contadorConjunto1Menor << std::endl;
    std::cout << "Quantidade de arquivos onde o conjunto 2 teve média menor que o conjunto 1: " << contadorConjunto2Menor << std::endl;

    return 0;
}