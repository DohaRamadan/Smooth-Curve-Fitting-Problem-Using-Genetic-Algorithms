#include <iostream>
#include <vector>
#include <cmath>
#define x first
#define y second
using namespace std;

struct Individual{
    vector<float> coefficients;
    double fitness = 0.0;
};

int pointsNum, D;
const int POP_SIZE = 4, MAX_GENERATION = 1000;
vector<Individual> currentGeneration, offsprings;
vector<pair<int, int>> points;

vector<Individual> InitializePopulation(){
    vector<Individual> newGeneration;
    while(newGeneration.size() < POP_SIZE){
        Individual v;
        v.coefficients.resize(D + 1);
        for (int i = 0; i < D+1; i++)
        {
            // min + ( std::rand() % ( max - min + 1 ) )
            float coff = -10 + (rand() % ( 20 + 1));
            v.coefficients[i] = coff;
        }
        newGeneration.push_back(v);
    }
    return newGeneration;
}


void evaluateFitness(){
    for (int i = 0; i < POP_SIZE; i++)
    {
        double SE = 0.0, MSE = 0.0;
        // J is the point (x, y)
        for(int j=0; j<pointsNum; j++){
            double predictedY = 0.0;
            // K is the coeffeints [0.1, 2.5, 10.5]
            for(int k=0; k<D+1; k++){
                double xp = pow(points[j].x, k);
                predictedY += (currentGeneration[i].coefficients[k] * xp);
            }
            double error = pow(predictedY - points[j].y, 2)  ;
            SE += error;
        }
        MSE = SE / pointsNum;
        currentGeneration[i].fitness = 1.00 / MSE ;
    }

}

pair<Individual, Individual> crossover(Individual parent1, Individual parent2){
    int r1 = 0 + (rand() % (D - 1));
    int r2 = r1+1 + (rand() % (D - 1));
    Individual child1, child2;
    child1.coefficients.resize(D+1);
    child2.coefficients.resize(D+1);
    for (int i = 0; i <= min(r1, r2); i++)
    {
        child1.coefficients[i] = parent1.coefficients[i];
        child2.coefficients[i] = parent2.coefficients[i];
    }

    for (int i = min(r1, r2)+1; i <= max(r1, r2); i++)
    {
        child1.coefficients[i] = parent2.coefficients[i];
        child2.coefficients[i] = parent1.coefficients[i];
    }
    
    for (int i = max(r1, r2) + 1; i < D + 1; i++)
    {
        child1.coefficients[i] = parent1.coefficients[i];
        child2.coefficients[i] = parent2.coefficients[i];
    }

    return {child1, child2};
}


Individual GA(){
    currentGeneration.clear();
    currentGeneration.resize(POP_SIZE);
    currentGeneration = InitializePopulation();
    for (int i = 0; i < MAX_GENERATION; i++)
    {
        evaluateFitness();
    }
}

int main(){
    freopen("curve_fitting_input.txt", "r", stdin);
    freopen("curve_fitting_output.txt", "w", stdout);
    int t;
    cin >> t;
    while(t--){
        cin >> pointsNum >> D;
        points.resize(pointsNum);
        for (int i = 0; i < pointsNum; i++)
        {
            cin >> points[i].x >> points[i].y;
        }
        GA();
    }
}