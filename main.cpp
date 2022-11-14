#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#define x first
#define y second
using namespace std;

struct Individual{
    vector<float> coefficients;
    double fitness = 0.0;
};

int pointsNum, D, t;
float b = 2;
const int POP_SIZE = 4, MAX_GENERATION = 1000, PM = 0.1, PC = 0.6;
vector<Individual> currentGeneration, offsprings;
vector<pair<int, int>> points;

vector<Individual> InitializePopulation(){
    vector<Individual> newGeneration;
    /*while(newGeneration.size() < POP_SIZE){
        Individual v;
        v.coefficients.resize(D + 1);
        for (int i = 0; i < D+1; i++)
        {
            // min + ( std::rand() % ( max - min + 1 ) )
            float coff = -10 + (rand() % ( 20 + 1));
            v.coefficients[i] = coff;
        }
        newGeneration.push_back(v);
    }*/
    Individual vect1;
    vect1.coefficients = {1.95, 8.16, -2};
    Individual vect2;
    vect2.coefficients = {4.26, -7.4, -2.5};
    Individual vect3;
    vect3.coefficients = {3.36, -0.3, -6.2};
    Individual vect4;
    vect4.coefficients = {0.23, 0.12, 4.62};

    newGeneration.push_back(vect1);
    newGeneration.push_back(vect2);
    newGeneration.push_back(vect3);
    newGeneration.push_back(vect4);

    return newGeneration;
}

Individual tournamentSelection(vector<Individual> population, int k){
    int mating_pool[k]; // Indices for the individual participating in the tournament
    for(int i = 0; i < k; i++)
        mating_pool[i] = 0 + (rand() % population.size()); // Randomize the indices
    double max_score = -1; // To find the highest score among contestants
    int max_score_index; // To find the index of the contestant with the highest score
    for(int i = 0; i < k; i++)
    {
        if(max_score < population[mating_pool[i]].fitness) {
            max_score = population[mating_pool[i]].fitness;
            max_score_index = mating_pool[i];
        }
    }
    return population[max_score_index];  // Return contestant with the highest score
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

void Mutation(){
    for (int i = 0; i < POP_SIZE; ++i) {
        float mn = *min_element(offsprings[i].coefficients.begin(), offsprings[i].coefficients.end());
        float mx = *max_element(offsprings[i].coefficients.begin(), offsprings[i].coefficients.end());
        for (int j = 0; j < D + 1; ++j) {
            float r = (float)rand() / (float)RAND_MAX;
            if(r <= PM){
                float dl = offsprings[i].coefficients[j] - mn;
                float du = mx - offsprings[i].coefficients[j];
                float xi = offsprings[i].coefficients[j];
                float r1 = (float)rand() / (float)RAND_MAX;
                float y;
                if(r1 <= 0.5){
                    y = dl;
                }else{
                    y = du;
                }
                float p = pow (1 - t / MAX_GENERATION, b);
                float d = y * (1 - pow(r, p));
                if(y == dl){
                    offsprings[i].coefficients[j] = xi - d;
                }else if(y == du){
                    offsprings[i].coefficients[j] = xi + d;
                }
            }
        }
    }

}

Individual GA(){
    currentGeneration.clear();
    currentGeneration.resize(POP_SIZE);
    currentGeneration = InitializePopulation();
    for (int i = 0; i < MAX_GENERATION; i++)
    {
        t = i;
        evaluateFitness();
    }
}

int main(){
    //freopen("curve_fitting_input.txt", "r", stdin);
    //freopen("curve_fitting_output.txt", "w", stdout);
    int t;
    cin >> t;
    while(t--){
        cin >> pointsNum >> D;
        points.resize(pointsNum);
        for (int i = 0; i < pointsNum; i++)
        {
            cin >> points[i].x >> points[i].y;
        }
        // GA();
        vector<Individual> population = InitializePopulation();
        currentGeneration = population;
        evaluateFitness();
        Individual max = tournamentSelection(currentGeneration, 3);
    }


}