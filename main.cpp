#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <iomanip>

using namespace std;

struct Individual{
    vector<float> coefficients;
    double fitness = 0.0;
};

const int POP_SIZE = 100, MAX_GENERATION = 1000;
const float PM = 0.1, PC = 0.6;
vector<Individual> currentGeneration, newGeneration;
vector<pair<float, float>> points;
int pointsNum, D, currGen = 0, k=2;
float b = 2;

bool sortByFitness(Individual& a, Individual& b){
    return a.fitness > b.fitness;
}

vector<Individual> InitializePopulation(){
    vector<Individual> Generation;
    while(Generation.size() < POP_SIZE){
        Individual v;
        v.coefficients.resize(D + 1);
        for (int i = 0; i < D+1; i++)
        {
            // min + ( std::rand() % ( max - min + 1 ) )
            float coff = -10 + (rand() % ( 20 + 1));
            v.coefficients[i] = coff;
        }
        Generation.push_back(v);
    }
    return Generation;
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
                double xp = pow(points[j].first, k);
                predictedY += (currentGeneration[i].coefficients[k] * xp);
            }
            double error = pow(predictedY - points[j].second, 2)  ;
            SE += error;
        }
        MSE = SE / pointsNum;
        currentGeneration[i].fitness = 1.00 / MSE ;
    }

}

Individual tournamentSelection(){
    int *mating_pool = new int[k]; // Indices for the individual participating in the tournament
    for(int i = 0; i < k; i++)
        mating_pool[i] = 0 + (rand() % currentGeneration.size()); // Randomize the indices
    double max_score = -1; // To find the highest score among contestants
    int max_score_index; // To find the index of the contestant with the highest score
    for(int i = 0; i < k; i++)
    {
        if(max_score < currentGeneration[mating_pool[i]].fitness) {
            max_score = currentGeneration[mating_pool[i]].fitness;
            max_score_index = mating_pool[i];
        }
    }
    return currentGeneration[max_score_index];  // Return contestant with the highest score
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
    for (int i = 0; i < newGeneration.size(); ++i) {
        float mn = *min_element(newGeneration[i].coefficients.begin(), newGeneration[i].coefficients.end());
        float mx = *max_element(newGeneration[i].coefficients.begin(), newGeneration[i].coefficients.end());
        for (int j = 0; j < D + 1; ++j) {
            float r = (float)rand() / (float)RAND_MAX;
            if(r <= PM){
                float dl = newGeneration[i].coefficients[j] - mn;
                float du = mx - newGeneration[i].coefficients[j];
                float xi = newGeneration[i].coefficients[j];
                float r1 = (float)rand() / (float)RAND_MAX;
                float y;
                if(r1 <= 0.5){
                    y = dl;
                }else{
                    y = du;
                }
                float r2 = (float)rand() / (float)RAND_MAX;
                float p = pow (1 - currGen / MAX_GENERATION, b);
                float d = y * (1 - pow(r2, p));
                if(y == dl){
                    newGeneration[i].coefficients[j] = xi - d;
                }else if(y == du){
                    newGeneration[i].coefficients[j] = xi + d;
                }
            }
        }
    }

}

void Replacement(){
    newGeneration.clear();
    double percent = (10.00 / 100.00);
    int best = ceil((double)percent * POP_SIZE);
    while(newGeneration.size() < (POP_SIZE - best)){
        Individual parent1 = tournamentSelection();
        Individual parent2 = tournamentSelection();
        float r = (float)rand() / (float)RAND_MAX;
        if(r <= PC){
            // perform crossover
            pair<Individual, Individual> springs = crossover(parent1, parent2);
            newGeneration.push_back(springs.first);
            if(newGeneration.size() < (POP_SIZE - best))  newGeneration.push_back(springs.second);
        }else{
            newGeneration.push_back(parent1);
            if(newGeneration.size() < (POP_SIZE - best)) newGeneration.push_back(parent2);
        }
    }
    Mutation();
    std::sort(currentGeneration.begin(), currentGeneration.end(), sortByFitness);
    for (int i = 0; i < best; i++) {
        newGeneration.push_back(currentGeneration[i]);
    }
}

Individual GA(){
    currentGeneration.clear();
    currentGeneration.resize(POP_SIZE);
    currentGeneration = InitializePopulation();
    for (int i = 0; i < MAX_GENERATION; i++)
    {
        currGen = i;
        evaluateFitness();
        Replacement();
        currentGeneration = newGeneration;
    }
    std::sort(currentGeneration.begin(), currentGeneration.end(), sortByFitness);
    return currentGeneration[0];
}

int main(){
    freopen("curve_fitting_input.txt", "r", stdin);
    freopen("curve_fitting_output.txt", "w", stdout);
    int testcases;
    cin >> testcases;
    for(int i=0; i < testcases; i++){
        cin >> pointsNum >> D;
        points.clear();
        points.resize(pointsNum);
        for (int j = 0; j < pointsNum; j++)
        {
            cin >> points[j].first >> points[j].second;
        }
        Individual bestIndividual = GA();
        cout << "------- test case " << i+1 << " -------\n";
        cout << "coefficients of the polynomial function of degree " << D << " : \n";
        cout << "{ ";
        for (int j = 0; j < bestIndividual.coefficients.size(); ++j) {
            if(j == bestIndividual.coefficients.size() - 1) {
                cout << fixed << setprecision(6) <<  bestIndividual.coefficients[j] << " }\n";
            }else{
                cout << fixed <<  setprecision(6) << bestIndividual.coefficients[j] << " , ";
            }
        }
        cout << "Fitness: " << fixed << setprecision(6) << bestIndividual.fitness << "\n";
        cout << "MSE: " << fixed << setprecision(6) << 1.00 / bestIndividual.fitness << "\n";
    }
}