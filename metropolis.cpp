#include <iostream>
#include <list>
#include <fstream>
#include <sstream>
#include <map>
#include <bitset>
#include<cmath>

using namespace std;

#include "data.h"

// Create a bitset for integers up to 2**128
bitset<n> bitset128(__int128_t Op)
{
    bitset<n> hi{ static_cast<unsigned long long>(Op >> 64) },
        lo{ static_cast<unsigned long long>(Op) },
        bits{ (hi << 64) | lo }; // Create binary representation

    return bits;
}

map<unsigned int, list<Interaction>> write_interactions_metropolis(double J, string file)
{
    string line, line2;
    unsigned int node; // First node involved with an interaction
    map<unsigned int, list<Interaction>> I_list; // Iterator for list of interactions
    Interaction I;

    __int128_t Op = 1;
    Op <<= (n - 1);

    bool flag = false;

    // Open file to read network
    ifstream myfile(file.c_str());
    if (myfile.is_open())
    {
        while (getline(myfile, line))
        {
            stringstream ss(line);
            // Create interaction
            I.Op = 0;
            I.g = J;
            // Read nodes involved in interaction
            while (getline(ss, line2, '\t'))
            {
                node = stoi(line2) - 1;
                I.Op += (Op >> node);
            }
            // Add interaction to list corresponding to 2nd node
            I_list[node].push_back(I);

            //bitset<n> bits = bitset128(I.Op);
            //cout << bits << endl;
        }
        myfile.close();
    }
    return I_list;
}

double delta_energy(__int128_t state, list<Interaction> edges, int spin)
{
    //signed int sign = (signed)bitset128((1 << spin) & state).count() * 2 - 1; // Get sign of spin;

    __int128_t Op;
    signed int sisj;
    double diff_E = 0;

    for (auto& edge : edges)
    {
        // Apply operator to state
        Op = state & edge.Op;
        // Calculate s_i * s_j
        bitset<n> bits = bitset128(Op);
        sisj = 2 * ( (1 + (signed)bits.count()) % 2 ) - 1;
        // Add to sum
        diff_E += 2.0 * edge.g * (double)sisj; // Difference in energy
    }
    return diff_E;
}

void sample_data_metropolis(double J, string input_file, string output_filename, unsigned int N = 1000)
{
    map<unsigned int, list<Interaction>> IA;
    list<Interaction> edges;
    int spin;
    __int128_t Op, state = 0, one = 1;
    double eps, energy, diff_E;

    // INTERACTIONS:
    IA = write_interactions_metropolis(J, input_file);

    // INITIAL ENERGY
    for (auto& elem_map : IA)
    {
        for (auto& edge : elem_map.second)
        {
            // Apply operator to state
            Op = state & edge.Op;
            // Create binary representation
            bitset<n> bits = bitset128(Op);
            // Add to energy
            energy += edge.g * (2 * ((1 + (signed)bits.count()) % 2) - 1);
        }
    }

    // WINDUP-PERIOD
    for (int interval = 0; interval < 100 * n; interval++)
    {
        // Select a node
        spin = rand() % n;
        // Draw uniform number
        eps = (double)rand() / RAND_MAX;
        // Select interactions that involve the selected node
        edges = IA[spin];
        // Calculate the difference in energy
        diff_E = delta_energy(state, edges, spin);
        // Accept/reject new state
        if (diff_E < 0) { state ^= (one << (n - 1 - spin)); state = ~state; }
        else if (eps < exp(-diff_E)) { state ^= (one << (n - 1 - spin)); state = ~state; }
    }

    // OUTPUT FILE:
    fstream file(output_filename.c_str(), ios::out);

    for (int i = 0; i < N; i++)
    {
        for (int interval = 0; interval < 100 * n; interval++)
        {
            // Select a node
            spin = rand() % n;
            // Draw uniform number
            eps = (double)rand() / RAND_MAX;
            // Select interactions that involve the selected node
            edges = IA[spin];
            // Calculate the difference in energy
            diff_E = delta_energy(state, edges, spin);
            // Accept/reject new state
            if (diff_E < 0) { state ^= (one << (n - 1 - spin)); state = ~state; }
            else if (eps < exp(-diff_E)) { state ^= (one << (n - 1 - spin)); state = ~state; }
        }
        bitset<n> bits = bitset128(state);
        
        file << bits << endl;
    }

    file.close();
}