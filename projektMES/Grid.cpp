#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <cmath>

#include "Grid.h"

using namespace std;

void readDataFromFile(const string& filename, GlobalData& globalData, Grid& grid)
{
    fstream file(filename);

    if (!file.is_open())
    {
        cerr << "The file cannot be opened." << endl;
        return;
    }

    string parameterName;

    while (file >> parameterName)
    {
        // GLOBAL DATA
        if (parameterName == "SimulationTime")
            file >> globalData.simulationTime;
        else if (parameterName == "SimulationStepTime")
            file >> globalData.simulationStepTime;
        else if (parameterName == "Conductivity")
            file >> globalData.conductivity;
        else if (parameterName == "Alfa")
            file >> globalData.alfa;
        else if (parameterName == "Tot")
            file >> globalData.tot;
        else if (parameterName == "InitialTemp")
            file >> globalData.initialTemp;
        else if (parameterName == "Density")
            file >> globalData.density;
        else if (parameterName == "SpecificHeat")
            file >> globalData.specificHeat;
        // NODES NUMBER
        else if (parameterName == "Nodes")
        {
            string parameterName2;
            file >> parameterName2;
            if (parameterName2 == "number")
                file >> globalData.nodesNumber;
        }
        // ELEMENTS NUMBER
        else if (parameterName == "Elements")
        {
            string parameterName2;
            file >> parameterName2;
            if (parameterName2 == "number")
                file >> globalData.elementsNumber;
        }
        // NODES
        else if (parameterName == "*Node")
        {
            grid.nodes = new Node[globalData.nodesNumber];
            for (int i = 0; i < globalData.nodesNumber; i++)
            {
                int id;
                double x, y;
                char comma;  // ',' separator
                if (file >> id >> comma >> x >> comma >> y)
                {
                    grid.nodes[i].x = x;
                    grid.nodes[i].y = y;
                }
                else
                {
                    cerr << "Error reading Node data at line " << (i + 1) << endl;
                    return;
                }
            }
        }
        // ELEMENTS
        else if (parameterName == "*Element,")
        {
            string parameterName2;
            file >> parameterName2;

            grid.elements = new Element[globalData.elementsNumber];
            for (int i = 0; i < globalData.elementsNumber; i++)
            {
                int id;
                int first, second, third, fourth;
                char comma;  // ',' separator
                if (file >> id >> comma >> first >> comma >> second >> comma >> third >> comma >> fourth)
                {
                    grid.elements[i].id[0] = first;
                    grid.elements[i].id[1] = second;
                    grid.elements[i].id[2] = third;
                    grid.elements[i].id[3] = fourth;
                }
                else
                {
                    cerr << "Error reading Element data at line " << (i + 1) << endl;
                    return;
                }
            }
        }
        // BC
        else if (parameterName == "*BC")
        {
            for (int i = 0; i < globalData.nodesNumber; i++)
            {
                grid.nodes[i].BC = 0;
            }

            for (int i = 0; i < globalData.nodesNumber; i++)
            {
                int id;
                char comma;  // ',' separator
                if (file >> id)
                {
                    grid.nodes[id - 1].BC = 1;

                    if (!(file >> comma) || comma != ',')
                    {
                        return;
                    }
                }
                else
                {
                    cerr << "Error reading BC data at line " << (i + 1) << endl;
                    return;
                }
            }
        }
    }

    file.close();
}

void printGridData(const GlobalData& globalData, const Grid& grid)
{
    cout << "SimulationTime: " << globalData.simulationTime << endl;
    cout << "SimulationStepTime: " << globalData.simulationStepTime << endl;
    cout << "Conductivity: " << globalData.conductivity << endl;
    cout << "Alfa: " << globalData.alfa << endl;
    cout << "Tot: " << globalData.tot << endl;
    cout << "InitialTemp: " << globalData.initialTemp << endl;
    cout << "Density: " << globalData.density << endl;
    cout << "SpecificHeat: " << globalData.specificHeat << endl;
    cout << "Nodes number: " << globalData.nodesNumber << endl;
    cout << "Elements number: " << globalData.elementsNumber << endl;

    cout << fixed << setprecision(11);

    for (int i = 0; i < globalData.nodesNumber; i++)
    {
        cout << "Node " << i + 1 << ": x = " << grid.nodes[i].x << ", y = " << grid.nodes[i].y << ", BC = " << grid.nodes[i].BC << endl;
    }

    for (int i = 0; i < globalData.elementsNumber; i++)
    {
        cout << "Element " << i + 1 << ":  " << grid.elements[i].id[0] << ", " << grid.elements[i].id[1] << ", " << grid.elements[i].id[2] << ", " << grid.elements[i].id[3] << endl;
    }

    cout << defaultfloat;
}