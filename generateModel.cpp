
#include "generateModel.h"
#include "parameters.h"
#include "constants.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

using namespace std;

parameters readInputFile(string );
void generateKnapsackModel(string );
void generateAssignmentModel(string );

/*================================
 Global input file variables
 ================================*/
int solverType;
int pointLimit;
double timeLimit;
double boundTolerance;
double territoryLength;
int numofObjective;
char inputType;
string problemType;
int numofKnapsack;
int numofItem;
int numofJob;
char runType;
string fileName;
int numofInstance;

float EPSILON = OBJ_EPSILON;


parameters generateModel(string path)
{
    parameters params = readInputFile(path);
    
    if (inputType=='M')
    {
        
    }
    else if (problemType=="KP")
    {
        generateKnapsackModel(path);
    }
    else if (problemType=="AP")
    {
        generateAssignmentModel(path);
    }
    
    return(params);
}

parameters readInputFile(string path)
{
    ifstream InputFile;
    string mainFile = path + "MainFile.txt";
    InputFile.open(mainFile.c_str());
    
    InputFile >> solverType;
    InputFile >> numofObjective;
    InputFile >> timeLimit;
    InputFile >> pointLimit;
    if (pointLimit == ALL_ND_INDICATOR) {
        pointLimit = POINT_LIMIT;
    }
    InputFile >> boundTolerance;
    InputFile >> territoryLength;
    InputFile >> inputType;
    InputFile >> problemType;
    InputFile >> numofKnapsack;
    InputFile >> numofItem;
    InputFile >> numofJob;
    InputFile >> fileName;
    InputFile.close();
    
    return parameters(solverType, numofObjective, timeLimit, pointLimit, boundTolerance, territoryLength);
}

void generateKnapsackModel(string path)
{
    ifstream InputFile;
    ofstream OutputFile;
    InputFile.open((path+fileName).c_str());
    OutputFile.open((path+"model.lp").c_str());
    
    int i,j,k;
    vector<vector<int>> objective;
    vector<vector<int>> weight;
    vector<float> capacity;
    vector<int> dummy;
    for(i=1;i<=numofItem;i++) dummy.push_back(0);
    
    
    /*================================================
     Reading input parameter values
     ==============================================*/
    
    
    // take objective function coefficients
    
    for(j=0;j<numofObjective;j++)
    {
        objective.push_back(dummy);
        for(i=0;i<numofItem;i++)
        {
            InputFile >> objective[j][i];
        }
    }
    
    // take weight coefficients
    
    for(k=0;k<numofKnapsack;k++)
    {
        weight.push_back(dummy);
        for(i=0;i<numofItem;i++)
        {
            InputFile >> weight[k][i];
        }
        capacity.push_back(0);
        InputFile >> capacity[k];
    }
    
    /*================================================
     Creating the model file
     ==============================================*/
    
    OutputFile << "MAX" << endl;
    if(solverType==2 || solverType==3){ // SBA or TDA solver
        for(j=0;j<numofObjective;j++)
        {
            OutputFile << EPSILON << "z" << j+1 << " + ";
        }
        OutputFile << "z" << numofObjective+1 << endl;
    } else{
        for(j=0;j<numofObjective-1;j++)
        {
            OutputFile << EPSILON << "z" << j+1 << " + ";
        }
        OutputFile << "z" << numofObjective << endl;
    }
    
    OutputFile << "SUBJECT TO" << endl;
    
    for(k=0;k<numofKnapsack;k++)
    {
        for(i=0;i<numofItem-1;i++)
        {
            if (weight[k][i]<0) OutputFile << weight[k][i] << "x" << i+1 ;
            else OutputFile << "+" << weight[k][i] << "x" << i+1 ;
        }
        
        if (weight[k][numofItem-1]<0) OutputFile << weight[k][numofItem-1] << "x" << numofItem << " < " << capacity[k] << endl;
        else OutputFile << "+" << weight[k][numofItem-1] << "x" << numofItem << " < " << capacity[k] << endl;
    }
    
    for(j=0;j<numofObjective;j++)
    {
        for(i=0;i<numofItem-1;i++)
        {
            if (objective[j][i]<0) OutputFile << objective[j][i] << "x" << i+1 ;
            else OutputFile << "+" << objective[j][i] << "x" << i+1 ;
        }
        
        if (objective[j][numofItem-1]<0) OutputFile << objective[j][numofItem-1] << "x" << numofItem << " - " << "z" << j+1 << " = " << "0" << endl;
        else OutputFile << "+" <<objective[j][numofItem-1] << "x" << numofItem << " - " << "z" << j+1 << " = " << "0" << endl;
    }
    
    OutputFile << "BOUNDS" << endl;
    
    for(j=0;j<numofObjective;j++)
    {
        OutputFile << "z" << j+1 << " free" <<endl;
    }
    
    OutputFile << "BINARIES" << endl;
    
    for(i=0;i<numofItem;i++)
    {
        OutputFile << "x" << i+1 << endl;
    }
    
    OutputFile << "END";
    
    OutputFile.close();
    InputFile.close();
}

void generateAssignmentModel(string path)
{
    ifstream InputFile;
    ofstream OutputFile;
    InputFile.open((path+fileName).c_str());
    OutputFile.open((path+"model.lp").c_str());
    
    int i,j,k;
    vector<vector<vector<int>>> objective;
    vector<int> dummy;
    vector<vector<int>> dummy2;
    for(i=1;i<=numofJob;i++) dummy.push_back(0);
    for(i=1;i<=numofJob;i++){ for(k=1;k<=numofJob;k++) dummy2.push_back(dummy);}
    
    
    /*================================================
     Reading input parameter values
     ==============================================*/
    for(j=0;j<numofObjective;j++)
    {
        objective.push_back(dummy2);
        for(i=0;i<numofJob;i++)
        {
            for(k=0;k<numofJob;k++)
            {
                InputFile >> objective[j][i][k];
            }
            
        }
    }
    
    /*================================================
     Creating the model file
     ==============================================*/
    
    OutputFile << "MAX" << endl;
    if(solverType==2 || solverType==3){ // SBA or TDA solver
        for(j=0;j<numofObjective;j++)
        {
            OutputFile << EPSILON << "z" << j+1 << " + ";
        }
        OutputFile << "z" << numofObjective+1 << endl;
    } else{
        for(j=0;j<numofObjective-1;j++)
        {
            OutputFile << EPSILON << "z" << j+1 << " + ";
        }
        OutputFile << "z" << numofObjective << endl;
    }
    
    OutputFile << "SUBJECT TO" << endl;
    
    for(i=0;i<numofJob;i++)
    {
        for(k=0;k<numofJob-1;k++)
        {
            OutputFile << "x" << i+1 << "_" << k+1 << " + ";
        }
        
        OutputFile << "x" << i+1 << "_" << numofJob << " = " << "1" << endl;
    }
    
    for(k=0;k<numofJob;k++)
    {
        for(i=0;i<numofJob-1;i++)
        {
            OutputFile << "x" << i+1 << "_" << k+1 << " + ";
        }
        
        OutputFile << "x" << numofJob << "_" << k+1 << " = " << "1" << endl;
    }
    
    for(j=0;j<numofObjective;j++)
    {
        for(i=0;i<numofJob-1;i++)
        {
            for(k=0;k<numofJob;k++)
            {
                if (objective[j][i][k]<0) OutputFile << objective[j][i][k] << "x" << i+1 << "_" << k+1 ;
                else OutputFile << "+" << objective[j][i][k] << "x" << i+1 << "_" << k+1 ;
            }
        }
        
        for(k=0;k<numofJob-1;k++)
        {
            if (objective[j][numofJob-1][k]<0) OutputFile << objective[j][numofJob-1][k] << "x" << numofJob << "_" << k+1 ;
            else OutputFile << "+" << objective[j][numofJob-1][k] << "x" << numofJob << "_" << k+1 ;
        }
        
        if (objective[j][numofJob-1][numofJob-1]<0) OutputFile << objective[j][numofJob-1][numofJob-1] << "x" << numofJob << "_" << numofJob << " - " << "z" << j+1 << " = " << "0" << endl;
        else OutputFile << "+" << objective[j][numofJob-1][numofJob-1] << "x" << numofJob << "_" << numofJob << " - " << "z" << j+1 << " = " << "0" << endl;
    }
    
    OutputFile << "BOUNDS" << endl;
    
    for(j=0;j<numofObjective;j++)
    {
        OutputFile << "z" << j+1 << " free" <<endl;
    }
    
    OutputFile << "BINARIES" << endl;
    
    for(i=0;i<numofJob;i++)
    {
        for(k=0;k<numofJob;k++)
        {
            OutputFile << "x" << i+1 << "_" << k+1 << endl;
        }
    }
    
    OutputFile << "END";
    
    OutputFile.close();
    InputFile.close();
}
