//
//  main.cpp
//  MO_solver
//
//  Created by Gokhan Ceyhan on 2/21/16.
//  Copyright (c) 2016 Gokhan Ceyhan. All rights reserved.
//

#include <iostream>
#include "exact.h"
#include "sba.h"
#include "tda.h"
#include "generateModel.h"
#include "parameters.h"

using namespace std;

int main(int argc, const char * argv[]) {
    
    // read the program settings
    string path = argv[1];
    parameters params = generateModel(path);
    int solverType = params.getSolverType();
    int numofObj = params.getNumOfObj();
    double timeLimit = params.getTimeLimit();
    int pointLimit = params.getPointLimit();
    double bound_tolerance = params.getBoundTolerance();
    double delta = params.getDelta();
    
    switch(solverType){
        case 1:
            // generate all nondominated points
            cout << "EXACT" << endl;
            exact(numofObj, path, timeLimit, bound_tolerance);
            break;
        case 2:
            // generate worst covered representative points
            cout << "SBA" << endl;
            sba(numofObj+1, path, timeLimit, pointLimit, bound_tolerance);
            break;
        case 3:
            // generate representative points achieving the desired coverage gap
            cout << "TDA" << endl;
            tda(numofObj+1, path, timeLimit, bound_tolerance, delta);
            break;
    }
    
    
    return 0;
}



