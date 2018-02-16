//
//  parameters.h
//  MO_solver
//
//  Created by Gokhan Ceyhan on 10/17/17.
//  Copyright © 2017 Gokhan Ceyhan. All rights reserved.
//

#ifndef parameters_h
#define parameters_h

using namespace std;

class parameters {
    
private:
    int numOfObj;
    int solverType; // 1: nMOCO-S, 2: sba, 3: tda
    double timeLimit; // in secs
    int pointLimit;
    double delta; // territory length
    
public:
    inline parameters(int solverType, int numOfObj, double timeLimit, int pointLimit, double delta){
        this->solverType = solverType;
        this->numOfObj = numOfObj;
        this->timeLimit = timeLimit;
        this->pointLimit = pointLimit;
        this->delta = delta;
    }
    
    inline int getNumOfObj(){
        return this->numOfObj;
    }
    
    inline int getSolverType(){
        return this->solverType;
    }
    
    inline double getTimeLimit(){
        return this->timeLimit;
    }
    
    inline int getPointLimit(){
        return this->pointLimit;
    }
    
    inline double getDelta(){
        return this->delta;
    }
};

#endif /* parameters_h */
