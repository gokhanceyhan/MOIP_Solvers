//
//  constants.h
//  MO_solver
//
//  Created by Gokhan Ceyhan on 2/15/18.
//  Copyright © 2018 Gokhan Ceyhan. All rights reserved.
//

#ifndef constants_h
#define constants_h

static const double OBJ_EPSILON = 1e-4;
static const double BOUND_TOLERANCE = 1e-4; // should be less than the simplex feasibility tolerance (default = 1e-6)
static const int ALL_ND_INDICATOR = -1; // it means that there is no point limit specified by the user.
static const int POINT_LIMIT = 1e+6; // maximm number of nd points to be generated.

// cplex parameters
static const double MIP_RELGAP = 1e-10; // mip relative gap
static const double SIMPLEX_OPTGAP = 1e-9; // simplex optimality gap
static const double BARRIER_CONV = 1e-9; // barrier complementarity tolerance
#endif /* constants_h */
