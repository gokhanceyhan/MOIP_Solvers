//
//  constants.h
//  MO_solver
//
//  Created by Gokhan Ceyhan on 2/15/18.
//  Copyright © 2018 Gokhan Ceyhan. All rights reserved.
//

#ifndef constants_h
#define constants_h

static const double OBJ_EPSILON = 0.0001;
static const int ALL_ND_INDICATOR = -1; // it means that there is no point limit specified by the user.
static const int POINT_LIMIT = 1000000; // maximm number of nd points to be generated.

// cplex parameters
static const double MIP_RELGAP = 1e-10; // mip relative gap
#endif /* constants_h */
