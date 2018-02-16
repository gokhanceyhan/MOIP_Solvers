//
//  tda.cpp
//  MO_solver
//
//  Created by Gokhan Ceyhan on 1/13/18.
//  Copyright © 2018 Gokhan Ceyhan. All rights reserved.
//

#include "tda.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <string>
#include <limits>
#include <cstdio>
#include <ctime>
#include "tree.hh"
#include "subspace.h"
#include "ilcplex/cplex.h"
#include "math.h"
#include "constants.h"

using namespace std;

//======= auxilary functions===============================================
void tda_form_subtree(tree<subspace>::sibling_iterator, int);
void tda_update_subtree(tree<subspace>::sibling_iterator, int);
void tda_set_lastBound(tree<subspace>::sibling_iterator);
void tda_solve_subspace (tree<subspace>::leaf_iterator); // finds the nondominated point with max z_m value in the subspace
int tda_search_subspace (tree<subspace>::leaf_iterator, vector<float>& p); // search the minimum volume subspaces that contain the given subspace
vector<double> tda_calculateCriteriaScalingCoeffs();
void tda_generatePayoffTable();
double tda_convertBoundTolerance(double ); // given a coverage gap threshold in [0,1]^m space, find an approximate threshold value in the original space
vector<float> tda_solveSingleObjectiveProblem(int);
vector<float> tda_solveForTheInitialSolution();
double tda_calculateCoverageGap(int); // calculates the coverage gap of the n th point when only first n-1 points are representative points
double find_search_space_volume(tree<subspace>::leaf_iterator); // calculates the volume of the subspace
int search_next_subspace_with_rule(); // picks and solves the next subspace. returns 0 if there is no feasible subspace, return 1 if a new nd point is found.
int search_next_subspace(); // tries subspaces in the order of leaf nodes


//======= general global variables==========================================
int tda_num_obj; // number of objectives
int tda_num_point; // number of nondominated points generated
int tda_num_model_solved=0; // number of times cplex is called
string tda_solver_status="Success"; // the status of the multi-objective solver
vector<float> tda_idealPoint;
vector<float> tda_nadirPoint; // approximate nadir unless complete nondominated point set is known.

tree<subspace> tda_spaceT ; // tree of generated subspaces
tree<subspace>::iterator tda_top,tda_root;
vector<vector<float>> tda_points; // set of nondominated points read from the input file
vector<vector<float>> tda_gen_points; // generated nondominated points
vector<float> tda_new_point; // last nondominated point generated
float tda_float_minus_infinity=-numeric_limits<float>::infinity( ); // set to minus infinity for type "float".
vector<float> tda_global_lower_bound; // global lower bound vector at the beginning of the algorithm
double tda_epsilon; // the amount by which the bounds are shifted.
double delta; // amount of shift on the criterion values of a nd point
string tda_path; //  path to the current folder

// ==== cplex and math model variables ======================================
CPXENVptr tda_env=NULL;
CPXLPptr tda_prob=NULL;
int tda_cp_status;
int tda_cp_opt; // stores the status of the model solved (for optimality 101 or 102 must be returned)
int tda_num_mathmodel_col;  // number of columns in the math model
int tda_num_mathmodel_row;  // number of rows in the math model
double *tda_mathmodel_sol;  // points to the array storing the criterion values of the optiaml solution
double *tda_rhs; // points to the array storing the initial right hand side values for bound constraints
char *tda_mathmodel_sense; // points to the array storing the type of ineqaulities for bound constraints

//======= input/output variables ============================================
ofstream tda_file_ND_points;
clock_t tda_start;
double tda_cpu_time;

int tda(int m, string path, double timeLimit, double scaledDelta)
{
    tda_start=clock();
    tda_num_obj=m;
    tda_epsilon = BOUND_TOLERANCE;
    tda_path = path;
    
    int flag,stop_flag; // zero if there is no new nondominated point
    int flag_insert; // used to indicate whether the new point is inserted at the first level
    subspace s;
    vector<float> b,initialP;
    
    tree<subspace>::iterator Iter,Pnode,Cnode; //position of the root of the tree
    tree<subspace>::sibling_iterator sibIter,sibIter_2;
    tree<subspace>::leaf_iterator leafIter;
    
    tda_file_ND_points.open((tda_path+"Output_nd_points.txt").c_str());
    
    //==========================================================================
    //          CPLEX INITIALIZATION and GENERATE FIRST NONDOMINATED POINT
    //==========================================================================
    
    tda_env = CPXopenCPLEX(&tda_cp_status);
    tda_cp_status = CPXsetdefaults (tda_env);
    tda_cp_status = CPXsetintparam(tda_env, CPX_PARAM_SCRIND, 0);
    tda_cp_status = CPXsetintparam(tda_env, CPXPARAM_MIP_Display, 2);
    tda_cp_status = CPXsetintparam(tda_env, CPX_PARAM_CLOCKTYPE, 2); // wall clock time
    tda_cp_status = CPXsetdblparam(tda_env, CPX_PARAM_TILIM, timeLimit); // sets the given time limit
    tda_cp_status = CPXsetdblparam (tda_env, CPX_PARAM_EPGAP, MIP_RELGAP);
    tda_cp_status = CPXsetdblparam (tda_env, CPX_PARAM_EPOPT, SIMPLEX_OPTGAP);
    tda_cp_status = CPXsetdblparam (tda_env, CPX_PARAM_BAREPCOMP, BARRIER_CONV);
    
    tda_prob = CPXcreateprob (tda_env, &tda_cp_status, "mathmodel"); // create problem in CPLEX
    tda_cp_status = CPXreadcopyprob (tda_env, tda_prob, (tda_path+"model.lp").c_str(), NULL); // read "model.lp" into CPLEX
    tda_num_mathmodel_col = CPXgetnumcols (tda_env, tda_prob); // get the number of variables in the math model
    tda_num_mathmodel_row = CPXgetnumrows (tda_env, tda_prob); // get the number of constraints in the math model
    
    // calculate criteria scaling coefficients
    vector<double> criteriaScalCoeffs = tda_calculateCriteriaScalingCoeffs();
    
    // update criterion constraint coefficients for scaling
    int j=0;
    for(int r=tda_num_mathmodel_row-tda_num_obj+1;r<tda_num_mathmodel_row;r++){
        tda_cp_status=CPXchgcoef (tda_env, tda_prob, r, j, -criteriaScalCoeffs[j]);
        j++;
    }
    
    // generate payoff table
    tda_generatePayoffTable();
    
    // set the delta value in the original space
    delta = tda_convertBoundTolerance(scaledDelta);
    
    if (tda_solver_status != "FailedToSolveSingleObjProblem")
    {
        // add bound constraints
        tda_mathmodel_sense = new char [tda_num_obj-1];
        tda_rhs=new double [tda_num_obj-1];
        for (int j=0; j<tda_num_obj-1; j++) {tda_mathmodel_sense[j]='G';    tda_rhs[j]=-CPX_INFBOUND;}
        tda_cp_status = CPXnewrows (tda_env, tda_prob, tda_num_obj-1, tda_rhs , tda_mathmodel_sense, NULL, NULL);
        for (int j=0; j<tda_num_obj-1; j++){
            tda_cp_status=CPXchgcoef (tda_env, tda_prob, tda_num_mathmodel_row + j, j, 1);
            tda_cp_status=CPXchgcoef (tda_env, tda_prob, tda_num_mathmodel_row + j, tda_num_obj-1, -1);
        }
        
        // tda_cp_status = CPXwriteprob (tda_env, tda_prob, (tda_path+"coverage_prob.lp").c_str(), NULL);
        initialP = tda_solveForTheInitialSolution();
        
        //==========================================================================
        //                    initialize the tree
        //==========================================================================
        
        for (int j=1; j<=tda_num_obj; j++) {
            tda_global_lower_bound.push_back(tda_float_minus_infinity);
        }
        
        s.set_bound(tda_global_lower_bound);
        s.set_point(tda_global_lower_bound);
        s.set_sub_status(-1);
        s.set_mod_status(-1);
        
        tda_top=tda_spaceT.begin();
        tda_root=tda_spaceT.insert(tda_top, s); // insert the global lower bound vector to the root
        Pnode=tda_root;
        
        for (int j=0; j<tda_num_obj-2; j++) {
            Cnode=tda_spaceT.append_child(Pnode,s);
            Pnode=Cnode;
        }
        
        tda_gen_points.push_back(initialP);
        // update nadir point
        for(int j=0; j<tda_num_obj-1;j++) if(initialP[j]<tda_nadirPoint[j]) tda_nadirPoint[j]=initialP[j];
        tda_num_point=1;
        (*Cnode).set_point(initialP);
        (*Cnode).set_mod_status(1);
        (*Cnode).set_sub_status(1);
        
        tda_new_point=initialP;
        
        flag=1;
    }
    else
    {
        cout << "Failed to solve single objective problem.\n";
        flag = 0;
    }
    
    // check the stopping conditions
    if(( clock() - tda_start ) / (double) CLOCKS_PER_SEC >= timeLimit)
    {
        flag = 0;
        tda_solver_status = "AbortTimeLimit";
    }
    
    //==========================================================================
    // Generation loop
    //==========================================================================
    
    while (flag) {
        
        b=tda_global_lower_bound;
        b[0]=tda_new_point[0]+tda_epsilon;
        
        s.set_bound(b);
        
        //==================================================================================================================================
        // traverse the first level and append the child to the appropriate place. In this way, siblings will be always sorted at each level.
        //==================================================================================================================================
        sibIter=tda_spaceT.begin(tda_root);
        
        flag_insert=0;
        
        while (sibIter!=tda_spaceT.end(tda_root)) { // first level traversing====================
            
            
            if (tda_spaceT.number_of_siblings(sibIter)<=0) {
                
                tda_spaceT.insert_after(sibIter, s);
                sibIter++;
                tda_form_subtree(sibIter,1);
                flag_insert=1;
                break;
            }
            
            else if ((*sibIter).get_bound()[0]>=b[0]) {
                
                flag_insert=1;
                
                if ((*sibIter).get_bound()[0]!=b[0]) {
                    
                    tda_spaceT.insert(sibIter, s);
                    sibIter--;
                    tda_form_subtree(sibIter,1);
                    
                }
                
                break;
            }
            
            sibIter++;
            
        }
        
        if (flag_insert==0) {
            sibIter--;
            sibIter=tda_spaceT.insert_after(sibIter, s);
            tda_form_subtree(sibIter,1);
        }
        
        //==================================================================================================================================
        //update the subtrees of all siblings on the left of the new node
        //==================================================================================================================================
        sibIter_2=tda_spaceT.begin(tda_root);
        
        while (sibIter_2!=sibIter){
            tda_update_subtree(sibIter_2,1);
            sibIter_2++;
        }
        
        //==================================================================================================================================
        // Determine the subspace search order. Pick a subspace and find a nondominated point there if exists.
        // Otherwise, skip to the next subspace in the search order.
        //==================================================================================================================================
        
        tda_new_point=tda_global_lower_bound;
        stop_flag = search_next_subspace();
        
        if (stop_flag==0) {
            flag = 0; // no more nondominated point
        } else
        {
            // add the new nondominated point to the list
            tda_gen_points.push_back(tda_new_point);
            // update nadir point
            for(int j=0; j<tda_num_obj-1;j++) if(tda_new_point[j]<tda_nadirPoint[j]) tda_nadirPoint[j]=tda_new_point[j];
            tda_num_point++;
        }
        
        // check the stopping conditions
        double elapsedTime = ( clock() - tda_start ) / (double) CLOCKS_PER_SEC;
        if (elapsedTime >= timeLimit)
        {
            flag = 0; // time limit exceeded
            tda_solver_status = "AbortTimeLimit";
        }
        
        cout << "point: " << tda_num_point << " model solved: " << tda_num_model_solved << " elapsed time " << elapsedTime <<  endl;
        
        //==================================================================================================================================
        //display the tree
        //==================================================================================================================================
        
        /*cout << "tree size: " << tda_spaceT.size() << "\n";
         
         Iter=tda_spaceT.begin(tda_root);
         while (Iter!=tda_spaceT.end(tda_root)) {
         
         if (tda_spaceT.depth(Iter)==tda_num_obj-2)
         {
         cout << "bound: ";
         for (int i=0; i<tda_num_obj-1; i++) {
         cout << (*Iter).get_bound()[i] << " ";
         }
         
         cout << "sub_status: " << (*Iter).get_sub_status() << " ";
         cout << "mod_status: " << (*Iter).get_mod_status() << " ";
         cout << "point: ";
         for (int i=0; i<tda_num_obj; i++) {
         cout << (*Iter).get_point()[i] << " ";
         }
         }
         
         cout << endl;
         
         Iter++;
         }
         
         cout << "new point: ";
         for (int i=0; i<tda_num_obj; i++) {
         
         cout << tda_new_point[i] << " ";
         }
         cout << "\n";
         */
        //========================================================================
        
    }
    
    tda_cpu_time = ( clock() - tda_start ) / (double) CLOCKS_PER_SEC;
    
    // display the generated points
    //cout << "Nondominated point list: \n";
    tda_file_ND_points << "#Solver type:" << endl << "rMOCO-S_tda" << endl;
    tda_file_ND_points << "#Solver status:" << endl << tda_solver_status << endl;
    if (tda_solver_status != "FailedToSolveSingleObjProblem") {
        tda_file_ND_points << "#Number of nondominated points:" << endl << tda_num_point << endl;
        tda_file_ND_points << "#Number of models solved:" << endl << tda_num_model_solved << endl;
        tda_file_ND_points << "#Elapsed time (seconds):" << endl << tda_cpu_time << endl;
        tda_file_ND_points << "#Ideal point:" << endl;
        for (int j=0; j<tda_num_obj-1; j++) tda_file_ND_points << tda_idealPoint[j]*criteriaScalCoeffs[j] << " ";
        tda_file_ND_points << endl;
        tda_file_ND_points << "#Incumbent nadir point:" << endl;
        for (int j=0; j<tda_num_obj-1; j++) tda_file_ND_points << tda_nadirPoint[j]*criteriaScalCoeffs[j] << " ";
        tda_file_ND_points << endl;
        
        // calculate maximum efficient range
        double scalingCoeff = 0.0;
        for(int j=0;j<tda_num_obj-1;j++){
            if(tda_idealPoint[j]-tda_nadirPoint[j]>scalingCoeff)
                scalingCoeff = tda_idealPoint[j]-tda_nadirPoint[j];
        }
        
        tda_file_ND_points << "#The set of nondominated points:" << endl;
        for (int i=0; i < tda_gen_points.size(); i++) {
            for (int j=0; j<tda_num_obj-1; j++) {
                tda_file_ND_points << tda_gen_points[i][j]*criteriaScalCoeffs[j]-delta << " ";
            }
            tda_file_ND_points << endl;
        }
    } else {
        tda_file_ND_points << "#Cplex error code:" << endl << tda_cp_status << endl;
        tda_file_ND_points << "#Cplex optimization status:" << endl << tda_cp_opt << endl;
    }
    
    tda_file_ND_points.close();
    CPXfreeprob(tda_env, &tda_prob);
    return 1;
    
}

void tda_solve_subspace (tree<subspace>::leaf_iterator sub)
{
    vector<float> ndp;
    vector<float> bound;
    int status=-1;
    
    
    for (int i=1; i<=tda_num_obj; i++) {
        ndp.push_back(tda_float_minus_infinity);
    }
    
    (*sub).set_mod_status(0);
    
    status=tda_search_subspace(sub,ndp);
    
    if (status==-1)
    {
        tda_num_model_solved++;
        (*sub).set_mod_status(1);
        bound=(*sub).get_bound();
        //==========================================================================
        //                           MODEL SOLUTION
        //==========================================================================
        
        for (int j=0; j<tda_num_obj-1; j++)
            tda_cp_status=CPXchgcoef (tda_env, tda_prob, tda_num_mathmodel_row + j, -1, bound[j]);
        
        // tda_cp_status = CPXwriteprob (tda_env, tda_prob, (tda_path+"myprob.lp").c_str(), NULL);
        
        tda_cp_status = CPXmipopt (tda_env, tda_prob); // solve the problem
        tda_cp_opt=CPXgetstat(tda_env,tda_prob);
        //cout << tda_cp_opt << endl;
        if (tda_cp_opt==101 || tda_cp_opt==102)
        {
            tda_cp_status = CPXgetx (tda_env, tda_prob, tda_mathmodel_sol, 0, tda_num_obj-1); // copies the optimal obj values to the "tda_mathmodel_sol"
            status=1;
            for (int j=0; j<tda_num_obj; j++) ndp[j]=(float)tda_mathmodel_sol[j]+delta;
            (*sub).set_mod_status(1);
            
        }
        else
            status=0;
        
        //===========================================================================
    }
    
    // ======= display ============
    /*cout << "RESULT FOR THIS SUBSPACE:" << endl;
     cout << "status:" << status << endl;
     cout << "point: " ;
     for (int i=0; i<tda_num_obj; i++) {
     cout << ndp[i] << " ";
     }
     cout << endl;*/
    
    (*sub).set_sub_status(status);
    (*sub).set_point(ndp);
    
    
}

void tda_update_subtree(tree<subspace>::sibling_iterator p_it, int d)
{
    vector<float> b;
    int flag,last_sib;
    tree<subspace>::sibling_iterator s_it,s_newP,next_sib;
    subspace s;
    
    s.set_bound(tda_global_lower_bound);
    s.set_point(tda_global_lower_bound);
    s.set_sub_status(-1);
    s.set_mod_status(-1);
    
    if (d<tda_num_obj-2) {
        b=(*p_it).get_bound();
        
        flag=0;
        for (int i=0; i<d; i++) {
            if (b[i] > tda_new_point[i]) { // which means the new_point is not in this subspace, no need to update
                flag=1;
                break;
            }
        }
        b[d]=tda_new_point[d]+tda_epsilon;
        s.set_bound(b);
        
        if (flag==0) { // update the children
            
            s_it=tda_spaceT.begin(p_it); // points to the first child
            
            if (tda_spaceT.number_of_siblings(s_it)<=0) {
                tda_spaceT.insert_after(s_it, s);
                s_it++;
                tda_form_subtree(s_it, d+1);
                s_it--;
                tda_update_subtree(s_it, d+1);
            }
            else {
                // Step 1. Insert the node to the appropriate place.
                last_sib=1;
                int insertion_index = 0;
                while (tda_spaceT.index(s_it)<tda_spaceT.number_of_children(p_it)-1) {
                    
                    next_sib=tda_spaceT.next_sibling(s_it);
                    
                    // check if we need to insert a new node after this one
                    if ( (*s_it).get_bound()[d] <= b[d])
                    {
                        if ((*next_sib).get_bound()[d] >= b[d])
                        {
                            last_sib=0; // the new node will not be the last sibling
                            insertion_index = tda_spaceT.index(s_it);
                            
                            if (b[d]!=(*s_it).get_bound()[d] && b[d]!=(*next_sib).get_bound()[d]) {
                                s_newP=tda_spaceT.insert_after(s_it, s);
                                tda_form_subtree(s_newP, d+1);
                            }
                            break;
                        }
                    }
                    s_it++;
                }
                if (last_sib==1) { // if the new node is to be inserted to the end
                    insertion_index = tda_spaceT.index(s_it);
                    s_newP=tda_spaceT.insert_after(s_it, s);
                    tda_form_subtree(s_newP, d+1);
                }
                
                // Step 2. Update all subtrees at the left-hand side of the inserted node
                s_it=tda_spaceT.begin(p_it);
                while (tda_spaceT.index(s_it) <= insertion_index) {
                    tda_update_subtree(s_it, d+1);
                    s_it++;
                }
            }
        }
    }
    else // this means that it is a leaf node
    {
        b=(*p_it).get_bound();
        flag=0;
        for (int i=0; i<tda_num_obj; i++) {
            if (b[i] > tda_new_point[i]) { // which means the new_point is not in this subspace, no need to solve this leaf node
                flag=1;
                break;
            }
        }
        
        if (flag==0) // update the leaf
        {
            tda_set_lastBound(p_it);
        }
        
    }
}

void tda_form_subtree(tree<subspace>::sibling_iterator p_it, int d)
{
    vector<float> b,p;
    tree<subspace>::sibling_iterator sibStart_it,sibEnd_it;
    tree<subspace>::iterator s_it;
    subspace s;
    
    s.set_bound(tda_global_lower_bound);
    s.set_point(tda_global_lower_bound);
    s.set_sub_status(-1);
    s.set_mod_status(-1);
    
    if (d < tda_num_obj-2) { // NOT a leaf node
        
        b=tda_global_lower_bound;
        
        // copy the first d bound of the parent node
        for (int i=0; i < d; i++) {
            b[i]=(*p_it).get_bound()[i];
        }
        
        if (tda_spaceT.number_of_siblings(p_it) > tda_spaceT.index(p_it)) // if parent is not the last sibling
        {
            // copy the left subtree
            sibStart_it = tda_spaceT.begin(tda_spaceT.previous_sibling(p_it));
            sibEnd_it = tda_spaceT.end(tda_spaceT.previous_sibling(p_it));
            tda_spaceT.append_children(p_it, sibStart_it, sibEnd_it);
            
            s_it = tda_spaceT.begin(p_it);
            // update the bounds and leaf node info
            while (s_it!=tda_spaceT.end(p_it)){
                for (int i=d; i < tda_num_obj-2; i++) b[i]=(*s_it).get_bound()[i];
                (*s_it).set_bound(b);
                if(tda_spaceT.depth(s_it)==(tda_num_obj-2)){
                    tda_set_lastBound(s_it);
                }
                s_it++;
            }
        }
        else // if parent is the last sibling
        {
            for (int i=d; i < tda_num_obj-2; i++)
                b[i]=tda_float_minus_infinity;
            s.set_bound(b);
            s_it=tda_spaceT.append_child(p_it,s);
            tda_form_subtree(s_it, d+1);
        }
    }
    else { // this means that it is a leaf node
        
        // we need to update the bound in the z_(m-1) objective
        tda_set_lastBound(p_it);
    }
    
}


void tda_set_lastBound(tree<subspace>::sibling_iterator it)
{
    int flag;
    vector<float> b;
    float max_z_m=tda_global_lower_bound[tda_num_obj-2];
    
    
    b=(*it).get_bound();
    
    
    for (int j=0; j<tda_gen_points.size(); j++) {
        
        flag=1;
        for (int k=0; k<tda_num_obj-2; k++) {
            if (b[k] > tda_gen_points[j][k]) {
                flag=0;
                break;
            }
        }
        
        if ( (flag==1) && (tda_gen_points[j][tda_num_obj-2] > max_z_m)) {
            
            max_z_m=tda_gen_points[j][tda_num_obj-2];
        }
        
    }
    
    b[tda_num_obj-2]=max_z_m+tda_epsilon;
    (*it).set_bound(b);
    (*it).set_sub_status(-1);
    (*it).set_mod_status(-1);
    
}

int tda_search_subspace (tree<subspace>::leaf_iterator sub, vector<float>& p)
{
    /* return 0 if the subspace is infeasible and we do not need to solve the model
     return -1 if we need to solve the model
     return 1 if optimal solution is guaranteed in this subspace and we do not need to solve the model
     */
    
    int flag, status;
    tree<subspace>::sibling_iterator sibIter,prev_sib;
    tree<subspace>::leaf_iterator leafIter;
    
    sibIter=sub;
    status=-1;
    
    // display the current subspace===================================
    /*cout << endl << endl <<"========CURRENT NODE=======" << endl;
     cout << "bound: ";
     for (int i=0; i<tda_num_obj-1; i++) {
     cout << (*sibIter).get_bound()[i] << " ";
     }
     cout << endl;*/
    //================================================================
    
    //cout << "========SCANNED NODES=======" << endl;
    
    if (tda_spaceT.depth(sibIter) == 1)
    {
        if (tda_spaceT.index(sibIter)!= 0)
        {
            prev_sib=tda_spaceT.previous_sibling(sibIter); // step back to the previous sibling
            
            // display the subspaces================================================
            /*cout << "bound: ";
             for (int i=0; i<tda_num_obj-1; i++) {
             cout << (*prev_sib).get_bound()[i] << " ";
             }
             
             cout << "sub_status: " << (*prev_sib).get_sub_status() << " ";
             
             if ((*prev_sib).get_sub_status()==1)
             {
             for (int i=0; i<tda_num_obj; i++) cout << (*prev_sib).get_point()[i] << " ";
             }
             
             cout << endl;*/
            
            //=======================================================================
            
            // compare with the nodes
            
            // bounds of the current subspace must be greater than or equal to each bound of the compared subspace
            flag=1;
            for (int i=0; i<tda_num_obj-1; i++)
            {
                if((*sub).get_bound()[i] < (*prev_sib).get_bound()[i] )
                {
                    flag=0; // mathematical model needs to be solved.
                    break;
                }
            }
            
            if (flag==1)
            {
                
                if ((*prev_sib).get_sub_status()==0)
                {
                    status=0; // current subspace is also infeasible.
                }
                else
                {
                    for (int i=0; i<tda_num_obj-1; i++)
                    {
                        if((*sub).get_bound()[i] > (*prev_sib).get_point()[i] )
                        {
                            flag=0; // mathematical model needs to be solved.
                            break;
                        }
                    }
                    
                    if (flag==1) { status=1; p=(*prev_sib).get_point();}
                }
            }
        }
    }
    else
    {
        
        while(tda_spaceT.depth(sibIter) != 1) // go up the tree until reaching to the first level of the tree
        {
            sibIter=tda_spaceT.parent(sibIter);
        }
        
        
        if (tda_spaceT.index(sibIter)!= 0)
        {
            
            prev_sib=tda_spaceT.previous_sibling(sibIter); // step back to the previous sibling
            
            leafIter=tda_spaceT.begin_leaf(prev_sib);
            
            // traversing previous largest subtree
            
            while (leafIter!=tda_spaceT.end_leaf(prev_sib) && status==-1)
            {
                
                // display the subspaces================================================
                /*cout << "bound: ";
                 for (int i=0; i<tda_num_obj-1; i++) {
                 cout << (*leafIter).get_bound()[i] << " ";
                 }
                 
                 cout << "sub_status: " << (*leafIter).get_sub_status() << " ";
                 cout << "mod_status: " << (*leafIter).get_mod_status() << " ";
                 
                 if ((*leafIter).get_sub_status()==1)
                 {
                 for (int i=0; i<tda_num_obj; i++) cout << (*leafIter).get_point()[i] << " ";
                 }
                 
                 cout << endl;*/
                
                //=======================================================================
                
                // compare the nodes
                
                // bounds of the current subspace must be greater than or equal to each bound of the compared subspace
                flag=1;
                for (int i=0; i<tda_num_obj-1; i++)
                {
                    if((*sub).get_bound()[i] < (*leafIter).get_bound()[i] )
                    {
                        flag=0; // mathematical model needs to be solved.
                        break;
                    }
                }
                
                if (flag==1)
                {
                    
                    if ((*leafIter).get_sub_status()==0)
                    {
                        status=0; // current subspace is also infeasible.
                    }
                    else
                    {
                        for (int i=0; i<tda_num_obj-1; i++)
                        {
                            if((*sub).get_bound()[i] > (*leafIter).get_point()[i] )
                            {
                                flag=0; // mathematical model needs to be solved.
                                break;
                            }
                        }
                        
                        if (flag==1) { status=1; p=(*leafIter).get_point();}
                    }
                }
                
                leafIter++;
                
            }
        }
        
        // traversing current largest subtree
        leafIter=tda_spaceT.begin_leaf(sibIter);
        
        while (leafIter!=sub && status==-1)
        {
            
            // display the subspaces======================================
            /*cout << "bound: ";
             for (int i=0; i<tda_num_obj-1; i++) {
             cout << (*leafIter).get_bound()[i] << " ";
             }
             
             cout << "sub_status: " << (*leafIter).get_sub_status() << " ";
             cout << "mod_status: " << (*leafIter).get_mod_status() << " ";
             
             if ((*leafIter).get_sub_status()==1)
             {
             for (int i=0; i<tda_num_obj; i++) cout << (*leafIter).get_point()[i] << " ";
             }
             
             cout << endl;*/
            
            //=============================================================
            
            // compare the nodes
            
            // bounds of the current subspace must be greater than or equal to each bound of the compared subspace
            flag=1;
            for (int i=0; i<tda_num_obj-1; i++)
            {
                if((*sub).get_bound()[i] < (*leafIter).get_bound()[i] )
                {
                    flag=0; // mathematical model needs to be solved.
                    break;
                }
            }
            
            if (flag==1)
            {
                
                if ((*leafIter).get_sub_status()==0)
                {
                    status=0; // current subspace is also infeasible.
                }
                else
                {
                    for (int i=0; i<tda_num_obj-1; i++)
                    {
                        if((*sub).get_bound()[i] > (*leafIter).get_point()[i] )
                        {
                            flag=0; // mathematical model needs to be solved.
                            break;
                        }
                    }
                    
                    if (flag==1)
                    {
                        status=1; p=(*leafIter).get_point();
                    }
                }
            }
            
            leafIter++;
            
        }
    }
    return (status);
    
}

vector<double> tda_calculateCriteriaScalingCoeffs(){
    
    vector<double> coeffs(tda_num_obj);
    double coeff;
    
    int j = 0;
    for(int r=tda_num_mathmodel_row-tda_num_obj+1;r<tda_num_mathmodel_row;r++){
        coeffs[j]=0.0;
        for(int c=tda_num_obj;c<tda_num_mathmodel_col;c++){
            CPXgetcoef (tda_env, tda_prob, r, c, &coeff);
            coeffs[j]+=pow(coeff,2);
        }
        coeffs[j]=pow(coeffs[j],0.5);
        j++;
    }
    coeffs[tda_num_obj-1]=1.0;
    
    for (int j=0; j<tda_num_obj; j++) {
        coeffs[j]=1.0;
    }
    
    return coeffs;
}

double tda_convertBoundTolerance(double scaledDelta){
    // call this method after tda_generatePayoffTable() method so that an incumbent nadir point is available.
    // assumes that the scaledDelta argument of the program is in the [0,1]^m space.
    double maxEfficientRange = 0.0;
    for(int i=0;i<tda_num_obj-1;i++){
        if(tda_idealPoint[i]-tda_nadirPoint[i] > maxEfficientRange)
            maxEfficientRange = tda_idealPoint[i]-tda_nadirPoint[i];
    }
    
    return maxEfficientRange * scaledDelta;
}

void tda_generatePayoffTable(){
    for(int j=0; j<tda_num_obj-1;j++){
        tda_idealPoint.push_back(-numeric_limits<float>::infinity( ));
        tda_nadirPoint.push_back(numeric_limits<float>::infinity( ));
    }
    vector<float> point;
    for(int j=0; j<tda_num_obj-1;j++){
        point = tda_solveSingleObjectiveProblem(j);
        if(tda_solver_status == "FailedToSolveSingleObjProblem") break;
        for(int k=0; k<tda_num_obj-1;k++){
            if(point[k]>tda_idealPoint[k]) tda_idealPoint[k]=point[k];
            if(point[k]<tda_nadirPoint[k]) tda_nadirPoint[k]=point[k];
        }
        tda_cp_status=CPXchgcoef (tda_env, tda_prob, -1, j, OBJ_EPSILON);
    }
    tda_cp_status=CPXchgcoef (tda_env, tda_prob, -1, tda_num_obj-1, 1.0);
}

vector<float> tda_solveSingleObjectiveProblem(int objIndex){
    vector<float> point(tda_num_obj-1);
    tda_cp_status=CPXchgcoef (tda_env, tda_prob, -1, tda_num_obj-1, 0.0);
    tda_cp_status=CPXchgcoef (tda_env, tda_prob, -1, objIndex, 1.0);
    tda_mathmodel_sol = new double [tda_num_obj];
    tda_cp_status = CPXmipopt (tda_env, tda_prob); // solve the problem
    tda_cp_opt=CPXgetstat(tda_env,tda_prob);
    
    if (tda_cp_opt==101 || tda_cp_opt==102)
    {
        tda_cp_status = CPXgetx (tda_env, tda_prob, tda_mathmodel_sol, 0, tda_num_obj-1); // copies the optimal obj values to the "tda_mathmodel_sol"
        for (int j=0; j<tda_num_obj; j++) point[j]= (float)tda_mathmodel_sol[j];
    } else {
        tda_solver_status = "FailedToSolveSingleObjProblem";
    }
    return point;
}

vector<float> tda_solveForTheInitialSolution(){
    // call this function after tda_generatePayoffTable() method and constructing the coverega constraints
    vector<float> point(tda_num_obj-1);
    
    // modify the model to obtain achievement scalarizing program
    tda_cp_status=CPXchgcoef (tda_env, tda_prob, -1, tda_num_obj-1, -1);
    for (int j=0; j<tda_num_obj-1; j++){
        tda_cp_status=CPXchgcoef (tda_env, tda_prob, tda_num_mathmodel_row + j, tda_num_obj-1, 1);
        tda_cp_status=CPXchgcoef (tda_env, tda_prob, tda_num_mathmodel_row + j, -1, tda_idealPoint[j]);
    }
    // tda_cp_status = CPXwriteprob (tda_env, tda_prob, (tda_path+"asp.lp").c_str(), NULL);
    
    tda_mathmodel_sol = new double [tda_num_obj];
    tda_cp_status = CPXmipopt (tda_env, tda_prob); // solve the problem
    tda_cp_opt=CPXgetstat(tda_env,tda_prob);
    
    if (tda_cp_opt==101 || tda_cp_opt==102)
    {
        tda_cp_status = CPXgetx (tda_env, tda_prob, tda_mathmodel_sol, 0, tda_num_obj-1); // copies the optimal obj values to the "tda_mathmodel_sol"
        for (int j=0; j<tda_num_obj; j++) point[j]= (float)tda_mathmodel_sol[j];
    } else {
        tda_solver_status = "FailedToSolveSingleObjProblem";
    }
    
    // reset the model
    tda_cp_status=CPXchgcoef (tda_env, tda_prob, -1, tda_num_obj-1, 1);
    for (int j=0; j<tda_num_obj-1; j++){
        tda_cp_status=CPXchgcoef (tda_env, tda_prob, tda_num_mathmodel_row + j, tda_num_obj-1, -1);
    }
    
    return point;
}

double tda_calculateCoverageGap(int n){
    double alpha_min = 1.0;
    double alpha = 0.0;
    double scalingCoeff = 0.0;
    for(int j=0;j<tda_num_obj-1;j++){
        if(tda_idealPoint[j]-tda_nadirPoint[j]>scalingCoeff)
            scalingCoeff = tda_idealPoint[j]-tda_nadirPoint[j];
    }
    
    for(int i=0;i<n;i++){
        alpha = 0.0;
        for(int j=0;j<tda_num_obj-1;j++){
            if((tda_gen_points[n][j]-tda_gen_points[i][j])/scalingCoeff>alpha)
                alpha = (tda_gen_points[n][j]-tda_gen_points[i][j])/scalingCoeff;
        }
        if (alpha < alpha_min) alpha_min = alpha;
    }
    return alpha_min;
}

int search_next_subspace_with_rule(){
    
    tree<subspace>::leaf_iterator leafIter, p;
    double max_vol,vol;
    
    while (1) {
        leafIter=tda_spaceT.begin_leaf();
        max_vol = -1.0;
        while (leafIter!=tda_spaceT.end_leaf()) {
            if ((*leafIter).get_sub_status()==-1) {
                vol = find_search_space_volume(leafIter);
                if(vol > max_vol){
                    max_vol = vol;
                    p=leafIter;
                }
            }
            leafIter++;
        }
        if (max_vol>-1.0) {
            tda_solve_subspace(p);
        } else {
            return 0;
        }
        
        if ((*p).get_sub_status()==1) {
            tda_new_point=(*p).get_point();
            return 1;
        }
    }
    return 1;
}

int search_next_subspace(){
    tree<subspace>::leaf_iterator leafIter;
    
    leafIter=tda_spaceT.begin_leaf();
    
    while (leafIter!=tda_spaceT.end_leaf()) {
        if ((*leafIter).get_sub_status()==-1) {
            tda_solve_subspace(leafIter);
        }
        if ((*leafIter).get_sub_status()==1) {
            tda_new_point=(*leafIter).get_point();
            return 1;
        }
        leafIter++;
    }
    
    return 0;
}

double find_search_space_volume(tree<subspace>::leaf_iterator sub){
    double v = tda_idealPoint[0]-(*sub).get_bound()[0];
    for (int i=1; i<tda_num_obj; i++) {
        if((tda_idealPoint[i]-(*sub).get_bound()[i])<v)
            v = (tda_idealPoint[i]-(*sub).get_bound()[i]);
    }
    return v;
}
