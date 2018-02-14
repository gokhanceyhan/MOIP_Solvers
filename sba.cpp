//
//  sba.cpp
//  MO_solver
//
//  Created by Gokhan Ceyhan on 9/26/17.
//  Copyright © 2017 Gokhan Ceyhan. All rights reserved.
//


#include "sba.h"
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

using namespace std;

//======= auxilary functions===============================================
void sba_form_subtree(tree<subspace>::sibling_iterator, int);
void sba_update_subtree(tree<subspace>::sibling_iterator, int);
void sba_set_lastBound(tree<subspace>::sibling_iterator);
void sba_solve_subspace (tree<subspace>::leaf_iterator); // finds the nondominated point with max z_m value in the subspace
int sba_search_subspace (tree<subspace>::leaf_iterator, vector<float>& p); // search the minimum volume subspaces that contain the given subspace
vector<double> calculateCriteriaScalingCoeffs();
void generatePayoffTable();
vector<float> solveSingleObjectiveProblem(int);
vector<float> solveForTheInitialSolution();
double calculateCoverageGap(int); // calculates the coverage gap of the n th point when only first n-1 points are representative points

//======= general global variables==========================================
int sba_num_obj; // number of objectives
int sba_num_point; // number of nondominated points generated
int sba_num_model_solved=0; // number of times cplex is called
string sba_solver_status="Success"; // the status of the multi-objective solver
vector<float> sba_idealPoint;
vector<float> sba_nadirPoint; // approximate nadir unless complete nondominated point set is known.

tree<subspace> sba_spaceT ; // tree of generated subspaces
tree<subspace>::iterator sba_top,sba_root;
vector<vector<float>> sba_points; // set of nondominated points read from the input file
vector<vector<float>> sba_gen_points; // generated nondominated points
vector<float> sba_new_point; // last nondominated point generated
float sba_float_minus_infinity=-numeric_limits<float>::infinity( ); // set to minus infinity for type "float".
vector<float> sba_global_lower_bound; // global lower bound vector at the beginning of the algorithm
double sba_epsilon; // the amount by which the bounds are shifted.
string sba_path; //  path to the current folder

// ==== cplex and math model variables ======================================
CPXENVptr sba_env=NULL;
CPXLPptr sba_prob=NULL;
int sba_cp_status;
int sba_cp_opt; // stores the status of the model solved (for optimality 101 or 102 must be returned)
int sba_num_mathmodel_col;  // number of columns in the math model
int sba_num_mathmodel_row;  // number of rows in the math model
double *sba_mathmodel_sol;  // points to the array storing the criterion values of the optiaml solution
double *sba_rhs; // points to the array storing the initial right hand side values for bound constraints
char *sba_mathmodel_sense; // points to the array storing the type of ineqaulities for bound constraints

//======= input/output variables ============================================
ofstream sba_file_ND_points;
clock_t sba_start;
double sba_cpu_time;

int sba(int m, string path, double timeLimit, int pointLimit, double bound_tolerance)
{
    sba_start=clock();
    sba_num_obj=m;
    sba_epsilon = bound_tolerance;
    sba_path = path;
    
    int flag,stop_flag; // zero if there is no new nondominated point
    int flag_insert; // used to indicate whether the new point is inserted at the first level
    subspace s;
    vector<float> b,initialP;
    
    tree<subspace>::iterator Iter,Pnode,Cnode; //position of the root of the tree
    tree<subspace>::sibling_iterator sibIter,sibIter_2;
    tree<subspace>::leaf_iterator leafIter;
    
    sba_file_ND_points.open((sba_path+"Output_nd_points.txt").c_str());
    
    //==========================================================================
    //          CPLEX INITIALIZATION and GENERATE FIRST NONDOMINATED POINT
    //==========================================================================
    
    sba_env = CPXopenCPLEX(&sba_cp_status);
    sba_cp_status = CPXsetdefaults (sba_env);
    sba_cp_status = CPXsetintparam(sba_env, CPX_PARAM_SCRIND, 0);
    sba_cp_status = CPXsetintparam(sba_env, CPXPARAM_MIP_Display, 2);
    sba_cp_status = CPXsetintparam(sba_env, CPX_PARAM_CLOCKTYPE, 2); // wall clock time
    sba_cp_status = CPXsetdblparam(sba_env, CPX_PARAM_TILIM, timeLimit); // sets the given time limit
    
    sba_prob = CPXcreateprob (sba_env, &sba_cp_status, "mathmodel"); // create problem in CPLEX
    sba_cp_status = CPXreadcopyprob (sba_env, sba_prob, (sba_path+"model.lp").c_str(), NULL); // read "model.lp" into CPLEX
    sba_num_mathmodel_col = CPXgetnumcols (sba_env, sba_prob); // get the number of variables in the math model
    sba_num_mathmodel_row = CPXgetnumrows (sba_env, sba_prob); // get the number of constraints in the math model
    
    // calculate criteria scaling coefficients
    vector<double> criteriaScalCoeffs = calculateCriteriaScalingCoeffs();
    
    // update criterion constraint coefficients for scaling
    int j=0;
    for(int r=sba_num_mathmodel_row-sba_num_obj+1;r<sba_num_mathmodel_row;r++){
        sba_cp_status=CPXchgcoef (sba_env, sba_prob, r, j, -criteriaScalCoeffs[j]);
        j++;
    }
    
    // generate payoff table
    generatePayoffTable();
    
    if (sba_solver_status != "ProblemInfeasible")
    {
        // add bound constraints
        sba_mathmodel_sense = new char [sba_num_obj-1];
        sba_rhs=new double [sba_num_obj-1];
        for (int j=0; j<sba_num_obj-1; j++) {sba_mathmodel_sense[j]='G';    sba_rhs[j]=-CPX_INFBOUND;}
        sba_cp_status = CPXnewrows (sba_env, sba_prob, sba_num_obj-1, sba_rhs , sba_mathmodel_sense, NULL, NULL);
        for (int j=0; j<sba_num_obj-1; j++){
            sba_cp_status=CPXchgcoef (sba_env, sba_prob, sba_num_mathmodel_row + j, j, 1);
            sba_cp_status=CPXchgcoef (sba_env, sba_prob, sba_num_mathmodel_row + j, sba_num_obj-1, -1);
        }
        
        // sba_cp_status = CPXwriteprob (sba_env, sba_prob, (sba_path+"coverage_prob.lp").c_str(), NULL);
        initialP = solveForTheInitialSolution();
        
        //==========================================================================
        //					initialize the tree
        //==========================================================================
        
        for (int j=1; j<=sba_num_obj; j++) {
            sba_global_lower_bound.push_back(sba_float_minus_infinity);
        }
        
        s.set_bound(sba_global_lower_bound);
        s.set_point(sba_global_lower_bound);
        s.set_sub_status(-1);
        s.set_mod_status(-1);
        
        sba_top=sba_spaceT.begin();
        sba_root=sba_spaceT.insert(sba_top, s); // insert the global lower bound vector to the root
        Pnode=sba_root;
        
        for (int j=0; j<sba_num_obj-2; j++) {
            Cnode=sba_spaceT.append_child(Pnode,s);
            Pnode=Cnode;
        }
        
        sba_gen_points.push_back(initialP);
        // update nadir point
        for(int j=0; j<sba_num_obj-1;j++) if(initialP[j]<sba_nadirPoint[j]) sba_nadirPoint[j]=initialP[j];
        sba_num_point=1;
        (*Cnode).set_point(initialP);
        (*Cnode).set_mod_status(1);
        (*Cnode).set_sub_status(1);
        
        sba_new_point=initialP;
        
        flag=1;
    }
    else
    {
        cout << "Model is infeasible\n";
        flag = 0;
        sba_solver_status = "ProblemInfeasible";
    }
    
    // check the stopping conditions
    if(( clock() - sba_start ) / (double) CLOCKS_PER_SEC >= timeLimit)
    {
        flag = 0;
        sba_solver_status = "AbortTimeLimit";
    } else if(sba_num_point>=pointLimit)
    {
        flag=0;
        sba_solver_status = "AbortPointLimit";
    }
    
    //==========================================================================
    // Generation loop
    //==========================================================================
    
    while (flag) {
        
        b=sba_global_lower_bound;
        b[0]=sba_new_point[0]+sba_epsilon;
        
        s.set_bound(b);
        
        //==================================================================================================================================
        // traverse the first level and append the child to the appropriate place. In this way, siblings will be always sorted at each level.
        //==================================================================================================================================
        sibIter=sba_spaceT.begin(sba_root);
        
        flag_insert=0;
        
        while (sibIter!=sba_spaceT.end(sba_root)) { // first level traversing====================
            
            
            if (sba_spaceT.number_of_siblings(sibIter)<=0) {
                
                sba_spaceT.insert_after(sibIter, s);
                sibIter++;
                sba_form_subtree(sibIter,1);
                flag_insert=1;
                break;
            }
            
            else if ((*sibIter).get_bound()[0]>=b[0]) {
                
                flag_insert=1;
                
                if ((*sibIter).get_bound()[0]!=b[0]) {
                    
                    sba_spaceT.insert(sibIter, s);
                    sibIter--;
                    sba_form_subtree(sibIter,1);
                    
                }
                
                break;
            }
            
            sibIter++;
            
        }
        
        if (flag_insert==0) {
            sibIter--;
            sibIter=sba_spaceT.insert_after(sibIter, s);
            sba_form_subtree(sibIter,1);
        }
        
        //==================================================================================================================================
        //update the subtrees of all siblings on the left of the new node
        //==================================================================================================================================
        sibIter_2=sba_spaceT.begin(sba_root);
        
        while (sibIter_2!=sibIter){
            sba_update_subtree(sibIter_2,1);
            sibIter_2++;
        }
        
        //==================================================================================================================================
        // pick the new generated point which has the max z_m value
        //==================================================================================================================================
        
        sba_new_point=sba_global_lower_bound;
        leafIter=sba_spaceT.begin_leaf();
        
        stop_flag=0;
        while (leafIter!=sba_spaceT.end_leaf()) {
            
            if ((*leafIter).get_sub_status()==1  && (*leafIter).get_point()[sba_num_obj-1]>sba_new_point[sba_num_obj-1]) {
                sba_new_point=(*leafIter).get_point();
                stop_flag=1;
            }
            leafIter++;
        }
        
        if (stop_flag==0) {
            flag = 0; // no more nondominated point
        }
        else
        {
            // add the new nondominated point to the list
            sba_gen_points.push_back(sba_new_point);
            // update nadir point
            for(int j=0; j<sba_num_obj-1;j++) if(sba_new_point[j]<sba_nadirPoint[j]) sba_nadirPoint[j]=sba_new_point[j];
            sba_num_point++;
        }
        
        // check the stopping conditions
        double elapsedTime = ( clock() - sba_start ) / (double) CLOCKS_PER_SEC;
        if ( elapsedTime >= timeLimit)
        {
            flag = 0; // time limit exceeded
            sba_solver_status = "AbortTimeLimit";
        } else if(sba_num_point>=pointLimit)
        {
            flag=0; // target number of points achieved.
            sba_solver_status = "AbortPointLimit";
        }
        
        cout << "point: " << sba_num_point << " model solved: " << sba_num_model_solved << " elapsed time " << elapsedTime <<  endl;
        
        //==================================================================================================================================
        //display the tree
        //==================================================================================================================================
        
        /*cout << "tree size: " << sba_spaceT.size() << "\n";
         
         Iter=sba_spaceT.begin(sba_root);
         while (Iter!=sba_spaceT.end(sba_root)) {
         
         if (sba_spaceT.depth(Iter)==sba_num_obj-2)
         {
         cout << "bound: ";
         for (int i=0; i<sba_num_obj-1; i++) {
         cout << (*Iter).get_bound()[i] << " ";
         }
         
         cout << "sub_status: " << (*Iter).get_sub_status() << " ";
         cout << "mod_status: " << (*Iter).get_mod_status() << " ";
         cout << "point: ";
         for (int i=0; i<sba_num_obj; i++) {
         cout << (*Iter).get_point()[i] << " ";
         }
         }
         
         cout << endl;
         
         Iter++;
         }
         
         cout << "new point: ";
         for (int i=0; i<sba_num_obj; i++) {
         
         cout << sba_new_point[i] << " ";
         }
         cout << "\n";
         */
        //========================================================================
        
    }
    
    sba_cpu_time = ( clock() - sba_start ) / (double) CLOCKS_PER_SEC;
    
    // display the generated points
    //cout << "Nondominated point list: \n";
    sba_file_ND_points << "#Solver status:" << endl << sba_solver_status << endl;
    sba_file_ND_points << "#Number of nondominated points:" << endl << sba_num_point << endl;
    sba_file_ND_points << "#Number of models solved:" << endl << sba_num_model_solved << endl;
    sba_file_ND_points << "#Elapsed time (seconds):" << endl << sba_cpu_time << endl;
    sba_file_ND_points << "#Ideal point:" << endl;
    for (int j=0; j<sba_num_obj-1; j++) sba_file_ND_points << sba_idealPoint[j]*criteriaScalCoeffs[j] << " ";
    sba_file_ND_points << endl;
    sba_file_ND_points << "#Incumbent nadir point:" << endl;
    for (int j=0; j<sba_num_obj-1; j++) sba_file_ND_points << sba_nadirPoint[j]*criteriaScalCoeffs[j] << " ";
    sba_file_ND_points << endl;
    
    // calculate maximum efficient range
    double scalingCoeff = 0.0;
    for(int j=0;j<sba_num_obj-1;j++){
        if(sba_idealPoint[j]-sba_nadirPoint[j]>scalingCoeff)
            scalingCoeff = sba_idealPoint[j]-sba_nadirPoint[j];
    }
    
    sba_file_ND_points << "#The set of nondominated points:" << endl;
    for (int i=0; i < sba_gen_points.size(); i++) {
        for (int j=0; j<sba_num_obj-1; j++) {
            sba_file_ND_points << sba_gen_points[i][j]*criteriaScalCoeffs[j] << " ";
        }
        if (i!=0) sba_file_ND_points << sba_gen_points[i][sba_num_obj-1]/scalingCoeff;
        sba_file_ND_points << endl;
    }
    sba_file_ND_points.close();
    CPXfreeprob(sba_env, &sba_prob);
    return 1;
    
}

void sba_solve_subspace (tree<subspace>::leaf_iterator sub)
{
    vector<float> ndp;
    vector<float> bound;
    int status=-1;
    
    
    for (int i=1; i<=sba_num_obj; i++) {
        ndp.push_back(sba_float_minus_infinity);
    }
    
    (*sub).set_mod_status(0);
    
    status=sba_search_subspace(sub,ndp);
    
    if (status==-1)
    {
        sba_num_model_solved++;
        (*sub).set_mod_status(1);
        bound=(*sub).get_bound();
        //==========================================================================
        //                           MODEL SOLUTION
        //==========================================================================
        
        for (int j=0; j<sba_num_obj-1; j++)
            sba_cp_status=CPXchgcoef (sba_env, sba_prob, sba_num_mathmodel_row + j, -1, bound[j]);
        
        // sba_cp_status = CPXwriteprob (sba_env, sba_prob, (sba_path+"myprob.lp").c_str(), NULL);
        
        sba_cp_status = CPXmipopt (sba_env, sba_prob); // solve the problem
        sba_cp_opt=CPXgetstat(sba_env,sba_prob);
        //cout << sba_cp_opt << endl;
        if (sba_cp_opt==101 || sba_cp_opt==102)
        {
            sba_cp_status = CPXgetx (sba_env, sba_prob, sba_mathmodel_sol, 0, sba_num_obj-1); // copies the optimal obj values to the "sba_mathmodel_sol"
            status=1;
            for (int j=0; j<sba_num_obj; j++) ndp[j]=(float)sba_mathmodel_sol[j];
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
     for (int i=0; i<sba_num_obj; i++) {
     cout << ndp[i] << " ";
     }
     cout << endl;*/
    
    (*sub).set_sub_status(status);
    (*sub).set_point(ndp);
    
    
}

void sba_update_subtree(tree<subspace>::sibling_iterator p_it, int d)
{
    vector<float> b;
    int flag,last_sib;
    tree<subspace>::sibling_iterator s_it,s_newP,next_sib;
    subspace s;
    
    s.set_bound(sba_global_lower_bound);
    s.set_point(sba_global_lower_bound);
    s.set_sub_status(-1);
    s.set_mod_status(-1);
    
    if (d<sba_num_obj-2) {
        b=(*p_it).get_bound();
        
        flag=0;
        for (int i=0; i<d; i++) {
            if (b[i] > sba_new_point[i]) { // which means the new_point is not in this subspace, no need to update
                flag=1;
                break;
            }
        }
        b[d]=sba_new_point[d]+sba_epsilon;
        s.set_bound(b);
        
        if (flag==0) { // update the children
            
            s_it=sba_spaceT.begin(p_it); // points to the first child
            
            if (sba_spaceT.number_of_siblings(s_it)<=0) {
                sba_spaceT.insert_after(s_it, s);
                s_it++;
                sba_form_subtree(s_it, d+1);
                s_it--;
                sba_update_subtree(s_it, d+1);
            }
            else {
                // Step 1. Insert the node to the appropriate place.
                last_sib=1;
                int insertion_index = 0;
                while (sba_spaceT.index(s_it)<sba_spaceT.number_of_children(p_it)-1) {
                    
                    next_sib=sba_spaceT.next_sibling(s_it);
                    
                    // check if we need to insert a new node after this one
                    if ( (*s_it).get_bound()[d] <= b[d])
                    {
                        if ((*next_sib).get_bound()[d] >= b[d])
                        {
                            last_sib=0; // the new node will not be the last sibling
                            insertion_index = sba_spaceT.index(s_it);
                            
                            if (b[d]!=(*s_it).get_bound()[d] && b[d]!=(*next_sib).get_bound()[d]) {
                                s_newP=sba_spaceT.insert_after(s_it, s);
                                sba_form_subtree(s_newP, d+1);
                            }
                            break;
                        }
                    }
                    s_it++;
                }
                if (last_sib==1) { // if the new node is to be inserted to the end
                    insertion_index = sba_spaceT.index(s_it);
                    s_newP=sba_spaceT.insert_after(s_it, s);
                    sba_form_subtree(s_newP, d+1);
                }
                
                // Step 2. Update all subtrees at the left-hand side of the inserted node
                s_it=sba_spaceT.begin(p_it);
                while (sba_spaceT.index(s_it) <= insertion_index) {
                    sba_update_subtree(s_it, d+1);
                    s_it++;
                }
            }
        }
    }
    else // this means that it is a leaf node
    {
        b=(*p_it).get_bound();
        flag=0;
        for (int i=0; i<sba_num_obj; i++) {
            if (b[i] > sba_new_point[i]) { // which means the new_point is not in this subspace, no need to solve this leaf node
                flag=1;
                break;
            }
        }
        
        if (flag==0) // solve the leaf
        {
            sba_set_lastBound(p_it);
            sba_solve_subspace(p_it);
        }
        
    }
}

void sba_form_subtree(tree<subspace>::sibling_iterator p_it, int d)
{
    vector<float> b,p;
    tree<subspace>::sibling_iterator sibStart_it,sibEnd_it;
    tree<subspace>::iterator s_it;
    subspace s;
    
    s.set_bound(sba_global_lower_bound);
    s.set_point(sba_global_lower_bound);
    s.set_sub_status(-1);
    s.set_mod_status(-1);
    
    if (d < sba_num_obj-2) { // NOT a leaf node
        
        b=sba_global_lower_bound;
        
        // copy the first d bound of the parent node
        for (int i=0; i < d; i++) {
            b[i]=(*p_it).get_bound()[i];
        }
        
        if (sba_spaceT.number_of_siblings(p_it) > sba_spaceT.index(p_it)) // if parent is not the last sibling
        {
            // copy the left subtree
            sibStart_it = sba_spaceT.begin(sba_spaceT.previous_sibling(p_it));
            sibEnd_it = sba_spaceT.end(sba_spaceT.previous_sibling(p_it));
            sba_spaceT.append_children(p_it, sibStart_it, sibEnd_it);
            
            s_it = sba_spaceT.begin(p_it);
            // update the bounds and leaf node info
            while (s_it!=sba_spaceT.end(p_it)){
                for (int i=d; i < sba_num_obj-2; i++) b[i]=(*s_it).get_bound()[i];
                (*s_it).set_bound(b);
                if(sba_spaceT.depth(s_it)==(sba_num_obj-2)){
                    sba_set_lastBound(s_it);
                    sba_solve_subspace(s_it);
                }
                s_it++;
            }
        }
        else // if parent is the last sibling
        {
            for (int i=d; i < sba_num_obj-2; i++)
                b[i]=sba_float_minus_infinity;
            s.set_bound(b);
            s_it=sba_spaceT.append_child(p_it,s);
            sba_form_subtree(s_it, d+1);
        }
    }
    else { // this means that it is a leaf node
        
        // we need to update the bound in the z_(m-1) objective
        sba_set_lastBound(p_it);
        
        // find the nondominated point having the maximum z_m value in this subspace
        sba_solve_subspace(p_it);
        
    }
    
}


void sba_set_lastBound(tree<subspace>::sibling_iterator it)
{
    int flag;
    vector<float> b;
    float max_z_m=sba_global_lower_bound[sba_num_obj-2];
    
    
    b=(*it).get_bound();
    
    
    for (int j=0; j<sba_gen_points.size(); j++) {
        
        flag=1;
        for (int k=0; k<sba_num_obj-2; k++) {
            if (b[k] > sba_gen_points[j][k]) {
                flag=0;
                break;
            }
        }
        
        if ( (flag==1) && (sba_gen_points[j][sba_num_obj-2] > max_z_m)) {
            
            max_z_m=sba_gen_points[j][sba_num_obj-2];
        }
        
    }
    
    b[sba_num_obj-2]=max_z_m+sba_epsilon;
    (*it).set_bound(b);
    
}

int sba_search_subspace (tree<subspace>::leaf_iterator sub, vector<float>& p)
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
     for (int i=0; i<sba_num_obj-1; i++) {
     cout << (*sibIter).get_bound()[i] << " ";
     }
     cout << endl;*/
    //================================================================
    
    //cout << "========SCANNED NODES=======" << endl;
    
    if (sba_spaceT.depth(sibIter) == 1)
    {
        if (sba_spaceT.index(sibIter)!= 0)
        {
            prev_sib=sba_spaceT.previous_sibling(sibIter); // step back to the previous sibling
            
            // display the subspaces================================================
            /*cout << "bound: ";
             for (int i=0; i<sba_num_obj-1; i++) {
             cout << (*prev_sib).get_bound()[i] << " ";
             }
             
             cout << "sub_status: " << (*prev_sib).get_sub_status() << " ";
             
             if ((*prev_sib).get_sub_status()==1)
             {
             for (int i=0; i<sba_num_obj; i++) cout << (*prev_sib).get_point()[i] << " ";
             }
             
             cout << endl;*/
            
            //=======================================================================
            
            // compare with the nodes
            
            // bounds of the current subspace must be greater than or equal to each bound of the compared subspace
            flag=1;
            for (int i=0; i<sba_num_obj-1; i++)
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
                    for (int i=0; i<sba_num_obj-1; i++)
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
        
        while(sba_spaceT.depth(sibIter) != 1) // go up the tree until reaching to the first level of the tree
        {
            sibIter=sba_spaceT.parent(sibIter);
        }
        
        
        if (sba_spaceT.index(sibIter)!= 0)
        {
            
            prev_sib=sba_spaceT.previous_sibling(sibIter); // step back to the previous sibling
            
            leafIter=sba_spaceT.begin_leaf(prev_sib);
            
            // traversing previous largest subtree
            
            while (leafIter!=sba_spaceT.end_leaf(prev_sib) && status==-1)
            {
                
                // display the subspaces================================================
                /*cout << "bound: ";
                 for (int i=0; i<sba_num_obj-1; i++) {
                 cout << (*leafIter).get_bound()[i] << " ";
                 }
                 
                 cout << "sub_status: " << (*leafIter).get_sub_status() << " ";
                 cout << "mod_status: " << (*leafIter).get_mod_status() << " ";
                 
                 if ((*leafIter).get_sub_status()==1)
                 {
                 for (int i=0; i<sba_num_obj; i++) cout << (*leafIter).get_point()[i] << " ";
                 }
                 
                 cout << endl;*/
                
                //=======================================================================
                
                // compare the nodes
                
                // bounds of the current subspace must be greater than or equal to each bound of the compared subspace
                flag=1;
                for (int i=0; i<sba_num_obj-1; i++)
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
                        for (int i=0; i<sba_num_obj-1; i++)
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
        leafIter=sba_spaceT.begin_leaf(sibIter);
        
        while (leafIter!=sub && status==-1)
        {
            
            // display the subspaces======================================
            /*cout << "bound: ";
             for (int i=0; i<sba_num_obj-1; i++) {
             cout << (*leafIter).get_bound()[i] << " ";
             }
             
             cout << "sub_status: " << (*leafIter).get_sub_status() << " ";
             cout << "mod_status: " << (*leafIter).get_mod_status() << " ";
             
             if ((*leafIter).get_sub_status()==1)
             {
             for (int i=0; i<sba_num_obj; i++) cout << (*leafIter).get_point()[i] << " ";
             }
             
             cout << endl;*/
            
            //=============================================================
            
            // compare the nodes
            
            // bounds of the current subspace must be greater than or equal to each bound of the compared subspace
            flag=1;
            for (int i=0; i<sba_num_obj-1; i++)
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
                    for (int i=0; i<sba_num_obj-1; i++)
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

vector<double> calculateCriteriaScalingCoeffs(){
    
    vector<double> coeffs(sba_num_obj);
    double coeff;
    
    int j = 0;
    for(int r=sba_num_mathmodel_row-sba_num_obj+1;r<sba_num_mathmodel_row;r++){
        coeffs[j]=0.0;
        for(int c=sba_num_obj;c<sba_num_mathmodel_col;c++){
            CPXgetcoef (sba_env, sba_prob, r, c, &coeff);
            coeffs[j]+=pow(coeff,2);
        }
        coeffs[j]=pow(coeffs[j],0.5);
        j++;
    }
    coeffs[sba_num_obj-1]=1.0;
    
    return coeffs;
}

void generatePayoffTable(){
    for(int j=0; j<sba_num_obj-1;j++){
        sba_idealPoint.push_back(-numeric_limits<float>::infinity( ));
        sba_nadirPoint.push_back(numeric_limits<float>::infinity( ));
    }
    vector<float> point;
    for(int j=0; j<sba_num_obj-1;j++){
        point = solveSingleObjectiveProblem(j);
        if(sba_solver_status == "ProblemInfeasible") break;
        for(int k=0; k<sba_num_obj-1;k++){
            if(point[k]>sba_idealPoint[k]) sba_idealPoint[k]=point[k];
            if(point[k]<sba_nadirPoint[k]) sba_nadirPoint[k]=point[k];
        }
        sba_cp_status=CPXchgcoef (sba_env, sba_prob, -1, j, 0.0001);
    }
    sba_cp_status=CPXchgcoef (sba_env, sba_prob, -1, sba_num_obj-1, 1.0);
}

vector<float> solveSingleObjectiveProblem(int objIndex){
    vector<float> point(sba_num_obj-1);
    sba_cp_status=CPXchgcoef (sba_env, sba_prob, -1, sba_num_obj-1, 0.0);
    sba_cp_status=CPXchgcoef (sba_env, sba_prob, -1, objIndex, 1.0);
    sba_mathmodel_sol = new double [sba_num_obj];
    sba_cp_status = CPXmipopt (sba_env, sba_prob); // solve the problem
    sba_cp_opt=CPXgetstat(sba_env,sba_prob);
    
    if (sba_cp_opt==101 || sba_cp_opt==102)
    {
        sba_cp_status = CPXgetx (sba_env, sba_prob, sba_mathmodel_sol, 0, sba_num_obj-1); // copies the optimal obj values to the "sba_mathmodel_sol"
        for (int j=0; j<sba_num_obj; j++) point[j]= (float)sba_mathmodel_sol[j];
    } else {
        sba_solver_status = "ProblemInfeasible";
        sba_cp_status = CPXwriteprob (sba_env, sba_prob, (sba_path+"infeasible_prob.lp").c_str(), NULL);
    }
    return point;
}

vector<float> solveForTheInitialSolution(){
    // call this function after finding the ideal point and constructing the coverega constraints
    vector<float> point(sba_num_obj-1);
    
    // modify the model to obtain achievement scalarizing program
    sba_cp_status=CPXchgcoef (sba_env, sba_prob, -1, sba_num_obj-1, -1);
    for (int j=0; j<sba_num_obj-1; j++){
        sba_cp_status=CPXchgcoef (sba_env, sba_prob, sba_num_mathmodel_row + j, sba_num_obj-1, 1);
        sba_cp_status=CPXchgcoef (sba_env, sba_prob, sba_num_mathmodel_row + j, -1, sba_idealPoint[j]);
    }
    // sba_cp_status = CPXwriteprob (sba_env, sba_prob, (sba_path+"asp.lp").c_str(), NULL);
    
    sba_mathmodel_sol = new double [sba_num_obj];
    sba_cp_status = CPXmipopt (sba_env, sba_prob); // solve the problem
    sba_num_model_solved++;
    sba_cp_opt=CPXgetstat(sba_env,sba_prob);
    
    if (sba_cp_opt==101 || sba_cp_opt==102)
    {
        sba_cp_status = CPXgetx (sba_env, sba_prob, sba_mathmodel_sol, 0, sba_num_obj-1); // copies the optimal obj values to the "sba_mathmodel_sol"
        for (int j=0; j<sba_num_obj; j++) point[j]= (float)sba_mathmodel_sol[j];
    } else {
        sba_solver_status = "ProblemInfeasible";
    }
    
    // reset the model
    sba_cp_status=CPXchgcoef (sba_env, sba_prob, -1, sba_num_obj-1, 1);
    for (int j=0; j<sba_num_obj-1; j++){
        sba_cp_status=CPXchgcoef (sba_env, sba_prob, sba_num_mathmodel_row + j, sba_num_obj-1, -1);
    }
    
    return point;
}

double calculateCoverageGap(int n){
    double alpha_min = 1.0;
    double alpha = 0.0;
    double scalingCoeff = 0.0;
    for(int j=0;j<sba_num_obj-1;j++){
        if(sba_idealPoint[j]-sba_nadirPoint[j]>scalingCoeff)
            scalingCoeff = sba_idealPoint[j]-sba_nadirPoint[j];
    }
    
    for(int i=0;i<n;i++){
        alpha = 0.0;
        for(int j=0;j<sba_num_obj-1;j++){
            if((sba_gen_points[n][j]-sba_gen_points[i][j])/scalingCoeff>alpha)
                alpha = (sba_gen_points[n][j]-sba_gen_points[i][j])/scalingCoeff;
        }
        if (alpha < alpha_min) alpha_min = alpha;
    }
    return alpha_min;
}
