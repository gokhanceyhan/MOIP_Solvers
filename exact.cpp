//
//  exact.cpp
//  MO_solver
//
//  Created by Gokhan Ceyhan on 9/26/17.
//  Copyright © 2017 Gokhan Ceyhan. All rights reserved.
//


#include "exact.h"
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
#include "constants.h"
#include "math.h"

using namespace std;

//======= auxilary functions===============================================
void form_subtree(tree<subspace>::sibling_iterator, int);
void update_subtree(tree<subspace>::sibling_iterator, int);
void set_lastBound(tree<subspace>::sibling_iterator);
void solve_subspace (tree<subspace>::leaf_iterator); // finds the nondominated point with max z_m value in the subspace
int search_subspace (tree<subspace>::leaf_iterator, vector<float>& p); // search the minimum volume subspaces that contain the given subspace
vector<double> calculateCriteriaScalingCoeffs(); // scales the criteron vectors by dividing them to their Euclidean norms

//======= general global variables==========================================
int num_obj; // number of objectives
int num_point; // number of nondominated points generated
int num_model_solved=0; // number of times cplex is called
string solver_status="Success"; // the status of the multi-objective solver
vector<float> idealPoint;
vector<float> nadirPoint; // approximate nadir unless complete nondominated point set is known.

tree<subspace> spaceT ; // tree of generated subspaces
tree<subspace>::iterator top,root;
vector<vector<float>> points; // set of nondominated points read from the input file
vector<vector<float>> gen_points; // generated nondominated points
vector<float> new_point; // last nondominated point generated
float float_minus_infinity=-numeric_limits<float>::infinity( ); // set to minus infinity for type "float".
vector<float> global_lower_bound; // global lower bound vector at the beginning of the algorithm
double epsilon; // the amount by which the bounds are shifted.



// ==== cplex and math model variables ======================================
CPXENVptr env=NULL;
CPXLPptr prob=NULL;
int cp_status;
int cp_opt; // stores the status of the model solved (for optimality 101 or 102 must be returned)
int num_mathmodel_col;  // number of columns in the math model
int num_mathmodel_row;  // number of rows in the math model
double *mathmodel_sol;  // points to the array storing the criterion values of the optiaml solution
double *rhs; // points to the array storing the initial right hand side values for bound constraints
char *mathmodel_sense; // points to the array storing the type of ineqaulities for bound constraints

//======= input/output variables ============================================
ofstream file_ND_points;
clock_t start;
double cpu_time;

int exact(int m, string path, double timeLimit)
{
    start=clock();
    num_obj=m;
    epsilon = BOUND_TOLERANCE;
    
    int flag,stop_flag; // zero if there is no new nondominated point
    int flag_insert; // used to indicate whether the new point is inserted at the first level
    subspace s;
    vector<float> b,initialP;
    
    tree<subspace>::iterator Iter,Pnode,Cnode; //position of the root of the tree
    tree<subspace>::sibling_iterator sibIter,sibIter_2;
    tree<subspace>::leaf_iterator leafIter;
    
    
    
    file_ND_points.open((path+"Output_nd_points.txt").c_str());
    
    //==========================================================================
    //          CPLEX INITIALIZATION and GENERATE FIRST NONDOMINATED POINT
    //==========================================================================
    
    env = CPXopenCPLEX(&cp_status);
    cp_status = CPXsetdefaults (env);
    cp_status = CPXsetintparam(env, CPX_PARAM_SCRIND, 0);
    cp_status = CPXsetintparam(env, CPXPARAM_MIP_Display, 0);
    cp_status = CPXsetintparam(env, CPX_PARAM_CLOCKTYPE, 2); // wall clock time
    cp_status = CPXsetdblparam(env, CPX_PARAM_TILIM, timeLimit); // sets the given time limit
    cp_status = CPXsetdblparam (env, CPX_PARAM_EPGAP, MIP_RELGAP);
    cp_status = CPXsetdblparam (env, CPX_PARAM_EPOPT, SIMPLEX_OPTGAP);
    cp_status = CPXsetdblparam (env, CPX_PARAM_BAREPCOMP, BARRIER_CONV);
    
    prob = CPXcreateprob (env, &cp_status, "mathmodel"); // create problem in CPLEX
    cp_status = CPXreadcopyprob (env, prob, (path+"model.lp").c_str(), NULL); // read "model.lp" into CPLEX
    num_mathmodel_col = CPXgetnumcols (env, prob); // get the number of variables in the math model
    num_mathmodel_row = CPXgetnumrows (env, prob); // get the number of constraints in the math model
    
    // calculate criteria scaling coefficients
    vector<double> criteriaScalCoeffs = calculateCriteriaScalingCoeffs();
    
    // update criterion constraint coefficients for scaling
    int j=0;
    for(int r=num_mathmodel_row-num_obj;r<num_mathmodel_row;r++){
        cp_status=CPXchgcoef (env, prob, r, j, -criteriaScalCoeffs[j]);
        j++;
    }
    
    // initialize ideal and nadir point vectors
    for(int j=0; j<num_obj;j++){
        idealPoint.push_back(-numeric_limits<float>::infinity( ));
        nadirPoint.push_back(numeric_limits<float>::infinity( ));
    }
    
    // add bound constraints
    mathmodel_sense = new char [num_obj-1];
    rhs=new double [num_obj-1];
    for (int j=0; j<num_obj-1; j++) {mathmodel_sense[j]='G';	rhs[j]=-CPX_INFBOUND;}
    cp_status = CPXnewrows (env, prob, num_obj-1, rhs , mathmodel_sense, NULL, NULL);
    for (int j=0; j<num_obj-1; j++)
        cp_status=CPXchgcoef (env, prob, num_mathmodel_row + j, j, 1);
    
    cp_status = CPXwriteprob (env, prob, (path+"myprob.lp").c_str(), NULL);
    
    mathmodel_sol = new double [num_obj];
    cp_status = CPXmipopt (env, prob); // solve the problem
    cp_opt=CPXgetstat(env,prob);
    
    if (cp_opt==101 || cp_opt==102)
    {
        
        cp_status = CPXgetx (env, prob, mathmodel_sol, 0, num_obj-1); // copies the optimal obj values to the "mathmodel_sol"
        
        for (int j=0; j<num_obj; j++)
        {
            initialP.push_back(0);
        }
        
        for (int j=0; j<num_obj; j++)
        {
            initialP[j]= (float)mathmodel_sol[j];
        }
        
        
        //==========================================================================
        //					initialize the tree
        //==========================================================================
        
        for (int j=1; j<=num_obj; j++) {
            global_lower_bound.push_back(float_minus_infinity);
        }
        
        s.set_bound(global_lower_bound);
        s.set_point(global_lower_bound);
        s.set_sub_status(-1);
        s.set_mod_status(-1);
        
        top=spaceT.begin();
        root=spaceT.insert(top, s); // insert the global lower bound vector to the root
        Pnode=root;
        
        for (int j=0; j<num_obj-2; j++) {
            Cnode=spaceT.append_child(Pnode,s);
            Pnode=Cnode;
        }
        
        gen_points.push_back(initialP);
        // update ideal and nadir points
        for(int j=0; j<num_obj;j++) if(initialP[j]>idealPoint[j]) idealPoint[j]=initialP[j];
        for(int j=0; j<num_obj;j++) if(initialP[j]<nadirPoint[j]) nadirPoint[j]=initialP[j];
        num_point=1;
        (*Cnode).set_point(initialP);
        (*Cnode).set_mod_status(1);
        (*Cnode).set_sub_status(1);
        
        new_point=initialP;
        
        flag=1;
    }
    else
    {
        cout << "Failed to solve single objective problem.\n";
        flag = 0;
        solver_status = "FailedToSolveSingleObjProblem";
    }
    
    // check the stopping conditions
    if(( clock() - start ) / (double) CLOCKS_PER_SEC >= timeLimit)
    {
        flag = 0;
        solver_status = "AbortTimeLimit";
    }
    
    //==========================================================================
    // Generation loop
    //==========================================================================
    
    while (flag) {
        
        b=global_lower_bound;
        b[0]=new_point[0]+epsilon;
        
        s.set_bound(b);
        
        //==================================================================================================================================
        // traverse the first level and append the child to the appropriate place. In this way, siblings will be always sorted at each level.
        //==================================================================================================================================
        sibIter=spaceT.begin(root);
        
        flag_insert=0;
        
        while (sibIter!=spaceT.end(root)) { // first level traversing====================
            
            if (spaceT.number_of_siblings(sibIter)<=0) {
                spaceT.insert_after(sibIter, s);
                sibIter++;
                form_subtree(sibIter,1);
                flag_insert=1;
                break;
            }
            
            else if ((*sibIter).get_bound()[0]>=b[0]) {
                flag_insert=1;
                if ((*sibIter).get_bound()[0]!=b[0]) {
                    spaceT.insert(sibIter, s);
                    sibIter--;
                    form_subtree(sibIter,1);
                }
                break;
            }
            sibIter++;
            
        }
        
        if (flag_insert==0) {
            sibIter--;
            sibIter=spaceT.insert_after(sibIter, s);
            form_subtree(sibIter,1);
        }
        
        //==================================================================================================================================
        //update the subtrees of all siblings on the left of the new node
        //==================================================================================================================================
        sibIter_2=spaceT.begin(root);
        
        while (sibIter_2!=sibIter){
            update_subtree(sibIter_2,1);
            sibIter_2++;
        }
        
        //==================================================================================================================================
        // pick the new generated point which has the max z_m value
        //==================================================================================================================================
        
        new_point=global_lower_bound;
        leafIter=spaceT.begin_leaf();
        
        stop_flag=0;
        while (leafIter!=spaceT.end_leaf()) {
            
            if ((*leafIter).get_sub_status()==1  && (*leafIter).get_point()[num_obj-1]>new_point[num_obj-1]) {
                new_point=(*leafIter).get_point();
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
            gen_points.push_back(new_point);
            // update ideal and nadir points
            for(int j=0; j<num_obj;j++) if(new_point[j]>idealPoint[j]) idealPoint[j]=new_point[j];
            for(int j=0; j<num_obj;j++) if(new_point[j]<nadirPoint[j]) nadirPoint[j]=new_point[j];
            num_point++;
        }
        
        // check the stopping conditions
        double elapsedTime = ( clock() - start ) / (double) CLOCKS_PER_SEC;
        if (elapsedTime >= timeLimit)
        {
            flag = 0; // time limit exceeded
            solver_status = "AbortTimeLimit";
        }
        
        cout << "point: " << num_point << " model solved: " << num_model_solved << " elapsed time " << elapsedTime <<  endl;
        
        //==================================================================================================================================
        //display the tree
        //==================================================================================================================================
        
        /*cout << "tree size: " << spaceT.size() << "\n";
         
         Iter=spaceT.begin(root);
         while (Iter!=spaceT.end(root)) {
         
         if (spaceT.depth(Iter)==num_obj-2)
         {
         cout << "bound: ";
         for (int i=0; i<num_obj-1; i++) {
         cout << (*Iter).get_bound()[i] << " ";
         }
         
         cout << "sub_status: " << (*Iter).get_sub_status() << " ";
         cout << "mod_status: " << (*Iter).get_mod_status() << " ";
         cout << "point: ";
         for (int i=0; i<num_obj; i++) {
         cout << (*Iter).get_point()[i] << " ";
         }
         }
         
         cout << endl;
         
         Iter++;
         }
         
         cout << "new point: ";
         for (int i=0; i<num_obj; i++) {
         
         cout << new_point[i] << " ";
         }
         cout << "\n";
         */
        //========================================================================
        
    }
    
    cpu_time = ( clock() - start ) / (double) CLOCKS_PER_SEC;
    
    // display the generated points
    //cout << "Nondominated point list: \n";
    file_ND_points << "#Solver type:" << endl << "nMOCO-S" << endl;
    file_ND_points << "#Solver status:" << endl << solver_status << endl;
    if (solver_status != "FailedToSolveSingleObjProblem") {
        file_ND_points << "#Number of nondominated points:" << endl << num_point << endl;
        file_ND_points << "#Number of models solved:" << endl << num_model_solved << endl;
        file_ND_points << "#Elapsed time (seconds):" << endl << cpu_time << endl;
        file_ND_points << "#Incumbent ideal point:" << endl;
        for (int j=0; j<num_obj; j++) file_ND_points << idealPoint[j]*criteriaScalCoeffs[j] << " ";
        file_ND_points << endl;
        file_ND_points << "#Incumbent nadir point:" << endl;
        for (int j=0; j<num_obj; j++) file_ND_points << nadirPoint[j]*criteriaScalCoeffs[j] << " ";
        file_ND_points << endl;
        file_ND_points << "#The set of nondominated points:" << endl;
        for (int i=0; i < gen_points.size(); i++) {
            for (int j=0; j<num_obj; j++) {
                file_ND_points << gen_points[i][j]*criteriaScalCoeffs[j] << " ";
            }
            file_ND_points << "\n";
        }
    } else {
        file_ND_points << "#Cplex error code:" << endl << cp_status << endl;
        file_ND_points << "#Cplex optimization status:" << endl << cp_opt << endl;
    }
    
    file_ND_points.close();
    CPXfreeprob(env, &prob);
    return 1;
    
}

void solve_subspace (tree<subspace>::leaf_iterator sub)
{
    vector<float> ndp;
    vector<float> bound;
    int status=-1;
    
    
    for (int i=1; i<=num_obj; i++) {
        ndp.push_back(float_minus_infinity);
    }
    
    (*sub).set_mod_status(0);
    
    status=search_subspace(sub,ndp);
    
    if (status==-1)
    {
        num_model_solved++;
        (*sub).set_mod_status(1);
        bound=(*sub).get_bound();
        //==========================================================================
        //                           MODEL SOLUTION
        //==========================================================================
        
        for (int j=0; j<num_obj-1; j++)
            cp_status=CPXchgcoef (env, prob, num_mathmodel_row + j, -1, bound[j]);
        
        //cp_status = CPXwriteprob (env, prob, "myprob.lp", NULL);
        
        cp_status = CPXmipopt (env, prob); // solve the problem
        cp_opt=CPXgetstat(env,prob);
        
        if (cp_opt==101 || cp_opt==102)
        {
            cp_status = CPXgetx (env, prob, mathmodel_sol, 0, num_obj-1); // copies the optimal obj values to the "mathmodel_sol"
            status=1;
            for (int j=0; j<num_obj; j++) ndp[j]=(float)mathmodel_sol[j];
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
     for (int i=0; i<num_obj; i++) {
     cout << ndp[i] << " ";
     }
     cout << endl;*/
    
    (*sub).set_sub_status(status);
    (*sub).set_point(ndp);
    
    
}

void update_subtree(tree<subspace>::sibling_iterator p_it, int d)
{
    vector<float> b;
    int flag,last_sib;
    tree<subspace>::sibling_iterator s_it,s_newP,next_sib;
    subspace s;
    
    s.set_bound(global_lower_bound);
    s.set_point(global_lower_bound);
    s.set_sub_status(-1);
    s.set_mod_status(-1);
    
    if (d<num_obj-2) {
        b=(*p_it).get_bound();
        
        flag=0;
        for (int i=0; i<d; i++) {
            if (b[i] > new_point[i]) { // which means the new_point is not in this subspace, no need to update
                flag=1;
                break;
            }
        }
        b[d]=new_point[d]+epsilon;
        s.set_bound(b);
        
        if (flag==0) { // update the children
            
            s_it=spaceT.begin(p_it); // points to the first child
            
            if (spaceT.number_of_siblings(s_it)<=0) {
                spaceT.insert_after(s_it, s);
                s_it++;
                form_subtree(s_it, d+1);
                s_it--;
                update_subtree(s_it, d+1);
            }
            else {
                // Step 1. Insert the node to the appropriate place.
                last_sib=1;
                int insertion_index = 0;
                while (spaceT.index(s_it)<spaceT.number_of_children(p_it)-1) {
                    
                    next_sib=spaceT.next_sibling(s_it);
                    
                    // check if we need to insert a new node after this one
                    if ( (*s_it).get_bound()[d] <= b[d])
                    {
                        if ((*next_sib).get_bound()[d] >= b[d])
                        {
                            last_sib=0; // the new node will not be the last sibling
                            insertion_index = spaceT.index(s_it);
                            
                            if (b[d]!=(*s_it).get_bound()[d] && b[d]!=(*next_sib).get_bound()[d]) {
                                s_newP=spaceT.insert_after(s_it, s);
                                form_subtree(s_newP, d+1);
                            }
                            break;
                        }
                    }
                    s_it++;
                }
                if (last_sib==1) { // if the new node is to be inserted to the end
                    insertion_index = spaceT.index(s_it);
                    s_newP=spaceT.insert_after(s_it, s);
                    form_subtree(s_newP, d+1);
                }
                
                // Step 2. Update all subtrees at the left-hand side of the inserted node
                s_it=spaceT.begin(p_it);
                while (spaceT.index(s_it) <= insertion_index) {
                    update_subtree(s_it, d+1);
                    s_it++;
                }
            }
        }
    }
    else // this means that it is a leaf node
    {
        b=(*p_it).get_bound();
        flag=0;
        for (int i=0; i<num_obj; i++) {
            if (b[i] > new_point[i]) { // which means the new_point is not in this subspace, no need to solve this leaf node
                flag=1;
                break;
            }
        }
        
        if (flag==0) // solve the leaf
        {
            set_lastBound(p_it);
            solve_subspace(p_it);
        }
        
    }
}

void form_subtree(tree<subspace>::sibling_iterator p_it, int d)
{
    vector<float> b,p;
    tree<subspace>::sibling_iterator sibStart_it,sibEnd_it;
    tree<subspace>::iterator s_it;
    subspace s;
    
    s.set_bound(global_lower_bound);
    s.set_point(global_lower_bound);
    s.set_sub_status(-1);
    s.set_mod_status(-1);
    
    if (d < num_obj-2) { // NOT a leaf node
        
        b=global_lower_bound;
        
        // copy the first d bound of the parent node
        for (int i=0; i < d; i++) {
            b[i]=(*p_it).get_bound()[i];
        }
        
        if (spaceT.number_of_siblings(p_it) > spaceT.index(p_it)) // if parent is not the last sibling
        {
            // copy the left subtree
            sibStart_it = spaceT.begin(spaceT.previous_sibling(p_it));
            sibEnd_it = spaceT.end(spaceT.previous_sibling(p_it));
            spaceT.append_children(p_it, sibStart_it, sibEnd_it);
            
            s_it = spaceT.begin(p_it);
            // update the bounds and leaf node info
            while (s_it!=spaceT.end(p_it)){
                for (int i=d; i < num_obj-2; i++) b[i]=(*s_it).get_bound()[i];
                (*s_it).set_bound(b);
                if(spaceT.depth(s_it)==(num_obj-2)){
                    set_lastBound(s_it);
                    solve_subspace(s_it);
                }
                s_it++;
            }
        }
        else // if parent is the last sibling
        {
            for (int i=d; i < num_obj-2; i++)
                b[i]=float_minus_infinity;
            s.set_bound(b);
            s_it=spaceT.append_child(p_it,s);
            form_subtree(s_it, d+1);
        }
    }
    else { // this means that it is a leaf node
        
        // we need to update the bound in the z_(m-1) objective
        set_lastBound(p_it);
        
        // find the nondominated point having the maximum z_m value in this subspace
        solve_subspace(p_it);
        
    }
    
}

void set_lastBound(tree<subspace>::sibling_iterator it)
{
    int flag;
    vector<float> b;
    float max_z_m=float_minus_infinity;
    
    
    b=(*it).get_bound();
    
    
    for (int j=0; j<gen_points.size(); j++) {
        
        flag=1;
        for (int k=0; k<num_obj-2; k++) {
            if (b[k] > gen_points[j][k]) {
                flag=0;
                break;
            }
        }
        
        if ( (flag==1) && (gen_points[j][num_obj-2] > max_z_m)) {
            
            max_z_m=gen_points[j][num_obj-2];
        }
        
    }
    
    b[num_obj-2]=max_z_m+epsilon;
    (*it).set_bound(b);
    
}

int search_subspace (tree<subspace>::leaf_iterator sub, vector<float>& p)
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
     for (int i=0; i<num_obj-1; i++) {
     cout << (*sibIter).get_bound()[i] << " ";
     }
     cout << endl;*/
    //================================================================
    
    //cout << "========SCANNED NODES=======" << endl;
    
    if (spaceT.depth(sibIter) == 1)
    {
        if (spaceT.index(sibIter)!= 0)
        {
            prev_sib=spaceT.previous_sibling(sibIter); // step back to the previous sibling
            
            // display the subspaces================================================
            /*cout << "bound: ";
             for (int i=0; i<num_obj-1; i++) {
             cout << (*prev_sib).get_bound()[i] << " ";
             }
             
             cout << "sub_status: " << (*prev_sib).get_sub_status() << " ";
             
             if ((*prev_sib).get_sub_status()==1)
             {
             for (int i=0; i<num_obj; i++) cout << (*prev_sib).get_point()[i] << " ";
             }
             
             cout << endl;*/
            
            //=======================================================================
            
            // compare with the nodes
            
            // bounds of the current subspace must be greater than or equal to each bound of the compared subspace
            flag=1;
            for (int i=0; i<num_obj-1; i++)
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
                    for (int i=0; i<num_obj-1; i++)
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
        
        while(spaceT.depth(sibIter) != 1) // go up the tree until reaching to the first level of the tree
        {
            sibIter=spaceT.parent(sibIter);
        }
        
        
        if (spaceT.index(sibIter)!= 0)
        {
            
            prev_sib=spaceT.previous_sibling(sibIter); // step back to the previous sibling
            
            leafIter=spaceT.begin_leaf(prev_sib);
            
            // traversing previous largest subtree
            
            while (leafIter!=spaceT.end_leaf(prev_sib) && status==-1)
            {
                
                // display the subspaces================================================
                /*cout << "bound: ";
                 for (int i=0; i<num_obj-1; i++) {
                 cout << (*leafIter).get_bound()[i] << " ";
                 }
                 
                 cout << "sub_status: " << (*leafIter).get_sub_status() << " ";
                 cout << "mod_status: " << (*leafIter).get_mod_status() << " ";
                 
                 if ((*leafIter).get_sub_status()==1)
                 {
                 for (int i=0; i<num_obj; i++) cout << (*leafIter).get_point()[i] << " ";
                 }
                 
                 cout << endl;*/
                
                //=======================================================================
                
                // compare the nodes
                
                // bounds of the current subspace must be greater than or equal to each bound of the compared subspace
                flag=1;
                for (int i=0; i<num_obj-1; i++)
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
                        for (int i=0; i<num_obj-1; i++)
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
        leafIter=spaceT.begin_leaf(sibIter);
        
        while (leafIter!=sub && status==-1)
        {
            
            // display the subspaces======================================
            /*cout << "bound: ";
             for (int i=0; i<num_obj-1; i++) {
             cout << (*leafIter).get_bound()[i] << " ";
             }
             
             cout << "sub_status: " << (*leafIter).get_sub_status() << " ";
             cout << "mod_status: " << (*leafIter).get_mod_status() << " ";
             
             if ((*leafIter).get_sub_status()==1)
             {
             for (int i=0; i<num_obj; i++) cout << (*leafIter).get_point()[i] << " ";
             }
             
             cout << endl;*/
            
            //=============================================================
            
            // compare the nodes
            
            // bounds of the current subspace must be greater than or equal to each bound of the compared subspace
            flag=1;
            for (int i=0; i<num_obj-1; i++)
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
                    for (int i=0; i<num_obj-1; i++)
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
    
    vector<double> coeffs(num_obj);
    double coeff;
    
    int j = 0;
    for(int r=num_mathmodel_row-num_obj;r<num_mathmodel_row;r++){
        coeffs[j]=0.0;
        for(int c=num_obj;c<num_mathmodel_col;c++){
            CPXgetcoef (env, prob, r, c, &coeff);
            coeffs[j]+=pow(coeff,2);
        }
        coeffs[j]=pow(coeffs[j],0.5);
        j++;
    }
    
    return coeffs;
}
