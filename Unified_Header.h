#include <iostream>
#include <stdio.h>
#include <vector>
#include <math.h>
#include <string.h>
#include <dlib/matrix.h>

#define _USE_MATH_DEFINES
#include <cmath>
#include "../../matplotlib-cpp/matplotlibcpp.h"

namespace plt = matplotlibcpp;

using namespace std;

#ifndef USERFUNCTION
#define USERFUNCTION

#include "snoptProblem.hh"
#ifdef __cplusplus
extern "C" {
#endif
    int Default_Init_Pr_( integer    *Status, integer *n,    doublereal x[],
        integer    *needF,  integer *neF,  doublereal F[],
        integer    *needG,  integer *neG,  doublereal G[],
        char       *cu,     integer *lencu,
        integer    iu[],    integer *leniu,
        doublereal ru[],    integer *lenru );

    int Seed_Conf_Optimization_Pr_fn_(integer    *Status, integer *n,    doublereal x[],
        integer    *needF,  integer *neF,  doublereal F[],
        integer    *needG,  integer *neG,  doublereal G[],
        char       *cu,     integer *lencu,
        integer    iu[],    integer *leniu,
        doublereal ru[],    integer *lenru );

    int Nodes_Optimization_Pr_fn(integer    *Status, integer *n,    doublereal x[],
        integer    *needF,  integer *neF,  doublereal F[],
        integer    *needG,  integer *neG,  doublereal G[],
        char       *cu,     integer *lencu,
        integer    iu[],    integer *leniu,
        doublereal ru[],    integer *lenru );

#ifdef __cplusplus
}
#endif

#endif

class Robot_StateNDot {
public:
	double rIx, rIy, theta, q1, q2, q3, q4, q5, q6, q7, q8, q9, q10;
	double rIxdot, rIydot, thetadot, q1dot, q2dot, q3dot, q4dot, q5dot, q6dot, q7dot, q8dot, q9dot, q10dot;
    Robot_StateNDot();
    Robot_StateNDot(std::vector<double> &);
};
typedef struct Tree_Node *Tree_Node_Ptr;    // This is a pointer to the Tree_Node
struct Tree_Node
{   Robot_StateNDot Node_StateNDot;     double KE;    int Node_Index;
    // This means the index of the current nodes in the whole node array
    Tree_Node_Ptr Parent_Node;          std::vector<Tree_Node_Ptr> Children_Nodes;
    std::vector<double> sigma;
    dlib::matrix<double,12,1> End_Effector_Pos, End_Effector_Vel;
};
class Unified_Structure_P{
public:
    Tree_Node Node_i;
    Tree_Node Node_i_child;
    std::vector<double> sigma_tran, sigma_goal;
    int Opt_Val_No, ObjNConst_No;
};
extern Unified_Structure_P Structure_P;
extern std::vector<Tree_Node_Ptr> All_Nodes, Children_Nodes, Frontier_Nodes;

/**
 * Functions that have been successfully tested!
 * Description
 */

 void Envi_Map_Defi();
 void Envi_Map_Normal_Cal(dlib::matrix<double> &Envi_Map);
 void Add_Node2Tree(Tree_Node &Current_Node);
 int Minimum_Index(std::vector<double> &Given_vec);
 Tree_Node Pop_Node();
 void ObjNConstraint_ValNType_Update(dlib::matrix<double> &Matrix_result, std::vector<double> &ObjNConstraint_Val, std::vector<double> &ObjNConstraint_Type, int Constraint_Type);
 std::vector<double> StateNDot2StateVec(const Robot_StateNDot &Robot_StateNDot_i);
 Robot_StateNDot StateVec2StateNDot(std::vector<double> &StateVec);
 std::vector<double> Default_Init(const std::vector<double> &sigma_i);
 std::vector<double> Default_Init_Opt(std::vector<double> &Robot_State_Init);
 void Default_Init_Pr_ObjNConstraint(std::vector<double> &Opt_Seed, std::vector<double> &ObjNConstraint_Val, std::vector<double> &ObjNConstraint_Type);
 dlib::matrix<double> Middle_Joint_Obs_Dist_Fn(Robot_StateNDot &StateNDot_Init_i);
 dlib::matrix<double> ONES_VECTOR_fn(int Dim);
 void Obs_Dist_Fn(std::vector<double> &r_Pos, double &Obs_Dist, int &Obs_Dist_Index, const char* name);

 void Robot_Plot_fn(const Robot_StateNDot &StateNDot_Init_i);
 void Robot_Plot_fn(const Robot_StateNDot &StateNDot_Init_i, std::string &name);
 dlib::matrix<double> D_q_fn(const Robot_StateNDot &Robot_StateNDot_i);
 dlib::matrix<double> B_q_fn();
 dlib::matrix<double> C_q_qdot_fn(const Robot_StateNDot &Robot_StateNDot_i);
 dlib::matrix<double> Jac_Full_fn(const Robot_StateNDot &Robot_StateNDot_i);
 dlib::matrix<double> Jacdot_qdot_fn(const Robot_StateNDot &Robot_StateNDot_i);
 std::vector<double> Ang_Pos_fn(const Robot_StateNDot &Robot_StateNDot_i, const char* s);
 std::vector<double> Ang_Vel_fn(const Robot_StateNDot &Robot_StateNDot_i, const char* s);
 double Kinetic_Energy_fn(Robot_StateNDot &Robot_StateNDot_i);
 std::vector<double> Vec_Minus(std::vector<double> &vec1, std::vector<double> &vec2);
 void End_Effector_Obs_Dist_Fn(dlib::matrix<double,12,1> &End_Effector_Pos, dlib::matrix<double,6,1> &End_Effector_Dist, std::vector<int> &End_Effector_Obs);
 void End_Effector_PosNVel(Robot_StateNDot &StateNDot_Init_i, dlib::matrix<double,12,1> &End_Effector_Pos, dlib::matrix<double,12,1> &End_Effector_Vel);
 void Node_UpdateNCon(Tree_Node &Node_i, Robot_StateNDot &Node_StateNDot_i, std::vector<double> &sigma);
 double Nodes_Optimization_ObjNConstraint(std::vector<double> &Opt_Seed, std::vector<double> &ObjNConstraint_Val, std::vector<double> &ObjNConstraint_Type);
 void Contact_Force_Feasibility_fn(dlib::matrix<double> &Contact_Force_k, dlib::matrix<int> &End_Effector_Obs_k, std::vector<double> &ObjNConstraint_Val, std::vector<double> &ObjNConstraint_Type);
 dlib::matrix<double> Contact_Force_Complem_Matrix_fn(std::vector<double> &sigma);
 dlib::matrix<double> Contact_Status_fn(std::vector<double> &sigma_i);
 dlib::matrix<double> Contact_Acc_Constraint(dlib::matrix<double> &Robotstate, dlib::matrix<double> &Robotstatedot, const dlib::matrix<double> &Sigma_Matrix);
 dlib::matrix<double> Diag_Matrix_fn(std::vector<double> &diag_vec);
 void Dynamics_Matrices(const Robot_StateNDot &Node_StateNDot, dlib::matrix<double> &D_q, dlib::matrix<double> &B_q, dlib::matrix<double> &C_q_qdot, dlib::matrix<double> &Jac_Full);
 dlib::matrix<double> State_Ctrl_CF_2_Statedot(dlib::matrix<double> &Robotstate_k, dlib::matrix<double> &Contact_Force_k, dlib::matrix<double> &Control_k);
 double Traj_Variation(dlib::matrix<double> &StateNDot_Traj);

 double Kinetic_Energy_End_Frame(dlib::matrix<double> &StateNDot_Traj);
 double KE_Variation_fn(dlib::matrix<double> &StateNDot_Traj);
 void Contact_Force_Proj(dlib::matrix<double> &Contact_Force_k, std::vector<double> &Normal_Force_vec, std::vector<double> &Tange_Force_vec, dlib::matrix<int> & End_Effector_Obs_k);
 int Nodes_Optimization_fn(Tree_Node &Node_i, Tree_Node &Node_i_child, std::vector<double> &Opt_Soln_Output);
 std::vector<double> Nodes_Optimization_Inner_Opt(Tree_Node &Node_i, Tree_Node &Node_i_child);
 double ObjNConstraint_Violation(const std::vector<double> &ObjNConstraint_Val, const std::vector<double> &ObjNConstraint_Type);
 void Opt_Soln_Write2Txt(Tree_Node &Node_i,Tree_Node &Node_i_child, std::vector<double> &Opt_Soln);
 dlib::matrix<double> Node_Expansion_fn(const Tree_Node &Node_i, int &Adjacent_Number);
 std::vector<double> End_RobotNDot_Extract(std::vector<double> &Opt_Soln, std::vector<double> &sigma, std::vector<double>&sigma_i_child);
 Robot_StateNDot Impact_Mapping_fn(dlib::matrix<double> &End_State_Dlib, double &Impulse_Mag, std::vector<double> &sigma_i_child);
 std::vector<double> Full_Row_Rank_Index(std::vector<double> &sigma_i_child);
 std::vector<double> Opt_Soln_Load();

 std::vector<double> Seed_Guess_Gene(Tree_Node &Node_i, Tree_Node &Node_i_child);
 void Opt_Seed_Unzip(std::vector<double> &Opt_Seed, double &T_tot, dlib::matrix<double> & StateNDot_Traj, dlib::matrix<double> & Ctrl_Traj, dlib::matrix<double> & Contact_Force_Traj, dlib::matrix<double> & Contact_Force_Mid_Traj);
 void Opt_Seed_Zip(std::vector<double> &Opt_Seed, dlib::matrix<double> & StateNDot_Traj, dlib::matrix<double> & Ctrl_Traj, dlib::matrix<double> & Contact_Force_Traj, dlib::matrix<double> & Contact_Force_Mid_Traj);
 std::vector<double> Seed_Guess_Gene_Robotstate(Tree_Node &Node_i, Tree_Node &Node_i_child);
 void Seed_Conf_Optimization_ObjNConstraint(std::vector<double> &Opt_Seed, std::vector<double> &ObjNConstraint_Val, std::vector<double> &ObjNConstraint_Type);
 std::vector<double> Time_Seed_Queue_fn(double Time_Interval, int Total_Num);

 std::vector<double> Sigma2Pos(std::vector<double> &sigma, int EqOrIneq);
 std::vector<double> Sigma2Vel(std::vector<double> &sigma);
 int Sigma_TransNGoal(std::vector<double> & sigma_i, std::vector<double> & sigma_i_child);
 Robot_StateNDot DlibRobotstate2StateNDot(dlib::matrix<double> &DlibRobotstate);
 dlib::matrix<double> Eqn_Maint_Matrix_fn(std::vector<double> &sigma_i, std::vector<double> &sigma_i_child);
 std::vector<double> CubicSpline_Coeff_fn(double T, double x_init, double x_end, double xdot_init, double xdot_end);
 dlib::matrix<double> Dynamics_RHS_Matrix_fn(dlib::matrix<double> &Jac_Full, dlib::matrix<double> &B_q);
 void StateNAcc2ContactForceNTorque(dlib::matrix<double> &State_k, dlib::matrix<double> &Acc_k, dlib::matrix<double> &Lamda_k, dlib::matrix<double> &u_k);
 dlib::matrix<double> StateNAccNTorque2ContactForce(dlib::matrix<double> &State_k, dlib::matrix<double> &Acc_k, dlib::matrix<double> &u_k);
 void State_Traj_Interpolater(std::vector<double>&Init_Robotstate, std::vector<double>&End_Robotstate, dlib::matrix<double> &StateNDot_Traj, dlib::matrix<double> &Acc_Traj, dlib::matrix<double> &StateNdot_Mid_Traj, dlib::matrix<double> &Acc_Mid_Traj, int Choice_Flag);

 void CubicSpline_Interpolater(dlib::matrix<double> &StateNDot_Traj, dlib::matrix<double> &Acc_Traj, dlib::matrix<double> &StateNdot_Mid_Traj, dlib::matrix<double> &Acc_Mid_Traj);
 double CubicSpline_Eval(double T, std::vector<double> & CubicSpline_Coeff_vec, double s, char name);
 std::vector<double> QuadraticSpline_Coeff_fn(double T, double x_init, double x_end, double xdot_init);
 double QuadraticSpline_Eval(double T, std::vector<double> & State_Traj_Coeff, double s, char name);
 int Sigma_Change(std::vector<double> &sigma_i, std::vector<double> &sigma_i_child);
