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
void Robot_Plot_fn(const Robot_StateNDot &StateNDot_Init_i);
void Robot_Plot_fn(const Robot_StateNDot &StateNDot_Init_i, std::string &name);
dlib::matrix<double> D_q_fn(const Robot_StateNDot &Robot_StateNDot_i);
dlib::matrix<double>  B_q_fn();
dlib::matrix<double>  C_q_qdot_fn(const Robot_StateNDot &Robot_StateNDot_i);
dlib::matrix<double>  Jac_Full_fn(const Robot_StateNDot &Robot_StateNDot_i);
dlib::matrix<double> Jacdot_qdot_fn(const Robot_StateNDot &Robot_StateNDot_i);
double Kinetic_Energy_fn(Robot_StateNDot &Robot_StateNDot_i);
std::vector<double> StateNDot2StateVec(const Robot_StateNDot &Robot_StateNDot_i);
Robot_StateNDot StateVec2StateNDot(std::vector<double> &StateVec);

void Default_Init_Pr_ObjNConstraint(std::vector<double> &Opt_Seed, std::vector<double> &ObjNConstraint_Val, std::vector<double> &ObjNConstraint_Type);
void ObjNConstraint_ValNType_Update(dlib::matrix<double> &Matrix_result, std::vector<double> &ObjNConstraint_Val, std::vector<double> &ObjNConstraint_Type, int Constraint_Type);
void End_Effector_Obs_Dist_Fn(dlib::matrix<double,12,1> &End_Effector_Pos, dlib::matrix<double,6,1> &End_Effector_Dist, std::vector<int> &End_Effector_Obs);

std::vector<double> Ang_Pos_fn(const Robot_StateNDot &Robot_StateNDot_i, const char* s);
std::vector<double> Ang_Vel_fn(const Robot_StateNDot &Robot_StateNDot_i, const char* s);

std::vector<double> Default_Init(const std::vector<double> &sigma_i);
void End_Effector_PosNVel(Robot_StateNDot &StateNDot_Init_i, dlib::matrix<double,12,1> &End_Effector_Pos, dlib::matrix<double,12,1> &End_Effector_Vel);
dlib::matrix<double> Middle_Joint_Obs_Dist_Fn(Robot_StateNDot &StateNDot_Init_i);
std::vector<double> Vec_Minus(std::vector<double> &vec1, std::vector<double> &vec2);
dlib::matrix<double> ONES_VECTOR_fn(int Dim);

void Obs_Dist_Fn(std::vector<double> &r_Pos, double &Obs_Dist, int &Obs_Dist_Index, const char* name);
dlib::matrix<double> Diag_Matrix_fn(std::vector<double> &diag_vec);

std::vector<double> Sigma2Pos(std::vector<double> &sigma, int EqOrIneq);
std::vector<double> Sigma2Vel(std::vector<double> &sigma);
dlib::matrix<double> Diag_Matrix_fn(std::vector<double> &diag_vec);
void Seed_Conf_Optimization_ObjNConstraint(std::vector<double> &Opt_Seed, std::vector<double> &ObjNConstraint_Val, std::vector<double> &ObjNConstraint_Type);
std::vector<double> Seed_Guess_Gene_Robotstate(Tree_Node &Node_i, Tree_Node &Node_i_child);
std::vector<double> Seed_Guess_Gene(Tree_Node &Node_i, Tree_Node &Node_i_child);
dlib::matrix<double> Eqn_Maint_Matrix_fn(std::vector<double> &sigma_i, std::vector<double> &sigma_i_child);
Tree_Node Pop_Node();
void Node_UpdateNCon(Tree_Node &Node_i, Robot_StateNDot &Node_StateNDot_i, std::vector<double> &sigma);

int Nodes_Optimization_fn(Tree_Node &Node_i, Tree_Node &Node_i_child, std::vector<double> &Opt_Soln_Output);

void Dynamics_Matrices(const Robot_StateNDot &Node_StateNDot, dlib::matrix<double> &D_q, dlib::matrix<double> &B_q, dlib::matrix<double> &C_q_qdot, dlib::matrix<double> &Jac_Full);
std::vector<double> CubicSpline_Coeff_fn(double T, double x_init, double x_end, double xdot_init, double xdot_end);
std::vector<double> CubicSpline_PosVelAcc4(double T, double a, double b, double c, double d, double s);
void Pos_Vel_Acc_VelfromPos_fromStateNdot_Coeff(double T, dlib::matrix<double> &StateNDot_Coeff, int Grid_Ind, double s, std::vector<double> &Robot_Config,  std::vector<double> &Robot_Vel, dlib::matrix<double> &Robot_Acc, std::vector<double> &Robot_VelfromPos);
dlib::matrix<double> Dynamics_RHS_Matrix_fn(dlib::matrix<double> &Jac_Full, dlib::matrix<double> &B_q);

std::vector<double> PosNVel2StateVec(std::vector<double> & Pos, std::vector<double> & Vel);
std::vector<double> CubicSpline_PosVelAcc8(double T, double x_a, double x_b, double x_c, double x_d, double xdot_a, double xdot_b, double xdot_c, double xdot_d, double s);
void Ctrl_Contact_Force_Coeff_fn(dlib::matrix<double> &Ctrl_Traj, dlib::matrix<double> &Contact_Force_Traj, dlib::matrix<double> &Ctrl_Coeff, dlib::matrix<double> &Contact_Force_Coeff);

void Opt_Seed_Zip(std::vector<double> &Opt_Seed, dlib::matrix<double> & StateNDot_Coeff, dlib::matrix<double> & Ctrl_Coeff, dlib::matrix<double> & Contact_Force_Coeff);
void Opt_Seed_Unzip(std::vector<double> &Opt_Seed, double &T_tot, dlib::matrix<double> & StateNDot_Traj, dlib::matrix<double> & Ctrl_Traj, dlib::matrix<double> & Contact_Force_Traj);

void Sigma_TransNGoal(std::vector<double> & sigma_i, std::vector<double> & sigma_i_child, int &Self_Opt_Flag, int &Crit_Grid);
void CtrlNContact_ForcefromCtrlNContact_Force_Coeff(dlib::matrix<double> &Ctrl_Coeff,dlib::matrix<double> &Contact_Force_Coeff, int Grid_Ind, double s, dlib::matrix<double> &Ctrl_i,  dlib::matrix<double> &Contact_Force_i);
dlib::matrix<double> StateVec2DlibMatrix_fn(const std::vector<double> &StateVec);

void Real_ObjNConstraint_Stage(Tree_Node &Node_i, Tree_Node &Node_i_child, std::vector<double>& sigma, double T, std::vector<double> &Pos, std::vector<double> &Vel, dlib::matrix<double> &Acc, dlib::matrix<double> &Ctrl, dlib::matrix<double> &Contact_Force, std::vector<double> &ObjNConstraint_Val, std::vector<double> &ObjNConstraint_Type);
void Dynamics_Constraint(std::vector<double> &Pos, std::vector<double> &Vel, dlib::matrix<double> &Acc, dlib::matrix<double> &Ctrl, dlib::matrix<double> &Contact_Force, std::vector<double> &ObjNConstraint_Val, std::vector<double> &ObjNConstraint_Type);
void Distance_Velocity_Constraint(std::vector<double>& sigma, std::vector<double> &Pos, std::vector<double> &Vel, std::vector<double> &ObjNConstraint_Val, std::vector<double> &ObjNConstraint_Type);
void Contact_Force_Complem_Constraint(dlib::matrix<double> &Contact_Force, std::vector<double> &sigma, std::vector<double> &ObjNConstraint_Val, std::vector<double> &ObjNConstraint_Type);
void Contact_Force_Feasibility_Constraint(std::vector<double> &Pos, std::vector<double> &Vel, dlib::matrix<double> &Contact_Force, std::vector<double> &ObjNConstraint_Val, std::vector<double> &ObjNConstraint_Type);
void Contact_Maintenance_Constraint(Tree_Node &Node_i, Tree_Node &Node_i_child, std::vector<double> &Pos, std::vector<double> &Vel, std::vector<double> &ObjNConstraint_Val, std::vector<double> &ObjNConstraint_Type);
double KE_Variation_fn(dlib::matrix<double> &StateNDot_Traj);

dlib::matrix<double> StateNDot_ref_fn(std::vector<double> &Robot_Config_i, std::vector<double> &Robot_Velocity_i);
std::vector<double> Opt_Soln_Load();
dlib::matrix<double> Quadratic_Minus(dlib::matrix<double> &Mat_A, dlib::matrix<double> &Mat_B);
void Quadratic_Angular_Sum_Cal(std::vector<double> &Robot_Vel,double &Quadratic_Angular_Sum);

Robot_StateNDot DlibRobotstate2StateNDot(dlib::matrix<double> &DlibRobotstate);
double CubicSpline_Evaluation_fn(const std::vector<double> &CubicSpline_Coeff, double s);
double CubicSpline_1stOrder_Evaluation_fn(const std::vector<double> &CubicSpline_Coeff, double s, double T);

void Robot_StateNDot_MidNAcc(double T, const Robot_StateNDot &Robot_StateNDot_Front, const Robot_StateNDot &Robot_StateNDot_Back, const dlib::matrix<double> &Ctrl_Front, const dlib::matrix<double> &Ctrl_Back, const dlib::matrix<double> &Contact_Force_Front, const dlib::matrix<double> &Contact_Force_Back, Robot_StateNDot &Robot_StateNDot_Mid, dlib::matrix<double> &Robotstate_Mid_Acc,std::vector<double> &ObjNConstraint_Val, std::vector<double> &ObjNConstraint_Type);
double Traj_Variation(dlib::matrix<double> &StateNDot_Traj);
double ObjNConstraint_Violation(const std::vector<double> &ObjNConstraint_Val, const std::vector<double> &ObjNConstraint_Type);
std::vector<double> Nodes_Optimization_Inner_Opt(Tree_Node &Node_i, Tree_Node &Node_i_child);

std::vector<double> Time_Seed_Queue_fn(double Time_Interval, int Total_Num);

dlib::matrix<double> Node_Expansion_fn(const Tree_Node &Node_i, int &Adjacent_Number);
std::vector<double> End_RobotNDot_Extract(std::vector<double> &Opt_Soln, std::vector<double> &sigma, std::vector<double>&sigma_i_child);

void Opt_Soln_Write2Txt(Tree_Node &Node_i,Tree_Node &Node_i_child, std::vector<double> &Opt_Soln);

std::vector<double> Final_State_Opt(std::vector<double> &Opt_Soln, Tree_Node &Node_i, Tree_Node &Node_i_child);

std::vector<double> Default_Init_Opt(std::vector<double> &Robot_State_Init);

void End_Effector_Upper_Vel(Robot_StateNDot &StateNDot_Init_i, dlib::matrix<double> &End_Effector_Upper_Normal_Speed);

void Contact_Force_Proj(dlib::matrix<double> &StateNDot_Traj, dlib::matrix<double> &Contact_Force_Traj, std::vector<double> &Normal_Force_vec, std::vector<double> &Tange_Force_vec, int Grid_Index, std::vector<double> &Full_Normal);

double Objective_Function_Cal(dlib::matrix<double> &StateNDot_Traj, int Opt_Type_Flag, std::vector<double> &sigma_i_child);

Robot_StateNDot Impact_Mapping_fn(dlib::matrix<double> &StateNDot_Traj, double &Impulse_Mag, std::vector<double> &sigma_i_child);

std::vector<double> Full_Row_Rank_Index(std::vector<double> &sigma_i_child);
void Integration_Consistent(dlib::matrix<double> &Robostate_Dlib_Front, dlib::matrix<double> &Robostate_Dlib_Back, std::vector<double> &ObjNConstraint_Val, std::vector<double> &ObjNConstraint_Type);

double Kinetic_Energy_End_Frame(dlib::matrix<double> &StateNDot_Traj);

double Velocity_Projection(std::vector<double> &Pos_A, std::vector<double> &Pos_B, std::vector<double> &Vel_B);

double Variation_Cal(dlib::matrix<double> &Dlib_Matrix);

double Quadratic_Sum(dlib::matrix<double> &Dlib_Matrix);

int Sigma_Change(std::vector<double> &sigma_i, std::vector<double> &sigma_i_child);
