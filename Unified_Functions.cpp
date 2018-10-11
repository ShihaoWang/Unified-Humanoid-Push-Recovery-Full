#include <stdio.h>
#include <ctime>
#include <cstdlib>
#include <string.h>
#include <iostream>
#include "snopt.hh"
#include "Unified_Header.h"

#include <fstream>
#include <cmath>
#include <dlib/matrix.h>
#include "snoptProblem.hh"
#include <algorithm>

double Inf = 1.1e20;										double PI = 3.1415926535897932384626;									   							double pi = PI;

// There are three types of variables in this optimization problem: robot state, control torques and contact forces
double rIxlow = -Inf;                  	double rIxupp = Inf;					double rIylow = -Inf;                  	double rIyupp = Inf;
double thetalow = -pi;                 	double thetaupp = pi;					double q1low = -2.3562;                	double q1upp = 0.733;
double q2low = 0.0;                    	double q2upp = 2.618;					double q3low = -1.3;                   	double q3upp = 0.733;
double q4low = -2.3562;                	double q4upp = 0.733;					double q5low = 0.0;                    	double q5upp = 2.618;
double q6low = -1.3;                   	double q6upp = 0.733;					double q7low = -3.14;                  	double q7upp = 1.047;
double q8low = -2.391;                 	double q8upp = 0.0;						double q9low = -3.14;                  	double q9upp = 1.047;
double q10low = -2.391;                	double q10upp = 0.0;
double AngRateMag = 3.0;			   				double AngRateLow = -AngRateMag;     	double AngRateHgh = AngRateMag;

double rIxdotlow = -Inf;               	double rIxdotupp = Inf;					double rIydotlow = -Inf;               	double rIydotupp = Inf;
double thetadotlow = -Inf;             	double thetadotupp = Inf;				double q1dotlow = AngRateLow;          	double q1dotupp = AngRateHgh;
double q2dotlow = AngRateLow;          	double q2dotupp = AngRateHgh;			double q3dotlow = AngRateLow;          	double q3dotupp = AngRateHgh;
double q4dotlow = AngRateLow;          	double q4dotupp = AngRateHgh;			double q5dotlow = AngRateLow;          	double q5dotupp = AngRateHgh;
double q6dotlow = AngRateLow;          	double q6dotupp = AngRateHgh;			double q7dotlow = AngRateLow;          	double q7dotupp = AngRateHgh;
double q8dotlow = AngRateLow;          	double q8dotupp = AngRateHgh;			double q9dotlow = AngRateLow;          	double q9dotupp = AngRateHgh;
double q10dotlow = AngRateLow;         	double q10dotupp = AngRateHgh;

double tau1_max = 100;             			double tau2_max = 100;					double tau3_max = 100;					double tau4_max = 100;
double tau5_max = 100;             			double tau6_max = 100;					double tau7_max = 60;              		double tau8_max = 50;
double tau9_max = 60;             			double tau10_max = 50;


// actually in the new method to formulate the optimization problem, the acceleration is left unbounded.
double Acc_max = 10;

dlib::matrix<double> xlow_vec;							dlib::matrix<double> xupp_vec;
dlib::matrix<double> ctrl_low_vec;					dlib::matrix<double> ctrl_upp_vec;
dlib::matrix<double>  Envi_Map;							dlib::matrix<double> Envi_Map_Normal, Envi_Map_Tange; // The Normal and tangential vector of the plane

/**
* Some global values are defined
* Description
*/

int Grids = 8;			double mu = 0.35;

int Variable_Num = (26 + 10 + 12) * Grids + 1;   //

double h_k; 															      // This value will be adaptively changed to formulate an optimal solution the meaning of this variable is the time step between two grids
std::vector<Tree_Node_Ptr> All_Nodes;						// All nodes are here!
std::vector<Tree_Node_Ptr> Children_Nodes;			// All children nodes!
std::vector<Tree_Node_Ptr> Frontier_Nodes;			// Only Frontier ndoes!
std::vector<double> Frontier_Nodes_Cost;		    // The kinetic energy of each nodes

// Used for the manifold projection
dlib::matrix<double> Jac_Act_Trans_Manifold;
std::vector<double> Pote_Robotstate(26), Sigma_Act_Manifold(6);

void Envi_Map_Defi()
{
  // 	This function is used to define the environment map for the simulation
  // Whenever this function gets called, it will return the array with the
  // environment obstacle information
  //
  // This is the default flat ground
  // This map is defined in a polyline manner with the first value denoting
  // the line length and the second value denoting the relative angle
  Envi_Map = dlib::ones_matrix<double>(2,4);
  Envi_Map(0,0) = -5.0;					Envi_Map(0,1) = 0.0;				Envi_Map(0,2) = 4.975;				Envi_Map(0,3) = 0.0;
  Envi_Map(1,0) = 4.975;				Envi_Map(1,1) = 0.0;				Envi_Map(1,2) = 4.975;				Envi_Map(1,3) = 3.0;
}

void Envi_Map_Normal_Cal(dlib::matrix<double> &Envi_Map)
{
  // This function is used to calculate the surface normal vector and tangential vector
  int NumOfObs = Envi_Map.nr();
  int Dimension = Envi_Map.nc()/2;
  Envi_Map_Normal = dlib::zeros_matrix<double>(NumOfObs,Dimension);
  Envi_Map_Tange = dlib::zeros_matrix<double>(NumOfObs,Dimension);
  double Envi_Map_Edge_A_x, Envi_Map_Edge_A_y, Envi_Map_Edge_B_x, Envi_Map_Edge_B_y;
  double Slope_Angle;
  for (int i = 0; i < NumOfObs; i++)
  {
    Envi_Map_Edge_A_x = Envi_Map(i, 0);						Envi_Map_Edge_A_y = Envi_Map(i, 1);
    Envi_Map_Edge_B_x = Envi_Map(i, 2);						Envi_Map_Edge_B_y = Envi_Map(i, 3);
    Slope_Angle = atan2(Envi_Map_Edge_B_y - Envi_Map_Edge_A_y, Envi_Map_Edge_B_x - Envi_Map_Edge_A_x);
    Envi_Map_Tange(i,0) = cos(Slope_Angle);					Envi_Map_Tange(i,1) = sin(Slope_Angle);
    Envi_Map_Normal(i,0) = cos(Slope_Angle + pi/2.0);		Envi_Map_Normal(i,1) = sin(Slope_Angle + pi/2.0);
  }
  return;
}

void Add_Node2Tree(Tree_Node &Current_Node)
{
  // This function will add the current node to the All_Nodes vector
  All_Nodes.push_back(&Current_Node);
  Frontier_Nodes.push_back(&Current_Node);
  // This function can only run after the node has been updated with Node_UpdateNCon()
  Frontier_Nodes_Cost.push_back(Current_Node.KE);
}

int Minimum_Index(std::vector<double> &Given_vec)
{ // Given a double vector, this function will return the index of the minimum value
  int index = 0;
  for(int i = 1; i < Given_vec.size(); i++)
  {
    if(Given_vec[i] < Given_vec[index])
    index = i;
  }
  return index;
}

Tree_Node Pop_Node()
{
  // This function will pop the node outfrom the current Frontier according to the kinetic energy
  // Another thing is that this function will operate on the global variables
  int Min_Ind = Minimum_Index(Frontier_Nodes_Cost);
  Tree_Node Current_Node = *Frontier_Nodes[Min_Ind];
  Frontier_Nodes.erase(Frontier_Nodes.begin()+Min_Ind);
  Frontier_Nodes_Cost.erase(Frontier_Nodes_Cost.begin()+Min_Ind);
  return Current_Node;
}

void ObjNConstraint_ValNType_Update(dlib::matrix<double> &Matrix_result, std::vector<double> &ObjNConstraint_Val, std::vector<double> &ObjNConstraint_Type, int Constraint_Type)
{
  // This function is used to update the value of ObjNConstraint_Val and Type
  for (int i = 0; i < Matrix_result.nr(); i++)
  {
    ObjNConstraint_Val.push_back(Matrix_result(i));
    ObjNConstraint_Type.push_back(Constraint_Type);
  }
}

std::vector<double> StateNDot2StateVec(const Robot_StateNDot &Robot_StateNDot_i)
{
  std::vector<double> StateVec(26);
  StateVec[0] = Robot_StateNDot_i.rIx;              StateVec[1] = Robot_StateNDot_i.rIy;
  StateVec[2] = Robot_StateNDot_i.theta;            StateVec[3] = Robot_StateNDot_i.q1;
  StateVec[4] = Robot_StateNDot_i.q2;               StateVec[5] = Robot_StateNDot_i.q3;
  StateVec[6] = Robot_StateNDot_i.q4;               StateVec[7] = Robot_StateNDot_i.q5;
  StateVec[8] = Robot_StateNDot_i.q6;               StateVec[9] = Robot_StateNDot_i.q7;
  StateVec[10] = Robot_StateNDot_i.q8;              StateVec[11] = Robot_StateNDot_i.q9;
  StateVec[12] = Robot_StateNDot_i.q10;

  StateVec[0+13] = Robot_StateNDot_i.rIxdot;        StateVec[1+13] = Robot_StateNDot_i.rIydot;
  StateVec[2+13] = Robot_StateNDot_i.thetadot;      StateVec[3+13] = Robot_StateNDot_i.q1dot;
  StateVec[4+13] = Robot_StateNDot_i.q2dot;         StateVec[5+13] = Robot_StateNDot_i.q3dot;
  StateVec[6+13] = Robot_StateNDot_i.q4dot;         StateVec[7+13] = Robot_StateNDot_i.q5dot;
  StateVec[8+13] = Robot_StateNDot_i.q6dot;         StateVec[9+13] = Robot_StateNDot_i.q7dot;
  StateVec[10+13] = Robot_StateNDot_i.q8dot;        StateVec[11+13] = Robot_StateNDot_i.q9dot;
  StateVec[12+13] = Robot_StateNDot_i.q10dot;
  return StateVec;
}

Robot_StateNDot StateVec2StateNDot(std::vector<double> &StateVec)
{
  Robot_StateNDot Robot_StateNDot_i;
  Robot_StateNDot_i.rIx = StateVec[0];              Robot_StateNDot_i.rIy = StateVec[1];
  Robot_StateNDot_i.theta = StateVec[2];            Robot_StateNDot_i.q1 = StateVec[3];
  Robot_StateNDot_i.q2 = StateVec[4];               Robot_StateNDot_i.q3 = StateVec[5];
  Robot_StateNDot_i.q4 = StateVec[6];               Robot_StateNDot_i.q5 = StateVec[7];
  Robot_StateNDot_i.q6 = StateVec[8];               Robot_StateNDot_i.q7 = StateVec[9];
  Robot_StateNDot_i.q8 = StateVec[10];              Robot_StateNDot_i.q9 = StateVec[11];
  Robot_StateNDot_i.q10 = StateVec[12];

  Robot_StateNDot_i.rIxdot = StateVec[13];          Robot_StateNDot_i.rIydot = StateVec[14];
  Robot_StateNDot_i.thetadot = StateVec[15];        Robot_StateNDot_i.q1dot = StateVec[16];
  Robot_StateNDot_i.q2dot = StateVec[17];           Robot_StateNDot_i.q3dot = StateVec[18];
  Robot_StateNDot_i.q4dot = StateVec[19];           Robot_StateNDot_i.q5dot = StateVec[20];
  Robot_StateNDot_i.q6dot = StateVec[21];           Robot_StateNDot_i.q7dot = StateVec[22];
  Robot_StateNDot_i.q8dot = StateVec[23];           Robot_StateNDot_i.q9dot = StateVec[24];
  Robot_StateNDot_i.q10dot = StateVec[25];
  return Robot_StateNDot_i;
}

std::vector<double> Default_Init(const std::vector<double> &sigma_i)
{
  // This function is used to initialize the whole optimization process
  // First, is to substitute the map info into the Envi_Map matrix
  // Second, is to give the proper bounds to the variables to be optimized
  // Thrid, is to generate a kinematically feasible initial robot state

  Envi_Map_Defi();
  // Then it is to compute the normal and tangential vector of the map
  Envi_Map_Normal_Cal(Envi_Map);
  xlow_vec = dlib::zeros_matrix<double>(26,1);					xupp_vec = dlib::zeros_matrix<double>(26,1);
  ctrl_low_vec = dlib::zeros_matrix<double>(10,1);			ctrl_upp_vec = dlib::matrix<double>(10,1) ;
  xlow_vec(0) = rIxlow; 						xupp_vec(0) = rIxupp;
  xlow_vec(1) = rIylow; 						xupp_vec(1) = rIyupp;
  xlow_vec(2) = thetalow; 					xupp_vec(2) = thetaupp;
  xlow_vec(3) = q1low; 							xupp_vec(3) = q1upp;
  xlow_vec(4) = q2low; 							xupp_vec(4) = q2upp;
  xlow_vec(5) = q3low; 							xupp_vec(5) = q3upp;
  xlow_vec(6) = q4low; 							xupp_vec(6) = q4upp;
  xlow_vec(7) = q5low; 							xupp_vec(7) = q5upp;
  xlow_vec(8) = q6low; 							xupp_vec(8) = q6upp;
  xlow_vec(9) = q7low; 							xupp_vec(9) = q7upp;
  xlow_vec(10) = q8low; 						xupp_vec(10) = q8upp;
  xlow_vec(11) = q9low; 						xupp_vec(11) = q9upp;
  xlow_vec(12) = q10low; 						xupp_vec(12) = q10upp;
  xlow_vec(0+13) = rIxdotlow; 			xupp_vec(0+13) = rIxdotupp;
  xlow_vec(1+13) = rIydotlow; 			xupp_vec(1+13) = rIydotupp;
  xlow_vec(2+13) = thetadotlow; 		xupp_vec(2+13) = thetadotupp;
  xlow_vec(3+13) = q1dotlow; 				xupp_vec(3+13) = q1dotupp;
  xlow_vec(4+13) = q2dotlow; 				xupp_vec(4+13) = q2dotupp;
  xlow_vec(5+13) = q3dotlow; 				xupp_vec(5+13) = q3dotupp;
  xlow_vec(6+13) = q4dotlow; 				xupp_vec(6+13) = q4dotupp;
  xlow_vec(7+13) = q5dotlow; 				xupp_vec(7+13) = q5dotupp;
  xlow_vec(8+13) = q6dotlow; 				xupp_vec(8+13) = q6dotupp;
  xlow_vec(9+13) = q7dotlow; 				xupp_vec(9+13) = q7dotupp;
  xlow_vec(10+13) = q8dotlow; 			xupp_vec(10+13) = q8dotupp;
  xlow_vec(11+13) = q9dotlow; 			xupp_vec(11+13) = q9dotupp;
  xlow_vec(12+13) = q10dotlow; 			xupp_vec(12+13) = q10dotupp;

  ctrl_low_vec(0) = -tau1_max;			ctrl_upp_vec(0) = -ctrl_low_vec(0);
  ctrl_low_vec(1) = -tau2_max;			ctrl_upp_vec(1) = -ctrl_low_vec(1);
  ctrl_low_vec(2) = -tau3_max;			ctrl_upp_vec(2) = -ctrl_low_vec(2);
  ctrl_low_vec(3) = -tau4_max;			ctrl_upp_vec(3) = -ctrl_low_vec(3);
  ctrl_low_vec(4) = -tau5_max;			ctrl_upp_vec(4) = -ctrl_low_vec(4);
  ctrl_low_vec(5) = -tau6_max;			ctrl_upp_vec(5) = -ctrl_low_vec(5);
  ctrl_low_vec(6) = -tau7_max;			ctrl_upp_vec(6) = -ctrl_low_vec(6);
  ctrl_low_vec(7) = -tau8_max;			ctrl_upp_vec(7) = -ctrl_low_vec(7);
  ctrl_low_vec(8) = -tau9_max;			ctrl_upp_vec(8) = -ctrl_low_vec(8);
  ctrl_low_vec(9) = -tau10_max;			ctrl_upp_vec(9) = -ctrl_low_vec(9);

  vector<double> Robot_State_Init;
  ifstream Initial_Robot_State_File;              // This is to read the initial angle and angular velocities
  Initial_Robot_State_File.open("robot_angle_init.txt");
  if(Initial_Robot_State_File.is_open())
  {
    double data_each_line = 0.0;
    while(Initial_Robot_State_File>>data_each_line)
    {
      Robot_State_Init.push_back(data_each_line);
    }
    Initial_Robot_State_File.close();
  }
  else
  {
    printf("Unable to open robot_angle_init.txt file!\n");
  }
  Initial_Robot_State_File.open("robot_velocity_init.txt");
  if(Initial_Robot_State_File.is_open())
  {
    double data_each_line = 0.0;
    while(Initial_Robot_State_File>>data_each_line)
    {
      Robot_State_Init.push_back(data_each_line);
    }
    Initial_Robot_State_File.close();
  }
  else
  {
    printf("Unable to open robot_velocity_init.txt file!\n");
  }

  // Here the robot initial state has been read into the Robot_State_Init vector
  Tree_Node RootNode;
  RootNode.Node_StateNDot = StateVec2StateNDot(Robot_State_Init);
  RootNode.sigma = sigma_i;
  Structure_P.Node_i = RootNode;
  // Robot_Plot_fn(RootNode.Node_StateNDot);
  // // If the default configuration would like to be viewed
  // Robot_StateNDot Robot_StateNDot_init(Robot_State_Init);
  // std::string input_name = "init_given";
  // Robot_Plot_fn(Robot_StateNDot_init,input_name);

  Robot_State_Init = Default_Init_Opt(Robot_State_Init);
  return Robot_State_Init;
}

std::vector<double> Default_Init_Opt(std::vector<double> &Robot_State_Init)
{
  snoptProblem Default_Init_Pr;                     // This is the name of the Optimization problem
  // Allocate and initialize
  std:vector<double> ObjNConstraint_Val, ObjNConstraint_Type;
  integer n = Robot_State_Init.size();
  Default_Init_Pr_ObjNConstraint(Robot_State_Init, ObjNConstraint_Val, ObjNConstraint_Type);
  integer neF = ObjNConstraint_Val.size();     							  // 1 objective function
  integer lenA  =  n * neF;                         // This is the number of nonzero elements in the linear part A    F(x) = f(x)+Ax

  integer *iAfun = new integer[lenA];              integer *jAvar = new integer[lenA];				doublereal *A  = new doublereal[lenA];

  integer lenG   = lenA;							 integer *iGfun = new integer[lenG];				integer *jGvar = new integer[lenG];

  doublereal *x      = new doublereal[n];			doublereal *xlow   = new doublereal[n];				doublereal *xupp   = new doublereal[n];
  doublereal *xmul   = new doublereal[n];			integer    *xstate = new    integer[n];

  doublereal *F      = new doublereal[neF];		doublereal *Flow   = new doublereal[neF];			doublereal *Fupp   = new doublereal[neF];
  doublereal *Fmul   = new doublereal[neF];		integer    *Fstate = new integer[neF];

  integer nxnames = 1;							integer nFnames = 1;						char *xnames = new char[nxnames*8];						char *Fnames = new char[nFnames*8];

  integer    ObjRow = 0;							doublereal ObjAdd = 0;

  // Set the upper and lower bounds.
  for (int i = 0; i < n; i++) {
    xlow[i] = xlow_vec(i);
    xupp[i] = xupp_vec(i);
    xstate[i] = 0.0;
    x[i] = Robot_State_Init[i];  	// Initial guess
  }

  for(int i = 0; i<neF; i++)
  {
    // The lower bound is the same
    Flow[i] = 0.0;
    if(ObjNConstraint_Type[i]>0)	// Inequality constraint
    {
      Fupp[i] = Inf;
    }
    else
    {
      Fupp[i] = 0.0;
    }
  }
  // Load the data for ToyProb ...
  Default_Init_Pr.setPrintFile  ( "Default_Init_Pr.out" );
  Default_Init_Pr.setProblemSize( n, neF );
  Default_Init_Pr.setObjective  ( ObjRow, ObjAdd );
  Default_Init_Pr.setA          ( lenA, iAfun, jAvar, A );
  Default_Init_Pr.setG          ( lenG, iGfun, jGvar );
  Default_Init_Pr.setX          ( x, xlow, xupp, xmul, xstate );
  Default_Init_Pr.setF          ( F, Flow, Fupp, Fmul, Fstate );
  Default_Init_Pr.setXNames     ( xnames, nxnames );
  Default_Init_Pr.setFNames     ( Fnames, nFnames );
  Default_Init_Pr.setProbName   ( "Default_Init_Pr" );
  Default_Init_Pr.setUserFun    ( Default_Init_Pr_);
  // snopta will compute the Jacobian by finite-differences.
  // The user has the option of calling  snJac  to define the
  // coordinate arrays (iAfun,jAvar,A) and (iGfun, jGvar).
  Default_Init_Pr.computeJac    ();
  Default_Init_Pr.setIntParameter( "Derivative option", 0 );
  Default_Init_Pr.setIntParameter( "Major print level", 0 );
  Default_Init_Pr.setIntParameter( "Minor print level", 0 );
  integer Cold = 0, Basis = 1, Warm = 2;
  Default_Init_Pr.solve( Cold );

  // Take the value out from x
  for (int i = 0; i < n; i++)
  {
    Robot_State_Init[i] = x[i];
  }
  // Robot_StateNDot Init_Opt_vec(Robot_State_Init);
  // Robot_Plot_fn(Init_Opt_vec);

  delete []iAfun;  delete []jAvar;  delete []A;
  delete []iGfun;  delete []jGvar;

  delete []x;      delete []xlow;   delete []xupp;
  delete []xmul;   delete []xstate;

  delete []F;		 delete []Flow;	  delete []Fupp;
  delete []Fmul;	 delete []Fstate;

  delete []xnames; delete []Fnames;

  return Robot_State_Init;
}

int Default_Init_Pr_(integer    *Status, integer *n,    doublereal x[],
  integer    *needF,  integer *neF,  doublereal F[],
  integer    *needG,  integer *neG,  doublereal G[],
  char       *cu,     integer *lencu,
  integer    iu[],    integer *leniu,
  doublereal ru[],    integer *lenru )
  {
    std::vector<double> ObjNConstraint_Val,ObjNConstraint_Type;
    // Initial guess of the robot configurations
    std::vector<double> Robot_State_Init = StateNDot2StateVec(Structure_P.Node_i.Node_StateNDot);
    std::vector<double> Robot_State_Opt;
    for (int i = 0; i < 26; i++)
    {
      Robot_State_Opt.push_back(x[i]);
    }
    Default_Init_Pr_ObjNConstraint(Robot_State_Opt, ObjNConstraint_Val,ObjNConstraint_Type);
    for (int i = 0; i < ObjNConstraint_Val.size(); i++)
    {
      F[i] = ObjNConstraint_Val[i];
    }
    return 0;
  }

void Default_Init_Pr_ObjNConstraint(std::vector<double> &Opt_Seed, std::vector<double> &ObjNConstraint_Val, std::vector<double> &ObjNConstraint_Type)
{
  Robot_StateNDot StateNDot_Init_i(Opt_Seed);		dlib::matrix<double,12,1> End_Effector_Pos, End_Effector_Vel;
  End_Effector_PosNVel(StateNDot_Init_i, End_Effector_Pos, End_Effector_Vel);
  std::vector<double> Robostate_ref = StateNDot2StateVec(Structure_P.Node_i.Node_StateNDot);
  std::vector<double> Robostate_offset = Vec_Minus(Opt_Seed, Robostate_ref);
  double Robostate_offset_val = 0.0;
  for (int i = 0; i < Robostate_offset.size()/2; i++)
  {
    Robostate_offset_val = Robostate_offset_val + Robostate_offset[i] * Robostate_offset[i];
  }
  std::vector<double> sigma = Structure_P.Node_i.sigma;
  ObjNConstraint_Val.push_back(Robostate_offset_val);
  ObjNConstraint_Type.push_back(1);

  dlib::matrix<double,6,1> End_Effector_Dist;
  std::vector<int> End_Effector_Obs(6);

  End_Effector_Obs_Dist_Fn(End_Effector_Pos, End_Effector_Dist, End_Effector_Obs);

  dlib::matrix<double> Eqn_Pos_Matrix, Ineqn_Pos_Matrix, Eqn_Vel_Matrix, Matrix_result;
  std::vector<double> sigma_temp;
  sigma_temp = Sigma2Pos(sigma, 0);
  Eqn_Pos_Matrix = Diag_Matrix_fn(sigma_temp);
  sigma_temp = Sigma2Pos(sigma, 1);
  Ineqn_Pos_Matrix = Diag_Matrix_fn(sigma_temp);
  sigma_temp = Sigma2Vel(sigma);
  Eqn_Vel_Matrix = Diag_Matrix_fn(sigma_temp);
  // cout<<Eqn_Pos_Matrix<<endl;				cout<<Ineqn_Pos_Matrix<<endl;			cout<<Eqn_Vel_Matrix<<endl;

  // 1. Active constraints have to be satisfied: Position and Velocity
  Matrix_result = Eqn_Pos_Matrix * End_Effector_Dist;
  ObjNConstraint_ValNType_Update(Matrix_result, ObjNConstraint_Val, ObjNConstraint_Type, 0);

  // cout<<Matrix_result<<endl;

  Matrix_result = Eqn_Vel_Matrix * End_Effector_Vel;
  ObjNConstraint_ValNType_Update(Matrix_result, ObjNConstraint_Val, ObjNConstraint_Type, 0);

  // cout<<Matrix_result<<endl;

  double mini = 0.025;

  // 2. Inactive constraints have to be strictly away from the obstacle
  dlib::matrix<double> ones_vector, temp_matrix;
  ones_vector = ONES_VECTOR_fn(6);
  Matrix_result = Ineqn_Pos_Matrix * (End_Effector_Dist - ones_vector * mini);
  ObjNConstraint_ValNType_Update(Matrix_result, ObjNConstraint_Val, ObjNConstraint_Type, 1);

  double rAx = End_Effector_Pos(0);			  double rAy = End_Effector_Pos(1);
  double rBx = End_Effector_Pos(2);			  double rBy = End_Effector_Pos(3);
  double rCx = End_Effector_Pos(4);			  double rCy = End_Effector_Pos(5);
  double rDx = End_Effector_Pos(6);			  double rDy = End_Effector_Pos(7);
  double rEx = End_Effector_Pos(8);			  double rEy = End_Effector_Pos(9);
  double rFx = End_Effector_Pos(10);			double rFy = End_Effector_Pos(11);

  // 3. Middle joints have to be strictly away from the obs
  temp_matrix = Middle_Joint_Obs_Dist_Fn(StateNDot_Init_i);
  Matrix_result = temp_matrix - ones_vector * mini;
  ObjNConstraint_ValNType_Update(Matrix_result, ObjNConstraint_Val, ObjNConstraint_Type, 1);

  // 4. There may also be a constraint on the kinetic energy of the robot
  double KE_Init = Kinetic_Energy_fn(StateNDot_Init_i);
  ObjNConstraint_Val.push_back(KE_Init - 40);
  ObjNConstraint_Type.push_back(0);

  return;
}

dlib::matrix<double> Middle_Joint_Obs_Dist_Fn(Robot_StateNDot &StateNDot_Init_i)
{
  dlib::matrix<double> Middle_Joint_Obs_Dist_Matrix;
  Middle_Joint_Obs_Dist_Matrix = dlib::ones_matrix<double>(6,1);
  std::vector<double> rH = Ang_Pos_fn(StateNDot_Init_i, "rH");
  std::vector<double> rK = Ang_Pos_fn(StateNDot_Init_i, "rK");
  std::vector<double> rM = Ang_Pos_fn(StateNDot_Init_i, "rM");
  std::vector<double> rN = Ang_Pos_fn(StateNDot_Init_i, "rN");
  std::vector<double> rI = Ang_Pos_fn(StateNDot_Init_i, "rI");
  std::vector<double> rT = Ang_Pos_fn(StateNDot_Init_i, "rT");
  int Obs_Dist_Index;		double rH_Obs, rK_Obs, rM_Obs, rN_Obs, rI_Obs, rT_Obs;

  Obs_Dist_Fn(rH, rH_Obs, Obs_Dist_Index, "None");
  Obs_Dist_Fn(rK, rK_Obs, Obs_Dist_Index, "None");
  Obs_Dist_Fn(rM, rM_Obs, Obs_Dist_Index, "None");
  Obs_Dist_Fn(rN, rN_Obs, Obs_Dist_Index, "None");
  Obs_Dist_Fn(rI, rI_Obs, Obs_Dist_Index, "None");
  Obs_Dist_Fn(rT, rT_Obs, Obs_Dist_Index, "None");

  Middle_Joint_Obs_Dist_Matrix(0) = rH_Obs;
  Middle_Joint_Obs_Dist_Matrix(1) = rK_Obs;
  Middle_Joint_Obs_Dist_Matrix(2) = rM_Obs;
  Middle_Joint_Obs_Dist_Matrix(3) = rN_Obs;
  Middle_Joint_Obs_Dist_Matrix(4) = rI_Obs;
  Middle_Joint_Obs_Dist_Matrix(5) = rT_Obs;
  return Middle_Joint_Obs_Dist_Matrix;
}

dlib::matrix<double> ONES_VECTOR_fn(int Dim)
{
  dlib::matrix<double> Ones_vector;
  Ones_vector = dlib::ones_matrix<double>(Dim,1);
  return Ones_vector;
}

void Obs_Dist_Fn(std::vector<double> &r_Pos, double &Obs_Dist, int &Obs_Dist_Index, const char* name)
{
  // 	This function is used to calculate the relative distance between the robot end effector and the nearby environment
  double Temp_Edge_x, Temp_Edge_y, Temp_offset_x, Temp_offset_y, Normal_vector_i_x, Normal_vector_i_y;
  std::vector<double> Obs_Dist_vec;
  for (int i = 0; i < Envi_Map.nr(); i++)
  {
    Temp_Edge_x = Envi_Map(i,0);
    Temp_Edge_y = Envi_Map(i,1);
    Temp_offset_x = r_Pos[0] - Temp_Edge_x;
    Temp_offset_y = r_Pos[1] - Temp_Edge_y;
    Normal_vector_i_x = Envi_Map_Normal(i,0);
    Normal_vector_i_y = Envi_Map_Normal(i,1);
    Obs_Dist_vec.push_back(Temp_offset_x * Normal_vector_i_x + Temp_offset_y * Normal_vector_i_y);
    if(strcmp(name,"floor")==0)
    {
      // In this case, we are talking about the robot foot, it can only make contact with the floor
      // So the iteration will terminate after the first run, the code assumption is that the ground floor will be put into the first row of the map
      break;
    }
    if(strcmp(name,"wall")==0)
    {
      Obs_Dist_vec[0] = 10;      // This is used to make sure that only the second row will be used
    }
  }
  Obs_Dist_Index = Minimum_Index(Obs_Dist_vec);
  Obs_Dist = Obs_Dist_vec[Obs_Dist_Index];
  return;
}

Robot_StateNDot::Robot_StateNDot()
{
  // A default constructor
  rIx = 0;			rIy = 0.7230;			theta = -0.0900;
  q1 = 0.3768;		q2 = 0.0045;			q3 = -0.2913;			q4 = -1.0015;			q5 = 0.1500;
  q6 = 0.2698;		q7 = -0.6600;			q8 = -0.6251;			q9 = 0.6900;			q10 = -0.2951;
  rIxdot = 0.2000;	rIydot = -0.0605;		thetadot = -0.2100;
  q1dot = -0.1239;	q2dot = 1.3108;			q3dot = -0.9768;		q4dot = -1.4999;		q5dot = 2.0000;
  q6dot = -1.2999;	q7dot = 1.0000;			q8dot = -2.0000;		q9dot = -1.5708;		q10dot = -1.5000;
}

Robot_StateNDot::Robot_StateNDot(std::vector<double> &Robot_AngleNRate)
{
  // An evaluated constructor
  rIx = Robot_AngleNRate[0];		rIy = Robot_AngleNRate[1];		theta = Robot_AngleNRate[2];
  q1 = Robot_AngleNRate[3];		q2 = Robot_AngleNRate[4];		q3 = Robot_AngleNRate[5];		q4 = Robot_AngleNRate[6];		q5 = Robot_AngleNRate[7];
  q6 = Robot_AngleNRate[8];		q7 = Robot_AngleNRate[9];		q8 = Robot_AngleNRate[10];		q9 = Robot_AngleNRate[11];		q10 = Robot_AngleNRate[12];

  rIxdot = Robot_AngleNRate[13];	rIydot = Robot_AngleNRate[14];	thetadot = Robot_AngleNRate[15];
  q1dot = Robot_AngleNRate[16];	q2dot = Robot_AngleNRate[17];	q3dot = Robot_AngleNRate[18];	q4dot = Robot_AngleNRate[19];	q5dot = Robot_AngleNRate[20];
  q6dot = Robot_AngleNRate[21];	q7dot = Robot_AngleNRate[22];	q8dot = Robot_AngleNRate[23];	q9dot = Robot_AngleNRate[24];	q10dot = Robot_AngleNRate[25];
}

void Robot_Plot_fn(const Robot_StateNDot &StateNDot_Init_i)
{
  // This function is used to plot the robot configuration
  std::vector<double> rA = Ang_Pos_fn(StateNDot_Init_i, "rA");
  std::vector<double> rB = Ang_Pos_fn(StateNDot_Init_i, "rB");
  std::vector<double> rC = Ang_Pos_fn(StateNDot_Init_i, "rC");
  std::vector<double> rD = Ang_Pos_fn(StateNDot_Init_i, "rD");
  std::vector<double> rE = Ang_Pos_fn(StateNDot_Init_i, "rE");
  std::vector<double> rF = Ang_Pos_fn(StateNDot_Init_i, "rF");
  std::vector<double> rG = Ang_Pos_fn(StateNDot_Init_i, "rG");
  std::vector<double> rH = Ang_Pos_fn(StateNDot_Init_i, "rH");
  std::vector<double> rI = Ang_Pos_fn(StateNDot_Init_i, "rI");
  std::vector<double> rJ = Ang_Pos_fn(StateNDot_Init_i, "rJ");
  std::vector<double> rK = Ang_Pos_fn(StateNDot_Init_i, "rK");
  std::vector<double> rL = Ang_Pos_fn(StateNDot_Init_i, "rL");
  std::vector<double> rM = Ang_Pos_fn(StateNDot_Init_i, "rM");
  std::vector<double> rN = Ang_Pos_fn(StateNDot_Init_i, "rN");
  std::vector<double> rT = Ang_Pos_fn(StateNDot_Init_i, "rT");

  // One thing that I would like to do is to change the coor of the plot in the matplotlib

  std::vector<double> x(2), y(2);
  // AB
  x.at(0) = rA[0];	x.at(1) = rB[0];	y.at(0) = rA[1];	y.at(1) = rB[1];	plt::plot(x,y);
  // CD
  x.at(0) = rC[0];	x.at(1) = rD[0];	y.at(0) = rC[1];	y.at(1) = rD[1];	plt::plot(x,y);
  // CJ
  x.at(0) = rC[0];	x.at(1) = rJ[0];	y.at(0) = rC[1];	y.at(1) = rJ[1];	plt::plot(x,y);
  // JD
  x.at(0) = rJ[0];	x.at(1) = rD[0];	y.at(0) = rJ[1];	y.at(1) = rD[1];	plt::plot(x,y);
  // KJ
  x.at(0) = rK[0];	x.at(1) = rJ[0];	y.at(0) = rK[1];	y.at(1) = rJ[1];	plt::plot(x,y);
  // IK
  x.at(0) = rI[0];	x.at(1) = rK[0];	y.at(0) = rI[1];	y.at(1) = rK[1];	plt::plot(x,y);
  // BG
  x.at(0) = rB[0];	x.at(1) = rG[0];	y.at(0) = rB[1];	y.at(1) = rG[1];	plt::plot(x,y);
  // AG
  x.at(0) = rA[0];	x.at(1) = rG[0];	y.at(0) = rA[1];	y.at(1) = rG[1];	plt::plot(x,y);
  // GH
  x.at(0) = rG[0];	x.at(1) = rH[0];	y.at(0) = rG[1];	y.at(1) = rH[1];	plt::plot(x,y);
  // HI
  x.at(0) = rH[0];	x.at(1) = rI[0];	y.at(0) = rH[1];	y.at(1) = rI[1];	plt::plot(x,y);
  // IT
  x.at(0) = rI[0];	x.at(1) = rT[0];	y.at(0) = rI[1];	y.at(1) = rT[1];	plt::plot(x,y);
  // LM
  x.at(0) = rL[0];	x.at(1) = rM[0];	y.at(0) = rL[1];	y.at(1) = rM[1];	plt::plot(x,y);
  // ME
  x.at(0) = rM[0];	x.at(1) = rE[0];	y.at(0) = rM[1];	y.at(1) = rE[1];	plt::plot(x,y);
  // LN
  x.at(0) = rL[0];	x.at(1) = rN[0];	y.at(0) = rL[1];	y.at(1) = rN[1];	plt::plot(x,y);
  // NF
  x.at(0) = rN[0];	x.at(1) = rF[0];	y.at(0) = rN[1];	y.at(1) = rF[1];	plt::plot(x,y);
  const char* filename = "./init.png";
  std::cout<<"Saving result to "<<filename<<std::endl;
  plt::save(filename);
}

void Robot_Plot_fn(const Robot_StateNDot &StateNDot_Init_i, std::string &name)
{
  // This function is used to plot the robot configuration
  std::vector<double> rA = Ang_Pos_fn(StateNDot_Init_i, "rA");
  std::vector<double> rB = Ang_Pos_fn(StateNDot_Init_i, "rB");
  std::vector<double> rC = Ang_Pos_fn(StateNDot_Init_i, "rC");
  std::vector<double> rD = Ang_Pos_fn(StateNDot_Init_i, "rD");
  std::vector<double> rE = Ang_Pos_fn(StateNDot_Init_i, "rE");
  std::vector<double> rF = Ang_Pos_fn(StateNDot_Init_i, "rF");
  std::vector<double> rG = Ang_Pos_fn(StateNDot_Init_i, "rG");
  std::vector<double> rH = Ang_Pos_fn(StateNDot_Init_i, "rH");
  std::vector<double> rI = Ang_Pos_fn(StateNDot_Init_i, "rI");
  std::vector<double> rJ = Ang_Pos_fn(StateNDot_Init_i, "rJ");
  std::vector<double> rK = Ang_Pos_fn(StateNDot_Init_i, "rK");
  std::vector<double> rL = Ang_Pos_fn(StateNDot_Init_i, "rL");
  std::vector<double> rM = Ang_Pos_fn(StateNDot_Init_i, "rM");
  std::vector<double> rN = Ang_Pos_fn(StateNDot_Init_i, "rN");
  std::vector<double> rT = Ang_Pos_fn(StateNDot_Init_i, "rT");

  std::vector<double> x(2), y(2);
  // AB
  x.at(0) = rA[0];	x.at(1) = rB[0];	y.at(0) = rA[1];	y.at(1) = rB[1];	plt::plot(x,y);
  // CD
  x.at(0) = rC[0];	x.at(1) = rD[0];	y.at(0) = rC[1];	y.at(1) = rD[1];	plt::plot(x,y);
  // CJ
  x.at(0) = rC[0];	x.at(1) = rJ[0];	y.at(0) = rC[1];	y.at(1) = rJ[1];	plt::plot(x,y);
  // JD
  x.at(0) = rJ[0];	x.at(1) = rD[0];	y.at(0) = rJ[1];	y.at(1) = rD[1];	plt::plot(x,y);
  // KJ
  x.at(0) = rK[0];	x.at(1) = rJ[0];	y.at(0) = rK[1];	y.at(1) = rJ[1];	plt::plot(x,y);
  // IK
  x.at(0) = rI[0];	x.at(1) = rK[0];	y.at(0) = rI[1];	y.at(1) = rK[1];	plt::plot(x,y);
  // BG
  x.at(0) = rB[0];	x.at(1) = rG[0];	y.at(0) = rB[1];	y.at(1) = rG[1];	plt::plot(x,y);
  // AG
  x.at(0) = rA[0];	x.at(1) = rG[0];	y.at(0) = rA[1];	y.at(1) = rG[1];	plt::plot(x,y);
  // GH
  x.at(0) = rG[0];	x.at(1) = rH[0];	y.at(0) = rG[1];	y.at(1) = rH[1];	plt::plot(x,y);
  // HI
  x.at(0) = rH[0];	x.at(1) = rI[0];	y.at(0) = rH[1];	y.at(1) = rI[1];	plt::plot(x,y);
  // IT
  x.at(0) = rI[0];	x.at(1) = rT[0];	y.at(0) = rI[1];	y.at(1) = rT[1];	plt::plot(x,y);
  // LM
  x.at(0) = rL[0];	x.at(1) = rM[0];	y.at(0) = rL[1];	y.at(1) = rM[1];	plt::plot(x,y);
  // ME
  x.at(0) = rM[0];	x.at(1) = rE[0];	y.at(0) = rM[1];	y.at(1) = rE[1];	plt::plot(x,y);
  // LN
  x.at(0) = rL[0];	x.at(1) = rN[0];	y.at(0) = rL[1];	y.at(1) = rN[1];	plt::plot(x,y);
  // NF
  x.at(0) = rN[0];	x.at(1) = rF[0];	y.at(0) = rN[1];	y.at(1) = rF[1];	plt::plot(x,y);

  std::string pre_filename = "./";
  std::string post_filename = ".png";
  std::string filename = pre_filename + name + post_filename;

  std::cout<<"Saving result to "<<filename<<std::endl;
  plt::save(filename);
}

dlib::matrix<double> D_q_fn(const Robot_StateNDot &Robot_StateNDot_i)
{
  double rIx = Robot_StateNDot_i.rIx;						double rIy = Robot_StateNDot_i.rIy;
  double theta = Robot_StateNDot_i.theta;				double q1 = Robot_StateNDot_i.q1;
  double q2 = Robot_StateNDot_i.q2;							double q3 = Robot_StateNDot_i.q3;
  double q4 = Robot_StateNDot_i.q4;							double q5 = Robot_StateNDot_i.q5;
  double q6 = Robot_StateNDot_i.q6;							double q7 = Robot_StateNDot_i.q7;
  double q8 = Robot_StateNDot_i.q8;							double q9 = Robot_StateNDot_i.q9;
  double q10 = Robot_StateNDot_i.q10;

  dlib::matrix<double>  T;									T = dlib::zeros_matrix<double>(13,13);
  T(0,0) = 1.13E2/2.0;
  T(0,2) = cos(q1+q2+theta)*(-9.1E1/1.0E2)-cos(q4+q5+theta)*(9.1E1/1.0E2)-cos(q7+q8+theta)*(9.9E1/1.6E2)-cos(q9+q10+theta)*(9.9E1/1.6E2)-cos(q1+q2+q3+theta)*(2.0/2.5E1)-cos(q4+q5+q6+theta)*(2.0/2.5E1)+sin(q1+q2+q3+theta)*(2.0/2.5E1)+sin(q4+q5+q6+theta)*(2.0/2.5E1)-cos(q1+theta)*(3.9E1/2.0E1)-cos(q4+theta)*(3.9E1/2.0E1)-cos(q7+theta)*(1.9E1/1.6E1)-cos(q9+theta)*(1.9E1/1.6E1)+cos(theta)*(2.97E2/2.0E1)-sqrt(4.1E1)*cos(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/5.0E1)-sqrt(4.1E1)*cos(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/5.0E1);
  T(0,3) = cos(q1+q2+theta)*(-9.1E1/1.0E2)-cos(q1+q2+q3+theta)*(2.0/2.5E1)+sin(q1+q2+q3+theta)*(2.0/2.5E1)-cos(q1+theta)*(3.9E1/2.0E1)-sqrt(4.1E1)*cos(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/5.0E1);
  T(0,4) = cos(q1+q2+theta)*(-9.1E1/1.0E2)-cos(q1+q2+q3+theta)*(2.0/2.5E1)+sin(q1+q2+q3+theta)*(2.0/2.5E1)-sqrt(4.1E1)*cos(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/5.0E1);
  T(0,5) = cos(q1+q2+q3+theta)*(-2.0/2.5E1)+sin(q1+q2+q3+theta)*(2.0/2.5E1)-sqrt(4.1E1)*cos(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/5.0E1);
  T(0,6) = cos(q4+q5+theta)*(-9.1E1/1.0E2)-cos(q4+q5+q6+theta)*(2.0/2.5E1)+sin(q4+q5+q6+theta)*(2.0/2.5E1)-cos(q4+theta)*(3.9E1/2.0E1)-sqrt(4.1E1)*cos(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/5.0E1);
  T(0,7) = cos(q4+q5+theta)*(-9.1E1/1.0E2)-cos(q4+q5+q6+theta)*(2.0/2.5E1)+sin(q4+q5+q6+theta)*(2.0/2.5E1)-sqrt(4.1E1)*cos(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/5.0E1);
  T(0,8) = cos(q4+q5+q6+theta)*(-2.0/2.5E1)+sin(q4+q5+q6+theta)*(2.0/2.5E1)-sqrt(4.1E1)*cos(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/5.0E1);
  T(0,9) = cos(q7+q8+theta)*(-9.9E1/1.6E2)-cos(q7+theta)*(1.9E1/1.6E1);
  T(0,10) = cos(q7+q8+theta)*(-9.9E1/1.6E2);
  T(0,11) = cos(q9+q10+theta)*(-9.9E1/1.6E2)-cos(q9+theta)*(1.9E1/1.6E1);
  T(0,12) = cos(q9+q10+theta)*(-9.9E1/1.6E2);
  T(1,1) = 1.13E2/2.0;
  T(1,2) = sin(q1+q2+theta)*(9.1E1/1.0E2)+sin(q4+q5+theta)*(9.1E1/1.0E2)+sin(q7+q8+theta)*(9.9E1/1.6E2)+sin(q9+q10+theta)*(9.9E1/1.6E2)+cos(q1+q2+q3+theta)*(2.0/2.5E1)+cos(q4+q5+q6+theta)*(2.0/2.5E1)+sin(q1+q2+q3+theta)*(2.0/2.5E1)+sin(q4+q5+q6+theta)*(2.0/2.5E1)+sin(q1+theta)*(3.9E1/2.0E1)+sin(q4+theta)*(3.9E1/2.0E1)+sin(q7+theta)*(1.9E1/1.6E1)+sin(q9+theta)*(1.9E1/1.6E1)-sin(theta)*(2.97E2/2.0E1)+sqrt(4.1E1)*sin(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/5.0E1)+sqrt(4.1E1)*sin(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/5.0E1);
  T(1,3) = sin(q1+q2+theta)*(9.1E1/1.0E2)+cos(q1+q2+q3+theta)*(2.0/2.5E1)+sin(q1+q2+q3+theta)*(2.0/2.5E1)+sin(q1+theta)*(3.9E1/2.0E1)+sqrt(4.1E1)*sin(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/5.0E1);
  T(1,4) = sin(q1+q2+theta)*(9.1E1/1.0E2)+cos(q1+q2+q3+theta)*(2.0/2.5E1)+sin(q1+q2+q3+theta)*(2.0/2.5E1)+sqrt(4.1E1)*sin(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/5.0E1);
  T(1,5) = cos(q1+q2+q3+theta)*(2.0/2.5E1)+sin(q1+q2+q3+theta)*(2.0/2.5E1)+sqrt(4.1E1)*sin(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/5.0E1);
  T(1,6) = sin(q4+q5+theta)*(9.1E1/1.0E2)+cos(q4+q5+q6+theta)*(2.0/2.5E1)+sin(q4+q5+q6+theta)*(2.0/2.5E1)+sin(q4+theta)*(3.9E1/2.0E1)+sqrt(4.1E1)*sin(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/5.0E1);
  T(1,7) = sin(q4+q5+theta)*(9.1E1/1.0E2)+cos(q4+q5+q6+theta)*(2.0/2.5E1)+sin(q4+q5+q6+theta)*(2.0/2.5E1)+sqrt(4.1E1)*sin(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/5.0E1);
  T(1,8) = cos(q4+q5+q6+theta)*(2.0/2.5E1)+sin(q4+q5+q6+theta)*(2.0/2.5E1)+sqrt(4.1E1)*sin(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/5.0E1);
  T(1,9) = sin(q7+q8+theta)*(9.9E1/1.6E2)+sin(q7+theta)*(1.9E1/1.6E1);
  T(1,10) = sin(q7+q8+theta)*(9.9E1/1.6E2);
  T(1,11) = sin(q9+q10+theta)*(9.9E1/1.6E2)+sin(q9+theta)*(1.9E1/1.6E1);
  T(1,12) = sin(q9+q10+theta)*(9.9E1/1.6E2);
  T(2,0) = cos(q1+q2+theta)*(-9.1E1/1.0E2)-cos(q4+q5+theta)*(9.1E1/1.0E2)-cos(q7+q8+theta)*(9.9E1/1.6E2)-cos(q9+q10+theta)*(9.9E1/1.6E2)-cos(q1+q2+q3+theta)*(2.0/2.5E1)-cos(q4+q5+q6+theta)*(2.0/2.5E1)+sin(q1+q2+q3+theta)*(2.0/2.5E1)+sin(q4+q5+q6+theta)*(2.0/2.5E1)-cos(q1+theta)*(3.9E1/2.0E1)-cos(q4+theta)*(3.9E1/2.0E1)-cos(q7+theta)*(1.9E1/1.6E1)-cos(q9+theta)*(1.9E1/1.6E1)+cos(theta)*(2.97E2/2.0E1)-sqrt(4.1E1)*cos(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/5.0E1)-sqrt(4.1E1)*cos(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/5.0E1);
  T(2,1) = sin(q1+q2+theta)*(9.1E1/1.0E2)+sin(q4+q5+theta)*(9.1E1/1.0E2)+sin(q7+q8+theta)*(9.9E1/1.6E2)+sin(q9+q10+theta)*(9.9E1/1.6E2)+cos(q1+q2+q3+theta)*(2.0/2.5E1)+cos(q4+q5+q6+theta)*(2.0/2.5E1)+sin(q1+q2+q3+theta)*(2.0/2.5E1)+sin(q4+q5+q6+theta)*(2.0/2.5E1)+sin(q1+theta)*(3.9E1/2.0E1)+sin(q4+theta)*(3.9E1/2.0E1)+sin(q7+theta)*(1.9E1/1.6E1)+sin(q9+theta)*(1.9E1/1.6E1)-sin(theta)*(2.97E2/2.0E1)+sqrt(4.1E1)*sin(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/5.0E1)+sqrt(4.1E1)*sin(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/5.0E1);
  T(2,2) = cos(q2+q3)*(1.3E1/2.5E2)+cos(q5+q6)*(1.3E1/2.5E2)-cos(q7+q8)*6.80625E-1-cos(q9+q10)*6.80625E-1-sin(q2+q3)*(1.3E1/2.5E2)-sin(q5+q6)*(1.3E1/2.5E2)+cos(q2)*5.915E-1+cos(q3)*(1.3E1/2.5E2)+cos(q5)*5.915E-1+cos(q6)*(1.3E1/2.5E2)-cos(q7)*(2.09E2/1.6E2)+cos(q8)*(9.9E1/3.2E2)-cos(q9)*(2.09E2/1.6E2)+cos(q10)*(9.9E1/3.2E2)-sin(q3)*(1.3E1/2.5E2)-sin(q6)*(1.3E1/2.5E2)+sqrt(4.1E1)*cos(q3-8.960553845713439E-1)*(1.3E1/1.0E3)+sqrt(4.1E1)*cos(q6-8.960553845713439E-1)*(1.3E1/1.0E3)+sqrt(4.1E1)*cos(8.960553845713439E-1)*(1.0/2.5E2)-sqrt(4.1E1)*sin(8.960553845713439E-1)*(1.0/2.5E2)+sqrt(4.1E1)*cos(q2+q3-8.960553845713439E-1)*(1.3E1/1.0E3)+sqrt(4.1E1)*cos(q5+q6-8.960553845713439E-1)*(1.3E1/1.0E3)+1.171739583333333E1;
  T(2,3) = cos(q2+q3)*(1.3E1/2.5E2)-sin(q2+q3)*(1.3E1/2.5E2)+cos(q2)*5.915E-1+cos(q3)*(1.3E1/2.5E2)-sin(q3)*(1.3E1/2.5E2)+sqrt(4.1E1)*cos(q3-8.960553845713439E-1)*(1.3E1/1.0E3)+sqrt(4.1E1)*cos(8.960553845713439E-1)*(1.0/5.0E2)-sqrt(4.1E1)*sin(8.960553845713439E-1)*(1.0/5.0E2)+sqrt(4.1E1)*cos(q2+q3-8.960553845713439E-1)*(1.3E1/1.0E3)+1.409583333333333;
  T(2,4) = cos(q2+q3)*(1.3E1/5.0E2)-sin(q2+q3)*(1.3E1/5.0E2)+cos(q2)*2.9575E-1+cos(q3)*(1.3E1/2.5E2)-sin(q3)*(1.3E1/2.5E2)+sqrt(4.1E1)*cos(q3-8.960553845713439E-1)*(1.3E1/1.0E3)+sqrt(4.1E1)*cos(8.960553845713439E-1)*(1.0/5.0E2)-sqrt(4.1E1)*sin(8.960553845713439E-1)*(1.0/5.0E2)+sqrt(4.1E1)*cos(q2+q3-8.960553845713439E-1)*6.5E-3+5.29375E-1;
  T(2,5) = cos(q2+q3)*(1.3E1/5.0E2)-sin(q2+q3)*(1.3E1/5.0E2)+cos(q3)*(1.3E1/5.0E2)-sin(q3)*(1.3E1/5.0E2)+sqrt(4.1E1)*cos(q3-8.960553845713439E-1)*6.5E-3+sqrt(4.1E1)*cos(8.960553845713439E-1)*(1.0/5.0E2)-sqrt(4.1E1)*sin(8.960553845713439E-1)*(1.0/5.0E2)+sqrt(4.1E1)*cos(q2+q3-8.960553845713439E-1)*6.5E-3+8.575E-2;
  T(2,6) = cos(q5+q6)*(1.3E1/2.5E2)-sin(q5+q6)*(1.3E1/2.5E2)+cos(q5)*5.915E-1+cos(q6)*(1.3E1/2.5E2)-sin(q6)*(1.3E1/2.5E2)+sqrt(4.1E1)*cos(q6-8.960553845713439E-1)*(1.3E1/1.0E3)+sqrt(4.1E1)*cos(8.960553845713439E-1)*(1.0/5.0E2)-sqrt(4.1E1)*sin(8.960553845713439E-1)*(1.0/5.0E2)+sqrt(4.1E1)*cos(q5+q6-8.960553845713439E-1)*(1.3E1/1.0E3)+1.409583333333333;
  T(2,7) = cos(q5+q6)*(1.3E1/5.0E2)-sin(q5+q6)*(1.3E1/5.0E2)+cos(q5)*2.9575E-1+cos(q6)*(1.3E1/2.5E2)-sin(q6)*(1.3E1/2.5E2)+sqrt(4.1E1)*cos(q6-8.960553845713439E-1)*(1.3E1/1.0E3)+sqrt(4.1E1)*cos(8.960553845713439E-1)*(1.0/5.0E2)-sqrt(4.1E1)*sin(8.960553845713439E-1)*(1.0/5.0E2)+sqrt(4.1E1)*cos(q5+q6-8.960553845713439E-1)*6.5E-3+5.29375E-1;
  T(2,8) = cos(q5+q6)*(1.3E1/5.0E2)-sin(q5+q6)*(1.3E1/5.0E2)+cos(q6)*(1.3E1/5.0E2)-sin(q6)*(1.3E1/5.0E2)+sqrt(4.1E1)*cos(q6-8.960553845713439E-1)*6.5E-3+sqrt(4.1E1)*cos(8.960553845713439E-1)*(1.0/5.0E2)-sqrt(4.1E1)*sin(8.960553845713439E-1)*(1.0/5.0E2)+sqrt(4.1E1)*cos(q5+q6-8.960553845713439E-1)*6.5E-3+8.575E-2;
  T(2,9) = cos(q7+q8)*(-3.403125E-1)-cos(q7)*(2.09E2/3.2E2)+cos(q8)*(9.9E1/3.2E2)+1.045989583333333;
  T(2,10) = cos(q7+q8)*(-3.403125E-1)+cos(q8)*(9.9E1/6.4E2)+6.0328125E-1;
  T(2,11) = cos(q9+q10)*(-3.403125E-1)-cos(q9)*(2.09E2/3.2E2)+cos(q10)*(9.9E1/3.2E2)+1.045989583333333;
  T(2,12) = cos(q9+q10)*(-3.403125E-1)+cos(q10)*(9.9E1/6.4E2)+6.0328125E-1;
  T(3,0) = cos(q1+q2+theta)*(-9.1E1/1.0E2)-cos(q1+q2+q3+theta)*(2.0/2.5E1)+sin(q1+q2+q3+theta)*(2.0/2.5E1)-cos(q1+theta)*(3.9E1/2.0E1)-sqrt(4.1E1)*cos(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/5.0E1);
  T(3,1) = sin(q1+q2+theta)*(9.1E1/1.0E2)+cos(q1+q2+q3+theta)*(2.0/2.5E1)+sin(q1+q2+q3+theta)*(2.0/2.5E1)+sin(q1+theta)*(3.9E1/2.0E1)+sqrt(4.1E1)*sin(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/5.0E1);
  T(3,2) = cos(q2+q3)*(1.3E1/2.5E2)-sin(q2+q3)*(1.3E1/2.5E2)+cos(q2)*5.915E-1+cos(q3)*(1.3E1/2.5E2)-sin(q3)*(1.3E1/2.5E2)+sqrt(4.1E1)*cos(q3-8.960553845713439E-1)*(1.3E1/1.0E3)+sqrt(4.1E1)*cos(8.960553845713439E-1)*(1.0/5.0E2)-sqrt(4.1E1)*sin(8.960553845713439E-1)*(1.0/5.0E2)+sqrt(4.1E1)*cos(q2+q3-8.960553845713439E-1)*(1.3E1/1.0E3)+1.409583333333333;
  T(3,3) = cos(q2+q3)*(1.3E1/2.5E2)-sin(q2+q3)*(1.3E1/2.5E2)+cos(q2)*5.915E-1+cos(q3)*(1.3E1/2.5E2)-sin(q3)*(1.3E1/2.5E2)+sqrt(4.1E1)*cos(q3-8.960553845713439E-1)*(1.3E1/1.0E3)+sqrt(4.1E1)*cos(8.960553845713439E-1)*(1.0/5.0E2)-sqrt(4.1E1)*sin(8.960553845713439E-1)*(1.0/5.0E2)+sqrt(4.1E1)*cos(q2+q3-8.960553845713439E-1)*(1.3E1/1.0E3)+1.409583333333333;
  T(3,4) = cos(q2+q3)*(1.3E1/5.0E2)-sin(q2+q3)*(1.3E1/5.0E2)+cos(q2)*2.9575E-1+cos(q3)*(1.3E1/2.5E2)-sin(q3)*(1.3E1/2.5E2)+sqrt(4.1E1)*cos(q3-8.960553845713439E-1)*(1.3E1/1.0E3)+sqrt(4.1E1)*cos(8.960553845713439E-1)*(1.0/5.0E2)-sqrt(4.1E1)*sin(8.960553845713439E-1)*(1.0/5.0E2)+sqrt(4.1E1)*cos(q2+q3-8.960553845713439E-1)*6.5E-3+5.29375E-1;
  T(3,5) = cos(q2+q3)*(1.3E1/5.0E2)-sin(q2+q3)*(1.3E1/5.0E2)+cos(q3)*(1.3E1/5.0E2)-sin(q3)*(1.3E1/5.0E2)+sqrt(4.1E1)*cos(q3-8.960553845713439E-1)*6.5E-3+sqrt(4.1E1)*cos(8.960553845713439E-1)*(1.0/5.0E2)-sqrt(4.1E1)*sin(8.960553845713439E-1)*(1.0/5.0E2)+sqrt(4.1E1)*cos(q2+q3-8.960553845713439E-1)*6.5E-3+8.575E-2;
  T(4,0) = cos(q1+q2+theta)*(-9.1E1/1.0E2)-cos(q1+q2+q3+theta)*(2.0/2.5E1)+sin(q1+q2+q3+theta)*(2.0/2.5E1)-sqrt(4.1E1)*cos(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/5.0E1);
  T(4,1) = sin(q1+q2+theta)*(9.1E1/1.0E2)+cos(q1+q2+q3+theta)*(2.0/2.5E1)+sin(q1+q2+q3+theta)*(2.0/2.5E1)+sqrt(4.1E1)*sin(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/5.0E1);
  T(4,2) = cos(q2+q3)*(1.3E1/5.0E2)-sin(q2+q3)*(1.3E1/5.0E2)+cos(q2)*2.9575E-1+cos(q3)*(1.3E1/2.5E2)-sin(q3)*(1.3E1/2.5E2)+sqrt(4.1E1)*cos(q3-8.960553845713439E-1)*(1.3E1/1.0E3)+sqrt(4.1E1)*cos(8.960553845713439E-1)*(1.0/5.0E2)-sqrt(4.1E1)*sin(8.960553845713439E-1)*(1.0/5.0E2)+sqrt(4.1E1)*cos(q2+q3-8.960553845713439E-1)*6.5E-3+5.29375E-1;
  T(4,3) = cos(q2+q3)*(1.3E1/5.0E2)-sin(q2+q3)*(1.3E1/5.0E2)+cos(q2)*2.9575E-1+cos(q3)*(1.3E1/2.5E2)-sin(q3)*(1.3E1/2.5E2)+sqrt(4.1E1)*cos(q3-8.960553845713439E-1)*(1.3E1/1.0E3)+sqrt(4.1E1)*cos(8.960553845713439E-1)*(1.0/5.0E2)-sqrt(4.1E1)*sin(8.960553845713439E-1)*(1.0/5.0E2)+sqrt(4.1E1)*cos(q2+q3-8.960553845713439E-1)*6.5E-3+5.29375E-1;
  T(4,4) = cos(q3)*(1.3E1/2.5E2)-sin(q3)*(1.3E1/2.5E2)+sqrt(4.1E1)*cos(q3-8.960553845713439E-1)*(1.3E1/1.0E3)+sqrt(4.1E1)*cos(8.960553845713439E-1)*(1.0/5.0E2)-sqrt(4.1E1)*sin(8.960553845713439E-1)*(1.0/5.0E2)+5.29375E-1;
  T(4,5) = cos(q3)*(1.3E1/5.0E2)-sin(q3)*(1.3E1/5.0E2)+sqrt(4.1E1)*cos(q3-8.960553845713439E-1)*6.5E-3+sqrt(4.1E1)*cos(8.960553845713439E-1)*(1.0/5.0E2)-sqrt(4.1E1)*sin(8.960553845713439E-1)*(1.0/5.0E2)+8.575E-2;
  T(5,0) = cos(q1+q2+q3+theta)*(-2.0/2.5E1)+sin(q1+q2+q3+theta)*(2.0/2.5E1)-sqrt(4.1E1)*cos(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/5.0E1);
  T(5,1) = cos(q1+q2+q3+theta)*(2.0/2.5E1)+sin(q1+q2+q3+theta)*(2.0/2.5E1)+sqrt(4.1E1)*sin(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/5.0E1);
  T(5,2) = cos(q2+q3)*(1.3E1/5.0E2)-sin(q2+q3)*(1.3E1/5.0E2)+cos(q3)*(1.3E1/5.0E2)-sin(q3)*(1.3E1/5.0E2)+sqrt(4.1E1)*cos(q3-8.960553845713439E-1)*6.5E-3+sqrt(4.1E1)*cos(8.960553845713439E-1)*(1.0/5.0E2)-sqrt(4.1E1)*sin(8.960553845713439E-1)*(1.0/5.0E2)+sqrt(4.1E1)*cos(q2+q3-8.960553845713439E-1)*6.5E-3+8.575E-2;
  T(5,3) = cos(q2+q3)*(1.3E1/5.0E2)-sin(q2+q3)*(1.3E1/5.0E2)+cos(q3)*(1.3E1/5.0E2)-sin(q3)*(1.3E1/5.0E2)+sqrt(4.1E1)*cos(q3-8.960553845713439E-1)*6.5E-3+sqrt(4.1E1)*cos(8.960553845713439E-1)*(1.0/5.0E2)-sqrt(4.1E1)*sin(8.960553845713439E-1)*(1.0/5.0E2)+sqrt(4.1E1)*cos(q2+q3-8.960553845713439E-1)*6.5E-3+8.575E-2;
  T(5,4) = cos(q3)*(1.3E1/5.0E2)-sin(q3)*(1.3E1/5.0E2)+sqrt(4.1E1)*cos(q3-8.960553845713439E-1)*6.5E-3+sqrt(4.1E1)*cos(8.960553845713439E-1)*(1.0/5.0E2)-sqrt(4.1E1)*sin(8.960553845713439E-1)*(1.0/5.0E2)+8.575E-2;
  T(5,5) = sqrt(4.1E1)*cos(8.960553845713439E-1)*(1.0/5.0E2)-sqrt(4.1E1)*sin(8.960553845713439E-1)*(1.0/5.0E2)+8.575E-2;
  T(6,0) = cos(q4+q5+theta)*(-9.1E1/1.0E2)-cos(q4+q5+q6+theta)*(2.0/2.5E1)+sin(q4+q5+q6+theta)*(2.0/2.5E1)-cos(q4+theta)*(3.9E1/2.0E1)-sqrt(4.1E1)*cos(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/5.0E1);
  T(6,1) = sin(q4+q5+theta)*(9.1E1/1.0E2)+cos(q4+q5+q6+theta)*(2.0/2.5E1)+sin(q4+q5+q6+theta)*(2.0/2.5E1)+sin(q4+theta)*(3.9E1/2.0E1)+sqrt(4.1E1)*sin(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/5.0E1);
  T(6,2) = cos(q5+q6)*(1.3E1/2.5E2)-sin(q5+q6)*(1.3E1/2.5E2)+cos(q5)*5.915E-1+cos(q6)*(1.3E1/2.5E2)-sin(q6)*(1.3E1/2.5E2)+sqrt(4.1E1)*cos(q6-8.960553845713439E-1)*(1.3E1/1.0E3)+sqrt(4.1E1)*cos(8.960553845713439E-1)*(1.0/5.0E2)-sqrt(4.1E1)*sin(8.960553845713439E-1)*(1.0/5.0E2)+sqrt(4.1E1)*cos(q5+q6-8.960553845713439E-1)*(1.3E1/1.0E3)+1.409583333333333;
  T(6,6) = cos(q5+q6)*(1.3E1/2.5E2)-sin(q5+q6)*(1.3E1/2.5E2)+cos(q5)*5.915E-1+cos(q6)*(1.3E1/2.5E2)-sin(q6)*(1.3E1/2.5E2)+sqrt(4.1E1)*cos(q6-8.960553845713439E-1)*(1.3E1/1.0E3)+sqrt(4.1E1)*cos(8.960553845713439E-1)*(1.0/5.0E2)-sqrt(4.1E1)*sin(8.960553845713439E-1)*(1.0/5.0E2)+sqrt(4.1E1)*cos(q5+q6-8.960553845713439E-1)*(1.3E1/1.0E3)+1.409583333333333;
  T(6,7) = cos(q5+q6)*(1.3E1/5.0E2)-sin(q5+q6)*(1.3E1/5.0E2)+cos(q5)*2.9575E-1+cos(q6)*(1.3E1/2.5E2)-sin(q6)*(1.3E1/2.5E2)+sqrt(4.1E1)*cos(q6-8.960553845713439E-1)*(1.3E1/1.0E3)+sqrt(4.1E1)*cos(8.960553845713439E-1)*(1.0/5.0E2)-sqrt(4.1E1)*sin(8.960553845713439E-1)*(1.0/5.0E2)+sqrt(4.1E1)*cos(q5+q6-8.960553845713439E-1)*6.5E-3+5.29375E-1;
  T(6,8) = cos(q5+q6)*(1.3E1/5.0E2)-sin(q5+q6)*(1.3E1/5.0E2)+cos(q6)*(1.3E1/5.0E2)-sin(q6)*(1.3E1/5.0E2)+sqrt(4.1E1)*cos(q6-8.960553845713439E-1)*6.5E-3+sqrt(4.1E1)*cos(8.960553845713439E-1)*(1.0/5.0E2)-sqrt(4.1E1)*sin(8.960553845713439E-1)*(1.0/5.0E2)+sqrt(4.1E1)*cos(q5+q6-8.960553845713439E-1)*6.5E-3+8.575E-2;
  T(7,0) = cos(q4+q5+theta)*(-9.1E1/1.0E2)-cos(q4+q5+q6+theta)*(2.0/2.5E1)+sin(q4+q5+q6+theta)*(2.0/2.5E1)-sqrt(4.1E1)*cos(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/5.0E1);
  T(7,1) = sin(q4+q5+theta)*(9.1E1/1.0E2)+cos(q4+q5+q6+theta)*(2.0/2.5E1)+sin(q4+q5+q6+theta)*(2.0/2.5E1)+sqrt(4.1E1)*sin(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/5.0E1);
  T(7,2) = cos(q5+q6)*(1.3E1/5.0E2)-sin(q5+q6)*(1.3E1/5.0E2)+cos(q5)*2.9575E-1+cos(q6)*(1.3E1/2.5E2)-sin(q6)*(1.3E1/2.5E2)+sqrt(4.1E1)*cos(q6-8.960553845713439E-1)*(1.3E1/1.0E3)+sqrt(4.1E1)*cos(8.960553845713439E-1)*(1.0/5.0E2)-sqrt(4.1E1)*sin(8.960553845713439E-1)*(1.0/5.0E2)+sqrt(4.1E1)*cos(q5+q6-8.960553845713439E-1)*6.5E-3+5.29375E-1;
  T(7,6) = cos(q5+q6)*(1.3E1/5.0E2)-sin(q5+q6)*(1.3E1/5.0E2)+cos(q5)*2.9575E-1+cos(q6)*(1.3E1/2.5E2)-sin(q6)*(1.3E1/2.5E2)+sqrt(4.1E1)*cos(q6-8.960553845713439E-1)*(1.3E1/1.0E3)+sqrt(4.1E1)*cos(8.960553845713439E-1)*(1.0/5.0E2)-sqrt(4.1E1)*sin(8.960553845713439E-1)*(1.0/5.0E2)+sqrt(4.1E1)*cos(q5+q6-8.960553845713439E-1)*6.5E-3+5.29375E-1;
  T(7,7) = cos(q6)*(1.3E1/2.5E2)-sin(q6)*(1.3E1/2.5E2)+sqrt(4.1E1)*cos(q6-8.960553845713439E-1)*(1.3E1/1.0E3)+sqrt(4.1E1)*cos(8.960553845713439E-1)*(1.0/5.0E2)-sqrt(4.1E1)*sin(8.960553845713439E-1)*(1.0/5.0E2)+5.29375E-1;
  T(7,8) = cos(q6)*(1.3E1/5.0E2)-sin(q6)*(1.3E1/5.0E2)+sqrt(4.1E1)*cos(q6-8.960553845713439E-1)*6.5E-3+sqrt(4.1E1)*cos(8.960553845713439E-1)*(1.0/5.0E2)-sqrt(4.1E1)*sin(8.960553845713439E-1)*(1.0/5.0E2)+8.575E-2;
  T(8,0) = cos(q4+q5+q6+theta)*(-2.0/2.5E1)+sin(q4+q5+q6+theta)*(2.0/2.5E1)-sqrt(4.1E1)*cos(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/5.0E1);
  T(8,1) = cos(q4+q5+q6+theta)*(2.0/2.5E1)+sin(q4+q5+q6+theta)*(2.0/2.5E1)+sqrt(4.1E1)*sin(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/5.0E1);
  T(8,2) = cos(q5+q6)*(1.3E1/5.0E2)-sin(q5+q6)*(1.3E1/5.0E2)+cos(q6)*(1.3E1/5.0E2)-sin(q6)*(1.3E1/5.0E2)+sqrt(4.1E1)*cos(q6-8.960553845713439E-1)*6.5E-3+sqrt(4.1E1)*cos(8.960553845713439E-1)*(1.0/5.0E2)-sqrt(4.1E1)*sin(8.960553845713439E-1)*(1.0/5.0E2)+sqrt(4.1E1)*cos(q5+q6-8.960553845713439E-1)*6.5E-3+8.575E-2;
  T(8,6) = cos(q5+q6)*(1.3E1/5.0E2)-sin(q5+q6)*(1.3E1/5.0E2)+cos(q6)*(1.3E1/5.0E2)-sin(q6)*(1.3E1/5.0E2)+sqrt(4.1E1)*cos(q6-8.960553845713439E-1)*6.5E-3+sqrt(4.1E1)*cos(8.960553845713439E-1)*(1.0/5.0E2)-sqrt(4.1E1)*sin(8.960553845713439E-1)*(1.0/5.0E2)+sqrt(4.1E1)*cos(q5+q6-8.960553845713439E-1)*6.5E-3+8.575E-2;
  T(8,7) = cos(q6)*(1.3E1/5.0E2)-sin(q6)*(1.3E1/5.0E2)+sqrt(4.1E1)*cos(q6-8.960553845713439E-1)*6.5E-3+sqrt(4.1E1)*cos(8.960553845713439E-1)*(1.0/5.0E2)-sqrt(4.1E1)*sin(8.960553845713439E-1)*(1.0/5.0E2)+8.575E-2;
  T(8,8) = sqrt(4.1E1)*cos(8.960553845713439E-1)*(1.0/5.0E2)-sqrt(4.1E1)*sin(8.960553845713439E-1)*(1.0/5.0E2)+8.575E-2;
  T(9,0) = cos(q7+q8+theta)*(-9.9E1/1.6E2)-cos(q7+theta)*(1.9E1/1.6E1);
  T(9,1) = sin(q7+q8+theta)*(9.9E1/1.6E2)+sin(q7+theta)*(1.9E1/1.6E1);
  T(9,2) = cos(q7+q8)*(-3.403125E-1)-cos(q7)*(2.09E2/3.2E2)+cos(q8)*(9.9E1/3.2E2)+1.045989583333333;
  T(9,9) = cos(q8)*(9.9E1/3.2E2)+1.045989583333333;
  T(9,10) = cos(q8)*(9.9E1/6.4E2)+6.0328125E-1;
  T(10,0) = cos(q7+q8+theta)*(-9.9E1/1.6E2);
  T(10,1) = sin(q7+q8+theta)*(9.9E1/1.6E2);
  T(10,2) = cos(q7+q8)*(-3.403125E-1)+cos(q8)*(9.9E1/6.4E2)+6.0328125E-1;
  T(10,9) = cos(q8)*(9.9E1/6.4E2)+6.0328125E-1;
  T(10,10) = 6.0328125E-1;
  T(11,0) = cos(q9+q10+theta)*(-9.9E1/1.6E2)-cos(q9+theta)*(1.9E1/1.6E1);
  T(11,1) = sin(q9+q10+theta)*(9.9E1/1.6E2)+sin(q9+theta)*(1.9E1/1.6E1);
  T(11,2) = cos(q9+q10)*(-3.403125E-1)-cos(q9)*(2.09E2/3.2E2)+cos(q10)*(9.9E1/3.2E2)+1.045989583333333;
  T(11,11) = cos(q10)*(9.9E1/3.2E2)+1.045989583333333;
  T(11,12) = cos(q10)*(9.9E1/6.4E2)+6.0328125E-1;
  T(12,0) = cos(q9+q10+theta)*(-9.9E1/1.6E2);
  T(12,1) = sin(q9+q10+theta)*(9.9E1/1.6E2);
  T(12,2) = cos(q9+q10)*(-3.403125E-1)+cos(q10)*(9.9E1/6.4E2)+6.0328125E-1;
  T(12,11) = cos(q10)*(9.9E1/6.4E2)+6.0328125E-1;
  T(12,12) = 6.0328125E-1;

  return T;
}

dlib::matrix<double> B_q_fn()
{
  dlib::matrix<double> B_q;
  B_q = dlib::zeros_matrix<double>(13,10);
  B_q(3,0) = 1.0;
  B_q(4,1) = 1.0;
  B_q(5,2) = 1.0;
  B_q(6,3) = 1.0;
  B_q(7,4) = 1.0;
  B_q(8,5) = 1.0;
  B_q(9,6) = 1.0;
  B_q(10,7) = 1.0;
  B_q(11,8) = 1.0;
  B_q(12,9) = 1.0;
  return B_q;
}

dlib::matrix<double> C_q_qdot_fn(const Robot_StateNDot &Robot_StateNDot_i)
{
  double rIx = Robot_StateNDot_i.rIx;					double rIy = Robot_StateNDot_i.rIy;						double theta = Robot_StateNDot_i.theta;
  double q1 = Robot_StateNDot_i.q1;					double q2 = Robot_StateNDot_i.q2;						double q3 = Robot_StateNDot_i.q3;
  double q4 = Robot_StateNDot_i.q4;					double q5 = Robot_StateNDot_i.q5;						double q6 = Robot_StateNDot_i.q6;
  double q7 = Robot_StateNDot_i.q7;					double q8 = Robot_StateNDot_i.q8;						double q9 = Robot_StateNDot_i.q9;
  double q10 = Robot_StateNDot_i.q10;

  double rIxdot = Robot_StateNDot_i.rIxdot;			double rIydot = Robot_StateNDot_i.rIydot;				double thetadot = Robot_StateNDot_i.thetadot;
  double q1dot = Robot_StateNDot_i.q1dot;				double q2dot = Robot_StateNDot_i.q2dot;					double q3dot = Robot_StateNDot_i.q3dot;
  double q4dot = Robot_StateNDot_i.q4dot;				double q5dot = Robot_StateNDot_i.q5dot;					double q6dot = Robot_StateNDot_i.q6dot;
  double q7dot = Robot_StateNDot_i.q7dot;				double q8dot = Robot_StateNDot_i.q8dot;					double q9dot = Robot_StateNDot_i.q9dot;
  double q10dot = Robot_StateNDot_i.q10dot;

  dlib::matrix<double> T;
  T = dlib::zeros_matrix<double>(13,1);
  T(0,0) = (thetadot*thetadot)*sin(theta)*(-2.97E2/2.0E1)+(q10dot*q10dot)*sin(q9+q10+theta)*(9.9E1/1.6E2)+(q1dot*q1dot)*sin(q1+q2+theta)*(9.1E1/1.0E2)+(q2dot*q2dot)*sin(q1+q2+theta)*(9.1E1/1.0E2)+(q4dot*q4dot)*sin(q4+q5+theta)*(9.1E1/1.0E2)+(q5dot*q5dot)*sin(q4+q5+theta)*(9.1E1/1.0E2)+(q7dot*q7dot)*sin(q7+q8+theta)*(9.9E1/1.6E2)+(q8dot*q8dot)*sin(q7+q8+theta)*(9.9E1/1.6E2)+(q9dot*q9dot)*sin(q9+q10+theta)*(9.9E1/1.6E2)+(thetadot*thetadot)*sin(q1+q2+theta)*(9.1E1/1.0E2)+(thetadot*thetadot)*sin(q4+q5+theta)*(9.1E1/1.0E2)+(thetadot*thetadot)*sin(q7+q8+theta)*(9.9E1/1.6E2)+(thetadot*thetadot)*sin(q9+q10+theta)*(9.9E1/1.6E2)+(q1dot*q1dot)*cos(q1+q2+q3+theta)*(2.0/2.5E1)+(q2dot*q2dot)*cos(q1+q2+q3+theta)*(2.0/2.5E1)+(q3dot*q3dot)*cos(q1+q2+q3+theta)*(2.0/2.5E1)+(q4dot*q4dot)*cos(q4+q5+q6+theta)*(2.0/2.5E1)+(q5dot*q5dot)*cos(q4+q5+q6+theta)*(2.0/2.5E1)+(q6dot*q6dot)*cos(q4+q5+q6+theta)*(2.0/2.5E1)+(thetadot*thetadot)*cos(q1+q2+q3+theta)*(2.0/2.5E1)+(thetadot*thetadot)*cos(q4+q5+q6+theta)*(2.0/2.5E1)+(q1dot*q1dot)*sin(q1+q2+q3+theta)*(2.0/2.5E1)+(q2dot*q2dot)*sin(q1+q2+q3+theta)*(2.0/2.5E1)+(q3dot*q3dot)*sin(q1+q2+q3+theta)*(2.0/2.5E1)+(q4dot*q4dot)*sin(q4+q5+q6+theta)*(2.0/2.5E1)+(q5dot*q5dot)*sin(q4+q5+q6+theta)*(2.0/2.5E1)+(q6dot*q6dot)*sin(q4+q5+q6+theta)*(2.0/2.5E1)+(thetadot*thetadot)*sin(q1+q2+q3+theta)*(2.0/2.5E1)+(thetadot*thetadot)*sin(q4+q5+q6+theta)*(2.0/2.5E1)+(q1dot*q1dot)*sin(q1+theta)*(3.9E1/2.0E1)+(q4dot*q4dot)*sin(q4+theta)*(3.9E1/2.0E1)+(q7dot*q7dot)*sin(q7+theta)*(1.9E1/1.6E1)+(q9dot*q9dot)*sin(q9+theta)*(1.9E1/1.6E1)+(thetadot*thetadot)*sin(q1+theta)*(3.9E1/2.0E1)+(thetadot*thetadot)*sin(q4+theta)*(3.9E1/2.0E1)+(thetadot*thetadot)*sin(q7+theta)*(1.9E1/1.6E1)+(thetadot*thetadot)*sin(q9+theta)*(1.9E1/1.6E1)+sqrt(4.1E1)*(q1dot*q1dot)*sin(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/5.0E1)+sqrt(4.1E1)*(q2dot*q2dot)*sin(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/5.0E1)+sqrt(4.1E1)*(q3dot*q3dot)*sin(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/5.0E1)+sqrt(4.1E1)*(q4dot*q4dot)*sin(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/5.0E1)+sqrt(4.1E1)*(q5dot*q5dot)*sin(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/5.0E1)+sqrt(4.1E1)*(q6dot*q6dot)*sin(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/5.0E1)+sqrt(4.1E1)*(thetadot*thetadot)*sin(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/5.0E1)+sqrt(4.1E1)*(thetadot*thetadot)*sin(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/5.0E1)+q10dot*q9dot*sin(q9+q10+theta)*(9.9E1/8.0E1)+q1dot*q2dot*sin(q1+q2+theta)*(9.1E1/5.0E1)+q4dot*q5dot*sin(q4+q5+theta)*(9.1E1/5.0E1)+q7dot*q8dot*sin(q7+q8+theta)*(9.9E1/8.0E1)+q10dot*thetadot*sin(q9+q10+theta)*(9.9E1/8.0E1)+q1dot*thetadot*sin(q1+q2+theta)*(9.1E1/5.0E1)+q2dot*thetadot*sin(q1+q2+theta)*(9.1E1/5.0E1)+q4dot*thetadot*sin(q4+q5+theta)*(9.1E1/5.0E1)+q5dot*thetadot*sin(q4+q5+theta)*(9.1E1/5.0E1)+q7dot*thetadot*sin(q7+q8+theta)*(9.9E1/8.0E1)+q8dot*thetadot*sin(q7+q8+theta)*(9.9E1/8.0E1)+q9dot*thetadot*sin(q9+q10+theta)*(9.9E1/8.0E1)+q1dot*q2dot*cos(q1+q2+q3+theta)*(4.0/2.5E1)+q1dot*q3dot*cos(q1+q2+q3+theta)*(4.0/2.5E1)+q2dot*q3dot*cos(q1+q2+q3+theta)*(4.0/2.5E1)+q4dot*q5dot*cos(q4+q5+q6+theta)*(4.0/2.5E1)+q4dot*q6dot*cos(q4+q5+q6+theta)*(4.0/2.5E1)+q5dot*q6dot*cos(q4+q5+q6+theta)*(4.0/2.5E1)+q1dot*thetadot*cos(q1+q2+q3+theta)*(4.0/2.5E1)+q2dot*thetadot*cos(q1+q2+q3+theta)*(4.0/2.5E1)+q3dot*thetadot*cos(q1+q2+q3+theta)*(4.0/2.5E1)+q4dot*thetadot*cos(q4+q5+q6+theta)*(4.0/2.5E1)+q5dot*thetadot*cos(q4+q5+q6+theta)*(4.0/2.5E1)+q6dot*thetadot*cos(q4+q5+q6+theta)*(4.0/2.5E1)+q1dot*q2dot*sin(q1+q2+q3+theta)*(4.0/2.5E1)+q1dot*q3dot*sin(q1+q2+q3+theta)*(4.0/2.5E1)+q2dot*q3dot*sin(q1+q2+q3+theta)*(4.0/2.5E1)+q4dot*q5dot*sin(q4+q5+q6+theta)*(4.0/2.5E1)+q4dot*q6dot*sin(q4+q5+q6+theta)*(4.0/2.5E1)+q5dot*q6dot*sin(q4+q5+q6+theta)*(4.0/2.5E1)+q1dot*thetadot*sin(q1+q2+q3+theta)*(4.0/2.5E1)+q2dot*thetadot*sin(q1+q2+q3+theta)*(4.0/2.5E1)+q3dot*thetadot*sin(q1+q2+q3+theta)*(4.0/2.5E1)+q4dot*thetadot*sin(q4+q5+q6+theta)*(4.0/2.5E1)+q5dot*thetadot*sin(q4+q5+q6+theta)*(4.0/2.5E1)+q6dot*thetadot*sin(q4+q5+q6+theta)*(4.0/2.5E1)+q1dot*thetadot*sin(q1+theta)*(3.9E1/1.0E1)+q4dot*thetadot*sin(q4+theta)*(3.9E1/1.0E1)+q7dot*thetadot*sin(q7+theta)*(1.9E1/8.0)+q9dot*thetadot*sin(q9+theta)*(1.9E1/8.0)+sqrt(4.1E1)*q1dot*q2dot*sin(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/2.5E1)+sqrt(4.1E1)*q1dot*q3dot*sin(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/2.5E1)+sqrt(4.1E1)*q2dot*q3dot*sin(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/2.5E1)+sqrt(4.1E1)*q4dot*q5dot*sin(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/2.5E1)+sqrt(4.1E1)*q4dot*q6dot*sin(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/2.5E1)+sqrt(4.1E1)*q5dot*q6dot*sin(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/2.5E1)+sqrt(4.1E1)*q1dot*thetadot*sin(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/2.5E1)+sqrt(4.1E1)*q2dot*thetadot*sin(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/2.5E1)+sqrt(4.1E1)*q3dot*thetadot*sin(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/2.5E1)+sqrt(4.1E1)*q4dot*thetadot*sin(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/2.5E1)+sqrt(4.1E1)*q5dot*thetadot*sin(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/2.5E1)+sqrt(4.1E1)*q6dot*thetadot*sin(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/2.5E1);
  T(1,0) = (thetadot*thetadot)*cos(theta)*(-2.97E2/2.0E1)+(q10dot*q10dot)*cos(q9+q10+theta)*(9.9E1/1.6E2)+(q1dot*q1dot)*cos(q1+q2+theta)*(9.1E1/1.0E2)+(q2dot*q2dot)*cos(q1+q2+theta)*(9.1E1/1.0E2)+(q4dot*q4dot)*cos(q4+q5+theta)*(9.1E1/1.0E2)+(q5dot*q5dot)*cos(q4+q5+theta)*(9.1E1/1.0E2)+(q7dot*q7dot)*cos(q7+q8+theta)*(9.9E1/1.6E2)+(q8dot*q8dot)*cos(q7+q8+theta)*(9.9E1/1.6E2)+(q9dot*q9dot)*cos(q9+q10+theta)*(9.9E1/1.6E2)+(thetadot*thetadot)*cos(q1+q2+theta)*(9.1E1/1.0E2)+(thetadot*thetadot)*cos(q4+q5+theta)*(9.1E1/1.0E2)+(thetadot*thetadot)*cos(q7+q8+theta)*(9.9E1/1.6E2)+(thetadot*thetadot)*cos(q9+q10+theta)*(9.9E1/1.6E2)+(q1dot*q1dot)*cos(q1+q2+q3+theta)*(2.0/2.5E1)+(q2dot*q2dot)*cos(q1+q2+q3+theta)*(2.0/2.5E1)+(q3dot*q3dot)*cos(q1+q2+q3+theta)*(2.0/2.5E1)+(q4dot*q4dot)*cos(q4+q5+q6+theta)*(2.0/2.5E1)+(q5dot*q5dot)*cos(q4+q5+q6+theta)*(2.0/2.5E1)+(q6dot*q6dot)*cos(q4+q5+q6+theta)*(2.0/2.5E1)+(thetadot*thetadot)*cos(q1+q2+q3+theta)*(2.0/2.5E1)+(thetadot*thetadot)*cos(q4+q5+q6+theta)*(2.0/2.5E1)-(q1dot*q1dot)*sin(q1+q2+q3+theta)*(2.0/2.5E1)-(q2dot*q2dot)*sin(q1+q2+q3+theta)*(2.0/2.5E1)-(q3dot*q3dot)*sin(q1+q2+q3+theta)*(2.0/2.5E1)-(q4dot*q4dot)*sin(q4+q5+q6+theta)*(2.0/2.5E1)-(q5dot*q5dot)*sin(q4+q5+q6+theta)*(2.0/2.5E1)-(q6dot*q6dot)*sin(q4+q5+q6+theta)*(2.0/2.5E1)-(thetadot*thetadot)*sin(q1+q2+q3+theta)*(2.0/2.5E1)-(thetadot*thetadot)*sin(q4+q5+q6+theta)*(2.0/2.5E1)+(q1dot*q1dot)*cos(q1+theta)*(3.9E1/2.0E1)+(q4dot*q4dot)*cos(q4+theta)*(3.9E1/2.0E1)+(q7dot*q7dot)*cos(q7+theta)*(1.9E1/1.6E1)+(q9dot*q9dot)*cos(q9+theta)*(1.9E1/1.6E1)+(thetadot*thetadot)*cos(q1+theta)*(3.9E1/2.0E1)+(thetadot*thetadot)*cos(q4+theta)*(3.9E1/2.0E1)+(thetadot*thetadot)*cos(q7+theta)*(1.9E1/1.6E1)+(thetadot*thetadot)*cos(q9+theta)*(1.9E1/1.6E1)+sqrt(4.1E1)*(q1dot*q1dot)*cos(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/5.0E1)+sqrt(4.1E1)*(q2dot*q2dot)*cos(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/5.0E1)+sqrt(4.1E1)*(q3dot*q3dot)*cos(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/5.0E1)+sqrt(4.1E1)*(q4dot*q4dot)*cos(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/5.0E1)+sqrt(4.1E1)*(q5dot*q5dot)*cos(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/5.0E1)+sqrt(4.1E1)*(q6dot*q6dot)*cos(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/5.0E1)+sqrt(4.1E1)*(thetadot*thetadot)*cos(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/5.0E1)+sqrt(4.1E1)*(thetadot*thetadot)*cos(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/5.0E1)+q10dot*q9dot*cos(q9+q10+theta)*(9.9E1/8.0E1)+q1dot*q2dot*cos(q1+q2+theta)*(9.1E1/5.0E1)+q4dot*q5dot*cos(q4+q5+theta)*(9.1E1/5.0E1)+q7dot*q8dot*cos(q7+q8+theta)*(9.9E1/8.0E1)+q10dot*thetadot*cos(q9+q10+theta)*(9.9E1/8.0E1)+q1dot*thetadot*cos(q1+q2+theta)*(9.1E1/5.0E1)+q2dot*thetadot*cos(q1+q2+theta)*(9.1E1/5.0E1)+q4dot*thetadot*cos(q4+q5+theta)*(9.1E1/5.0E1)+q5dot*thetadot*cos(q4+q5+theta)*(9.1E1/5.0E1)+q7dot*thetadot*cos(q7+q8+theta)*(9.9E1/8.0E1)+q8dot*thetadot*cos(q7+q8+theta)*(9.9E1/8.0E1)+q9dot*thetadot*cos(q9+q10+theta)*(9.9E1/8.0E1)+q1dot*q2dot*cos(q1+q2+q3+theta)*(4.0/2.5E1)+q1dot*q3dot*cos(q1+q2+q3+theta)*(4.0/2.5E1)+q2dot*q3dot*cos(q1+q2+q3+theta)*(4.0/2.5E1)+q4dot*q5dot*cos(q4+q5+q6+theta)*(4.0/2.5E1)+q4dot*q6dot*cos(q4+q5+q6+theta)*(4.0/2.5E1)+q5dot*q6dot*cos(q4+q5+q6+theta)*(4.0/2.5E1)+q1dot*thetadot*cos(q1+q2+q3+theta)*(4.0/2.5E1)+q2dot*thetadot*cos(q1+q2+q3+theta)*(4.0/2.5E1)+q3dot*thetadot*cos(q1+q2+q3+theta)*(4.0/2.5E1)+q4dot*thetadot*cos(q4+q5+q6+theta)*(4.0/2.5E1)+q5dot*thetadot*cos(q4+q5+q6+theta)*(4.0/2.5E1)+q6dot*thetadot*cos(q4+q5+q6+theta)*(4.0/2.5E1)-q1dot*q2dot*sin(q1+q2+q3+theta)*(4.0/2.5E1)-q1dot*q3dot*sin(q1+q2+q3+theta)*(4.0/2.5E1)-q2dot*q3dot*sin(q1+q2+q3+theta)*(4.0/2.5E1)-q4dot*q5dot*sin(q4+q5+q6+theta)*(4.0/2.5E1)-q4dot*q6dot*sin(q4+q5+q6+theta)*(4.0/2.5E1)-q5dot*q6dot*sin(q4+q5+q6+theta)*(4.0/2.5E1)-q1dot*thetadot*sin(q1+q2+q3+theta)*(4.0/2.5E1)-q2dot*thetadot*sin(q1+q2+q3+theta)*(4.0/2.5E1)-q3dot*thetadot*sin(q1+q2+q3+theta)*(4.0/2.5E1)-q4dot*thetadot*sin(q4+q5+q6+theta)*(4.0/2.5E1)-q5dot*thetadot*sin(q4+q5+q6+theta)*(4.0/2.5E1)-q6dot*thetadot*sin(q4+q5+q6+theta)*(4.0/2.5E1)+q1dot*thetadot*cos(q1+theta)*(3.9E1/1.0E1)+q4dot*thetadot*cos(q4+theta)*(3.9E1/1.0E1)+q7dot*thetadot*cos(q7+theta)*(1.9E1/8.0)+q9dot*thetadot*cos(q9+theta)*(1.9E1/8.0)+sqrt(4.1E1)*q1dot*q2dot*cos(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/2.5E1)+sqrt(4.1E1)*q1dot*q3dot*cos(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/2.5E1)+sqrt(4.1E1)*q2dot*q3dot*cos(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/2.5E1)+sqrt(4.1E1)*q4dot*q5dot*cos(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/2.5E1)+sqrt(4.1E1)*q4dot*q6dot*cos(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/2.5E1)+sqrt(4.1E1)*q5dot*q6dot*cos(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/2.5E1)+sqrt(4.1E1)*q1dot*thetadot*cos(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/2.5E1)+sqrt(4.1E1)*q2dot*thetadot*cos(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/2.5E1)+sqrt(4.1E1)*q3dot*thetadot*cos(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/2.5E1)+sqrt(4.1E1)*q4dot*thetadot*cos(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/2.5E1)+sqrt(4.1E1)*q5dot*thetadot*cos(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/2.5E1)+sqrt(4.1E1)*q6dot*thetadot*cos(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/2.5E1)+5.54265E2;
  T(2,0) = sin(q1+q2+theta)*8.9271+sin(q4+q5+theta)*8.9271+sin(q7+q8+theta)*6.0699375+sin(q9+q10+theta)*6.0699375+sin(q1+theta)*1.91295E1+sin(q4+theta)*1.91295E1+sin(q7+theta)*1.1649375E1+sin(q9+theta)*1.1649375E1-sin(theta)*1.456785E2-(q3dot*q3dot)*cos(q3)*(1.3E1/5.0E2)-(q6dot*q6dot)*cos(q6)*(1.3E1/5.0E2)+sqrt(2.0)*sin(q1+q2+q3+theta+3.141592653589793*(1.0/4.0))*7.848E-1+sqrt(2.0)*sin(q4+q5+q6+theta+3.141592653589793*(1.0/4.0))*7.848E-1-(q10dot*q10dot)*sin(q10)*(9.9E1/6.4E2)-(q2dot*q2dot)*sin(q2)*2.9575E-1-(q3dot*q3dot)*sin(q3)*(1.3E1/5.0E2)-(q5dot*q5dot)*sin(q5)*2.9575E-1-(q6dot*q6dot)*sin(q6)*(1.3E1/5.0E2)+(q7dot*q7dot)*sin(q7)*(2.09E2/3.2E2)-(q8dot*q8dot)*sin(q8)*(9.9E1/6.4E2)+(q9dot*q9dot)*sin(q9)*(2.09E2/3.2E2)+sqrt(4.1E1)*sin(q1+q2+q3+theta-8.960553845713439E-1)*1.962E-1+sqrt(4.1E1)*sin(q4+q5+q6+theta-8.960553845713439E-1)*1.962E-1-(q2dot*q2dot)*cos(q2+q3)*(1.3E1/5.0E2)-(q3dot*q3dot)*cos(q2+q3)*(1.3E1/5.0E2)-(q5dot*q5dot)*cos(q5+q6)*(1.3E1/5.0E2)-(q6dot*q6dot)*cos(q5+q6)*(1.3E1/5.0E2)+(q10dot*q10dot)*sin(q9+q10)*3.403125E-1-(q2dot*q2dot)*sin(q2+q3)*(1.3E1/5.0E2)-(q3dot*q3dot)*sin(q2+q3)*(1.3E1/5.0E2)-(q5dot*q5dot)*sin(q5+q6)*(1.3E1/5.0E2)-(q6dot*q6dot)*sin(q5+q6)*(1.3E1/5.0E2)+(q7dot*q7dot)*sin(q7+q8)*3.403125E-1+(q8dot*q8dot)*sin(q7+q8)*3.403125E-1+(q9dot*q9dot)*sin(q9+q10)*3.403125E-1-q1dot*q3dot*cos(q3)*(1.3E1/2.5E2)-q2dot*q3dot*cos(q3)*(1.3E1/2.5E2)-q4dot*q6dot*cos(q6)*(1.3E1/2.5E2)-q5dot*q6dot*cos(q6)*(1.3E1/2.5E2)-q3dot*thetadot*cos(q3)*(1.3E1/2.5E2)-q6dot*thetadot*cos(q6)*(1.3E1/2.5E2)-q10dot*q9dot*sin(q10)*(9.9E1/3.2E2)-q1dot*q2dot*sin(q2)*5.915E-1-q1dot*q3dot*sin(q3)*(1.3E1/2.5E2)-q2dot*q3dot*sin(q3)*(1.3E1/2.5E2)-q4dot*q5dot*sin(q5)*5.915E-1-q4dot*q6dot*sin(q6)*(1.3E1/2.5E2)-q5dot*q6dot*sin(q6)*(1.3E1/2.5E2)-q7dot*q8dot*sin(q8)*(9.9E1/3.2E2)-q10dot*thetadot*sin(q10)*(9.9E1/3.2E2)-q2dot*thetadot*sin(q2)*5.915E-1-q3dot*thetadot*sin(q3)*(1.3E1/2.5E2)-q5dot*thetadot*sin(q5)*5.915E-1-q6dot*thetadot*sin(q6)*(1.3E1/2.5E2)+q7dot*thetadot*sin(q7)*(2.09E2/1.6E2)-q8dot*thetadot*sin(q8)*(9.9E1/3.2E2)+q9dot*thetadot*sin(q9)*(2.09E2/1.6E2)-sqrt(4.1E1)*(q2dot*q2dot)*sin(q2+q3-8.960553845713439E-1)*6.5E-3-sqrt(4.1E1)*(q3dot*q3dot)*sin(q2+q3-8.960553845713439E-1)*6.5E-3-sqrt(4.1E1)*(q5dot*q5dot)*sin(q5+q6-8.960553845713439E-1)*6.5E-3-sqrt(4.1E1)*(q6dot*q6dot)*sin(q5+q6-8.960553845713439E-1)*6.5E-3-q1dot*q2dot*cos(q2+q3)*(1.3E1/2.5E2)-q1dot*q3dot*cos(q2+q3)*(1.3E1/2.5E2)-q2dot*q3dot*cos(q2+q3)*(1.3E1/2.5E2)-q4dot*q5dot*cos(q5+q6)*(1.3E1/2.5E2)-q4dot*q6dot*cos(q5+q6)*(1.3E1/2.5E2)-q5dot*q6dot*cos(q5+q6)*(1.3E1/2.5E2)-q2dot*thetadot*cos(q2+q3)*(1.3E1/2.5E2)-q3dot*thetadot*cos(q2+q3)*(1.3E1/2.5E2)-q5dot*thetadot*cos(q5+q6)*(1.3E1/2.5E2)-q6dot*thetadot*cos(q5+q6)*(1.3E1/2.5E2)+q10dot*q9dot*sin(q9+q10)*6.80625E-1-q1dot*q2dot*sin(q2+q3)*(1.3E1/2.5E2)-q1dot*q3dot*sin(q2+q3)*(1.3E1/2.5E2)-q2dot*q3dot*sin(q2+q3)*(1.3E1/2.5E2)-q4dot*q5dot*sin(q5+q6)*(1.3E1/2.5E2)-q4dot*q6dot*sin(q5+q6)*(1.3E1/2.5E2)-q5dot*q6dot*sin(q5+q6)*(1.3E1/2.5E2)+q7dot*q8dot*sin(q7+q8)*6.80625E-1+q10dot*thetadot*sin(q9+q10)*6.80625E-1-q2dot*thetadot*sin(q2+q3)*(1.3E1/2.5E2)-q3dot*thetadot*sin(q2+q3)*(1.3E1/2.5E2)-q5dot*thetadot*sin(q5+q6)*(1.3E1/2.5E2)-q6dot*thetadot*sin(q5+q6)*(1.3E1/2.5E2)+q7dot*thetadot*sin(q7+q8)*6.80625E-1+q8dot*thetadot*sin(q7+q8)*6.80625E-1+q9dot*thetadot*sin(q9+q10)*6.80625E-1-sqrt(4.1E1)*(q3dot*q3dot)*sin(q3-8.960553845713439E-1)*6.5E-3-sqrt(4.1E1)*(q6dot*q6dot)*sin(q6-8.960553845713439E-1)*6.5E-3-sqrt(4.1E1)*q1dot*q3dot*sin(q3-8.960553845713439E-1)*(1.3E1/1.0E3)-sqrt(4.1E1)*q2dot*q3dot*sin(q3-8.960553845713439E-1)*(1.3E1/1.0E3)-sqrt(4.1E1)*q4dot*q6dot*sin(q6-8.960553845713439E-1)*(1.3E1/1.0E3)-sqrt(4.1E1)*q5dot*q6dot*sin(q6-8.960553845713439E-1)*(1.3E1/1.0E3)-sqrt(4.1E1)*q3dot*thetadot*sin(q3-8.960553845713439E-1)*(1.3E1/1.0E3)-sqrt(4.1E1)*q6dot*thetadot*sin(q6-8.960553845713439E-1)*(1.3E1/1.0E3)-sqrt(4.1E1)*q1dot*q2dot*sin(q2+q3-8.960553845713439E-1)*(1.3E1/1.0E3)-sqrt(4.1E1)*q1dot*q3dot*sin(q2+q3-8.960553845713439E-1)*(1.3E1/1.0E3)-sqrt(4.1E1)*q2dot*q3dot*sin(q2+q3-8.960553845713439E-1)*(1.3E1/1.0E3)-sqrt(4.1E1)*q4dot*q5dot*sin(q5+q6-8.960553845713439E-1)*(1.3E1/1.0E3)-sqrt(4.1E1)*q4dot*q6dot*sin(q5+q6-8.960553845713439E-1)*(1.3E1/1.0E3)-sqrt(4.1E1)*q5dot*q6dot*sin(q5+q6-8.960553845713439E-1)*(1.3E1/1.0E3)-sqrt(4.1E1)*q2dot*thetadot*sin(q2+q3-8.960553845713439E-1)*(1.3E1/1.0E3)-sqrt(4.1E1)*q3dot*thetadot*sin(q2+q3-8.960553845713439E-1)*(1.3E1/1.0E3)-sqrt(4.1E1)*q5dot*thetadot*sin(q5+q6-8.960553845713439E-1)*(1.3E1/1.0E3)-sqrt(4.1E1)*q6dot*thetadot*sin(q5+q6-8.960553845713439E-1)*(1.3E1/1.0E3);
  T(3,0) = sin(q1+q2+theta)*8.9271+sin(q1+theta)*1.91295E1-(q3dot*q3dot)*cos(q3)*(1.3E1/5.0E2)+sqrt(2.0)*sin(q1+q2+q3+theta+3.141592653589793*(1.0/4.0))*7.848E-1-(q2dot*q2dot)*sin(q2)*2.9575E-1-(q3dot*q3dot)*sin(q3)*(1.3E1/5.0E2)+sqrt(4.1E1)*sin(q1+q2+q3+theta-8.960553845713439E-1)*1.962E-1-(q2dot*q2dot)*cos(q2+q3)*(1.3E1/5.0E2)-(q3dot*q3dot)*cos(q2+q3)*(1.3E1/5.0E2)-(q2dot*q2dot)*sin(q2+q3)*(1.3E1/5.0E2)-(q3dot*q3dot)*sin(q2+q3)*(1.3E1/5.0E2)-q1dot*q3dot*cos(q3)*(1.3E1/2.5E2)-q2dot*q3dot*cos(q3)*(1.3E1/2.5E2)-q3dot*thetadot*cos(q3)*(1.3E1/2.5E2)-q1dot*q2dot*sin(q2)*5.915E-1-q1dot*q3dot*sin(q3)*(1.3E1/2.5E2)-q2dot*q3dot*sin(q3)*(1.3E1/2.5E2)-q2dot*thetadot*sin(q2)*5.915E-1-q3dot*thetadot*sin(q3)*(1.3E1/2.5E2)-sqrt(4.1E1)*(q2dot*q2dot)*sin(q2+q3-8.960553845713439E-1)*6.5E-3-sqrt(4.1E1)*(q3dot*q3dot)*sin(q2+q3-8.960553845713439E-1)*6.5E-3-q1dot*q2dot*cos(q2+q3)*(1.3E1/2.5E2)-q1dot*q3dot*cos(q2+q3)*(1.3E1/2.5E2)-q2dot*q3dot*cos(q2+q3)*(1.3E1/2.5E2)-q2dot*thetadot*cos(q2+q3)*(1.3E1/2.5E2)-q3dot*thetadot*cos(q2+q3)*(1.3E1/2.5E2)-q1dot*q2dot*sin(q2+q3)*(1.3E1/2.5E2)-q1dot*q3dot*sin(q2+q3)*(1.3E1/2.5E2)-q2dot*q3dot*sin(q2+q3)*(1.3E1/2.5E2)-q2dot*thetadot*sin(q2+q3)*(1.3E1/2.5E2)-q3dot*thetadot*sin(q2+q3)*(1.3E1/2.5E2)-sqrt(4.1E1)*(q3dot*q3dot)*sin(q3-8.960553845713439E-1)*6.5E-3-sqrt(4.1E1)*q1dot*q3dot*sin(q3-8.960553845713439E-1)*(1.3E1/1.0E3)-sqrt(4.1E1)*q2dot*q3dot*sin(q3-8.960553845713439E-1)*(1.3E1/1.0E3)-sqrt(4.1E1)*q3dot*thetadot*sin(q3-8.960553845713439E-1)*(1.3E1/1.0E3)-sqrt(4.1E1)*q1dot*q2dot*sin(q2+q3-8.960553845713439E-1)*(1.3E1/1.0E3)-sqrt(4.1E1)*q1dot*q3dot*sin(q2+q3-8.960553845713439E-1)*(1.3E1/1.0E3)-sqrt(4.1E1)*q2dot*q3dot*sin(q2+q3-8.960553845713439E-1)*(1.3E1/1.0E3)-sqrt(4.1E1)*q2dot*thetadot*sin(q2+q3-8.960553845713439E-1)*(1.3E1/1.0E3)-sqrt(4.1E1)*q3dot*thetadot*sin(q2+q3-8.960553845713439E-1)*(1.3E1/1.0E3);
  T(4,0) = sin(q1+q2+theta)*8.9271-(q3dot*q3dot)*cos(q3)*(1.3E1/5.0E2)+sqrt(2.0)*sin(q1+q2+q3+theta+3.141592653589793*(1.0/4.0))*7.848E-1+(q1dot*q1dot)*sin(q2)*2.9575E-1-(q3dot*q3dot)*sin(q3)*(1.3E1/5.0E2)+(thetadot*thetadot)*sin(q2)*2.9575E-1+sqrt(4.1E1)*sin(q1+q2+q3+theta-8.960553845713439E-1)*1.962E-1+(q1dot*q1dot)*cos(q2+q3)*(1.3E1/5.0E2)+(thetadot*thetadot)*cos(q2+q3)*(1.3E1/5.0E2)+(q1dot*q1dot)*sin(q2+q3)*(1.3E1/5.0E2)+(thetadot*thetadot)*sin(q2+q3)*(1.3E1/5.0E2)-q1dot*q3dot*cos(q3)*(1.3E1/2.5E2)-q2dot*q3dot*cos(q3)*(1.3E1/2.5E2)-q3dot*thetadot*cos(q3)*(1.3E1/2.5E2)-q1dot*q3dot*sin(q3)*(1.3E1/2.5E2)-q2dot*q3dot*sin(q3)*(1.3E1/2.5E2)+q1dot*thetadot*sin(q2)*5.915E-1-q3dot*thetadot*sin(q3)*(1.3E1/2.5E2)+sqrt(4.1E1)*(q1dot*q1dot)*sin(q2+q3-8.960553845713439E-1)*6.5E-3+sqrt(4.1E1)*(thetadot*thetadot)*sin(q2+q3-8.960553845713439E-1)*6.5E-3+q1dot*thetadot*cos(q2+q3)*(1.3E1/2.5E2)+q1dot*thetadot*sin(q2+q3)*(1.3E1/2.5E2)-sqrt(4.1E1)*(q3dot*q3dot)*sin(q3-8.960553845713439E-1)*6.5E-3-sqrt(4.1E1)*q1dot*q3dot*sin(q3-8.960553845713439E-1)*(1.3E1/1.0E3)-sqrt(4.1E1)*q2dot*q3dot*sin(q3-8.960553845713439E-1)*(1.3E1/1.0E3)-sqrt(4.1E1)*q3dot*thetadot*sin(q3-8.960553845713439E-1)*(1.3E1/1.0E3)+sqrt(4.1E1)*q1dot*thetadot*sin(q2+q3-8.960553845713439E-1)*(1.3E1/1.0E3);
  T(5,0) = (q1dot*q1dot)*cos(q3)*(1.3E1/5.0E2)+(q2dot*q2dot)*cos(q3)*(1.3E1/5.0E2)+(thetadot*thetadot)*cos(q3)*(1.3E1/5.0E2)+sqrt(2.0)*sin(q1+q2+q3+theta+3.141592653589793*(1.0/4.0))*7.848E-1+(q1dot*q1dot)*sin(q3)*(1.3E1/5.0E2)+(q2dot*q2dot)*sin(q3)*(1.3E1/5.0E2)+(thetadot*thetadot)*sin(q3)*(1.3E1/5.0E2)+sqrt(4.1E1)*sin(q1+q2+q3+theta-8.960553845713439E-1)*1.962E-1+(q1dot*q1dot)*cos(q2+q3)*(1.3E1/5.0E2)+(thetadot*thetadot)*cos(q2+q3)*(1.3E1/5.0E2)+(q1dot*q1dot)*sin(q2+q3)*(1.3E1/5.0E2)+(thetadot*thetadot)*sin(q2+q3)*(1.3E1/5.0E2)+sqrt(4.1E1)*(thetadot*thetadot)*sin(q3-8.960553845713439E-1)*6.5E-3+q1dot*q2dot*cos(q3)*(1.3E1/2.5E2)+q1dot*thetadot*cos(q3)*(1.3E1/2.5E2)+q2dot*thetadot*cos(q3)*(1.3E1/2.5E2)+q1dot*q2dot*sin(q3)*(1.3E1/2.5E2)+q1dot*thetadot*sin(q3)*(1.3E1/2.5E2)+q2dot*thetadot*sin(q3)*(1.3E1/2.5E2)+sqrt(4.1E1)*(q1dot*q1dot)*sin(q2+q3-8.960553845713439E-1)*6.5E-3+sqrt(4.1E1)*(thetadot*thetadot)*sin(q2+q3-8.960553845713439E-1)*6.5E-3+q1dot*thetadot*cos(q2+q3)*(1.3E1/2.5E2)+q1dot*thetadot*sin(q2+q3)*(1.3E1/2.5E2)+sqrt(4.1E1)*(q1dot*q1dot)*sin(q3-8.960553845713439E-1)*6.5E-3+sqrt(4.1E1)*(q2dot*q2dot)*sin(q3-8.960553845713439E-1)*6.5E-3+sqrt(4.1E1)*q1dot*q2dot*sin(q3-8.960553845713439E-1)*(1.3E1/1.0E3)+sqrt(4.1E1)*q1dot*thetadot*sin(q3-8.960553845713439E-1)*(1.3E1/1.0E3)+sqrt(4.1E1)*q2dot*thetadot*sin(q3-8.960553845713439E-1)*(1.3E1/1.0E3)+sqrt(4.1E1)*q1dot*thetadot*sin(q2+q3-8.960553845713439E-1)*(1.3E1/1.0E3);
  T(6,0) = sin(q4+q5+theta)*8.9271+sin(q4+theta)*1.91295E1-(q6dot*q6dot)*cos(q6)*(1.3E1/5.0E2)+sqrt(2.0)*sin(q4+q5+q6+theta+3.141592653589793*(1.0/4.0))*7.848E-1-(q5dot*q5dot)*sin(q5)*2.9575E-1-(q6dot*q6dot)*sin(q6)*(1.3E1/5.0E2)+sqrt(4.1E1)*sin(q4+q5+q6+theta-8.960553845713439E-1)*1.962E-1-(q5dot*q5dot)*cos(q5+q6)*(1.3E1/5.0E2)-(q6dot*q6dot)*cos(q5+q6)*(1.3E1/5.0E2)-(q5dot*q5dot)*sin(q5+q6)*(1.3E1/5.0E2)-(q6dot*q6dot)*sin(q5+q6)*(1.3E1/5.0E2)-q4dot*q6dot*cos(q6)*(1.3E1/2.5E2)-q5dot*q6dot*cos(q6)*(1.3E1/2.5E2)-q6dot*thetadot*cos(q6)*(1.3E1/2.5E2)-q4dot*q5dot*sin(q5)*5.915E-1-q4dot*q6dot*sin(q6)*(1.3E1/2.5E2)-q5dot*q6dot*sin(q6)*(1.3E1/2.5E2)-q5dot*thetadot*sin(q5)*5.915E-1-q6dot*thetadot*sin(q6)*(1.3E1/2.5E2)-sqrt(4.1E1)*(q5dot*q5dot)*sin(q5+q6-8.960553845713439E-1)*6.5E-3-sqrt(4.1E1)*(q6dot*q6dot)*sin(q5+q6-8.960553845713439E-1)*6.5E-3-q4dot*q5dot*cos(q5+q6)*(1.3E1/2.5E2)-q4dot*q6dot*cos(q5+q6)*(1.3E1/2.5E2)-q5dot*q6dot*cos(q5+q6)*(1.3E1/2.5E2)-q5dot*thetadot*cos(q5+q6)*(1.3E1/2.5E2)-q6dot*thetadot*cos(q5+q6)*(1.3E1/2.5E2)-q4dot*q5dot*sin(q5+q6)*(1.3E1/2.5E2)-q4dot*q6dot*sin(q5+q6)*(1.3E1/2.5E2)-q5dot*q6dot*sin(q5+q6)*(1.3E1/2.5E2)-q5dot*thetadot*sin(q5+q6)*(1.3E1/2.5E2)-q6dot*thetadot*sin(q5+q6)*(1.3E1/2.5E2)-sqrt(4.1E1)*(q6dot*q6dot)*sin(q6-8.960553845713439E-1)*6.5E-3-sqrt(4.1E1)*q4dot*q6dot*sin(q6-8.960553845713439E-1)*(1.3E1/1.0E3)-sqrt(4.1E1)*q5dot*q6dot*sin(q6-8.960553845713439E-1)*(1.3E1/1.0E3)-sqrt(4.1E1)*q6dot*thetadot*sin(q6-8.960553845713439E-1)*(1.3E1/1.0E3)-sqrt(4.1E1)*q4dot*q5dot*sin(q5+q6-8.960553845713439E-1)*(1.3E1/1.0E3)-sqrt(4.1E1)*q4dot*q6dot*sin(q5+q6-8.960553845713439E-1)*(1.3E1/1.0E3)-sqrt(4.1E1)*q5dot*q6dot*sin(q5+q6-8.960553845713439E-1)*(1.3E1/1.0E3)-sqrt(4.1E1)*q5dot*thetadot*sin(q5+q6-8.960553845713439E-1)*(1.3E1/1.0E3)-sqrt(4.1E1)*q6dot*thetadot*sin(q5+q6-8.960553845713439E-1)*(1.3E1/1.0E3);
  T(7,0) = sin(q4+q5+theta)*8.9271-(q6dot*q6dot)*cos(q6)*(1.3E1/5.0E2)+sqrt(2.0)*sin(q4+q5+q6+theta+3.141592653589793*(1.0/4.0))*7.848E-1+(q4dot*q4dot)*sin(q5)*2.9575E-1-(q6dot*q6dot)*sin(q6)*(1.3E1/5.0E2)+(thetadot*thetadot)*sin(q5)*2.9575E-1+sqrt(4.1E1)*sin(q4+q5+q6+theta-8.960553845713439E-1)*1.962E-1+(q4dot*q4dot)*cos(q5+q6)*(1.3E1/5.0E2)+(thetadot*thetadot)*cos(q5+q6)*(1.3E1/5.0E2)+(q4dot*q4dot)*sin(q5+q6)*(1.3E1/5.0E2)+(thetadot*thetadot)*sin(q5+q6)*(1.3E1/5.0E2)-q4dot*q6dot*cos(q6)*(1.3E1/2.5E2)-q5dot*q6dot*cos(q6)*(1.3E1/2.5E2)-q6dot*thetadot*cos(q6)*(1.3E1/2.5E2)-q4dot*q6dot*sin(q6)*(1.3E1/2.5E2)-q5dot*q6dot*sin(q6)*(1.3E1/2.5E2)+q4dot*thetadot*sin(q5)*5.915E-1-q6dot*thetadot*sin(q6)*(1.3E1/2.5E2)+sqrt(4.1E1)*(q4dot*q4dot)*sin(q5+q6-8.960553845713439E-1)*6.5E-3+sqrt(4.1E1)*(thetadot*thetadot)*sin(q5+q6-8.960553845713439E-1)*6.5E-3+q4dot*thetadot*cos(q5+q6)*(1.3E1/2.5E2)+q4dot*thetadot*sin(q5+q6)*(1.3E1/2.5E2)-sqrt(4.1E1)*(q6dot*q6dot)*sin(q6-8.960553845713439E-1)*6.5E-3-sqrt(4.1E1)*q4dot*q6dot*sin(q6-8.960553845713439E-1)*(1.3E1/1.0E3)-sqrt(4.1E1)*q5dot*q6dot*sin(q6-8.960553845713439E-1)*(1.3E1/1.0E3)-sqrt(4.1E1)*q6dot*thetadot*sin(q6-8.960553845713439E-1)*(1.3E1/1.0E3)+sqrt(4.1E1)*q4dot*thetadot*sin(q5+q6-8.960553845713439E-1)*(1.3E1/1.0E3);
  T(8,0) = (q4dot*q4dot)*cos(q6)*(1.3E1/5.0E2)+(q5dot*q5dot)*cos(q6)*(1.3E1/5.0E2)+(thetadot*thetadot)*cos(q6)*(1.3E1/5.0E2)+sqrt(2.0)*sin(q4+q5+q6+theta+3.141592653589793*(1.0/4.0))*7.848E-1+(q4dot*q4dot)*sin(q6)*(1.3E1/5.0E2)+(q5dot*q5dot)*sin(q6)*(1.3E1/5.0E2)+(thetadot*thetadot)*sin(q6)*(1.3E1/5.0E2)+sqrt(4.1E1)*sin(q4+q5+q6+theta-8.960553845713439E-1)*1.962E-1+(q4dot*q4dot)*cos(q5+q6)*(1.3E1/5.0E2)+(thetadot*thetadot)*cos(q5+q6)*(1.3E1/5.0E2)+(q4dot*q4dot)*sin(q5+q6)*(1.3E1/5.0E2)+(thetadot*thetadot)*sin(q5+q6)*(1.3E1/5.0E2)+sqrt(4.1E1)*(thetadot*thetadot)*sin(q6-8.960553845713439E-1)*6.5E-3+q4dot*q5dot*cos(q6)*(1.3E1/2.5E2)+q4dot*thetadot*cos(q6)*(1.3E1/2.5E2)+q5dot*thetadot*cos(q6)*(1.3E1/2.5E2)+q4dot*q5dot*sin(q6)*(1.3E1/2.5E2)+q4dot*thetadot*sin(q6)*(1.3E1/2.5E2)+q5dot*thetadot*sin(q6)*(1.3E1/2.5E2)+sqrt(4.1E1)*(q4dot*q4dot)*sin(q5+q6-8.960553845713439E-1)*6.5E-3+sqrt(4.1E1)*(thetadot*thetadot)*sin(q5+q6-8.960553845713439E-1)*6.5E-3+q4dot*thetadot*cos(q5+q6)*(1.3E1/2.5E2)+q4dot*thetadot*sin(q5+q6)*(1.3E1/2.5E2)+sqrt(4.1E1)*(q4dot*q4dot)*sin(q6-8.960553845713439E-1)*6.5E-3+sqrt(4.1E1)*(q5dot*q5dot)*sin(q6-8.960553845713439E-1)*6.5E-3+sqrt(4.1E1)*q4dot*q5dot*sin(q6-8.960553845713439E-1)*(1.3E1/1.0E3)+sqrt(4.1E1)*q4dot*thetadot*sin(q6-8.960553845713439E-1)*(1.3E1/1.0E3)+sqrt(4.1E1)*q5dot*thetadot*sin(q6-8.960553845713439E-1)*(1.3E1/1.0E3)+sqrt(4.1E1)*q4dot*thetadot*sin(q5+q6-8.960553845713439E-1)*(1.3E1/1.0E3);
  T(9,0) = sin(q7+q8+theta)*6.0699375+sin(q7+theta)*1.1649375E1-(q8dot*q8dot)*sin(q8)*(9.9E1/6.4E2)-(thetadot*thetadot)*sin(q7)*(2.09E2/3.2E2)-(thetadot*thetadot)*sin(q7+q8)*3.403125E-1-q7dot*q8dot*sin(q8)*(9.9E1/3.2E2)-q8dot*thetadot*sin(q8)*(9.9E1/3.2E2);
  T(10,0) = sin(q7+q8+theta)*6.0699375+(q7dot*q7dot)*sin(q8)*(9.9E1/6.4E2)+(thetadot*thetadot)*sin(q8)*(9.9E1/6.4E2)-(thetadot*thetadot)*sin(q7+q8)*3.403125E-1+q7dot*thetadot*sin(q8)*(9.9E1/3.2E2);
  T(11,0) = sin(q9+q10+theta)*6.0699375+sin(q9+theta)*1.1649375E1-(q10dot*q10dot)*sin(q10)*(9.9E1/6.4E2)-(thetadot*thetadot)*sin(q9)*(2.09E2/3.2E2)-(thetadot*thetadot)*sin(q9+q10)*3.403125E-1-q10dot*q9dot*sin(q10)*(9.9E1/3.2E2)-q10dot*thetadot*sin(q10)*(9.9E1/3.2E2);
  T(12,0) = sin(q9+q10+theta)*6.0699375+(q9dot*q9dot)*sin(q10)*(9.9E1/6.4E2)+(thetadot*thetadot)*sin(q10)*(9.9E1/6.4E2)-(thetadot*thetadot)*sin(q9+q10)*3.403125E-1+q9dot*thetadot*sin(q10)*(9.9E1/3.2E2);

  return T;
}

dlib::matrix<double> Jac_Full_fn(const Robot_StateNDot &Robot_StateNDot_i)
{
  double rIx = Robot_StateNDot_i.rIx;
  double rIy = Robot_StateNDot_i.rIy;
  double theta = Robot_StateNDot_i.theta;
  double q1 = Robot_StateNDot_i.q1;
  double q2 = Robot_StateNDot_i.q2;
  double q3 = Robot_StateNDot_i.q3;
  double q4 = Robot_StateNDot_i.q4;
  double q5 = Robot_StateNDot_i.q5;
  double q6 = Robot_StateNDot_i.q6;
  double q7 = Robot_StateNDot_i.q7;
  double q8 = Robot_StateNDot_i.q8;
  double q9 = Robot_StateNDot_i.q9;
  double q10 = Robot_StateNDot_i.q10;

  dlib::matrix<double> T;
  T = dlib::zeros_matrix<double>(12,13);
  T(0,0) = 1.0;
  T(0,2) = sin(q1+theta+3.141592653589793*(1.0/2.0))*(-1.3E1/4.0E1)-sin(q1+q2+theta+3.141592653589793*(1.0/2.0))*(1.3E1/4.0E1)-sqrt(4.1E1)*sin(q1+q2+q3+theta+3.141592653589793*(1.0/2.0)-8.960553845713439E-1)*(1.0/4.0E1);
  T(0,3) = sin(q1+theta+3.141592653589793*(1.0/2.0))*(-1.3E1/4.0E1)-sin(q1+q2+theta+3.141592653589793*(1.0/2.0))*(1.3E1/4.0E1)-sqrt(4.1E1)*sin(q1+q2+q3+theta+3.141592653589793*(1.0/2.0)-8.960553845713439E-1)*(1.0/4.0E1);
  T(0,4) = sin(q1+q2+theta+3.141592653589793*(1.0/2.0))*(-1.3E1/4.0E1)-sqrt(4.1E1)*sin(q1+q2+q3+theta+3.141592653589793*(1.0/2.0)-8.960553845713439E-1)*(1.0/4.0E1);
  T(0,5) = sqrt(4.1E1)*sin(q1+q2+q3+theta+3.141592653589793*(1.0/2.0)-8.960553845713439E-1)*(-1.0/4.0E1);
  T(1,1) = 1.0;
  T(1,2) = cos(q1+theta+3.141592653589793*(1.0/2.0))*(-1.3E1/4.0E1)-cos(q1+q2+theta+3.141592653589793*(1.0/2.0))*(1.3E1/4.0E1)-sqrt(4.1E1)*cos(q1+q2+q3+theta+3.141592653589793*(1.0/2.0)-8.960553845713439E-1)*(1.0/4.0E1);
  T(1,3) = cos(q1+theta+3.141592653589793*(1.0/2.0))*(-1.3E1/4.0E1)-cos(q1+q2+theta+3.141592653589793*(1.0/2.0))*(1.3E1/4.0E1)-sqrt(4.1E1)*cos(q1+q2+q3+theta+3.141592653589793*(1.0/2.0)-8.960553845713439E-1)*(1.0/4.0E1);
  T(1,4) = cos(q1+q2+theta+3.141592653589793*(1.0/2.0))*(-1.3E1/4.0E1)-sqrt(4.1E1)*cos(q1+q2+q3+theta+3.141592653589793*(1.0/2.0)-8.960553845713439E-1)*(1.0/4.0E1);
  T(1,5) = sqrt(4.1E1)*cos(q1+q2+q3+theta+3.141592653589793*(1.0/2.0)-8.960553845713439E-1)*(-1.0/4.0E1);
  T(2,0) = 1.0;
  T(2,2) = sin(q1+theta+3.141592653589793*(1.0/2.0))*(-1.3E1/4.0E1)-sin(q1+q2+theta+3.141592653589793*(1.0/2.0))*(1.3E1/4.0E1)-sqrt(2.0)*sin(q1+q2+q3+theta-3.141592653589793*(5.0/4.0))*(1.0/1.0E1);
  T(2,3) = sin(q1+theta+3.141592653589793*(1.0/2.0))*(-1.3E1/4.0E1)-sin(q1+q2+theta+3.141592653589793*(1.0/2.0))*(1.3E1/4.0E1)-sqrt(2.0)*sin(q1+q2+q3+theta-3.141592653589793*(5.0/4.0))*(1.0/1.0E1);
  T(2,4) = sin(q1+q2+theta+3.141592653589793*(1.0/2.0))*(-1.3E1/4.0E1)-sqrt(2.0)*sin(q1+q2+q3+theta-3.141592653589793*(5.0/4.0))*(1.0/1.0E1);
  T(2,5) = sqrt(2.0)*sin(q1+q2+q3+theta-3.141592653589793*(5.0/4.0))*(-1.0/1.0E1);
  T(3,1) = 1.0;
  T(3,2) = cos(q1+theta+3.141592653589793*(1.0/2.0))*(-1.3E1/4.0E1)-cos(q1+q2+theta+3.141592653589793*(1.0/2.0))*(1.3E1/4.0E1)-sqrt(2.0)*cos(q1+q2+q3+theta-3.141592653589793*(5.0/4.0))*(1.0/1.0E1);
  T(3,3) = cos(q1+theta+3.141592653589793*(1.0/2.0))*(-1.3E1/4.0E1)-cos(q1+q2+theta+3.141592653589793*(1.0/2.0))*(1.3E1/4.0E1)-sqrt(2.0)*cos(q1+q2+q3+theta-3.141592653589793*(5.0/4.0))*(1.0/1.0E1);
  T(3,4) = cos(q1+q2+theta+3.141592653589793*(1.0/2.0))*(-1.3E1/4.0E1)-sqrt(2.0)*cos(q1+q2+q3+theta-3.141592653589793*(5.0/4.0))*(1.0/1.0E1);
  T(3,5) = sqrt(2.0)*cos(q1+q2+q3+theta-3.141592653589793*(5.0/4.0))*(-1.0/1.0E1);
  T(4,0) = 1.0;
  T(4,2) = sin(q4+theta+3.141592653589793*(1.0/2.0))*(-1.3E1/4.0E1)-sin(q4+q5+theta+3.141592653589793*(1.0/2.0))*(1.3E1/4.0E1)-sqrt(4.1E1)*sin(q4+q5+q6+theta+3.141592653589793*(1.0/2.0)-8.960553845713439E-1)*(1.0/4.0E1);
  T(4,6) = sin(q4+theta+3.141592653589793*(1.0/2.0))*(-1.3E1/4.0E1)-sin(q4+q5+theta+3.141592653589793*(1.0/2.0))*(1.3E1/4.0E1)-sqrt(4.1E1)*sin(q4+q5+q6+theta+3.141592653589793*(1.0/2.0)-8.960553845713439E-1)*(1.0/4.0E1);
  T(4,7) = sin(q4+q5+theta+3.141592653589793*(1.0/2.0))*(-1.3E1/4.0E1)-sqrt(4.1E1)*sin(q4+q5+q6+theta+3.141592653589793*(1.0/2.0)-8.960553845713439E-1)*(1.0/4.0E1);
  T(4,8) = sqrt(4.1E1)*sin(q4+q5+q6+theta+3.141592653589793*(1.0/2.0)-8.960553845713439E-1)*(-1.0/4.0E1);
  T(5,1) = 1.0;
  T(5,2) = cos(q4+theta+3.141592653589793*(1.0/2.0))*(-1.3E1/4.0E1)-cos(q4+q5+theta+3.141592653589793*(1.0/2.0))*(1.3E1/4.0E1)-sqrt(4.1E1)*cos(q4+q5+q6+theta+3.141592653589793*(1.0/2.0)-8.960553845713439E-1)*(1.0/4.0E1);
  T(5,6) = cos(q4+theta+3.141592653589793*(1.0/2.0))*(-1.3E1/4.0E1)-cos(q4+q5+theta+3.141592653589793*(1.0/2.0))*(1.3E1/4.0E1)-sqrt(4.1E1)*cos(q4+q5+q6+theta+3.141592653589793*(1.0/2.0)-8.960553845713439E-1)*(1.0/4.0E1);
  T(5,7) = cos(q4+q5+theta+3.141592653589793*(1.0/2.0))*(-1.3E1/4.0E1)-sqrt(4.1E1)*cos(q4+q5+q6+theta+3.141592653589793*(1.0/2.0)-8.960553845713439E-1)*(1.0/4.0E1);
  T(5,8) = sqrt(4.1E1)*cos(q4+q5+q6+theta+3.141592653589793*(1.0/2.0)-8.960553845713439E-1)*(-1.0/4.0E1);
  T(6,0) = 1.0;
  T(6,2) = sin(q4+theta+3.141592653589793*(1.0/2.0))*(-1.3E1/4.0E1)-sin(q4+q5+theta+3.141592653589793*(1.0/2.0))*(1.3E1/4.0E1)-sqrt(2.0)*sin(q4+q5+q6+theta+3.141592653589793*(3.0/4.0))*(1.0/1.0E1);
  T(6,6) = sin(q4+theta+3.141592653589793*(1.0/2.0))*(-1.3E1/4.0E1)-sin(q4+q5+theta+3.141592653589793*(1.0/2.0))*(1.3E1/4.0E1)-sqrt(2.0)*sin(q4+q5+q6+theta+3.141592653589793*(3.0/4.0))*(1.0/1.0E1);
  T(6,7) = sin(q4+q5+theta+3.141592653589793*(1.0/2.0))*(-1.3E1/4.0E1)-sqrt(2.0)*sin(q4+q5+q6+theta+3.141592653589793*(3.0/4.0))*(1.0/1.0E1);
  T(6,8) = sqrt(2.0)*sin(q4+q5+q6+theta+3.141592653589793*(3.0/4.0))*(-1.0/1.0E1);
  T(7,1) = 1.0;
  T(7,2) = cos(q4+theta+3.141592653589793*(1.0/2.0))*(-1.3E1/4.0E1)-cos(q4+q5+theta+3.141592653589793*(1.0/2.0))*(1.3E1/4.0E1)-sqrt(2.0)*cos(q4+q5+q6+theta+3.141592653589793*(3.0/4.0))*(1.0/1.0E1);
  T(7,6) = cos(q4+theta+3.141592653589793*(1.0/2.0))*(-1.3E1/4.0E1)-cos(q4+q5+theta+3.141592653589793*(1.0/2.0))*(1.3E1/4.0E1)-sqrt(2.0)*cos(q4+q5+q6+theta+3.141592653589793*(3.0/4.0))*(1.0/1.0E1);
  T(7,7) = cos(q4+q5+theta+3.141592653589793*(1.0/2.0))*(-1.3E1/4.0E1)-sqrt(2.0)*cos(q4+q5+q6+theta+3.141592653589793*(3.0/4.0))*(1.0/1.0E1);
  T(7,8) = sqrt(2.0)*cos(q4+q5+q6+theta+3.141592653589793*(3.0/4.0))*(-1.0/1.0E1);
  T(8,0) = 1.0;
  T(8,2) = sin(q7+theta+3.141592653589793*(1.0/2.0))*(-1.0/4.0)-sin(q7+q8+theta+3.141592653589793*(1.0/2.0))*(9.0/2.0E1)-sin(theta-3.141592653589793*(1.0/2.0))*(1.1E1/2.0E1);
  T(8,9) = sin(q7+theta+3.141592653589793*(1.0/2.0))*(-1.0/4.0)-sin(q7+q8+theta+3.141592653589793*(1.0/2.0))*(9.0/2.0E1);
  T(8,10) = sin(q7+q8+theta+3.141592653589793*(1.0/2.0))*(-9.0/2.0E1);
  T(9,1) = 1.0;
  T(9,2) = cos(q7+theta+3.141592653589793*(1.0/2.0))*(-1.0/4.0)-cos(q7+q8+theta+3.141592653589793*(1.0/2.0))*(9.0/2.0E1)-cos(theta-3.141592653589793*(1.0/2.0))*(1.1E1/2.0E1);
  T(9,9) = cos(q7+theta+3.141592653589793*(1.0/2.0))*(-1.0/4.0)-cos(q7+q8+theta+3.141592653589793*(1.0/2.0))*(9.0/2.0E1);
  T(9,10) = cos(q7+q8+theta+3.141592653589793*(1.0/2.0))*(-9.0/2.0E1);
  T(10,0) = 1.0;
  T(10,2) = sin(q9+theta+3.141592653589793*(1.0/2.0))*(-1.0/4.0)-sin(q9+q10+theta+3.141592653589793*(1.0/2.0))*(9.0/2.0E1)-sin(theta-3.141592653589793*(1.0/2.0))*(1.1E1/2.0E1);
  T(10,11) = sin(q9+theta+3.141592653589793*(1.0/2.0))*(-1.0/4.0)-sin(q9+q10+theta+3.141592653589793*(1.0/2.0))*(9.0/2.0E1);
  T(10,12) = sin(q9+q10+theta+3.141592653589793*(1.0/2.0))*(-9.0/2.0E1);
  T(11,1) = 1.0;
  T(11,2) = cos(q9+theta+3.141592653589793*(1.0/2.0))*(-1.0/4.0)-cos(q9+q10+theta+3.141592653589793*(1.0/2.0))*(9.0/2.0E1)-cos(theta-3.141592653589793*(1.0/2.0))*(1.1E1/2.0E1);
  T(11,11) = cos(q9+theta+3.141592653589793*(1.0/2.0))*(-1.0/4.0)-cos(q9+q10+theta+3.141592653589793*(1.0/2.0))*(9.0/2.0E1);
  T(11,12) = cos(q9+q10+theta+3.141592653589793*(1.0/2.0))*(-9.0/2.0E1);
  return T;
}

// dlib::matrix<double> Jacdot_qdot_fn(const Robot_StateNDot &Robot_StateNDot_i)
// {
//   double rIx = Robot_StateNDot_i.rIx;					double rIy = Robot_StateNDot_i.rIy;						double theta = Robot_StateNDot_i.theta;
//   double q1 = Robot_StateNDot_i.q1;					  double q2 = Robot_StateNDot_i.q2;						double q3 = Robot_StateNDot_i.q3;
//   double q4 = Robot_StateNDot_i.q4;					  double q5 = Robot_StateNDot_i.q5;						double q6 = Robot_StateNDot_i.q6;
//   double q7 = Robot_StateNDot_i.q7;					  double q8 = Robot_StateNDot_i.q8;						double q9 = Robot_StateNDot_i.q9;
//   double q10 = Robot_StateNDot_i.q10;
//
//   double rIxdot = Robot_StateNDot_i.rIxdot;			double rIydot = Robot_StateNDot_i.rIydot;				double thetadot = Robot_StateNDot_i.thetadot;
//   double q1dot = Robot_StateNDot_i.q1dot;				double q2dot = Robot_StateNDot_i.q2dot;					double q3dot = Robot_StateNDot_i.q3dot;
//   double q4dot = Robot_StateNDot_i.q4dot;				double q5dot = Robot_StateNDot_i.q5dot;					double q6dot = Robot_StateNDot_i.q6dot;
//   double q7dot = Robot_StateNDot_i.q7dot;				double q8dot = Robot_StateNDot_i.q8dot;					double q9dot = Robot_StateNDot_i.q9dot;
//   double q10dot = Robot_StateNDot_i.q10dot;
//
//   dlib::matrix<double> T;
//   T = dlib::zeros_matrix<double>(13,1);
//
//   T(0,0) = (q1dot*q1dot)*sin(q1+q2+theta)*(1.3E1/4.0E1)+(q2dot*q2dot)*sin(q1+q2+theta)*(1.3E1/4.0E1)+(thetadot*thetadot)*sin(q1+q2+theta)*(1.3E1/4.0E1)+(q1dot*q1dot)*sin(q1+theta)*(1.3E1/4.0E1)+(thetadot*thetadot)*sin(q1+theta)*(1.3E1/4.0E1)+sqrt(4.1E1)*(q1dot*q1dot)*sin(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/4.0E1)+sqrt(4.1E1)*(q2dot*q2dot)*sin(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/4.0E1)+sqrt(4.1E1)*(q3dot*q3dot)*sin(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/4.0E1)+sqrt(4.1E1)*(thetadot*thetadot)*sin(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/4.0E1)+q1dot*q2dot*sin(q1+q2+theta)*(1.3E1/2.0E1)+q1dot*thetadot*sin(q1+q2+theta)*(1.3E1/2.0E1)+q2dot*thetadot*sin(q1+q2+theta)*(1.3E1/2.0E1)+q1dot*thetadot*sin(q1+theta)*(1.3E1/2.0E1)+sqrt(4.1E1)*q1dot*q2dot*sin(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/2.0E1)+sqrt(4.1E1)*q1dot*q3dot*sin(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/2.0E1)+sqrt(4.1E1)*q2dot*q3dot*sin(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/2.0E1)+sqrt(4.1E1)*q1dot*thetadot*sin(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/2.0E1)+sqrt(4.1E1)*q2dot*thetadot*sin(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/2.0E1)+sqrt(4.1E1)*q3dot*thetadot*sin(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/2.0E1);
//   T(1,0) = (q1dot*q1dot)*cos(q1+q2+theta)*(1.3E1/4.0E1)+(q2dot*q2dot)*cos(q1+q2+theta)*(1.3E1/4.0E1)+(thetadot*thetadot)*cos(q1+q2+theta)*(1.3E1/4.0E1)+(q1dot*q1dot)*cos(q1+theta)*(1.3E1/4.0E1)+(thetadot*thetadot)*cos(q1+theta)*(1.3E1/4.0E1)+sqrt(4.1E1)*(q1dot*q1dot)*cos(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/4.0E1)+sqrt(4.1E1)*(q2dot*q2dot)*cos(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/4.0E1)+sqrt(4.1E1)*(q3dot*q3dot)*cos(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/4.0E1)+sqrt(4.1E1)*(thetadot*thetadot)*cos(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/4.0E1)+q1dot*q2dot*cos(q1+q2+theta)*(1.3E1/2.0E1)+q1dot*thetadot*cos(q1+q2+theta)*(1.3E1/2.0E1)+q2dot*thetadot*cos(q1+q2+theta)*(1.3E1/2.0E1)+q1dot*thetadot*cos(q1+theta)*(1.3E1/2.0E1)+sqrt(4.1E1)*q1dot*q2dot*cos(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/2.0E1)+sqrt(4.1E1)*q1dot*q3dot*cos(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/2.0E1)+sqrt(4.1E1)*q2dot*q3dot*cos(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/2.0E1)+sqrt(4.1E1)*q1dot*thetadot*cos(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/2.0E1)+sqrt(4.1E1)*q2dot*thetadot*cos(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/2.0E1)+sqrt(4.1E1)*q3dot*thetadot*cos(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/2.0E1);
//   T(2,0) = (q1dot*q1dot)*sin(q1+q2+theta)*(1.3E1/4.0E1)+(q2dot*q2dot)*sin(q1+q2+theta)*(1.3E1/4.0E1)+(thetadot*thetadot)*sin(q1+q2+theta)*(1.3E1/4.0E1)+(q1dot*q1dot)*cos(q1+q2+q3+theta)*(1.0/1.0E1)+(q2dot*q2dot)*cos(q1+q2+q3+theta)*(1.0/1.0E1)+(q3dot*q3dot)*cos(q1+q2+q3+theta)*(1.0/1.0E1)+(thetadot*thetadot)*cos(q1+q2+q3+theta)*(1.0/1.0E1)+(q1dot*q1dot)*sin(q1+q2+q3+theta)*(1.0/1.0E1)+(q2dot*q2dot)*sin(q1+q2+q3+theta)*(1.0/1.0E1)+(q3dot*q3dot)*sin(q1+q2+q3+theta)*(1.0/1.0E1)+(thetadot*thetadot)*sin(q1+q2+q3+theta)*(1.0/1.0E1)+(q1dot*q1dot)*sin(q1+theta)*(1.3E1/4.0E1)+(thetadot*thetadot)*sin(q1+theta)*(1.3E1/4.0E1)+q1dot*q2dot*sin(q1+q2+theta)*(1.3E1/2.0E1)+q1dot*thetadot*sin(q1+q2+theta)*(1.3E1/2.0E1)+q2dot*thetadot*sin(q1+q2+theta)*(1.3E1/2.0E1)+q1dot*q2dot*cos(q1+q2+q3+theta)*(1.0/5.0)+q1dot*q3dot*cos(q1+q2+q3+theta)*(1.0/5.0)+q2dot*q3dot*cos(q1+q2+q3+theta)*(1.0/5.0)+q1dot*thetadot*cos(q1+q2+q3+theta)*(1.0/5.0)+q2dot*thetadot*cos(q1+q2+q3+theta)*(1.0/5.0)+q3dot*thetadot*cos(q1+q2+q3+theta)*(1.0/5.0)+q1dot*q2dot*sin(q1+q2+q3+theta)*(1.0/5.0)+q1dot*q3dot*sin(q1+q2+q3+theta)*(1.0/5.0)+q2dot*q3dot*sin(q1+q2+q3+theta)*(1.0/5.0)+q1dot*thetadot*sin(q1+q2+q3+theta)*(1.0/5.0)+q2dot*thetadot*sin(q1+q2+q3+theta)*(1.0/5.0)+q3dot*thetadot*sin(q1+q2+q3+theta)*(1.0/5.0)+q1dot*thetadot*sin(q1+theta)*(1.3E1/2.0E1);
//   T(3,0) = (q1dot*q1dot)*cos(q1+q2+theta)*(1.3E1/4.0E1)+(q2dot*q2dot)*cos(q1+q2+theta)*(1.3E1/4.0E1)+(thetadot*thetadot)*cos(q1+q2+theta)*(1.3E1/4.0E1)+(q1dot*q1dot)*cos(q1+q2+q3+theta)*(1.0/1.0E1)+(q2dot*q2dot)*cos(q1+q2+q3+theta)*(1.0/1.0E1)+(q3dot*q3dot)*cos(q1+q2+q3+theta)*(1.0/1.0E1)+(thetadot*thetadot)*cos(q1+q2+q3+theta)*(1.0/1.0E1)-(q1dot*q1dot)*sin(q1+q2+q3+theta)*(1.0/1.0E1)-(q2dot*q2dot)*sin(q1+q2+q3+theta)*(1.0/1.0E1)-(q3dot*q3dot)*sin(q1+q2+q3+theta)*(1.0/1.0E1)-(thetadot*thetadot)*sin(q1+q2+q3+theta)*(1.0/1.0E1)+(q1dot*q1dot)*cos(q1+theta)*(1.3E1/4.0E1)+(thetadot*thetadot)*cos(q1+theta)*(1.3E1/4.0E1)+q1dot*q2dot*cos(q1+q2+theta)*(1.3E1/2.0E1)+q1dot*thetadot*cos(q1+q2+theta)*(1.3E1/2.0E1)+q2dot*thetadot*cos(q1+q2+theta)*(1.3E1/2.0E1)+q1dot*q2dot*cos(q1+q2+q3+theta)*(1.0/5.0)+q1dot*q3dot*cos(q1+q2+q3+theta)*(1.0/5.0)+q2dot*q3dot*cos(q1+q2+q3+theta)*(1.0/5.0)+q1dot*thetadot*cos(q1+q2+q3+theta)*(1.0/5.0)+q2dot*thetadot*cos(q1+q2+q3+theta)*(1.0/5.0)+q3dot*thetadot*cos(q1+q2+q3+theta)*(1.0/5.0)-q1dot*q2dot*sin(q1+q2+q3+theta)*(1.0/5.0)-q1dot*q3dot*sin(q1+q2+q3+theta)*(1.0/5.0)-q2dot*q3dot*sin(q1+q2+q3+theta)*(1.0/5.0)-q1dot*thetadot*sin(q1+q2+q3+theta)*(1.0/5.0)-q2dot*thetadot*sin(q1+q2+q3+theta)*(1.0/5.0)-q3dot*thetadot*sin(q1+q2+q3+theta)*(1.0/5.0)+q1dot*thetadot*cos(q1+theta)*(1.3E1/2.0E1);
//   T(4,0) = (q4dot*q4dot)*sin(q4+q5+theta)*(1.3E1/4.0E1)+(q5dot*q5dot)*sin(q4+q5+theta)*(1.3E1/4.0E1)+(thetadot*thetadot)*sin(q4+q5+theta)*(1.3E1/4.0E1)+(q4dot*q4dot)*sin(q4+theta)*(1.3E1/4.0E1)+(thetadot*thetadot)*sin(q4+theta)*(1.3E1/4.0E1)+sqrt(4.1E1)*(q4dot*q4dot)*sin(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/4.0E1)+sqrt(4.1E1)*(q5dot*q5dot)*sin(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/4.0E1)+sqrt(4.1E1)*(q6dot*q6dot)*sin(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/4.0E1)+sqrt(4.1E1)*(thetadot*thetadot)*sin(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/4.0E1)+q4dot*q5dot*sin(q4+q5+theta)*(1.3E1/2.0E1)+q4dot*thetadot*sin(q4+q5+theta)*(1.3E1/2.0E1)+q5dot*thetadot*sin(q4+q5+theta)*(1.3E1/2.0E1)+q4dot*thetadot*sin(q4+theta)*(1.3E1/2.0E1)+sqrt(4.1E1)*q4dot*q5dot*sin(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/2.0E1)+sqrt(4.1E1)*q4dot*q6dot*sin(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/2.0E1)+sqrt(4.1E1)*q5dot*q6dot*sin(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/2.0E1)+sqrt(4.1E1)*q4dot*thetadot*sin(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/2.0E1)+sqrt(4.1E1)*q5dot*thetadot*sin(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/2.0E1)+sqrt(4.1E1)*q6dot*thetadot*sin(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/2.0E1);
//   T(5,0) = (q4dot*q4dot)*cos(q4+q5+theta)*(1.3E1/4.0E1)+(q5dot*q5dot)*cos(q4+q5+theta)*(1.3E1/4.0E1)+(thetadot*thetadot)*cos(q4+q5+theta)*(1.3E1/4.0E1)+(q4dot*q4dot)*cos(q4+theta)*(1.3E1/4.0E1)+(thetadot*thetadot)*cos(q4+theta)*(1.3E1/4.0E1)+sqrt(4.1E1)*(q4dot*q4dot)*cos(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/4.0E1)+sqrt(4.1E1)*(q5dot*q5dot)*cos(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/4.0E1)+sqrt(4.1E1)*(q6dot*q6dot)*cos(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/4.0E1)+sqrt(4.1E1)*(thetadot*thetadot)*cos(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/4.0E1)+q4dot*q5dot*cos(q4+q5+theta)*(1.3E1/2.0E1)+q4dot*thetadot*cos(q4+q5+theta)*(1.3E1/2.0E1)+q5dot*thetadot*cos(q4+q5+theta)*(1.3E1/2.0E1)+q4dot*thetadot*cos(q4+theta)*(1.3E1/2.0E1)+sqrt(4.1E1)*q4dot*q5dot*cos(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/2.0E1)+sqrt(4.1E1)*q4dot*q6dot*cos(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/2.0E1)+sqrt(4.1E1)*q5dot*q6dot*cos(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/2.0E1)+sqrt(4.1E1)*q4dot*thetadot*cos(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/2.0E1)+sqrt(4.1E1)*q5dot*thetadot*cos(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/2.0E1)+sqrt(4.1E1)*q6dot*thetadot*cos(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/2.0E1);
//   T(6,0) = (q4dot*q4dot)*sin(q4+q5+theta)*(1.3E1/4.0E1)+(q5dot*q5dot)*sin(q4+q5+theta)*(1.3E1/4.0E1)+(thetadot*thetadot)*sin(q4+q5+theta)*(1.3E1/4.0E1)+(q4dot*q4dot)*cos(q4+q5+q6+theta)*(1.0/1.0E1)+(q5dot*q5dot)*cos(q4+q5+q6+theta)*(1.0/1.0E1)+(q6dot*q6dot)*cos(q4+q5+q6+theta)*(1.0/1.0E1)+(thetadot*thetadot)*cos(q4+q5+q6+theta)*(1.0/1.0E1)+(q4dot*q4dot)*sin(q4+q5+q6+theta)*(1.0/1.0E1)+(q5dot*q5dot)*sin(q4+q5+q6+theta)*(1.0/1.0E1)+(q6dot*q6dot)*sin(q4+q5+q6+theta)*(1.0/1.0E1)+(thetadot*thetadot)*sin(q4+q5+q6+theta)*(1.0/1.0E1)+(q4dot*q4dot)*sin(q4+theta)*(1.3E1/4.0E1)+(thetadot*thetadot)*sin(q4+theta)*(1.3E1/4.0E1)+q4dot*q5dot*sin(q4+q5+theta)*(1.3E1/2.0E1)+q4dot*thetadot*sin(q4+q5+theta)*(1.3E1/2.0E1)+q5dot*thetadot*sin(q4+q5+theta)*(1.3E1/2.0E1)+q4dot*q5dot*cos(q4+q5+q6+theta)*(1.0/5.0)+q4dot*q6dot*cos(q4+q5+q6+theta)*(1.0/5.0)+q5dot*q6dot*cos(q4+q5+q6+theta)*(1.0/5.0)+q4dot*thetadot*cos(q4+q5+q6+theta)*(1.0/5.0)+q5dot*thetadot*cos(q4+q5+q6+theta)*(1.0/5.0)+q6dot*thetadot*cos(q4+q5+q6+theta)*(1.0/5.0)+q4dot*q5dot*sin(q4+q5+q6+theta)*(1.0/5.0)+q4dot*q6dot*sin(q4+q5+q6+theta)*(1.0/5.0)+q5dot*q6dot*sin(q4+q5+q6+theta)*(1.0/5.0)+q4dot*thetadot*sin(q4+q5+q6+theta)*(1.0/5.0)+q5dot*thetadot*sin(q4+q5+q6+theta)*(1.0/5.0)+q6dot*thetadot*sin(q4+q5+q6+theta)*(1.0/5.0)+q4dot*thetadot*sin(q4+theta)*(1.3E1/2.0E1);
//   T(7,0) = (q4dot*q4dot)*cos(q4+q5+theta)*(1.3E1/4.0E1)+(q5dot*q5dot)*cos(q4+q5+theta)*(1.3E1/4.0E1)+(thetadot*thetadot)*cos(q4+q5+theta)*(1.3E1/4.0E1)+(q4dot*q4dot)*cos(q4+q5+q6+theta)*(1.0/1.0E1)+(q5dot*q5dot)*cos(q4+q5+q6+theta)*(1.0/1.0E1)+(q6dot*q6dot)*cos(q4+q5+q6+theta)*(1.0/1.0E1)+(thetadot*thetadot)*cos(q4+q5+q6+theta)*(1.0/1.0E1)-(q4dot*q4dot)*sin(q4+q5+q6+theta)*(1.0/1.0E1)-(q5dot*q5dot)*sin(q4+q5+q6+theta)*(1.0/1.0E1)-(q6dot*q6dot)*sin(q4+q5+q6+theta)*(1.0/1.0E1)-(thetadot*thetadot)*sin(q4+q5+q6+theta)*(1.0/1.0E1)+(q4dot*q4dot)*cos(q4+theta)*(1.3E1/4.0E1)+(thetadot*thetadot)*cos(q4+theta)*(1.3E1/4.0E1)+q4dot*q5dot*cos(q4+q5+theta)*(1.3E1/2.0E1)+q4dot*thetadot*cos(q4+q5+theta)*(1.3E1/2.0E1)+q5dot*thetadot*cos(q4+q5+theta)*(1.3E1/2.0E1)+q4dot*q5dot*cos(q4+q5+q6+theta)*(1.0/5.0)+q4dot*q6dot*cos(q4+q5+q6+theta)*(1.0/5.0)+q5dot*q6dot*cos(q4+q5+q6+theta)*(1.0/5.0)+q4dot*thetadot*cos(q4+q5+q6+theta)*(1.0/5.0)+q5dot*thetadot*cos(q4+q5+q6+theta)*(1.0/5.0)+q6dot*thetadot*cos(q4+q5+q6+theta)*(1.0/5.0)-q4dot*q5dot*sin(q4+q5+q6+theta)*(1.0/5.0)-q4dot*q6dot*sin(q4+q5+q6+theta)*(1.0/5.0)-q5dot*q6dot*sin(q4+q5+q6+theta)*(1.0/5.0)-q4dot*thetadot*sin(q4+q5+q6+theta)*(1.0/5.0)-q5dot*thetadot*sin(q4+q5+q6+theta)*(1.0/5.0)-q6dot*thetadot*sin(q4+q5+q6+theta)*(1.0/5.0)+q4dot*thetadot*cos(q4+theta)*(1.3E1/2.0E1);
//   T(8,0) = (thetadot*thetadot)*sin(theta)*(-1.1E1/2.0E1)+(q7dot*q7dot)*sin(q7+q8+theta)*(9.0/2.0E1)+(q8dot*q8dot)*sin(q7+q8+theta)*(9.0/2.0E1)+(thetadot*thetadot)*sin(q7+q8+theta)*(9.0/2.0E1)+(q7dot*q7dot)*sin(q7+theta)*(1.0/4.0)+(thetadot*thetadot)*sin(q7+theta)*(1.0/4.0)+q7dot*q8dot*sin(q7+q8+theta)*(9.0/1.0E1)+q7dot*thetadot*sin(q7+q8+theta)*(9.0/1.0E1)+q8dot*thetadot*sin(q7+q8+theta)*(9.0/1.0E1)+q7dot*thetadot*sin(q7+theta)*(1.0/2.0);
//   T(9,0) = (thetadot*thetadot)*cos(theta)*(-1.1E1/2.0E1)+(q7dot*q7dot)*cos(q7+q8+theta)*(9.0/2.0E1)+(q8dot*q8dot)*cos(q7+q8+theta)*(9.0/2.0E1)+(thetadot*thetadot)*cos(q7+q8+theta)*(9.0/2.0E1)+(q7dot*q7dot)*cos(q7+theta)*(1.0/4.0)+(thetadot*thetadot)*cos(q7+theta)*(1.0/4.0)+q7dot*q8dot*cos(q7+q8+theta)*(9.0/1.0E1)+q7dot*thetadot*cos(q7+q8+theta)*(9.0/1.0E1)+q8dot*thetadot*cos(q7+q8+theta)*(9.0/1.0E1)+q7dot*thetadot*cos(q7+theta)*(1.0/2.0);
//   T(10,0) = (thetadot*thetadot)*sin(theta)*(-1.1E1/2.0E1)+(q10dot*q10dot)*sin(q9+q10+theta)*(9.0/2.0E1)+(q9dot*q9dot)*sin(q9+q10+theta)*(9.0/2.0E1)+(thetadot*thetadot)*sin(q9+q10+theta)*(9.0/2.0E1)+(q9dot*q9dot)*sin(q9+theta)*(1.0/4.0)+(thetadot*thetadot)*sin(q9+theta)*(1.0/4.0)+q10dot*q9dot*sin(q9+q10+theta)*(9.0/1.0E1)+q10dot*thetadot*sin(q9+q10+theta)*(9.0/1.0E1)+q9dot*thetadot*sin(q9+q10+theta)*(9.0/1.0E1)+q9dot*thetadot*sin(q9+theta)*(1.0/2.0);
//   T(11,0) = (thetadot*thetadot)*cos(theta)*(-1.1E1/2.0E1)+(q10dot*q10dot)*cos(q9+q10+theta)*(9.0/2.0E1)+(q9dot*q9dot)*cos(q9+q10+theta)*(9.0/2.0E1)+(thetadot*thetadot)*cos(q9+q10+theta)*(9.0/2.0E1)+(q9dot*q9dot)*cos(q9+theta)*(1.0/4.0)+(thetadot*thetadot)*cos(q9+theta)*(1.0/4.0)+q10dot*q9dot*cos(q9+q10+theta)*(9.0/1.0E1)+q10dot*thetadot*cos(q9+q10+theta)*(9.0/1.0E1)+q9dot*thetadot*cos(q9+q10+theta)*(9.0/1.0E1)+q9dot*thetadot*cos(q9+theta)*(1.0/2.0);
//
//   return T;
// }

std::vector<double> Ang_Pos_fn(const Robot_StateNDot &Robot_StateNDot_i, const char* s)
{

  std::vector<double> r_pos(2);
  r_pos[0] = 0.0;
  r_pos[1] = 0.0;

  double rIx = Robot_StateNDot_i.rIx;							double rIy = Robot_StateNDot_i.rIy;							double theta = Robot_StateNDot_i.theta;
  double q1 = Robot_StateNDot_i.q1;							double q2 = Robot_StateNDot_i.q2;							double q3 = Robot_StateNDot_i.q3;
  double q4 = Robot_StateNDot_i.q4;							double q5 = Robot_StateNDot_i.q5;							double q6 = Robot_StateNDot_i.q6;
  double q7 = Robot_StateNDot_i.q7;							double q8 = Robot_StateNDot_i.q8;							double q9 = Robot_StateNDot_i.q9;
  double q10 = Robot_StateNDot_i.q10;

  if (strcmp(s,"rA")==0)
  {
    r_pos[0] = rIx+cos(q1+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+cos(q1+q2+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+sqrt(4.1E1)*cos(q1+q2+q3+theta+PI*(1.0/2.0)-8.960553845713439E-1)*(1.0/4.0E1);
    r_pos[1] = rIy-sin(q1+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)-sin(q1+q2+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)-sqrt(4.1E1)*sin(q1+q2+q3+theta+PI*(1.0/2.0)-8.960553845713439E-1)*(1.0/4.0E1);
  }
  if (strcmp(s,"rB")==0)
  {
    r_pos[0] = rIx+cos(q1+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+cos(q1+q2+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+sqrt(2.0)*cos(q1+q2+q3+theta-PI*(5.0/4.0))*(1.0/1.0E1);
    r_pos[1] = rIy-sin(q1+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)-sin(q1+q2+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)-sqrt(2.0)*sin(q1+q2+q3+theta-PI*(5.0/4.0))*(1.0/1.0E1);
  }
  if (strcmp(s,"rC")==0)
  {
    r_pos[0] = rIx+cos(q4+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+cos(q4+q5+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+sqrt(4.1E1)*cos(q4+q5+q6+theta+PI*(1.0/2.0)-8.960553845713439E-1)*(1.0/4.0E1);
    r_pos[1] = rIy-sin(q4+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)-sin(q4+q5+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)-sqrt(4.1E1)*sin(q4+q5+q6+theta+PI*(1.0/2.0)-8.960553845713439E-1)*(1.0/4.0E1);
  }
  if (strcmp(s,"rD")==0)
  {
    r_pos[0] = rIx+cos(q4+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+cos(q4+q5+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+sqrt(2.0)*cos(q4+q5+q6+theta+PI*(3.0/4.0))*(1.0/1.0E1);
    r_pos[1] = rIy-sin(q4+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)-sin(q4+q5+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)-sqrt(2.0)*sin(q4+q5+q6+theta+PI*(3.0/4.0))*(1.0/1.0E1);
  }
  if (strcmp(s,"rE")==0)
  {
    r_pos[0] = rIx+cos(q7+theta+PI*(1.0/2.0))*(1.0/4.0)+cos(q7+q8+theta+PI*(1.0/2.0))*(9.0/2.0E1)+cos(theta-PI*(1.0/2.0))*(1.1E1/2.0E1);
    r_pos[1] = rIy-sin(q7+theta+PI*(1.0/2.0))*(1.0/4.0)-sin(q7+q8+theta+PI*(1.0/2.0))*(9.0/2.0E1)-sin(theta-PI*(1.0/2.0))*(1.1E1/2.0E1);
  }
  if (strcmp(s,"rF")==0)
  {
    r_pos[0] = rIx+cos(q9+theta+PI*(1.0/2.0))*(1.0/4.0)+cos(q9+q10+theta+PI*(1.0/2.0))*(9.0/2.0E1)+cos(theta-PI*(1.0/2.0))*(1.1E1/2.0E1);
    r_pos[1] = rIy-sin(q9+theta+PI*(1.0/2.0))*(1.0/4.0)-sin(q9+q10+theta+PI*(1.0/2.0))*(9.0/2.0E1)-sin(theta-PI*(1.0/2.0))*(1.1E1/2.0E1);
  }
  if (strcmp(s,"rG")==0)
  {
    r_pos[0] = rIx+cos(q1+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+cos(q1+q2+theta+PI*(1.0/2.0))*(1.3E1/4.0E1);
    r_pos[1] = rIy-sin(q1+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)-sin(q1+q2+theta+PI*(1.0/2.0))*(1.3E1/4.0E1);
  }
  if (strcmp(s,"rH")==0)
  {
    r_pos[0] = rIx+cos(q1+theta+PI*(1.0/2.0))*(1.3E1/4.0E1);
    r_pos[1] = rIy-sin(q1+theta+PI*(1.0/2.0))*(1.3E1/4.0E1);
  }
  if (strcmp(s,"rI")==0)
  {
    r_pos[0] = rIx;
    r_pos[1] = rIy;
  }
  if (strcmp(s,"rJ")==0)
  {
    r_pos[0] = rIx+cos(q4+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+cos(q4+q5+theta+PI*(1.0/2.0))*(1.3E1/4.0E1);
    r_pos[1] = rIy-sin(q4+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)-sin(q4+q5+theta+PI*(1.0/2.0))*(1.3E1/4.0E1);
  }
  if (strcmp(s,"rK")==0)
  {
    r_pos[0] = rIx+cos(q4+theta+PI*(1.0/2.0))*(1.3E1/4.0E1);
    r_pos[1] = rIy-sin(q4+theta+PI*(1.0/2.0))*(1.3E1/4.0E1);
  }
  if (strcmp(s,"rL")==0)
  {
    r_pos[0] = rIx+cos(theta-PI*(1.0/2.0))*(1.1E1/2.0E1);
    r_pos[1] = rIy-sin(theta-PI*(1.0/2.0))*(1.1E1/2.0E1);
  }
  if (strcmp(s,"rM")==0)
  {
    r_pos[0] = rIx+cos(q7+theta+PI*(1.0/2.0))*(1.0/4.0)+cos(theta-PI*(1.0/2.0))*(1.1E1/2.0E1);
    r_pos[1] = rIy-sin(q7+theta+PI*(1.0/2.0))*(1.0/4.0)-sin(theta-PI*(1.0/2.0))*(1.1E1/2.0E1);
  }
  if (strcmp(s,"rN")==0)
  {
    r_pos[0] = rIx+cos(q9+theta+PI*(1.0/2.0))*(1.0/4.0)+cos(theta-PI*(1.0/2.0))*(1.1E1/2.0E1);
    r_pos[1] = rIy-sin(q9+theta+PI*(1.0/2.0))*(1.0/4.0)-sin(theta-PI*(1.0/2.0))*(1.1E1/2.0E1);
  }
  if (strcmp(s,"rT")==0)
  {
    r_pos[0] = rIx+cos(theta-PI*(1.0/2.0))*(7.0/1.0E1);
    r_pos[1] = rIy-sin(theta-PI*(1.0/2.0))*(7.0/1.0E1);
  }
  if (strcmp(s,"rCOM")==0)
  {
    r_pos[0] = rIx-sin(q1+q2+theta)*1.678966789667897E-2-sin(q4+q5+theta)*1.678966789667897E-2-sin(q7+q8+theta)*8.717712177121771E-3-sin(q9+q10+theta)*8.717712177121771E-3-sin(q1+theta)*3.597785977859779E-2-sin(q4+theta)*3.597785977859779E-2-sin(q7+theta)*1.775830258302583E-2-sin(q9+theta)*1.775830258302583E-2+sin(theta)*2.506457564575646E-1-sqrt(2.0)*sin(q1+q2+q3+theta+PI*(1.0/4.0))*1.476014760147601E-3-sqrt(2.0)*sin(q4+q5+q6+theta+PI*(1.0/4.0))*1.476014760147601E-3-sqrt(4.1E1)*sin(q1+q2+q3+theta-8.960553845713439E-1)*3.690036900369004E-4-sqrt(4.1E1)*sin(q4+q5+q6+theta-8.960553845713439E-1)*3.690036900369004E-4;
    r_pos[1] = rIy-cos(q1+q2+theta)*1.678966789667897E-2-cos(q4+q5+theta)*1.678966789667897E-2-cos(q7+q8+theta)*8.717712177121771E-3-cos(q9+q10+theta)*8.717712177121771E-3-cos(q1+theta)*3.597785977859779E-2-cos(q4+theta)*3.597785977859779E-2-cos(q7+theta)*1.775830258302583E-2-cos(q9+theta)*1.775830258302583E-2+cos(theta)*2.506457564575646E-1-sqrt(4.1E1)*cos(q1+q2+q3+theta-8.960553845713439E-1)*3.690036900369004E-4-sqrt(4.1E1)*cos(q4+q5+q6+theta-8.960553845713439E-1)*3.690036900369004E-4-sqrt(2.0)*cos(q1+q2+q3+theta+PI*(1.0/4.0))*1.476014760147601E-3-sqrt(2.0)*cos(q4+q5+q6+theta+PI*(1.0/4.0))*1.476014760147601E-3;
  }
  return r_pos;
}

std::vector<double> Ang_Vel_fn(const Robot_StateNDot &Robot_StateNDot_i, const char* s)
{
  double PI = 3.1415926535897932384626;
  std::vector<double> T(2);

  T[0] = 0.0;
  T[1] = 0.0;

  // double rIx = Robot_StateNDot_i.rIx;
  // double rIy = Robot_StateNDot_i.rIy;
  double theta = Robot_StateNDot_i.theta;					double q1 = Robot_StateNDot_i.q1;				double q2 = Robot_StateNDot_i.q2;
  double q3 = Robot_StateNDot_i.q3;						double q4 = Robot_StateNDot_i.q4;				double q5 = Robot_StateNDot_i.q5;
  double q6 = Robot_StateNDot_i.q6;						double q7 = Robot_StateNDot_i.q7;				double q8 = Robot_StateNDot_i.q8;
  double q9 = Robot_StateNDot_i.q9;						double q10 = Robot_StateNDot_i.q10;

  double rIxdot = Robot_StateNDot_i.rIxdot;				double rIydot = Robot_StateNDot_i.rIydot;		double thetadot = Robot_StateNDot_i.thetadot;
  double q1dot = Robot_StateNDot_i.q1dot;					double q2dot = Robot_StateNDot_i.q2dot;			double q3dot = Robot_StateNDot_i.q3dot;
  double q4dot = Robot_StateNDot_i.q4dot;					double q5dot = Robot_StateNDot_i.q5dot;			double q6dot = Robot_StateNDot_i.q6dot;
  double q7dot = Robot_StateNDot_i.q7dot;					double q8dot = Robot_StateNDot_i.q8dot;			double q9dot = Robot_StateNDot_i.q9dot;
  double q10dot = Robot_StateNDot_i.q10dot;

  if (strcmp(s,"vA")==0)
  {
    T[0] = rIxdot-q1dot*(sin(q1+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+sin(q1+q2+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+sqrt(4.1E1)*sin(q1+q2+q3+theta+PI*(1.0/2.0)-8.960553845713439E-1)*(1.0/4.0E1))-thetadot*(sin(q1+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+sin(q1+q2+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+sqrt(4.1E1)*sin(q1+q2+q3+theta+PI*(1.0/2.0)-8.960553845713439E-1)*(1.0/4.0E1))-q2dot*(sin(q1+q2+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+sqrt(4.1E1)*sin(q1+q2+q3+theta+PI*(1.0/2.0)-8.960553845713439E-1)*(1.0/4.0E1))-sqrt(4.1E1)*q3dot*sin(q1+q2+q3+theta+PI*(1.0/2.0)-8.960553845713439E-1)*(1.0/4.0E1);
    T[1] = rIydot-q1dot*(cos(q1+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+cos(q1+q2+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+sqrt(4.1E1)*cos(q1+q2+q3+theta+PI*(1.0/2.0)-8.960553845713439E-1)*(1.0/4.0E1))-thetadot*(cos(q1+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+cos(q1+q2+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+sqrt(4.1E1)*cos(q1+q2+q3+theta+PI*(1.0/2.0)-8.960553845713439E-1)*(1.0/4.0E1))-q2dot*(cos(q1+q2+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+sqrt(4.1E1)*cos(q1+q2+q3+theta+PI*(1.0/2.0)-8.960553845713439E-1)*(1.0/4.0E1))-sqrt(4.1E1)*q3dot*cos(q1+q2+q3+theta+PI*(1.0/2.0)-8.960553845713439E-1)*(1.0/4.0E1);
  }
  if (strcmp(s,"vB")==0)
  {
    T[0] = rIxdot-q1dot*(sin(q1+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+sin(q1+q2+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+sqrt(2.0)*sin(q1+q2+q3+theta-PI*(5.0/4.0))*(1.0/1.0E1))-thetadot*(sin(q1+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+sin(q1+q2+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+sqrt(2.0)*sin(q1+q2+q3+theta-PI*(5.0/4.0))*(1.0/1.0E1))-q2dot*(sin(q1+q2+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+sqrt(2.0)*sin(q1+q2+q3+theta-PI*(5.0/4.0))*(1.0/1.0E1))-sqrt(2.0)*q3dot*sin(q1+q2+q3+theta-PI*(5.0/4.0))*(1.0/1.0E1);
    T[1] = rIydot-q1dot*(cos(q1+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+cos(q1+q2+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+sqrt(2.0)*cos(q1+q2+q3+theta-PI*(5.0/4.0))*(1.0/1.0E1))-thetadot*(cos(q1+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+cos(q1+q2+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+sqrt(2.0)*cos(q1+q2+q3+theta-PI*(5.0/4.0))*(1.0/1.0E1))-q2dot*(cos(q1+q2+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+sqrt(2.0)*cos(q1+q2+q3+theta-PI*(5.0/4.0))*(1.0/1.0E1))-sqrt(2.0)*q3dot*cos(q1+q2+q3+theta-PI*(5.0/4.0))*(1.0/1.0E1);
  }
  if (strcmp(s,"vC")==0)
  {
    T[0] = rIxdot-q4dot*(sin(q4+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+sin(q4+q5+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+sqrt(4.1E1)*sin(q4+q5+q6+theta+PI*(1.0/2.0)-8.960553845713439E-1)*(1.0/4.0E1))-thetadot*(sin(q4+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+sin(q4+q5+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+sqrt(4.1E1)*sin(q4+q5+q6+theta+PI*(1.0/2.0)-8.960553845713439E-1)*(1.0/4.0E1))-q5dot*(sin(q4+q5+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+sqrt(4.1E1)*sin(q4+q5+q6+theta+PI*(1.0/2.0)-8.960553845713439E-1)*(1.0/4.0E1))-sqrt(4.1E1)*q6dot*sin(q4+q5+q6+theta+PI*(1.0/2.0)-8.960553845713439E-1)*(1.0/4.0E1);
    T[1] = rIydot-q4dot*(cos(q4+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+cos(q4+q5+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+sqrt(4.1E1)*cos(q4+q5+q6+theta+PI*(1.0/2.0)-8.960553845713439E-1)*(1.0/4.0E1))-thetadot*(cos(q4+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+cos(q4+q5+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+sqrt(4.1E1)*cos(q4+q5+q6+theta+PI*(1.0/2.0)-8.960553845713439E-1)*(1.0/4.0E1))-q5dot*(cos(q4+q5+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+sqrt(4.1E1)*cos(q4+q5+q6+theta+PI*(1.0/2.0)-8.960553845713439E-1)*(1.0/4.0E1))-sqrt(4.1E1)*q6dot*cos(q4+q5+q6+theta+PI*(1.0/2.0)-8.960553845713439E-1)*(1.0/4.0E1);
  }
  if (strcmp(s,"vD")==0)
  {
    T[0] = rIxdot-q4dot*(sin(q4+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+sin(q4+q5+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+sqrt(2.0)*sin(q4+q5+q6+theta+PI*(3.0/4.0))*(1.0/1.0E1))-thetadot*(sin(q4+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+sin(q4+q5+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+sqrt(2.0)*sin(q4+q5+q6+theta+PI*(3.0/4.0))*(1.0/1.0E1))-q5dot*(sin(q4+q5+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+sqrt(2.0)*sin(q4+q5+q6+theta+PI*(3.0/4.0))*(1.0/1.0E1))-sqrt(2.0)*q6dot*sin(q4+q5+q6+theta+PI*(3.0/4.0))*(1.0/1.0E1);
    T[1] = rIydot-q4dot*(cos(q4+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+cos(q4+q5+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+sqrt(2.0)*cos(q4+q5+q6+theta+PI*(3.0/4.0))*(1.0/1.0E1))-thetadot*(cos(q4+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+cos(q4+q5+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+sqrt(2.0)*cos(q4+q5+q6+theta+PI*(3.0/4.0))*(1.0/1.0E1))-q5dot*(cos(q4+q5+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+sqrt(2.0)*cos(q4+q5+q6+theta+PI*(3.0/4.0))*(1.0/1.0E1))-sqrt(2.0)*q6dot*cos(q4+q5+q6+theta+PI*(3.0/4.0))*(1.0/1.0E1);
  }
  if (strcmp(s,"vE")==0)
  {
    T[0] = rIxdot-q7dot*(sin(q7+theta+PI*(1.0/2.0))*(1.0/4.0)+sin(q7+q8+theta+PI*(1.0/2.0))*(9.0/2.0E1))-thetadot*(sin(q7+theta+PI*(1.0/2.0))*(1.0/4.0)+sin(q7+q8+theta+PI*(1.0/2.0))*(9.0/2.0E1)+sin(theta-PI*(1.0/2.0))*(1.1E1/2.0E1))-q8dot*sin(q7+q8+theta+PI*(1.0/2.0))*(9.0/2.0E1);
    T[1] = rIydot-q7dot*(cos(q7+theta+PI*(1.0/2.0))*(1.0/4.0)+cos(q7+q8+theta+PI*(1.0/2.0))*(9.0/2.0E1))-thetadot*(cos(q7+theta+PI*(1.0/2.0))*(1.0/4.0)+cos(q7+q8+theta+PI*(1.0/2.0))*(9.0/2.0E1)+cos(theta-PI*(1.0/2.0))*(1.1E1/2.0E1))-q8dot*cos(q7+q8+theta+PI*(1.0/2.0))*(9.0/2.0E1);
  }
  if(strcmp(s,"vF")==0)
  {
    T[0] = rIxdot-q9dot*(sin(q9+theta+PI*(1.0/2.0))*(1.0/4.0)+sin(q9+q10+theta+PI*(1.0/2.0))*(9.0/2.0E1))-thetadot*(sin(q9+theta+PI*(1.0/2.0))*(1.0/4.0)+sin(q9+q10+theta+PI*(1.0/2.0))*(9.0/2.0E1)+sin(theta-PI*(1.0/2.0))*(1.1E1/2.0E1))-q10dot*sin(q9+q10+theta+PI*(1.0/2.0))*(9.0/2.0E1);
    T[1] = rIydot-q9dot*(cos(q9+theta+PI*(1.0/2.0))*(1.0/4.0)+cos(q9+q10+theta+PI*(1.0/2.0))*(9.0/2.0E1))-thetadot*(cos(q9+theta+PI*(1.0/2.0))*(1.0/4.0)+cos(q9+q10+theta+PI*(1.0/2.0))*(9.0/2.0E1)+cos(theta-PI*(1.0/2.0))*(1.1E1/2.0E1))-q10dot*cos(q9+q10+theta+PI*(1.0/2.0))*(9.0/2.0E1);
  }
  if(strcmp(s,"vI")==0)
  {
    T[0] = rIxdot;
    T[1] = rIydot;
  }
  if(strcmp(s,"vM")==0)
  {
    T[0] = rIxdot-q7dot*sin(q7+theta+3.141592653589793*(1.0/2.0))*(1.0/4.0)-thetadot*(sin(q7+theta+3.141592653589793*(1.0/2.0))*(1.0/4.0)+sin(theta-3.141592653589793*(1.0/2.0))*(1.1E1/2.0E1));
    T[1] = rIydot-q7dot*cos(q7+theta+3.141592653589793*(1.0/2.0))*(1.0/4.0)-thetadot*(cos(q7+theta+3.141592653589793*(1.0/2.0))*(1.0/4.0)+cos(theta-3.141592653589793*(1.0/2.0))*(1.1E1/2.0E1));
  }
  if(strcmp(s,"vN")==0)
  {
    T[0] = rIxdot-q9dot*sin(q9+theta+3.141592653589793*(1.0/2.0))*(1.0/4.0)-thetadot*(sin(q9+theta+3.141592653589793*(1.0/2.0))*(1.0/4.0)+sin(theta-3.141592653589793*(1.0/2.0))*(1.1E1/2.0E1));
    T[1] = rIydot-q9dot*cos(q9+theta+3.141592653589793*(1.0/2.0))*(1.0/4.0)-thetadot*(cos(q9+theta+3.141592653589793*(1.0/2.0))*(1.0/4.0)+cos(theta-3.141592653589793*(1.0/2.0))*(1.1E1/2.0E1));
  }
  if(strcmp(s,"vG")==0)
  {
    T[0] = rIxdot-q1dot*(sin(q1+theta+3.141592653589793*(1.0/2.0))*(1.3E1/4.0E1)+sin(q1+q2+theta+3.141592653589793*(1.0/2.0))*(1.3E1/4.0E1))-thetadot*(sin(q1+theta+3.141592653589793*(1.0/2.0))*(1.3E1/4.0E1)+sin(q1+q2+theta+3.141592653589793*(1.0/2.0))*(1.3E1/4.0E1))-q2dot*sin(q1+q2+theta+3.141592653589793*(1.0/2.0))*(1.3E1/4.0E1);
    T[1] = rIydot-q1dot*(cos(q1+theta+3.141592653589793*(1.0/2.0))*(1.3E1/4.0E1)+cos(q1+q2+theta+3.141592653589793*(1.0/2.0))*(1.3E1/4.0E1))-thetadot*(cos(q1+theta+3.141592653589793*(1.0/2.0))*(1.3E1/4.0E1)+cos(q1+q2+theta+3.141592653589793*(1.0/2.0))*(1.3E1/4.0E1))-q2dot*cos(q1+q2+theta+3.141592653589793*(1.0/2.0))*(1.3E1/4.0E1);
  }
  if(strcmp(s,"vJ")==0)
  {
    T[0] = rIxdot-q4dot*(sin(q4+theta+3.141592653589793*(1.0/2.0))*(1.3E1/4.0E1)+sin(q4+q5+theta+3.141592653589793*(1.0/2.0))*(1.3E1/4.0E1))-thetadot*(sin(q4+theta+3.141592653589793*(1.0/2.0))*(1.3E1/4.0E1)+sin(q4+q5+theta+3.141592653589793*(1.0/2.0))*(1.3E1/4.0E1))-q5dot*sin(q4+q5+theta+3.141592653589793*(1.0/2.0))*(1.3E1/4.0E1);
    T[1] = rIydot-q4dot*(cos(q4+theta+3.141592653589793*(1.0/2.0))*(1.3E1/4.0E1)+cos(q4+q5+theta+3.141592653589793*(1.0/2.0))*(1.3E1/4.0E1))-thetadot*(cos(q4+theta+3.141592653589793*(1.0/2.0))*(1.3E1/4.0E1)+cos(q4+q5+theta+3.141592653589793*(1.0/2.0))*(1.3E1/4.0E1))-q5dot*cos(q4+q5+theta+3.141592653589793*(1.0/2.0))*(1.3E1/4.0E1);
  }
  if(strcmp(s,"vL")==0)
  {
    T[0] = rIxdot-thetadot*sin(theta-3.141592653589793*(1.0/2.0))*(1.1E1/2.0E1);
    T[1] = rIydot-thetadot*cos(theta-3.141592653589793*(1.0/2.0))*(1.1E1/2.0E1);
  }
  if(strcmp(s,"vCOM")==0)
  {
    T[0] = rIxdot-q1dot*cos(q1+q2+q3+theta)*1.415929203539823E-3-q2dot*cos(q1+q2+q3+theta)*1.415929203539823E-3-q3dot*cos(q1+q2+q3+theta)*1.415929203539823E-3-q4dot*cos(q4+q5+q6+theta)*1.415929203539823E-3-q5dot*cos(q4+q5+q6+theta)*1.415929203539823E-3-q6dot*cos(q4+q5+q6+theta)*1.415929203539823E-3-thetadot*cos(q1+q2+q3+theta)*1.415929203539823E-3-thetadot*cos(q4+q5+q6+theta)*1.415929203539823E-3+q1dot*sin(q1+q2+q3+theta)*1.415929203539823E-3+q2dot*sin(q1+q2+q3+theta)*1.415929203539823E-3+q3dot*sin(q1+q2+q3+theta)*1.415929203539823E-3+q4dot*sin(q4+q5+q6+theta)*1.415929203539823E-3+q5dot*sin(q4+q5+q6+theta)*1.415929203539823E-3+q6dot*sin(q4+q5+q6+theta)*1.415929203539823E-3+thetadot*sin(q1+q2+q3+theta)*1.415929203539823E-3+thetadot*sin(q4+q5+q6+theta)*1.415929203539823E-3-q1dot*cos(q1+theta)*3.451327433628319E-2-q4dot*cos(q4+theta)*3.451327433628319E-2-q7dot*cos(q7+theta)*(1.9E1/9.04E2)-q9dot*cos(q9+theta)*(1.9E1/9.04E2)-thetadot*cos(q1+theta)*3.451327433628319E-2-thetadot*cos(q4+theta)*3.451327433628319E-2-thetadot*cos(q7+theta)*(1.9E1/9.04E2)-thetadot*cos(q9+theta)*(1.9E1/9.04E2)+thetadot*cos(theta)*2.628318584070796E-1-q10dot*cos(q9+q10+theta)*1.095132743362832E-2-q1dot*cos(q1+q2+theta)*1.610619469026549E-2-q2dot*cos(q1+q2+theta)*1.610619469026549E-2-q4dot*cos(q4+q5+theta)*1.610619469026549E-2-q5dot*cos(q4+q5+theta)*1.610619469026549E-2-q7dot*cos(q7+q8+theta)*1.095132743362832E-2-q8dot*cos(q7+q8+theta)*1.095132743362832E-2-q9dot*cos(q9+q10+theta)*1.095132743362832E-2-thetadot*cos(q1+q2+theta)*1.610619469026549E-2-thetadot*cos(q4+q5+theta)*1.610619469026549E-2-thetadot*cos(q7+q8+theta)*1.095132743362832E-2-thetadot*cos(q9+q10+theta)*1.095132743362832E-2-sqrt(4.1E1)*q1dot*cos(q1+q2+q3+theta-8.960553845713439E-1)*3.539823008849558E-4-sqrt(4.1E1)*q2dot*cos(q1+q2+q3+theta-8.960553845713439E-1)*3.539823008849558E-4-sqrt(4.1E1)*q3dot*cos(q1+q2+q3+theta-8.960553845713439E-1)*3.539823008849558E-4-sqrt(4.1E1)*q4dot*cos(q4+q5+q6+theta-8.960553845713439E-1)*3.539823008849558E-4-sqrt(4.1E1)*q5dot*cos(q4+q5+q6+theta-8.960553845713439E-1)*3.539823008849558E-4-sqrt(4.1E1)*q6dot*cos(q4+q5+q6+theta-8.960553845713439E-1)*3.539823008849558E-4-sqrt(4.1E1)*thetadot*cos(q1+q2+q3+theta-8.960553845713439E-1)*3.539823008849558E-4-sqrt(4.1E1)*thetadot*cos(q4+q5+q6+theta-8.960553845713439E-1)*3.539823008849558E-4;
    T[1] = rIydot+q1dot*cos(q1+q2+q3+theta)*1.415929203539823E-3+q2dot*cos(q1+q2+q3+theta)*1.415929203539823E-3+q3dot*cos(q1+q2+q3+theta)*1.415929203539823E-3+q4dot*cos(q4+q5+q6+theta)*1.415929203539823E-3+q5dot*cos(q4+q5+q6+theta)*1.415929203539823E-3+q6dot*cos(q4+q5+q6+theta)*1.415929203539823E-3+thetadot*cos(q1+q2+q3+theta)*1.415929203539823E-3+thetadot*cos(q4+q5+q6+theta)*1.415929203539823E-3+q1dot*sin(q1+q2+q3+theta)*1.415929203539823E-3+q2dot*sin(q1+q2+q3+theta)*1.415929203539823E-3+q3dot*sin(q1+q2+q3+theta)*1.415929203539823E-3+q4dot*sin(q4+q5+q6+theta)*1.415929203539823E-3+q5dot*sin(q4+q5+q6+theta)*1.415929203539823E-3+q6dot*sin(q4+q5+q6+theta)*1.415929203539823E-3+thetadot*sin(q1+q2+q3+theta)*1.415929203539823E-3+thetadot*sin(q4+q5+q6+theta)*1.415929203539823E-3+q1dot*sin(q1+theta)*3.451327433628319E-2+q4dot*sin(q4+theta)*3.451327433628319E-2+q7dot*sin(q7+theta)*(1.9E1/9.04E2)+q9dot*sin(q9+theta)*(1.9E1/9.04E2)+thetadot*sin(q1+theta)*3.451327433628319E-2+thetadot*sin(q4+theta)*3.451327433628319E-2+thetadot*sin(q7+theta)*(1.9E1/9.04E2)+thetadot*sin(q9+theta)*(1.9E1/9.04E2)-thetadot*sin(theta)*2.628318584070796E-1+q10dot*sin(q9+q10+theta)*1.095132743362832E-2+q1dot*sin(q1+q2+theta)*1.610619469026549E-2+q2dot*sin(q1+q2+theta)*1.610619469026549E-2+q4dot*sin(q4+q5+theta)*1.610619469026549E-2+q5dot*sin(q4+q5+theta)*1.610619469026549E-2+q7dot*sin(q7+q8+theta)*1.095132743362832E-2+q8dot*sin(q7+q8+theta)*1.095132743362832E-2+q9dot*sin(q9+q10+theta)*1.095132743362832E-2+thetadot*sin(q1+q2+theta)*1.610619469026549E-2+thetadot*sin(q4+q5+theta)*1.610619469026549E-2+thetadot*sin(q7+q8+theta)*1.095132743362832E-2+thetadot*sin(q9+q10+theta)*1.095132743362832E-2+sqrt(4.1E1)*q1dot*sin(q1+q2+q3+theta-8.960553845713439E-1)*3.539823008849558E-4+sqrt(4.1E1)*q2dot*sin(q1+q2+q3+theta-8.960553845713439E-1)*3.539823008849558E-4+sqrt(4.1E1)*q3dot*sin(q1+q2+q3+theta-8.960553845713439E-1)*3.539823008849558E-4+sqrt(4.1E1)*q4dot*sin(q4+q5+q6+theta-8.960553845713439E-1)*3.539823008849558E-4+sqrt(4.1E1)*q5dot*sin(q4+q5+q6+theta-8.960553845713439E-1)*3.539823008849558E-4+sqrt(4.1E1)*q6dot*sin(q4+q5+q6+theta-8.960553845713439E-1)*3.539823008849558E-4+sqrt(4.1E1)*thetadot*sin(q1+q2+q3+theta-8.960553845713439E-1)*3.539823008849558E-4+sqrt(4.1E1)*thetadot*sin(q4+q5+q6+theta-8.960553845713439E-1)*3.539823008849558E-4;
  }
  return T;
}

double Kinetic_Energy_fn(Robot_StateNDot &Robot_StateNDot_i)
{
  std::vector<double> Robotstate_vector = StateNDot2StateVec(Robot_StateNDot_i);
  dlib::matrix<double> D_q, B_q, C_q_qdot, Jac_Full;
  Dynamics_Matrices(Robot_StateNDot_i, D_q, B_q, C_q_qdot, Jac_Full);
  dlib::matrix<double> Robotstatedot_Dlib, Robotstatedot_Dlib_Trans; Robotstatedot_Dlib = dlib::zeros_matrix<double>(13,1);
  for (int i = 0; i < 13; i++) {Robotstatedot_Dlib(i) = Robotstate_vector[13+i];}
  Robotstatedot_Dlib_Trans = dlib::trans(Robotstatedot_Dlib);
  double T;
  T = 0.5 * Robotstatedot_Dlib_Trans * D_q * Robotstatedot_Dlib;
  return T;
}

std::vector<double> Vec_Minus(std::vector<double> &vec1, std::vector<double> &vec2)
{
  int Vec_Length = vec1.size();

  std::vector<double> vec_minus;

  for (int i = 0; i < Vec_Length; i++)
  {
    vec_minus.push_back(vec1[i] - vec2[i]);
  }
  return vec_minus;
}

void End_Effector_Obs_Dist_Fn(dlib::matrix<double,12,1> &End_Effector_Pos, dlib::matrix<double,6,1> &End_Effector_Dist, std::vector<int> &End_Effector_Obs)
{
  std::vector<double> r_Pos(2);		double Obs_Dist_i;			int Obs_Dist_Index;
  for (int i = 0; i < 6; i++)
  {
    r_Pos[0] = End_Effector_Pos(2*i);
    r_Pos[1] = End_Effector_Pos(2*i+1);
    if(i<4)
    {
      Obs_Dist_Fn(r_Pos, Obs_Dist_i, Obs_Dist_Index, "floor");
    }
    else
    {
      // This is to directly choose the specific environmental feature for the hand contacts

      // Obs_Dist_Fn(r_Pos, Obs_Dist_i, Obs_Dist_Index, "None");
      Obs_Dist_Fn(r_Pos, Obs_Dist_i, Obs_Dist_Index, "floor");
      // Obs_Dist_Fn(r_Pos, Obs_Dist_i, Obs_Dist_Index, "wall");
    }
    End_Effector_Dist(i) = Obs_Dist_i;
    End_Effector_Obs[i] = Obs_Dist_Index;
  }
  return;
}

void End_Effector_PosNVel(Robot_StateNDot &StateNDot_Init_i, dlib::matrix<double,12,1> &End_Effector_Pos, dlib::matrix<double,12,1> &End_Effector_Vel)
{
  std::vector<double> rA = Ang_Pos_fn(StateNDot_Init_i, "rA");        std::vector<double> rB = Ang_Pos_fn(StateNDot_Init_i, "rB");
  std::vector<double> rC = Ang_Pos_fn(StateNDot_Init_i, "rC");        std::vector<double> rD = Ang_Pos_fn(StateNDot_Init_i, "rD");
  std::vector<double> rE = Ang_Pos_fn(StateNDot_Init_i, "rE");        std::vector<double> rF = Ang_Pos_fn(StateNDot_Init_i, "rF");

  End_Effector_Pos(0) = rA[0];                 	End_Effector_Pos(1) = rA[1];
  End_Effector_Pos(2) = rB[0];                 	End_Effector_Pos(3) = rB[1];
  End_Effector_Pos(4) = rC[0];                 	End_Effector_Pos(5) = rC[1];
  End_Effector_Pos(6) = rD[0];                 	End_Effector_Pos(7) = rD[1];
  End_Effector_Pos(8) = rE[0];                 	End_Effector_Pos(9) = rE[1];
  End_Effector_Pos(10) = rF[0];                 End_Effector_Pos(11) = rF[1];

  std::vector<double> vA = Ang_Vel_fn(StateNDot_Init_i, "vA");        std::vector<double> vB = Ang_Vel_fn(StateNDot_Init_i, "vB");
  std::vector<double> vC = Ang_Vel_fn(StateNDot_Init_i, "vC");        std::vector<double> vD = Ang_Vel_fn(StateNDot_Init_i, "vD");
  std::vector<double> vE = Ang_Vel_fn(StateNDot_Init_i, "vE");        std::vector<double> vF = Ang_Vel_fn(StateNDot_Init_i, "vF");

  End_Effector_Vel(0) = vA[0];                 	End_Effector_Vel(1) = vA[1];
  End_Effector_Vel(2) = vB[0];                 	End_Effector_Vel(3) = vB[1];
  End_Effector_Vel(4) = vC[0];                 	End_Effector_Vel(5) = vC[1];
  End_Effector_Vel(6) = vD[0];                 	End_Effector_Vel(7) = vD[1];
  End_Effector_Vel(8) = vE[0];                 	End_Effector_Vel(9) = vE[1];
  End_Effector_Vel(10) = vF[0];                	End_Effector_Vel(11) = vF[1];

  return;
}

void Node_UpdateNCon(Tree_Node &Node_i, Robot_StateNDot &Node_StateNDot_i, std::vector<double> &sigma)
{
  // This function is used to substitute the attribute and add it to the All_Nodes
  Node_i.Node_StateNDot = Node_StateNDot_i;
  Node_i.sigma = sigma;
  Node_i.KE = Kinetic_Energy_fn(Node_StateNDot_i);
  Node_i.Node_Index = All_Nodes.size();
  dlib::matrix<double,12,1> End_Effector_Pos, End_Effector_Vel;
  End_Effector_PosNVel(Node_StateNDot_i, End_Effector_Pos, End_Effector_Vel);
  Node_i.End_Effector_Pos = End_Effector_Pos;
  Node_i.End_Effector_Vel = End_Effector_Vel;
  All_Nodes.push_back(&Node_i);
  Frontier_Nodes.push_back(&Node_i);
  Frontier_Nodes_Cost.push_back(Node_i.KE);
}

double Nodes_Optimization_ObjNConstraint(std::vector<double> &Opt_Seed, std::vector<double> &ObjNConstraint_Val, std::vector<double> &ObjNConstraint_Type)
{
  // This function is the main optimization constraint function
  double T_tot, T, KE_End;;										dlib::matrix<double> StateNDot_Traj, Ctrl_Traj, Contact_Force_Traj;

  Opt_Seed_Unzip(Opt_Seed, T_tot, StateNDot_Traj, Ctrl_Traj, Contact_Force_Traj);

  T = T_tot/(Grids - 1) * 1.0;

  Tree_Node Node_i, Node_i_child;							Node_i = Structure_P.Node_i;						Node_i_child = Structure_P.Node_i_child;

  std::vector<double> Robotstate_Initial_Vec = StateNDot2StateVec(Node_i.Node_StateNDot);

  std::vector<double> sigma_i = Node_i.sigma;
  std::vector<double> sigma_i_child = Node_i_child.sigma;

  int Opt_Type_Flag = Sigma_TransNGoal(sigma_i, sigma_i_child);

  double KE_ref = 0.001;			// A sufficiently trivial value

  // Objective function initialization
  ObjNConstraint_Val.push_back(0);
  ObjNConstraint_Type.push_back(1);

  // 1. The first constraint is to make sure that the initial condition matches the given initial condition
  dlib::matrix<double> StateNDot_Traj_1st_colm, Matrix_result;
  StateNDot_Traj_1st_colm = dlib::colm(StateNDot_Traj, 0);
  // The first 26 constrait is set to be zero
  for (int i = 0; i < 26; i++)
  {
    ObjNConstraint_Val.push_back(StateNDot_Traj_1st_colm(i));
    ObjNConstraint_Type.push_back(1);
  }

  dlib::matrix<double> Contact_Active_Matrix = Contact_Status_fn(sigma_i);
  // Then is to find the active contact status where I am thinking that a trick might work which is when both feet are in contact...... anyway save taht for later

  // THe constraint will be added term by term
  dlib::matrix<double> Robostate_Dlib_Front,  Robostate_Dlib_Mid,   Robostate_Dlib_Back;
  dlib::matrix<double> Ctrl_Front,            Ctrl_Mid,             Ctrl_Back;
  dlib::matrix<double> Contact_Force_Front,   Contact_Force_Mid,    Contact_Force_Back;
  dlib::matrix<double> Robostatedot_Front,    Robostatedot_Mid,     Robostatedot_Back;

  // 2. The second constraint is Dynamics constraints
  // 2.a The numerical integration constraint at the knot points
  // 2.b The collocation constraint at the mid point where additional contact force will be introduced
  for (int i = 0; i < Grids-1; i++)
  {
    // Get the robot state, ctrl, and contact force at the front and end of each segment
    Robostate_Dlib_Front =  dlib::colm(StateNDot_Traj, i);	     Ctrl_Front = dlib::colm(Ctrl_Traj,i);       Contact_Force_Front = dlib::colm(Contact_Force_Traj,i);
    Robostate_Dlib_Back =   dlib::colm(StateNDot_Traj, i+1);	   Ctrl_Back = dlib::colm(Ctrl_Traj,i+1);      Contact_Force_Back = dlib::colm(Contact_Force_Traj,i+1);

    Robostatedot_Front =  State_Ctrl_CF_2_Statedot(Robostate_Dlib_Front,  Contact_Force_Front,  Ctrl_Front);
    Robostatedot_Back =   State_Ctrl_CF_2_Statedot(Robostate_Dlib_Back,   Contact_Force_Back,   Ctrl_Back);

    // First is to compute the robot state at the middle and the derivative at the middle
    Robostate_Dlib_Mid = 0.5 * (Robostatedot_Front + Robostate_Dlib_Back) + T/8.0 * (Robostatedot_Front - Robostatedot_Back);
    Ctrl_Mid = 0.5 * (Ctrl_Front + Ctrl_Back);
    Contact_Force_Mid = 0.5 * (Contact_Force_Front + Contact_Force_Back);
    Robostatedot_Mid =    State_Ctrl_CF_2_Statedot(Robostate_Dlib_Mid,    Contact_Force_Mid,    Ctrl_Mid);

    // Second is to impose the numerical integration constraint
    Matrix_result = Robostate_Dlib_Back - Robostate_Dlib_Front - T/6.0 * (Robostatedot_Front + 4 * Robostatedot_Mid + Robostatedot_Back);
    ObjNConstraint_ValNType_Update(Matrix_result, ObjNConstraint_Val, ObjNConstraint_Type, 0);

    // // Third is the manifold constraint on the acceleration since it is convenient to impose this constraint here
    // // Here the constraint on the last part will be taken care of individually
    // Matrix_result = Contact_Acc_Constraint(Robostate_Dlib_Front, Robostatedot_Front, Contact_Active_Matrix);
    // ObjNConstraint_ValNType_Update(Matrix_result, ObjNConstraint_Val, ObjNConstraint_Type, 0);

    // Last is the constraint on the change of state based on the integration of velocity and acceleration bounds
    for (int ss = 3; ss < 13; ss++)
    {
      // This is the bounds on the position change according to the velocity
      double xdot_max = max(Robostatedot_Front(ss), Robostatedot_Back(ss));
      double xdot_min = min(Robostatedot_Front(ss), Robostatedot_Back(ss));
      double Pos_Act_Change = Robostate_Dlib_Back(ss) - Robostate_Dlib_Front(ss);   // The actual change of the position
      // The change of the
      ObjNConstraint_Val.push_back(xdot_max*T - Pos_Act_Change);
      ObjNConstraint_Type.push_back(1);
      ObjNConstraint_Val.push_back(Pos_Act_Change - xdot_min*T);
      ObjNConstraint_Type.push_back(1);
    }
    // for (int aa = 0; aa < 13; aa++)
    // {
    //   // This is the bound on the velocity change according to the acceleration
    //   double xddot_max = max(Robostatedot_Front(aa + 13), Robostatedot_Back(aa + 13));
    //   double xddot_min = min(Robostatedot_Front(aa + 13), Robostatedot_Back(aa + 13));
    //   double Vel_Act_Change = Robostate_Dlib_Back(aa + 13) - Robostate_Dlib_Front(aa + 13);   // The actual change of the position
    //   // The change of the
    //   ObjNConstraint_Val.push_back(xddot_max*T - Vel_Act_Change);
    //   ObjNConstraint_Type.push_back(1);
    //   ObjNConstraint_Val.push_back(Vel_Act_Change - xddot_min*T);
    //   ObjNConstraint_Type.push_back(1);
    // }
  }

  // 3. Contact manifold constraints on position and velocty

  dlib::matrix<double,12,1> End_Effector_Pos, End_Effector_Vel;			            // Reservered for the constraint on the position and velocity

  dlib::matrix<double> End_Effector_Pos_ref;										                End_Effector_Pos_ref = Node_i.End_Effector_Pos;

  dlib::matrix<double> Robotstate_k;                                            Robot_StateNDot Robot_StateNDot_k;

  dlib::matrix<double> Robotstate_End, Contact_Force_End, Ctrl_End, Robostatedot_End;             Robot_StateNDot Robot_StateNDot_End;

  dlib::matrix<double> Pos_Maint_Constrait;

  for (int i = 0; i < Grids-1; i++)
  {
    // At the knot points, the position, velocity and acceleration will also be constrained
    // The default idea is that at least one contact point will be on the ground. As a result, the position constraint will be a comparison between the reference and the actual state

    // a. Get the end effector position at each state

    Robotstate_k = dlib::colm(StateNDot_Traj, i);
    Robot_StateNDot_k = DlibRobotstate2StateNDot(Robotstate_k);
    End_Effector_PosNVel(Robot_StateNDot_k, End_Effector_Pos, End_Effector_Vel);

    // Constraint on the position
    Pos_Maint_Constrait = End_Effector_Pos - End_Effector_Pos_ref;
    Matrix_result = Contact_Active_Matrix * Pos_Maint_Constrait;
    ObjNConstraint_ValNType_Update(Matrix_result, ObjNConstraint_Val, ObjNConstraint_Type, 0);

    // Constraint on the velocity
    Matrix_result = Contact_Active_Matrix * End_Effector_Vel;
    ObjNConstraint_ValNType_Update(Matrix_result, ObjNConstraint_Val, ObjNConstraint_Type, 0);

    // Constraint on the acceleration: this part has already been constrained in the previous dynamics constraint part
  }
  // However, the last grid is a special case which needs to be considered in a different manner
  Robotstate_End = dlib::colm(StateNDot_Traj, Grids - 1);
  Robot_StateNDot_End = DlibRobotstate2StateNDot(Robotstate_End);
  End_Effector_PosNVel(Robot_StateNDot_End, End_Effector_Pos, End_Effector_Vel);

  if (Opt_Type_Flag==0)
  {
    // This is the self-stabilizing process, so the contact constraint will be added in both the position, and velocity

    Pos_Maint_Constrait = End_Effector_Pos - End_Effector_Pos_ref;
    Matrix_result = Contact_Active_Matrix * Pos_Maint_Constrait;
    ObjNConstraint_ValNType_Update(Matrix_result, ObjNConstraint_Val, ObjNConstraint_Type, 0);

    // Constraint on the velocity
    Matrix_result = Contact_Active_Matrix * End_Effector_Vel;
    ObjNConstraint_ValNType_Update(Matrix_result, ObjNConstraint_Val, ObjNConstraint_Type, 0);

  }
  else
  {
    if(Opt_Type_Flag == -1)
    {
      // This is the contact removing case, where the position and velocity constraint will be satisfied while the acceleration is not satisfied by default
      Pos_Maint_Constrait = End_Effector_Pos - End_Effector_Pos_ref;
      Matrix_result = Contact_Active_Matrix * Pos_Maint_Constrait;
      ObjNConstraint_ValNType_Update(Matrix_result, ObjNConstraint_Val, ObjNConstraint_Type, 0);

      // Constraint on the velocity
      Matrix_result = Contact_Active_Matrix * End_Effector_Vel;
      ObjNConstraint_ValNType_Update(Matrix_result, ObjNConstraint_Val, ObjNConstraint_Type, 0);
    }
    else
    {
      // This is the impact mapping case, where the position constraint is the only constraint
      Pos_Maint_Constrait = End_Effector_Pos - End_Effector_Pos_ref;
      Matrix_result = Contact_Active_Matrix * Pos_Maint_Constrait;
      ObjNConstraint_ValNType_Update(Matrix_result, ObjNConstraint_Val, ObjNConstraint_Type, 0);
    }
  }

  // 3. Complementarity constraints: Distance
  //								   3.a All contact points have to remain nonnegative distance
  //								   3.b Certian distance should be zero when a contact is added

  dlib::matrix<double,6,1>  End_Effector_Dist;                                  // Reservered for the inequality constraint always remain nonnegative
  dlib::matrix<int> End_Effector_Obs_Index_Matrix = dlib::ones_matrix<int>(6, Grids);
  std::vector<int> End_Effector_Obs(6);
  dlib::matrix<double> ones_vector, temp_matrix;
  ones_vector = ONES_VECTOR_fn(6);
  double mini = 0.025;
  for (int i = 0; i < Grids; i++)
  {
    Robotstate_k = dlib::colm(StateNDot_Traj, i);
    Robot_StateNDot_k = DlibRobotstate2StateNDot(Robotstate_k);
    End_Effector_PosNVel(Robot_StateNDot_k, End_Effector_Pos, End_Effector_Vel);
    End_Effector_Obs_Dist_Fn(End_Effector_Pos, End_Effector_Dist, End_Effector_Obs);

    Matrix_result = End_Effector_Dist;
    ObjNConstraint_ValNType_Update(Matrix_result, ObjNConstraint_Val, ObjNConstraint_Type, 1);

    for (int j = 0; j < 6; j++)
    {
      End_Effector_Obs_Index_Matrix(j,i) = End_Effector_Obs[j];
    }
    // The middle joint can also be considered if possible
    temp_matrix = Middle_Joint_Obs_Dist_Fn(Robot_StateNDot_k);
    ones_vector(5) = 10;
    Matrix_result = temp_matrix - ones_vector * mini;
    ObjNConstraint_ValNType_Update(Matrix_result, ObjNConstraint_Val, ObjNConstraint_Type, 1);
  }

  // Impact mapping constraint
  if (Opt_Type_Flag == 1)
  {
    dlib::matrix<double> Dist_value = dlib::ones_matrix<double>(6,1);
    Dist_value(0) = End_Effector_Dist(0) * sigma_i_child[0];
    Dist_value(1) = End_Effector_Dist(1) * sigma_i_child[0];
    Dist_value(2) = End_Effector_Dist(2) * sigma_i_child[1];
    Dist_value(3) = End_Effector_Dist(3) * sigma_i_child[1];
    Dist_value(4) = End_Effector_Dist(4) * sigma_i_child[2];
    Dist_value(5) = End_Effector_Dist(5) * sigma_i_child[3];
    ObjNConstraint_ValNType_Update(Dist_value, ObjNConstraint_Val, ObjNConstraint_Type, 0);
  }

  // 4. Contact Force: Complementarity
  dlib::matrix<double> Contact_Force_Complem_Trans_Matrix = Contact_Force_Complem_Matrix_fn(sigma_i);
  dlib::matrix<double> Contact_Force_Complem_End_Matrix = Contact_Force_Complem_Matrix_fn(sigma_i_child);
  dlib::matrix<double> Contact_Force_Complem_Matrix;

  dlib::matrix<double> Contact_Force_k;

  for (int i = 0; i < Grids; i++)
  {
    if(i<Grids-1)
    {
      Contact_Force_Complem_Matrix = Contact_Force_Complem_Trans_Matrix;
    }
    else
    {
      Contact_Force_Complem_Matrix = Contact_Force_Complem_End_Matrix;
    }

    Contact_Force_k = dlib::colm(Contact_Force_Traj,i);
    Matrix_result = Contact_Force_Complem_Matrix * Contact_Force_k;
    ObjNConstraint_ValNType_Update(Matrix_result, ObjNConstraint_Val, ObjNConstraint_Type, 0);
  }

  // 5. Contact Force: Feasibility
  //                                           a. Normal force should be positive and
  //											                     b. The friction cone constraint has to be satisfied
  //											                     c. The horizontal force has to lie on the same side to avoid force oscilation (* This constraint may not be so important)

  dlib::matrix<int> End_Effector_Obs_k;

  for (int i = 0; i < Grids; i++)
  {
    Contact_Force_k = dlib::colm(Contact_Force_Traj,i);
    End_Effector_Obs_k = dlib::colm(End_Effector_Obs_Index_Matrix, i);
    Contact_Force_Feasibility_fn(Contact_Force_k, End_Effector_Obs_k, ObjNConstraint_Val, ObjNConstraint_Type);
  }

  // The cost function
  KE_End = Kinetic_Energy_End_Frame(StateNDot_Traj);
  ObjNConstraint_Val[0] = KE_End;

  // Here all constraints have been formulated and the cost function is written as follows.
  dlib::matrix<double> Robostate_Dlib_End;
  if(Opt_Type_Flag == 0)
  {
    // One constraint on the acceleration for inertial shaping to be zero
    Robostate_Dlib_End = dlib::colm(StateNDot_Traj, Grids-1);                     Robot_StateNDot_End = DlibRobotstate2StateNDot(Robostate_Dlib_End);

    // One more constraint to be added is the acceleration to be trivial value
    dlib::matrix<double> End_State_Traj = dlib::colm(StateNDot_Traj,Grids-1);			Robot_StateNDot End_StateNDot = DlibRobotstate2StateNDot(End_State_Traj);
    dlib::matrix<double> End_Ctrl_Traj = dlib::colm(Ctrl_Traj, Grids-1);				  dlib::matrix<double> End_Contact_Force_Traj = dlib::colm(Contact_Force_Traj, Grids-1);

    dlib::matrix<double> D_q_End, B_q_End, C_q_qdot_End, Jac_Full_End;            Dynamics_Matrices(End_StateNDot, D_q_End, B_q_End, C_q_qdot_End, Jac_Full_End);

    dlib::matrix<double> Trivial_Acc_Matrix, Jac_Full_Trans_End;						      Jac_Full_Trans_End = dlib::trans(Jac_Full_End);
    Trivial_Acc_Matrix = Jac_Full_Trans_End * End_Contact_Force_Traj + B_q_End * End_Ctrl_Traj - C_q_qdot_End;
    ObjNConstraint_ValNType_Update(Trivial_Acc_Matrix, ObjNConstraint_Val, ObjNConstraint_Type, 0);

    ObjNConstraint_Val.push_back(KE_ref - KE_End);
    ObjNConstraint_Type.push_back(1);
  }

  double State_Var = Traj_Variation(StateNDot_Traj);

  return State_Var;
}

void Contact_Force_Feasibility_fn(dlib::matrix<double> &Contact_Force_k, dlib::matrix<int> &End_Effector_Obs_k, std::vector<double> &ObjNConstraint_Val, std::vector<double> &ObjNConstraint_Type)
{
  // This function is used to impose the contraint on the contact force feasibility

  std::vector<double> Normal_Force_vec(6), 	Tange_Force_vec(6);

  Contact_Force_Proj(Contact_Force_k, Normal_Force_vec, Tange_Force_vec, End_Effector_Obs_k);

  for (int nn = 0; nn < 6; nn++)
  {
    // a. Normal force should be positive
    ObjNConstraint_Val.push_back(Normal_Force_vec[nn]);
    ObjNConstraint_Type.push_back(1);

    // b. Friction cone constraint
    ObjNConstraint_Val.push_back(Normal_Force_vec[nn] * Normal_Force_vec[nn] * mu * mu - Tange_Force_vec[nn] * Tange_Force_vec[nn]);
    ObjNConstraint_Type.push_back(1);

    // // c. Horizontal force change limitation?
    // ObjNConstraint_Val.push_back(Tange_Force_Front_vec[nn] * Tange_Force_Back_vec[nn]);
    // ObjNConstraint_Type.push_back(1);
  }
}

dlib::matrix<double> Contact_Force_Complem_Matrix_fn(std::vector<double> &sigma)
{
  std::vector<double> sigma_temp;       dlib::matrix<double> Contact_Force_Complem_Matrix;
  sigma_temp.push_back(!sigma[0]);			sigma_temp.push_back(!sigma[0]);			sigma_temp.push_back(!sigma[0]);				sigma_temp.push_back(!sigma[0]);
  sigma_temp.push_back(!sigma[1]);			sigma_temp.push_back(!sigma[1]);			sigma_temp.push_back(!sigma[1]);				sigma_temp.push_back(!sigma[1]);
  sigma_temp.push_back(!sigma[2]);			sigma_temp.push_back(!sigma[2]);			sigma_temp.push_back(!sigma[3]);				sigma_temp.push_back(!sigma[3]);
  Contact_Force_Complem_Matrix = Diag_Matrix_fn(sigma_temp);
  return Contact_Force_Complem_Matrix;
}

dlib::matrix<double> Contact_Status_fn(std::vector<double> &sigma_i)
{
  // This function is used to generaete the contact status function based on the Node_i and Node_i_child sigma
  // The contact status function only considers about the contact status during the integration period not the end which is the critial index

  // THis constraint should be satisfied in all the intermediate states: knot points
  std::vector<double> sigma_maint(12);

  sigma_maint[0] = sigma_i[0];							sigma_maint[1] = sigma_i[0];
  sigma_maint[2] = sigma_i[0];							sigma_maint[3] = sigma_i[0];
  sigma_maint[4] = sigma_i[1];							sigma_maint[5] = sigma_i[1];
  sigma_maint[6] = sigma_i[1];							sigma_maint[7] = sigma_i[1];
  sigma_maint[8] = sigma_i[2];							sigma_maint[9] = sigma_i[2];
  sigma_maint[10] = sigma_i[3];							sigma_maint[11] = sigma_i[3];

  dlib::matrix<double> Contact_Status_Matrix = Diag_Matrix_fn(sigma_maint);
  return Contact_Status_Matrix;
}

// dlib::matrix<double> Contact_Acc_Constraint(dlib::matrix<double> &Robotstate, dlib::matrix<double> &Robotstatedot, const dlib::matrix<double> &Sigma_Matrix)
// {
//   // This function is used to calculate the constraint on the acceleration at the contact points the outpu will be a dlib::matrix
//   Robot_StateNDot Robot_StateNDot_k = DlibRobotstate2StateNDot(Robotstate);
//   dlib::matrix<double> Jacdot_qdot_k = Jacdot_qdot_fn(Robot_StateNDot_k);
//   dlib::matrix<double> Jac_Full_k = Jac_Full_fn(Robot_StateNDot_k);
//   dlib::matrix<double> qddot = dlib::ones_matrix<double>(13,1);
//
//   for (int i = 0; i < 13; i++)
//   {
//     qddot(i) = Robotstatedot(i + 13);
//   }
//   dlib::matrix<double> Matrix_result = Jacdot_qdot_k + Jac_Full_k * qddot;
//
//   Matrix_result = Sigma_Matrix * Matrix_result;
//
//   return Matrix_result;
// }

dlib::matrix<double> Diag_Matrix_fn(std::vector<double> &diag_vec)
{
  // This function is used to generate a diagonal matrix within diag_vec to its diagonal elements
  int dim = diag_vec.size();
  dlib::matrix<double> Diag_Matrix;
  Diag_Matrix = dlib::zeros_matrix<double>(dim,dim);
  for (int i = 0; i < dim; i++)
  {
    Diag_Matrix(i,i) = diag_vec[i];
  }
  return Diag_Matrix;
}

void Dynamics_Matrices(const Robot_StateNDot &Node_StateNDot, dlib::matrix<double> &D_q, dlib::matrix<double> &B_q, dlib::matrix<double> &C_q_qdot, dlib::matrix<double> &Jac_Full)
{
  D_q = D_q_fn(Node_StateNDot);
  B_q = B_q_fn();
  C_q_qdot = C_q_qdot_fn(Node_StateNDot);
  Jac_Full = Jac_Full_fn(Node_StateNDot);
}

dlib::matrix<double> State_Ctrl_CF_2_Statedot(dlib::matrix<double> &Robotstate_k, dlib::matrix<double> &Contact_Force_k, dlib::matrix<double> &Control_k)
{
  // This function is used to calculate the statedot based on the information of the state, control and contact force

  dlib::matrix<double> D_q, B_q, 	C_q_qdot, Jac_Full, Jac_Full_Trans;
  dlib::matrix<double> Dynamics_LHS, Dynamics_RHS;      // The default way to write it is Aq'' = J^(T)Lambda + B * u - C(q,q')
  Robot_StateNDot Robot_StateNDot_k = DlibRobotstate2StateNDot(Robotstate_k);
  Dynamics_Matrices(Robot_StateNDot_k, D_q, B_q, C_q_qdot, Jac_Full);
  Jac_Full_Trans = dlib::trans(Jac_Full);

  // Here the update is for the whole system state. As a result, both the velocity and the position will be updated
  Dynamics_LHS = D_q;
  Dynamics_RHS = Jac_Full_Trans * Contact_Force_k + B_q * Control_k - C_q_qdot;
  dlib::matrix<double> qddot = dlib::inv(Dynamics_LHS) * Dynamics_RHS;

  // The last step is to give the value back to form the statedot vector
  dlib::matrix<double> Statedot = dlib::zeros_matrix<double>(26,1);    // Here 26 is 13 + 13
  for (int i = 0; i < 13; i++)
  {
    Statedot(i) = Robotstate_k(i+13);
  }
  for (int i = 0; i < 13; i++)
  {
    Statedot(i+13) = qddot(i);
  }
  return Statedot;
}

double Traj_Variation(dlib::matrix<double> &StateNDot_Traj)
{
  // This function is used to calcualte the variation of the state trajectory
  dlib::matrix<double> StateNDot_Traj_Front, StateNDot_Traj_Back, Matrix_result;
  double Traj_Variation_Val = 0.0;
  for (int i = 0; i < StateNDot_Traj.nc()-1; i++)
  {
    StateNDot_Traj_Front = dlib::colm(StateNDot_Traj, i);
    StateNDot_Traj_Back = dlib::colm(StateNDot_Traj, i+1);
    Matrix_result = StateNDot_Traj_Front - StateNDot_Traj_Back;
    for (int j = 0; j < Matrix_result.nr()/2; j++)
    {
      Traj_Variation_Val = Traj_Variation_Val +  Matrix_result(j) * Matrix_result(j);
    }
  }
  return Traj_Variation_Val;
}

double Kinetic_Energy_End_Frame(dlib::matrix<double> &StateNDot_Traj)
{
  dlib::matrix<double> End_Robotstate;
  End_Robotstate = dlib::colm(StateNDot_Traj, Grids-1);
  Robot_StateNDot Robot_StateNDot_i = DlibRobotstate2StateNDot(End_Robotstate);
  double KE_i = Kinetic_Energy_fn(Robot_StateNDot_i);
  return KE_i;
}

double KE_Variation_fn(dlib::matrix<double> &StateNDot_Traj)
{
  // This function should describe the change of the kinetic energy instead of the sum
  double KE_Variation = 0.0;
  std::vector<double> KE_tot;
  dlib::matrix<double> Robostate_Dlib_i;
  Robot_StateNDot Robot_StateNDot_i;
  for (int i = 0; i < StateNDot_Traj.nc(); i++)
  {
    Robostate_Dlib_i = dlib::colm(StateNDot_Traj, i);
    Robot_StateNDot_i = DlibRobotstate2StateNDot(Robostate_Dlib_i);
    KE_tot.push_back(Kinetic_Energy_fn(Robot_StateNDot_i));
  }
  for (int i = 0; i < KE_tot.size()-1; i++)
  {
    KE_Variation = KE_Variation + (KE_tot[i+1] - KE_tot[i]) * (KE_tot[i+1] - KE_tot[i]);
  }
  return KE_Variation;
}

void Contact_Force_Proj(dlib::matrix<double> &Contact_Force_k, std::vector<double> &Normal_Force_vec, std::vector<double> &Tange_Force_vec, dlib::matrix<int> & End_Effector_Obs_k)
{
  // This function is used to project the contact force in to the normal and tangential direction
  // The default is 6 contact forces
  double Contact_Force_i_x, Contact_Force_i_y, Contact_Force_i_Normal, Contact_Force_i_Tange;
  for (int j = 0; j < 6; j++)
  {
    Contact_Force_i_x = Contact_Force_k(2*j);											Contact_Force_i_y = Contact_Force_k(2*j+1);
    Contact_Force_i_Normal = Contact_Force_i_x * Envi_Map_Normal(End_Effector_Obs_k(j),0) + Contact_Force_i_y * Envi_Map_Normal(End_Effector_Obs_k(j),1);
    Contact_Force_i_Tange = Contact_Force_i_x * Envi_Map_Tange(End_Effector_Obs_k(j),0)  + Contact_Force_i_y * Envi_Map_Tange(End_Effector_Obs_k(j),1);
    Normal_Force_vec[j] = Contact_Force_i_Normal;
    Tange_Force_vec[j] = Contact_Force_i_Tange;
  }
}

int Nodes_Optimization_fn(Tree_Node &Node_i, Tree_Node &Node_i_child, std::vector<double> &Opt_Soln_Output)
{
  // This function will optimize the joint trajectories to minimize the robot kinetic energy while maintaining a smooth transition fashion
  // However, the constraint will be set to use the direct collocation method

  int Opt_Flag = 0;                Structure_P.Node_i = Node_i;		Structure_P.Node_i_child = Node_i_child;

  double Time_Interval = 0.02;		 int Total_Num = 41;

  std::vector<double> Time_Queue = Time_Seed_Queue_fn(Time_Interval, Total_Num);       // Here Time_Queue will output the queue of the h_k but the optimization variable is h_k_Tot

  std::vector<double> Opt_Soln_Tot, State_Variation_Tot;

  int Feasible_Num = 0;

  for (int j = 0; j < Time_Queue.size(); j++)
  {
    std::vector<double> Opt_Soln, ObjNConstraint_Val, ObjNConstraint_Type;

    h_k = Time_Queue[j];								// I forget what is the meaning of my Time_Seed but I would say that is the duration between two knots

    cout<<"--------------------------------- Time Guess Iteration: "<<j<<" --------------------------------- "<<endl;

    cout<<"The Guess for the time duaration h_k is :" <<h_k<<" s"<<endl;

    Opt_Soln = Nodes_Optimization_Inner_Opt(Node_i, Node_i_child);

    double State_Variation = Nodes_Optimization_ObjNConstraint(Opt_Soln, ObjNConstraint_Val, ObjNConstraint_Type);
    State_Variation_Tot.push_back(State_Variation);

    double ObjNConstraint_Violation_Val = ObjNConstraint_Violation(ObjNConstraint_Val, ObjNConstraint_Type);
    cout<<endl;
    cout<<"-----------------------------ObjNConstraint_Violation_Val: "<<ObjNConstraint_Violation_Val<<"-------------------------------"<<endl;

    if(ObjNConstraint_Violation_Val<0.1)
    {
      Opt_Flag = 1;

      Feasible_Num = Feasible_Num + 1;

      // However, due to the small constraint violation, the robot may not 100% satisfy the holonomic constraint so here we would like to check the robot final state
      // Opt_Soln =  Final_State_Opt(Opt_Soln, Node_i, Node_i_child);

      for (int i = 0; i < Opt_Soln.size(); i++)
      {
        Opt_Soln_Tot.push_back(Opt_Soln[i]);
      }
      time_t now = time(0);
      // convert now to string form
      char* dt = ctime(&now);

      ofstream output_file;
      std::string pre_filename = "From_Node_";
      std::string Node_i_name = to_string(Node_i.Node_Index);
      std::string mid_filename = "_Expansion_At_Feasible_";
      std::string Node_i_Child_name = to_string(Feasible_Num);
      std::string post_filename = ".txt";

      std::string filename = pre_filename + Node_i_name + mid_filename + Node_i_Child_name + dt + post_filename;
      output_file.open(filename, std::ofstream::out);
      for (int i = 0; i < Opt_Soln.size(); i++)
      {
        output_file<<Opt_Soln[i]<<endl;
      }
      output_file.close();
    }
    else
    {
      // The differentiation between the pure seed guess and the failure guess needs to be conducted
      // Then this is a failure attempt
      if(ObjNConstraint_Violation_Val<10)
      {
        time_t now = time(0);
        // convert now to string form
        char* dt = ctime(&now);
        ofstream output_file;
        std::string pre_filename = "From_Node_";
        std::string Node_i_name = to_string(Node_i.Node_Index);
        std::string mid_filename = "_Expansion_Quasi_At_";
        std::string post_filename = ".txt";

        std::string filename = pre_filename + Node_i_name + mid_filename + dt + post_filename;
        output_file.open(filename, std::ofstream::out);
        for (int i = 0; i < Opt_Soln.size(); i++)
        {
          output_file<<Opt_Soln[i]<<endl;
        }
        output_file.close();
      }
    }
  }

  // Which means that there exits feasible solutions

  if(Opt_Flag == 1)
  {
    // Now it is time to select the best one from all feasible solutions. The best solution is the one whose state drift is the minimum

    int Opt_Soln_Index = Minimum_Index(State_Variation_Tot);

    time_t now = time(0);
    // convert now to string form
    char* dt = ctime(&now);

    ofstream output_file;
    std::string pre_filename = "From_Node_";
    std::string Node_i_name = to_string(Node_i.Node_Index);
    std::string mid_filename = "_Expansion_At_";
    std::string Node_i_Child_name = dt;
    std::string post_filename = "_Opt_Soln.txt";

    std::string filename = pre_filename + Node_i_name + mid_filename + Node_i_Child_name + post_filename;

    output_file.open(filename, std::ofstream::out);

    int Opt_Soln_Index_Start = Variable_Num * (Opt_Soln_Index - 1);


    for (int i = 0; i < Variable_Num; i++)
    {
      output_file<<Opt_Soln_Tot[i + Opt_Soln_Index_Start]<<endl;
      Opt_Soln_Output[i] = Opt_Soln_Tot[i + Opt_Soln_Index_Start];
    }
    output_file.close();

  }
  return Opt_Flag;
}

std::vector<double> Nodes_Optimization_Inner_Opt(Tree_Node &Node_i, Tree_Node &Node_i_child)
{
  std::vector<double> Opt_Seed = Seed_Guess_Gene(Node_i, Node_i_child);

  std::vector<double> ObjNConstraint_Val, ObjNConstraint_Type;

  Nodes_Optimization_ObjNConstraint(Opt_Seed, ObjNConstraint_Val, ObjNConstraint_Type);

  snoptProblem Nodes_Optimization_Pr;                     // This is the name of the Optimization problem for the robot configuration

  integer n = Opt_Seed.size();
  integer neF = ObjNConstraint_Val.size();
  integer lenA  =  n * neF;

  Structure_P.Opt_Val_No = Opt_Seed.size();			     Structure_P.ObjNConst_No = ObjNConstraint_Val.size();

  integer *iAfun = new integer[lenA];              	  integer *jAvar = new integer[lenA];					   doublereal *A  = new doublereal[lenA];
  integer lenG   = lenA;								              integer *iGfun = new integer[lenG];					   integer *jGvar = new integer[lenG];
  doublereal *x      = new doublereal[n];				      doublereal *xlow   = new doublereal[n];				 doublereal *xupp   = new doublereal[n];
  doublereal *xmul   = new doublereal[n];				      integer    *xstate = new    integer[n];

  doublereal *F      = new doublereal[neF];			      doublereal *Flow   = new doublereal[neF];			 doublereal *Fupp   = new doublereal[neF];
  doublereal *Fmul   = new doublereal[neF];			      integer    *Fstate = new integer[neF];

  integer nxnames = 1;								                integer nFnames = 1;					char *xnames = new char[nxnames*8];					char *Fnames = new char[nFnames*8];
  integer    ObjRow = 0;								              doublereal ObjAdd = 0;

  for (int i = 0; i < n; i++)
  {
    xlow[i] =-Inf;
    xupp[i] = Inf;
  }

  xlow[0] = 0.25;		xupp[0] = 3.5;        // The first value is the time

  int Index_Count = 1;

  for (int i = 0; i < Grids; i++)
  {
    for (int j = 0; j < 26; j++)
    {
      xlow[Index_Count] = xlow_vec(j);
      xupp[Index_Count] = xupp_vec(j);
      Index_Count = Index_Count + 1;
    }
  }

  for (int i = 0; i < Grids; i++)
  {
    for (int j = 0; j < 10; j++)
    {
      xlow[Index_Count] = ctrl_low_vec(j);
      xupp[Index_Count] = ctrl_upp_vec(j);
      Index_Count = Index_Count + 1;
    }
  }

  for (int i = 0; i < n; i++)
  {
    xstate[i] = 0.0;
    x[i] = Opt_Seed[i];  	// Initial guess
  }

  for(int i = 0; i<neF; i++)
  {
    if (ObjNConstraint_Type[i]>1.0)
    {
      // This seems like a constraint on the acceleration. However, I would say curently the compact version is not to consider the qddot
      Flow[i] = -Acc_max;
      Fupp[i] = Acc_max;
    }
    else
    {
      Flow[i] = 0.0;
      if(ObjNConstraint_Type[i]>0.0)
      {
        Fupp[i] = Inf;
      }
      else
      {
        Fupp[i] = 0.0;
      }
    }
  }

  // Here the linear initial condition matching is conducted as bounds
  std::vector<double> Init_Config = StateNDot2StateVec(Node_i.Node_StateNDot);

  for (int i = 0; i < 26; i++)
  {
    Flow[i+1] = Init_Config[i];
    Fupp[i+1] = Init_Config[i];
  }

  // Load the data for ToyProb ...
  Nodes_Optimization_Pr.setPrintFile  ( "Nodes_Optimization_Pr.out" );
  Nodes_Optimization_Pr.setProblemSize( n, neF );
  Nodes_Optimization_Pr.setObjective  ( ObjRow, ObjAdd );
  Nodes_Optimization_Pr.setA          ( lenA, iAfun, jAvar, A );
  Nodes_Optimization_Pr.setG          ( lenG, iGfun, jGvar );
  Nodes_Optimization_Pr.setX          ( x, xlow, xupp, xmul, xstate );
  Nodes_Optimization_Pr.setF          ( F, Flow, Fupp, Fmul, Fstate );
  Nodes_Optimization_Pr.setXNames     ( xnames, nxnames );
  Nodes_Optimization_Pr.setFNames     ( Fnames, nFnames );
  Nodes_Optimization_Pr.setProbName   ( "Nodes_Optimization_Pr" );
  Nodes_Optimization_Pr.setUserFun    ( Nodes_Optimization_Pr_fn);
  // snopta will compute the Jacobian by finite-differences.
  // The user has the option of calling  snJac  to define the
  // coordinate arrays (iAfun,jAvar,A) and (iGfun, jGvar).
  Nodes_Optimization_Pr.computeJac    ();
  Nodes_Optimization_Pr.setIntParameter( "Derivative option", 0 );
  Nodes_Optimization_Pr.setIntParameter("Major iterations limit", 250);
  Nodes_Optimization_Pr.setIntParameter("Minor iterations limit", 50000);
  Nodes_Optimization_Pr.setIntParameter("Iterations limit", 2000000);
  Nodes_Optimization_Pr.setIntParameter("Minor print level", 1);
  // Nodes_Optimization_Pr.setSpecsFile   ( ".\Nodes_Optimization_Pr.spc" );

  // string QPsolver_opt = "QN";
  // Nodes_Optimization_Pr.setRealParameter("QPsolver", QPsolver_opt);


  integer Cold = 0, Basis = 1, Warm = 2;
  Nodes_Optimization_Pr.solve( Cold );

  std::vector<double> Opt_Soln;

  for (int i = 0; i < n; i++)
  {
    Opt_Soln.push_back(x[i]);
  }
  delete []iAfun;  delete []jAvar;  delete []A;
  delete []iGfun;  delete []jGvar;

  delete []x;      delete []xlow;   delete []xupp;
  delete []xmul;   delete []xstate;

  delete []F;      delete []Flow;   delete []Fupp;
  delete []Fmul;   delete []Fstate;

  delete []xnames; delete []Fnames;

  return Opt_Soln;
}

int Nodes_Optimization_Pr_fn(integer    *Status, integer *n,    doublereal x[],
     integer    *needF,  integer *neF,  doublereal F[],
     integer    *needG,  integer *neG,  doublereal G[],
     char       *cu,     integer *lencu,
	 integer    iu[],    integer *leniu,
	 doublereal ru[],    integer *lenru )
{
  std::vector<double> Opt_Seed, ObjNConstraint_Val, ObjNConstraint_Type;
  const int Opt_Val_No = Structure_P.Opt_Val_No ;
  const int ObjNConst_No = Structure_P.ObjNConst_No;
  for (int i = 0; i < Opt_Val_No; i++)
  {
    Opt_Seed.push_back(x[i]);
  }
  Nodes_Optimization_ObjNConstraint(Opt_Seed, ObjNConstraint_Val, ObjNConstraint_Type);
  for (int i = 0; i < ObjNConst_No; i++)
  {
    F[i] = ObjNConstraint_Val[i];
  }
  return 0;
}

double ObjNConstraint_Violation(const std::vector<double> &ObjNConstraint_Val, const std::vector<double> &ObjNConstraint_Type)
{
  double ObjNConstraint_Violation_Val = 0.0;

  for (int i = 0; i < ObjNConstraint_Val.size(); i++)
  {
    if(ObjNConstraint_Type[i]==0)
    {
      if(abs(ObjNConstraint_Val[i])>ObjNConstraint_Violation_Val)
      {
        ObjNConstraint_Violation_Val = abs(ObjNConstraint_Val[i]);
      }
    }
  }
  return ObjNConstraint_Violation_Val;
}

void Opt_Soln_Write2Txt(Tree_Node &Node_i,Tree_Node &Node_i_child, std::vector<double> &Opt_Soln)
{
  ofstream output_file;
  std::string pre_filename = "From_Node_";
  std::string Node_i_name = to_string(Node_i.Node_Index);
  std::string mid_filename = "_To_Node_";
  std::string Node_i_Child_name = to_string(Node_i_child.Node_Index);
  std::string post_filename = "_Opt_Soln.txt";
  std::string filename = pre_filename + Node_i_name + mid_filename + Node_i_Child_name + post_filename;
  output_file.open(filename, std::ofstream::out);
  for (int i = 0; i < Opt_Soln.size(); i++)
  {
    output_file<<Opt_Soln[i]<<endl;
  }
  output_file.close();
}

dlib::matrix<double> Node_Expansion_fn(const Tree_Node &Node_i, int &Adjacent_Number)
{
  // Checked! Oct.4th.2018 9:46PM
  // This function is used to conduct the node expansion for a given parent node
  // The basic consideration is not to have the flying-in-air phase

  std::vector<double> sigma_i = Node_i.sigma;
  dlib::matrix<double> Nodes_Sigma_Matrix;
  Nodes_Sigma_Matrix = dlib::zeros_matrix<double>(4,4);
  double sigma_sum;
  Adjacent_Number = 0;
  for (int i = 0; i < 4; i++)
  {
    std::vector<double> sigma_t = sigma_i;
    sigma_t[i] = !sigma_t[i];
    sigma_sum = sigma_t[0] + sigma_t[1] + sigma_t[2] + sigma_t[3];
    if(sigma_sum==0)
    {
      continue;
    }
    else
    {
      Nodes_Sigma_Matrix(Adjacent_Number, 0) = sigma_t[0];
      Nodes_Sigma_Matrix(Adjacent_Number, 1) = sigma_t[1];
      Nodes_Sigma_Matrix(Adjacent_Number, 2) = sigma_t[2];
      Nodes_Sigma_Matrix(Adjacent_Number, 3) = sigma_t[3];
      Adjacent_Number = Adjacent_Number + 1;
    }
  }
  // cout<<Nodes_Sigma_Matrix<<endl;
  return Nodes_Sigma_Matrix;
}

std::vector<double> End_RobotNDot_Extract(std::vector<double> &Opt_Soln, std::vector<double> &sigma_i, std::vector<double>&sigma_i_child)
{
  // This function is used to extract the end robot state out from the Opt_Soln while conducting Impact Mapping if there exists
  double T_tot;
  dlib::matrix<double> StateNDot_Traj, Ctrl_Traj, Contact_Force_Traj;
  Opt_Seed_Unzip(Opt_Soln, T_tot, StateNDot_Traj, Ctrl_Traj, Contact_Force_Traj);
  int Opt_Type_Flag = Sigma_TransNGoal(sigma_i, sigma_i_child);

  std::vector<double> Robotstate_vec;
  if(Opt_Type_Flag ==1)
  {
    // In this case, there is impact mapping needed to be conducted
    Robot_StateNDot Robot_End_State;
    double Impulse_Mag;
    Robot_End_State = Impact_Mapping_fn(StateNDot_Traj, Impulse_Mag, sigma_i_child);
    Robotstate_vec = StateNDot2StateVec(Robot_End_State);
  }
  else
  {
    dlib::matrix<double> End_RobotstateDlib;
    End_RobotstateDlib = dlib::colm(StateNDot_Traj, Grids-1);
    for (int i = 0; i < End_RobotstateDlib.nr(); i++)
    {
      Robotstate_vec.push_back(End_RobotstateDlib(i));
    }
  }
  return Robotstate_vec;
}

Robot_StateNDot Impact_Mapping_fn(dlib::matrix<double> &End_State_Dlib, double &Impulse_Mag, std::vector<double> &sigma_i_child)
{
  // Since the impact mapping only happens at the end of the state so there is no need to take in all the state trajectories into consideration
  Robot_StateNDot End_StateNDot = DlibRobotstate2StateNDot(End_State_Dlib);
  std::vector<double> End_Statevec = StateNDot2StateVec(End_StateNDot);
  dlib::matrix<double> D_q, B_q, C_q_qdot, Jac_Full;
  Dynamics_Matrices(End_StateNDot, D_q, B_q, C_q_qdot, Jac_Full);

  // Here sigma_i_child is used to select the full row rank Jacobian matrix
  std::vector<double> Jac_Act_Index = Full_Row_Rank_Index(sigma_i_child);
  const int Jac_Row_Number = Jac_Act_Index.size();
  dlib::matrix<double> Jac_Act;					Jac_Act = dlib::zeros_matrix<double>(Jac_Row_Number,13);
  for (int j = 0; j < Jac_Row_Number; j++)
  {
    for (int k = 0; k < 13; k++)
    {
      Jac_Act(j,k) = Jac_Full(Jac_Act_Index[j],k);
    }
  }
  dlib::matrix<double> Jac_Act_Trans = dlib::trans(Jac_Act);
  dlib::matrix<double> D_q_Inv = dlib::inv(D_q);
  dlib::matrix<double> Pre_Impact_Vel, Post_Impact_Vel;
  Pre_Impact_Vel = dlib::zeros_matrix<double>(13,1);
  Post_Impact_Vel = Pre_Impact_Vel;
  for (int i = 0; i < 13; i++)
  {
    Pre_Impact_Vel(i) = End_Statevec[i+13];
  }
  dlib::matrix<double> Impulse_Lamda;
  // cout<<D_q<<endl;			cout<<B_q<<endl;			cout<<C_q_qdot<<endl;			cout<<Jac_Act<<endl;
  Impulse_Lamda = -dlib::pinv(Jac_Act * D_q_Inv * Jac_Act_Trans) * Jac_Act * Pre_Impact_Vel;
  // cout<<Impulse_Lamda<<endl;
  Post_Impact_Vel = D_q_Inv * Jac_Act_Trans * Impulse_Lamda + Pre_Impact_Vel;
  // cout<<Post_Impact_Vel<<endl;
  std::vector<double> Robotstate_End(26);
  for (int i = 0; i < 13; i++)
  {
    Robotstate_End[i] = End_Statevec[i];
    Robotstate_End[i+13] = Post_Impact_Vel(i);
  }
  Robot_StateNDot Robot_StateNDot_End;
  Robot_StateNDot_End = StateVec2StateNDot(Robotstate_End);

  Impulse_Mag = 0.0;
  for (int i = 0; i < Impulse_Lamda.nr(); i++)
  {
    Impulse_Mag = Impulse_Mag + Impulse_Lamda(i) * Impulse_Lamda(i);
  }
  return Robot_StateNDot_End;
}

std::vector<double> Full_Row_Rank_Index(std::vector<double> &sigma_i_child)
{
  // This function is used select the full row rank matrix out from the Full jacobian matrix
  std::vector<double> Row_Rank_Index_Vec;
  if(sigma_i_child[0]==1)
  {
    Row_Rank_Index_Vec.push_back(0);
    Row_Rank_Index_Vec.push_back(1);
    Row_Rank_Index_Vec.push_back(3);
  }
  if(sigma_i_child[1]==1)
  {
    Row_Rank_Index_Vec.push_back(4);
    Row_Rank_Index_Vec.push_back(5);
    Row_Rank_Index_Vec.push_back(7);
  }
  if(sigma_i_child[2]==1)
  {
    Row_Rank_Index_Vec.push_back(8);
    Row_Rank_Index_Vec.push_back(9);
  }
  if(sigma_i_child[3]==1)
  {
    Row_Rank_Index_Vec.push_back(10);
    Row_Rank_Index_Vec.push_back(11);
  }
  return Row_Rank_Index_Vec;
}

std::vector<double> Opt_Soln_Load()
{
  // This function is used to load the computed optimal solution for data analysis
  std::vector<double> Opt_Seed;
  ifstream Opt_Soln_File;              // This is to read the initial angle and angular velocities
  Opt_Soln_File.open("4.txt");
  if(Opt_Soln_File.is_open())
  {
    double data_each_line = 0.0;

    while(Opt_Soln_File>>data_each_line)
    {
      Opt_Seed.push_back(data_each_line);
    }
    Opt_Soln_File.close();
  }
  return Opt_Seed;
}

std::vector<double> Seed_Guess_Gene(Tree_Node &Node_i, Tree_Node &Node_i_child)
{
  // This function will generate the spline coefficients needed for the further optimization
  // Here Time_Seed is a global variable whose value is changed in the outer loop

  double T = h_k;

  // The first step is to generate a feasible configuration that can satisfy the contact mode in the node i child
  std::vector<double> Init_Config = StateNDot2StateVec(Node_i.Node_StateNDot);
  std::vector<double> Seed_Config = Seed_Guess_Gene_Robotstate(Node_i, Node_i_child);
  const int StateNDot_len = Init_Config.size();			const int Control_len = 10;			const int Contact_Force_len = 12;
  std::vector<double> sigma_i = Node_i.sigma;
  std::vector<double> sigma_i_child = Node_i_child.sigma;
  std::vector<double> sigma;

  dlib::matrix<double> StateNDot_Traj, Acc_Traj, StateNdot_Mid_Traj, Acc_Mid_Traj;
  StateNDot_Traj = dlib::zeros_matrix<double>(StateNDot_len, Grids);				// each variable traj will be in a row fashion
  Acc_Traj = dlib::zeros_matrix<double>(13, Grids);				                  // each variable traj will be in a row fashion
  StateNdot_Mid_Traj = dlib::zeros_matrix<double>(StateNDot_len, Grids-1);
  Acc_Mid_Traj = dlib::zeros_matrix<double>(13, Grids-1);

  dlib::matrix<double> Ctrl_Traj, Contact_Force_Traj;
  Ctrl_Traj = dlib::zeros_matrix<double>(Control_len, Grids);
  Contact_Force_Traj = dlib::zeros_matrix<double>(Contact_Force_len, Grids);

  // Then there are actually three methods to initialize the seed for optimizations
  // The difference lies in the way that the state is treated.

  // Method 1: The whole states are chosen in a linear fashion. Then the acceleration is found by treating each segment between grids to be cubic spline.
  // Method 2: Each whole state is chosen in a cubic spline.    Then the acceleration is found by evaluating this cubic spline at each knot pointd.
  // Method 3: The position is chosen in a quadratic fashion while the velocity is a linear function. Then the acceleration is a constant value which match (v0 - vN)/T_tot

  // Method 4: The position is chosen in  a quadratic fashion while there are some manifold projection has to be conducted to ensure the initial point and the end point on the manifold

  // No matter what method is chosen, the output from this step are four things: States/Accelerations at grids and States/Accelerations at mid-points

  // State_Traj_Interpolater(Init_Config, Seed_Config, StateNDot_Traj, Acc_Traj, StateNdot_Mid_Traj, Acc_Mid_Traj,1);
  // State_Traj_Interpolater(Init_Config, Seed_Config, StateNDot_Traj, Acc_Traj, StateNdot_Mid_Traj, Acc_Mid_Traj,2);
  // State_Traj_Interpolater(Init_Config, Seed_Config, StateNDot_Traj, Acc_Traj, StateNdot_Mid_Traj, Acc_Mid_Traj,3);

  // Here the contact sigma has to be specified hreer the sigma should match the previuos sigma idea

  int Opt_Type_Flag = Sigma_TransNGoal(sigma_i, sigma_i_child);
  if(Opt_Type_Flag == 1)
  {
    sigma = sigma_i_child;
  }
  else
  {
    sigma = sigma_i;
  }

  // Sigma_Act_Manifold substitution
  Sigma_Act_Manifold[0] = sigma[0];
  Sigma_Act_Manifold[1] = sigma[0];
  Sigma_Act_Manifold[2] = sigma[1];
  Sigma_Act_Manifold[3] = sigma[1];
  Sigma_Act_Manifold[4] = sigma[2];
  Sigma_Act_Manifold[5] = sigma[3];

  State_Traj_Interpolater_Manifold(Init_Config, sigma, StateNDot_Traj, Acc_Traj, StateNdot_Mid_Traj, Acc_Mid_Traj);

  cout<<StateNDot_Traj<<endl;
  cout<<Acc_Traj<<endl;


  // Second is to initialize: Control, and Contact Force
  // Since the acceleration is directly given in advance. This makes the initialization easier.

  dlib::matrix<double> State_k, Acc_k, Lamda_k, u_k;
  Lamda_k = dlib::zeros_matrix<double>(12,1);
  u_k = dlib::zeros_matrix<double>(10,1);

  // Contact Force Traj and Ctrl Traj
  for (int i = 0; i < Grids; i++)
  {
    State_k = dlib::colm(StateNDot_Traj, i);

    Acc_k = dlib::colm(Acc_Traj, i);

    StateNAcc2ContactForceNTorque(State_k, Acc_k, Lamda_k, u_k);

    // This is the evaluation of the contact force
    for (int j = 0; j < 12; j++)
    {
      Contact_Force_Traj(j,i) = Lamda_k(j);
    }

    // This is the evaluation of the control torque
    for (int j = 0; j < 10; j++)
    {
      Ctrl_Traj(j,i) = u_k(j);
    }
  }

  // The final task is to pile them into a single vector
  std::vector<double> Opt_Seed;
  Opt_Seed.push_back(T * (Grids - 1));
  Opt_Seed_Zip(Opt_Seed, StateNDot_Traj, Ctrl_Traj, Contact_Force_Traj);

  // cout<<StateNDot_Traj<<endl;				cout<<Ctrl_Traj<<endl;        cout<<Contact_Force_Traj<<endl;

  return Opt_Seed;
}
void StateNAcc2ContactForceNTorque(dlib::matrix<double> &State_k, dlib::matrix<double> &Acc_k, dlib::matrix<double> &Lamda_k, dlib::matrix<double> &u_k)
{
  // This function is sued to initialize the contact force and torque vector based on the information of State and Acc

  Robot_StateNDot Robot_StateNDot_k = DlibRobotstate2StateNDot(State_k);
  dlib::matrix<double> D_q_k, B_q_k, C_q_qdot_k, Jac_Full_k;
  Dynamics_Matrices(Robot_StateNDot_k, D_q_k, B_q_k, C_q_qdot_k, Jac_Full_k);

  dlib::matrix<double> Dynamics_LHS = D_q_k * Acc_k + C_q_qdot_k;
  dlib::matrix<double> Dynamics_RHS_Matrix = Dynamics_RHS_Matrix_fn(Jac_Full_k, B_q_k);
  dlib::matrix<double> Dynamics_RHS = dlib::pinv(Dynamics_RHS_Matrix) * Dynamics_LHS;
  for (int j = 0; j < Dynamics_RHS.nr(); j++) {
    if(j<12)
    {
      Lamda_k(j) = Dynamics_RHS(j);
    }
    else
    {
      u_k(j - 12) = Dynamics_RHS(j);
    }
  }
}
dlib::matrix<double> StateNAccNTorque2ContactForce(dlib::matrix<double> &State_k, dlib::matrix<double> &Acc_k, dlib::matrix<double> &u_k)
{
  // This function is used to initialize the contact force
  Robot_StateNDot Robot_StateNDot_k = DlibRobotstate2StateNDot(State_k);

  dlib::matrix<double> D_q_k, B_q_k, C_q_qdot_k, Jac_Full_k;

  Dynamics_Matrices(Robot_StateNDot_k, D_q_k, B_q_k, C_q_qdot_k, Jac_Full_k);

  dlib::matrix<double> B_LHS = D_q_k * Acc_k + C_q_qdot_k - B_q_k * u_k;

  dlib::matrix<double> Lamda_k = dlib::pinv(Jac_Full_k) * B_LHS;

  return Lamda_k;
}

void Seed_Config_Manifold_Pr_ObjNConstraint(std::vector<double> &Opt_Seed, std::vector<double> &ObjNConstraint_Val, std::vector<double> &ObjNConstraint_Type)
{
  // This function is used to formulate the constraint for the optimization problem

  dlib::matrix<int> End_Effector_Obs_Index_Matrix = dlib::ones_matrix<int>(6, Grids);

  dlib::matrix<double> temp_matrix;       std::vector<double> Seed_Config_on_Manifold;

  for (int i = 0; i < 13; i++)
  {
    Seed_Config_on_Manifold.push_back(Opt_Seed[i]);
  }
  // Pote_Robotstate is a global variable

  // Pure feasibility test with the requirement to minimize the difference between the Pote and Seed

  std::vector<double> Robostate_offset = Vec_Minus(Pote_Robotstate, Seed_Config_on_Manifold);
  double State_Diff = 0;
  for (int i = 0; i < Robostate_offset.size(); i++)
  {
    State_Diff = State_Diff + Robostate_offset[i] * Robostate_offset[i];
  }

  ObjNConstraint_Val.push_back(State_Diff);
  ObjNConstraint_Type.push_back(1);

  // Second is to make sure that Seed_Config_on_Manifold is on the manifold

  Robot_StateNDot Robot_StateNDot_Manifold = StateVec2StateNDot(Seed_Config_on_Manifold);
  dlib::matrix<double,12,1> End_Effector_Pos, End_Effector_Vel;
  dlib::matrix<double,6,1> End_Effector_Dist; std::vector<int> End_Effector_Obs(6);

  End_Effector_PosNVel(Robot_StateNDot_Manifold, End_Effector_Pos, End_Effector_Vel);
  End_Effector_Obs_Dist_Fn(End_Effector_Pos, End_Effector_Dist, End_Effector_Obs);

  // Here the End_Effector_Dist should multiply with the active sigma matrix
  for (int i = 0; i < 6; i++)
  {
    ObjNConstraint_Val.push_back(Sigma_Act_Manifold[i] * End_Effector_Dist(i));
    ObjNConstraint_Type.push_back(0);
  }
  return;
}

int Seed_Config_Manifold_Pr_(integer    *Status, integer *n,    doublereal x[],
     integer    *needF,  integer *neF,  doublereal F[],
     integer    *needG,  integer *neG,  doublereal G[],
     char       *cu,     integer *lencu,
	 integer    iu[],    integer *leniu,
	 doublereal ru[],    integer *lenru )
{
  std::vector<double> Opt_Seed, ObjNConstraint_Val, ObjNConstraint_Type;

  for (int i = 0; i < *n; i++)
  {
    Opt_Seed.push_back(x[i]);
    // cout<<x[i]<<endl;
  }

  Seed_Config_Manifold_Pr_ObjNConstraint(Opt_Seed, ObjNConstraint_Val, ObjNConstraint_Type);

  for (int i = 0; i < ObjNConstraint_Val.size(); i++)
  {
    F[i] = ObjNConstraint_Val[i];
  }
    return 0;
  }

void State_Traj_Interpolater_Manifold(std::vector<double>&Init_Robotstate, std::vector<double> &sigma, dlib::matrix<double> &StateNDot_Traj, dlib::matrix<double> &Acc_Traj, dlib::matrix<double> &StateNdot_Mid_Traj, dlib::matrix<double> &Acc_Mid_Traj)
{
  // This method is based on a parabola assumption of the state spline while the initial point and the end point should lie on the manifold on pos and velocity

  double Pos_k_Init, Pos_k_End, Vel_k_Init, Acc_k;             // The Position, Velocity and Acceleration for variable k

  // One thing needs to be noticed is that in this case, the final state does not need to be known in advance since we will directly make use of the T_tot to compute the final state
  // Here in this initialization of optimization seed, the acceleration is determined by the h_k inst4ad of the whole h_k

  for (int i = 0; i < 13; i++)
  {
    Pos_k_Init = Init_Robotstate[i];
    Vel_k_Init = Init_Robotstate[i+13];
    Acc_k = Vel_k_Init/h_k;                  // Since the assumption is that the acceleration remains constant.
    // Then it is the computation of the final robot configuration
    Pos_k_End = Pos_k_Init + Vel_k_Init * h_k + 1/2 * Acc_k * h_k * h_k;
    Pote_Robotstate[i] = Pos_k_End;           // The parabola state at the end
    Pote_Robotstate[i+13] = 0;                // The robot velocity at the end is set to 0 for the sake of kinetic energy minimization
  }

  // However, it is known in adavance that this Pote_Robotstate will not lie on the manifold so the idea is to enforce this contact constraint
  // This projection can be conducted by an orthogonal projection in the direction of the transpose of the Jacobian matrix

  // It is very likely that Pote_Robotstate is not on the manifold so the next job is to project this point back into the manifold.

  std::vector<double> Jac_Act_Index = Full_Row_Rank_Index(sigma);

  const int Jac_Row_Number = Jac_Act_Index.size();

  dlib::matrix<double> Jac_Act = dlib::zeros_matrix<double>(Jac_Row_Number,13);

  Robot_StateNDot Pote_RobotstateNDot = StateVec2StateNDot(Pote_Robotstate);

  dlib::matrix<double> Jac_Full = Jac_Full_fn(Pote_RobotstateNDot);

  for (int j = 0; j < Jac_Row_Number; j++)
  {
    for (int k = 0; k < 13; k++)
    {
      Jac_Act(j,k) = Jac_Full(Jac_Act_Index[j],k);
    }
  }

  Jac_Act_Trans_Manifold = dlib::trans(Jac_Act);

  // Here we have the descent direction! Then the job is to find the projection and the value. What we would like to do is to conduct an optimization to find the seed config

  snoptProblem Seed_Config_Manifold_Pr;                     // This is the name of the Optimization problem
  // Allocate and initialize
  std:vector<double> ObjNConstraint_Val, ObjNConstraint_Type;

  integer n = 13;

  // This is time to initialize the guess for this seed configuration manifold optimization problem
  std::vector<double> seed_config_opt_var(n);
  for (int i = 0; i < n; i++)
  {
    seed_config_opt_var[i] = Pote_Robotstate[i];
  }

  Seed_Config_Manifold_Pr_ObjNConstraint(seed_config_opt_var, ObjNConstraint_Val, ObjNConstraint_Type);

  integer neF = ObjNConstraint_Val.size();     			// 1 objective function
  integer lenA  =  n * neF;                         // This is the number of nonzero elements in the linear part A    F(x) = f(x)+Ax

  integer *iAfun = new integer[lenA];               integer *jAvar = new integer[lenA];				    doublereal *A  = new doublereal[lenA];

  integer lenG   = lenA;							              integer *iGfun = new integer[lenG];				    integer *jGvar = new integer[lenG];

  doublereal *x      = new doublereal[n];			      doublereal *xlow   = new doublereal[n];				doublereal *xupp   = new doublereal[n];
  doublereal *xmul   = new doublereal[n];			      integer    *xstate = new    integer[n];

  doublereal *F      = new doublereal[neF];		      doublereal *Flow   = new doublereal[neF];			doublereal *Fupp   = new doublereal[neF];
  doublereal *Fmul   = new doublereal[neF];		      integer    *Fstate = new integer[neF];

  integer nxnames = 1;							integer nFnames = 1;						char *xnames = new char[nxnames*8];						char *Fnames = new char[nFnames*8];

  integer    ObjRow = 0;							doublereal ObjAdd = 0;

  // Set the upper and lower bounds.
  for (int i = 0; i < n; i++) {
    xlow[i] = xlow_vec(i);
    xupp[i] = xupp_vec(i);
    xstate[i] = 0.0;
    x[i] = seed_config_opt_var[i];  	// Initial guess
  }

  for(int i = 0; i<neF; i++)
  {
    // The lower bound is the same
    Flow[i] = 0.0;
    if(ObjNConstraint_Type[i]>0)	// Inequality constraint
    {
      Fupp[i] = Inf;
    }
    else
    {
      Fupp[i] = 0.0;
    }
  }

  // Load the data for ToyProb ...
  Seed_Config_Manifold_Pr.setPrintFile  ( "Seed_Config_Manifold_Pr.out" );
  Seed_Config_Manifold_Pr.setProblemSize( n, neF );
  Seed_Config_Manifold_Pr.setObjective  ( ObjRow, ObjAdd );
  Seed_Config_Manifold_Pr.setA          ( lenA, iAfun, jAvar, A );
  Seed_Config_Manifold_Pr.setG          ( lenG, iGfun, jGvar );
  Seed_Config_Manifold_Pr.setX          ( x, xlow, xupp, xmul, xstate );
  Seed_Config_Manifold_Pr.setF          ( F, Flow, Fupp, Fmul, Fstate );
  Seed_Config_Manifold_Pr.setXNames     ( xnames, nxnames );
  Seed_Config_Manifold_Pr.setFNames     ( Fnames, nFnames );
  Seed_Config_Manifold_Pr.setProbName   ( "Seed_Config_Manifold_Pr" );
  Seed_Config_Manifold_Pr.setUserFun    ( Seed_Config_Manifold_Pr_);
  // snopta will compute the Jacobian by finite-differences.
  // The user has the option of calling  snJac  to define the
  // coordinate arrays (iAfun,jAvar,A) and (iGfun, jGvar).
  Seed_Config_Manifold_Pr.computeJac    ();
  Seed_Config_Manifold_Pr.setIntParameter( "Derivative option", 0 );
  Seed_Config_Manifold_Pr.setIntParameter( "Major print level", 0 );
  Seed_Config_Manifold_Pr.setIntParameter( "Minor print level", 0 );
  integer Cold = 0, Basis = 1, Warm = 2;
  Seed_Config_Manifold_Pr.solve( Cold );

  // Take the value out from x
  for (int i = 0; i < 13; i++)
  {
    Pote_Robotstate[i] = x[i];
  }

  // cin.get();
  Robot_StateNDot Init_Opt_vec(Pote_Robotstate);
  Robot_Plot_fn(Init_Opt_vec);

  delete []iAfun;  delete []jAvar;  delete []A;
  delete []iGfun;  delete []jGvar;

  delete []x;      delete []xlow;   delete []xupp;
  delete []xmul;   delete []xstate;

  delete []F;		   delete []Flow;	  delete []Fupp;
  delete []Fmul;	 delete []Fstate;

  delete []xnames; delete []Fnames;

  // Here the seed_configuration have been achieved to lie on the constraint manifold

  std::vector<double> State_Traj_Coeff_i;

  double Pos_i_Init, Pos_i_Goal, Vel_i_Init;

  for (int i = 0; i < 13; i++)
  {
    Pos_i_Init = Init_Robotstate[i];				Pos_i_Goal = Pote_Robotstate[i];         Vel_i_Init = Init_Robotstate[i+13];

    State_Traj_Coeff_i = QuadraticSpline_Coeff_fn(h_k, Pos_i_Init, Pos_i_Goal, Vel_i_Init); // 3 by 1 vector: a, b, c

    // After this step, we hacve the cubic spline for the whole trajectories. Then it is time to discretize them
    // for (int qq = 0; qq < 3; qq++)
    // {
    //   cout<<State_Traj_Coeff_i[qq]<<endl;
    // }

    double ds = 1.0/(Grids * 1.0 - 1.0), s_j;

    for (int j = 0; j < Grids; j++)
    {
      s_j = ds * (j);

      StateNDot_Traj(i, j) = QuadraticSpline_Eval(h_k, State_Traj_Coeff_i, s_j, 'P');

      StateNDot_Traj(i + 13,j) = QuadraticSpline_Eval(h_k, State_Traj_Coeff_i, s_j, 'V');

      Acc_Traj(i, j) = QuadraticSpline_Eval(h_k, State_Traj_Coeff_i, s_j, 'A');

      // Now it is the information at the middle point
      if(j < Grids-1)
      {
        double s_j_mid = s_j + 0.5 * ds;

        StateNdot_Mid_Traj(i, j) = QuadraticSpline_Eval(h_k, State_Traj_Coeff_i, s_j_mid, 'P');

        StateNdot_Mid_Traj(i, j+13) = QuadraticSpline_Eval(h_k, State_Traj_Coeff_i, s_j_mid, 'V');

        Acc_Mid_Traj(i, j) = QuadraticSpline_Eval(h_k, State_Traj_Coeff_i, s_j_mid, 'A');
      }
    }
  }

  cout<<StateNDot_Traj<<endl;
  cout<<Acc_Traj<<endl;


}

void State_Traj_Interpolater(std::vector<double>&Init_Robotstate, std::vector<double>&End_Robotstate, dlib::matrix<double> &StateNDot_Traj, dlib::matrix<double> &Acc_Traj, dlib::matrix<double> &StateNdot_Mid_Traj, dlib::matrix<double> &Acc_Mid_Traj, int Choice_Flag)
{
  // This function is used to interpolate the State Trajectory and get ready the acceleration

  // Based on the different methods, there are three possible ways to address these problem

  if(Choice_Flag == 1)
  {
    // Method 1: The whole states are chosen in a linear fashion. Then the acceleration is found by treating each segment between grids to be cubic spline.

    dlib::matrix<double> Robot_State_Interpol_i;

    for (int i = 0; i < 26; i++)
    {
      Robot_State_Interpol_i = dlib::linspace(Init_Robotstate[i], End_Robotstate[i], Grids);

      for (int j = 0; j < Grids; j++)
      {
        StateNDot_Traj(i,j) = Robot_State_Interpol_i(j);
      }
    }
    CubicSpline_Interpolater(StateNDot_Traj, Acc_Traj, StateNdot_Mid_Traj, Acc_Mid_Traj);
  }
  else
  {
    if(Choice_Flag == 2)
    {
      // Method 2: Each state trajectory will be an individual cubic spline

      double T_tot = h_k * (Grids - 1);

      std::vector<double> State_Traj_Coeff_i;

      double Pos_i_Init, Pos_i_Goal, Vel_i_Init, Vel_i_Goal;

      for (int i = 0; i < 13; i++)
      {
        Pos_i_Init = Init_Robotstate[i];				Pos_i_Goal = End_Robotstate[i];

        Vel_i_Init = Init_Robotstate[i+13];			Vel_i_Goal = End_Robotstate[i+13];

        State_Traj_Coeff_i = CubicSpline_Coeff_fn(T_tot, Pos_i_Init, Pos_i_Goal, Vel_i_Init, Vel_i_Goal); // 4 by 1 vector: a, b, c, d

        // After this step, we hacve the cubic spline for the whole trajectories. Then it is time to discretize them

        double ds = 1.0/(Grids * 1.0 - 1.0), s_j;

        for (int j = 0; j < Grids; j++)
        {
          s_j = ds * (j);

          StateNDot_Traj(i, j) = CubicSpline_Eval(T_tot, State_Traj_Coeff_i, s_j, 'P');

          StateNDot_Traj(i + 13,j) = CubicSpline_Eval(T_tot, State_Traj_Coeff_i, s_j, 'V');

          Acc_Traj(i, j) = CubicSpline_Eval(T_tot, State_Traj_Coeff_i, s_j, 'A');

          // Now it is the information at the middle point
          if(j < Grids-1)
          {
            double s_j_mid = s_j + 0.5 * ds;

            StateNdot_Mid_Traj(i, j) = CubicSpline_Eval(T_tot, State_Traj_Coeff_i, s_j_mid, 'P');

            StateNdot_Mid_Traj(i, j+13) = CubicSpline_Eval(T_tot, State_Traj_Coeff_i, s_j_mid, 'V');

            Acc_Mid_Traj(i, j) = CubicSpline_Eval(T_tot, State_Traj_Coeff_i, s_j_mid, 'A');

          }
        }
      }
    }
    else
    {
      // Method 3: Each state trajectory is a parabola where the acceleration is set to be a constant value during the whole duration

      double T_tot = h_k * (Grids - 1);

      std::vector<double> State_Traj_Coeff_i;

      double Pos_i_Init, Pos_i_Goal, Vel_i_Init;

      for (int i = 0; i < 13; i++)
      {
        Pos_i_Init = Init_Robotstate[i];				Pos_i_Goal = End_Robotstate[i];         Vel_i_Init = Init_Robotstate[i+13];

        State_Traj_Coeff_i = QuadraticSpline_Coeff_fn(T_tot, Pos_i_Init, Pos_i_Goal, Vel_i_Init); // 4 by 1 vector: a, b, c, d

        // After this step, we hacve the cubic spline for the whole trajectories. Then it is time to discretize them

        double ds = 1.0/(Grids * 1.0 - 1.0), s_j;

        for (int j = 0; j < Grids; j++)
        {
          s_j = ds * (j);

          StateNDot_Traj(i, j) = QuadraticSpline_Eval(T_tot, State_Traj_Coeff_i, s_j, 'P');

          StateNDot_Traj(i + 13,j) = QuadraticSpline_Eval(T_tot, State_Traj_Coeff_i, s_j, 'V');

          Acc_Traj(i, j) = QuadraticSpline_Eval(T_tot, State_Traj_Coeff_i, s_j, 'A');

          // Now it is the information at the middle point
          if(j < Grids-1)
          {
            double s_j_mid = s_j + 0.5 * ds;

            StateNdot_Mid_Traj(i, j) = QuadraticSpline_Eval(T_tot, State_Traj_Coeff_i, s_j_mid, 'P');

            StateNdot_Mid_Traj(i, j+13) = QuadraticSpline_Eval(T_tot, State_Traj_Coeff_i, s_j_mid, 'V');

            Acc_Mid_Traj(i, j) = QuadraticSpline_Eval(T_tot, State_Traj_Coeff_i, s_j_mid, 'A');
          }
        }
      }
    }
  }
}

double QuadraticSpline_Eval(double T, std::vector<double> & State_Traj_Coeff, double s, char name)
{
  // This function is used to calcualte the position, velocity and acceleration given the spline coefficients
  // Here T is the duration constant, s is the path variable
  // y(s) = a*s^2 + b * s + c
  double a = State_Traj_Coeff[0];
  double b = State_Traj_Coeff[1];
  double c = State_Traj_Coeff[2];

  if (name == 'P')
  {
    double Pos = a*s*s + b*s + c;
    return Pos;
  }
  else
  {
    if (name == 'V')
    {
      double Vel = (2 * a * s + b)/T;
      return Vel;
    }
    else
    {
      double Acc = (2 * a)/(T*T);
      return Acc;
    }
  }
}

void CubicSpline_Interpolater(dlib::matrix<double> &StateNDot_Traj, dlib::matrix<double> &Acc_Traj, dlib::matrix<double> &StateNdot_Mid_Traj, dlib::matrix<double> &Acc_Mid_Traj)
{
  // This function is used to interpolate the given state trajectory in a cubic spline manner
  // Which means that each sequential grid points

  double x_init, q_mid, qdot_mid, x_end, xdot_init, xdot_end, xddot_init, qddot_mid, xddot_end;

  std::vector<double> CubicSpline_Coeff_Pos_vec, CubicSpline_Coeff_Vel_vec;

  for (int i = 0; i < Grids-1; i++)
  {
    for (int j = 0; j < 13; j++)
    {
      // Position first
      x_init = StateNDot_Traj(j,i);
      x_end = StateNDot_Traj(j,i+1);
      xdot_init = StateNDot_Traj(j + 13, i);
      xdot_end = StateNDot_Traj(j + 13, i+1);
      CubicSpline_Coeff_Pos_vec = CubicSpline_Coeff_fn(h_k, x_init, x_end, xdot_init, xdot_end);        // 4 by 1 vector: a, b, c, d

      // Velocity second
      xddot_init = CubicSpline_Eval(h_k, CubicSpline_Coeff_Pos_vec, 0, 'A');
      xddot_end = CubicSpline_Eval(h_k, CubicSpline_Coeff_Pos_vec, 1, 'A');

      Acc_Traj(j,i) = xddot_init;
      Acc_Traj(j,i+1) = xddot_end;

      // Now it is the middle point state and acceleration

      // Position and Velocity

      q_mid = CubicSpline_Eval(h_k, CubicSpline_Coeff_Pos_vec, 0.5, 'P');
      qdot_mid = CubicSpline_Eval(h_k, CubicSpline_Coeff_Pos_vec, 0.5, 'V');
      qddot_mid = CubicSpline_Eval(h_k, CubicSpline_Coeff_Pos_vec, 0.5, 'A');

      StateNdot_Mid_Traj(j,i) = q_mid;
      StateNdot_Mid_Traj(j,i + 13) = qdot_mid;
      Acc_Mid_Traj(j,i) = qddot_mid;
    }
  }
}

void Opt_Seed_Unzip(std::vector<double> &Opt_Seed, double &T_tot, dlib::matrix<double> & StateNDot_Traj, dlib::matrix<double> & Ctrl_Traj, dlib::matrix<double> & Contact_Force_Traj)
{
  const int NumOfStateNDot_Traj = 26;
  const int NumOfCtrl_Traj = 10;
  const int NumOfContactForce_Traj = 12;

  T_tot = Opt_Seed[0];			int Opt_Seed_Index = 1;

  StateNDot_Traj = dlib::zeros_matrix<double>(NumOfStateNDot_Traj, Grids);
  Ctrl_Traj = dlib::zeros_matrix<double>(NumOfCtrl_Traj, Grids);
  Contact_Force_Traj = dlib::zeros_matrix<double>(NumOfContactForce_Traj, Grids);

  // 1. Retrieve the StateNDot_Traj matrix
  for (int i = 0; i < Grids; i++)
  {
    for (int j = 0; j < NumOfStateNDot_Traj; j++)
    {
      StateNDot_Traj(j,i) = Opt_Seed[Opt_Seed_Index];
      Opt_Seed_Index = Opt_Seed_Index + 1;
    }
  }
  // cout<<StateNDot_Traj<<endl;
  // 2. Retrieve the control matrix
  for (int i = 0; i < Grids; i++)
  {
    for (int j = 0; j < NumOfCtrl_Traj; j++)
    {
      Ctrl_Traj(j,i) = Opt_Seed[Opt_Seed_Index];
      Opt_Seed_Index = Opt_Seed_Index + 1;
    }
  }
  // cout<<Ctrl_Traj<<endl;
  // 3. Retrieve the contact force matrix
  for (int i = 0; i < Grids; i++)
  {
    for (int j = 0; j < NumOfContactForce_Traj; j++)
    {
      Contact_Force_Traj(j,i) = Opt_Seed[Opt_Seed_Index];
      Opt_Seed_Index = Opt_Seed_Index + 1;
    }
  }
  // cout<<Contact_Force_Traj<<endl;
  return ;
}

void Opt_Seed_Zip(std::vector<double> &Opt_Seed, dlib::matrix<double> & StateNDot_Traj, dlib::matrix<double> & Ctrl_Traj, dlib::matrix<double> & Contact_Force_Traj)
{
  // This function is used to stack the coefficient matrices into a column vector
  for (int i = 0; i < StateNDot_Traj.nc(); i++)
  {
    for (int j = 0; j < StateNDot_Traj.nr(); j++)
    {
      Opt_Seed.push_back(StateNDot_Traj(j,i));
    }
  }
  for (int i = 0; i < Ctrl_Traj.nc(); i++)
  {
    for (int j = 0; j < Ctrl_Traj.nr(); j++)
    {
      Opt_Seed.push_back(Ctrl_Traj(j,i));
    }
  }
  for (int i = 0; i < Contact_Force_Traj.nc(); i++)
  {
    for (int j = 0; j < Contact_Force_Traj.nr(); j++)
    {
      Opt_Seed.push_back(Contact_Force_Traj(j,i));
    }
  }
}

std::vector<double> Seed_Guess_Gene_Robotstate(Tree_Node &Node_i, Tree_Node &Node_i_child)
{
  // This function is used to generate a configuration to initialize the optimization
  std::vector<double> Robot_State_Seed = StateNDot2StateVec(Node_i.Node_StateNDot);
  std::vector<double> ObjNConstraint_Val, ObjNConstraint_Type;
  Seed_Conf_Optimization_ObjNConstraint(Robot_State_Seed, ObjNConstraint_Val, ObjNConstraint_Type);
  snoptProblem Seed_Conf_Optimization_Pr;                     // This is the name of the Optimization problem for the robot configuration

  integer n = Robot_State_Seed.size();
  integer neF = ObjNConstraint_Val.size();
  integer lenA  =  n * neF;

  integer *iAfun = new integer[lenA];              integer *jAvar = new integer[lenA];					doublereal *A  = new doublereal[lenA];

  integer lenG   = lenA;							integer *iGfun = new integer[lenG];						integer *jGvar = new integer[lenG];

  doublereal *x      = new doublereal[n];			doublereal *xlow   = new doublereal[n];					doublereal *xupp   = new doublereal[n];
  doublereal *xmul   = new doublereal[n];			integer    *xstate = new    integer[n];

  doublereal *F      = new doublereal[neF];		doublereal *Flow   = new doublereal[neF];				doublereal *Fupp   = new doublereal[neF];
  doublereal *Fmul   = new doublereal[neF];		integer    *Fstate = new integer[neF];

  integer nxnames = 1;							integer nFnames = 1;				char *xnames = new char[nxnames*8];					char *Fnames = new char[nFnames*8];

  integer    ObjRow = 0;							doublereal ObjAdd = 0;

  for (int i = 0; i < n; i++) {
    xlow[i] = xlow_vec(i);						xupp[i] = xupp_vec(i);				xstate[i] = 0.0;							x[i] = Robot_State_Seed[i];  	// Initial guess
  }

  for(int i = 0; i<neF; i++){
    // The lower bound is the same
    Flow[i] = 0.0;
    if(ObjNConstraint_Type[i]>0)	// Inequality constraint
    {	Fupp[i] = Inf;}
    else{
      Fupp[i] = 0.0;}
    }

    // Load the data for ToyProb ...
    Seed_Conf_Optimization_Pr.setPrintFile  ( "Seed_Conf_Optimization_Pr.out" );
    Seed_Conf_Optimization_Pr.setProblemSize( n, neF );
    Seed_Conf_Optimization_Pr.setObjective  ( ObjRow, ObjAdd );
    Seed_Conf_Optimization_Pr.setA          ( lenA, iAfun, jAvar, A );
    Seed_Conf_Optimization_Pr.setG          ( lenG, iGfun, jGvar );
    Seed_Conf_Optimization_Pr.setX          ( x, xlow, xupp, xmul, xstate );
    Seed_Conf_Optimization_Pr.setF          ( F, Flow, Fupp, Fmul, Fstate );
    Seed_Conf_Optimization_Pr.setXNames     ( xnames, nxnames );
    Seed_Conf_Optimization_Pr.setFNames     ( Fnames, nFnames );
    Seed_Conf_Optimization_Pr.setProbName   ( "Seed_Conf_Optimization_Pr" );
    Seed_Conf_Optimization_Pr.setUserFun    ( Seed_Conf_Optimization_Pr_fn_);
    // snopta will compute the Jacobian by finite-differences.
    // The user has the option of calling  snJac  to define the
    // coordinate arrays (iAfun,jAvar,A) and (iGfun, jGvar).
    Seed_Conf_Optimization_Pr.computeJac    ();
    Seed_Conf_Optimization_Pr.setIntParameter( "Derivative option", 0 );
    Seed_Conf_Optimization_Pr.setIntParameter( "Print level", 0 );
    Seed_Conf_Optimization_Pr.setIntParameter( "Major print level", 0 );
    Seed_Conf_Optimization_Pr.setIntParameter( "Minor print level", 0 );

    integer Cold = 0, Basis = 1, Warm = 2;
    Seed_Conf_Optimization_Pr.solve( Cold );
    for (int i = 0; i < 26; i++){Robot_State_Seed[i] = x[i];}
    Robot_StateNDot Robot_StateNDot_Seed(Robot_State_Seed);
    std::string input_name = "Seed Configuration in Opt";
    Robot_Plot_fn(Robot_StateNDot_Seed, input_name);
    delete []iAfun;  delete []jAvar;  delete []A;		delete []iGfun;  delete []jGvar;
    delete []x;      delete []xlow;   delete []xupp;	delete []xmul;   delete []xstate;
    delete []F;      delete []Flow;   delete []Fupp;	delete []Fmul;   delete []Fstate;
    delete []xnames; delete []Fnames;
    return Robot_State_Seed;
}

int Seed_Conf_Optimization_Pr_fn_(integer    *Status, integer *n,    doublereal x[],
     integer    *needF,  integer *neF,  doublereal F[],
     integer    *needG,  integer *neG,  doublereal G[],
     char       *cu,     integer *lencu,
	 integer    iu[],    integer *leniu,
	 doublereal ru[],    integer *lenru )
{
  std::vector<double> Opt_Seed, ObjNConstraint_Val, ObjNConstraint_Type;
  for (int i = 0; i < 26; i++)
  {
    Opt_Seed.push_back(x[i]);
  }
  Seed_Conf_Optimization_ObjNConstraint(Opt_Seed, ObjNConstraint_Val, ObjNConstraint_Type);
  for (int i = 0; i < ObjNConstraint_Val.size(); i++)
  {
    F[i] = ObjNConstraint_Val[i];
  }
    return 0;
  }

void Seed_Conf_Optimization_ObjNConstraint(std::vector<double> &Opt_Seed, std::vector<double> &ObjNConstraint_Val, std::vector<double> &ObjNConstraint_Type)
{
  Robot_StateNDot StateNDot_Init_i(Opt_Seed);		dlib::matrix<double,12,1> End_Effector_Pos, End_Effector_Vel, End_Effector_Pos_ref, End_Effector_Vel_ref;

  End_Effector_PosNVel(StateNDot_Init_i, End_Effector_Pos, End_Effector_Vel);

  End_Effector_Pos_ref = Structure_P.Node_i.End_Effector_Pos;                              End_Effector_Vel_ref = Structure_P.Node_i.End_Effector_Vel;

  Robot_StateNDot StateNDot_Init_ref = Structure_P.Node_i.Node_StateNDot;

  std::vector<double> StateVec_Init_ref = StateNDot2StateVec(StateNDot_Init_ref);

  std::vector<double> sigma_i = Structure_P.Node_i.sigma;                                  std::vector<double> sigma_i_child = Structure_P.Node_i_child.sigma;

  std::vector<double> rCOM_ref = Ang_Pos_fn(Structure_P.Node_i.Node_StateNDot, "rCOM");
  std::vector<double> rA_ref = Ang_Pos_fn(Structure_P.Node_i.Node_StateNDot, "rA");
  std::vector<double> rB_ref = Ang_Pos_fn(Structure_P.Node_i.Node_StateNDot, "rB");
  std::vector<double> rC_ref = Ang_Pos_fn(Structure_P.Node_i.Node_StateNDot, "rC");
  std::vector<double> rD_ref = Ang_Pos_fn(Structure_P.Node_i.Node_StateNDot, "rD");
  std::vector<double> rE_ref = Ang_Pos_fn(Structure_P.Node_i.Node_StateNDot, "rE");
  std::vector<double> rF_ref = Ang_Pos_fn(Structure_P.Node_i.Node_StateNDot, "rF");
  std::vector<double> rT_ref = Ang_Pos_fn(Structure_P.Node_i.Node_StateNDot, "rT");

  std::vector<double> rCOM_opt = Ang_Pos_fn(StateNDot_Init_i, "rCOM");
  std::vector<double> rA_opt = Ang_Pos_fn(StateNDot_Init_i, "rA");
  std::vector<double> rB_opt = Ang_Pos_fn(StateNDot_Init_i, "rB");
  std::vector<double> rC_opt = Ang_Pos_fn(StateNDot_Init_i, "rC");
  std::vector<double> rD_opt = Ang_Pos_fn(StateNDot_Init_i, "rD");
  std::vector<double> rE_opt = Ang_Pos_fn(StateNDot_Init_i, "rE");
  std::vector<double> rF_opt = Ang_Pos_fn(StateNDot_Init_i, "rF");
  std::vector<double> rT_opt = Ang_Pos_fn(StateNDot_Init_i, "rT");

  std::vector<double> vCOM_ref = Ang_Vel_fn(StateNDot_Init_i,"vCOM");

  ObjNConstraint_Val.push_back(0);     ObjNConstraint_Type.push_back(1);

  int Opt_Type_Flag = Sigma_TransNGoal(sigma_i, sigma_i_child);

  double sigma_diff = 0.0; double Obj_val = 0.0;

  double mini = 0.025;

  for (int i = 0; i < 4; i++)
  {
    sigma_diff = sigma_diff + (sigma_i_child[i] - sigma_i[i]);
  }
  if(abs(sigma_diff)>0)
  {
    // In this case, the robot changes its contact point.

    Obj_val = 0.0;

    // In this case, there must a contact modification during this whole process. As a result, the cost function is the quadratic sum of the state change.

    for (int i = 0; i < StateVec_Init_ref.size(); i++)
    {
      Obj_val = Obj_val + (StateVec_Init_ref[i] - Opt_Seed[i]) * (StateVec_Init_ref[i] - Opt_Seed[i]);
    }

    if(sigma_diff<0)
    {
      // Contact break/retract this part is a little tough to handle. However, in our problem, we only consider the foot contact retract

      int Contact_Index = Sigma_Change(sigma_i, sigma_i_child);

      if(Contact_Index<2)
      {
        // In this case, the robot will lift up left/right foot
        if(Contact_Index == 0)
        {
          // In this case, the robot will lift AB foot so we would like to move the COM to CD
          ObjNConstraint_Val.push_back(rCOM_opt[0] - rD_ref[0] - mini);
          ObjNConstraint_Type.push_back(1);
          ObjNConstraint_Val.push_back(rC_ref[0] - rCOM_opt[0] - mini);
          ObjNConstraint_Type.push_back(1);
        }
        else
        {
          // In this case, the robot will lift CD foot so we would like to move the COM to AB
          ObjNConstraint_Val.push_back(rCOM_opt[0] - rB_ref[0] - mini);
          ObjNConstraint_Type.push_back(1);
          ObjNConstraint_Val.push_back(rA_ref[0] - rCOM_opt[0] - mini);
          ObjNConstraint_Type.push_back(1);
        }
      }
    }
    else
    {
      // Making contact
      int Contact_Index = Sigma_Change(sigma_i, sigma_i_child);

      // The next step is to figure out which foot to step and which direction to go

      double vCOM_sign = vCOM_ref[0]/abs(vCOM_ref[0]);

      if(Contact_Index<2)
      {
        // Making contact with foot
        if(Contact_Index > 0)
        {
          // Making contact with the right foot
          if(vCOM_sign>0)
          {
            // Robot moves forward
            ObjNConstraint_Val.push_back(rD_opt[0] - rA_ref[0] - mini);
            ObjNConstraint_Type.push_back(0);
          }
          else
          {
            // Robot moves backward
            ObjNConstraint_Val.push_back(rB_ref[0] - rC_ref[0] - mini);
            ObjNConstraint_Type.push_back(1);
          }
        }
        else
        {
          // Making contact with the left foot
          if(vCOM_sign>0)
          {
            ObjNConstraint_Val.push_back(rB_opt[0] - rC_ref[0] - mini);
            ObjNConstraint_Type.push_back(1);
          }
          else
          {
            ObjNConstraint_Val.push_back(rD_ref[0] - rA_opt[0] - mini);
            ObjNConstraint_Type.push_back(1);
          }
          // ObjNConstraint_Val.push_back((rT_opt[0] - rT_ref[0]) * (rT_opt[0] - rT_ref[0]));
          // ObjNConstraint_Type.push_back(0);
        }
      }
    }
  }
  else
  {
    // This is used for the self-stasbilization process so kinetic energy is the only cost
    Obj_val = Kinetic_Energy_fn(StateNDot_Init_i);
  }

  ObjNConstraint_Val[0] = Obj_val;

  // Here some heuristic inequality constraint have been added. Now the job is to figure out the equality constraint

  // double Hand_max = max(rE_ref[0], rF_ref[0]);
  // Hand_max = max(Hand_max, rD_ref[0]);
  // Hand_max = max(Hand_max, rC_ref[0]);
  // Hand_max = max(Hand_max, rA_ref[0]);
  // Hand_max = max(Hand_max, rT_opt[0]);

  dlib::matrix<double,6,1> End_Effector_Dist;

  std::vector<int> End_Effector_Obs(6);

  End_Effector_Obs_Dist_Fn(End_Effector_Pos, End_Effector_Dist, End_Effector_Obs);

  dlib::matrix<double> Eqn_Pos_Matrix, Ineqn_Pos_Matrix, Eqn_Vel_Matrix, Eqn_Maint_Matrix, Matrix_result;

  std::vector<double> sigma_temp, sigma_real;

  if(Opt_Type_Flag == -1)
  {
    // In this case, sigma_i has to be maintained.
    sigma_real = sigma_i;
  }
  else
  {
    // In other cases, sigma_i_child has to be maintained.
    sigma_real = sigma_i_child;
  }

  sigma_temp = Sigma2Pos(sigma_real, 0);			Eqn_Pos_Matrix = Diag_Matrix_fn(sigma_temp);
  sigma_temp = Sigma2Pos(sigma_real, 1);			Ineqn_Pos_Matrix = Diag_Matrix_fn(sigma_temp);
  sigma_temp = Sigma2Vel(sigma_real);				  Eqn_Vel_Matrix = Diag_Matrix_fn(sigma_temp);

  // 1. Active constraints have to be satisfied: Position and Velocity
  Matrix_result = Eqn_Pos_Matrix * End_Effector_Dist;
  ObjNConstraint_ValNType_Update(Matrix_result, ObjNConstraint_Val, ObjNConstraint_Type, 0);

  Matrix_result = Eqn_Vel_Matrix * End_Effector_Vel;
  ObjNConstraint_ValNType_Update(Matrix_result, ObjNConstraint_Val, ObjNConstraint_Type, 0);

  // 2. Inactive constraints have to be strictly away from the obstacle
  dlib::matrix<double> ones_vector, temp_matrix;
  ones_vector = ONES_VECTOR_fn(6);
  Matrix_result = Ineqn_Pos_Matrix * (End_Effector_Dist - ones_vector * mini);
  ObjNConstraint_ValNType_Update(Matrix_result, ObjNConstraint_Val, ObjNConstraint_Type, 1);

  // 3. Middle joints have to be strictly away from the obs
  temp_matrix = Middle_Joint_Obs_Dist_Fn(StateNDot_Init_i);
  Matrix_result = temp_matrix - ones_vector * mini;
  ObjNConstraint_ValNType_Update(Matrix_result, ObjNConstraint_Val, ObjNConstraint_Type, 1);

  // 4. One more constraint to be added is to maintain the active unchanged constraint
  if(Opt_Type_Flag == -1)
  {
    Eqn_Maint_Matrix = Eqn_Maint_Matrix_fn(sigma_i, sigma_i);
  }
  else
  {
    Eqn_Maint_Matrix = Eqn_Maint_Matrix_fn(sigma_i, sigma_i_child);
  }
  Matrix_result = Eqn_Maint_Matrix * (End_Effector_Pos_ref - End_Effector_Pos);
  ObjNConstraint_ValNType_Update(Matrix_result, ObjNConstraint_Val, ObjNConstraint_Type, 0);
  return;
}

std::vector<double> Time_Seed_Queue_fn(double Time_Interval, int Total_Num)
{
  // This function is used to generate the Time_Seed_Queue
  double Time_Center = 0.5;
  int One_Side_Num = (Total_Num - 1)/2;
  std::vector<double> Time_Seed_Queue;
  Time_Seed_Queue.push_back(Time_Center);
  for (int i = 1; i < One_Side_Num; i++)
  {
    Time_Seed_Queue.push_back(Time_Center + Time_Interval * i);
    Time_Seed_Queue.push_back(Time_Center - Time_Interval * i);
  }
  // std::vector<double> Time_Seed_Queue;
  // Time_Seed_Queue.push_back(Time_Center);
  // for (int i = 0; i < Total_Num; i++) {
  // 	Time_Seed_Queue.push_back(Time_Center + (i+1) * Time_Interval);
  // }
  return Time_Seed_Queue;
}

std::vector<double> Sigma2Pos(std::vector<double> &sigma, int EqOrIneq)
{
  std::vector<double> sigma_pos;
  if(EqOrIneq==0)
  {
    sigma_pos.push_back(sigma[0]);
    sigma_pos.push_back(sigma[0]);
    sigma_pos.push_back(sigma[1]);
    sigma_pos.push_back(sigma[1]);
    sigma_pos.push_back(sigma[2]);
    sigma_pos.push_back(sigma[3]);
  }
  else
  {
    sigma_pos.push_back(!sigma[0]);
    sigma_pos.push_back(!sigma[0]);
    sigma_pos.push_back(!sigma[1]);
    sigma_pos.push_back(!sigma[1]);
    sigma_pos.push_back(!sigma[2]);
    sigma_pos.push_back(!sigma[3]);
  }
  return sigma_pos;
}

std::vector<double> Sigma2Vel(std::vector<double> &sigma)
{
  std::vector<double> sigma_pos;
  sigma_pos.push_back(sigma[0]);
  sigma_pos.push_back(sigma[0]);
  sigma_pos.push_back(sigma[0]);
  sigma_pos.push_back(sigma[0]);
  sigma_pos.push_back(sigma[1]);
  sigma_pos.push_back(sigma[1]);
  sigma_pos.push_back(sigma[1]);
  sigma_pos.push_back(sigma[1]);
  sigma_pos.push_back(sigma[2]);
  sigma_pos.push_back(sigma[2]);
  sigma_pos.push_back(sigma[3]);
  sigma_pos.push_back(sigma[3]);
  return sigma_pos;
}

int Sigma_TransNGoal(std::vector<double> & sigma_i, std::vector<double> & sigma_i_child)
{
  // The role of this function is actually very important because making contact is easier than retracting contact since retracting contact requires an accumulation of
  // the necessary kinetic energy of certain body parts before a retract can be made
  double sigma_result = 0.0;
  int Opt_Type_Flag = 0;
  for (int i = 0; i < sigma_i.size(); i++) {
    sigma_result = sigma_result + sigma_i_child[i] - sigma_i[i];
  }
  if(sigma_result>0)
  {
    // In this case, it is making the contact so the Critical frame is the ending frame
    Opt_Type_Flag = 1;
  }
  else
  {
    if (sigma_result ==0)
    {
      // In this case, it is the self-opt process
      Opt_Type_Flag = 0;
    }
    else
    {
      // In this case, it is the retracting contact process: the whole process is divided into two subprocesses, energy accumulation, energy release
      Opt_Type_Flag = -1;
    }
  }
  return Opt_Type_Flag;
}

Robot_StateNDot DlibRobotstate2StateNDot(dlib::matrix<double> &DlibRobotstate)
{
  // This function is used to convert the dlib matrix robot state to Robot_StateNDot type
  std::vector<double> Robot_StateNDot_vec;
  for (int i = 0; i < DlibRobotstate.nr(); i++)
  {
    Robot_StateNDot_vec.push_back(DlibRobotstate(i));
  }
  Robot_StateNDot Robot_StateNDot_i(Robot_StateNDot_vec);

  return Robot_StateNDot_i;
}

dlib::matrix<double> Eqn_Maint_Matrix_fn(std::vector<double> &sigma_i, std::vector<double> &sigma_i_child)
{
  // This function is used to generate the contact maintenance matrix
  std::vector<double> sigma_maint;
  sigma_maint.push_back(sigma_i[0] * sigma_i_child[0]);
  sigma_maint.push_back(sigma_i[0] * sigma_i_child[0]);
  sigma_maint.push_back(sigma_i[0] * sigma_i_child[0]);
  sigma_maint.push_back(sigma_i[0] * sigma_i_child[0]);
  sigma_maint.push_back(sigma_i[1] * sigma_i_child[1]);
  sigma_maint.push_back(sigma_i[1] * sigma_i_child[1]);
  sigma_maint.push_back(sigma_i[1] * sigma_i_child[1]);
  sigma_maint.push_back(sigma_i[1] * sigma_i_child[1]);
  sigma_maint.push_back(sigma_i[2] * sigma_i_child[2]);
  sigma_maint.push_back(sigma_i[2] * sigma_i_child[2]);
  sigma_maint.push_back(sigma_i[3] * sigma_i_child[3]);
  sigma_maint.push_back(sigma_i[3] * sigma_i_child[3]);

  dlib::matrix<double> Eqn_Maint_Matrix;
  Eqn_Maint_Matrix = Diag_Matrix_fn(sigma_maint);
  return Eqn_Maint_Matrix;
}

std::vector<double> QuadraticSpline_Coeff_fn(double T, double x_init, double x_end, double xdot_init)
{
  // This function is used to calcualte the coefficients for the parabola spline
  // The cubic spline is expressed to be : y(s) = a*s^2 + b*s + c
  std::vector<double> QuadraticSpline_Coeff_vec(3);     double a, b, c;
  b = xdot_init * T;        c = x_init;         a = x_end - b - c;

  QuadraticSpline_Coeff_vec[0] = a;
  QuadraticSpline_Coeff_vec[1] = b;
  QuadraticSpline_Coeff_vec[2] = c;

  return QuadraticSpline_Coeff_vec;
}

std::vector<double> CubicSpline_Coeff_fn(double T, double x_init, double x_end, double xdot_init, double xdot_end)
{
  // This function is used to calcualte the coefficients for the cubic spline
  // The cubic spline is expressed to be : y(s) = a*s^3 + b*s^2 + c*s + d
  std::vector<double> CubicSpline_Coeff_vec(4);     double a, b, c, d;
  a = 2 * x_init - 2 * x_end + T * xdot_end + T * xdot_init;
  b = 3 * x_end - 3 * x_init - T * xdot_end - 2 * T * xdot_init;
  c = T * xdot_init;
  d = x_init;

  CubicSpline_Coeff_vec[0] = a;               CubicSpline_Coeff_vec[1] = b;
  CubicSpline_Coeff_vec[2] = c;               CubicSpline_Coeff_vec[3] = d;

  return CubicSpline_Coeff_vec;
}

double CubicSpline_Eval(double T, std::vector<double> & CubicSpline_Coeff_vec, double s, char name)
{
  //# This function is used to calcualte the position, velocity and acceleration given the spline coefficients
  //# Here T is the duration constant, s is the path variable
  double a, b, c, d;
  a = CubicSpline_Coeff_vec[0];             b = CubicSpline_Coeff_vec[1];
  c = CubicSpline_Coeff_vec[2];             d = CubicSpline_Coeff_vec[3];

  if (name == 'P')
  {
    double Pos = a*s*s*s + b*s*s + c*s + d;
    return Pos;
  }
  else
  {
    if (name == 'V')
    {
      double Vel = (3*a*s*s + 2*b*s + c)/T;
      return Vel;
    }
    else
    {
      double Acc = (6*a*s+2*b)/(T*T);
      return Acc;
    }
  }
}

dlib::matrix<double> Dynamics_RHS_Matrix_fn(dlib::matrix<double> &Jac_Full, dlib::matrix<double> &B_q)
{
  dlib::matrix<double> Jac_Full_Trans, Dynamics_RHS_Matrix;

  Jac_Full_Trans = dlib::trans(Jac_Full);

  const int Dynamics_RHS_Matrix_Row = Jac_Full_Trans.nr();
  const int Dynamics_RHS_Matrix_Col = Jac_Full_Trans.nc() + B_q.nc();
  Dynamics_RHS_Matrix = dlib::zeros_matrix<double>(Dynamics_RHS_Matrix_Row, Dynamics_RHS_Matrix_Col);
  for (int i = 0; i < Dynamics_RHS_Matrix_Col; i++)
  {
    if (i<Jac_Full_Trans.nc())
    {
      dlib::set_colm(Dynamics_RHS_Matrix, i) = dlib::colm(Jac_Full_Trans,i);
    }
    else
    {
      dlib::set_colm(Dynamics_RHS_Matrix, i) = dlib::colm(B_q,i-Jac_Full_Trans.nc());
    }
  }
  return Dynamics_RHS_Matrix;
}

int Sigma_Change(std::vector<double> &sigma_i, std::vector<double> &sigma_i_child)
{
  // This function is only used when a contact change is made for sure
  int Sigma_Index = 0;
  double Sigma_Result = 0;
  for (int i = 0; i < 4; i++)
  {
    Sigma_Result = sigma_i[i] - sigma_i_child[i];
    if(abs(Sigma_Result)==1)
    {
      Sigma_Index = i;
      break;
    }
  }
  return Sigma_Index;
}
