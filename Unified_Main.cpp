#include <stdio.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <numeric>

#include "Unified_Header.h"
#include "snfilewrapper.hh"
#include "snopt.hh"
#include "snoptProblem.hh"

using namespace std;
Unified_Structure_P Structure_P;
extern int Variable_Num;

int main( int argc, char **argv)
{
	// This function is used to generate the contact transition tree
	std::vector<double> sigma_init;

	// sigma_init is used initialize the robot contact status at the initial time
	ifstream Sigma_Init_File;              // This is to read the initial angle and angular velocities
	Sigma_Init_File.open("sigma_init.txt");
	if(Sigma_Init_File.is_open())
	{
		double data_each_line = 0.0;

		while(Sigma_Init_File>>data_each_line)
		{
			sigma_init.push_back(data_each_line);
		}
		Sigma_Init_File.close();
	}
	else
	{
		printf("Unable to open sigma_init.txt file!\n");
	}
	std::vector<double> Robot_Init_vec = Default_Init(sigma_init);
	Robot_StateNDot StateNDot_Init_Opt(Robot_Init_vec);
	std::string input_name = "Initial State Plot at Beginning";
	Robot_Plot_fn(StateNDot_Init_Opt,input_name);
	double Init_KE = Kinetic_Energy_fn(StateNDot_Init_Opt);
	cout<<"-----------------------------------Initial Kinetic Energy: "<<Init_KE<<"-----------------------------------"<<endl;
	// After the robot state initialization, the next job is to conduct the multi-contact staiblization strategy: the root node initialization
	Tree_Node Root_Node;
	Node_UpdateNCon(Root_Node, StateNDot_Init_Opt, sigma_init);			// Node_UpdateNCon can only be used when the state and sigma at a given node is known
	Tree_Node Node_i;
	int Opt_Flag = 0;
	int Iteration = 0;
	while(Frontier_Nodes.size()>0)
	{
		/**
		* For the current node, first is the Node_Self_Opt to optimize a motion while maintain the current mode
		* if this does not work, then expand the current node into the adjacent nodes then do the Nodes_Connectivity_Opt
		*/
		std::vector<double> Opt_Soln(Variable_Num);
		Node_i = Pop_Node();
		// if(Iteration>0)
		// {
			Opt_Flag = Nodes_Optimization_fn(Node_i, Node_i, Opt_Soln);
		// }
		if(Opt_Flag==1)
		{
			// Optimal solution has been found
			Opt_Soln_Write2Txt(Node_i, Node_i, Opt_Soln);
			break;
		}
		else
		{
			// Here it is the node expansion
			int Adjacent_Number;	int Nodes_Opt_Flag;
			dlib::matrix<double> Nodes_Sigma_Matrix = Node_Expansion_fn(Node_i, Adjacent_Number);
			cout<<Nodes_Sigma_Matrix<<endl;
			for (int i = 0; i < Adjacent_Number; i++)
			{
				Tree_Node Node_i_child;
				Node_i_child.sigma.push_back(Nodes_Sigma_Matrix(i,0));
				Node_i_child.sigma.push_back(Nodes_Sigma_Matrix(i,1));
				Node_i_child.sigma.push_back(Nodes_Sigma_Matrix(i,2));
				Node_i_child.sigma.push_back(Nodes_Sigma_Matrix(i,3));
				Nodes_Opt_Flag = Nodes_Optimization_fn(Node_i, Node_i_child, Opt_Soln);
				if(Nodes_Opt_Flag==1)
				{
					std::vector<double> Robot_StateNDot_vec_i;
					Robot_StateNDot_vec_i = End_RobotNDot_Extract(Opt_Soln, Node_i.sigma, Node_i_child.sigma);
					Robot_StateNDot Robot_StateNDot_Child_i(Robot_StateNDot_vec_i);
					Node_UpdateNCon(Node_i_child, Robot_StateNDot_Child_i, Node_i_child.sigma);
					Node_i_child.Parent_Node = &Node_i;
					Node_i.Children_Nodes.push_back(&Node_i_child);
					Opt_Soln_Write2Txt(Node_i, Node_i_child, Opt_Soln);
				}
			}
		}
		Iteration = Iteration + 1;
	}
	return 0;
}
