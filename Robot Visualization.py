#! /usr/bin/env python

import sys, os, time
import numpy as np
import matplotlib.pyplot as plt
import logging
import pickle as pkl
from klampt import WorldModel
from klampt import vis
from klampt.math import vectorops
from klampt.model.trajectory import Trajectory
import ipdb

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

Grids = 0
pi = 3.1415926535897932384626
Aux_Link_Ind = [1, 3, 5, 6, 7, 11, 12, 13, 17, 18, 19, 20, 21, 23, 24, 26, 27, 28, 30, 31, 33, 34, 35]  # This is the auxilliary joints causing lateral motion
Act_Link_Ind = [0, 2, 4, 8, 9, 10, 14, 15, 16, 22, 25, 29, 32]                                          # This is the active joints to be considered
Ctrl_Link_Ind = Act_Link_Ind[3:]
Tot_Link_No = len(Aux_Link_Ind) + len(Act_Link_Ind)
Local_Extremeties = [0.1, 0, -0.1, -0.15 , 0, -0.1, 0.1, 0, -0.1, -0.15 , 0, -0.1, 0, 0, -0.22, 0, 0, -0.205]         # This 6 * 3 vector describes the local coordinate of the contact extremeties in their local coordinate
End_Effector_Ind = [11, 11, 17, 17, 27, 34]                       # The link index of the end effectors

class MyGLPlugin(vis.GLPluginInterface):
    def __init__(self, world):
        vis.GLPluginInterface.__init__(self)
        self.world = world
        self.quit = False
        self.starp = False

    def mousefunc(self, button, state, x, y):
        print("mouse",button,state,x,y)
        if button==2:
            if state==0:
                print("Click list...",[o.getName() for o in self.click_world(x,y)])
            return True
        return False

    def motionfunc(self, x, y, dx, dy):
        return False

    def keyboardfunc(self, c, x, y):
        print("Pressed",c)
        if c == 'q':
            self.quit = True
            return True
        if c == 's':
            self.starp = not self.starp
            return True
        return False

    def click_world(self, x, y):
        """Helper: returns a list of world objects sorted in order of
        increasing distance."""
        #get the viewport ray
        (s, d) = self.click_ray(x, y)

        #run the collision tests
        collided = []
        for g in self.collider.geomList:
            (hit, pt) = g[1].rayCast(s, d)
            if hit:
                dist = vectorops.dot(vectorops.sub(pt, s), d)
                collided.append((dist,g[0]))
        return [g[1] for g in sorted(collided)]


def getPath(data):
    """Load configurations and return a traj, i.e. list of list

    :param data: a dictionary structure from parsing the solution. It has 't', 'x', 'u', etc
    :return t, traj, force: time stamp, configuration, and support force

    """
    t = data['t']
    x = data['x'][:, :7]  # only use q
    force = data['p']
    x_offset = 0
    # parse into traj
    traj = []
    for i in range(len(t)):
        qT, q1, q2, p1, p2, qx, qy = x[i]
        q = np.zeros(10)
        q[0] = -qx  # x-loc
        q[2] = qy  # z-loc
        q[4] = qT
        q[6] = q1
        q[7] = q2
        q[8] = p1
        q[9] = p2
        traj.append(q.tolist())
    return t, traj, force


def copy_copy(q, switched):
    qq = copy.copy(q)
    if switched:
        qq[6], qq[8] = qq[8], qq[6]
        qq[7], qq[9] = qq[9], qq[7]
    return qq

def Dimension_Recovery(low_dim_obj):
    high_dim_obj = np.zeros(Tot_Link_No)
    for i in range(0,len(Act_Link_Ind)):
        high_dim_obj[Act_Link_Ind[i]] = low_dim_obj[i]
    return high_dim_obj

def main():
    #creates a world and loads all the items on the command line
    world = WorldModel()
    res = world.readFile("./HRP2_Robot.xml")
    if not res:
        raise RuntimeError("Unable to load model")
    robot = world.robot(0)
    # coordinates.setWorldModel(world)
    plugin = MyGLPlugin(world)
    vis.pushPlugin(plugin)
    #add the world to the visualizer
    vis.add("world", world)
    vis.add("robot",world.robot(0))
    vis.show()
    # ipdb.set_trace()

    Robotstate_Traj, Contact_Force_Traj = Traj_Loader()
    h = 0.0068               # 0.0043: vert
                            # 0.0033: flat
    playspeed = 2.5
    norm = 500
    while vis.shown():
        # This is the main plot program
        for i in range(0, Robotstate_Traj.shape[1]):
            vis.lock()
            Robotstate_Traj_i = Robotstate_Traj[:,i]
            Robotstate_Traj_Full_i = Dimension_Recovery(Robotstate_Traj_i)
            # Now it is the plot of the contact force at the contact extremities
            robot.setConfig(Robotstate_Traj_Full_i)
            Contact_Force_Traj_i = Contact_Force_Traj[:,i]
            left_ft, rght_ft, left_hd, rght_hd, left_ft_end, rght_ft_end, left_hd_end, rght_hd_end = Contact_Force_vec(robot, Contact_Force_Traj_i, norm)
            vis.add("left foot force", Trajectory([0, 1], [left_ft, left_ft_end]))
            vis.add("right foot force", Trajectory([0, 1], [rght_ft, rght_ft_end]))
            vis.add("left hand force", Trajectory([0, 1], [left_hd, left_hd_end]))
            vis.add("right hand force", Trajectory([0, 1], [rght_hd, rght_hd_end]))
            COMPos_start = robot.getCom()
            COMPos_end = COMPos_start
            COMPos_end[2] = COMPos_end[2] + 100
            vis.add("Center of Mass",  Trajectory([0, 1], [COMPos_start, COMPos_end]))

            # ipdb.set_trace()

            vis.unlock()
            time.sleep(h * playspeed)
            # vis.hide('left_ft force')
            # vis.hide('rght_ft force')
            # vis.hide('left_hd force')
            # vis.hide('rght_hd force')
def Contact_Force_vec(robot, Contact_Force_Traj_i, norm):
    length = 0.5
    End_Effector_Pos = get_End_Effector_Pos(robot)          # 18 by 1
    left_ft_x = 0.5 * End_Effector_Pos[0] + 0.5 * End_Effector_Pos[3]
    left_ft_y = 0.5 * End_Effector_Pos[1] + 0.5 * End_Effector_Pos[4]
    left_ft_z = 0.5 * End_Effector_Pos[2] + 0.5 * End_Effector_Pos[5]
    rght_ft_x = 0.5 * End_Effector_Pos[6] + 0.5 * End_Effector_Pos[9]
    rght_ft_y = 0.5 * End_Effector_Pos[7] + 0.5 * End_Effector_Pos[10]
    rght_ft_z = 0.5 * End_Effector_Pos[8] + 0.5 * End_Effector_Pos[11]
    left_hd_x = End_Effector_Pos[12]
    left_hd_y = End_Effector_Pos[13]
    left_hd_z = End_Effector_Pos[14]
    rght_hd_x = End_Effector_Pos[15]
    rght_hd_y = End_Effector_Pos[16]
    rght_hd_z = End_Effector_Pos[17]

    left_ft = [left_ft_x, left_ft_y, left_ft_z]
    rght_ft = [rght_ft_x, rght_ft_y, rght_ft_z]
    left_hd = [left_hd_x, left_hd_y, left_hd_z]
    rght_hd = [rght_hd_x, rght_hd_y, rght_hd_z]

    left_ft_x_off = Contact_Force_Traj_i[0]/norm * length
    left_ft_z_off = Contact_Force_Traj_i[1]/norm * length
    rght_ft_x_off = Contact_Force_Traj_i[2]/norm * length
    rght_ft_z_off = Contact_Force_Traj_i[3]/norm * length
    left_hd_x_off = Contact_Force_Traj_i[4]/norm * length
    left_hd_z_off = Contact_Force_Traj_i[5]/norm * length
    rght_hd_x_off = Contact_Force_Traj_i[6]/norm * length
    rght_hd_z_off = Contact_Force_Traj_i[7]/norm * length

    left_ft_end = [left_ft_x + left_ft_x_off, left_ft_y, left_ft_z + left_ft_z_off]
    rght_ft_end = [rght_ft_x + rght_ft_x_off, rght_ft_y, rght_ft_z + rght_ft_z_off]
    left_hd_end = [left_hd_x + left_hd_x_off, left_hd_y, left_hd_z + left_hd_z_off]
    rght_hd_end = [rght_hd_x + rght_hd_x_off, rght_hd_y, rght_hd_z + rght_hd_z_off]

    return left_ft, rght_ft, left_hd, rght_hd, left_ft_end, rght_ft_end, left_hd_end, rght_hd_end


def get_End_Effector_Pos(hrp2_robot):
    End_Effector_Pos_Array = np.array([0, 0, 0])
    End_Link_No_Index = -1
    for End_Effector_Link_Index in End_Effector_Ind:
        End_Link_No_Index = End_Link_No_Index + 1
        End_Link_i = hrp2_robot.link(End_Effector_Link_Index)
        End_Link_i_Extre_Loc = Local_Extremeties[End_Link_No_Index*3:End_Link_No_Index*3+3]
        End_Link_i_Extre_Pos = End_Link_i.getWorldPosition(End_Link_i_Extre_Loc)
        End_Effector_Pos_Array = np.append(End_Effector_Pos_Array, End_Link_i_Extre_Pos)
    return End_Effector_Pos_Array[3:]

def Robot_ConfigNVel_Update(robot, x):
    OptConfig_Low = x[0:len(x)/2]
    OptVelocity_Low = x[len(x)/2:]
    OptConfig_High = Dimension_Recovery(OptConfig_Low)
    OptVelocity_High = Dimension_Recovery(OptVelocity_Low)
    robot.setConfig(OptConfig_High)
    robot.setVelocity(OptVelocity_High)

def KE_fn(robot, dataArray):
    #First we have to set the robot to be the corresponding configuration and angular velocities
    Robot_ConfigNVel_Update(robot, dataArray)
    D_q = np.asarray(robot.getMassMatrix())    # Here D_q is the effective inertia matrix
    qdot_i = dataArray[13:None]
    qdot_i = np.reshape(qdot_i,[13,1])
    qdot_i = Dimension_Recovery(qdot_i)
    qdot_i_trans = np.transpose(qdot_i)
    T = 0.5 * qdot_i_trans.dot(D_q.dot(qdot_i))
    return T

def Traj_Loader():
    # This function will load the robotstate and contact force trajectories
    Robotstate_Traj = np.array([])
    # with open("./Exp Data/Exp 4/State_Exp3_.txt",'r') as robot_soln_file:
    # with open("./Exp Data/Highlights/Vert_Wall/State_Vert_.txt",'r') as robot_soln_file:
    # with open("./Exp Data/Highlights/Flat_Gnd/State_flat_.txt",'r') as robot_soln_file:
    with open("./Exp Data/Evenly/State_Flat_150_.txt",'r') as robot_soln_file:
        for line in robot_soln_file:
            currentline = line.split(",")
            currentline = [x.replace("\r\n","") for x in currentline]
            currentline = map(float, currentline)
            currentline = np.array([currentline])
            Robotstate_Traj = np.append(Robotstate_Traj, currentline)
    Robotstate_Traj = np.reshape(Robotstate_Traj, (13, Robotstate_Traj.shape[0]/13))

    Contact_Force_Traj = np.array([])
    # with open("./Exp Data/Exp 4/Contact_Force_Exp3_.txt",'r') as robot_soln_file:
    # with open("./Exp Data/Highlights/Vert_Wall/Contact_Force_Vert_.txt",'r') as robot_soln_file:
    # with open("./Exp Data/Highlights/Flat_Gnd/Contact_Force_flat_.txt",'r') as robot_soln_file:
    with open("./Exp Data/Evenly/Contact_Force_Flat_150_.txt",'r') as robot_soln_file:
        for line in robot_soln_file:
            currentline = line.split(",")
            currentline = [x.replace("\r\n","") for x in currentline]
            currentline = map(float, currentline)
            currentline = np.array([currentline])
            Contact_Force_Traj = np.append(Contact_Force_Traj, currentline)
    Contact_Force_Traj = np.reshape(Contact_Force_Traj, (8, Contact_Force_Traj.shape[0]/8))
    return Robotstate_Traj, Contact_Force_Traj
def Path_Loader():
    # This function is used to read the Opt_Soln txt file
    global Grids
    with open("Opt_Soln3.txt",'r') as robot_soln_file:
        robot_soln_i = robot_soln_file.readlines()
        robot_soln = [x.strip() for x in robot_soln_i]
        robot_soln = np.array(robot_soln, dtype = float)
    Grids = (robot_soln.shape[0]-1)/48
    T_tot, StateNDot_Traj, Ctrl_Traj, Contact_Force_Traj = Seed_Guess_Unzip(robot_soln)
    return T_tot, StateNDot_Traj, Ctrl_Traj, Contact_Force_Traj

def Seed_Guess_Unzip(Opt_Seed):
    # This function is used to unzip the Opt_Seed back into the coefficients
    T = Opt_Seed[0]
    # First is to retrieve the stateNdot coefficients from the Opt_Seed: 1 to 1 + (Grids - 1) * 26 * 4
    Var_Init = 1
    Var_Length = (Grids) * 26
    Var_End = Var_Init + Var_Length
    StateNDot_Traj = Opt_Seed[Var_Init:Var_End]
    StateNDot_Traj = np.reshape(StateNDot_Traj, (Grids, 26)).transpose()

    Var_Init = Var_End
    Var_Length = (Grids) * 10
    Var_End = Var_Init + Var_Length
    Ctrl_Traj = Opt_Seed[Var_Init:Var_End]
    Ctrl_Traj = np.reshape(Ctrl_Traj, (Grids, 10)).transpose()

    Var_Init = Var_End
    Var_Length = (Grids) * 12
    Var_End = Var_Init + Var_Length
    Contact_Force_Traj = Opt_Seed[Var_Init:Var_End]
    Contact_Force_Traj = np.reshape(Contact_Force_Traj, (Grids, 12)).transpose()
    # ipdb.set_trace()
    return T, StateNDot_Traj, Ctrl_Traj, Contact_Force_Traj

def showTwoStep(plugin, world, robot, ground1, ground2, ground3, ground4, link_foot, link_foot_other):
    """Read the data file and select trajectories to show, they are composed of two steps"""
    data = np.load('data/twoStepBunch.npz')['output']
    ndata = len(data)
    N = 20
    force_len = 0.5
    vis.show()
    for j in range(ndata):
        i = j
        print('Entering traj %d' % i)
        l0, l1, h0, h1, vel, phase0, phase1 = data[i]
        # set those grounds
        ground1.setConfig([0, 0, 0, 0, 0, 0])
        ground2.setConfig([l0, 0, h0, 0, 0, 0])
        ground3.setConfig([l0 + l1, 0, h0 + h1, 0, 0, 0])
        ground4.setConfig([2 * l0 + l1, 0, 2 * h0 + h1, 0, 0, 0])
        while vis.shown() and not plugin.quit:
            if plugin.nextone:  # check if we want next one
                plugin.nextone = False
                break
            # show phase0
            t, q, force = getTraj(phase0)
            h_ = t[1] - t[0]
            if h_ < 0:
                break
            nPoint = len(t)
            for k in range(nPoint):
                vis.lock()
                useq = copy_copy(q[k], False)
                robot.setConfig(useq)
                footpos = link_foot.getWorldPosition([0, 0, 0.5])
                support = np.array([force[k, 0], 0, force[k, 1]])
                use_support = support / np.linalg.norm(support) * force_len
                force_end = vectorops.add(footpos, use_support.tolist())
                vis.add("force", Trajectory([0, 1], [footpos, force_end]))
                vis.unlock()
                time.sleep(h_ * 5.0)
                vis.remove('force')
            # phase 1
            t, q, force = getTraj(phase1)
            h_ = t[1] - t[0]
            if h_ < 0:
                break
            print('h_ = %f' % h_)
            for k in range(nPoint):
                vis.lock()
                useq = copy_copy(q[k], True)
                robot.setConfig(useq)
                footpos = link_foot_other.getWorldPosition([0, 0, 0.5])
                support = np.array([force[k, 0], 0, force[k, 1]])
                use_support = support / np.linalg.norm(support) * force_len
                force_end = vectorops.add(footpos, use_support.tolist())
                vis.add("force", Trajectory([0, 1], [footpos, force_end]))
                vis.unlock()
                time.sleep(h_ * 5.0)
                vis.remove('force')
    while vis.shown() and not plugin.quit:
        vis.lock()
        vis.unlock()
        time.sleep(0.05)
    print("Ending visualization.")
    vis.kill()

if __name__ == '__main__':
    main()
