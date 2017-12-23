import numpy as np
from numpy import pi as PI
import sys

import PyKDL
from PyKDL import *



fetch_chain = Chain()
fetch_chain.addSegment(Segment(Joint(Joint.NoJoint), Frame(Rotation.RPY(0, 0, 0), Vector(0.0, 0.0, 0.0))))
fetch_chain.addSegment(Segment(Joint(Joint.NoJoint), Frame(Rotation.RPY(0, 0, 0), Vector(0, -0.19, 0.3))))
fetch_chain.addSegment(Segment(Joint(Joint.NoJoint), Frame(Rotation.RPY(0, 0, 0), Vector(0, 0, 0))))
fetch_chain.addSegment(Segment(Joint(Joint.RotY), Frame(Rotation.RPY(0, 0, 0), Vector(0.0, 0.0, 0.0))))
fetch_chain.addSegment(Segment(Joint(Joint.RotZ), Frame(Rotation.RPY(0, 0, 0), Vector(0.0, 0.0, 0.0))))
fetch_chain.addSegment(Segment(Joint(Joint.RotX), Frame(Rotation.RPY(0, 0, 0), Vector(0.19, 0, 0))))
fetch_chain.addSegment(Segment(Joint(Joint.RotY), Frame(Rotation.RPY(0, 0, 0), Vector(0, 0, 0))))
fetch_chain.addSegment(Segment(Joint(Joint.RotX), Frame(Rotation.RPY(0, 0, 0), Vector(0.28, 0, 0))))
fetch_chain.addSegment(Segment(Joint(Joint.RotY), Frame(Rotation.RPY(0, 0, 0), Vector(0.06, 0.03, 0))))
fetch_chain.addSegment(Segment(Joint(Joint.NoJoint), Frame(Rotation.RPY(0, 0, 0), Vector(0.05, 0.0, 0.0))))

print("use fk and ik in module derectly ! ")

q0 = JntArray(6)
q0_list = [0.7, -0.4, -1.5, -1.4, 1.0, -0.2]
for i in range(6):
	q0[i] = q0_list[i]
p0 = Frame()

JntToCart(fetch_chain, q0, p0)
print(p0.p[0], p0.p[1], p0.p[2])

q1 = JntArray(6)
CartToJnt(fetch_chain, q0, p0, q1)
for i in range(6):
	print(q1[i])




print("use fk and ik in solver ! ")

jac_solver = ChainJntToJacSolver(fetch_chain)

ik_v_kdl = ChainIkSolverVel_pinv(fetch_chain)
fk_kdl = ChainFkSolverPos_recursive(fetch_chain)

q_min = JntArray(6)
q_max = JntArray(6)
for i in range(6):
    q_min[i] = -3
    q_max[i] = 3


ik_p_kdl = ChainIkSolverPos_NR_JL(fetch_chain, q_min, q_max, fk_kdl, ik_v_kdl)


p0 = Frame()
q1 = JntArray(6)
fk_kdl.JntToCart(q0, p0, -1)


p1 = Frame(Rotation.RPY(-0.05, 0.1, 1.1), Vector(0.3, -0.05, 0.07))
ik_p_kdl.CartToJnt(q0, p0, q1)
for i in range(6):
	print(q1[i])
