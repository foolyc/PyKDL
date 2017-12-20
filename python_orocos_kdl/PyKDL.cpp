#include <string>
#include <boost/python.hpp>
#include <boost/python/wrapper.hpp>
# include <boost/python/detail/wrapper_base.hpp>


#include <kdl/frames.hpp>
#include <kdl/frames_io.hpp>
#include <kdl/joint.hpp>
#include <kdl/kinfam_io.hpp>
#include <kdl/segment.hpp>
#include <kdl/kinfam_io.hpp>
#include <kdl/chain.hpp>
#include <kdl/jntarray.hpp>
#include <kdl/jacobian.hpp>

#include <kdl/chainfksolver.hpp>
#include <kdl/chainiksolver.hpp>
#include <kdl/chainfksolverpos_recursive.hpp>
#include <kdl/solveri.hpp>
#include <kdl/chainiksolvervel_pinv.hpp>
#include <kdl/chainiksolverpos_nr.hpp>
#include <kdl/chainiksolverpos_nr_jl.hpp>

#include <kdl/chainjnttojacsolver.hpp>


#include <sstream>

using namespace boost::python;
using namespace KDL;



int    (ChainFkSolverPos_recursive::*JntToCart_1)(const JntArray &, Frame &, int)              = &ChainFkSolverPos_recursive::JntToCart;

void JntToCart(const Chain& chain, const JntArray& q_in, Frame& p_out)
{
	ChainFkSolverPos_recursive fk_solver = ChainFkSolverPos_recursive(chain);
	fk_solver.JntToCart(q_in, p_out);
}


int CartToJnt(const Chain& chain, const JntArray& q_init, Frame& F_dest, JntArray& q)
{
	ChainFkSolverPos_recursive fksolver(chain);//Forward position solver
	ChainIkSolverVel_pinv iksolverv(chain);//Inverse velocity solver
	ChainIkSolverPos_NR iksolverpos(chain,fksolver,iksolverv,100,1e-6);//Maximum 100 iterations, stop at accuracy 1e-6
  
	int ret = iksolverpos.CartToJnt(q_init,F_dest,q);
	return ret;

}



// expose index operator

void Vector_setitem(Vector& v, int index, double value)
{
    if (index >= 0 && index < 3)
    {
        v[index] = value;
    }
    else {
        PyErr_SetString(PyExc_IndexError, "index out of range");
        throw_error_already_set();
    }
}

double Vector_getitem(Vector& v, int index)
{
    if (index >= 0 && index < 3)
    {
        return v[index];
    }
    else {
        PyErr_SetString(PyExc_IndexError, "index out of range");
        throw_error_already_set();
    }
}


void JntArray_setitem(JntArray& q, int index, double value)
{
    if (index >= 0 && index < (q.rows())*(q.columns()))
    {
        q(index) = value;
    }
    else {
        PyErr_SetString(PyExc_IndexError, "index out of range");
        throw_error_already_set();
    }
}

double JntArray_getitem(JntArray& q, int index)
{
    if (index >= 0 && index < (q.rows())*(q.columns()))
    {
        return q(index);
    }
    else {
        PyErr_SetString(PyExc_IndexError, "index out of range");
        throw_error_already_set();
    }
}


void Rotation_setitem(Rotation& r, int index, double value)
{

    if (index >= 0 && index < 9)
    {
    	int i = int(index/3);
		int j= index - 3*i;
        r(i, j) = value;
    }
    else {
        PyErr_SetString(PyExc_IndexError, "index out of range");
        throw_error_already_set();
    }
}

double Rotation_getitem(Rotation& r, int index)
{
   if (index >= 0 && index < 9)
    {
    	int i = int(index/3);
		int j= index - 3*i;
        return r(i, j);
    }
    else {
        PyErr_SetString(PyExc_IndexError, "index out of range");
        throw_error_already_set();
    }
}

void Jacobian_setitem(Jacobian& jac, int index, double value)
{

    if (index >= 0 && index < jac.rows()*jac.columns())
    {
    	int i = int(index/jac.columns());
		int j= index - jac.columns()*i;
        jac(i, j) = value;
    }
    else {
        PyErr_SetString(PyExc_IndexError, "index out of range");
        throw_error_already_set();
    }
}

double Jacobian_getitem(Jacobian& jac, int index)
{
   if (index >= 0 && index < jac.rows()*jac.columns())
    {
    	int i = int(index/jac.columns());
		int j= index - jac.columns()*i;
        return jac(i, j);
    }
    else {
        PyErr_SetString(PyExc_IndexError, "index out of range");
        throw_error_already_set();
    }
}



//Wrapper
BOOST_PYTHON_MODULE(PyKDL)
{
	Py_Initialize();

	class_<Vector>("Vector")
		.def(init<double, double, double>())
		.def("__getitem__", &Vector_getitem)
		.def("__setitem__", &Vector_setitem)

		.def(self + self)
		.def(self - self)
		.def(self * double())
		.def(self / double())

	;

	class_<Rotation>("Rotation")
		.def(init<double, double, double, double, double, double, double, double, double>())
		.def(init<Vector, Vector, Vector>())
		.def("RPY", &Rotation::RPY)
		.def("Quaternion", &Rotation::Quaternion)
		.def("EulerZYZ", &Rotation::EulerZYZ)
		.def("Rot", &Rotation::Rot)
		.def("Rot2", &Rotation::Rot2)
		.def("__getitem__", &Rotation_getitem)
		.def("__setitem__", &Rotation_setitem)
		.def("GetRPY", &Rotation::GetRPY)
		.def("GetEulerZYZ", &Rotation::GetEulerZYZ)
		.def("GetQuaternion", &Rotation::GetQuaternion)
		.def("GetRot", &Rotation::GetRot)
		.def("GetRotAngle", &Rotation::GetRotAngle)
		.def(self * self)

	;

	class_<Jacobian>("Jacobian")
		.def("rows", &Jacobian::rows)
		.def("columns", &Jacobian::columns)
		.def("resize", &Jacobian::resize)
		.def("__getitem__", &Jacobian_getitem)
		.def("__setitem__", &Jacobian_setitem)
	;


	class_<Frame>("Frame")
		.def(init<Rotation, Vector>())
		.def(init<Rotation>())
		.def(init<Vector>())
		.def_readonly("p", &Frame::p)
		.def_readonly("M", &Frame::M)
	;



	class_<Joint>("Joint")
		.def(init<Joint::JointType>())
		.def_readonly("RotAxis", Joint::RotAxis)
		.def_readonly("RotX", Joint::RotX)
		.def_readonly("RotY", Joint::RotY)
		.def_readonly("RotZ", Joint::RotZ)
		.def_readonly("TransAxis", Joint::TransAxis)
		.def_readonly("TransX", Joint::TransX)
		.def_readonly("TransY", Joint::TransY)
		.def_readonly("TransZ", Joint::TransZ)
		.def_readonly("NoJoint", Joint::None)
	;


	enum_<Joint::JointType>("JointType")
    	.value("RotAxis", Joint::RotAxis)
    	.value("RotX", Joint::RotX)
    	.value("RotY", Joint::RotY)
    	.value("RotZ", Joint::RotZ)
    	.value("TransAxis", Joint::TransAxis)
    	.value("TransX", Joint::TransX)
    	.value("TransY", Joint::TransY)
    	.value("TransZ", Joint::TransZ)
    	.value("NoJoint", Joint::None)
	;


	class_<Segment>("Segment")
		.def(init<Joint, Frame>())

	;

	class_<Chain>("Chain")
		.def("addSegment", &Chain::addSegment)
		.def("addChain", &Chain::addChain)
		.def("getNrOfJoints", &Chain::getNrOfJoints)
		.def("getNrOfSegments", &Chain::getNrOfSegments)
	;

	class_<JntArray>("JntArray")
		.def(init<int>())
		.def("rows", &JntArray::rows)
		.def("columns", &JntArray::columns)
		.def("__getitem__", &JntArray_getitem)
		.def("__setitem__", &JntArray_setitem)
	;

	class_<SolverI, boost::noncopyable>("SolverI", no_init);
	
	class_<ChainFkSolverPos, bases<SolverI>, boost::noncopyable>("ChainFkSolverPos", no_init);

	class_<ChainFkSolverPos_recursive, bases<ChainFkSolverPos>, boost::noncopyable >("ChainFkSolverPos_recursive", init<Chain>())
		.def("JntToCart", JntToCart_1)
	;

	class_<ChainIkSolverPos, bases<SolverI>, boost::noncopyable>("ChainIkSolverPos", no_init);

	class_<ChainIkSolverVel, bases<SolverI>, boost::noncopyable>("ChainIkSolverVel", no_init);

	class_<ChainIkSolverVel_pinv, bases<ChainIkSolverVel>, boost::noncopyable >("ChainIkSolverVel_pinv", init<Chain>());

	class_<ChainIkSolverPos_NR_JL, bases<ChainIkSolverPos>, boost::noncopyable >("ChainIkSolverPos_NR_JL", init<const Chain&, const JntArray, const JntArray, ChainFkSolverPos &, ChainIkSolverVel &, unsigned int, double>())
		// .def(init<const Chain&, ChainFkSolverPos &, ChainIkSolverVel &, unsigned int, double>())
		.def("CartToJnt", &ChainIkSolverPos_NR_JL::CartToJnt)
	;

	class_<ChainJntToJacSolver, bases<SolverI>, boost::noncopyable>("ChainJntToJacSolver", init<Chain>())
		.def("JntToJac", &ChainJntToJacSolver::JntToJac)
	;



	def("JntToCart",JntToCart);

	def("CartToJnt",CartToJnt);
}


int main()
{
	PyImport_AppendInittab("PyKDL", &PyInit_PyKDL); // Add example to built-in.
	Py_Initialize(); // Start interpreter.
	try
	{
		object main = import("__main__");
    	object main_namespace = main.attr("__dict__");
    	scope scope(main); // Force main scope
    	main_namespace["PyKDL"] = import("PyKDL");
	}
	catch (const error_already_set&)
  	{
		PyErr_Print();
  	}

}