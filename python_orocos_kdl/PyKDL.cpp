#include <string>
#include <boost/python.hpp>
#include <boost/python/wrapper.hpp>
#include <boost/python/detail/wrapper_base.hpp>
#include "boost/python/suite/indexing/vector_indexing_suite.hpp"


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
#include "kdl/chainfksolvervel_recursive.hpp"
#include <kdl/solveri.hpp>
#include <kdl/chainiksolvervel_pinv.hpp>
#include <kdl/chainiksolverpos_nr.hpp>
#include <kdl/chainiksolverpos_nr_jl.hpp>

#include <kdl/chainjnttojacsolver.hpp>

#include "kdl/chainiksolverpos_lma.hpp"
#include "kdl/chainiksolvervel_wdls.hpp"
#include "kdl/chainiksolvervel_pinv_nso.hpp"
#include "kdl/chainiksolvervel_pinv_givens.hpp"
#include "kdl/chainidsolver.hpp"
#include "kdl/chainidsolver_recursive_newton_euler.hpp"

#include <Eigen/Dense>
#include <sstream>

namespace bp = boost::python;
typedef double ScalarType;
typedef Eigen::Matrix<ScalarType,Eigen::Dynamic,Eigen::Dynamic> MatrixXq;
typedef Eigen::Matrix<ScalarType,Eigen::Dynamic,1> VectorXq;


double epsilon = 0.0000001;



struct ChainFkSolverPos_wrapper : KDL::ChainFkSolverPos, bp::wrapper< KDL::ChainFkSolverPos > {

    ChainFkSolverPos_wrapper()
            : KDL::ChainFkSolverPos()
            , bp::wrapper< KDL::ChainFkSolverPos >(){
        // null constructor

    }

    virtual int JntToCart( ::KDL::JntArray const & q_in, ::KDL::Frame & p_out, int segmentNr=-1 ){
        bp::override func_JntToCart = this->get_override( "JntToCart" );
        return func_JntToCart( boost::ref(q_in), boost::ref(p_out), segmentNr );
    }

    virtual int getError(  ) const  {
        if( bp::override func_getError = this->get_override( "getError" ) )
            return func_getError(  );
        else{
            return this->KDL::SolverI::getError(  );
        }
    }

    int default_getError(  ) const  {
        return KDL::SolverI::getError( );
    }

    virtual char const * strError( int const error ) const  {
        if( bp::override func_strError = this->get_override( "strError" ) )
            return func_strError( error );
        else{
            return this->KDL::SolverI::strError( error );
        }
    }

    char const * default_strError( int const error ) const  {
        return KDL::SolverI::strError( error );
    }

};

struct ChainFkSolverPos_recursive_wrapper : KDL::ChainFkSolverPos_recursive, bp::wrapper< KDL::ChainFkSolverPos_recursive > {

    ChainFkSolverPos_recursive_wrapper(::KDL::Chain const & chain )
            : KDL::ChainFkSolverPos_recursive( boost::ref(chain) )
            , bp::wrapper< KDL::ChainFkSolverPos_recursive >(){
        // constructor

    }

    virtual int JntToCart( ::KDL::JntArray const & q_in, ::KDL::Frame & p_out, int segmentNr=-1 ) {
        if( bp::override func_JntToCart = this->get_override( "JntToCart" ) )
            return func_JntToCart( boost::ref(q_in), boost::ref(p_out), segmentNr );
        else{
            return this->KDL::ChainFkSolverPos_recursive::JntToCart( boost::ref(q_in), boost::ref(p_out), segmentNr );
        }
    }

    int default_JntToCart( ::KDL::JntArray const & q_in, ::KDL::Frame & p_out, int segmentNr=-1 ) {
        return KDL::ChainFkSolverPos_recursive::JntToCart( boost::ref(q_in), boost::ref(p_out), segmentNr );
    }

    virtual int getError(  ) const  {
        if( bp::override func_getError = this->get_override( "getError" ) )
            return func_getError(  );
        else{
            return this->KDL::SolverI::getError(  );
        }
    }

    int default_getError(  ) const  {
        return KDL::SolverI::getError( );
    }

    virtual char const * strError( int const error ) const  {
        if( bp::override func_strError = this->get_override( "strError" ) )
            return func_strError( error );
        else{
            return this->KDL::SolverI::strError( error );
        }
    }

    char const * default_strError( int const error ) const  {
        return KDL::SolverI::strError( error );
    }

};

struct ChainFkSolverVel_wrapper : KDL::ChainFkSolverVel, bp::wrapper< KDL::ChainFkSolverVel > {

    ChainFkSolverVel_wrapper()
            : KDL::ChainFkSolverVel()
            , bp::wrapper< KDL::ChainFkSolverVel >(){
        // null constructor

    }

    virtual int JntToCart( ::KDL::JntArrayVel const & q_in, ::KDL::FrameVel & out, int segmentNr=-1 ){
        bp::override func_JntToCart = this->get_override( "JntToCart" );
        return func_JntToCart( boost::ref(q_in), boost::ref(out), segmentNr );
    }

    virtual int getError(  ) const  {
        if( bp::override func_getError = this->get_override( "getError" ) )
            return func_getError(  );
        else{
            return this->KDL::SolverI::getError(  );
        }
    }

    int default_getError(  ) const  {
        return KDL::SolverI::getError( );
    }

    virtual char const * strError( int const error ) const  {
        if( bp::override func_strError = this->get_override( "strError" ) )
            return func_strError( error );
        else{
            return this->KDL::SolverI::strError( error );
        }
    }

    char const * default_strError( int const error ) const  {
        return KDL::SolverI::strError( error );
    }

};

struct ChainFkSolverVel_recursive_wrapper : KDL::ChainFkSolverVel_recursive, bp::wrapper< KDL::ChainFkSolverVel_recursive > {

    ChainFkSolverVel_recursive_wrapper(::KDL::Chain const & chain )
            : KDL::ChainFkSolverVel_recursive( boost::ref(chain) )
            , bp::wrapper< KDL::ChainFkSolverVel_recursive >(){
        // constructor

    }

    virtual int JntToCart( ::KDL::JntArrayVel const & q_in, ::KDL::FrameVel & out, int segmentNr=-1 ) {
        if( bp::override func_JntToCart = this->get_override( "JntToCart" ) )
            return func_JntToCart( boost::ref(q_in), boost::ref(out), segmentNr );
        else{
            return this->KDL::ChainFkSolverVel_recursive::JntToCart( boost::ref(q_in), boost::ref(out), segmentNr );
        }
    }

    int default_JntToCart( ::KDL::JntArrayVel const & q_in, ::KDL::FrameVel & out, int segmentNr=-1 ) {
        return KDL::ChainFkSolverVel_recursive::JntToCart( boost::ref(q_in), boost::ref(out), segmentNr );
    }

    virtual int getError(  ) const  {
        if( bp::override func_getError = this->get_override( "getError" ) )
            return func_getError(  );
        else{
            return this->KDL::SolverI::getError(  );
        }
    }

    int default_getError(  ) const  {
        return KDL::SolverI::getError( );
    }

    virtual char const * strError( int const error ) const  {
        if( bp::override func_strError = this->get_override( "strError" ) )
            return func_strError( error );
        else{
            return this->KDL::SolverI::strError( error );
        }
    }

    char const * default_strError( int const error ) const  {
        return KDL::SolverI::strError( error );
    }

};

struct ChainIdSolver_wrapper : KDL::ChainIdSolver, bp::wrapper< KDL::ChainIdSolver > {

    ChainIdSolver_wrapper()
            : KDL::ChainIdSolver()
            , bp::wrapper< KDL::ChainIdSolver >(){
        // null constructor

    }

    virtual int CartToJnt( ::KDL::JntArray const & q, ::KDL::JntArray const & q_dot, ::KDL::JntArray const & q_dotdot, ::KDL::Wrenches const & f_ext, ::KDL::JntArray & torques ){
        bp::override func_CartToJnt = this->get_override( "CartToJnt" );
        return func_CartToJnt( boost::ref(q), boost::ref(q_dot), boost::ref(q_dotdot), boost::ref(f_ext), boost::ref(torques) );
    }

};

struct ChainIdSolver_RNE_wrapper : KDL::ChainIdSolver_RNE, bp::wrapper< KDL::ChainIdSolver_RNE > {

    ChainIdSolver_RNE_wrapper(KDL::ChainIdSolver_RNE const & arg )
            : KDL::ChainIdSolver_RNE( arg )
            , bp::wrapper< KDL::ChainIdSolver_RNE >(){
        // copy constructor

    }

    ChainIdSolver_RNE_wrapper(::KDL::Chain const & chain, ::KDL::Vector grav )
            : KDL::ChainIdSolver_RNE( boost::ref(chain), grav )
            , bp::wrapper< KDL::ChainIdSolver_RNE >(){
        // constructor

    }

    virtual int CartToJnt( ::KDL::JntArray const & q, ::KDL::JntArray const & q_dot, ::KDL::JntArray const & q_dotdot, ::KDL::Wrenches const & f_ext, ::KDL::JntArray & torques ) {
        if( bp::override func_CartToJnt = this->get_override( "CartToJnt" ) )
            return func_CartToJnt( boost::ref(q), boost::ref(q_dot), boost::ref(q_dotdot), boost::ref(f_ext), boost::ref(torques) );
        else{
            return this->KDL::ChainIdSolver_RNE::CartToJnt( boost::ref(q), boost::ref(q_dot), boost::ref(q_dotdot), boost::ref(f_ext), boost::ref(torques) );
        }
    }

    int default_CartToJnt( ::KDL::JntArray const & q, ::KDL::JntArray const & q_dot, ::KDL::JntArray const & q_dotdot, ::KDL::Wrenches const & f_ext, ::KDL::JntArray & torques ) {
        return KDL::ChainIdSolver_RNE::CartToJnt( boost::ref(q), boost::ref(q_dot), boost::ref(q_dotdot), boost::ref(f_ext), boost::ref(torques) );
    }

};


struct ChainIkSolverPos_wrapper : KDL::ChainIkSolverPos, bp::wrapper< KDL::ChainIkSolverPos > {

    ChainIkSolverPos_wrapper()
            : KDL::ChainIkSolverPos()
            , bp::wrapper< KDL::ChainIkSolverPos >(){
        // null constructor

    }

    virtual int CartToJnt( ::KDL::JntArray const & q_init, ::KDL::Frame const & p_in, ::KDL::JntArray & q_out ){
        bp::override func_CartToJnt = this->get_override( "CartToJnt" );
        return func_CartToJnt( boost::ref(q_init), boost::ref(p_in), boost::ref(q_out) );
    }

    virtual int getError(  ) const  {
        if( bp::override func_getError = this->get_override( "getError" ) )
            return func_getError(  );
        else{
            return this->KDL::SolverI::getError(  );
        }
    }

    int default_getError(  ) const  {
        return KDL::SolverI::getError( );
    }

    virtual char const * strError( int const error ) const  {
        if( bp::override func_strError = this->get_override( "strError" ) )
            return func_strError( error );
        else{
            return this->KDL::SolverI::strError( error );
        }
    }

    char const * default_strError( int const error ) const  {
        return KDL::SolverI::strError( error );
    }

};

struct ChainIkSolverPos_LMA_wrapper : KDL::ChainIkSolverPos_LMA, bp::wrapper< KDL::ChainIkSolverPos_LMA > {

    ChainIkSolverPos_LMA_wrapper(::KDL::Chain const & _chain, ::Eigen::Matrix< double, 6, 1, 0, 6, 1 > const & _L, double _eps=1.0000000000000001E-5, int _maxiter=500, double _eps_joints=1.0000000000000001E-15 )
            : KDL::ChainIkSolverPos_LMA( boost::ref(_chain), boost::ref(_L), _eps, _maxiter, _eps_joints )
            , bp::wrapper< KDL::ChainIkSolverPos_LMA >(){
        // constructor

    }

    ChainIkSolverPos_LMA_wrapper(::KDL::Chain const & _chain, double _eps=1.0000000000000001E-5, int _maxiter=500, double _eps_joints=1.0000000000000001E-15 )
            : KDL::ChainIkSolverPos_LMA( boost::ref(_chain), _eps, _maxiter, _eps_joints )
            , bp::wrapper< KDL::ChainIkSolverPos_LMA >(){
        // constructor

    }

    virtual int CartToJnt( ::KDL::JntArray const & q_init, ::KDL::Frame const & T_base_goal, ::KDL::JntArray & q_out ) {
        if( bp::override func_CartToJnt = this->get_override( "CartToJnt" ) )
            return func_CartToJnt( boost::ref(q_init), boost::ref(T_base_goal), boost::ref(q_out) );
        else{
            return this->KDL::ChainIkSolverPos_LMA::CartToJnt( boost::ref(q_init), boost::ref(T_base_goal), boost::ref(q_out) );
        }
    }

    int default_CartToJnt( ::KDL::JntArray const & q_init, ::KDL::Frame const & T_base_goal, ::KDL::JntArray & q_out ) {
        return KDL::ChainIkSolverPos_LMA::CartToJnt( boost::ref(q_init), boost::ref(T_base_goal), boost::ref(q_out) );
    }

    virtual int getError(  ) const  {
        if( bp::override func_getError = this->get_override( "getError" ) )
            return func_getError(  );
        else{
            return this->KDL::SolverI::getError(  );
        }
    }

    int default_getError(  ) const  {
        return KDL::SolverI::getError( );
    }

    virtual char const * strError( int const error ) const  {
        if( bp::override func_strError = this->get_override( "strError" ) )
            return func_strError( error );
        else{
            return this->KDL::SolverI::strError( error );
        }
    }

    char const * default_strError( int const error ) const  {
        return KDL::SolverI::strError( error );
    }

};

struct ChainIkSolverPos_NR_wrapper : KDL::ChainIkSolverPos_NR, bp::wrapper< KDL::ChainIkSolverPos_NR > {

    ChainIkSolverPos_NR_wrapper(::KDL::Chain const & chain, ::KDL::ChainFkSolverPos & fksolver, ::KDL::ChainIkSolverVel & iksolver, unsigned int maxiter=100, double eps=9.9999999999999995E-7 )
            : KDL::ChainIkSolverPos_NR( boost::ref(chain), boost::ref(fksolver), boost::ref(iksolver), maxiter, eps )
            , bp::wrapper< KDL::ChainIkSolverPos_NR >(){
        // constructor

    }

    virtual int CartToJnt( ::KDL::JntArray const & q_init, ::KDL::Frame const & p_in, ::KDL::JntArray & q_out ) {
        if( bp::override func_CartToJnt = this->get_override( "CartToJnt" ) )
            return func_CartToJnt( boost::ref(q_init), boost::ref(p_in), boost::ref(q_out) );
        else{
            return this->KDL::ChainIkSolverPos_NR::CartToJnt( boost::ref(q_init), boost::ref(p_in), boost::ref(q_out) );
        }
    }

    int default_CartToJnt( ::KDL::JntArray const & q_init, ::KDL::Frame const & p_in, ::KDL::JntArray & q_out ) {
        return KDL::ChainIkSolverPos_NR::CartToJnt( boost::ref(q_init), boost::ref(p_in), boost::ref(q_out) );
    }

    virtual char const * strError( int const error ) const  {
        if( bp::override func_strError = this->get_override( "strError" ) )
            return func_strError( error );
        else{
            return this->KDL::ChainIkSolverPos_NR::strError( error );
        }
    }

    char const * default_strError( int const error ) const  {
        return KDL::ChainIkSolverPos_NR::strError( error );
    }

    virtual int getError(  ) const  {
        if( bp::override func_getError = this->get_override( "getError" ) )
            return func_getError(  );
        else{
            return this->KDL::SolverI::getError(  );
        }
    }

    int default_getError(  ) const  {
        return KDL::SolverI::getError( );
    }

};

struct ChainIkSolverPos_NR_JL_wrapper : KDL::ChainIkSolverPos_NR_JL, bp::wrapper< KDL::ChainIkSolverPos_NR_JL > {

    ChainIkSolverPos_NR_JL_wrapper(::KDL::Chain const & chain, ::KDL::JntArray const & q_min, ::KDL::JntArray const & q_max, ::KDL::ChainFkSolverPos & fksolver, ::KDL::ChainIkSolverVel & iksolver, unsigned int maxiter=100, double eps=9.9999999999999995E-7 )
            : KDL::ChainIkSolverPos_NR_JL( boost::ref(chain), boost::ref(q_min), boost::ref(q_max), boost::ref(fksolver), boost::ref(iksolver), maxiter, eps )
            , bp::wrapper< KDL::ChainIkSolverPos_NR_JL >(){
        // constructor

    }

    virtual int CartToJnt( ::KDL::JntArray const & q_init, ::KDL::Frame const & p_in, ::KDL::JntArray & q_out ) {
        if( bp::override func_CartToJnt = this->get_override( "CartToJnt" ) )
            return func_CartToJnt( boost::ref(q_init), boost::ref(p_in), boost::ref(q_out) );
        else{
            return this->KDL::ChainIkSolverPos_NR_JL::CartToJnt( boost::ref(q_init), boost::ref(p_in), boost::ref(q_out) );
        }
    }

    int default_CartToJnt( ::KDL::JntArray const & q_init, ::KDL::Frame const & p_in, ::KDL::JntArray & q_out ) {
        return KDL::ChainIkSolverPos_NR_JL::CartToJnt( boost::ref(q_init), boost::ref(p_in), boost::ref(q_out) );
    }

    virtual int getError(  ) const  {
        if( bp::override func_getError = this->get_override( "getError" ) )
            return func_getError(  );
        else{
            return this->KDL::SolverI::getError(  );
        }
    }

    int default_getError(  ) const  {
        return KDL::SolverI::getError( );
    }

    virtual char const * strError( int const error ) const  {
        if( bp::override func_strError = this->get_override( "strError" ) )
            return func_strError( error );
        else{
            return this->KDL::SolverI::strError( error );
        }
    }

    char const * default_strError( int const error ) const  {
        return KDL::SolverI::strError( error );
    }

};

struct ChainIkSolverVel_wrapper : KDL::ChainIkSolverVel, bp::wrapper< KDL::ChainIkSolverVel > {

    ChainIkSolverVel_wrapper()
            : KDL::ChainIkSolverVel()
            , bp::wrapper< KDL::ChainIkSolverVel >(){
        // null constructor

    }

    virtual int CartToJnt( ::KDL::JntArray const & q_in, ::KDL::Twist const & v_in, ::KDL::JntArray & qdot_out ){
        bp::override func_CartToJnt = this->get_override( "CartToJnt" );
        return func_CartToJnt( boost::ref(q_in), boost::ref(v_in), boost::ref(qdot_out) );
    }

    virtual int CartToJnt( ::KDL::JntArray const & q_init, ::KDL::FrameVel const & v_in, ::KDL::JntArrayVel & q_out ){
        bp::override func_CartToJnt = this->get_override( "CartToJnt" );
        return func_CartToJnt( boost::ref(q_init), boost::ref(v_in), boost::ref(q_out) );
    }

    virtual int getError(  ) const  {
        if( bp::override func_getError = this->get_override( "getError" ) )
            return func_getError(  );
        else{
            return this->KDL::SolverI::getError(  );
        }
    }

    int default_getError(  ) const  {
        return KDL::SolverI::getError( );
    }

    virtual char const * strError( int const error ) const  {
        if( bp::override func_strError = this->get_override( "strError" ) )
            return func_strError( error );
        else{
            return this->KDL::SolverI::strError( error );
        }
    }

    char const * default_strError( int const error ) const  {
        return KDL::SolverI::strError( error );
    }

};

struct ChainIkSolverVel_pinv_wrapper : KDL::ChainIkSolverVel_pinv, bp::wrapper< KDL::ChainIkSolverVel_pinv > {

    ChainIkSolverVel_pinv_wrapper(::KDL::Chain const & chain, double eps=1.0000000000000001E-5, int maxiter=150 )
            : KDL::ChainIkSolverVel_pinv( boost::ref(chain), eps, maxiter )
            , bp::wrapper< KDL::ChainIkSolverVel_pinv >(){
        // constructor

    }

    virtual int CartToJnt( ::KDL::JntArray const & q_in, ::KDL::Twist const & v_in, ::KDL::JntArray & qdot_out ) {
        if( bp::override func_CartToJnt = this->get_override( "CartToJnt" ) )
            return func_CartToJnt( boost::ref(q_in), boost::ref(v_in), boost::ref(qdot_out) );
        else{
            return this->KDL::ChainIkSolverVel_pinv::CartToJnt( boost::ref(q_in), boost::ref(v_in), boost::ref(qdot_out) );
        }
    }

    int default_CartToJnt( ::KDL::JntArray const & q_in, ::KDL::Twist const & v_in, ::KDL::JntArray & qdot_out ) {
        return KDL::ChainIkSolverVel_pinv::CartToJnt( boost::ref(q_in), boost::ref(v_in), boost::ref(qdot_out) );
    }

    virtual int CartToJnt( ::KDL::JntArray const & q_init, ::KDL::FrameVel const & v_in, ::KDL::JntArrayVel & q_out ) {
        if( bp::override func_CartToJnt = this->get_override( "CartToJnt" ) )
            return func_CartToJnt( boost::ref(q_init), boost::ref(v_in), boost::ref(q_out) );
        else{
            return this->KDL::ChainIkSolverVel_pinv::CartToJnt( boost::ref(q_init), boost::ref(v_in), boost::ref(q_out) );
        }
    }

    int default_CartToJnt( ::KDL::JntArray const & q_init, ::KDL::FrameVel const & v_in, ::KDL::JntArrayVel & q_out ) {
        return KDL::ChainIkSolverVel_pinv::CartToJnt( boost::ref(q_init), boost::ref(v_in), boost::ref(q_out) );
    }

    virtual char const * strError( int const error ) const  {
        if( bp::override func_strError = this->get_override( "strError" ) )
            return func_strError( error );
        else{
            return this->KDL::ChainIkSolverVel_pinv::strError( error );
        }
    }

    char const * default_strError( int const error ) const  {
        return KDL::ChainIkSolverVel_pinv::strError( error );
    }

    virtual int getError(  ) const  {
        if( bp::override func_getError = this->get_override( "getError" ) )
            return func_getError(  );
        else{
            return this->KDL::SolverI::getError(  );
        }
    }

    int default_getError(  ) const  {
        return KDL::SolverI::getError( );
    }

};

struct ChainIkSolverVel_pinv_givens_wrapper : KDL::ChainIkSolverVel_pinv_givens, bp::wrapper< KDL::ChainIkSolverVel_pinv_givens > {

    ChainIkSolverVel_pinv_givens_wrapper(::KDL::Chain const & chain )
            : KDL::ChainIkSolverVel_pinv_givens( boost::ref(chain) )
            , bp::wrapper< KDL::ChainIkSolverVel_pinv_givens >(){
        // constructor

    }

    virtual int CartToJnt( ::KDL::JntArray const & q_in, ::KDL::Twist const & v_in, ::KDL::JntArray & qdot_out ) {
        if( bp::override func_CartToJnt = this->get_override( "CartToJnt" ) )
            return func_CartToJnt( boost::ref(q_in), boost::ref(v_in), boost::ref(qdot_out) );
        else{
            return this->KDL::ChainIkSolverVel_pinv_givens::CartToJnt( boost::ref(q_in), boost::ref(v_in), boost::ref(qdot_out) );
        }
    }

    int default_CartToJnt( ::KDL::JntArray const & q_in, ::KDL::Twist const & v_in, ::KDL::JntArray & qdot_out ) {
        return KDL::ChainIkSolverVel_pinv_givens::CartToJnt( boost::ref(q_in), boost::ref(v_in), boost::ref(qdot_out) );
    }

    virtual int CartToJnt( ::KDL::JntArray const & q_init, ::KDL::FrameVel const & v_in, ::KDL::JntArrayVel & q_out ) {
        if( bp::override func_CartToJnt = this->get_override( "CartToJnt" ) )
            return func_CartToJnt( boost::ref(q_init), boost::ref(v_in), boost::ref(q_out) );
        else{
            return this->KDL::ChainIkSolverVel_pinv_givens::CartToJnt( boost::ref(q_init), boost::ref(v_in), boost::ref(q_out) );
        }
    }

    int default_CartToJnt( ::KDL::JntArray const & q_init, ::KDL::FrameVel const & v_in, ::KDL::JntArrayVel & q_out ) {
        return KDL::ChainIkSolverVel_pinv_givens::CartToJnt( boost::ref(q_init), boost::ref(v_in), boost::ref(q_out) );
    }

    virtual int getError(  ) const  {
        if( bp::override func_getError = this->get_override( "getError" ) )
            return func_getError(  );
        else{
            return this->KDL::SolverI::getError(  );
        }
    }

    int default_getError(  ) const  {
        return KDL::SolverI::getError( );
    }

    virtual char const * strError( int const error ) const  {
        if( bp::override func_strError = this->get_override( "strError" ) )
            return func_strError( error );
        else{
            return this->KDL::SolverI::strError( error );
        }
    }

    char const * default_strError( int const error ) const  {
        return KDL::SolverI::strError( error );
    }

};

struct ChainIkSolverVel_pinv_nso_wrapper : KDL::ChainIkSolverVel_pinv_nso, bp::wrapper< KDL::ChainIkSolverVel_pinv_nso > {

    ChainIkSolverVel_pinv_nso_wrapper(::KDL::Chain const & chain, ::KDL::JntArray opt_pos, ::KDL::JntArray weights, double eps=1.0000000000000001E-5, int maxiter=150, double alpha=0.25 )
            : KDL::ChainIkSolverVel_pinv_nso( boost::ref(chain), opt_pos, weights, eps, maxiter, alpha )
            , bp::wrapper< KDL::ChainIkSolverVel_pinv_nso >(){
        // constructor

    }

    ChainIkSolverVel_pinv_nso_wrapper(::KDL::Chain const & chain, double eps=1.0000000000000001E-5, int maxiter=150, double alpha=0.25 )
            : KDL::ChainIkSolverVel_pinv_nso( boost::ref(chain), eps, maxiter, alpha )
            , bp::wrapper< KDL::ChainIkSolverVel_pinv_nso >(){
        // constructor

    }

    virtual int CartToJnt( ::KDL::JntArray const & q_in, ::KDL::Twist const & v_in, ::KDL::JntArray & qdot_out ) {
        if( bp::override func_CartToJnt = this->get_override( "CartToJnt" ) )
            return func_CartToJnt( boost::ref(q_in), boost::ref(v_in), boost::ref(qdot_out) );
        else{
            return this->KDL::ChainIkSolverVel_pinv_nso::CartToJnt( boost::ref(q_in), boost::ref(v_in), boost::ref(qdot_out) );
        }
    }

    int default_CartToJnt( ::KDL::JntArray const & q_in, ::KDL::Twist const & v_in, ::KDL::JntArray & qdot_out ) {
        return KDL::ChainIkSolverVel_pinv_nso::CartToJnt( boost::ref(q_in), boost::ref(v_in), boost::ref(qdot_out) );
    }

    virtual int CartToJnt( ::KDL::JntArray const & q_init, ::KDL::FrameVel const & v_in, ::KDL::JntArrayVel & q_out ) {
        if( bp::override func_CartToJnt = this->get_override( "CartToJnt" ) )
            return func_CartToJnt( boost::ref(q_init), boost::ref(v_in), boost::ref(q_out) );
        else{
            return this->KDL::ChainIkSolverVel_pinv_nso::CartToJnt( boost::ref(q_init), boost::ref(v_in), boost::ref(q_out) );
        }
    }

    int default_CartToJnt( ::KDL::JntArray const & q_init, ::KDL::FrameVel const & v_in, ::KDL::JntArrayVel & q_out ) {
        return KDL::ChainIkSolverVel_pinv_nso::CartToJnt( boost::ref(q_init), boost::ref(v_in), boost::ref(q_out) );
    }

    virtual int setAlpha( double const alpha ) {
        if( bp::override func_setAlpha = this->get_override( "setAlpha" ) )
            return func_setAlpha( alpha );
        else{
            return this->KDL::ChainIkSolverVel_pinv_nso::setAlpha( alpha );
        }
    }

    int default_setAlpha( double const alpha ) {
        return KDL::ChainIkSolverVel_pinv_nso::setAlpha( alpha );
    }

    virtual int setOptPos( ::KDL::JntArray const & opt_pos ) {
        if( bp::override func_setOptPos = this->get_override( "setOptPos" ) )
            return func_setOptPos( boost::ref(opt_pos) );
        else{
            return this->KDL::ChainIkSolverVel_pinv_nso::setOptPos( boost::ref(opt_pos) );
        }
    }

    int default_setOptPos( ::KDL::JntArray const & opt_pos ) {
        return KDL::ChainIkSolverVel_pinv_nso::setOptPos( boost::ref(opt_pos) );
    }

    virtual int setWeights( ::KDL::JntArray const & weights ) {
        if( bp::override func_setWeights = this->get_override( "setWeights" ) )
            return func_setWeights( boost::ref(weights) );
        else{
            return this->KDL::ChainIkSolverVel_pinv_nso::setWeights( boost::ref(weights) );
        }
    }

    int default_setWeights( ::KDL::JntArray const & weights ) {
        return KDL::ChainIkSolverVel_pinv_nso::setWeights( boost::ref(weights) );
    }

    virtual int getError(  ) const  {
        if( bp::override func_getError = this->get_override( "getError" ) )
            return func_getError(  );
        else{
            return this->KDL::SolverI::getError(  );
        }
    }

    int default_getError(  ) const  {
        return KDL::SolverI::getError( );
    }

    virtual char const * strError( int const error ) const  {
        if( bp::override func_strError = this->get_override( "strError" ) )
            return func_strError( error );
        else{
            return this->KDL::SolverI::strError( error );
        }
    }

    char const * default_strError( int const error ) const  {
        return KDL::SolverI::strError( error );
    }

};

struct ChainIkSolverVel_wdls_wrapper : KDL::ChainIkSolverVel_wdls, bp::wrapper< KDL::ChainIkSolverVel_wdls > {

    ChainIkSolverVel_wdls_wrapper(::KDL::Chain const & chain, double eps=1.0000000000000001E-5, int maxiter=150 )
            : KDL::ChainIkSolverVel_wdls( boost::ref(chain), eps, maxiter )
            , bp::wrapper< KDL::ChainIkSolverVel_wdls >(){
        // constructor

    }

    virtual int CartToJnt( ::KDL::JntArray const & q_in, ::KDL::Twist const & v_in, ::KDL::JntArray & qdot_out ) {
        if( bp::override func_CartToJnt = this->get_override( "CartToJnt" ) )
            return func_CartToJnt( boost::ref(q_in), boost::ref(v_in), boost::ref(qdot_out) );
        else{
            return this->KDL::ChainIkSolverVel_wdls::CartToJnt( boost::ref(q_in), boost::ref(v_in), boost::ref(qdot_out) );
        }
    }

    int default_CartToJnt( ::KDL::JntArray const & q_in, ::KDL::Twist const & v_in, ::KDL::JntArray & qdot_out ) {
        return KDL::ChainIkSolverVel_wdls::CartToJnt( boost::ref(q_in), boost::ref(v_in), boost::ref(qdot_out) );
    }

    virtual int CartToJnt( ::KDL::JntArray const & q_init, ::KDL::FrameVel const & v_in, ::KDL::JntArrayVel & q_out ) {
        if( bp::override func_CartToJnt = this->get_override( "CartToJnt" ) )
            return func_CartToJnt( boost::ref(q_init), boost::ref(v_in), boost::ref(q_out) );
        else{
            return this->KDL::ChainIkSolverVel_wdls::CartToJnt( boost::ref(q_init), boost::ref(v_in), boost::ref(q_out) );
        }
    }

    int default_CartToJnt( ::KDL::JntArray const & q_init, ::KDL::FrameVel const & v_in, ::KDL::JntArrayVel & q_out ) {
        return KDL::ChainIkSolverVel_wdls::CartToJnt( boost::ref(q_init), boost::ref(v_in), boost::ref(q_out) );
    }

    virtual char const * strError( int const error ) const  {
        if( bp::override func_strError = this->get_override( "strError" ) )
            return func_strError( error );
        else{
            return this->KDL::ChainIkSolverVel_wdls::strError( error );
        }
    }

    char const * default_strError( int const error ) const  {
        return KDL::ChainIkSolverVel_wdls::strError( error );
    }

    virtual int getError(  ) const  {
        if( bp::override func_getError = this->get_override( "getError" ) )
            return func_getError(  );
        else{
            return this->KDL::SolverI::getError(  );
        }
    }

    int default_getError(  ) const  {
        return KDL::SolverI::getError( );
    }

};

struct ChainJntToJacSolver_wrapper : KDL::ChainJntToJacSolver, bp::wrapper< KDL::ChainJntToJacSolver > {

    ChainJntToJacSolver_wrapper(::KDL::Chain const & chain )
            : KDL::ChainJntToJacSolver( boost::ref(chain) )
            , bp::wrapper< KDL::ChainJntToJacSolver >(){
        // constructor

    }

    virtual int JntToJac( ::KDL::JntArray const & q_in, ::KDL::Jacobian & jac, int segmentNR=-1 ) {
        if( bp::override func_JntToJac = this->get_override( "JntToJac" ) )
            return func_JntToJac( boost::ref(q_in), boost::ref(jac), segmentNR );
        else{
            return this->KDL::ChainJntToJacSolver::JntToJac( boost::ref(q_in), boost::ref(jac), segmentNR );
        }
    }

    int default_JntToJac( ::KDL::JntArray const & q_in, ::KDL::Jacobian & jac, int segmentNR=-1 ) {
        return KDL::ChainJntToJacSolver::JntToJac( boost::ref(q_in), boost::ref(jac), segmentNR );
    }

    virtual char const * strError( int const error ) const  {
        if( bp::override func_strError = this->get_override( "strError" ) )
            return func_strError( error );
        else{
            return this->KDL::ChainJntToJacSolver::strError( error );
        }
    }

    char const * default_strError( int const error ) const  {
        return KDL::ChainJntToJacSolver::strError( error );
    }

    virtual int getError(  ) const  {
        if( bp::override func_getError = this->get_override( "getError" ) )
            return func_getError(  );
        else{
            return this->KDL::SolverI::getError(  );
        }
    }

    int default_getError(  ) const  {
        return KDL::SolverI::getError( );
    }

};




struct Vector_wrapper : KDL::Vector, bp::wrapper< KDL::Vector > {

    Vector_wrapper( )
            : KDL::Vector( )
            , bp::wrapper< KDL::Vector >(){
        // null constructor

    }

    Vector_wrapper(double x, double y, double z )
            : KDL::Vector( x, y, z )
            , bp::wrapper< KDL::Vector >(){
        // constructor

    }

    Vector_wrapper(::KDL::Vector const & arg )
            : KDL::Vector( boost::ref(arg) )
            , bp::wrapper< KDL::Vector >(){
        // copy constructor

    }


};





struct Rotation_wrapper : KDL::Rotation, bp::wrapper< KDL::Rotation > {

    Rotation_wrapper(KDL::Rotation const & arg )
            : KDL::Rotation( arg )
            , bp::wrapper< KDL::Rotation >(){
        // copy constructor

    }

    Rotation_wrapper( )
            : KDL::Rotation( )
            , bp::wrapper< KDL::Rotation >(){
        // null constructor

    }

    Rotation_wrapper(double Xx, double Yx, double Zx, double Xy, double Yy, double Zy, double Xz, double Yz, double Zz )
            : KDL::Rotation( Xx, Yx, Zx, Xy, Yy, Zy, Xz, Yz, Zz )
            , bp::wrapper< KDL::Rotation >(){
        // constructor

    }

    Rotation_wrapper(::KDL::Vector const & x, ::KDL::Vector const & y, ::KDL::Vector const & z )
            : KDL::Rotation( boost::ref(x), boost::ref(y), boost::ref(z) )
            , bp::wrapper< KDL::Rotation >(){
        // constructor

    }

//    static pyplusplus::containers::static_sized::array_1_t< double, 9>
//    pyplusplus_data_wrapper( ::KDL::Rotation & inst ){
//        return pyplusplus::containers::static_sized::array_1_t< double, 9>( inst.data );
//    }

};









int    (KDL::ChainFkSolverPos_recursive::*JntToCart_1)(const KDL::JntArray &, KDL::Frame &, int)              = &KDL::ChainFkSolverPos_recursive::JntToCart;

void JntToCart(const KDL::Chain& chain, const KDL::JntArray& q_in, KDL::Frame& p_out)
{
    KDL::ChainFkSolverPos_recursive fk_solver = KDL::ChainFkSolverPos_recursive(chain);
	fk_solver.JntToCart(q_in, p_out);
}


int CartToJnt(const KDL::Chain& chain, const KDL::JntArray& q_init, KDL::Frame& F_dest, KDL::JntArray& q)
{
    KDL::ChainFkSolverPos_recursive fksolver(chain);//Forward position solver
    KDL::ChainIkSolverVel_pinv iksolverv(chain);//Inverse velocity solver
    KDL::ChainIkSolverPos_NR iksolverpos(chain,fksolver,iksolverv,100,1e-6);//Maximum 100 iterations, stop at accuracy 1e-6
  
	int ret = iksolverpos.CartToJnt(q_init,F_dest,q);
	return ret;

}



// expose index operator

void Vector_setitem(KDL::Vector& v, int index, double value)
{
    if (index >= 0 && index < 3)
    {
        v[index] = value;
    }
    else {
        PyErr_SetString(PyExc_IndexError, "index out of range");
        bp::throw_error_already_set();
    }
}

double Vector_getitem(KDL::Vector& v, int index)
{
    if (index >= 0 && index < 3)
    {
        return v[index];
    }
    else {
        PyErr_SetString(PyExc_IndexError, "index out of range");
        bp::throw_error_already_set();
    }
}


void JntArray_setitem(KDL::JntArray& q, int index, double value)
{
    if (index >= 0 && index < (q.rows())*(q.columns()))
    {
        q(index) = value;
    }
    else {
        PyErr_SetString(PyExc_IndexError, "index out of range");
        bp::throw_error_already_set();
    }
}

double JntArray_getitem(KDL::JntArray& q, int index)
{
    if (index >= 0 && index < (q.rows())*(q.columns()))
    {
        return q(index);
    }
    else {
        PyErr_SetString(PyExc_IndexError, "index out of range");
        bp::throw_error_already_set();
    }
}


void Rotation_setitem(KDL::Rotation& r, int index, double value)
{

    if (index >= 0 && index < 9)
    {
    	int i = int(index/3);
		int j= index - 3*i;
        r(i, j) = value;
    }
    else {
        PyErr_SetString(PyExc_IndexError, "index out of range");
        bp::throw_error_already_set();
    }
}

double Rotation_getitem(KDL::Rotation& r, int index)
{
   if (index >= 0 && index < 9)
    {
    	int i = int(index/3);
		int j= index - 3*i;
        return r(i, j);
    }
   else {
       PyErr_SetString(PyExc_IndexError, "index out of range");
       bp::throw_error_already_set();
   }
}

void Jacobian_setitem(KDL::Jacobian& jac, int index, double value)
{

    if (index >= 0 && index < jac.rows()*jac.columns())
    {
    	int i = int(index/jac.columns());
		int j= index - jac.columns()*i;
        jac(i, j) = value;
    }
    else {
        PyErr_SetString(PyExc_IndexError, "index out of range");
        bp::throw_error_already_set();
    }
}

double Jacobian_getitem(KDL::Jacobian& jac, int index)
{
   if (index >= 0 && index < jac.rows()*jac.columns())
    {
    	int i = int(index/jac.columns());
		int j= index - jac.columns()*i;
        return jac(i, j);
    }
   else {
       PyErr_SetString(PyExc_IndexError, "index out of range");
       bp::throw_error_already_set();
   }
}



//Wrapper
BOOST_PYTHON_MODULE(PyKDL)
{
	Py_Initialize();


    { //::KDL::Frame
        typedef bp::class_< KDL::Frame > Frame_exposer_t;
        Frame_exposer_t Frame_exposer = Frame_exposer_t( "Frame", bp::init< KDL::Rotation const &, KDL::Vector const & >(( bp::arg("R"), bp::arg("V") )) );
        bp::scope Frame_scope( Frame_exposer );
        Frame_exposer.def( bp::init< KDL::Vector const & >(( bp::arg("V") )) );
        bp::implicitly_convertible< KDL::Vector const &, KDL::Frame >();
        Frame_exposer.def( bp::init< KDL::Rotation const & >(( bp::arg("R") )) );
        bp::implicitly_convertible< KDL::Rotation const &, KDL::Frame >();
        Frame_exposer.def( bp::init< >() );
        Frame_exposer.def( bp::init< KDL::Frame const & >(( bp::arg("arg") )) );
        { //::KDL::Frame::DH

            typedef ::KDL::Frame ( *DH_function_type )( double,double,double,double );

            Frame_exposer.def(
                    "DH"
                    , DH_function_type( &::KDL::Frame::DH )
                    , ( bp::arg("a"), bp::arg("alpha"), bp::arg("d"), bp::arg("theta") ) );

        }
        { //::KDL::Frame::DH_Craig1989

            typedef ::KDL::Frame ( *DH_Craig1989_function_type )( double,double,double,double );

            Frame_exposer.def(
                    "DH_Craig1989"
                    , DH_Craig1989_function_type( &::KDL::Frame::DH_Craig1989 )
                    , ( bp::arg("a"), bp::arg("alpha"), bp::arg("d"), bp::arg("theta") ) );

        }
        { //::KDL::Frame::Identity

            typedef ::KDL::Frame ( *Identity_function_type )(  );

            Frame_exposer.def(
                    "Identity"
                    , Identity_function_type( &::KDL::Frame::Identity ) );

        }
         { //::KDL::Frame::Integrate

             typedef void ( ::KDL::Frame::*Integrate_function_type)( ::KDL::Twist const &,double ) ;

             Frame_exposer.def(
                 "Integrate"
                 , Integrate_function_type( &::KDL::Frame::Integrate )
                 , ( bp::arg("t_this"), bp::arg("frequency") ) );

         }
        { //::KDL::Frame::Inverse

            typedef ::KDL::Frame ( ::KDL::Frame::*Inverse_function_type)(  ) const;

            Frame_exposer.def(
                    "Inverse"
                    , Inverse_function_type( &::KDL::Frame::Inverse ) );

        }
        { //::KDL::Frame::Inverse

            typedef ::KDL::Vector ( ::KDL::Frame::*Inverse_function_type)( ::KDL::Vector const & ) const;

            Frame_exposer.def(
                    "Inverse"
                    , Inverse_function_type( &::KDL::Frame::Inverse )
                    , ( bp::arg("arg") ) );

        }
        { //::KDL::Frame::Inverse

            typedef ::KDL::Wrench ( ::KDL::Frame::*Inverse_function_type)( ::KDL::Wrench const & ) const;

            Frame_exposer.def(
                    "Inverse"
                    , Inverse_function_type( &::KDL::Frame::Inverse )
                    , ( bp::arg("arg") ) );

        }
        { //::KDL::Frame::Inverse

            typedef ::KDL::Twist ( ::KDL::Frame::*Inverse_function_type)( ::KDL::Twist const & ) const;

            Frame_exposer.def(
                    "Inverse"
                    , Inverse_function_type( &::KDL::Frame::Inverse )
                    , ( bp::arg("arg") ) );

        }

         { //::KDL::Frame::Make4x4

             typedef void ( ::KDL::Frame::*Make4x4_function_type)( double * ) ;

             Frame_exposer.def(
                 "Make4x4"
                 , Make4x4_function_type( &::KDL::Frame::Make4x4 )
                 , ( bp::arg("d") ) );

         }

         { //::KDL::Frame::operator()

             typedef double ( ::KDL::Frame::*__call___function_type)( int,int ) ;

             Frame_exposer.def(
                 "__call__"
                 , __call___function_type( &::KDL::Frame::operator() )
                 , ( bp::arg("i"), bp::arg("j") ) );

         }
         { //::KDL::Frame::operator()

             typedef double ( ::KDL::Frame::*__call___function_type)( int,int ) const;

             Frame_exposer.def(
                 "__call__"
                 , __call___function_type( &::KDL::Frame::operator() )
                 , ( bp::arg("i"), bp::arg("j") ) );

         }
        Frame_exposer.def( bp::self * bp::other< KDL::Vector >() );
        Frame_exposer.def( bp::self * bp::other< KDL::Wrench >() );
        Frame_exposer.def( bp::self * bp::other< KDL::Twist >() );
        { //::KDL::Frame::operator=

            typedef ::KDL::Frame & ( ::KDL::Frame::*assign_function_type)( ::KDL::Frame const & ) ;

            Frame_exposer.def(
                    "assign"
                    , assign_function_type( &::KDL::Frame::operator= )
                    , ( bp::arg("arg") )
                    , bp::return_self< >() );

        }
        Frame_exposer.def_readwrite( "M", &KDL::Frame::M );
        Frame_exposer.def_readwrite( "p", &KDL::Frame::p );
        Frame_exposer.staticmethod( "DH" );
        Frame_exposer.staticmethod( "DH_Craig1989" );
        Frame_exposer.staticmethod( "Identity" );
        Frame_exposer.def( bp::self != bp::self );
        Frame_exposer.def( bp::self * bp::other< KDL::RigidBodyInertia >() );
        Frame_exposer.def( bp::self * bp::self );
        Frame_exposer.def( bp::self == bp::self );
    }

    { //::KDL::Frame2
        typedef bp::class_< KDL::Frame2 > Frame2_exposer_t;
        Frame2_exposer_t Frame2_exposer = Frame2_exposer_t( "Frame2", bp::init< KDL::Rotation2 const &, KDL::Vector2 const & >(( bp::arg("R"), bp::arg("V") )) );
        bp::scope Frame2_scope( Frame2_exposer );
        Frame2_exposer.def( bp::init< KDL::Vector2 const & >(( bp::arg("V") )) );
        bp::implicitly_convertible< KDL::Vector2 const &, KDL::Frame2 >();
        Frame2_exposer.def( bp::init< KDL::Rotation2 const & >(( bp::arg("R") )) );
        bp::implicitly_convertible< KDL::Rotation2 const &, KDL::Frame2 >();
        Frame2_exposer.def( bp::init< >() );
        Frame2_exposer.def( bp::init< KDL::Frame2 const & >(( bp::arg("arg") )) );
        { //::KDL::Frame2::Identity

            typedef ::KDL::Frame2 ( *Identity_function_type )(  );

            Frame2_exposer.def(
                    "Identity"
                    , Identity_function_type( &::KDL::Frame2::Identity ) );

        }
        { //::KDL::Frame2::Inverse

            typedef ::KDL::Frame2 ( ::KDL::Frame2::*Inverse_function_type)(  ) const;

            Frame2_exposer.def(
                    "Inverse"
                    , Inverse_function_type( &::KDL::Frame2::Inverse ) );

        }
        { //::KDL::Frame2::Inverse

            typedef ::KDL::Vector2 ( ::KDL::Frame2::*Inverse_function_type)( ::KDL::Vector2 const & ) const;

            Frame2_exposer.def(
                    "Inverse"
                    , Inverse_function_type( &::KDL::Frame2::Inverse )
                    , ( bp::arg("arg") ) );

        }
        { //::KDL::Frame2::SetIdentity

            typedef void ( ::KDL::Frame2::*SetIdentity_function_type)(  ) ;

            Frame2_exposer.def(
                    "SetIdentity"
                    , SetIdentity_function_type( &::KDL::Frame2::SetIdentity ) );

        }
        { //::KDL::Frame2::SetInverse

            typedef void ( ::KDL::Frame2::*SetInverse_function_type)(  ) ;

            Frame2_exposer.def(
                    "SetInverse"
                    , SetInverse_function_type( &::KDL::Frame2::SetInverse ) );

        }
         { //::KDL::Frame2::operator()

             typedef double ( ::KDL::Frame2::*__call___function_type)( int,int ) ;

             Frame2_exposer.def(
                 "__call__"
                 , __call___function_type( &::KDL::Frame2::operator() )
                 , ( bp::arg("i"), bp::arg("j") ) );

         }
         { //::KDL::Frame2::operator()

             typedef double ( ::KDL::Frame2::*__call___function_type)( int,int ) const;

             Frame2_exposer.def(
                 "__call__"
                 , __call___function_type( &::KDL::Frame2::operator() )
                 , ( bp::arg("i"), bp::arg("j") ) );

         }
        Frame2_exposer.def( bp::self * bp::other< KDL::Vector2 >() );
        { //::KDL::Frame2::operator=

            typedef ::KDL::Frame2 & ( ::KDL::Frame2::*assign_function_type)( ::KDL::Frame2 const & ) ;

            Frame2_exposer.def(
                    "assign"
                    , assign_function_type( &::KDL::Frame2::operator= )
                    , ( bp::arg("arg") )
                    , bp::return_self< >() );

        }
        Frame2_exposer.def_readwrite( "M", &KDL::Frame2::M );
        Frame2_exposer.def_readwrite( "p", &KDL::Frame2::p );
        Frame2_exposer.staticmethod( "Identity" );
        Frame2_exposer.def( bp::self * bp::self );
    }


    { //::KDL::FrameAcc
        typedef bp::class_< KDL::FrameAcc > FrameAcc_exposer_t;
        FrameAcc_exposer_t FrameAcc_exposer = FrameAcc_exposer_t( "FrameAcc", bp::init< >() );
        bp::scope FrameAcc_scope( FrameAcc_exposer );
        FrameAcc_exposer.def( bp::init< KDL::Frame const & >(( bp::arg("_T") )) );
        bp::implicitly_convertible< KDL::Frame const &, KDL::FrameAcc >();
        FrameAcc_exposer.def( bp::init< KDL::Frame const &, KDL::Twist const &, KDL::Twist const & >(( bp::arg("_T"), bp::arg("_t"), bp::arg("_dt") )) );
        FrameAcc_exposer.def( bp::init< KDL::RotationAcc const &, KDL::VectorAcc const & >(( bp::arg("_M"), bp::arg("_p") )) );
        { //::KDL::FrameAcc::GetAccTwist

            typedef ::KDL::Twist ( ::KDL::FrameAcc::*GetAccTwist_function_type)(  ) const;

            FrameAcc_exposer.def(
                    "GetAccTwist"
                    , GetAccTwist_function_type( &::KDL::FrameAcc::GetAccTwist ) );

        }
        { //::KDL::FrameAcc::GetFrame

            typedef ::KDL::Frame ( ::KDL::FrameAcc::*GetFrame_function_type)(  ) const;

            FrameAcc_exposer.def(
                    "GetFrame"
                    , GetFrame_function_type( &::KDL::FrameAcc::GetFrame ) );

        }
        { //::KDL::FrameAcc::GetTwist

            typedef ::KDL::Twist ( ::KDL::FrameAcc::*GetTwist_function_type)(  ) const;

            FrameAcc_exposer.def(
                    "GetTwist"
                    , GetTwist_function_type( &::KDL::FrameAcc::GetTwist ) );

        }
        { //::KDL::FrameAcc::Identity

            typedef ::KDL::FrameAcc ( *Identity_function_type )(  );

            FrameAcc_exposer.def(
                    "Identity"
                    , Identity_function_type( &::KDL::FrameAcc::Identity ) );

        }
        { //::KDL::FrameAcc::Inverse

            typedef ::KDL::FrameAcc ( ::KDL::FrameAcc::*Inverse_function_type)(  ) const;

            FrameAcc_exposer.def(
                    "Inverse"
                    , Inverse_function_type( &::KDL::FrameAcc::Inverse ) );

        }
        { //::KDL::FrameAcc::Inverse

            typedef ::KDL::VectorAcc ( ::KDL::FrameAcc::*Inverse_function_type)( ::KDL::VectorAcc const & ) const;

            FrameAcc_exposer.def(
                    "Inverse"
                    , Inverse_function_type( &::KDL::FrameAcc::Inverse )
                    , ( bp::arg("arg") ) );

        }
        { //::KDL::FrameAcc::Inverse

            typedef ::KDL::VectorAcc ( ::KDL::FrameAcc::*Inverse_function_type)( ::KDL::Vector const & ) const;

            FrameAcc_exposer.def(
                    "Inverse"
                    , Inverse_function_type( &::KDL::FrameAcc::Inverse )
                    , ( bp::arg("arg") ) );

        }
        { //::KDL::FrameAcc::Inverse

            typedef ::KDL::TwistAcc ( ::KDL::FrameAcc::*Inverse_function_type)( ::KDL::TwistAcc const & ) const;

            FrameAcc_exposer.def(
                    "Inverse"
                    , Inverse_function_type( &::KDL::FrameAcc::Inverse )
                    , ( bp::arg("arg") ) );

        }
        { //::KDL::FrameAcc::Inverse

            typedef ::KDL::TwistAcc ( ::KDL::FrameAcc::*Inverse_function_type)( ::KDL::Twist const & ) const;

            FrameAcc_exposer.def(
                    "Inverse"
                    , Inverse_function_type( &::KDL::FrameAcc::Inverse )
                    , ( bp::arg("arg") ) );

        }
        FrameAcc_exposer.def( bp::self * bp::other< KDL::VectorAcc >() );
        FrameAcc_exposer.def( bp::self * bp::other< KDL::Vector >() );
        FrameAcc_exposer.def( bp::self * bp::other< KDL::TwistAcc >() );
        FrameAcc_exposer.def( bp::self * bp::other< KDL::Twist >() );
        { //::KDL::FrameAcc::operator=

            typedef ::KDL::FrameAcc & ( ::KDL::FrameAcc::*assign_function_type)( ::KDL::FrameAcc const & ) ;

            FrameAcc_exposer.def(
                    "assign"
                    , assign_function_type( &::KDL::FrameAcc::operator= )
                    , ( bp::arg("arg") )
                    , bp::return_self< >() );

        }
        { //::KDL::FrameAcc::operator=

            typedef ::KDL::FrameAcc & ( ::KDL::FrameAcc::*assign_function_type)( ::KDL::Frame const & ) ;

            FrameAcc_exposer.def(
                    "assign"
                    , assign_function_type( &::KDL::FrameAcc::operator= )
                    , ( bp::arg("arg") )
                    , bp::return_self< >() );

        }
        FrameAcc_exposer.def_readwrite( "M", &KDL::FrameAcc::M );
        FrameAcc_exposer.def_readwrite( "p", &KDL::FrameAcc::p );
        FrameAcc_exposer.staticmethod( "Identity" );
        FrameAcc_exposer.def( bp::self * bp::self );
        FrameAcc_exposer.def( bp::self * bp::other< KDL::Frame >() );
    }

    { //::KDL::FrameVel
        typedef bp::class_< KDL::FrameVel > FrameVel_exposer_t;
        FrameVel_exposer_t FrameVel_exposer = FrameVel_exposer_t( "FrameVel", bp::init< >() );
        bp::scope FrameVel_scope( FrameVel_exposer );
        FrameVel_exposer.def( bp::init< KDL::Frame const & >(( bp::arg("_T") )) );
        bp::implicitly_convertible< KDL::Frame const &, KDL::FrameVel >();
        FrameVel_exposer.def( bp::init< KDL::Frame const &, KDL::Twist const & >(( bp::arg("_T"), bp::arg("_t") )) );
        FrameVel_exposer.def( bp::init< KDL::RotationVel const &, KDL::VectorVel const & >(( bp::arg("_M"), bp::arg("_p") )) );
        { //::KDL::FrameVel::GetFrame

            typedef ::KDL::Frame ( ::KDL::FrameVel::*GetFrame_function_type)(  ) const;

            FrameVel_exposer.def(
                    "GetFrame"
                    , GetFrame_function_type( &::KDL::FrameVel::GetFrame ) );

        }
        { //::KDL::FrameVel::GetTwist

            typedef ::KDL::Twist ( ::KDL::FrameVel::*GetTwist_function_type)(  ) const;

            FrameVel_exposer.def(
                    "GetTwist"
                    , GetTwist_function_type( &::KDL::FrameVel::GetTwist ) );

        }
        { //::KDL::FrameVel::Identity

            typedef ::KDL::FrameVel ( *Identity_function_type )(  );

            FrameVel_exposer.def(
                    "Identity"
                    , Identity_function_type( &::KDL::FrameVel::Identity ) );

        }
        { //::KDL::FrameVel::Inverse

            typedef ::KDL::FrameVel ( ::KDL::FrameVel::*Inverse_function_type)(  ) const;

            FrameVel_exposer.def(
                    "Inverse"
                    , Inverse_function_type( &::KDL::FrameVel::Inverse ) );

        }
        { //::KDL::FrameVel::Inverse

            typedef ::KDL::VectorVel ( ::KDL::FrameVel::*Inverse_function_type)( ::KDL::VectorVel const & ) const;

            FrameVel_exposer.def(
                    "Inverse"
                    , Inverse_function_type( &::KDL::FrameVel::Inverse )
                    , ( bp::arg("arg") ) );

        }
        { //::KDL::FrameVel::Inverse

            typedef ::KDL::VectorVel ( ::KDL::FrameVel::*Inverse_function_type)( ::KDL::Vector const & ) const;

            FrameVel_exposer.def(
                    "Inverse"
                    , Inverse_function_type( &::KDL::FrameVel::Inverse )
                    , ( bp::arg("arg") ) );

        }
        { //::KDL::FrameVel::Inverse

            typedef ::KDL::TwistVel ( ::KDL::FrameVel::*Inverse_function_type)( ::KDL::TwistVel const & ) const;

            FrameVel_exposer.def(
                    "Inverse"
                    , Inverse_function_type( &::KDL::FrameVel::Inverse )
                    , ( bp::arg("arg") ) );

        }
        { //::KDL::FrameVel::Inverse

            typedef ::KDL::TwistVel ( ::KDL::FrameVel::*Inverse_function_type)( ::KDL::Twist const & ) const;

            FrameVel_exposer.def(
                    "Inverse"
                    , Inverse_function_type( &::KDL::FrameVel::Inverse )
                    , ( bp::arg("arg") ) );

        }
        { //::KDL::FrameVel::deriv

            typedef ::KDL::Twist ( ::KDL::FrameVel::*deriv_function_type)(  ) const;

            FrameVel_exposer.def(
                    "deriv"
                    , deriv_function_type( &::KDL::FrameVel::deriv ) );

        }
        FrameVel_exposer.def( bp::self * bp::other< KDL::VectorVel >() );
        FrameVel_exposer.def( bp::self * bp::other< KDL::Vector >() );
        FrameVel_exposer.def( bp::self * bp::other< KDL::TwistVel >() );
        FrameVel_exposer.def( bp::self * bp::other< KDL::Twist >() );
        { //::KDL::FrameVel::operator=

            typedef ::KDL::FrameVel & ( ::KDL::FrameVel::*assign_function_type)( ::KDL::Frame const & ) ;

            FrameVel_exposer.def(
                    "assign"
                    , assign_function_type( &::KDL::FrameVel::operator= )
                    , ( bp::arg("arg") )
                    , bp::return_self< >() );

        }
        { //::KDL::FrameVel::operator=

            typedef ::KDL::FrameVel & ( ::KDL::FrameVel::*assign_function_type)( ::KDL::FrameVel const & ) ;

            FrameVel_exposer.def(
                    "assign"
                    , assign_function_type( &::KDL::FrameVel::operator= )
                    , ( bp::arg("arg") )
                    , bp::return_self< >() );

        }
        { //::KDL::FrameVel::value

            typedef ::KDL::Frame ( ::KDL::FrameVel::*value_function_type)(  ) const;

            FrameVel_exposer.def(
                    "value"
                    , value_function_type( &::KDL::FrameVel::value ) );

        }
        FrameVel_exposer.def_readwrite( "M", &KDL::FrameVel::M );
        FrameVel_exposer.def_readwrite( "p", &KDL::FrameVel::p );
        FrameVel_exposer.staticmethod( "Identity" );
        FrameVel_exposer.def( bp::self * bp::self );
        FrameVel_exposer.def( bp::self * bp::other< KDL::Frame >() );
    }


    { //::KDL::Segment
        typedef bp::class_< KDL::Segment > Segment_exposer_t;
        Segment_exposer_t Segment_exposer = Segment_exposer_t( "Segment", bp::init<KDL::Joint, KDL::Frame >() );
        bp::scope Segment_scope( Segment_exposer );
//        Segment_exposer.def( bp::init< std::string const &, bp::optional< KDL::Joint const &, KDL::Frame const &, KDL::RigidBodyInertia const & > >(( bp::arg("name"), bp::arg("joint")=KDL::Joint(KDL::Joint::None), bp::arg("f_tip")=KDL::Frame::Identity(), bp::arg("I")=KDL::RigidBodyInertia::Zero() )) );
//        Segment_exposer.def( bp::init< bp::optional< KDL::Joint const &, KDL::Frame const &, KDL::RigidBodyInertia const & > >(( bp::arg("joint")=KDL::Joint(KDL::Joint::None), bp::arg("f_tip")=KDL::Frame::Identity(), bp::arg("I")=KDL::RigidBodyInertia::Zero() )) );
        Segment_exposer.def( bp::init< KDL::Segment const & >(( bp::arg("in") )) );
        { //::KDL::Segment::getFrameToTip

            typedef ::KDL::Frame ( ::KDL::Segment::*getFrameToTip_function_type)(  ) const;
            Segment_exposer.def(
                    "getFrameToTip"
                    , getFrameToTip_function_type( &::KDL::Segment::getFrameToTip ) );

        }
        { //::KDL::Segment::getInertia

            typedef ::KDL::RigidBodyInertia const & ( ::KDL::Segment::*getInertia_function_type)(  ) const;

            Segment_exposer.def(
                    "getInertia"
                    , getInertia_function_type( &::KDL::Segment::getInertia )
                    , bp::return_value_policy< bp::copy_const_reference >() );

        }
        { //::KDL::Segment::getJoint

            typedef ::KDL::Joint const & ( ::KDL::Segment::*getJoint_function_type)(  ) const;

            Segment_exposer.def(
                    "getJoint"
                    , getJoint_function_type( &::KDL::Segment::getJoint )
                    , bp::return_value_policy< bp::copy_const_reference >() );

        }
        { //::KDL::Segment::getName

            typedef ::std::string const & ( ::KDL::Segment::*getName_function_type)(  ) const;

            Segment_exposer.def(
                    "getName"
                    , getName_function_type( &::KDL::Segment::getName )
                    , bp::return_value_policy< bp::copy_const_reference >() );

        }
        { //::KDL::Segment::operator=

            typedef ::KDL::Segment & ( ::KDL::Segment::*assign_function_type)( ::KDL::Segment const & ) ;

            Segment_exposer.def(
                    "assign"
                    , assign_function_type( &::KDL::Segment::operator= )
                    , ( bp::arg("arg") )
                    , bp::return_self< >() );

        }
        { //::KDL::Segment::pose

            typedef ::KDL::Frame ( ::KDL::Segment::*pose_function_type)( double const & ) const;

            Segment_exposer.def(
                    "pose"
                    , pose_function_type( &::KDL::Segment::pose )
                    , ( bp::arg("q") ) );

        }
        { //::KDL::Segment::setInertia

            typedef void ( ::KDL::Segment::*setInertia_function_type)( ::KDL::RigidBodyInertia const & ) ;

            Segment_exposer.def(
                    "setInertia"
                    , setInertia_function_type( &::KDL::Segment::setInertia )
                    , ( bp::arg("Iin") ) );

        }
        { //::KDL::Segment::twist

            typedef ::KDL::Twist ( ::KDL::Segment::*twist_function_type)( double const &,double const & ) const;

            Segment_exposer.def(
                    "twist"
                    , twist_function_type( &::KDL::Segment::twist )
                    , ( bp::arg("q"), bp::arg("qdot") ) );

        }
    }



    { //::KDL::Rotation
        typedef bp::class_< Rotation_wrapper > Rotation_exposer_t;
        Rotation_exposer_t Rotation_exposer = Rotation_exposer_t( "Rotation", bp::init< >() );
        bp::scope Rotation_scope( Rotation_exposer );
        Rotation_exposer.def( bp::init< double, double, double, double, double, double, double, double, double >(( bp::arg("Xx"), bp::arg("Yx"), bp::arg("Zx"), bp::arg("Xy"), bp::arg("Yy"), bp::arg("Zy"), bp::arg("Xz"), bp::arg("Yz"), bp::arg("Zz") )) );
        Rotation_exposer.def( bp::init< KDL::Vector const &, KDL::Vector const &, KDL::Vector const & >(( bp::arg("x"), bp::arg("y"), bp::arg("z") )) );
        { //::KDL::Rotation::DoRotX

            typedef void ( ::KDL::Rotation::*DoRotX_function_type)( double ) ;

            Rotation_exposer.def(
                    "DoRotX"
                    , DoRotX_function_type( &::KDL::Rotation::DoRotX )
                    , ( bp::arg("angle") ) );

        }
        { //::KDL::Rotation::DoRotY

            typedef void ( ::KDL::Rotation::*DoRotY_function_type)( double ) ;

            Rotation_exposer.def(
                    "DoRotY"
                    , DoRotY_function_type( &::KDL::Rotation::DoRotY )
                    , ( bp::arg("angle") ) );

        }
        { //::KDL::Rotation::DoRotZ

            typedef void ( ::KDL::Rotation::*DoRotZ_function_type)( double ) ;

            Rotation_exposer.def(
                    "DoRotZ"
                    , DoRotZ_function_type( &::KDL::Rotation::DoRotZ )
                    , ( bp::arg("angle") ) );

        }
        { //::KDL::Rotation::EulerZYX

            typedef ::KDL::Rotation ( *EulerZYX_function_type )( double,double,double );

            Rotation_exposer.def(
                    "EulerZYX"
                    , EulerZYX_function_type( &::KDL::Rotation::EulerZYX )
                    , ( bp::arg("Alfa"), bp::arg("Beta"), bp::arg("Gamma") ) );

        }
        { //::KDL::Rotation::EulerZYZ

            typedef ::KDL::Rotation ( *EulerZYZ_function_type )( double,double,double );

            Rotation_exposer.def(
                    "EulerZYZ"
                    , EulerZYZ_function_type( &::KDL::Rotation::EulerZYZ )
                    , ( bp::arg("Alfa"), bp::arg("Beta"), bp::arg("Gamma") ) );

        }
        { //::KDL::Rotation::GetEulerZYX

            typedef void ( ::KDL::Rotation::*GetEulerZYX_function_type)( double &,double &,double & ) const;

            Rotation_exposer.def(
                    "GetEulerZYX"
                    , GetEulerZYX_function_type( &::KDL::Rotation::GetEulerZYX )
                    , ( bp::arg("Alfa"), bp::arg("Beta"), bp::arg("Gamma") ) );

        }
        { //::KDL::Rotation::GetEulerZYZ

            typedef void ( ::KDL::Rotation::*GetEulerZYZ_function_type)( double &,double &,double & ) const;

            Rotation_exposer.def(
                    "GetEulerZYZ"
                    , GetEulerZYZ_function_type( &::KDL::Rotation::GetEulerZYZ )
                    , ( bp::arg("alpha"), bp::arg("beta"), bp::arg("gamma") ) );

        }
        { //::KDL::Rotation::GetQuaternion

            typedef void ( ::KDL::Rotation::*GetQuaternion_function_type)( double &,double &,double &,double & ) const;

            Rotation_exposer.def(
                    "GetQuaternion"
                    , GetQuaternion_function_type( &::KDL::Rotation::GetQuaternion )
                    , ( bp::arg("x"), bp::arg("y"), bp::arg("z"), bp::arg("w") ) );

        }
        { //::KDL::Rotation::GetRPY

            typedef void ( ::KDL::Rotation::*GetRPY_function_type)( double &,double &,double & ) const;

            Rotation_exposer.def(
                    "GetRPY"
                    , GetRPY_function_type( &::KDL::Rotation::GetRPY )
                    , ( bp::arg("roll"), bp::arg("pitch"), bp::arg("yaw") ) );

        }
        { //::KDL::Rotation::GetRot

            typedef ::KDL::Vector ( ::KDL::Rotation::*GetRot_function_type)(  ) const;

            Rotation_exposer.def(
                    "GetRot"
                    , GetRot_function_type( &::KDL::Rotation::GetRot ) );

        }
        { //::KDL::Rotation::GetRotAngle

            typedef double ( ::KDL::Rotation::*GetRotAngle_function_type)( ::KDL::Vector &,double ) const;

            Rotation_exposer.def(
                    "GetRotAngle"
                    , GetRotAngle_function_type( &::KDL::Rotation::GetRotAngle )
                    , ( bp::arg("axis"), bp::arg("eps")=epsilon ) );

        }
        { //::KDL::Rotation::Identity

            typedef ::KDL::Rotation ( *Identity_function_type )(  );

            Rotation_exposer.def(
                    "Identity"
                    , Identity_function_type( &::KDL::Rotation::Identity ) );

        }
        { //::KDL::Rotation::Inverse

            typedef ::KDL::Rotation ( ::KDL::Rotation::*Inverse_function_type)(  ) const;

            Rotation_exposer.def(
                    "Inverse"
                    , Inverse_function_type( &::KDL::Rotation::Inverse ) );

        }
        { //::KDL::Rotation::Inverse

            typedef ::KDL::Vector ( ::KDL::Rotation::*Inverse_function_type)( ::KDL::Vector const & ) const;

            Rotation_exposer.def(
                    "Inverse"
                    , Inverse_function_type( &::KDL::Rotation::Inverse )
                    , ( bp::arg("v") ) );

        }
        { //::KDL::Rotation::Inverse

            typedef ::KDL::Wrench ( ::KDL::Rotation::*Inverse_function_type)( ::KDL::Wrench const & ) const;

            Rotation_exposer.def(
                    "Inverse"
                    , Inverse_function_type( &::KDL::Rotation::Inverse )
                    , ( bp::arg("arg") ) );

        }
        { //::KDL::Rotation::Inverse

            typedef ::KDL::Twist ( ::KDL::Rotation::*Inverse_function_type)( ::KDL::Twist const & ) const;

            Rotation_exposer.def(
                    "Inverse"
                    , Inverse_function_type( &::KDL::Rotation::Inverse )
                    , ( bp::arg("arg") ) );

        }
        { //::KDL::Rotation::Quaternion

            typedef ::KDL::Rotation ( *Quaternion_function_type )( double,double,double,double );

            Rotation_exposer.def(
                    "Quaternion"
                    , Quaternion_function_type( &::KDL::Rotation::Quaternion )
                    , ( bp::arg("x"), bp::arg("y"), bp::arg("z"), bp::arg("w") ) );

        }
        { //::KDL::Rotation::RPY

            typedef ::KDL::Rotation ( *RPY_function_type )( double,double,double );

            Rotation_exposer.def(
                    "RPY"
                    , RPY_function_type( &::KDL::Rotation::RPY )
                    , ( bp::arg("roll"), bp::arg("pitch"), bp::arg("yaw") ) );

        }
        { //::KDL::Rotation::Rot

            typedef ::KDL::Rotation ( *Rot_function_type )( ::KDL::Vector const &,double );

            Rotation_exposer.def(
                    "Rot"
                    , Rot_function_type( &::KDL::Rotation::Rot )
                    , ( bp::arg("rotvec"), bp::arg("angle") ) );

        }
        { //::KDL::Rotation::Rot2

            typedef ::KDL::Rotation ( *Rot2_function_type )( ::KDL::Vector const &,double );

            Rotation_exposer.def(
                    "Rot2"
                    , Rot2_function_type( &::KDL::Rotation::Rot2 )
                    , ( bp::arg("rotvec"), bp::arg("angle") ) );

        }
        { //::KDL::Rotation::RotX

            typedef ::KDL::Rotation ( *RotX_function_type )( double );

            Rotation_exposer.def(
                    "RotX"
                    , RotX_function_type( &::KDL::Rotation::RotX )
                    , ( bp::arg("angle") ) );

        }
        { //::KDL::Rotation::RotY

            typedef ::KDL::Rotation ( *RotY_function_type )( double );

            Rotation_exposer.def(
                    "RotY"
                    , RotY_function_type( &::KDL::Rotation::RotY )
                    , ( bp::arg("angle") ) );

        }
        { //::KDL::Rotation::RotZ

            typedef ::KDL::Rotation ( *RotZ_function_type )( double );

            Rotation_exposer.def(
                    "RotZ"
                    , RotZ_function_type( &::KDL::Rotation::RotZ )
                    , ( bp::arg("angle") ) );

        }
        { //::KDL::Rotation::SetInverse

            typedef void ( ::KDL::Rotation::*SetInverse_function_type)(  ) ;

            Rotation_exposer.def(
                    "SetInverse"
                    , SetInverse_function_type( &::KDL::Rotation::SetInverse ) );

        }
        { //::KDL::Rotation::UnitX

            typedef ::KDL::Vector ( ::KDL::Rotation::*UnitX_function_type)(  ) const;

            Rotation_exposer.def(
                    "UnitX"
                    , UnitX_function_type( &::KDL::Rotation::UnitX ) );

        }
        { //::KDL::Rotation::UnitX

            typedef void ( ::KDL::Rotation::*UnitX_function_type)( ::KDL::Vector const & ) ;

            Rotation_exposer.def(
                    "UnitX"
                    , UnitX_function_type( &::KDL::Rotation::UnitX )
                    , ( bp::arg("X") ) );

        }
        { //::KDL::Rotation::UnitY

            typedef ::KDL::Vector ( ::KDL::Rotation::*UnitY_function_type)(  ) const;

            Rotation_exposer.def(
                    "UnitY"
                    , UnitY_function_type( &::KDL::Rotation::UnitY ) );

        }
        { //::KDL::Rotation::UnitY

            typedef void ( ::KDL::Rotation::*UnitY_function_type)( ::KDL::Vector const & ) ;

            Rotation_exposer.def(
                    "UnitY"
                    , UnitY_function_type( &::KDL::Rotation::UnitY )
                    , ( bp::arg("X") ) );

        }
        { //::KDL::Rotation::UnitZ

            typedef ::KDL::Vector ( ::KDL::Rotation::*UnitZ_function_type)(  ) const;

            Rotation_exposer.def(
                    "UnitZ"
                    , UnitZ_function_type( &::KDL::Rotation::UnitZ ) );

        }
        { //::KDL::Rotation::UnitZ

            typedef void ( ::KDL::Rotation::*UnitZ_function_type)( ::KDL::Vector const & ) ;

            Rotation_exposer.def(
                    "UnitZ"
                    , UnitZ_function_type( &::KDL::Rotation::UnitZ )
                    , ( bp::arg("X") ) );

        }

        Rotation_exposer.def( bp::self * bp::other< KDL::Vector >() );
        Rotation_exposer.def( bp::self * bp::other< KDL::Twist >() );
        Rotation_exposer.def( bp::self * bp::other< KDL::Wrench >() );
        { //::KDL::Rotation::operator=

            typedef ::KDL::Rotation & ( ::KDL::Rotation::*assign_function_type)( ::KDL::Rotation const & ) ;

            Rotation_exposer.def(
                    "assign"
                    , assign_function_type( &::KDL::Rotation::operator= )
                    , ( bp::arg("arg") )
                    , bp::return_self< >() );

        }
//        { //KDL::Rotation::data [variable], type=double [9]
//
//            typedef pyplusplus::containers::static_sized::array_1_t< double, 9> ( *array_wrapper_creator )( ::KDL::Rotation & );
//
//            Rotation_exposer.add_property( "data"
//                    , bp::make_function( array_wrapper_creator(&Rotation_wrapper::pyplusplus_data_wrapper)
//                            , bp::with_custodian_and_ward_postcall< 0, 1 >() ) );
//        }
        Rotation_exposer.staticmethod( "EulerZYX" );
        Rotation_exposer.staticmethod( "EulerZYZ" );
        Rotation_exposer.staticmethod( "Identity" );
        Rotation_exposer.staticmethod( "Quaternion" );
        Rotation_exposer.staticmethod( "RPY" );
        Rotation_exposer.staticmethod( "Rot" );
        Rotation_exposer.staticmethod( "Rot2" );
        Rotation_exposer.staticmethod( "RotX" );
        Rotation_exposer.staticmethod( "RotY" );
        Rotation_exposer.staticmethod( "RotZ" );
        Rotation_exposer.def( bp::self != bp::self );
        Rotation_exposer.def( bp::self * bp::other< KDL::RigidBodyInertia >() );
        Rotation_exposer.def( bp::self * bp::self );
        Rotation_exposer.def( bp::self * bp::other< KDL::RotationAcc >() );
        Rotation_exposer.def( bp::self * bp::other< KDL::VectorAcc >() );
        Rotation_exposer.def( bp::self * bp::other< KDL::RotationVel >() );
        Rotation_exposer.def( bp::self * bp::other< KDL::VectorVel >() );
        Rotation_exposer.def( bp::self == bp::self );
    }






    { //::KDL::Jacobian
        typedef bp::class_< KDL::Jacobian > Jacobian_exposer_t;
        Jacobian_exposer_t Jacobian_exposer = Jacobian_exposer_t( "Jacobian", bp::init< >() );
        bp::scope Jacobian_scope( Jacobian_exposer );
        Jacobian_exposer.def( bp::init< unsigned int >(( bp::arg("nr_of_columns") )) );
        bp::implicitly_convertible< unsigned int, KDL::Jacobian >();
        Jacobian_exposer.def( bp::init< KDL::Jacobian const & >(( bp::arg("arg") )) );
        { //::KDL::Jacobian::changeBase

            typedef void ( ::KDL::Jacobian::*changeBase_function_type)( ::KDL::Rotation const & ) ;

            Jacobian_exposer.def(
                    "changeBase"
                    , changeBase_function_type( &::KDL::Jacobian::changeBase )
                    , ( bp::arg("rot") ) );

        }
        { //::KDL::Jacobian::changeRefFrame

            typedef void ( ::KDL::Jacobian::*changeRefFrame_function_type)( ::KDL::Frame const & ) ;

            Jacobian_exposer.def(
                    "changeRefFrame"
                    , changeRefFrame_function_type( &::KDL::Jacobian::changeRefFrame )
                    , ( bp::arg("frame") ) );

        }
        { //::KDL::Jacobian::changeRefPoint

            typedef void ( ::KDL::Jacobian::*changeRefPoint_function_type)( ::KDL::Vector const & ) ;

            Jacobian_exposer.def(
                    "changeRefPoint"
                    , changeRefPoint_function_type( &::KDL::Jacobian::changeRefPoint )
                    , ( bp::arg("base_AB") ) );

        }
        { //::KDL::Jacobian::columns

            typedef unsigned int ( ::KDL::Jacobian::*columns_function_type)(  ) const;

            Jacobian_exposer.def(
                    "columns"
                    , columns_function_type( &::KDL::Jacobian::columns ) );

        }
        { //::KDL::Jacobian::getColumn

            typedef ::KDL::Twist ( ::KDL::Jacobian::*getColumn_function_type)( unsigned int ) const;

            Jacobian_exposer.def(
                    "getColumn"
                    , getColumn_function_type( &::KDL::Jacobian::getColumn )
                    , ( bp::arg("i") ) );

        }
        Jacobian_exposer.def( bp::self != bp::self );

        { //::KDL::Jacobian::operator=

            typedef ::KDL::Jacobian & ( ::KDL::Jacobian::*assign_function_type)( ::KDL::Jacobian const & ) ;

            Jacobian_exposer.def(
                    "assign"
                    , assign_function_type( &::KDL::Jacobian::operator= )
                    , ( bp::arg("arg") )
                    , bp::return_self< >() );

        }
        Jacobian_exposer.def( bp::self == bp::self );
        { //::KDL::Jacobian::resize

            typedef void ( ::KDL::Jacobian::*resize_function_type)( unsigned int ) ;

            Jacobian_exposer.def(
                    "resize"
                    , resize_function_type( &::KDL::Jacobian::resize )
                    , ( bp::arg("newNrOfColumns") ) );

        }
        { //::KDL::Jacobian::rows

            typedef unsigned int ( ::KDL::Jacobian::*rows_function_type)(  ) const;

            Jacobian_exposer.def(
                    "rows"
                    , rows_function_type( &::KDL::Jacobian::rows ) );

        }
        { //::KDL::Jacobian::setColumn

            typedef void ( ::KDL::Jacobian::*setColumn_function_type)( unsigned int,::KDL::Twist const & ) ;

            Jacobian_exposer.def(
                    "setColumn"
                    , setColumn_function_type( &::KDL::Jacobian::setColumn )
                    , ( bp::arg("i"), bp::arg("t") ) );

        }
        Jacobian_exposer.def_readwrite( "data", &KDL::Jacobian::data );
    }


    { //::KDL::JntArray
        typedef bp::class_< KDL::JntArray > JntArray_exposer_t;
        JntArray_exposer_t JntArray_exposer = JntArray_exposer_t( "JntArray", bp::init< >() );
        bp::scope JntArray_scope( JntArray_exposer );
        JntArray_exposer.def( bp::init< unsigned int >(( bp::arg("size") )) );
        bp::implicitly_convertible< unsigned int, KDL::JntArray >();
        JntArray_exposer.def( bp::init< KDL::JntArray const & >(( bp::arg("arg") )) );



        JntArray_exposer.def("__getitem__", &JntArray_getitem);
        JntArray_exposer.def("__setitem__", &JntArray_setitem);

        { //::KDL::JntArray::columns

            typedef unsigned int ( ::KDL::JntArray::*columns_function_type)(  ) const;

            JntArray_exposer.def(
                    "columns"
                    , columns_function_type( &::KDL::JntArray::columns ) );

        }
         { //::KDL::JntArray::operator()

             typedef double ( ::KDL::JntArray::*__call___function_type)( unsigned int,unsigned int ) const;

             JntArray_exposer.def(
                 "__call__"
                 , __call___function_type( &::KDL::JntArray::operator() )
                 , ( bp::arg("i"), bp::arg("j")=(unsigned int)(0) ) );

         }
        { //::KDL::JntArray::operator=

            typedef ::KDL::JntArray & ( ::KDL::JntArray::*assign_function_type)( ::KDL::JntArray const & ) ;

            JntArray_exposer.def(
                    "assign"
                    , assign_function_type( &::KDL::JntArray::operator= )
                    , ( bp::arg("arg") )
                    , bp::return_self< >() );

        }
        { //::KDL::JntArray::resize

            typedef void ( ::KDL::JntArray::*resize_function_type)( unsigned int ) ;

            JntArray_exposer.def(
                    "resize"
                    , resize_function_type( &::KDL::JntArray::resize )
                    , ( bp::arg("newSize") ) );

        }
        { //::KDL::JntArray::rows

            typedef unsigned int ( ::KDL::JntArray::*rows_function_type)(  ) const;

            JntArray_exposer.def(
                    "rows"
                    , rows_function_type( &::KDL::JntArray::rows ) );

        }
        JntArray_exposer.def_readwrite( "data", &KDL::JntArray::data );
        JntArray_exposer.def( bp::self == bp::self );
    }

    { //::KDL::JntArrayAcc
        typedef bp::class_< KDL::JntArrayAcc > JntArrayAcc_exposer_t;
        JntArrayAcc_exposer_t JntArrayAcc_exposer = JntArrayAcc_exposer_t( "JntArrayAcc", bp::init< >() );
        bp::scope JntArrayAcc_scope( JntArrayAcc_exposer );
        JntArrayAcc_exposer.def( bp::init< unsigned int >(( bp::arg("size") )) );
        bp::implicitly_convertible< unsigned int, KDL::JntArrayAcc >();
        JntArrayAcc_exposer.def( bp::init< KDL::JntArray const &, KDL::JntArray const &, KDL::JntArray const & >(( bp::arg("q"), bp::arg("qdot"), bp::arg("qdotdot") )) );
        JntArrayAcc_exposer.def( bp::init< KDL::JntArray const &, KDL::JntArray const & >(( bp::arg("q"), bp::arg("qdot") )) );
        JntArrayAcc_exposer.def( bp::init< KDL::JntArray const & >(( bp::arg("q") )) );
        bp::implicitly_convertible< KDL::JntArray const &, KDL::JntArrayAcc >();
        { //::KDL::JntArrayAcc::dderiv

            typedef ::KDL::JntArray ( ::KDL::JntArrayAcc::*dderiv_function_type)(  ) const;

            JntArrayAcc_exposer.def(
                    "dderiv"
                    , dderiv_function_type( &::KDL::JntArrayAcc::dderiv ) );

        }
        { //::KDL::JntArrayAcc::deriv

            typedef ::KDL::JntArray ( ::KDL::JntArrayAcc::*deriv_function_type)(  ) const;

            JntArrayAcc_exposer.def(
                    "deriv"
                    , deriv_function_type( &::KDL::JntArrayAcc::deriv ) );

        }
        { //::KDL::JntArrayAcc::resize

            typedef void ( ::KDL::JntArrayAcc::*resize_function_type)( unsigned int ) ;

            JntArrayAcc_exposer.def(
                    "resize"
                    , resize_function_type( &::KDL::JntArrayAcc::resize )
                    , ( bp::arg("newSize") ) );

        }
        { //::KDL::JntArrayAcc::value

            typedef ::KDL::JntArray ( ::KDL::JntArrayAcc::*value_function_type)(  ) const;

            JntArrayAcc_exposer.def(
                    "value"
                    , value_function_type( &::KDL::JntArrayAcc::value ) );

        }
        JntArrayAcc_exposer.def_readwrite( "q", &KDL::JntArrayAcc::q );
        JntArrayAcc_exposer.def_readwrite( "qdot", &KDL::JntArrayAcc::qdot );
        JntArrayAcc_exposer.def_readwrite( "qdotdot", &KDL::JntArrayAcc::qdotdot );
    }

    { //::KDL::JntArrayVel
        typedef bp::class_< KDL::JntArrayVel > JntArrayVel_exposer_t;
        JntArrayVel_exposer_t JntArrayVel_exposer = JntArrayVel_exposer_t( "JntArrayVel", bp::init< >() );
        bp::scope JntArrayVel_scope( JntArrayVel_exposer );
        JntArrayVel_exposer.def( bp::init< unsigned int >(( bp::arg("size") )) );
        bp::implicitly_convertible< unsigned int, KDL::JntArrayVel >();
        JntArrayVel_exposer.def( bp::init< KDL::JntArray const &, KDL::JntArray const & >(( bp::arg("q"), bp::arg("qdot") )) );
        JntArrayVel_exposer.def( bp::init< KDL::JntArray const & >(( bp::arg("q") )) );
        bp::implicitly_convertible< KDL::JntArray const &, KDL::JntArrayVel >();
        { //::KDL::JntArrayVel::deriv

            typedef ::KDL::JntArray ( ::KDL::JntArrayVel::*deriv_function_type)(  ) const;

            JntArrayVel_exposer.def(
                    "deriv"
                    , deriv_function_type( &::KDL::JntArrayVel::deriv ) );

        }
        { //::KDL::JntArrayVel::resize

            typedef void ( ::KDL::JntArrayVel::*resize_function_type)( unsigned int ) ;

            JntArrayVel_exposer.def(
                    "resize"
                    , resize_function_type( &::KDL::JntArrayVel::resize )
                    , ( bp::arg("newSize") ) );

        }
        { //::KDL::JntArrayVel::value

            typedef ::KDL::JntArray ( ::KDL::JntArrayVel::*value_function_type)(  ) const;

            JntArrayVel_exposer.def(
                    "value"
                    , value_function_type( &::KDL::JntArrayVel::value ) );

        }
        JntArrayVel_exposer.def_readwrite( "q", &KDL::JntArrayVel::q );
        JntArrayVel_exposer.def_readwrite( "qdot", &KDL::JntArrayVel::qdot );
    }


    { //::KDL::Vector
        typedef bp::class_< Vector_wrapper > Vector_exposer_t;
        Vector_exposer_t Vector_exposer = Vector_exposer_t( "Vector", bp::init< >() );
        bp::scope Vector_scope( Vector_exposer );
        Vector_exposer.def( bp::init< double, double, double >(( bp::arg("x"), bp::arg("y"), bp::arg("z") )) );
        Vector_exposer.def( bp::init< KDL::Vector const & >(( bp::arg("arg") )) );

        Vector_exposer.def("__setitem__", &Vector_setitem);
        { //::KDL::Vector::Norm

            typedef double ( ::KDL::Vector::*Norm_function_type)(  ) const;

            Vector_exposer.def(
                    "Norm"
                    , Norm_function_type( &::KDL::Vector::Norm ) );

        }
        { //::KDL::Vector::Normalize

            typedef double ( ::KDL::Vector::*Normalize_function_type)( double ) ;

            Vector_exposer.def(
                    "Normalize"
                    , Normalize_function_type( &::KDL::Vector::Normalize )
                    , ( bp::arg("eps")=epsilon ) );

        }
        { //::KDL::Vector::ReverseSign

            typedef void ( ::KDL::Vector::*ReverseSign_function_type)(  ) ;

            Vector_exposer.def(
                    "ReverseSign"
                    , ReverseSign_function_type( &::KDL::Vector::ReverseSign ) );

        }
        { //::KDL::Vector::Set2DPlane

            typedef void ( ::KDL::Vector::*Set2DPlane_function_type)( ::KDL::Frame const &,::KDL::Vector2 const & ) ;

            Vector_exposer.def(
                    "Set2DPlane"
                    , Set2DPlane_function_type( &::KDL::Vector::Set2DPlane )
                    , ( bp::arg("F_someframe_XY"), bp::arg("v_XY") ) );

        }
        { //::KDL::Vector::Set2DXY

            typedef void ( ::KDL::Vector::*Set2DXY_function_type)( ::KDL::Vector2 const & ) ;

            Vector_exposer.def(
                    "Set2DXY"
                    , Set2DXY_function_type( &::KDL::Vector::Set2DXY )
                    , ( bp::arg("v") ) );

        }
        { //::KDL::Vector::Set2DYZ

            typedef void ( ::KDL::Vector::*Set2DYZ_function_type)( ::KDL::Vector2 const & ) ;

            Vector_exposer.def(
                    "Set2DYZ"
                    , Set2DYZ_function_type( &::KDL::Vector::Set2DYZ )
                    , ( bp::arg("v") ) );

        }
        { //::KDL::Vector::Set2DZX

            typedef void ( ::KDL::Vector::*Set2DZX_function_type)( ::KDL::Vector2 const & ) ;

            Vector_exposer.def(
                    "Set2DZX"
                    , Set2DZX_function_type( &::KDL::Vector::Set2DZX )
                    , ( bp::arg("v") ) );

        }
        { //::KDL::Vector::Zero

            typedef ::KDL::Vector ( *Zero_function_type )(  );

            Vector_exposer.def(
                    "Zero"
                    , Zero_function_type( &::KDL::Vector::Zero ) );

        }
         { //::KDL::Vector::operator()

             typedef double ( ::KDL::Vector::*__call___function_type)( int ) const;

             Vector_exposer.def(
                 "__call__"
                 , __call___function_type( &::KDL::Vector::operator() )
                 , ( bp::arg("index") ) );

         }

        Vector_exposer.def( bp::self += bp::self );
        Vector_exposer.def( bp::self -= bp::self );
        { //::KDL::Vector::operator=

            typedef ::KDL::Vector & ( ::KDL::Vector::*assign_function_type)( ::KDL::Vector const & ) ;

            Vector_exposer.def(
                    "assign"
                    , assign_function_type( &::KDL::Vector::operator= )
                    , ( bp::arg("arg") )
                    , bp::return_self< >() );

        }
        { //::KDL::Vector::operator[]

            typedef double ( ::KDL::Vector::*__getitem___function_type)( int ) const;

            Vector_exposer.def(
                    "__getitem__"
                    , __getitem___function_type( &::KDL::Vector::operator[] )
                    , ( bp::arg("index") ) );

        }
         { //::KDL::Vector::operator[]

             typedef double & ( ::KDL::Vector::*__getitem___function_type)( int ) ;

             Vector_exposer.def(
                 "__getitem__"
                 , __getitem___function_type( &::KDL::Vector::operator[] )
                 , ( bp::arg("index") )
                 , bp::return_value_policy< bp::copy_non_const_reference >() );

         }
        { //::KDL::Vector::x

            typedef double ( ::KDL::Vector::*x_function_type)(  ) const;

            Vector_exposer.def(
                    "x"
                    , x_function_type( &::KDL::Vector::x ) );

        }
        { //::KDL::Vector::x

            typedef void ( ::KDL::Vector::*x_function_type)( double ) ;

            Vector_exposer.def(
                    "x"
                    , x_function_type( &::KDL::Vector::x )
                    , ( bp::arg("arg0") ) );

        }
        { //::KDL::Vector::y

            typedef double ( ::KDL::Vector::*y_function_type)(  ) const;

            Vector_exposer.def(
                    "y"
                    , y_function_type( &::KDL::Vector::y ) );

        }
        { //::KDL::Vector::y

            typedef void ( ::KDL::Vector::*y_function_type)( double ) ;

            Vector_exposer.def(
                    "y"
                    , y_function_type( &::KDL::Vector::y )
                    , ( bp::arg("arg0") ) );

        }
        { //::KDL::Vector::z

            typedef double ( ::KDL::Vector::*z_function_type)(  ) const;

            Vector_exposer.def(
                    "z"
                    , z_function_type( &::KDL::Vector::z ) );

        }
        { //::KDL::Vector::z

            typedef void ( ::KDL::Vector::*z_function_type)( double ) ;

            Vector_exposer.def(
                    "z"
                    , z_function_type( &::KDL::Vector::z )
                    , ( bp::arg("arg0") ) );

        }
//        pyplusplus::containers::static_sized::register_array_1< double, 3 >( "__array_1_double_3" );
//        { //KDL::Vector::data [variable], type=double [3]
//
//            typedef pyplusplus::containers::static_sized::array_1_t< double, 3> ( *array_wrapper_creator )( ::KDL::Vector & );
//
//            Vector_exposer.add_property( "data"
//                    , bp::make_function( array_wrapper_creator(&Vector_wrapper::pyplusplus_data_wrapper)
//                            , bp::with_custodian_and_ward_postcall< 0, 1 >() ) );
//        }
        Vector_exposer.staticmethod( "Zero" );
        Vector_exposer.def( bp::self != bp::self );
        Vector_exposer.def( bp::self * bp::other< double >() );
        Vector_exposer.def( bp::other< double >() * bp::self );
        Vector_exposer.def( bp::self * bp::self );
        Vector_exposer.def( bp::self + bp::self );
        Vector_exposer.def( bp::self - bp::self );
        Vector_exposer.def( -bp::self );
        Vector_exposer.def( bp::self / bp::other< double >() );
        Vector_exposer.def( bp::self == bp::self );
    }



































    { //::KDL::Joint
		typedef bp::class_< KDL::Joint > Joint_exposer_t;
		Joint_exposer_t Joint_exposer = Joint_exposer_t( "Joint");
		bp::scope Joint_scope( Joint_exposer );
		bp::enum_< KDL::Joint::JointType>("JointType")
				.value("RotAxis", KDL::Joint::RotAxis)
				.value("RotX", KDL::Joint::RotX)
				.value("RotY", KDL::Joint::RotY)
				.value("RotZ", KDL::Joint::RotZ)
				.value("TransAxis", KDL::Joint::TransAxis)
				.value("TransX", KDL::Joint::TransX)
				.value("TransY", KDL::Joint::TransY)
				.value("TransZ", KDL::Joint::TransZ)
				.value("NoJoint", KDL::Joint::None)
				.export_values()
				;

        Joint_exposer.def( bp::init< bp::optional< KDL::Joint::JointType const &, double const &, double const &, double const &, double const &, double const & > >(( bp::arg("type")=::KDL::Joint::None, bp::arg("scale")=1, bp::arg("offset")=0, bp::arg("inertia")=0, bp::arg("damping")=0, bp::arg("stiffness")=0 )) );
        Joint_exposer.def( bp::init< std::string const &, KDL::Vector const &, KDL::Vector const &, KDL::Joint::JointType const &, bp::optional< double const &, double const &, double const &, double const &, double const & > >(( bp::arg("name"), bp::arg("_origin"), bp::arg("_axis"), bp::arg("type"), bp::arg("_scale")=1, bp::arg("_offset")=0, bp::arg("_inertia")=0, bp::arg("_damping")=0, bp::arg("_stiffness")=0 )) );
        Joint_exposer.def( bp::init< KDL::Vector const &, KDL::Vector const &, KDL::Joint::JointType const &, bp::optional< double const &, double const &, double const &, double const &, double const & > >(( bp::arg("_origin"), bp::arg("_axis"), bp::arg("type"), bp::arg("_scale")=1, bp::arg("_offset")=0, bp::arg("_inertia")=0, bp::arg("_damping")=0, bp::arg("_stiffness")=0 )) );

        { //::KDL::Joint::JointAxis

            typedef ::KDL::Vector ( ::KDL::Joint::*JointAxis_function_type)(  ) const;

            Joint_exposer.def(
                    "JointAxis"
                    , JointAxis_function_type( &::KDL::Joint::JointAxis ) );

        }
        { //::KDL::Joint::JointOrigin

            typedef ::KDL::Vector ( ::KDL::Joint::*JointOrigin_function_type)(  ) const;

            Joint_exposer.def(
                    "JointOrigin"
                    , JointOrigin_function_type( &::KDL::Joint::JointOrigin ) );

        }
        { //::KDL::Joint::getName

            typedef ::std::string const & ( ::KDL::Joint::*getName_function_type)(  ) const;

            Joint_exposer.def(
                    "getName"
                    , getName_function_type( &::KDL::Joint::getName )
                    , bp::return_value_policy< bp::copy_const_reference >() );

        }
        { //::KDL::Joint::getType

            typedef ::KDL::Joint::JointType const & ( ::KDL::Joint::*getType_function_type)(  ) const;

            Joint_exposer.def(
                    "getType"
                    , getType_function_type( &::KDL::Joint::getType )
                    , bp::return_value_policy< bp::copy_const_reference >() );

        }
        { //::KDL::Joint::getTypeName

            typedef ::std::string const ( ::KDL::Joint::*getTypeName_function_type)(  ) const;

            Joint_exposer.def(
                    "getTypeName"
                    , getTypeName_function_type( &::KDL::Joint::getTypeName ) );

        }
        { //::KDL::Joint::pose

            typedef ::KDL::Frame ( ::KDL::Joint::*pose_function_type)( double const & ) const;

            Joint_exposer.def(
                    "pose"
                    , pose_function_type( &::KDL::Joint::pose )
                    , ( bp::arg("q") ) );

        }
        { //::KDL::Joint::twist

            typedef ::KDL::Twist ( ::KDL::Joint::*twist_function_type)( double const & ) const;

            Joint_exposer.def(
                    "twist"
                    , twist_function_type( &::KDL::Joint::twist )
                    , ( bp::arg("qdot") ) );

        }

    }







    { //::KDL::Chain
        typedef bp::class_< KDL::Chain > Chain_exposer_t;
        Chain_exposer_t Chain_exposer = Chain_exposer_t( "Chain", bp::init< >() );
        bp::scope Chain_scope( Chain_exposer );
        Chain_exposer.def( bp::init< KDL::Chain const & >(( bp::arg("in") )) );
        { //::KDL::Chain::addChain

            typedef void ( ::KDL::Chain::*addChain_function_type)( ::KDL::Chain const & ) ;

            Chain_exposer.def(
                    "addChain"
                    , addChain_function_type( &::KDL::Chain::addChain )
                    , ( bp::arg("chain") ) );

        }
        { //::KDL::Chain::addSegment

            typedef void ( ::KDL::Chain::*addSegment_function_type)( ::KDL::Segment const & ) ;

            Chain_exposer.def(
                    "addSegment"
                    , addSegment_function_type( &::KDL::Chain::addSegment )
                    , ( bp::arg("segment") ) );

        }
        { //::KDL::Chain::getNrOfJoints

            typedef unsigned int ( ::KDL::Chain::*getNrOfJoints_function_type)(  ) const;

            Chain_exposer.def(
                    "getNrOfJoints"
                    , getNrOfJoints_function_type( &::KDL::Chain::getNrOfJoints ) );

        }
        { //::KDL::Chain::getNrOfSegments

            typedef unsigned int ( ::KDL::Chain::*getNrOfSegments_function_type)(  ) const;

            Chain_exposer.def(
                    "getNrOfSegments"
                    , getNrOfSegments_function_type( &::KDL::Chain::getNrOfSegments ) );

        }
        { //::KDL::Chain::getSegment

            typedef ::KDL::Segment const & ( ::KDL::Chain::*getSegment_function_type)( unsigned int ) const;

            Chain_exposer.def(
                    "getSegment"
                    , getSegment_function_type( &::KDL::Chain::getSegment )
                    , ( bp::arg("nr") )
                    , bp::return_value_policy< bp::copy_const_reference >() );

        }
        { //::KDL::Chain::operator=

            typedef ::KDL::Chain & ( ::KDL::Chain::*assign_function_type)( ::KDL::Chain const & ) ;

            Chain_exposer.def(
                    "assign"
                    , assign_function_type( &::KDL::Chain::operator= )
                    , ( bp::arg("arg") )
                    , bp::return_self< >() );

        }
        Chain_exposer.def_readwrite( "segments", &KDL::Chain::segments );
    }


	bp::class_<KDL::SolverI, boost::noncopyable>("SolverI", bp::no_init);

	bp::class_< ChainFkSolverPos_wrapper, bp::bases< KDL::SolverI >, boost::noncopyable >( "ChainFkSolverPos" )
			.def(
					"JntToCart"
					, bp::pure_virtual( (int ( ::KDL::ChainFkSolverPos::* )( ::KDL::JntArray const &,::KDL::Frame &,int ))(&::KDL::ChainFkSolverPos::JntToCart) )
					, ( bp::arg("q_in"), bp::arg("p_out"), bp::arg("segmentNr")=(int)(-1) ) )
			.def(
			"getError"
			, (int ( ::KDL::SolverI::* )(  )const)(&::KDL::SolverI::getError)
			, (int ( ChainFkSolverPos_wrapper::* )(  )const)(&ChainFkSolverPos_wrapper::default_getError) )
			.def(
			"strError"
			, (char const * ( ::KDL::SolverI::* )( int const )const)(&::KDL::SolverI::strError)
			, (char const * ( ChainFkSolverPos_wrapper::* )( int const )const)(&ChainFkSolverPos_wrapper::default_strError)
			, ( bp::arg("error") ) );



    { //::KDL::ChainFkSolverPos_recursive
        typedef bp::class_< ChainFkSolverPos_recursive_wrapper, bp::bases< KDL::ChainFkSolverPos >, boost::noncopyable > ChainFkSolverPos_recursive_exposer_t;
        ChainFkSolverPos_recursive_exposer_t ChainFkSolverPos_recursive_exposer = ChainFkSolverPos_recursive_exposer_t( "ChainFkSolverPos_recursive", bp::init< KDL::Chain const & >(( bp::arg("chain") )) );
        bp::scope ChainFkSolverPos_recursive_scope( ChainFkSolverPos_recursive_exposer );
        bp::implicitly_convertible< KDL::Chain const &, KDL::ChainFkSolverPos_recursive >();
        { //::KDL::ChainFkSolverPos_recursive::JntToCart

            typedef int ( ::KDL::ChainFkSolverPos_recursive::*JntToCart_function_type)( ::KDL::JntArray const &,::KDL::Frame &,int ) ;
            typedef int ( ChainFkSolverPos_recursive_wrapper::*default_JntToCart_function_type)( ::KDL::JntArray const &,::KDL::Frame &,int ) ;

            ChainFkSolverPos_recursive_exposer.def(
                    "JntToCart"
                    , JntToCart_function_type(&::KDL::ChainFkSolverPos_recursive::JntToCart)
                    , default_JntToCart_function_type(&ChainFkSolverPos_recursive_wrapper::default_JntToCart)
                    , ( bp::arg("q_in"), bp::arg("p_out"), bp::arg("segmentNr")=(int)(-1) ) );

        }
        { //::KDL::SolverI::getError

            typedef int ( ::KDL::SolverI::*getError_function_type)(  ) const;
            typedef int ( ChainFkSolverPos_recursive_wrapper::*default_getError_function_type)(  ) const;

            ChainFkSolverPos_recursive_exposer.def(
                    "getError"
                    , getError_function_type(&::KDL::SolverI::getError)
                    , default_getError_function_type(&ChainFkSolverPos_recursive_wrapper::default_getError) );

        }
        { //::KDL::SolverI::strError

            typedef char const * ( ::KDL::SolverI::*strError_function_type)( int const ) const;
            typedef char const * ( ChainFkSolverPos_recursive_wrapper::*default_strError_function_type)( int const ) const;

            ChainFkSolverPos_recursive_exposer.def(
                    "strError"
                    , strError_function_type(&::KDL::SolverI::strError)
                    , default_strError_function_type(&ChainFkSolverPos_recursive_wrapper::default_strError)
                    , ( bp::arg("error") ) );

        }
    }

    bp::class_< ChainFkSolverVel_wrapper, bp::bases< KDL::SolverI >, boost::noncopyable >( "ChainFkSolverVel" )
            .def(
                    "JntToCart"
                    , bp::pure_virtual( (int ( ::KDL::ChainFkSolverVel::* )( ::KDL::JntArrayVel const &,::KDL::FrameVel &,int ))(&::KDL::ChainFkSolverVel::JntToCart) )
    , ( bp::arg("q_in"), bp::arg("out"), bp::arg("segmentNr")=(int)(-1) ) )
    .def(
            "getError"
            , (int ( ::KDL::SolverI::* )(  )const)(&::KDL::SolverI::getError)
            , (int ( ChainFkSolverVel_wrapper::* )(  )const)(&ChainFkSolverVel_wrapper::default_getError) )
    .def(
            "strError"
            , (char const * ( ::KDL::SolverI::* )( int const )const)(&::KDL::SolverI::strError)
            , (char const * ( ChainFkSolverVel_wrapper::* )( int const )const)(&ChainFkSolverVel_wrapper::default_strError)
            , ( bp::arg("error") ) );

    { //::KDL::ChainFkSolverVel_recursive
        typedef bp::class_< ChainFkSolverVel_recursive_wrapper, bp::bases< KDL::ChainFkSolverVel >, boost::noncopyable > ChainFkSolverVel_recursive_exposer_t;
        ChainFkSolverVel_recursive_exposer_t ChainFkSolverVel_recursive_exposer = ChainFkSolverVel_recursive_exposer_t( "ChainFkSolverVel_recursive", bp::init< KDL::Chain const & >(( bp::arg("chain") )) );
        bp::scope ChainFkSolverVel_recursive_scope( ChainFkSolverVel_recursive_exposer );
        bp::implicitly_convertible< KDL::Chain const &, KDL::ChainFkSolverVel_recursive >();
        { //::KDL::ChainFkSolverVel_recursive::JntToCart

            typedef int ( ::KDL::ChainFkSolverVel_recursive::*JntToCart_function_type)( ::KDL::JntArrayVel const &,::KDL::FrameVel &,int ) ;
            typedef int ( ChainFkSolverVel_recursive_wrapper::*default_JntToCart_function_type)( ::KDL::JntArrayVel const &,::KDL::FrameVel &,int ) ;

            ChainFkSolverVel_recursive_exposer.def(
                    "JntToCart"
                    , JntToCart_function_type(&::KDL::ChainFkSolverVel_recursive::JntToCart)
                    , default_JntToCart_function_type(&ChainFkSolverVel_recursive_wrapper::default_JntToCart)
                    , ( bp::arg("q_in"), bp::arg("out"), bp::arg("segmentNr")=(int)(-1) ) );

        }
        { //::KDL::SolverI::getError

            typedef int ( ::KDL::SolverI::*getError_function_type)(  ) const;
            typedef int ( ChainFkSolverVel_recursive_wrapper::*default_getError_function_type)(  ) const;

            ChainFkSolverVel_recursive_exposer.def(
                    "getError"
                    , getError_function_type(&::KDL::SolverI::getError)
                    , default_getError_function_type(&ChainFkSolverVel_recursive_wrapper::default_getError) );

        }
        { //::KDL::SolverI::strError

            typedef char const * ( ::KDL::SolverI::*strError_function_type)( int const ) const;
            typedef char const * ( ChainFkSolverVel_recursive_wrapper::*default_strError_function_type)( int const ) const;

            ChainFkSolverVel_recursive_exposer.def(
                    "strError"
                    , strError_function_type(&::KDL::SolverI::strError)
                    , default_strError_function_type(&ChainFkSolverVel_recursive_wrapper::default_strError)
                    , ( bp::arg("error") ) );

        }
    }



    bp::class_< ChainIdSolver_wrapper, boost::noncopyable >( "ChainIdSolver" )
            .def(
                    "CartToJnt"
                    , bp::pure_virtual( (int ( ::KDL::ChainIdSolver::* )( ::KDL::JntArray const &,::KDL::JntArray const &,::KDL::JntArray const &,::KDL::Wrenches const &,::KDL::JntArray & ))(&::KDL::ChainIdSolver::CartToJnt) )
    , ( bp::arg("q"), bp::arg("q_dot"), bp::arg("q_dotdot"), bp::arg("f_ext"), bp::arg("torques") ) );

    bp::class_< ChainIdSolver_RNE_wrapper, bp::bases< KDL::ChainIdSolver > >( "ChainIdSolver_RNE", bp::init< KDL::Chain const &, KDL::Vector >(( bp::arg("chain"), bp::arg("grav") )) )
            .def(
                    "CartToJnt"
                    , (int ( ::KDL::ChainIdSolver_RNE::* )( ::KDL::JntArray const &,::KDL::JntArray const &,::KDL::JntArray const &,::KDL::Wrenches const &,::KDL::JntArray & ))(&::KDL::ChainIdSolver_RNE::CartToJnt)
            , (int ( ChainIdSolver_RNE_wrapper::* )( ::KDL::JntArray const &,::KDL::JntArray const &,::KDL::JntArray const &,::KDL::Wrenches const &,::KDL::JntArray & ))(&ChainIdSolver_RNE_wrapper::default_CartToJnt)
            , ( bp::arg("q"), bp::arg("q_dot"), bp::arg("q_dotdot"), bp::arg("f_ext"), bp::arg("torques") ) );

    bp::class_< ChainIkSolverPos_wrapper, bp::bases< KDL::SolverI >, boost::noncopyable >( "ChainIkSolverPos" )
            .def(
                    "CartToJnt"
                    , bp::pure_virtual( (int ( ::KDL::ChainIkSolverPos::* )( ::KDL::JntArray const &,::KDL::Frame const &,::KDL::JntArray & ))(&::KDL::ChainIkSolverPos::CartToJnt) )
    , ( bp::arg("q_init"), bp::arg("p_in"), bp::arg("q_out") ) )
    .def(
            "getError"
            , (int ( ::KDL::SolverI::* )(  )const)(&::KDL::SolverI::getError)
            , (int ( ChainIkSolverPos_wrapper::* )(  )const)(&ChainIkSolverPos_wrapper::default_getError) )
    .def(
            "strError"
            , (char const * ( ::KDL::SolverI::* )( int const )const)(&::KDL::SolverI::strError)
            , (char const * ( ChainIkSolverPos_wrapper::* )( int const )const)(&ChainIkSolverPos_wrapper::default_strError)
            , ( bp::arg("error") ) );


    { //::KDL::ChainIkSolverPos_LMA
        typedef bp::class_< ChainIkSolverPos_LMA_wrapper, bp::bases< KDL::ChainIkSolverPos >, boost::noncopyable > ChainIkSolverPos_LMA_exposer_t;
        ChainIkSolverPos_LMA_exposer_t ChainIkSolverPos_LMA_exposer = ChainIkSolverPos_LMA_exposer_t( "ChainIkSolverPos_LMA", bp::init< KDL::Chain const &, Eigen::Matrix< double, 6, 1, 0, 6, 1 > const &, bp::optional< double, int, double > >(( bp::arg("_chain"), bp::arg("_L"), bp::arg("_eps")=1.0000000000000001E-5, bp::arg("_maxiter")=(int)(500), bp::arg("_eps_joints")=1.0000000000000001E-15 )) );
        bp::scope ChainIkSolverPos_LMA_scope( ChainIkSolverPos_LMA_exposer );
        ChainIkSolverPos_LMA_exposer.def( bp::init< KDL::Chain const &, bp::optional< double, int, double > >(( bp::arg("_chain"), bp::arg("_eps")=1.0000000000000001E-5, bp::arg("_maxiter")=(int)(500), bp::arg("_eps_joints")=1.0000000000000001E-15 )) );
        bp::implicitly_convertible< KDL::Chain const &, KDL::ChainIkSolverPos_LMA >();
        { //::KDL::ChainIkSolverPos_LMA::CartToJnt

            typedef int ( ::KDL::ChainIkSolverPos_LMA::*CartToJnt_function_type)( ::KDL::JntArray const &,::KDL::Frame const &,::KDL::JntArray & ) ;
            typedef int ( ChainIkSolverPos_LMA_wrapper::*default_CartToJnt_function_type)( ::KDL::JntArray const &,::KDL::Frame const &,::KDL::JntArray & ) ;

            ChainIkSolverPos_LMA_exposer.def(
                    "CartToJnt"
                    , CartToJnt_function_type(&::KDL::ChainIkSolverPos_LMA::CartToJnt)
                    , default_CartToJnt_function_type(&ChainIkSolverPos_LMA_wrapper::default_CartToJnt)
                    , ( bp::arg("q_init"), bp::arg("T_base_goal"), bp::arg("q_out") ) );

        }
        { //::KDL::ChainIkSolverPos_LMA::compute_fwdpos

         typedef void ( KDL::ChainIkSolverPos_LMA::*compute_fwdpos_function_type)( const VectorXq & ) ;

         ChainIkSolverPos_LMA_exposer.def(
             "compute_fwdpos"
             , compute_fwdpos_function_type( &::KDL::ChainIkSolverPos_LMA::compute_fwdpos )
             , ( bp::arg("q") ) );

        }
        { //::KDL::ChainIkSolverPos_LMA::compute_jacobian

         typedef void ( ::KDL::ChainIkSolverPos_LMA::*compute_jacobian_function_type)(  const VectorXq & ) ;

         ChainIkSolverPos_LMA_exposer.def(
             "compute_jacobian"
             , compute_jacobian_function_type( &::KDL::ChainIkSolverPos_LMA::compute_jacobian )
             , ( bp::arg("q") ) );

        }
        { //::KDL::ChainIkSolverPos_LMA::display_jac

            typedef void ( ::KDL::ChainIkSolverPos_LMA::*display_jac_function_type)( ::KDL::JntArray const & ) ;

            ChainIkSolverPos_LMA_exposer.def(
                    "display_jac"
                    , display_jac_function_type( &::KDL::ChainIkSolverPos_LMA::display_jac )
                    , ( bp::arg("jval") ) );

        }
        ChainIkSolverPos_LMA_exposer.def_readwrite( "T_base_head", &KDL::ChainIkSolverPos_LMA::T_base_head );
        ChainIkSolverPos_LMA_exposer.def_readwrite( "display_information", &KDL::ChainIkSolverPos_LMA::display_information );
        ChainIkSolverPos_LMA_exposer.def_readwrite( "grad", &KDL::ChainIkSolverPos_LMA::grad );
        ChainIkSolverPos_LMA_exposer.def_readwrite( "jac", &KDL::ChainIkSolverPos_LMA::jac );
        ChainIkSolverPos_LMA_exposer.def_readwrite( "lastDifference", &KDL::ChainIkSolverPos_LMA::lastDifference );
        ChainIkSolverPos_LMA_exposer.def_readwrite( "lastNrOfIter", &KDL::ChainIkSolverPos_LMA::lastNrOfIter );
        ChainIkSolverPos_LMA_exposer.def_readwrite( "lastRotDiff", &KDL::ChainIkSolverPos_LMA::lastRotDiff );
        ChainIkSolverPos_LMA_exposer.def_readwrite( "lastSV", &KDL::ChainIkSolverPos_LMA::lastSV );
        ChainIkSolverPos_LMA_exposer.def_readwrite( "lastTransDiff", &KDL::ChainIkSolverPos_LMA::lastTransDiff );
        { //::KDL::SolverI::getError

            typedef int ( ::KDL::SolverI::*getError_function_type)(  ) const;
            typedef int ( ChainIkSolverPos_LMA_wrapper::*default_getError_function_type)(  ) const;

            ChainIkSolverPos_LMA_exposer.def(
                    "getError"
                    , getError_function_type(&::KDL::SolverI::getError)
                    , default_getError_function_type(&ChainIkSolverPos_LMA_wrapper::default_getError) );

        }
        { //::KDL::SolverI::strError

            typedef char const * ( ::KDL::SolverI::*strError_function_type)( int const ) const;
            typedef char const * ( ChainIkSolverPos_LMA_wrapper::*default_strError_function_type)( int const ) const;

            ChainIkSolverPos_LMA_exposer.def(
                    "strError"
                    , strError_function_type(&::KDL::SolverI::strError)
                    , default_strError_function_type(&ChainIkSolverPos_LMA_wrapper::default_strError)
                    , ( bp::arg("error") ) );

        }
    }

    bp::class_< ChainIkSolverPos_NR_wrapper, bp::bases< KDL::ChainIkSolverPos >, boost::noncopyable >( "ChainIkSolverPos_NR", bp::init< KDL::Chain const &, KDL::ChainFkSolverPos &, KDL::ChainIkSolverVel &, bp::optional< unsigned int, double > >(( bp::arg("chain"), bp::arg("fksolver"), bp::arg("iksolver"), bp::arg("maxiter")=(unsigned int)(100), bp::arg("eps")=9.9999999999999995E-7 )) )
            .def(
                    "CartToJnt"
                    , (int ( ::KDL::ChainIkSolverPos_NR::* )( ::KDL::JntArray const &,::KDL::Frame const &,::KDL::JntArray & ))(&::KDL::ChainIkSolverPos_NR::CartToJnt)
            , (int ( ChainIkSolverPos_NR_wrapper::* )( ::KDL::JntArray const &,::KDL::Frame const &,::KDL::JntArray & ))(&ChainIkSolverPos_NR_wrapper::default_CartToJnt)
            , ( bp::arg("q_init"), bp::arg("p_in"), bp::arg("q_out") ) )
    .def(
            "strError"
            , (char const * ( ::KDL::ChainIkSolverPos_NR::* )( int const )const)(&::KDL::ChainIkSolverPos_NR::strError)
            , (char const * ( ChainIkSolverPos_NR_wrapper::* )( int const )const)(&ChainIkSolverPos_NR_wrapper::default_strError)
            , ( bp::arg("error") ) )
    .def_readonly( "E_IKSOLVER_FAILED", -100 )
            .def(
                    "getError"
                    , (int ( ::KDL::SolverI::* )(  )const)(&::KDL::SolverI::getError)
            , (int ( ChainIkSolverPos_NR_wrapper::* )(  )const)(&ChainIkSolverPos_NR_wrapper::default_getError) );

    bp::class_< ChainIkSolverPos_NR_JL_wrapper, bp::bases< KDL::ChainIkSolverPos >, boost::noncopyable >( "ChainIkSolverPos_NR_JL", bp::init< KDL::Chain const &, KDL::JntArray const &, KDL::JntArray const &, KDL::ChainFkSolverPos &, KDL::ChainIkSolverVel &, bp::optional< unsigned int, double > >(( bp::arg("chain"), bp::arg("q_min"), bp::arg("q_max"), bp::arg("fksolver"), bp::arg("iksolver"), bp::arg("maxiter")=(unsigned int)(100), bp::arg("eps")=9.9999999999999995E-7 )) )
            .def(
                    "CartToJnt"
                    , (int ( ::KDL::ChainIkSolverPos_NR_JL::* )( ::KDL::JntArray const &,::KDL::Frame const &,::KDL::JntArray & ))(&::KDL::ChainIkSolverPos_NR_JL::CartToJnt)
            , (int ( ChainIkSolverPos_NR_JL_wrapper::* )( ::KDL::JntArray const &,::KDL::Frame const &,::KDL::JntArray & ))(&ChainIkSolverPos_NR_JL_wrapper::default_CartToJnt)
            , ( bp::arg("q_init"), bp::arg("p_in"), bp::arg("q_out") ) )
    .def(
            "getError"
            , (int ( ::KDL::SolverI::* )(  )const)(&::KDL::SolverI::getError)
            , (int ( ChainIkSolverPos_NR_JL_wrapper::* )(  )const)(&ChainIkSolverPos_NR_JL_wrapper::default_getError) )
    .def(
            "strError"
            , (char const * ( ::KDL::SolverI::* )( int const )const)(&::KDL::SolverI::strError)
            , (char const * ( ChainIkSolverPos_NR_JL_wrapper::* )( int const )const)(&ChainIkSolverPos_NR_JL_wrapper::default_strError)
            , ( bp::arg("error") ) );

    bp::class_< ChainIkSolverVel_wrapper, bp::bases< KDL::SolverI >, boost::noncopyable >( "ChainIkSolverVel" )
            .def(
                    "CartToJnt"
                    , bp::pure_virtual( (int ( ::KDL::ChainIkSolverVel::* )( ::KDL::JntArray const &,::KDL::Twist const &,::KDL::JntArray & ))(&::KDL::ChainIkSolverVel::CartToJnt) )
    , ( bp::arg("q_in"), bp::arg("v_in"), bp::arg("qdot_out") ) )
    .def(
            "CartToJnt"
            , bp::pure_virtual( (int ( ::KDL::ChainIkSolverVel::* )( ::KDL::JntArray const &,::KDL::FrameVel const &,::KDL::JntArrayVel & ))(&::KDL::ChainIkSolverVel::CartToJnt) )
    , ( bp::arg("q_init"), bp::arg("v_in"), bp::arg("q_out") ) )
    .def(
            "getError"
            , (int ( ::KDL::SolverI::* )(  )const)(&::KDL::SolverI::getError)
            , (int ( ChainIkSolverVel_wrapper::* )(  )const)(&ChainIkSolverVel_wrapper::default_getError) )
    .def(
            "strError"
            , (char const * ( ::KDL::SolverI::* )( int const )const)(&::KDL::SolverI::strError)
            , (char const * ( ChainIkSolverVel_wrapper::* )( int const )const)(&ChainIkSolverVel_wrapper::default_strError)
            , ( bp::arg("error") ) );


    { //::KDL::ChainIkSolverVel_pinv
        typedef bp::class_< ChainIkSolverVel_pinv_wrapper, bp::bases< KDL::ChainIkSolverVel >, boost::noncopyable > ChainIkSolverVel_pinv_exposer_t;
        ChainIkSolverVel_pinv_exposer_t ChainIkSolverVel_pinv_exposer = ChainIkSolverVel_pinv_exposer_t( "ChainIkSolverVel_pinv", bp::init< KDL::Chain const &, bp::optional< double, int > >(( bp::arg("chain"), bp::arg("eps")=1.0000000000000001E-5, bp::arg("maxiter")=(int)(150) )) );
        bp::scope ChainIkSolverVel_pinv_scope( ChainIkSolverVel_pinv_exposer );
        bp::implicitly_convertible< KDL::Chain const &, KDL::ChainIkSolverVel_pinv >();
        { //::KDL::ChainIkSolverVel_pinv::CartToJnt

            typedef int ( ::KDL::ChainIkSolverVel_pinv::*CartToJnt_function_type)( ::KDL::JntArray const &,::KDL::Twist const &,::KDL::JntArray & ) ;
            typedef int ( ChainIkSolverVel_pinv_wrapper::*default_CartToJnt_function_type)( ::KDL::JntArray const &,::KDL::Twist const &,::KDL::JntArray & ) ;

            ChainIkSolverVel_pinv_exposer.def(
                    "CartToJnt"
                    , CartToJnt_function_type(&::KDL::ChainIkSolverVel_pinv::CartToJnt)
                    , default_CartToJnt_function_type(&ChainIkSolverVel_pinv_wrapper::default_CartToJnt)
                    , ( bp::arg("q_in"), bp::arg("v_in"), bp::arg("qdot_out") ) );

        }
        { //::KDL::ChainIkSolverVel_pinv::CartToJnt

            typedef int ( ::KDL::ChainIkSolverVel_pinv::*CartToJnt_function_type)( ::KDL::JntArray const &,::KDL::FrameVel const &,::KDL::JntArrayVel & ) ;
            typedef int ( ChainIkSolverVel_pinv_wrapper::*default_CartToJnt_function_type)( ::KDL::JntArray const &,::KDL::FrameVel const &,::KDL::JntArrayVel & ) ;

            ChainIkSolverVel_pinv_exposer.def(
                    "CartToJnt"
                    , CartToJnt_function_type(&::KDL::ChainIkSolverVel_pinv::CartToJnt)
                    , default_CartToJnt_function_type(&ChainIkSolverVel_pinv_wrapper::default_CartToJnt)
                    , ( bp::arg("q_init"), bp::arg("v_in"), bp::arg("q_out") ) );

        }
        { //::KDL::ChainIkSolverVel_pinv::getNrZeroSigmas

            typedef unsigned int ( ::KDL::ChainIkSolverVel_pinv::*getNrZeroSigmas_function_type)(  ) const;

            ChainIkSolverVel_pinv_exposer.def(
                    "getNrZeroSigmas"
                    , getNrZeroSigmas_function_type( &::KDL::ChainIkSolverVel_pinv::getNrZeroSigmas ) );

        }
        { //::KDL::ChainIkSolverVel_pinv::getSVDResult

            typedef int ( ::KDL::ChainIkSolverVel_pinv::*getSVDResult_function_type)(  ) const;

            ChainIkSolverVel_pinv_exposer.def(
                    "getSVDResult"
                    , getSVDResult_function_type( &::KDL::ChainIkSolverVel_pinv::getSVDResult ) );

        }
        { //::KDL::ChainIkSolverVel_pinv::strError

            typedef char const * ( ::KDL::ChainIkSolverVel_pinv::*strError_function_type)( int const ) const;
            typedef char const * ( ChainIkSolverVel_pinv_wrapper::*default_strError_function_type)( int const ) const;

            ChainIkSolverVel_pinv_exposer.def(
                    "strError"
                    , strError_function_type(&::KDL::ChainIkSolverVel_pinv::strError)
                    , default_strError_function_type(&ChainIkSolverVel_pinv_wrapper::default_strError)
                    , ( bp::arg("error") ) );

        }
        ChainIkSolverVel_pinv_exposer.def_readonly( "E_CONVERGE_PINV_SINGULAR", +100 );
        ChainIkSolverVel_pinv_exposer.def_readonly( "E_SVD_FAILED", -100 );
        { //::KDL::SolverI::getError

            typedef int ( ::KDL::SolverI::*getError_function_type)(  ) const;
            typedef int ( ChainIkSolverVel_pinv_wrapper::*default_getError_function_type)(  ) const;

            ChainIkSolverVel_pinv_exposer.def(
                    "getError"
                    , getError_function_type(&::KDL::SolverI::getError)
                    , default_getError_function_type(&ChainIkSolverVel_pinv_wrapper::default_getError) );

        }
    }

    { //::KDL::ChainIkSolverVel_pinv_givens
        typedef bp::class_< ChainIkSolverVel_pinv_givens_wrapper, bp::bases< KDL::ChainIkSolverVel >, boost::noncopyable > ChainIkSolverVel_pinv_givens_exposer_t;
        ChainIkSolverVel_pinv_givens_exposer_t ChainIkSolverVel_pinv_givens_exposer = ChainIkSolverVel_pinv_givens_exposer_t( "ChainIkSolverVel_pinv_givens", bp::init< KDL::Chain const & >(( bp::arg("chain") )) );
        bp::scope ChainIkSolverVel_pinv_givens_scope( ChainIkSolverVel_pinv_givens_exposer );
        bp::implicitly_convertible< KDL::Chain const &, KDL::ChainIkSolverVel_pinv_givens >();
        { //::KDL::ChainIkSolverVel_pinv_givens::CartToJnt

            typedef int ( ::KDL::ChainIkSolverVel_pinv_givens::*CartToJnt_function_type)( ::KDL::JntArray const &,::KDL::Twist const &,::KDL::JntArray & ) ;
            typedef int ( ChainIkSolverVel_pinv_givens_wrapper::*default_CartToJnt_function_type)( ::KDL::JntArray const &,::KDL::Twist const &,::KDL::JntArray & ) ;

            ChainIkSolverVel_pinv_givens_exposer.def(
                    "CartToJnt"
                    , CartToJnt_function_type(&::KDL::ChainIkSolverVel_pinv_givens::CartToJnt)
                    , default_CartToJnt_function_type(&ChainIkSolverVel_pinv_givens_wrapper::default_CartToJnt)
                    , ( bp::arg("q_in"), bp::arg("v_in"), bp::arg("qdot_out") ) );

        }
        { //::KDL::ChainIkSolverVel_pinv_givens::CartToJnt

            typedef int ( ::KDL::ChainIkSolverVel_pinv_givens::*CartToJnt_function_type)( ::KDL::JntArray const &,::KDL::FrameVel const &,::KDL::JntArrayVel & ) ;
            typedef int ( ChainIkSolverVel_pinv_givens_wrapper::*default_CartToJnt_function_type)( ::KDL::JntArray const &,::KDL::FrameVel const &,::KDL::JntArrayVel & ) ;

            ChainIkSolverVel_pinv_givens_exposer.def(
                    "CartToJnt"
                    , CartToJnt_function_type(&::KDL::ChainIkSolverVel_pinv_givens::CartToJnt)
                    , default_CartToJnt_function_type(&ChainIkSolverVel_pinv_givens_wrapper::default_CartToJnt)
                    , ( bp::arg("q_init"), bp::arg("v_in"), bp::arg("q_out") ) );

        }
        { //::KDL::SolverI::getError

            typedef int ( ::KDL::SolverI::*getError_function_type)(  ) const;
            typedef int ( ChainIkSolverVel_pinv_givens_wrapper::*default_getError_function_type)(  ) const;

            ChainIkSolverVel_pinv_givens_exposer.def(
                    "getError"
                    , getError_function_type(&::KDL::SolverI::getError)
                    , default_getError_function_type(&ChainIkSolverVel_pinv_givens_wrapper::default_getError) );

        }
        { //::KDL::SolverI::strError

            typedef char const * ( ::KDL::SolverI::*strError_function_type)( int const ) const;
            typedef char const * ( ChainIkSolverVel_pinv_givens_wrapper::*default_strError_function_type)( int const ) const;

            ChainIkSolverVel_pinv_givens_exposer.def(
                    "strError"
                    , strError_function_type(&::KDL::SolverI::strError)
                    , default_strError_function_type(&ChainIkSolverVel_pinv_givens_wrapper::default_strError)
                    , ( bp::arg("error") ) );

        }
    }

    { //::KDL::ChainIkSolverVel_pinv_nso
        typedef bp::class_< ChainIkSolverVel_pinv_nso_wrapper, bp::bases< KDL::ChainIkSolverVel >, boost::noncopyable > ChainIkSolverVel_pinv_nso_exposer_t;
        ChainIkSolverVel_pinv_nso_exposer_t ChainIkSolverVel_pinv_nso_exposer = ChainIkSolverVel_pinv_nso_exposer_t( "ChainIkSolverVel_pinv_nso", bp::init< KDL::Chain const &, KDL::JntArray, KDL::JntArray, bp::optional< double, int, double > >(( bp::arg("chain"), bp::arg("opt_pos"), bp::arg("weights"), bp::arg("eps")=1.0000000000000001E-5, bp::arg("maxiter")=(int)(150), bp::arg("alpha")=0.25 )) );
        bp::scope ChainIkSolverVel_pinv_nso_scope( ChainIkSolverVel_pinv_nso_exposer );
        ChainIkSolverVel_pinv_nso_exposer.def( bp::init< KDL::Chain const &, bp::optional< double, int, double > >(( bp::arg("chain"), bp::arg("eps")=1.0000000000000001E-5, bp::arg("maxiter")=(int)(150), bp::arg("alpha")=0.25 )) );
        bp::implicitly_convertible< KDL::Chain const &, KDL::ChainIkSolverVel_pinv_nso >();
        { //::KDL::ChainIkSolverVel_pinv_nso::CartToJnt

            typedef int ( ::KDL::ChainIkSolverVel_pinv_nso::*CartToJnt_function_type)( ::KDL::JntArray const &,::KDL::Twist const &,::KDL::JntArray & ) ;
            typedef int ( ChainIkSolverVel_pinv_nso_wrapper::*default_CartToJnt_function_type)( ::KDL::JntArray const &,::KDL::Twist const &,::KDL::JntArray & ) ;

            ChainIkSolverVel_pinv_nso_exposer.def(
                    "CartToJnt"
                    , CartToJnt_function_type(&::KDL::ChainIkSolverVel_pinv_nso::CartToJnt)
                    , default_CartToJnt_function_type(&ChainIkSolverVel_pinv_nso_wrapper::default_CartToJnt)
                    , ( bp::arg("q_in"), bp::arg("v_in"), bp::arg("qdot_out") ) );

        }
        { //::KDL::ChainIkSolverVel_pinv_nso::CartToJnt

            typedef int ( ::KDL::ChainIkSolverVel_pinv_nso::*CartToJnt_function_type)( ::KDL::JntArray const &,::KDL::FrameVel const &,::KDL::JntArrayVel & ) ;
            typedef int ( ChainIkSolverVel_pinv_nso_wrapper::*default_CartToJnt_function_type)( ::KDL::JntArray const &,::KDL::FrameVel const &,::KDL::JntArrayVel & ) ;

            ChainIkSolverVel_pinv_nso_exposer.def(
                    "CartToJnt"
                    , CartToJnt_function_type(&::KDL::ChainIkSolverVel_pinv_nso::CartToJnt)
                    , default_CartToJnt_function_type(&ChainIkSolverVel_pinv_nso_wrapper::default_CartToJnt)
                    , ( bp::arg("q_init"), bp::arg("v_in"), bp::arg("q_out") ) );

        }
        { //::KDL::ChainIkSolverVel_pinv_nso::setAlpha

            typedef int ( ::KDL::ChainIkSolverVel_pinv_nso::*setAlpha_function_type)( double const ) ;
            typedef int ( ChainIkSolverVel_pinv_nso_wrapper::*default_setAlpha_function_type)( double const ) ;

            ChainIkSolverVel_pinv_nso_exposer.def(
                    "setAlpha"
                    , setAlpha_function_type(&::KDL::ChainIkSolverVel_pinv_nso::setAlpha)
                    , default_setAlpha_function_type(&ChainIkSolverVel_pinv_nso_wrapper::default_setAlpha)
                    , ( bp::arg("alpha") ) );

        }
        { //::KDL::ChainIkSolverVel_pinv_nso::setOptPos

            typedef int ( ::KDL::ChainIkSolverVel_pinv_nso::*setOptPos_function_type)( ::KDL::JntArray const & ) ;
            typedef int ( ChainIkSolverVel_pinv_nso_wrapper::*default_setOptPos_function_type)( ::KDL::JntArray const & ) ;

            ChainIkSolverVel_pinv_nso_exposer.def(
                    "setOptPos"
                    , setOptPos_function_type(&::KDL::ChainIkSolverVel_pinv_nso::setOptPos)
                    , default_setOptPos_function_type(&ChainIkSolverVel_pinv_nso_wrapper::default_setOptPos)
                    , ( bp::arg("opt_pos") ) );

        }
        { //::KDL::ChainIkSolverVel_pinv_nso::setWeights

            typedef int ( ::KDL::ChainIkSolverVel_pinv_nso::*setWeights_function_type)( ::KDL::JntArray const & ) ;
            typedef int ( ChainIkSolverVel_pinv_nso_wrapper::*default_setWeights_function_type)( ::KDL::JntArray const & ) ;

            ChainIkSolverVel_pinv_nso_exposer.def(
                    "setWeights"
                    , setWeights_function_type(&::KDL::ChainIkSolverVel_pinv_nso::setWeights)
                    , default_setWeights_function_type(&ChainIkSolverVel_pinv_nso_wrapper::default_setWeights)
                    , ( bp::arg("weights") ) );

        }
        { //::KDL::SolverI::getError

            typedef int ( ::KDL::SolverI::*getError_function_type)(  ) const;
            typedef int ( ChainIkSolverVel_pinv_nso_wrapper::*default_getError_function_type)(  ) const;

            ChainIkSolverVel_pinv_nso_exposer.def(
                    "getError"
                    , getError_function_type(&::KDL::SolverI::getError)
                    , default_getError_function_type(&ChainIkSolverVel_pinv_nso_wrapper::default_getError) );

        }
        { //::KDL::SolverI::strError

            typedef char const * ( ::KDL::SolverI::*strError_function_type)( int const ) const;
            typedef char const * ( ChainIkSolverVel_pinv_nso_wrapper::*default_strError_function_type)( int const ) const;

            ChainIkSolverVel_pinv_nso_exposer.def(
                    "strError"
                    , strError_function_type(&::KDL::SolverI::strError)
                    , default_strError_function_type(&ChainIkSolverVel_pinv_nso_wrapper::default_strError)
                    , ( bp::arg("error") ) );

        }
    }

    { //::KDL::ChainIkSolverVel_wdls
        typedef bp::class_< ChainIkSolverVel_wdls_wrapper, bp::bases< KDL::ChainIkSolverVel >, boost::noncopyable > ChainIkSolverVel_wdls_exposer_t;
        ChainIkSolverVel_wdls_exposer_t ChainIkSolverVel_wdls_exposer = ChainIkSolverVel_wdls_exposer_t( "ChainIkSolverVel_wdls", bp::init< KDL::Chain const &, bp::optional< double, int > >(( bp::arg("chain"), bp::arg("eps")=1.0000000000000001E-5, bp::arg("maxiter")=(int)(150) )) );
        bp::scope ChainIkSolverVel_wdls_scope( ChainIkSolverVel_wdls_exposer );
        bp::implicitly_convertible< KDL::Chain const &, KDL::ChainIkSolverVel_wdls >();
        { //::KDL::ChainIkSolverVel_wdls::CartToJnt

            typedef int ( ::KDL::ChainIkSolverVel_wdls::*CartToJnt_function_type)( ::KDL::JntArray const &,::KDL::Twist const &,::KDL::JntArray & ) ;
            typedef int ( ChainIkSolverVel_wdls_wrapper::*default_CartToJnt_function_type)( ::KDL::JntArray const &,::KDL::Twist const &,::KDL::JntArray & ) ;

            ChainIkSolverVel_wdls_exposer.def(
                    "CartToJnt"
                    , CartToJnt_function_type(&::KDL::ChainIkSolverVel_wdls::CartToJnt)
                    , default_CartToJnt_function_type(&ChainIkSolverVel_wdls_wrapper::default_CartToJnt)
                    , ( bp::arg("q_in"), bp::arg("v_in"), bp::arg("qdot_out") ) );

        }
        { //::KDL::ChainIkSolverVel_wdls::CartToJnt

            typedef int ( ::KDL::ChainIkSolverVel_wdls::*CartToJnt_function_type)( ::KDL::JntArray const &,::KDL::FrameVel const &,::KDL::JntArrayVel & ) ;
            typedef int ( ChainIkSolverVel_wdls_wrapper::*default_CartToJnt_function_type)( ::KDL::JntArray const &,::KDL::FrameVel const &,::KDL::JntArrayVel & ) ;

            ChainIkSolverVel_wdls_exposer.def(
                    "CartToJnt"
                    , CartToJnt_function_type(&::KDL::ChainIkSolverVel_wdls::CartToJnt)
                    , default_CartToJnt_function_type(&ChainIkSolverVel_wdls_wrapper::default_CartToJnt)
                    , ( bp::arg("q_init"), bp::arg("v_in"), bp::arg("q_out") ) );

        }
        { //::KDL::ChainIkSolverVel_wdls::getLambda

            typedef double ( ::KDL::ChainIkSolverVel_wdls::*getLambda_function_type)(  ) const;

            ChainIkSolverVel_wdls_exposer.def(
                    "getLambda"
                    , getLambda_function_type( &::KDL::ChainIkSolverVel_wdls::getLambda ) );

        }
        { //::KDL::ChainIkSolverVel_wdls::getLambdaScaled

            typedef double ( ::KDL::ChainIkSolverVel_wdls::*getLambdaScaled_function_type)(  ) const;

            ChainIkSolverVel_wdls_exposer.def(
                    "getLambdaScaled"
                    , getLambdaScaled_function_type( &::KDL::ChainIkSolverVel_wdls::getLambdaScaled ) );

        }
        { //::KDL::ChainIkSolverVel_wdls::getNrZeroSigmas

            typedef unsigned int ( ::KDL::ChainIkSolverVel_wdls::*getNrZeroSigmas_function_type)(  ) const;

            ChainIkSolverVel_wdls_exposer.def(
                    "getNrZeroSigmas"
                    , getNrZeroSigmas_function_type( &::KDL::ChainIkSolverVel_wdls::getNrZeroSigmas ) );

        }
        { //::KDL::ChainIkSolverVel_wdls::getSVDResult

            typedef int ( ::KDL::ChainIkSolverVel_wdls::*getSVDResult_function_type)(  ) const;

            ChainIkSolverVel_wdls_exposer.def(
                    "getSVDResult"
                    , getSVDResult_function_type( &::KDL::ChainIkSolverVel_wdls::getSVDResult ) );

        }
        { //::KDL::ChainIkSolverVel_wdls::getSigmaMin

            typedef double ( ::KDL::ChainIkSolverVel_wdls::*getSigmaMin_function_type)(  ) const;

            ChainIkSolverVel_wdls_exposer.def(
                    "getSigmaMin"
                    , getSigmaMin_function_type( &::KDL::ChainIkSolverVel_wdls::getSigmaMin ) );

        }
        { //::KDL::ChainIkSolverVel_wdls::setEps

            typedef void ( ::KDL::ChainIkSolverVel_wdls::*setEps_function_type)( double const ) ;

            ChainIkSolverVel_wdls_exposer.def(
                    "setEps"
                    , setEps_function_type( &::KDL::ChainIkSolverVel_wdls::setEps )
                    , ( bp::arg("eps_in") ) );

        }
        { //::KDL::ChainIkSolverVel_wdls::setLambda

            typedef void ( ::KDL::ChainIkSolverVel_wdls::*setLambda_function_type)( double const ) ;

            ChainIkSolverVel_wdls_exposer.def(
                    "setLambda"
                    , setLambda_function_type( &::KDL::ChainIkSolverVel_wdls::setLambda )
                    , ( bp::arg("lambda") ) );

        }
        { //::KDL::ChainIkSolverVel_wdls::setMaxIter

            typedef void ( ::KDL::ChainIkSolverVel_wdls::*setMaxIter_function_type)( int const ) ;

            ChainIkSolverVel_wdls_exposer.def(
                    "setMaxIter"
                    , setMaxIter_function_type( &::KDL::ChainIkSolverVel_wdls::setMaxIter )
                    , ( bp::arg("maxiter_in") ) );

        }
        { //::KDL::ChainIkSolverVel_wdls::setWeightJS

            typedef void ( ::KDL::ChainIkSolverVel_wdls::*setWeightJS_function_type)( ::Eigen::MatrixXd const & ) ;

            ChainIkSolverVel_wdls_exposer.def(
                    "setWeightJS"
                    , setWeightJS_function_type( &::KDL::ChainIkSolverVel_wdls::setWeightJS )
                    , ( bp::arg("Mq") ) );

        }
        { //::KDL::ChainIkSolverVel_wdls::setWeightTS

            typedef void ( ::KDL::ChainIkSolverVel_wdls::*setWeightTS_function_type)( ::Eigen::MatrixXd const & ) ;

            ChainIkSolverVel_wdls_exposer.def(
                    "setWeightTS"
                    , setWeightTS_function_type( &::KDL::ChainIkSolverVel_wdls::setWeightTS )
                    , ( bp::arg("Mx") ) );

        }
        { //::KDL::ChainIkSolverVel_wdls::strError

            typedef char const * ( ::KDL::ChainIkSolverVel_wdls::*strError_function_type)( int const ) const;
            typedef char const * ( ChainIkSolverVel_wdls_wrapper::*default_strError_function_type)( int const ) const;

            ChainIkSolverVel_wdls_exposer.def(
                    "strError"
                    , strError_function_type(&::KDL::ChainIkSolverVel_wdls::strError)
                    , default_strError_function_type(&ChainIkSolverVel_wdls_wrapper::default_strError)
                    , ( bp::arg("error") ) );

        }
        ChainIkSolverVel_wdls_exposer.def_readonly( "E_CONVERGE_PINV_SINGULAR", +100 );
        ChainIkSolverVel_wdls_exposer.def_readonly( "E_SVD_FAILED", -100 );
        { //::KDL::SolverI::getError

            typedef int ( ::KDL::SolverI::*getError_function_type)(  ) const;
            typedef int ( ChainIkSolverVel_wdls_wrapper::*default_getError_function_type)(  ) const;

            ChainIkSolverVel_wdls_exposer.def(
                    "getError"
                    , getError_function_type(&::KDL::SolverI::getError)
                    , default_getError_function_type(&ChainIkSolverVel_wdls_wrapper::default_getError) );

        }
    }


    { //::KDL::ChainJntToJacSolver
        typedef bp::class_< ChainJntToJacSolver_wrapper, bp::bases< KDL::SolverI >, boost::noncopyable > ChainJntToJacSolver_exposer_t;
        ChainJntToJacSolver_exposer_t ChainJntToJacSolver_exposer = ChainJntToJacSolver_exposer_t( "ChainJntToJacSolver", bp::init< KDL::Chain const & >(( bp::arg("chain") )) );
        bp::scope ChainJntToJacSolver_scope( ChainJntToJacSolver_exposer );
        bp::implicitly_convertible< KDL::Chain const &, KDL::ChainJntToJacSolver >();
        { //::KDL::ChainJntToJacSolver::JntToJac

            typedef int ( ::KDL::ChainJntToJacSolver::*JntToJac_function_type)( ::KDL::JntArray const &,::KDL::Jacobian &,int ) ;
            typedef int ( ChainJntToJacSolver_wrapper::*default_JntToJac_function_type)( ::KDL::JntArray const &,::KDL::Jacobian &,int ) ;

            ChainJntToJacSolver_exposer.def(
                    "JntToJac"
                    , JntToJac_function_type(&::KDL::ChainJntToJacSolver::JntToJac)
                    , default_JntToJac_function_type(&ChainJntToJacSolver_wrapper::default_JntToJac)
                    , ( bp::arg("q_in"), bp::arg("jac"), bp::arg("segmentNR")=(int)(-1) ) );

        }
        { //::KDL::ChainJntToJacSolver::setLockedJoints

            typedef int ( ::KDL::ChainJntToJacSolver::*setLockedJoints_function_type)( ::std::vector< bool > const ) ;

            ChainJntToJacSolver_exposer.def(
                    "setLockedJoints"
                    , setLockedJoints_function_type( &::KDL::ChainJntToJacSolver::setLockedJoints )
                    , ( bp::arg("locked_joints") ) );

        }
        { //::KDL::ChainJntToJacSolver::strError

            typedef char const * ( ::KDL::ChainJntToJacSolver::*strError_function_type)( int const ) const;
            typedef char const * ( ChainJntToJacSolver_wrapper::*default_strError_function_type)( int const ) const;

            ChainJntToJacSolver_exposer.def(
                    "strError"
                    , strError_function_type(&::KDL::ChainJntToJacSolver::strError)
                    , default_strError_function_type(&ChainJntToJacSolver_wrapper::default_strError)
                    , ( bp::arg("error") ) );

        }
        ChainJntToJacSolver_exposer.def_readonly( "E_JAC_FAILED", -100 );
        { //::KDL::SolverI::getError

            typedef int ( ::KDL::SolverI::*getError_function_type)(  ) const;
            typedef int ( ChainJntToJacSolver_wrapper::*default_getError_function_type)(  ) const;

            ChainJntToJacSolver_exposer.def(
                    "getError"
                    , getError_function_type(&::KDL::SolverI::getError)
                    , default_getError_function_type(&ChainJntToJacSolver_wrapper::default_getError) );

        }
    }






    bp::def("JntToCart",JntToCart);

    bp::def("CartToJnt",CartToJnt);
}


int main()
{
	PyImport_AppendInittab("PyKDL", &PyInit_PyKDL); // Add example to built-in.
	Py_Initialize(); // Start interpreter.
	try
	{
        bp::object main = bp::import("__main__");
        bp::object main_namespace = main.attr("__dict__");
        bp::scope scope(main); // Force main scope
    	main_namespace["PyKDL"] = bp::import("PyKDL");
	}
	catch (const bp::error_already_set&)
  	{
		PyErr_Print();
  	}

}