//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   test_StateERME.cc
 *  @author Jeremy Roberts
 *  @date   Aug 19, 2012
 *  @brief  Test of StateERME class.
 *  @note   Copyright (C) 2012 Jeremy Roberts.
 */
//---------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST           \
        FUNC(test_StateERME)

#include "linear_algebra/LinearAlgebraSetup.hh"
#include "utilities/TestDriver.hh"
#include "StateERME.hh"
#include <iostream>

using namespace erme;
using namespace detran_test;
using detran_utilities::soft_equiv;
using std::cout;
using std::endl;

int main(int argc, char *argv[])
{
  RUN(argc, argv);
}

//----------------------------------------------//
// TEST DEFINITIONS
//----------------------------------------------//

// Test of basic public interface
int test_StateERME(int argc, char *argv[])
{
	serment_comm::Comm::initialize(argc, argv);
	int ng = 1;
	if (serment_comm::Comm::size() > 2) ng = 2;
	serment_comm::Comm::setup_communicators(ng);
	linear_algebra::initialize(argc, argv);

	{

		// Create node list
		StateERME state(10);

		linear_algebra::Vector::SP_vector v;
		if (serment_comm::Comm::is_global())
		{
			v = new linear_algebra::Vector(10, 1.0);
		}

		state.set_k(1.12);
		state.set_lambda(1.13);
		TEST(soft_equiv(state.k(),      1.12));
		TEST(soft_equiv(state.lambda(), 1.13));

		state.update(v, 2.0, 3.0);

		TEST(soft_equiv(state.k(),      2.0));
		TEST(soft_equiv(state.lambda(), 3.0));

		std::cout << " local size = " << state.local_size() << std::endl;

	}

  linear_algebra::finalize();
  serment_comm::Comm::finalize();

  return 0;

}

//---------------------------------------------------------------------------//
//              end of test_StateERME.cc
//---------------------------------------------------------------------------//
