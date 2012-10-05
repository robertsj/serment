//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   xerme.cc
 *  @brief  xerme
 *  @author Jeremy Roberts
 *  @date   Oct 5, 2012
 */
//---------------------------------------------------------------------------//

#include "utilities/Timer.hh"
#include <iostream>
#include <ctime>

void print_welcome();

int main(int argc, char **argv)
{

  //START_PROFILER();

  // Print the welcome header.
  print_welcome();

  // Timer
  detran_utilities::Timer timer;

  // Start timer.
  timer.tic();

  // Stop timer
  timer.toc(true);

  //STOP_PROFILER();

  return 0;
}

void print_welcome()
{
  std::time_t t;
  std::time(&t);
  std::cout << "---SERMENT---" << std::endl;
  std::cout << " Run on: " << std::ctime(&t);
  std::cout << std::endl << std::endl;
}
//---------------------------------------------------------------------------//
//              end of file xerme.cc
//---------------------------------------------------------------------------//
