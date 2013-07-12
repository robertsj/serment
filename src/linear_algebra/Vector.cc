//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file   Vector.cc
 *  @brief  Vector member definitions
 *  @note   Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#include "Vector.hh"
#include "comm/Comm.hh"

namespace linear_algebra
{

//----------------------------------------------------------------------------//
Vector::Vector(const size_type m, const double val, const bool seq)
  : d_local_size(m)
  , d_is_assembled(false)
  , d_is_temporary(false)
  , d_is_sequential(seq)
{
  Require(m > 0);

  if (serment_comm::Comm::size() == 1) d_is_sequential = true;

  // Create and initialize the vector
  PetscErrorCode ierr;
  if (d_is_sequential)
    ierr = VecCreateSeq(PETSC_COMM_SELF, m, &d_V);
  else
    ierr = VecCreateMPI(PETSC_COMM_WORLD, m, PETSC_DETERMINE, &d_V);
  ierr = VecSet(d_V, val);

  // Set sizes and ranges
  ierr = setup();

  // Because we initialize to zero, we can assemble right away
  assemble();

  Ensure(!ierr);
  Ensure(d_upper_bound - d_lower_bound == d_local_size);
}

//----------------------------------------------------------------------------//
Vector::Vector(const Vector& V)
{
  Require(V.is_assembled());
  Require(!V.is_temporary());

  d_is_sequential = V.is_sequential();
  d_local_size    = V.local_size();
  d_global_size   = V.global_size();

  // Create the vector
  PetscErrorCode ierr;
  if (d_is_sequential)
    ierr = VecCreateSeq(PETSC_COMM_SELF, d_local_size, &d_V);
  else
    ierr = VecCreateMPI(PETSC_COMM_WORLD, d_local_size, PETSC_DETERMINE, &d_V);

  // Set sizes and ranges
  ierr = setup();

  // Because we initialize to zero, we can assemble right
  // away, since we might just be filling this from a multiplication.
  assemble();

  // Copy the incoming vector
  copy(V);

  Ensure(!ierr);
  Ensure(d_upper_bound - d_lower_bound == d_local_size);
}

//----------------------------------------------------------------------------//
Vector::Vector(Vector &V, const int m)
  : d_is_assembled(true)
  , d_is_temporary(false)
{
  Insist(m <= V.local_size(), "Can't wrap a smaller vector.");

  // While this *is* a temporary object, the Vec created does not
  // allocate memory (as we pass it that of the incoming V), and
  // hence, it must be destroyed.

  d_is_sequential = V.is_sequential();
  PetscErrorCode ierr = 0;
  if (serment_comm::Comm::size() == 1 || d_is_sequential)
  {
    ierr =
      VecCreateSeqWithArray(PETSC_COMM_WORLD, 1, m, &V[0], &d_V);
  }
  else
  {
    ierr =
      VecCreateMPIWithArray(PETSC_COMM_WORLD, 1, m, PETSC_DECIDE, &V[0], &d_V);
  }

  d_local_size = m;

  ierr = setup();

  Ensure(!ierr);
  Ensure(d_upper_bound - d_lower_bound == d_local_size);
}

//----------------------------------------------------------------------------//
Vector::Vector(Vec pv)
  : d_is_assembled(true)
  , d_is_temporary(true)
  , d_is_sequential(false)
{
  d_V = pv;
  int gs = 0, ls = 0;
  PetscErrorCode ierr = 0;
  ierr = VecGetSize(pv, &gs);
  ierr = VecGetLocalSize(pv, &ls);
  // if local size = global size, the vector can't live on >1 process
  if (gs == ls) d_is_sequential = true;
  d_local_size = ls;
  ierr = setup();
  Ensure(!ierr);
}

//----------------------------------------------------------------------------//
Vector::~Vector()
{
  if (!d_is_temporary) VecDestroy(&d_V);
}

//----------------------------------------------------------------------------//
void Vector::assemble()
{
  PetscErrorCode ierr = 0;
  if (1)//!d_is_assembled)
  {
    ierr = VecAssemblyBegin(d_V);
    ierr = VecAssemblyEnd(d_V);
    d_is_assembled = true;
  }
  Ensure(!ierr);
}

//----------------------------------------------------------------------------//
void Vector::display(const size_type    output,
                     const std::string &name) const
{
  Require(output < END_VECTOR_DISPLAY_TYPES);
  PetscErrorCode ierr = 0;
  if (output == STDOUT)
  {
    ierr = is_sequential() ? VecView(d_V, PETSC_VIEWER_STDOUT_SELF)
                           : VecView(d_V, PETSC_VIEWER_STDOUT_WORLD);
  }
  else
  {
    PetscViewer viewer;
    MPI_Comm com = is_sequential() ? PETSC_COMM_SELF : PETSC_COMM_WORLD;
    if (output == BINARY)
    {
      ierr = PetscViewerBinaryOpen(com, name.c_str(), FILE_MODE_WRITE, &viewer);
    }
    else if (output == ASCII)
    {
      ierr = PetscViewerASCIIOpen(com, name.c_str(), &viewer);
    }
    ierr = VecView(d_V, viewer);
    PetscViewerDestroy(&viewer);
  }
  Ensure(!ierr);
}

//----------------------------------------------------------------------------//
Vector::SP_vector Vector::collect_on_root(const size_type root)
{
  using serment_comm::Comm;
  int ierr = 0;
  SP_vector V_seq;

  if (Comm::rank() == root)
  {
    V_seq = new Vector(d_global_size, 0.0, true);
    for (int i = d_lower_bound; i < d_upper_bound; ++i)
      (*V_seq)[i]= (*this)[i-d_lower_bound];
  }

  // get all sizes
  vec_int sizes = ranges();

  // send chunks
  for (int p = 0; p < Comm::size(); ++p)
  {
    if (p == root) continue;
    if (Comm::rank() == p)
    {
      // send my chunk to root
      ierr = Comm::send(&d_local_size,  1,            root); Assert(!ierr);
      ierr = Comm::send(&d_array[0],    d_local_size, root); Assert(!ierr);
    }
    else if (Comm::rank() == root)
    {
      // receive chunk
      size_type chunk_size;
      Comm::receive(&chunk_size,             1,            p);
      Comm::receive(&(*V_seq)[0+sizes[p]],   chunk_size,   p);
    }
  }

  if (Comm::rank() == root)
    V_seq->assemble();

  Ensure(!ierr);
  return V_seq;
}

//----------------------------------------------------------------------------//
Vector::vec_int Vector::ranges() const
{
  vec_int v(serment_comm::Comm::size()+1, 0);
  const int *v_a;
  PetscErrorCode ierr = VecGetOwnershipRanges(d_V, &v_a);
  for (int i = 0; i < v.size(); ++i) v[i] = v_a[i];
  Ensure(!ierr);
  return v;
}

//----------------------------------------------------------------------------//
PetscErrorCode Vector::setup()
{
  // Get the global size.
  int gs;
  PetscErrorCode ierr = VecGetSize(d_V, &gs);
  Assert(gs > 0);
  d_global_size = gs;

  // Get the global ranges.
  int lb, ub;
  ierr = VecGetOwnershipRange(d_V, &lb, &ub);
  Assert(lb >= 0 and ub > 0);
  d_lower_bound = lb;
  d_upper_bound = ub;

  // Get access to and then "restore" the array, but not actually. PETSc
  // requires the restore call to ensure safe coding, but by passing null, we
  // get to keep it.  We'll code safely...
  ierr = VecGetArray(d_V, &d_array);
  ierr = VecRestoreArray(d_V, PETSC_NULL);

  return ierr;
}

} // end namespace linear_algebra

//----------------------------------------------------------------------------//
//              end of file Vector.cc
//----------------------------------------------------------------------------//
