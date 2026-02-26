/**
 ******************************************************************************
 * \file tardigrade_finite_element_utilities.cpp
 ******************************************************************************
 * The source file for utilities which can assist using the balance equations
 * in finite element codes. We here assume that unknown quantities are
 * computed at the evaluation point using interpolation functions (here called
 * interp) and are projected to the nodes using test functions (here called
 * test).
 ******************************************************************************
 */

#include "tardigrade_finite_element_utilities.h"
