// ================================================================================================
// 
// This file is part of M2S2 - Matrices for Mechanices of Solids and Structures
//
// Copyright(C) 2023 
//		Dorival Piedade Neto &
//		Rodrigo Ribeiro Paccola &
//		Rog�rio Carrazedo
// 
// This source code form is subject to the terms of the Apache License 2.0.
// If a copy of Apache License 2.0 was not distributed with this file, you can obtain one at
// https://www.apache.org/licenses/LICENSE-2.0
// 
// ================================================================================================
#pragma once

// Standard libraries
#include <array>
#include <vector>
#include <iostream>		// required by std::cout
#include <iomanip>		// Required by ios manipulations
#include <sstream>		// required by std::ostringstream
#include <cmath>		// required by std::sqrt / std::acos
#include <utility>      // required by std::move
#include <algorithm>    // required by std::copy() and std::assign()
#include <cassert>		// required by assert (programing checks)
#include <exception>	// required by exception (runtime errors)

// Coloring output
#define RESET "\x1B[0m"
#define RED "\x1B[31m"
#define YELLOW "\x1B[33m"
#define GREEN "\x1B[32m"
#define MAGENTA "\x1B[35m"
#define CYAN "\x1B[36m"

#define ERROR(x) RED "ERROR: " x "\n" RESET
#define WARN(x) YELLOW "WARNING: " x "\n" RESET
#define INFO(x) CYAN x RESET
#define INPUT(x) MAGENTA x RESET
#define OK(x) GREEN x RESET

constexpr double MV_PI = 3.14159265358979323846;
constexpr double mg_zero = 0.;
