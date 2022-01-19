#pragma once

#include <iostream>
#include <ctime>
#include <cmath>
#include <vector>
#include "Eigen/Dense"

using namespace std;
using namespace Eigen;

//this class defines a line segment in 2D Cartesian system
namespace DFN
{
	class Line
	{
	public:
		double xa; // x coordinate
		double ya;
		double xb;
		double yb;
		Line(std::vector<Vector3d> A);
		Line(Vector3d A, Vector3d B);

		double get_max_x(); //find max x coordinate value
		double get_min_x();
		double get_max_y();
		double get_min_y();
	};

	inline Line::Line(std::vector<Vector3d> A)
	{
		this->xa = A[0](0);
		this->ya = A[0](1);
		this->xb = A[1](0);
		this->yb = A[1](1);
	};

	inline Line::Line(Vector3d A, Vector3d B)
	{
		this->xa = A(0);
		this->ya = A(1);
		this->xb = B(0);
		this->yb = B(1);
	}

	inline double Line::get_max_x()
	{
		return xa > xb ? xa : xb;
	}

	inline double Line::get_min_x()
	{
		return xa > xb ? xb : xa;
	}

	inline double Line::get_max_y()
	{
		return ya > yb ? ya : yb;
	}

	inline double Line::get_min_y()
	{
		return ya > yb ? yb : ya;
	}
}; // namespace DFN
