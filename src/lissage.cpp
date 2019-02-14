#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <map>
#include <vector>
#include <queue>
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
typedef CGAL::Vector_3<Kernel> Vector3;
typedef Kernel::Point_3 Point_3;
typedef Polyhedron::Facet_iterator Facet_iterator;
typedef Polyhedron::Vertex_iterator Vertex_iterator;
typedef Polyhedron::Halfedge_iterator Halfedge_iterator;
typedef Polyhedron::Halfedge_around_facet_circulator Halfedge_facet_circulator;
typedef std::map<Polyhedron::Facet_handle, double> Facet_double_map;
typedef std::map<Polyhedron::Facet_handle, int> Facet_int_map;
typedef std::map<Polyhedron::Facet_handle, bool> Facet_bool_map;
typedef std::vector<double> Vect_double;
typedef std::vector<int> Vect_int;
typedef std::vector<Vect_double> Vect_Color;
typedef std::queue<Polyhedron::Facet_handle> Facet_Queue;

int main(int argc, char *argv[]) {
	return 0;
}
