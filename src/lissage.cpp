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
typedef Polyhedron::Facet_handle Facet;
typedef Polyhedron::Halfedge_around_facet_circulator Halfedge_facet_circulator;
typedef Polyhedron::Halfedge_around_vertex_circulator Halfedge_vertex_circulator;
typedef std::map<Polyhedron::Vertex_handle, Point_3> Vertex_coord_map;
typedef std::map<Polyhedron::Vertex_handle, double> Vertex_double_map;

double maxDist = 0;

void getCoords(Polyhedron &mesh, Vertex_coord_map &coords) {
	Vertex_iterator v_it = mesh.vertices_begin();
	while (v_it != mesh.vertices_end()) {
		coords[v_it] = v_it->point();
		++v_it;
	}
}

void areaComputing(Polyhedron &mesh, Vertex_double_map &map) {
	Vertex_iterator v_it = mesh.vertices_begin();
	double area;
	Point_3 p,q,r;
	while (v_it != mesh.vertices_end()) {
		area = 0.0;
		Halfedge_vertex_circulator h_it = v_it->vertex_begin();
		do {
			if (!(h_it->is_border())) {
				Facet f = h_it->facet();
				Halfedge_facet_circulator hf_it = f->facet_begin();
				p = hf_it->vertex()->point();
				++hf_it;
				q = hf_it->vertex()->point();
				++hf_it;
				r = hf_it->vertex()->point();
				area += sqrt(CGAL::squared_area(p,q,r));
			}
			++h_it;
		} while (h_it != v_it->vertex_begin());
		map[v_it] = area;
		++v_it;
	}
}

void laplacianSmoothing(Polyhedron &mesh, Vertex_coord_map &coords, Vertex_coord_map &newCoords) {
	Vertex_iterator v_it = mesh.vertices_begin();
	Vertex_double_map mapArea;
	double x, y, z;
	double sumArea;
	areaComputing(mesh, mapArea);
	while (v_it != mesh.vertices_end()) {
		x = 0;
		y = 0;
		z = 0;
		sumArea = 0.0;
		Halfedge_vertex_circulator h_it = v_it->vertex_begin();
		do {
			x += coords[h_it->opposite()->vertex()].x() * mapArea[v_it];
			y += coords[h_it->opposite()->vertex()].y() * mapArea[v_it];
			z += coords[h_it->opposite()->vertex()].z() * mapArea[v_it];
			sumArea += mapArea[v_it];
			++h_it;
		} while (h_it != v_it->vertex_begin());
		newCoords[v_it] = Point_3(x/sumArea, y/sumArea, z/sumArea);
		++v_it;
	}
}

void save(Polyhedron &mesh, Vertex_coord_map &coords, const char *filename) {
	std::fstream output(filename);
	output << "OFF" << std::endl; //ligne d'entête
	output << mesh.size_of_vertices() << " " << mesh.size_of_facets() << " 0" << std::endl; //infos sur le mesh
	for (Vertex_iterator v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it) {
		output << coords[v_it] << std::endl;
	}
	for (Facet_iterator it = mesh.facets_begin(); it != mesh.facets_end(); ++it) {
		Halfedge_facet_circulator j = it->facet_begin();
        CGAL_assertion(CGAL::circulator_size(j) >= 3);
        output << CGAL::circulator_size(j) << ' ';
        do {
            output << ' ' << std::distance(mesh.vertices_begin(), j->vertex());
        } while ( ++j != it->facet_begin());
        output << std::endl;
	}
}

void distanceBetweenCoords(Polyhedron &mesh, Vertex_coord_map &coords1, Vertex_coord_map &coords2, Vertex_double_map &distMap) {
	Vertex_iterator v_it = mesh.vertices_begin();
	double distance;
	while (v_it != mesh.vertices_end()) {
		distance = sqrt(CGAL::squared_distance(coords1[v_it], coords2[v_it]));
		if (distance > maxDist) maxDist = distance;
		distMap[v_it] = distance;
		++v_it;
	}
}

void saveColorDistance(Polyhedron &mesh, Vertex_coord_map &coords1, Vertex_coord_map &coords2, const char *filename) {
	std::fstream output(filename);
	output << "COFF" << std::endl; //ligne d'entête
	output << mesh.size_of_vertices() << " " << mesh.size_of_facets() << " 0" << std::endl; //infos sur le mesh
	Vertex_double_map distMap;
	distanceBetweenCoords(mesh, coords1, coords2, distMap);
	for (Vertex_iterator v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it) {
		output << coords2[v_it] << " 255 " << (1.0 - (distMap[v_it]/maxDist)) * 255 << ' ' << (1.0 - (distMap[v_it]/maxDist)) * 255 << " 1" << std::endl;
	}
	for (Facet_iterator it = mesh.facets_begin(); it != mesh.facets_end(); ++it) {
		Halfedge_facet_circulator j = it->facet_begin();
        CGAL_assertion(CGAL::circulator_size(j) >= 3);
        output << CGAL::circulator_size(j) << ' ';
        do {
            output << ' ' << std::distance(mesh.vertices_begin(), j->vertex());
        } while ( ++j != it->facet_begin());
        output << std::endl;
	}
}

int main(int argc, char *argv[]) {
	if (argc < 2) {
		std::cerr << "Il manque un paramètre au programme. Veuillez lui donner en entrée un nom de fichier au format off." << std::endl;
		return 1;
	}
	
	Polyhedron mesh;
	std::ifstream input(argv[1]);
	if (!input || !(input >> mesh) || mesh.is_empty()) {
		std::cerr << "Le fichier donné n'est pas un fichier off valide." << std::endl;
		return 1;
	}
	
	Vertex_coord_map coords0, coords1;
	getCoords(mesh, coords0);
	laplacianSmoothing(mesh, coords0, coords1);
	saveColorDistance(mesh, coords0, coords1, "resultMesh.off");
	return 0;
}
