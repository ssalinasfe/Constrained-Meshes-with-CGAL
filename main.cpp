#define BOOST_BIND_NO_PLACEHOLDERS

//Memory usage
//#include <malloc_count-0.7.1/malloc_count.h>
#include <malloc_count-0.7.1/malloc_count.c>

// CGAL
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/IO/Color.h>
#include <CGAL/Triangulation_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Constrained_triangulation_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
//#include <CGAL/draw_triangulation_2.h>
#include <CGAL/intersections.h>
#include <CGAL/Regular_triangulation_2.h>
#include <fstream>
#include <iterator>
#include <array>
#include <vector>
#include <iostream>
#include <cmath>
#include <string>
#include <sstream>
#include <random>       // std::default_random_engine
#include <algorithm>    // std::shuffle
#include <cstdlib> 
#include <cctype>
#include <iomanip>
#include <set>

struct ElementInfo
{
	int index;
	bool border = false;
};

//Used for all
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

//Used for the random triangulation
typedef CGAL::Triangulation_vertex_base_with_info_2<ElementInfo, K>    Vb;
typedef CGAL::Triangulation_face_base_with_info_2<ElementInfo,K> Fbb;//era fb
typedef CGAL::Constrained_triangulation_face_base_2<K,Fbb> Fb;  //added
typedef CGAL::Triangulation_data_structure_2<Vb,Fb> Tds;
typedef CGAL::Exact_predicates_tag Itag; // added

//Typedef for the constrained triangulation
typedef CGAL::Constrained_triangulation_2<K, Tds, Itag>  ConstrainedTriangulation;
typedef ConstrainedTriangulation::Constrained_edges_iterator Constrained_edges_iterator;
typedef ConstrainedTriangulation::Face_handle Face_handle;
typedef ConstrainedTriangulation::Edge Edge;
typedef ConstrainedTriangulation::Point Point;
typedef CGAL::Polygon_2<K> Polygon_2;

//typedefs for the constrained Delaunay Triangulation
typedef CGAL::Constrained_Delaunay_triangulation_2<K, Tds, Itag>  ConstrainedDelaunay;
typedef ConstrainedDelaunay::Constrained_edges_iterator Constrained_edges_iterator;
typedef ConstrainedDelaunay::Face_handle Face_handle;
typedef ConstrainedDelaunay::Edge Edge;
typedef ConstrainedDelaunay::Point Point;

//Regular triangulation
typedef CGAL::Regular_triangulation_vertex_base_2<K>                   Vbase;
typedef CGAL::Triangulation_vertex_base_with_info_2<ElementInfo, K,Vbase> VbR;
typedef CGAL::Regular_triangulation_face_base_2<K>                     FbR;
typedef CGAL::Triangulation_data_structure_2<VbR,FbR>                    TdsR;
typedef CGAL::Regular_triangulation_2<K, TdsR>                          Regular;
typedef K::Point_2                                                     Point;
typedef K::Weighted_point_2                                            Wpoint;


//typedefs for Cropped Voronoi
typedef CGAL::Delaunay_triangulation_2<K>  Delaunay;
typedef Delaunay::Edge_iterator  Edge_iterator;
typedef K::Segment_2 Segment_2;
typedef K::Point_2 Point_2;
typedef K::Iso_rectangle_2 Iso_rectangle_2;
typedef K::Segment_2 Segment_2;
typedef K::Ray_2 Ray_2;
typedef K::Line_2 Line_2;

#include <CGAL/Voronoi_diagram_2.h>
#include <CGAL/Delaunay_triangulation_adaptation_traits_2.h>
#include <CGAL/Delaunay_triangulation_adaptation_policies_2.h>
//#include <CGAL/draw_voronoi_diagram_2.h>

// typedefs for defining the adaptor
typedef CGAL::Delaunay_triangulation_2<K>                                    DT;
typedef CGAL::Delaunay_triangulation_adaptation_traits_2<DT>                 AT;
typedef CGAL::Delaunay_triangulation_caching_degeneracy_removal_policy_2<DT> AP;
typedef CGAL::Voronoi_diagram_2<DT,AT,AP>                                    VD;
// typedef for the result type of the point location
typedef AT::Site_2                    Site_2;


//Read node file in .node format and nodes in point vector
std::vector<Point> read_nodes_from_file(std::string name){
		std::string line;
		std::ifstream nodefile(name);
		double a1, a2, a3, a4;
		std::vector<Point> points;
		unsigned n_vertices = 0;
		std::string tmp;
		if (nodefile.is_open())
		{
			while (std::getline(nodefile, line)){ //add check boundary vertices flag
				std::istringstream(line) >> tmp;
				if (tmp[0] != '#' ) //check if first line is a comentary
				{
						n_vertices = std::stoi(tmp);
						points.reserve(n_vertices);
						break;
				}
			}
			while (std::getline(nodefile, line))
			{
				std::istringstream(line) >> tmp;
				if (tmp[0] != '#' ) //check if first line is a comentary
				{
					std::istringstream(line) >> a1 >> a2 >> a3 >> a4;
					points.push_back(Point(a2, a3));
				}
			}  
		}
		else 
				std::cout << "Unable to open node file"; 
		nodefile.close();
		return points;
}


//Print off file based in a Triangulation_2, the file name and the number of interior faces
void print_off(ConstrainedTriangulation Tr, std::string name, int n_faces){
		std::cout << "Printing off file" << std::endl;
		std::ofstream outfile(name + ".1.off");
		outfile << "OFF" << std::endl;
		outfile << Tr.number_of_vertices() << " " << n_faces << " " << 0 << std::endl;
		for(ConstrainedTriangulation::Vertex_handle v : Tr.finite_vertex_handles()){
				outfile << std::setprecision(15)<< v->point().x() <<" "<<v->point().y()<<" 0.0\n";
		}
		for(ConstrainedTriangulation::Face_handle f : Tr.finite_face_handles()){
			if(f->info().index!=-1)
				outfile << 3 << " " << f->vertex(0)->info().index << " " << f->vertex(1)->info().index << " " << f->vertex(2)->info().index << std::endl;
		}
		outfile.close();
		std::cout << "Off file printed as "<< name << ".1.off" << std::endl;
}

//Print a node file based in a Triangulation_2, the file name and the number of interior faces
void print_node_file(ConstrainedTriangulation Tr, std::string name){
		std::cout<<"Printing node file"<<std::endl;
		std::ofstream outfile(name + ".1.node");
		outfile << Tr.number_of_vertices() <<" 2 0 1"<< std::endl;
		for(ConstrainedTriangulation::Vertex_handle v : Tr.finite_vertex_handles()){
				outfile << v->info().index<<" "<< std::setprecision(15)<< v->point().x() << " " << v->point().y() << " "<<v->info().border << std::endl;
		}
		outfile.close();
		std::cout<<"Node file printed as "<<name + ".1.node"<<std::endl;
}

//Print triangle file based in a Triangulation_2, the file name and the number of interior faces
void print_ele_file(ConstrainedTriangulation Tr, std::string name, int n_faces){
		std::cout<<"Printing ele file"<<std::endl;
		std::ofstream outfile(name + ".1.ele");
		outfile << n_faces <<" 3 0"<< std::endl;
		for(ConstrainedTriangulation::Face_handle f : Tr.finite_face_handles()){
				if(f->info().index!=-1)
						outfile << f->info().index << " " << f->vertex(0)->info().index << " " << f->vertex(1)->info().index << " " << f->vertex(2)->info().index << std::endl;
		}
		outfile.close();
		std::cout<<"Ele file printed as "<<name + ".1.ele"<<std::endl;
}

//Print a adjacency file based in a Triangulation_2, the file name and the number of interior faces
void print_neigh_file(ConstrainedTriangulation Tr, std::string name, int n_faces){
	printf("Printing neigh file\n");
	std::ofstream outfile(name + ".1.neigh");
	outfile << n_faces <<" 3"<< std::endl;
	for(ConstrainedTriangulation::Face_handle f : Tr.finite_face_handles()){
			if(f->info().index!=-1){
					outfile << f->info().index << " " << f->neighbor(0)->info().index << " " << f->neighbor(1)->info().index << " " << f->neighbor(2)->info().index << std::endl;
			}
	}
	outfile.close();
	printf("Neigh file printed as %s.1.neigh\n", name.c_str());
}


bool check_inside(Point pt, std::vector<Point> boundary, K traits)
{
	switch(CGAL::bounded_side_2(boundary.begin(), boundary.end(), pt, traits)) {
		case CGAL::ON_BOUNDED_SIDE :
			return true;
		case CGAL::ON_BOUNDARY:
			return true;
		case CGAL::ON_UNBOUNDED_SIDE:
			return false;
	}
	return false;
}


//Print a constrained Triangulation TR with a constrained boundaryu as a node file, ele file, neigh file and off file for vizualization
void printConstrainedTriangulation(ConstrainedTriangulation Tr, std::vector<Point> boundary, std::string output_file){
	int index = 0;
	for(ConstrainedTriangulation::Vertex_handle v : Tr.finite_vertex_handles()) {
		v->info().index = index;
		index++;
	}

	for (Constrained_edges_iterator it = Tr.constrained_edges_begin(); it != Tr.constrained_edges_end(); ++it) {
		Edge e = *it;
		e.first->vertex( (e.second+1)%3 )->info().border=true;
		e.first->vertex( (e.second+2)%3 )->info().border=true;
	}

	index = 0;
	// Calcula mid point of each edge of the face, to check if the face is inside the boundary
	for(ConstrainedTriangulation::Face_handle f : Tr.all_face_handles()){
		if(!Tr.is_infinite(f)){
			auto v0 = f->vertex(0)->point();
			auto v1 = f->vertex(1)->point();
			auto v2 = f->vertex(2)->point();
			Point mid1 = Point((v0.x()+v1.x())/2,(v0.y()+v1.y())/2);
			Point mid2 = Point((v1.x()+v2.x())/2,(v1.y()+v2.y())/2);
			Point mid3 = Point((v2.x()+v0.x())/2,(v2.y()+v0.y())/2);
			if(check_inside(mid1, boundary, K()) && check_inside(mid2,boundary, K()) && check_inside(mid3, boundary, K())){
				f->info().index = index;
				index++;
			}
			else{
				f->info().index = -1;
			}
		}else{
			f->info().index = -1;
		}
	}

	//std::cout<<"Printing mesh with "<<index<<" faces and "<< Tr.number_of_vertices() <<" vertices"<<std::endl;
	//print_off(Tr, output_file, index);
	print_node_file(Tr, output_file);
	print_ele_file(Tr, output_file, index);
	print_neigh_file(Tr, output_file, index);
	//std::cout<<"Printed files in "<<output_file<<std::endl;
}

//Print a print triangulation with the .node, ele and neigh. files. It uses ConstrainedDelunay, but it is not a constrained triangulation
void printDelaunayTriangulation(ConstrainedDelaunay Tr, std::string output_file){
	int index = 0;
	for(ConstrainedDelaunay::Vertex_handle v : Tr.finite_vertex_handles()) {
		v->info().index = index;
		index++;
	}

	for (Constrained_edges_iterator it = Tr.constrained_edges_begin(); it != Tr.constrained_edges_end(); ++it) {
		Edge e = *it;
		e.first->vertex( (e.second+1)%3 )->info().border=true;
		e.first->vertex( (e.second+2)%3 )->info().border=true;
	}

	index = 0;
	for(ConstrainedDelaunay::Face_handle f : Tr.all_face_handles()){
		if(!Tr.is_infinite(f)){
			f->info().index = index;
			index++;
		}else{
			f->info().index = -1;
		}
	}

	//std::cout<<"Printing mesh with "<<index<<" faces and "<< Tr.number_of_vertices() <<" vertices"<<std::endl;
	//print_off(Tr, output_file, index);
	print_node_file(Tr, output_file);
	print_ele_file(Tr, output_file, index);
	print_neigh_file(Tr, output_file, index);
	//std::cout<<"Printed files in "<<output_file<<std::endl;
}


//Print a constrained Triangulation TR with a constrained boundaryu as a node file, ele file, neigh file and off file for vizualization
void printConstrainedTriangulation(ConstrainedDelaunay Tr, std::vector<Point> boundary, std::string output_file){
	int index = 0;
	for(ConstrainedDelaunay::Vertex_handle v : Tr.finite_vertex_handles()) {
		v->info().index = index;
		index++;
	}

	for (Constrained_edges_iterator it = Tr.constrained_edges_begin(); it != Tr.constrained_edges_end(); ++it) {
		Edge e = *it;
		e.first->vertex( (e.second+1)%3 )->info().border=true;
		e.first->vertex( (e.second+2)%3 )->info().border=true;
	}

	index = 0;
	for(ConstrainedDelaunay::Face_handle f : Tr.all_face_handles()){
		if(!Tr.is_infinite(f)){
			auto v0 = f->vertex(0)->point();
			auto v1 = f->vertex(1)->point();
			auto v2 = f->vertex(2)->point();
			Point mid1 = Point((v0.x()+v1.x())/2,(v0.y()+v1.y())/2);
			Point mid2 = Point((v1.x()+v2.x())/2,(v1.y()+v2.y())/2);
			Point mid3 = Point((v2.x()+v0.x())/2,(v2.y()+v0.y())/2);
			if(check_inside(mid1, boundary, K()) && check_inside(mid2,boundary, K()) && check_inside(mid3, boundary, K())){
				f->info().index = index;
				index++;
			}
			else{
				f->info().index = -1;
			}
		}else{
			f->info().index = -1;
		}
	}

	//std::cout<<"Printing mesh with "<<index<<" faces and "<< Tr.number_of_vertices() <<" vertices"<<std::endl;
	print_off(Tr, output_file, index);
	print_node_file(Tr, output_file);
	print_ele_file(Tr, output_file, index);
	print_neigh_file(Tr, output_file, index);
	//std::cout<<"Printed files in "<<output_file<<std::endl;
}


void printRegularOFF(Regular Rtr, std::vector<Point> non_hidden_points, std::string output_file){
	std::ofstream outfile(output_file + "_regular.off");
	outfile<<"OFF"<<std::endl;
	outfile<<Rtr.number_of_vertices()<<" "<<Rtr.number_of_faces()<<" "<<0<<std::endl;
	
  	for (auto v : non_hidden_points){
		outfile<<v.x()<<" "<<v.y()<<" "<<0<<std::endl;
	}

	for (auto it = Rtr.finite_faces_begin(); it != Rtr.finite_faces_end(); it++)
	{
		int v1 = std::distance(non_hidden_points.begin(), std::find(non_hidden_points.begin(), non_hidden_points.end(), it->vertex(0)->point()) );
		int v2 = std::distance(non_hidden_points.begin(), std::find(non_hidden_points.begin(), non_hidden_points.end(), it->vertex(1)->point()) );
		int v3 = std::distance(non_hidden_points.begin(), std::find(non_hidden_points.begin(), non_hidden_points.end(), it->vertex(2)->point()) );
		outfile<<"3 "<<v1<<" "<<v2<<" "<<v3<<std::endl;
	}
	outfile.close();
}

void printDelaunayOFF(Delaunay Tr, std::vector<Point> points, std::string output_file){
	std::ofstream outfile(output_file + "_delaunay.off");
	outfile<<"OFF"<<std::endl;
	outfile<<Tr.number_of_vertices()<<" "<<Tr.number_of_faces()<<" "<<0<<std::endl;

  	for (auto v : points){
		outfile<<v.x()<<" "<<v.y()<<" "<<0<<std::endl;
	}

	for (auto it = Tr.finite_faces_begin(); it != Tr.finite_faces_end(); it++)
	{
		int v1 = std::distance(points.begin(), std::find(points.begin(), points.end(), it->vertex(0)->point()) );
		int v2 = std::distance(points.begin(), std::find(points.begin(), points.end(), it->vertex(1)->point()) );
		int v3 = std::distance(points.begin(), std::find(points.begin(), points.end(), it->vertex(2)->point()) );
		outfile<<"3 "<<v1<<" "<<v2<<" "<<v3<<std::endl;
	}
	outfile.close();
}

void printRandomOFF(ConstrainedTriangulation Tr, std::vector<Point> points, std::string output_file){
	std::ofstream outfile(output_file + "_randomTR.off");
	outfile<<"OFF"<<std::endl;
	outfile<<Tr.number_of_vertices()<<" "<<Tr.number_of_faces()<<" "<<0<<std::endl;

  	for (auto v : points){
		outfile<<v.x()<<" "<<v.y()<<" "<<0<<std::endl;
	}

	for (auto it = Tr.finite_faces_begin(); it != Tr.finite_faces_end(); it++)
	{
		int v1 = std::distance(points.begin(), std::find(points.begin(), points.end(), it->vertex(0)->point()) );
		int v2 = std::distance(points.begin(), std::find(points.begin(), points.end(), it->vertex(1)->point()) );
		int v3 = std::distance(points.begin(), std::find(points.begin(), points.end(), it->vertex(2)->point()) );
		outfile<<"3 "<<v1<<" "<<v2<<" "<<v3<<std::endl;
	}
	outfile.close();
}

void printVoronoiOFF(std::vector<std::vector<Point>> voronoi_mesh, std::string output_file){
	std::set<Point> unique_points;

	//Generate a set with unique points in nlogn time
	for(auto region: voronoi_mesh){
		for(auto p : region){
			unique_points.insert(p);
		}
	}
	std::ofstream outfile(output_file + "_voronoi.off");
	outfile<<"OFF"<<std::endl;
	outfile<<unique_points.size()<<" "<<voronoi_mesh.size()<<" "<<0<<std::endl;

  	for (auto v : unique_points){
		outfile<<v.x()<<" "<<v.y()<<" "<<0<<std::endl;
	}

	//int index = 0;
	for(auto region : voronoi_mesh){
		outfile<<region.size()<<" ";
		for(auto point : region){
			int v1 = std::distance(unique_points.begin(), std::find(unique_points.begin(), unique_points.end(), point) );
			outfile<<v1<<" ";
		}
		//std::cout<<"printing region "<<index++<<std::endl;
		outfile<<std::endl;
	}

	outfile.close();


}

// Cut a Voronoi edge with a given boundary
//The boundary is assumed to be polygon in ccw, 
//boundary is a vector with the segments of the boundary and the boundary_points is a vector with the points of the boundary
//if an edge only containts its source point inside the boundary, the edge is cut and we store the source point, the intersection and the taget of the boundary point
//if an edge only containts its target point inside the boundary, the edge is cut and we only store its target point
//Function return a vector with the points of the cell in ccw
//if the boundary containts collinear points, the should be removed after generate the cell
std::vector<Point> cut_halfedge(VD &vd, VD::Halfedge_handle e, std::vector<Segment_2> &boundary, std::vector<Point> &boundary_points){
	//std::cout << "\t";
	std::vector<Point> constrained_region;
	bool flag = false;
	//if e is a segment
	if(e->is_segment()){
		bool flag = false;
		Segment_2 seg(e->source()->point(), e->target()->point());
		//Check intersection with the boundary
		for(auto borderseg : boundary){
			const auto res = CGAL::intersection(seg, borderseg);
			if (res)
			{
				const auto pt = boost::get<Point>(&*res);
				if(CGAL::area(borderseg.source(), borderseg.target(), seg.source()) > 0){ //if source is inside the border
					constrained_region.push_back(seg.source());
					constrained_region.push_back(*pt);
					constrained_region.push_back(borderseg.target()); //add the border point
					//std::cout<<seg.source()<<" "<<*pt<<" ";
				}else{ //target is inside the border
					//std::cout<<"c"<<*pt<<" ";
					constrained_region.push_back(*pt);
				}
				flag = true;
				break;
			}
		}
		if(!flag){ //if the segment does not intersect with the border, check if is inside the border
			if(check_inside(seg.source(), boundary_points, K()) && check_inside(seg.target(), boundary_points, K()))
				//std::cout << seg.source();
				constrained_region.push_back(seg.source());
		}
	//If e is a ray (segment to the infinity) and the source is inside the boundary
	}else if( e->has_source() ){
		if(check_inside(e->source()->point(), boundary_points, K())) { //if source is inside the border
			const CGAL::Object rdual = vd.dual().dual(e->dual());
			Ray_2 ray = CGAL::object_cast<K::Ray_2>(rdual);
			//Check intersection with the boundary
			for(auto borderseg : boundary){
				const auto res = CGAL::intersection(ray, borderseg);
				if (res)
				{
					const auto pt = boost::get<Point>(&*res);
					//std::cout<<"sou"<<e->source()->point()<<" ";
					//std::cout<<" inf"<<*pt<<" b"<<borderseg.target()<<" ";
					//As the source is inside the boundary
					//we stored the source of the ray, the intersection with the bondary and the target of the boundary
					constrained_region.push_back(e->source()->point());
					constrained_region.push_back(*pt);
					constrained_region.push_back(borderseg.target());
					break;
				}
			}			
		}
	//if e is a ray (segment to the infinity) and the target is inside the boundary
	}else{
		if(check_inside(e->target()->point(), boundary_points, K())){
			const CGAL::Object rdual = vd.dual().dual(e->dual());
			Ray_2 ray = CGAL::object_cast<K::Ray_2>(rdual);
			//Check intersection with the boundary
			for(auto borderseg : boundary){
				const auto res = CGAL::intersection(ray, borderseg);
				if (res)
				{
					const auto pt = boost::get<Point>(&*res);
					//we just store the intersection with the boundary
					constrained_region.push_back(*pt);
					//std::cout<<" inf"<<*pt<<" ";
					break;
				}
			}
		}
	}
//	std::cout << std::endl;
	return constrained_region;
}

std::vector<std::vector<Point>> cut_voronoi(VD &vd, std::vector<Segment_2> &boundary_segments, std::vector<Point> &boundary){
	std::vector<std::vector<Point>> voronoi_mesh;
	for(VD::Face_iterator fit = vd.faces_begin(); fit!=vd.faces_end();++fit){
		VD::Face f = *fit;
		VD::Ccb_halfedge_circulator ec_start = f.ccb(); //Start the circulator at the first edge of face
		VD::Ccb_halfedge_circulator ec = ec_start; 
		std::vector<Point> face;
		//std::cout<<"Region: "<<std::endl;
		do { //Loop through the edges of the face
        	std::vector<Point> boundary_vertices = cut_halfedge(vd, ec, boundary_segments, boundary); //Cut the region
			face.insert(face.end(), boundary_vertices.begin(), boundary_vertices.end());
      	}while ( ++ec != ec_start ); 
		// Remove colliniear points, this is to avoid the problem of degenerate faces do to the adittion of boundary points
		for (size_t i = 0; i < face.size(); i++)
		{
			Point p1 = face[i];
			Point p2 = face[(i+1)%face.size()];
			Point p3 = face[(i+2)%face.size()];
			if(CGAL::collinear(p1,p2,p3)){ 
				face.erase(face.begin() + (i+1)%face.size());
			}
		}
		if(face.size() > 0)
			voronoi_mesh.push_back(face);
	}
	return voronoi_mesh;
}

int main(int argc, char **argv) {


/*
	//Read nodes and output file name
	std::vector<Point> points = read_nodes_from_file(argv[1]); 
	std::string output_file = std::string(argv[2]);


	//std::cout<<"Read "<<points.size()<<" points from file "<<argv[1]<<std::endl;
	unsigned seed = 138;

	//Shuffle points
	std::shuffle( points.begin(), points.end(), std::default_random_engine(seed));
	std::vector<Point> boundary = { Point(0.0,0.0), Point(1.0,0.0), Point(1.0,1.0), Point(-1.0,1.0), Point(-1.0,-1.0), Point(0.0,-1.0) };
	//std::vector<Point> boundary = { Point(0.0,0.0), Point(10000.0,0.0), Point(10000.0,10000.0), Point(0.0,10000.0) };
	std::vector<Segment_2> boundary_segments;
	for(int i = 0; i < boundary.size(); i++){
		Segment_2 s(boundary[i], boundary[ (i + 1) % boundary.size()]);
		boundary_segments.push_back(s);
	}


	std::default_random_engine myRandomEngine(seed);
    // Initialize a uniform_int_distribution to produce values between -10 and 10
    //std::uniform_int_distribution<int> myUnifIntDist(-12, 12);
	std::uniform_int_distribution<int> myUnifIntDist(0, 10000000);
    
	std::vector<Wpoint> weighted_points(points.size());
	for(size_t i = 0; i < points.size(); i++){
		weighted_points[i] = Wpoint(points[i], myUnifIntDist(myRandomEngine));
	}
	
	//Regular triangulation generation
	Regular Rtr(weighted_points.begin(), weighted_points.end());

	//std::cout<<"number triangles "<< Rtr.number_of_faces()<<std::endl;
	CGAL::draw(Rtr);

	std::vector<Point> non_hidden_points;
  	for (Regular::Vertex_handle v : Rtr.finite_vertex_handles())
  		non_hidden_points.push_back(Point(v->point()));	
	printRegularOFF(Rtr, non_hidden_points, output_file);

	Delaunay T(non_hidden_points.begin(), non_hidden_points.end());
	printDelaunayOFF(T, non_hidden_points, output_file);
	CGAL::draw(T);
	
	ConstrainedTriangulation Tr_from_regular;
	Tr_from_regular.insert(non_hidden_points.begin(), non_hidden_points.end());
	printRandomOFF(Tr_from_regular, non_hidden_points, output_file);
	CGAL::draw(Tr_from_regular);


	// Random Constrained Triangulation generation
	ConstrainedTriangulation Tr;
	auto tb_randomTr = std::chrono::high_resolution_clock::now();
	Tr.insert(points.begin(),points.end());
	Tr.insert_constraint(boundary.begin(), boundary.end(), true);
	auto te_randomTr = std::chrono::high_resolution_clock::now();
	uint t_randomTr = std::chrono::duration_cast<std::chrono::milliseconds>(te_randomTr - tb_randomTr).count(); 
	//printConstrainedTriangulation(Tr, boundary, output_file);
	//CGAL::draw(Tr);


	//Delaunay triangulation generation without constraints
	Delaunay dt2;
	auto tb_delaunayTR = std::chrono::high_resolution_clock::now();
	dt2.insert(points.begin(),points.end());
	auto te_delaunayTR = std::chrono::high_resolution_clock::now();
	uint t_delaunayTR = std::chrono::duration_cast<std::chrono::milliseconds>(te_delaunayTR - tb_delaunayTR).count(); 
	//CGAL::draw(dt2);


	//Constrained Delaunay triangulation generation
	ConstrainedDelaunay CDT;
	auto tb_constrainedDelaunayTR = std::chrono::high_resolution_clock::now();
	CDT.insert(points.begin(),points.end());
	CDT.insert_constraint(boundary.begin(), boundary.end(), true);
	auto te_constrainedDelaunayTR = std::chrono::high_resolution_clock::now();
	uint t_constrainedDelaunayTR = std::chrono::duration_cast<std::chrono::milliseconds>(te_constrainedDelaunayTR - tb_constrainedDelaunayTR).count();
	//printConstrainedTriangulation(CDT, boundary, output_file);
	//CGAL::draw(CDT);

	//Voronoi diagram generation with constrains
	VD vd;
	auto tb_voronoiAdaptator = std::chrono::high_resolution_clock::now();
	vd.insert(points.begin(),points.end());
	auto te_voronoiAdaptator = std::chrono::high_resolution_clock::now();
	uint t_voronoiAdaptator = std::chrono::duration_cast<std::chrono::milliseconds>(te_voronoiAdaptator - tb_voronoiAdaptator).count();
	
	//Cut Voronoi faces and store it a vector of vectors
	auto tb_cutVoronoi = std::chrono::high_resolution_clock::now();
	std::vector<std::vector<Point>> voronoi_mesh = cut_voronoi(vd, boundary_segments, boundary);
	auto te_cutVoronoi = std::chrono::high_resolution_clock::now();
	uint t_cutVoronoi = std::chrono::duration_cast<std::chrono::milliseconds>(te_cutVoronoi - tb_cutVoronoi).count();
	
	//printVoronoiOFF(voronoi_mesh, output_file);

//	for (auto region : voronoi_mesh){
//		std::cout<<"Region: ";
//		for(auto point : region){
//			std::cout<<"("<<point<<") ";
//		}
//		std::cout<<std::endl;
//	}
//	CGAL::draw(vd);

	std::cout<<"Random Constrained Triangulation: "<<t_randomTr<<" ms"<<std::endl;
	std::cout<<"Delaunay Triangulation: "<<t_delaunayTR<<" ms"<<std::endl;
	std::cout<<"Constrained Delaunay Triangulation: "<<t_constrainedDelaunayTR<<" ms"<<std::endl;
	std::cout<<"Voronoi adaptator: "<<t_voronoiAdaptator<<" ms"<<std::endl;
	std::cout<<"Cut Voronoi: "<<t_cutVoronoi<<" ms"<<std::endl;
	std::cout<<"Constrained Voronoi: "<<t_voronoiAdaptator+t_cutVoronoi<<" ms"<<std::endl;

	std::cout<<t_randomTr<<" "<<t_delaunayTR<<" "<<t_constrainedDelaunayTR<<" "<<t_voronoiAdaptator<<" "<<t_cutVoronoi<<" "<<t_voronoiAdaptator+t_cutVoronoi<<std::endl;
*/

	//Read nodes and output file 
	std::cout<<"Reading file "<<argv[1]<<std::endl;
	std::vector<Point> points = read_nodes_from_file(argv[1]); 
	
	std::string output_file = std::string(argv[2]);

	std::cout<<"Read "<<points.size()<<" points from file "<<argv[1]<<std::endl;
	std::cout<<"Generating Delaunay triangulation"<<std::endl;


	int n_points = 0;
	//Constrained Delaunay triangulation generation
	ConstrainedDelaunay CDT;
	auto tb_constrainedDelaunayTR = std::chrono::high_resolution_clock::now();
	n_points = CDT.insert(points.begin(),points.end());
	auto te_constrainedDelaunayTR = std::chrono::high_resolution_clock::now();
	uint t_constrainedDelaunayTR = std::chrono::duration_cast<std::chrono::milliseconds>(te_constrainedDelaunayTR - tb_constrainedDelaunayTR).count();
	std::cout<<"Constrained Delaunay Triangulation with "<<n_points<<std::endl;
	std::cout<<"Generated Delaunay Triangulation in "<<t_constrainedDelaunayTR<<" ms"<<std::endl;
	long long mem_peak = malloc_count_peak();
	long long mem_current = malloc_count_current();
	std::cout<<"Memory peak: "<< (long long) mem_peak<<" bytes"<<std::endl;
	std::cout<<"Memory current: "<< (long long) mem_current<<" bytes"<<std::endl;

	printDelaunayTriangulation(CDT, output_file);
	
	std::ofstream file;
    file.open(output_file + "_triangulation_info.json");
	file << "{" << std::endl;
	file << "\"n_points\": " << n_points << "," << std::endl;
	file << "\"triangulation_time\": " << t_constrainedDelaunayTR << "," << std::endl;
	file << "\"memory_usage\": " << mem_peak << "," << std::endl;
	file << "\"memory_peak\": " << mem_current << std::endl;
	file << "}" << std::endl;

	return EXIT_SUCCESS;
}

