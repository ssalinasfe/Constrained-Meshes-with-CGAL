#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/IO/Color.h>
#include <CGAL/Triangulation_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Constrained_triangulation_2.h>
#include <CGAL/draw_triangulation_2.h>
#include <CGAL/intersections.h>
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

struct ElementInfo
{
  int index;
  bool border = false;
};


typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
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

//typedefs for Cropped Voronoi
typedef CGAL::Delaunay_triangulation_2<K>  Delaunay;
typedef K::Segment_2 Segment_2;
typedef K::Point_2 Point_2;
typedef K::Iso_rectangle_2 Iso_rectangle_2;
typedef K::Segment_2 Segment_2;
typedef K::Ray_2 Ray_2;
typedef K::Line_2 Line_2;

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
}

//Print a node file based in a Triangulation_2, the file name and the number of interior faces
void print_node_file(ConstrainedTriangulation Tr, std::string name){
    std::ofstream outfile(name + ".1.node");
    outfile << Tr.number_of_vertices() <<" 2 0 1"<< std::endl;
    for(ConstrainedTriangulation::Vertex_handle v : Tr.finite_vertex_handles()){
        outfile << v->info().index<<" "<< std::setprecision(15)<< v->point().x() << " " << v->point().y() << " "<<v->info().border << std::endl;
    }
    outfile.close();
}

//Print triangle file based in a Triangulation_2, the file name and the number of interior faces
void print_ele_file(ConstrainedTriangulation Tr, std::string name, int n_faces){
    std::ofstream outfile(name + ".1.ele");
    outfile << n_faces <<" 3 0"<< std::endl;
    for(ConstrainedTriangulation::Face_handle f : Tr.finite_face_handles()){
        if(f->info().index!=-1)
            outfile << f->info().index << " " << f->vertex(0)->info().index << " " << f->vertex(1)->info().index << " " << f->vertex(2)->info().index << std::endl;
    }
    outfile.close();
}

//Print a adjacency file based in a Triangulation_2, the file name and the number of interior faces
void print_neigh_file(ConstrainedTriangulation Tr, std::string name, int n_faces){
  std::ofstream outfile(name + ".1.neigh");
  outfile << n_faces <<" 3"<< std::endl;
  for(ConstrainedTriangulation::Face_handle f : Tr.finite_face_handles()){
      if(f->info().index!=-1){
          outfile << f->info().index << " " << f->neighbor(0)->info().index << " " << f->neighbor(1)->info().index << " " << f->neighbor(2)->info().index << std::endl;
      }
  }
  outfile.close();
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

  std::cout<<"Printing mesh with "<<index<<" faces and "<< Tr.number_of_vertices() <<" vertices"<<std::endl;
  print_off(Tr, output_file, index);
  print_node_file(Tr, output_file);
  print_ele_file(Tr, output_file, index);
  print_neigh_file(Tr, output_file, index);
  std::cout<<"Printed files in "<<output_file<<std::endl;
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

  std::cout<<"Printing mesh with "<<index<<" faces and "<< Tr.number_of_vertices() <<" vertices"<<std::endl;
  print_off(Tr, output_file, index);
  print_node_file(Tr, output_file);
  print_ele_file(Tr, output_file, index);
  print_neigh_file(Tr, output_file, index);
  std::cout<<"Printed files in "<<output_file<<std::endl;
}

//A class to recover Voronoi diagram from stream.
//Rays, lines and segments are cropped to a rectangle
//so that only segments are stored
struct Cropped_voronoi_from_delaunay{
	std::list<Segment_2> m_cropped_vd;
	//Iso_rectangle_2 m_bbox;
	std::vector<Segment_2> m_segments;
	std::vector<Point_2> m_points;

	Cropped_voronoi_from_delaunay(std::vector<Point> points){
		for(int i = 0; i < points.size(); i++){
			Segment_2 s(points[i], points[ (i + 1) % points.size()]);
			m_segments.push_back(s);
		}
		m_points = std::vector<Point_2>(points.begin(), points.end()); //borrar esto y revisar como checkear puinto en poligono usando segmento
	}

	void operator<<(const Ray_2& ray){ 
		//std::cout<<"Ray: "<<ray<<std::endl;
		for(auto borderseg : m_segments){
			const auto res = CGAL::intersection(ray, borderseg);
			if (res)
			{
				const auto pt = boost::get<Point>(&*res);
				std::cout<<"Ray: "<< *pt << std::endl;
				if(CGAL::area(borderseg.source(), borderseg.target(), ray.source()) > 0){ //if source is inside the border
					m_cropped_vd.push_back(Segment_2(*pt, ray.source()));
					std::cout<<"Added ray "<< Segment_2(*pt, ray.source()) << std::endl;
					break;
				}
			}
		}
	}
	void operator<<(const Line_2& line){ 
		std::cout<<"Line: "<<line<<std::endl;
		exit(0);
	}

	void operator<<(const Segment_2& seg){ 
		//std::cout<<"Segment: "<<seg<<std::endl;
		bool flag = false;
		for(auto borderseg : m_segments){
			const auto res = CGAL::intersection(seg, borderseg);
			if (res)
			{
				const auto pt = boost::get<Point>(&*res);
				//std::cout<<"seg: "<< *pt << std::endl;
				if(CGAL::area(borderseg.source(), borderseg.target(), seg.source()) > 0){ //if source is inside the border
					m_cropped_vd.push_back(Segment_2(*pt, seg.source()));
				}else{ //target is inside the border
					m_cropped_vd.push_back(Segment_2(*pt, seg.target()));
				}
				flag = true;
				break;
			}
		}
		if(!flag){ //if the segment does not intersect with the border, check if is inside the border
			if(check_inside(seg.source(), m_points, K()) && check_inside(seg.target(), m_points, K()))
				m_cropped_vd.push_back(seg);
		}
	}
};

int main(int argc, char **argv) {

  //Read nodes and output file name
  std::vector<Point> points = read_nodes_from_file(argv[1]); 
  std::string output_file = std::string(argv[2]);


  std::cout<<"Read "<<points.size()<<" points from file "<<argv[1]<<std::endl;
  unsigned seed = 138;

  //Shuffle points
  std::shuffle( points.begin(), points.end(), std::default_random_engine(seed));
  std::vector<Point> boundary = { Point(0.0,0.0), Point(1.0,0.0), Point(1.0,1.0), Point(-1.0,1.0), Point(-1.0,-1.0), Point(0.0,-1.0) };
  


  // Random Constrained Triangulation generation
  ConstrainedTriangulation Tr;
  Tr.insert(points.begin(),points.end());
  Tr.insert_constraint(boundary.begin(), boundary.end(), true);
  //printConstrainedTriangulation(Tr, boundary, output_file);
  //CGAL::draw(Tr);

  //Voronoi diagram generatio based on Delaunay triangulation
  Delaunay dt2;
  dt2.insert(points.begin(),points.end());
  Cropped_voronoi_from_delaunay vor(boundary);
  //extract the cropped Voronoi diagram
  dt2.draw_dual(vor);
  //for(auto seg : vor.m_cropped_vd){
//	std::cout<<"Segment: "<<seg<<std::endl;
  //}
  //print the cropped Voronoi diagram as segments
  //std::copy(vor.m_cropped_vd.begin(),vor.m_cropped_vd.end(), std::ostream_iterator<Segment_2>(std::cout,"\n"));


  //Constrained Delaunay triangulation generation
  ConstrainedDelaunay CDT;
  CDT.insert(points.begin(),points.end());
  CDT.insert_constraint(boundary.begin(), boundary.end(), true);
 // printConstrainedTriangulation(CDT, boundary, output_file);
  //CGAL::draw(CDT);



  return EXIT_SUCCESS;
}